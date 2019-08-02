"""

------------
OPTIMUS BIND 
------------

Version: 0.0.1

Description: 

Contributors: 

Contact:

"""
import sys
import os

# Biopython modules

import click
import logging
from pathlib import Path
from dotenv import find_dotenv, load_dotenv
import re

# Subprocess modules for calling on making the 

import subprocess
from subprocess import call

# Scipy Stack

import numpy as np
import pandas as pd
import scipy
import seaborn as sns
import matplotlib.pyplot as plt


def SKEMPItoPandas(SKEMPI_loc):
    '''
    Purpose:
        1. Loads SKEMPI CSV file.
        2. Calculates ddG
        3. For multiple measurements, keeps the median value
        4. Eliminates entries with mutations on both sides of the interface
    Input:
        SKEMPI_loc : Location of SKEMPI CSV file
    Output:
        SKEMPI_df : Pandas dataframe    
    '''
    import pandas as pd
    import numpy as np
    import re
	# fix this
    pd.options.mode.chained_assignment = None  # default='warn'

    # Constants
    R = 1.9872036e-3  # Ideal Gas Constant in kcal

    SKEMPI_df = pd.read_csv(SKEMPI_loc, sep=';')

    # Convert non numeric temperature comments to numeric values. Default is 298K 
    ConvertTemp = lambda x: int(re.search(r'\d+', x).group(0) or 298)
    BadTemps = SKEMPI_df.Temperature.str.isnumeric() == 0
    SKEMPI_df['Temperature'].loc[BadTemps] = SKEMPI_df['Temperature'].loc[BadTemps].map(ConvertTemp)
    SKEMPI_df['Temperature'] = pd.to_numeric(SKEMPI_df['Temperature'], errors='coerce')

    # Drop missing values
    SKEMPI_df.dropna(subset=['Affinity_wt_parsed'], inplace=True)
    SKEMPI_df.dropna(subset=['Affinity_mut_parsed'], inplace=True)

    # Calculate free energies
    SKEMPI_df['dgWT'] = -R*SKEMPI_df['Temperature']*np.log(SKEMPI_df['Affinity_wt_parsed'])
    SKEMPI_df['dgMut'] = -R*SKEMPI_df['Temperature']*np.log(SKEMPI_df['Affinity_mut_parsed'])
    SKEMPI_df['ddG'] = SKEMPI_df['dgWT']-SKEMPI_df['dgMut']

    # Create a key for unique mutations based on PDB and 
    SKEMPI_df['MutKey'] = SKEMPI_df['#Pdb']+'_'+SKEMPI_df['Mutation(s)_PDB']
    # Replace multiple measurements of the same mutation with the group mean
    # May consider grouping by experimental method as well
    SKEMPI_df['ddgMedian'] = SKEMPI_df.groupby('MutKey')['ddG'].transform('median')        
    SKEMPI_df = SKEMPI_df.drop_duplicates(subset=['MutKey', 'Temperature'], keep='first', inplace=False)

    # Flag multiple mutations in the same protein
    SKEMPI_df['NumMutations'] = SKEMPI_df['Mutation(s)_PDB'].str.count(',')+1 

    # Extract Chains and remove cross chain mutations. Chain is the second position in the mutation code
    SKEMPI_df['Prot1Chain'] = SKEMPI_df['#Pdb'].str.split('_').str[1]
    SKEMPI_df['Prot2Chain'] = SKEMPI_df['#Pdb'].str.split('_').str[2]
    SKEMPI_df['MutSplit'] = SKEMPI_df['Mutation(s)_PDB'].str.split(',')
    SKEMPI_df['MutCleanSplit'] = SKEMPI_df['Mutation(s)_cleaned'].str.split(',')

	# SYN added - Added a pdb name column to make it easier to identiy pdb when it comes to implementing
	# mutations
	
    NAME = [] 
    for pdbname in SKEMPI_df['#Pdb']:
        name = pdbname.split('_')[0]
        NAME.append(name)
    SKEMPI_df['NAME'] = NAME

    def ChainCheck(df):
        if df['NumMutations'] == 1:
            CrossChain = False
            return CrossChain
        else:
            Chain = df['MutSplit'][0][1]
            if Chain in df['Prot1Chain']:
                ChainSet = df['Prot1Chain']
            elif Chain in df['Prot2Chain']:
                ChainSet = df['Prot2Chain']
            for i in range(len(df['MutSplit'])):
                Chain = df['MutSplit'][i][1]
                if Chain in ChainSet:
                    CrossChain = False
                else:
                    CrossChain = True
                    break
        return CrossChain

    SKEMPI_df['CrossChain'] = SKEMPI_df.apply(ChainCheck, axis=1)
    SKEMPI_SingleSided = SKEMPI_df[SKEMPI_df.CrossChain == False]

    NumProteins = SKEMPI_SingleSided['#Pdb'].nunique()
    NumMutations = SKEMPI_SingleSided['#Pdb'].count()
    print("There are %s unique single sided mutations in %s proteins" % (NumMutations, NumProteins))             
    return SKEMPI_SingleSided

DataFrame = SKEMPItoPandas('skempi_v2.csv')


# -----------------------
# Calling perl for ialign
# -----------------------

# Ialign subprocess. Example input
# ialign.pl -w output 1lyl.pdb AC 12as.pdb AB

IalignPath = "/home/oohnohnoh1/Desktop/ACADEMIA/Papermaking/OPTIMUS_BIND/ialign/bin/ialign.pl"

# -------------
# Calling foldx
# -------------

# FoldX subprocess. Example input:
# foldx --command=Optimize --pdb=example.pdb

FoldxPath = "/home/oohnohnoh1/Desktop/ACADEMIA/Papermaking/OPTIMUS_BIND/FoldX/foldx"

#subprocess.call(IalignPath)

def callialign(folder, pdb1, chain1, pdb2, chain2, Ialignpath = None): # Default None for Ialignpath for now
	"""
	Purpose:
	
	Function to call ialign with subprocess

	Parameters
	----------
	folder: path/to/directory 
	    Path to directory where the WT and mutations are stored
	pdb1: pdb file name
	    The string of the pdb file
	chain1: chain name
	    Chain in pdb1 
	pdb2: pdb file name
	    The string of the pdb file
	chain2: chain name
	    Chain in pdb1 
	Ialignpath: placeholder
	    ---
	"""
	# TODO - add timeit 
	try:
		import os
		from os import path
	except ImportError:
		print ("..")
	IalignPath = "/home/oohnohnoh1/Desktop/ACADEMIA/Papermaking/OPTIMUS_BIND/ialign/bin/ialign.pl"
	print ("Running ialign...")
	outputFile = ">> {}.dat"
	try: 
		p = subprocess.Popen([IalignPath, "-w", folder, pdb1, chain1, pdb2, chain2, outputFile], stdout = subprocess.PIPE) # Ok this works
		stdout, stderr = p.communicate()
	except subprocess.CalledProcessError as e:
		print ("ERROR: Cannot read ialign properly. Please check the input file/chain/pdb files, or check the path to the ialign binary")
	else:
		print ("Running ialign sucessfully..")

		
def callfoldx(pdb, foldxpath = None): # Default None for Ialignpath for now
	"""
	Purpose:
	
	Function to call folx with subprocess
    
	Parameters
	----------
	pdb: pdb file name
	    The string of the pdb file
    foldxpath: path/to/binary
	    The path to the foldx binary
	"""
	# We need to make sure the rotabase.txt file is in the pdb folder
	try:
		import errno
		import os
		from os import path
	except ImportError:
		print ("..")
	FoldxPath = "/home/oohnohnoh1/Desktop/ACADEMIA/Papermaking/OPTIMUS_BIND/FoldX/foldx"
	print ("Running foldx ...")

	try:
		pdbstring = "--pdb={}".format(str(pdb))
		p = subprocess.Popen([FoldxPath, "--command=Optimize", pdbstring], stdout = subprocess.PIPE) # Need to check if this works -works
		stdout, stderr = p.communicate()
	except subprocess.CalledProcessError as e:
		print ("ERROR: Cannot run foldx properly. Please check the input file/chain/pdb files, or check the path to the foldx binary")             
	else:
		print ("Running foldx successfully")

	# -------------------------------------------------
	# Check if the optimized files have been produced |
	# -------------------------------------------------
	
	output_FX = "OP_{}.fxout".format(pdb.split('.')[0]) # FX name output
	output_PDB = "Optimized_{}.pdb".format(pdb.split('.')[0]) # Optimized PDB name output
	rotabase = "rotabase.txt" # rotabase.txt required for foldx

	# Exception conditionals for the rotabase file
	if os.path.isfile(rotabase) is True:
		print ("rotabase.txt is there - required for foldx")
	else:
		raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), rotabase)
	
	# Exception conditionals for the FX file
	if os.path.isfile(output_FX) is True:
		print ("{} has been produced successfully".format(output_FX))
	else:
		raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), output_FX)
	
	# Exception conditionals for the PDB file 
	if os.path.isfile(output_PDB) is True:
		print ("{} has been produced successfully".format(output_PDB))
	else:
		raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), output_PDB)

def function_to_read_ialign_profiles(outputdata, PDBname, foldPandasInput):
	"""
	Purpose:

	- IS-Score - brief explanation  
	- P-Value - ditto 
	- Z-Score - ditto 
	- Number of aligned residues - ditto 
	- Number of aligned contacts - ditto
	- RMSD - ditto
	- Seq Identity - ditto
	
	Parameters
	----------
	outputdata: pdb file name
	    Placeholder
	PDBname: path/to/binary
	    PLaceholder
	foldPandasInput: 
	    Placeholder
	"""
	
	Entries = ['IS-score', 'P-value', 'Z-score', 'Number of aligned residues', 'Number of aligned Contacts', 'RMSD', 'Seq Identity']
	ialignOutput = open(str(outputdata), "r")
	ialignLines = ialignOutput.readlines()
	datBlock = None
	output = []
	#print (ialignLines)
	# Concatenate the data block where we have the data
	for index, entry in enumerate(ialignLines):
		if 'IS-score' in entry:
			datBlock = ialignLines[index:index+4]
			break
	#print (datBlock)
	datBlock  = ','.join(datBlock)
	datBlock = datBlock.split(',')
	print (datBlock)
	for entry in datBlock:
		output.append([entry.split(' = ')[0].strip(), float(entry.split(' = ')[1])])
	print ("The WT vs mutation ialign results for {}.pdb is as follows:".format(PDBname))
	return output

		

			   
def GenerateMutations(DataFrame, PDB):
	"""
	Purpose:
	
	This function returns the mutated pdb protein files 
	from skempi_v2 database (https://life.bsc.es/pid/skempi2/). 

	Both single mutations and multiple comma separated mutations 
	are taken in to account. 

>	If there are multiple mutation indices for the same protein, 
	then this will generate multiple pdb files.

	Parameters
	----------
	DataFrame: pandas table 
	    The pandas table to read_csv
	PDB: str
	    The string of the pdb file
	"""
	try:
		from Bio.PDB.PDBIO import PDBIO
		from Bio.PDB.PDBParser import PDBParser
		from Bio.Data.IUPACData import protein_letters
		from Bio.SeqUtils.ProtParam import ProteinAnalysis
		from Bio.PDB.Polypeptide import PPBuilder # important 
		from Bio.PDB.Polypeptide import standard_aa_names # Standard amino acid names - https://biopython.org/DIST/docs/api/Bio.PDB.Polypeptide-module.html#standard_aa_names
		from Bio.PDB.Polypeptide import aa1 #  aa1 = 'ACDEFGHIKLMNPQRSTVWY'
		from Bio.PDB.Polypeptide import aa3 #  aa3 = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE',... 
	except ImportError:
		print ("ERROR: Need to check Biopython imports!")
		
	AminoAcidListDict = {} # Dictionary to assign alpabetical letters to amino acids
	for index, code in enumerate(standard_aa_names):
		AminoAcidListDict[aa1[index]] = aa3[index] # Building the mutation dictionary for each code
	parser = PDBParser(PERMISSIVE=1) # Standard PDB parser
	title = PDB.split('.')
	PDBList = set()	
	for pdb in DataFrame['#Pdb']:
		pdbname = pdb.split('_')[0]
		string = "{}.pdb".format(pdbname)
		PDBList.add(string)
		
	if PDB not in PDBList:
		raise Exception("The PDB is not in the SKEMPI list") # Not in the PDB list we expect - i.e. from the SKEMPI list 

	# Search for PDB mutations that contain the PDB string - e.g. the 1CSE mutations will have the format 1CSE_E_I
	# where it indicates the mutations were made in the 1CSE E and I chains 

    # Check that the mutation in the pandas column is comma separated or not
	MutationList = DataFrame.loc[(DataFrame['NAME'] == PDB.split('.')[0])] # This should get the PDB mutations 
	MutationList = MutationList.reset_index()
	# Make a dictionary (hash map) with the mutation name and the residue lists to change
	# This part below - WIP
	for index, entry in MutationList.iterrows():
		structure = parser.get_structure(str(title[0]),PDB) # reset structure each time 
		model = structure[0]  # Switch back to the unchanged one 
		for mut in entry['MutCleanSplit']:
			initAA, chain, loc, mutAA = re.findall('(\d+|.)', mut)
			# Check we are reading the right residue and index
			assert(model[chain][int(loc)].resname == AminoAcidListDict[str(initAA)]) # This will check that the model is the unmutated pdb
			print ("Mutating {} on index {} of chain {} to {}".format(AminoAcidListDict[str(initAA)], chain, loc , AminoAcidListDict[str(mutAA)]))
			model[chain][int(loc)].resname = AminoAcidListDict[str(mutAA)] # This command replaces the nonmutated species into the mutated one
			assert(model[chain][int(loc)].resname == AminoAcidListDict[str(mutAA)]) # This will check that the mutation was successful
		mutanttotalstring = '_'.join(entry['MutCleanSplit'])
		mutatedname = "{}_{}_{}.pdb".format(entry['#Pdb'], mutanttotalstring, index)
		io = PDBIO(structure)
		io.set_structure(model)
		io.save(mutatedname) # This should print out the name of protein, the mutaton list, and the index on the pandas file 
		print ("Produced new mutation PDB file {}".format(mutatedname)) # Printing out sign to say the pdb was produced
			   
		
def mapped_index(pdb, chain, index, basis='FASTA_index'):
	names = ["Residue", "Chain", "FASTA_index", "PDB_index"]
	map_data = pd.read_csv(f'PDBs/{pdb}.mapping', sep=r"\s*", 
                           header=None, names=names)
	map_data = map_data.loc[map_data['Chain'] == chain]
	map_data = map_data.set_index(basis)
	output = names[-2:][basis=='FASTA_index']
	return map_data[output][index]

	# Now that we have a dictionary, we need to make sure the mutations can be read whether a single
	# mutation or a comma separated mutation.
	
	#for index in range(0, len(MutationList)-1):
	#	MutationDict[MutationList['#Pdb'][index]] = MutationList['MutSplit'][index]
	# Parsing mutation output 
	# If we have multiple mutations that are comma separated,
	# then we will use this 
	#mutations = test_mutation.split(',')  # Thanks to the suggestion by Thomas. J. Card of the Optimus Prime project 
	#for mut in mutations:
	#	initAA, chain, loc, mutAA = re.findall('(\d+|.)', mut)
	

#if __name__ == '__main__':
#    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
#    logging.basicConfig(level=logging.INFO, format=log_fmt)#
#   # not used in this stub but often useful for finding various files
#    project_dir = Path(__file__).resolve().parents[2]

   # find .env automagically by walking up directories until it's found, then
   # load up the .env entries as environment variables
#    load_dotenv(find_dotenv())
#    main()


			   
#@click.command()
# fix this patchwork later
#@click.argument('input_filepath', type=click.Path(exists=True))
#@click.argument('output_filepath', type=click.Path())
#def main(): #main(input_filepath, output_filepath):
#    """ Runs data processing scripts to turn raw data from (../raw) into
#        cleaned data ready to be analyzed (saved in ../processed).
#    """
#    input_filepath = 'data/raw/'
#    SKEMPItoPandas('skempi_v2.csv')  # GENERALIZE!
#
#    logger = logging.getLogger(__name__)
#    logger.info('making final data set from raw data')
	

