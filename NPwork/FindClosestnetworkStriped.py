
"""
Code to generate the bond restraints where pyCGtools do not seem to be adequate. Not a c
criticism of pyCGtools, but that it is not perhaps designed to optimize NP structures

Author: Sang Young Noh

Date: 25/08/2019


For the CA CA bonds

CA CA         1    0.14000   392459.2 ; 7,(1986),230; BENZENE,PHE,TRP,TYR ; This is the amber aa - between aromatic AA
CA S          1    0.17500   189953.6 ; Au_cluster_ff ; This is the parameters for the amber aa - between aromatic AA and Sulfer



"""
## Module imports

import numpy as np 
import MDAnalysis
from MDAnalysis.analysis.distances import distance_array

universe = MDAnalysis.Universe("striped_CG.gro") ## Reading in the .gro file5

# Ligand attached P5 atoms
P5Atoms = universe.select_atoms("name P5")
P5ID = P5Atoms.atoms.ids

# ST atoms
STAtoms = universe.select_atoms("name ST")
STID = STAtoms.atoms.ids

# BEN ligands
PETAtoms = universe.select_atoms("resname PET")
PETID = PETAtoms.atoms.ids

# BEN ligands
SAtoms = universe.select_atoms("resname S")
SID = SAtoms.atoms.ids


P5AtomsPositionArray = P5Atoms.positions # Core atoms
STAtomsPositionArray = STAtoms.positions # Sulfur atoms
PETAtomsPositionArray = PETAtoms.positions # 
SAtomsPositionArray = SAtoms.positions 

# Print out whole index
# We follow the format of the NP_template of the previous MARTINI NP I made
# which has the following format:

"""
How to print out the atoms part of the itp file5

[ moleculetype ]
NP             1

[ atoms ]
; nr  type  resnr residue atom cgnr charge  mass
     1       P5        1   NP   P5        1    0.0000    30.9738
     2       C1        1   NP   C1        2    0.0000    12.0110
     3       C2        1   NP   C2        3    0.0000    12.0110
     4       Qa        1   NP   Qa        4    0.0000    12.0000
     5       P5        1   NP   P5        5    0.0000    30.9738
     6       C1        1   NP   C1        6    0.0000    12.0110
     7       C2        1   NP   C2        7    0.0000    12.0110
"""

BeadDict = {'C1' : 'C1', 'C2' : 'C2', 'C3' : 'C3', 'C4' : 'C4', 'Na': 'Na', 'SC1' : 'SC1', 'SC2' : 'SC2', 'SC3' : 'SC3'} # Dict to store name and type of the MARTINI beads used for making the CG NP

BeadMass = {'C1' : 30.0, 'C2' : 30.0, 'C3' : 30.0, 'C4' : 30.0, 'Na': 30.0 , 'SC1' : 30.0, 'SC2' : 30.0, 'SC3' : 30.0} # Dict to store name and mass of the MARTINI beads used for making the CG NP 

print ('[ moleculetype ]')
print ('NP	       1')


allAtoms = universe.select_atoms('all')
print ("[ atoms ]")
print ("; nr  type  resnr residue atom cgnr charge  mass")
for row in allAtoms:
	print ("     {}       {}        {}   {}   {}        {}    {}    {}".format(row.id, BeadDict[row.id], 1, 'NP', row.name, row.id, 0.000, BeadMass[row.id]))
	

# -------------------------------------------------------------
## We want to print out bond restraints in the following format:
# -------------------------------------------------------- #
# index, index2, bond type (?), distance, spring constant  #
# -------------------------------------------------------- #

# [ bonds ]
#1 785  1 0.20        5000			
#1 397	1 0.20        5000

# ----------------------------------------------------------------
##  As for angular restraints, we want them to be printed out as:
# -------------------------------------------------------- #
# index, index2, index3, bond type (?), angle, spring constant  #
# -------------------------------------------------------- #
#[ angles ]
#2 3 4 2 180 25
#6 7 8 2 180 25
#10 11 12 2 180 25
#14 15 16 2 180 25

# Some sanity checks 

assert (len(P5ID) == len(P5AtomsPositionArray))
#print (len(P5ID), len(P5AtomsPositionArray)) 

assert (len(STID) == len(STAtomsPositionArray))
#print (len(STID), len(STAtomsPositionArray)) 

assert (len(PETID) == len(PETAtomsPositionArray))
#print (len(PETID), len(PETAtomsPositionArray)) 



#print (STID, P5ID)

## Find the closest ST group to the Gold Core  
print ('[ bonds ]')
print('; i  j  func')
print ('; P5 - P5 ')

#print (P5AtomsPositionArray)
#print (P5ID)

for index, atom in enumerate(P5AtomsPositionArray):
	distarray = (distance_array(atom, P5AtomsPositionArray))
	for index2, entry in enumerate(distarray[0]):
		if index == index2:
			pass
		else:
			print (P5ID[index], P5ID[index2], 1, distarray[0][index2]/10, 5000) # Explanation required
			
print ('; ST - P5 ')
for index, atom in enumerate(STAtomsPositionArray):
		distarray = (distance_array(atom, P5AtomsPositionArray))
		print (STID[index], P5ID[np.argmin(distarray)], 1, np.amin(distarray), 5000) # Explanation required 		

print ('; PET ligands ')
for res in PETAtoms.residues:
	# First index is ST, second SC1, third SC2 and fourth SC4
	ligandResIDS = res.atoms.ids 
	#	print (ligandsResIDS)
	print (ligandResIDS[0], ligandResIDS[1], 1, 0.5000, 392459.2) #  ST - C1  
	print (ligandResIDS[1], ligandResIDS[2], 1, 0.2000, 392459.2) #  C1 - C2  
	print (ligandResIDS[2], ligandResIDS[3], 1, 0.2000, 392459.2) #  C2 - C3 
	print (ligandResIDS[3], ligandResIDS[1], 1, 0.2000, 392459.2) #  C3 - C1

print ('; S ligands ')
for res in SAtoms.residues:
	# zeroth index is ST, first RA, second R1, third R2 and fourth R3	
	ligandResIDS = res.atoms.ids 
	#	print (ligandsResIDS)
	print (ligandResIDS[0], ligandResIDS[1], 1, 0.5000, 392459.2) # ST - RA 
	print (ligandResIDS[1], ligandResIDS[2], 1, 0.2000, 392459.2) # RA - R1 
	print (ligandResIDS[2], ligandResIDS[3], 1, 0.2000, 392459.2) # R1 - R2 
	print (ligandResIDS[3], ligandResIDS[4], 1, 0.2000, 392459.2) # R2 - R3 
	print (ligandResIDS[4], ligandResIDS[1], 1, 0.2000, 392459.2) # R3 - R1

# Finding the closest ST section
#	closestSTArray = distance_array(atom, STAtomsPositionArray)
#print (closestSTArray)
## Need to find the indices of the STAtoms, so that I can link that on to the P5 atoms

## Print out the BEN ligands and allocate the bond and angle restrains

# print (PETAtoms.residues)

## Making Angular restraints for the hydrophilic ligands

print ('[ angles ]')
for res in SAtoms.residues:
	# zeroth index is ST, first RA, second R1, third R2 and fourth R3	
	ligandResIDS = res.atoms.ids 	
	print (ligandResIDS[0], ligandResIDS[1], ligandResIDS[2], 1, 120, 25) # ST - RA - R1  
	print (ligandResIDS[1], ligandResIDS[2], ligandResIDS[3], 1, 120, 25) # RA - R1 - R2 
	print (ligandResIDS[2], ligandResIDS[3], ligandResIDS[4], 1, 120, 25) # R1 - R2 - R3
	print (ligandResIDS[4], ligandResIDS[3], ligandResIDS[1], 1, 120, 25) # R2 - R3 - R1
	ligandRes = res.atoms
	#print(ligandRes)


for res in PETAtoms.residues:
	# zeroth index is ST, first RA, second R1, third R2 and fourth R3	
	ligandResIDS = res.atoms.ids
	print (ligandResIDS[1], ligandResIDS[2], ligandResIDS[3], 1, 120, 25) # ST - C1 - C2 
	print (ligandResIDS[1], ligandResIDS[2], ligandResIDS[3], 1, 120, 25) # C1 - C2 - C3 
	print (ligandResIDS[2], ligandResIDS[3], ligandResIDS[1], 1, 120, 25) # C2 -  - R3
	print (ligandResIDS[3], ligandResIDS[1], ligandResIDS[2], 1, 120, 25) # R2 - R3 - R1
	ligandRes = res.atoms
	#print(ligandRes)

	
