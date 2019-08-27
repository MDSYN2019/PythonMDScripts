
"""
Code to generate the bond restraints where pyCGtools do not seem to be adequate. Not a c
criticism of pyCGtools, but that it is not perhaps designed to optimize NP structures

Author: Sang Young Noh

Date: 25/08/2019

"""

## Module imports

import MDAnalysis
from MDAnalysis.analysis.distances import distance_array

universe = MDAnalysis.Universe("out.gro") ## Reading in the .gro file5
P5Atoms = universe.select_atoms("name P5")
P5ID = P5Atoms.atoms.ids

STAtoms = universe.select_atoms("name ST")
STID = STAtoms.atoms.ids

# BEN ligands
BENAtoms = universe.select_atoms("resname BEN")
BENID = BENAtoms.atoms.ids


P5AtomsPositionArray = P5Atoms.positions
STAtomsPositionArray = STAtoms.positions
BENAtomsPositionArray = BENAtoms.positions 

## We want to print out bond restraints in the following format:
# -------------------------------------------------------- #
# index, index2, bond type (?), distance, spring constant  #
# -------------------------------------------------------- #

# [ bonds ]
#1 785  1 0.20        5000			
#1 397	1 0.20        5000

# Check that the lengths are identical

assert (len(P5ID) == len(P5AtomsPositionArray))
print (len(P5ID), len(P5AtomsPositionArray)) 

assert (len(STID) == len(STAtomsPositionArray))
print (len(STID), len(STAtomsPositionArray)) 


assert (len(BENID) == len(BENAtomsPositionArray))
print (len(BENID), len(BENAtomsPositionArray)) 


for index, atom in enumerate(P5AtomsPositionArray):
	distarray = (distance_array(atom, P5AtomsPositionArray))
#	print (distarray[0])
	for index2, entry in enumerate(distarray[0]):
		if index == index2:
			pass
		else:
			print (index+1, index2+1, 1, distarray[0][index2], 5000)

	# Finding the closest ST section
#	closestSTArray = distance_array(atom, STAtomsPositionArray)
	#print (closestSTArray)
	## Need to find the indices of the STAtoms, so that I can link that on to the P5 atoms

## Print out the BEN ligands and allocate the bond and angle restrains

print (BENAtoms.residues)

for res in BENAtoms.residues:
	# First index is ST, second SC1, third SC2 and fourth SC4
	ligandResIDS = res.atoms.ids 
	#print (ligandsResIDS)
#	print (ligandResIDS[0], ligandResIDS[1], 1, 0.470, 1250)
#	print (ligandResIDS[1], ligandResIDS[2], 1, 0.470, 1250)
#	print (ligandResIDS[2], ligandResIDS[3], 1, 0.470, 1250)
		

#print (BENID, P5ID, STID)

print (len(P5ID), len(P5AtomsPositionArray)) 

print (len(STID), len(STAtomsPositionArray)) 

print (len(BENID), len(BENAtomsPositionArray)) 
