import numpy as np
import MDAnalysis
from MDAnalysis.analysis.distances import distance_array
from MDAnalysis.analysis.distances import between
#from MDAnalysis.tests.datafiles import PDB

universe = MDAnalysis.Universe("minimised.gro", "nvt_SMD_hydrophobic.trr") ## Reading in the .gro file
protein = universe.select_atoms("not resname DPPC W ION NP")
NP = universe.select_atoms("resname NP")

dist = 25

for dist in range(1, 70, 1):
	
	for ts in universe.trajectory:
		distance = abs(protein.center_of_mass()[0] - NP.center_of_mass()[0])
		distance = round(distance)
		#print(distance)

		#print (protein.center_of_mass(), NP.center_of_mass() )
		if distance == dist:
			print(dist, distance)
			system = universe.select_atoms("all")
			system.write("system_{}.gro".format(distance))
			break
