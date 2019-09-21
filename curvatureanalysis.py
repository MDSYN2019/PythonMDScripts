"""

Author: Sang Young Noh 

Date (last updated): 21/09/2019

Notes: Curvature analysis to interpret which energetic component is responsible for 
       the aggregation of nanoparticles in the mixed bilayers

Useful Links for this work: 
------------
https://buildmedia.readthedocs.org/media/pdf/pybilt/latest/pybilt.pdf

https://github.com/LoLab-VU/PyBILT/blob/master/jupyter_notebooks/bilayer_analyzer.ipynb


"""

import MDAnalysis as mda # MDAnalysis

import tqdm # tqdm 
from tqdm import tqdm # tqdm

# pybilt modules

from pybilt.bilayer_analyzer import BilayerAnalyzer
from pybilt.bilayer_analyzer import print_valid_analyses
from pybilt.com_trajectory import COMTraj
from pybilt.lipid_grid.lipid_grid import grid_curvature
from pybilt.lipid_grid.lipid_grid import LipidGrid2d
from pybilt.lipid_grid.lipid_grid_opt import LipidGrids

sel_string = "resname DPPC or resname DUPC"
ba = BilayerAnalyzer(
    structure='../replicate.gro',
    trajectory='../md.xtc',
    selection=sel_string,
)
ba.add_analysis('ndcorr ndcorr_1')
#ba.run_analysis()

#ndcorr_1Dat = ba.get_analysis_data('ndcorr_1')


# Separate Curvature analysis

u = mda.Universe('../replicate.gro','../md.xtc')
lipids = u.select_atoms("resname DPPC DUPC")
lipidindices = lipids.indices
G = [u.select_atoms("resid {}".format(index)) for index in lipidindices]
COM = [i.center_of_mass() for i in G]

bilayer = u.select_atoms("resname DPPC DUPC")
com_traj = COMTraj.COMTraj(u.trajectory, bilayer)
#for ts in tqdm(u.trajectory):
#	print (COM)
	#print(lipidindices.center_of_mass())

	#A = ts.positions
	#print (A.shape)

for frame_idx in tqdm(range(com_traj.nframes)):
	lgs = LipidGrids(com_traj.frame[frame_idx], com_traj.leaflets, [0,1], area_per_bin=5.0)
	apl_vals = lgs.curvature()
	print (apl_vals)
#    print("frame: ",frame_idx, " composite APL: ",apl_vals[0], " partial for DPPC: ", apl_vals[1]['DPPC'][0], " DUPC: ", apl_vals[1]['DUPC'][0])
