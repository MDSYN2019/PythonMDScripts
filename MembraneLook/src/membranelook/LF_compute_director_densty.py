"""
This program computes the extent of the orderphobic effect around each lipid tailgroup beads, of C4#,
where # is A or B. For the cholesterol equivalents, we use the C2 bead as the equivalents.
"""

# Module imports
from tqdm import tqdm
import argparse
import MDAnalysis
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as ml
import seaborn as sns

from MDAnalysis.analysis.leaflet import LeafletFinder
from numpy import (array, dot, arccos, clip)
from numpy.linalg import norm

def calc_vec(x1, x2, y1, y2, z1, z2):
    vec = np.array([x2 - x1, y2 - y1, z2 - z1])
    return vec

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def Distance_points(x, x_ref, y, y_ref):
    dist = np.sqrt(np.square(x-x_ref) + np.square(y-y_ref))
    return dist

def Calculate_director_density():
    pass

parser = argparse.ArgumentParser()
parser.add_argument("--echo", help = "echo the string you use here")
parser.add_argument("--square", help = "display a square of a given number", type = int)

args = parser.parse_args()
if args.square:
    print(args.square**2)
elif args.echo:
    print(args.echo)

u = MDAnalysis.Universe("min_bilayer_NP.gro", "traj_comp.xtc")

L = LeafletFinder(u,'name PO4')
leaflet0 = L.groups(0)
leaflet1 = L.groups(1)

l1_id = leaflet0.resids
l2_id = leaflet1.resids

for ids in l1_id:
    residue_id = "resid {}".format(ids)
    head_position = u.select_atoms(residue_id).select_atoms('name PO4')
    tail_position1 = u.select_atoms(residue_id).select_atoms('name C3A or name D3A')
    tail_position2 = u.select_atoms(residue_id).select_atoms('name C3B or name D3B')
    #print(head_position.positions, tail_position1.positions, tail_position2.positions) 
    print(calc_vec(tail_position1.positions[0][0], head_position.positions[0][0], tail_position1.positions[0][1], head_position.positions[0][1], tail_position1.positions[0][2], head_position.positions[0][2]), ids)
    #print (A)
for ids in l2_id:
    pass

atoms_tail = u.select_atoms('name C2A or name C2B or name D2A or name D2B')

Rgyr = []
NP = u.select_atoms("resname NP")
NP_positions = NP.center_of_geometry()
#print (NP_positions)

NP_array = []
lipid_array = []

X = np.array([])
Y = np.array([])
Z = np.array([])

L = LeafletFinder(u,'name PO4')
leaflet0 = L.groups(0)
leaflet1 = L.groups(1)#
list(leaflet0.select_atoms('name PO4').positions) # Head cartesian coordinate                                                                                            
list(leaflet0.select_atoms('name C3A or name D3A').positions) # Tail Cartesian coordinate for the A chain                                                                


#for ts in u.trajectory[300:400:1]:
#    L = LeafletFinder(u,'name PO4')
#    leaflet0 = L.groups(0)
#    leaflet1 = L.groups(1)
    

#np.save('NP.npy', NP_array)
#np.save('dat.npy', lipid_array)

#NP_dat = np.load("NP.npy")
#dat = np.load("dat.npy")
#avBin = []

