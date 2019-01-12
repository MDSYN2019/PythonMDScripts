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



def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def Distance_points(x, x_ref, y, y_ref):
    dist = np.sqrt(np.square(x-x_ref) + np.square(y-y_ref))
    return dist

def CalculatePhi6(totalCoordinatesArray):
    referenceVector = np.array([100.0, 100.0, 0.0])
    coordinate_vector = []
    for referenceCoordinates in tqdm(totalCoordinatesArray):
        val = 0
        closestArraycoordinates = []
        for hexagonCoordinates in totalCoordinatesArray:
            #if referenceCoordinates.all() != hexagonCoordinates.all(): # Ensure we do not double count the reference coordinates
               # print("ljnd")
            closestArraycoordinates.append([Distance_points(hexagonCoordinates[0], referenceCoordinates[0], hexagonCoordinates[1], referenceCoordinates[1]), hexagonCoordinates[0], hexagonCoordinates[1], hexagonCoordinates[2]]) 
            closestArraycoordinates.sort(key=lambda x: x[0])
        for coordinate in closestArraycoordinates[1:6]:
            # fix vector to that from 0 to that from the central point
            coordinate[1] = coordinate[1] - referenceCoordinates[0]
            coordinate[2] = coordinate[2] - referenceCoordinates[1]
            coordinate[3] = coordinate[3] - referenceCoordinates[2]
            c = dot(coordinate[1:], referenceVector)/norm(coordinate[1:])/norm(referenceVector)
            angle = arccos(clip(c, -1, 1))
            val += np.exp(1j* 6 * np.degrees(angle))
        referenceCoordinates = np.append(referenceCoordinates, np.abs(1/6 * val) * np.abs(1/6 * val) )
        coordinate_vector.append(referenceCoordinates.tolist())
        #print (coordinate_vector)
    return coordinate_vector

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

for ts in u.trajectory[300:400:1]:
    atoms_tail = u.select_atoms('name C2A or name C2B or name D2A or name D2B')
    NP = u.select_atoms("resname NP")
    NP_positions = NP.center_of_geometry()
    dat = CalculatePhi6(atoms_tail.positions)    
    NP_array.append(NP_positions)
    lipid_array.append(dat)
    print(NP_array)

    for line in dat:
        np.append(X,line[0]) # x 
        np.append(Y,line[1]) # y 
        np.append(Z,line[3]) # phi

np.save('NP.npy', NP_array)
np.save('dat.npy', lipid_array)

NP_dat = np.load("NP.npy")
dat = np.load("dat.npy")
avBin = []

for i in range(0,70,1):     
    phibin = [] 
    for index in range(0,len(dat)):
            for line in dat[index]:
                if Distance_points(line[0], NP_dat[index][0], line[1], NP_dat[index][1]) > i - 1 and Distance_points(line[0], NP_dat[index][0], line[1], NP_dat[index][1]) < i:
                    print (Distance_points(line[0], NP_dat[index][0], line[1], NP_dat[index][1]), i-0.5, i + 0.5, line[3]) 
                    phibin.append(line[3])
    try:
       avBin.append([i, sum(phibin)/len(phibin), np.std(phibin)/len(NP_dat)])
    except ZeroDivisionError:
       avBin.append([i,0,0])

x = []
y = []
err = [] 
for aa in avBin:
    x.append(aa[0])
    y.append(aa[1])
    err.append(aa[2])
    
np.savetxt('phi_SE.dat', np.c_[x,y,err])
#plt.plot(x, y)
#plt.plot(x, smooth(y,6), 'b-', lw=2)
#plt.plot(x, smooth(y,10), 'r-', lw=2)
#plt.plot(x, smooth(y,15), 'r-', lw=2)
#plt.show()


#dat = CalculatePhi6(atoms_tail.positions)
#np.savetxt('dat.npy', dat)

X = []
Y = []
Z = []

for line in dat[1]:
    X.append(line[0])
    Y.append(line[1])
    Z.append(line[3])
    
nx = 150
ny = 150
xi = np.linspace(min(X), max(X), nx)
yi = np.linspace(min(Y), max(Y), ny)
zi = ml.griddata(X, Y, Z, xi, yi, interp='linear')

plt.contourf(xi, yi, zi, 30)
plt.axis([0, 150, 0, 150])
plt.colorbar()

plt.scatter(X, Y, marker = 'o', c = 'b', s = 10, zorder = 10)                                                                                      
plt.show()

#f, (ax1, ax2) = plt.subplots(1, 2)
#ax1.contourf(xi, yi, zi, 30)

#ax1.axis([0, 150, 0, 150])
#ax1.colorbar()
