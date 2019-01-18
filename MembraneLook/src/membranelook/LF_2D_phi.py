import numpy as np
import matplotlib.pyplot as plt
from math import isnan
import matplotlib.mlab as ml
import seaborn as sns

NP_coord = np.load('NP_coordinates.npy')
psi = np.load('psi6s_upper_tail.npy')
coordinate = np.load('coordinate_upper_tail.npy')

psi_lower = np.load('psi6s_upper_tail.npy')
coordinate_lower = np.load('coordinate_upper_tail.npy')


def distance(x1,x2,y1,y2,z1,z2):
    dist = np.sqrt(np.power((x1-x2),2) + np.power((y1-y2),2) + np.power((z1-z2),2))
    return dist

Total_sim = {}
for index, NP_i in enumerate(NP_coord):
    result = {}
    NP_x = NP_coord[index][0]
    NP_y = NP_coord[index][1]
    NP_z = NP_coord[index][2]
    psi_here = psi[index] 
    coordinate_here = coordinate[index]
    min_r = 0.0
    max_r = 0.0 
    for r in range(2,200):
        min_r = r - 2.0 
        max_r = r + 2.0 
        ans = []
        for ind in range(0,len(coordinate_here)):
            count = 0 
            coordinate_x = coordinate_here[ind][0]
            coordinate_y = coordinate_here[ind][1]
            coordinate_z = coordinate_here[ind][2]
            NP_coordinate_dist = distance(NP_x, coordinate_x, NP_y, coordinate_y, NP_z, coordinate_z)
            if NP_coordinate_dist <= max_r:
               ans.append(np.abs(psi_here[ind]))
        result[r] = ans
    Total_sim[index] = result     

for snapshot in Total_sim.keys():
    new_result = {}
    clean_dict = {k: Total_sim[snapshot][k] for k in Total_sim[snapshot] if not len(Total_sim[snapshot][k]) == 0}
    #print (clean_dict, snapshot)
    X = [key for key in clean_dict.keys()]
    Y = [np.mean(clean_dict[key]) for key in clean_dict.keys()]
    plt.plot(X, Y, label = 'Snapshot {}'.format(snapshot))

#psi = np.load('psi6s_lower_tail.npy')
#coordinate = np.load('coordinate_lower_tail.npy')
#print (psi)

td_psi_upper = psi[1]
td_coordinate_upper = coordinate[1]

td_psi_lower = psi_lower[1]
td_coordinate_lower = coordinate_lower[1]

td_psi = np.concatenate((td_psi_upper,td_psi_lower))
td_coordinate = np.concatenate((td_coordinate_upper,td_coordinate_lower))
NP_coordinate = NP_coord[1]

x = td_coordinate[:,0] 
y = td_coordinate[:,1] 
np_x = NP_coordinate[0] 
np_y = NP_coordinate[1] 

X = list(td_coordinate[:,0])
Y = list(td_coordinate[:,1])
Z = list(np.power(np.abs(td_psi), 2))
nx = 150
ny = 150 
xi = np.linspace(min(X), max(X), nx)
yi = np.linspace(min(Y), max(Y), ny)
zi = ml.griddata(X, Y, Z, xi, yi, interp='linear')

plt.contour(xi, yi, zi, 15, linewidths = 0.2, colors = 'k')
plt.pcolormesh(xi, yi, zi)
plt.colorbar()
plt.title('2D order parameter')
plt.ylabel('ylabel')
plt.xlabel('xlabel')
plt.axis([0, 150, 0, 150])
#plt.legend()
plt.scatter(x, y, marker = 'o', c = 'b', s = 10)    
plt.scatter(np_x, np_y, marker = 'o', c = 'r', s = 5000)    

plt.show()
