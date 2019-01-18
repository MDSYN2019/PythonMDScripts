import numpy as np
import matplotlib.pyplot as plt
from math import isnan

NP_coord = np.load('NP_coordinates.npy')
psi = np.load('psi6s_lower_tail.npy')
coordinate = np.load('coordinate_lower_tail.npy')

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
            if NP_coordinate_dist <= max_r and NP_coordinate_dist >= min_r:
               ans.append(np.abs(psi_here[ind]))
        result[r] = ans
    Total_sim[index] = result     

for snapshot in Total_sim.keys():
    new_result = {}
    clean_dict = {k: Total_sim[snapshot][k] for k in Total_sim[snapshot] if not len(Total_sim[snapshot][k]) == 0}
    #print (clean_dict, snapshot)
    X = [key for key in clean_dict.keys()]
    Y = [np.mean(clean_dict[key]) for key in clean_dict.keys()]
    Y_std = [np.std(clean_dict[key]) for key in clean_dict.keys()]
    Y_std_E = [np.std(clean_dict[key])/np.sqrt(len(NP_coord)) for key in clean_dict.keys()]    
    #plt.plot(X,Y, label = 'Snapshot {}'.format(snapshot))

np.savetxt('phi_SE.txt', np.c_[X,Y,Y_std_E])

plt.title('2D order parameter')
plt.ylabel('ylabel')
plt.xlabel('xlabel')
#plt.axis([0, 60, 0, 1.0])
plt.legend()
plt.show()
