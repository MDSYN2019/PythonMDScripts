import numpy as np
import matplotlib.pyplot as plt

LC_1 = 0.9996 # No units required
LC_1_5 = 0.9990  # No units required 
LC_2_0 = 0.7600 # No units required 

L_1 = 2.8 # units -> nanometers
L_1_5 = 3.3 # units -> nanometers
L_2 = 3.8 #

Co = 0.47 # Goes from 0.47 to 0.54, units of nm^-1
kappa = 5.4
# TODO

def Energy(LC,L,CSP):
    E = (2*np.pi*(L**2)*kappa)(LC - (L*CSP))**2 + ((L*0.225)* np.sqrt(1- (LC/2)**2))
    return E

def phispont(phi, beta):
	CC = 0.47 - (beta*phi)
	return CC 
	
A = np.linspace(0,1,100)
B = []
for i in A:
	new = phispont(i, 1.0)
	B.append(Energy(0.9996, 2.8, new))
	
plt.scatter(A,B)
plt.show()
