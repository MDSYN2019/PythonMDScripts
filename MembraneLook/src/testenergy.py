
import numpy as np
import matplotlib.pyplot as plt

def Energy(LC,L,CSP):
    E = (LC - (L*CSP))**2 + ((L*0.225)* np.sqrt(1- (LC/2)**2))
    return E

A = np.linspace(-2,2,100)
B = []
for i in A:
    B.append(Energy(i,2.8, 0.54))

plt.scatter(A,B)
plt.show()
