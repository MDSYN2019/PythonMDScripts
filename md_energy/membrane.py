import numpy as np
import scipy.optimize as optim
import scipy.integrate as integ
import ipdb
import matplotlib.pyplot as plt
import pprint


class Membrane:
    'generic membrane object containing membrane geometry'

    def __init__(self, psi=np.zeros(10),
                 elementLength=np.ones(10), startZ=0, endZ=0):
        self.psi = psi
        self.elementLength = elementLength * np.ones(self.psi.size)
        self.length = np.append(0, self.elementLength.cumsum())
        self.startZ = startZ
        self.endZ = endZ
        self.xcoor = self.get_X()
        self.zcoor = self.get_Z()

    def print_geom(self):
        print("ElementLength ", self.elementLength)
        print("ElementPsi ", self.psi)
        print("startZ ", self.startZ)
        print("endZ ", self.endZ)
        print("xCoor", self.get_X())
        print("zCoor", self.get_Z())
        print("zEnd", self.get_Zend())
        print("end of out")

    def get_X(self):
        dX = self.elementLength * np.cos(self.psi)
        dX = np.append([0], dX)
        return dX.cumsum()

    def get_Z(self):
        dZ = - self.elementLength * np.sin(self.psi)
        dZ = np.append([0], dZ)
        return self.startZ + dZ.cumsum()

    def get_Zend(self):
        return self.get_Z()[-1]

    def get_diffZ(self):
        return self.endZ - self.get_Zend()

    def funZ_zero(self, psiAngles):
        psi2 = np.append(np.append(self.psi[0], psiAngles[::]), self.psi[-1])
        dZ2 = - self.elementLength * np.sin(psi2)
        dZ = np.append(0, dZ2)
        Z2 = self.startZ + dZ2.cumsum()
        return (self.endZ - Z2[-1])**2

    def set_starting_shape(self):
        x0 = optim.fmin_bfgs(self.funZ_zero, np.zeros(self.psi.size - 2))
        self.psi = np.append(np.append(self.psi[0], x0[::]), self.psi[-1])

    def get_dPsi(self):
        return np.append(np.diff(self.psi), 0)

    def get_c1(self):
        return np.append(self.get_dPsi() / self.elementLength, 0)

    def getLength(self):
        return self.length[-1]
