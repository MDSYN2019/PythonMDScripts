# E_{edge} = 2 \times sigma \times L \sqrt(L - (LC/2)^2)

import numpy as np
import matplotlib.pyplot as plt
import scipy

class curvatureEnergetics:
	def __self__(self, L, C, Csp, sigma, A, kappa):
		self.L = L
		self.C = C
		self.Csp = Csp
		self.A = A 
		self.sigma = sigma
		self.kappa = kappa
	def E_edge(self):
		E = (2 * self.sigma * np.pi * self.L) * np.sqrt(1 - (np.pow((self.L * self.C)/2),2))
		return E
	def E_bend(self):
		E = (0.5 * self.A * self.kappa) * np.pow(((self.L * self.C) - (self.L*self.Csp)), 2)
		return E
	
	

