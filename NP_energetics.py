# E_{edge} = 2 \times sigma \times L \sqrt(L - (LC/2)^2)

import numpy as np
import matplotlib.pyplot as plt
import scipy

class curvatureEnergetics:
	def __self__(self, L, C, Csp, sigma, A, kappa):
		self.L = L # Radii?
		self.C = C # Curvature 
		self.Csp = Csp # Spontaneous Curvature
		self.A = A  # Area 
		self.sigma = sigma # Line tension 
		self.kappa = kappa # Bilayer Rigidity 
	def E_edge(self):
		"""
		Equation from reference - 
		"""
		E = (2 * self.sigma * np.pi * self.L) * np.sqrt(1 - (np.pow((self.L * self.C)/2),2))
		return E
	def E_bend(self):
		"""
		Equation from reference - 
		"""
		E = (0.5 * self.A * self.kappa) * np.pow(((self.L * self.C) - (self.L*self.Csp)), 2)
		return E
	
	

