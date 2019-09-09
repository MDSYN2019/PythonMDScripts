"""
-------------------------------------------------------
Using optimization routines from scipy and statsmodels
------------------------------------------------------

"""

import scipy.linalg as la
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
import pandas as pd


def f(x):
	return x**3 -3*x + 1

x = np.linspace(-3,3,100)
plt.axhline(0, c = 'red')
plt.plot(x,f(x))

