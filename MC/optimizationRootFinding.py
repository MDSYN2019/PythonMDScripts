import matplotlib.pyplot as plt
import numpy as np
import scipy as scipy
from scipy.interpolate import interp1d


import os
import sys
import glob
import operator as op
import itertools as it
from functools import reduce, partial
from pandas import DataFrame, Series
import matplotlib.pyplot as plt
import seaborn as sns


"""
-----------------------------
Optimization and Root Finding
-----------------------------
"""

def f(x):
	return x**2 + 3*x**2- 3

x = np.linspace(-3.1,0,100)
