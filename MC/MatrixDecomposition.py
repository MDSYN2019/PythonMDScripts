import os
import sys
import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
plt.style.use('ggplot')

np.set_printoptions(suppress=True)

import scipy.linalg as la

A = np.array([[1,3,4], [2,1,3], [4,1,2]])
L = np.array([[1,0,0], [2,1,0], [4,11/5, 1]])
U = np.array([[1,3,4], [0, -5, -5], [0, 0,-3]])

print(L.dot(U))


