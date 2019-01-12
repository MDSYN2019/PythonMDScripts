# general import
import numpy as np
import scipy as sp
import scipy.constants as const
import ipdb
import copy
import multiprocessing

# own classes
import membrane
import energy
import parameters
import plotMembrane
import monteCarlo2
import variables

var = variables.Variables()
var.bend_mod = 9  # kbT
var.lipid_area = 72.4 * 1e-2  # angstrom2
var.stretch_mod = 23  # kbTnm-2
var.c0 = 0  # nm-1
nano_r = 2
var.r1 = nano_r
var.r2 = nano_r
nano_dist = 10
var.x1 = nano_dist / 2 + var.r1
var.x2 = -nano_dist / 2 - var.r1

# testing shape
startZ = var.hc0 + var.r1
endZ = 0
elementNum = 50  # numbers of elements
membraneLength = 20.0  # membrane length
elementPsi = np.zeros(elementNum)
elementPsi[0] = 0
elementPsi[-1] = 0
elementLength = membraneLength / elementNum

i = 0

m = membrane.Membrane(elementPsi, elementLength, startZ, endZ)
e0 = energy.energy(m, var)

# define parameters
num_rnd = 5
num_rnd1 = 10
i = 0
tempSA = 1e2
piChange = 1 * const.pi / 180
numProcessors = 1

noChange = 1
e1 = e0

# set energy for plots
eTotal = np.zeros(100)
eMin = np.zeros(100)
step = 0
noChange = 0
steps = np.zeros(100)

stopIteration = False

m1 = copy.deepcopy(m)

while (not stopIteration):
    step = step + 1

    energyMC = []
    membraneMC = []
    tasksMC = []

    # pool = multiprocessing.Pool(numProcessors)

    # for i in range(numProcessors):

    resultsMC = monteCarlo2.monteCarlo(copy.deepcopy(
        m1), num_rnd, num_rnd1, tempSA, piChange, e0, var, i)

    energyMC = resultsMC[0]
    membraneMC = resultsMC[1]


#            decision process MC

    emin = energyMC

    # ipdb.set_trace()

    if emin < e0:
        m = copy.deepcopy(membraneMC)
        m1 = copy.deepcopy(membraneMC)
        e0 = emin
        noChange = 0
    else:
        noChange = noChange + 1
        paramSA = np.exp((emin - e0) / tempSA)
        # print(e0, e1, e2)
        if paramSA > np.random.random():
            m1 = copy.deepcopy(membraneMC)
            e1 = emin

    if noChange > 10:
        tempSA = tempSA * 0.5
        piChange = piChange * 0.8
        m1 = copy.deepcopy(m)
        noChange = 0

    print(step, noChange, tempSA, piChange, e0, energyMC)
    eTotal = np.append(np.delete(eTotal, 0), energyMC)
    eMin = np.append(np.delete(eMin, 0), e0)
    steps = np.append(np.delete(steps, 0), step)

    if i > 10:
        stopIteration = True

    if (step % 10) == 0:
        # ipdb.set_trace()
        plotMembrane.plotMembrane(m, var, steps, eMin, eTotal, 0, 0)

    # if (step % 100) == 0:
        # ipdb.set_trace()

    if (np.sum(np.abs(eTotal - eTotal.mean())) / e0 < 1e-5):
        stopIteration = True
