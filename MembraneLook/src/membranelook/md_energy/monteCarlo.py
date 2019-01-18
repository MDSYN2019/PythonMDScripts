import numpy as np
import scipy as sp
import scipy.constants as const
import ipdb
import copy

# own classes
import membrane
import energy
import parameters
import plotMembrane


def monteCarlo(m, num_rnd, num_rnd1, tempSA, piChange, e0, c0, jobNo=00):

    np.random.seed()  # initialize by random

    print("   ", jobNo, energy.energy(m, c0))

    i = 0
    e1 = e0
    m1 = copy.deepcopy(m)
    m2 = copy.deepcopy(m)
    # set energy for plots
    stopIteration = False

    while (not stopIteration):
        i = i + 1
        # select n random elements
        rndAll = np.arange(1, m.psi.size - 1)
        np.random.shuffle(rndAll)
        rndIdx = rndAll[:num_rnd]
        rndIdx2 = rndAll[num_rnd:num_rnd + num_rnd1]
        rndIdx = np.int8(rndIdx)
        rndIdx2 = np.int8(rndIdx2)
        m2.psi[rndIdx] = piChange * (2 * np.random.random(num_rnd) - 1)
        localVal = energy.findLocalMin(m2, rndIdx2, c0)
        m2.psi[rndIdx2] = localVal
        e2 = energy.energy(m2, c0)

        # print(rndIdx, " - ", rndIdx2)
        # print(m2.psi[rndIdx])

#            decision process MC
        if e2 <= e0:
            m1 = copy.deepcopy(m2)
            m = copy.deepcopy(m2)
            e1 = e2.copy()
            e0 = e2.copy()
            noChange = 1
        else:
            paramSA = np.exp((e1 - e2) / tempSA)
            #print(e0, e1, e2)
            if paramSA > np.random.random():
                m1 = copy.deepcopy(m2)
                e1 = e2.copy()
                print("-", jobNo, "-")
            else:
                m2 = copy.deepcopy(m1)

        if i > 10:
            stopIteration = True

        # print("e1 = ", e1)

    return (e1, m)
