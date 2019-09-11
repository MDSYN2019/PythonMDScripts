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


def monteCarlo(m, num_rnd, num_rnd1, tempSA, piChange, e0, var, jobNo=00):

    np.random.seed()  # initialize by random

    m2 = copy.deepcopy(m)
    # set energy for plots

    # select n random elements
    rndAll = np.arange(1, m.psi.size - 1)
    np.random.shuffle(rndAll)
    rndIdx = rndAll[:num_rnd]
    rndIdx2 = rndAll[num_rnd:num_rnd + num_rnd1]
    rndIdx = np.int8(rndIdx)
    rndIdx.sort()
    rndIdx2 = np.int8(rndIdx2)
    rndIdx2.sort()
    m2.psi[rndIdx] = m2.psi[rndIdx] + piChange * \
        (2 * (np.random.random(num_rnd) - 0.5))

    localVal = energy.findLocalMin(m2, rndIdx2, var)
    m2.psi[rndIdx2] = localVal[1:]
    m2.startZ = localVal[0]

    # print('localVal',  localVal)
    # print('m2psi', m2.psi[rndIdx])

    # print(energy.funZ_zero(localVal, m2, rndIdx2))
    # if (abs(energy.funZ_zero(localVal, m2, rndIdx2)) > 1e-3):
    #    e2 = 1e5
    # else:
    e2 = energy.energy(m2, var)
    # ipdb.set_trace()
    # print(rndIdx, " - ", rndIdx2)
    # print(m2.psi[rndIdx])
    # print('e2', e2)

    # print("e1 = ", e1)

    return (e2, m2)
