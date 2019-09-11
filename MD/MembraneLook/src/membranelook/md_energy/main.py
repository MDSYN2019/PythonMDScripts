# general import
import numpy as np
import scipy as sp
import scipy.constants as const
import ipdb
import copy
import multiprocessing
import pickle

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
elementNum = 100  # numbers of elements
membraneLength = 30.0  # membrane length
elementPsi = np.zeros(elementNum)
elementPsi[0] = 0
elementPsi[-1] = 0
elementLength = membraneLength / elementNum

i = 0

# m = membrane.Membrane(elementPsi, elementLength, startZ, endZ)

# define parameters
num_rnd = 15
num_rnd1 = 15
i = 0
tempSA = 1e2
piChange = 5 * const.pi / 180
numProcessors = 3

noChange = 1

# set energy for plots
eTotal = np.zeros(100)
eMin = np.zeros(100)
step = 0
noChange = 0
steps = np.zeros(100)

stopIteration = False


for nano_r in np.array([1.5, 2.5, 3.5, 4.5]):
    for var.c0 in np.array([-0.2]):  # ([-0.1, -0.2, -0.3]):
        for var.stretch_mod in np.array([45.0]):  # ([30.0, 45.0, 60.0]):
            for var.bend_mod in np.array([10.0]):  # ([5.0, 10.0, 15.0]):

                res = []
                resE = []
                mid = membraneLength / 2
                var.r1 = nano_r
                var.r2 = nano_r
                # testing shape
                startZ = (var.hc0 + nano_r)
                endZ = 0
                elementPsi = np.zeros(elementNum)
                elementPsi[0] = 0
                elementPsi[-1] = 0  # const.pi / 2

                for nano_dist in np.append(np.append([0, 0.5, 1, 1.5, 1.8], np.linspace(2, 5, 31)), [5.5, 6, 6.5, 7, 7.5, 8, 10, 15, 16, 17]):
                    stopIteration = False
                    step = 0
                    noChange = 0
                    steps = np.zeros(100)
                    var.x1 = nano_dist / 2 + var.r1
                    var.x2 = -nano_dist / 2 - var.r1
                    tempSA = 1e2
                    piChange = 5 * const.pi / 180
                    m = membrane.Membrane(
                        elementPsi, elementLength, startZ, endZ)
                    e0 = energy.energy(m, var)
                    m1 = copy.deepcopy(m)
                    e1 = e0
                    while (not stopIteration):
                        step = step + 1

                        energyMC = []
                        membraneMC = []
                        tasksMC = []

                        pool = multiprocessing.Pool(numProcessors)

                        for i in range(numProcessors):
                            membraneMC.append(copy.deepcopy(m1))
                            energyMC.append(e0)
                            tasksMC.append((membraneMC[i], num_rnd, num_rnd1,
                                            tempSA, piChange, energyMC[i], var, i))

                        resultsMC = [pool.apply_async(
                            monteCarlo2.monteCarlo, t) for t in tasksMC]

                        pool.close()
                        pool.join()

                        processNum = 0
                        for result in resultsMC:
                            energyMC[processNum] = result.get()[0]
                            membraneMC[processNum] = result.get()[1]
                            processNum = processNum + 1

                        #            decision process MC

                        emin = min(energyMC)
                        idxMin = energyMC.index(min(energyMC))

                        # ipdb.set_trace()

                        if emin < e0:
                            m = copy.deepcopy(membraneMC[idxMin])
                            m1 = copy.deepcopy(membraneMC[idxMin])
                            e0 = emin
                            noChange = 0
                        else:
                            noChange = noChange + 1
                            paramSA = np.exp((emin - e0) / tempSA)
                            # print(e0, e1, e2)
                            if paramSA > np.random.random():
                                m1 = copy.deepcopy(membraneMC[idxMin])
                                e1 = emin

                        if noChange > 10:
                            tempSA = tempSA * 0.5
                            piChange = piChange * 0.8
                            m1 = copy.deepcopy(m)
                            noChange = 0

                        print(step, idxMin, noChange,
                              tempSA, piChange, e0, energyMC)
                        eTotal = np.append(
                            np.delete(eTotal, 0), min(energyMC))
                        eMin = np.append(np.delete(eMin, 0), e0)
                        steps = np.append(np.delete(steps, 0), step)

                        if (step % 10) == 0:
                            filenamePlot = "R_bend_{0}_stretch_{1}_c0_{2}_nano_r_{3}_dist_{4:05.2f}.png".format(
                                var.bend_mod,  var.stretch_mod, var.c0, nano_r, nano_dist)
                            # print(filenamePlot)
                            # plotMembrane.plotMembrane(
                            #    m, var, steps, eMin, eTotal, filenamePlot)

                        if (np.sum(np.abs(eTotal - eTotal.mean())) / e0 < 1e-3) or (steps[-1] > 2000):
                            stopIteration = True

                    energyDiff = energy.energyPerMolecule(
                        m, var, energy.energy(m, var) - energy.energyPlanar(m, var))

                    result = np.append(
                        [var.bend_mod, var.stretch_mod, var.c0, nano_r,  nano_dist, steps[-1], tempSA, m.startZ], m.psi)

                    resultE = np.array([var.bend_mod, var.stretch_mod, var.c0, nano_r,  nano_dist, energyDiff, energy.energy(
                        m, var), energy.energyBending(m, var), energy.energyContact(m, var), energy.energyPlanar(m, var)])

                    filenamePickle = "ResPickle/R_bend_{0}_stretch_{1}_c0_{2}_nano_r_{3}_dist_r_{4:05.2f}.pkl".format(
                        var.bend_mod,  var.stretch_mod, var.c0, nano_r, nano_dist)

                    pickle.dump(m, open(filenamePickle, 'wb'))
                    pickle.dump(var, open(filenamePickle, 'wb'))

                    if np.size(res) == 0:
                        res = result
                        resE = resultE
                    else:
                        res = np.vstack((res, result))
                        resE = np.vstack((resE, resultE))

                    filenamePlot = "R_bend_{0}_stretch_{1}_c0_{2}_nano_r_{3}_dist_{4:05.2f}.png".format(
                        var.bend_mod,  var.stretch_mod, var.c0, nano_r, nano_dist)
                    print(filenamePlot)
                    plotMembrane.plotMembraneOnly(m, var, filenamePlot)

                filename = "R_bend_{0}_stretch_{1}_c0_{2}_nano_r_{3}.gz".format(
                    var.bend_mod,  var.stretch_mod, var.c0, nano_r)
                np.savetxt(''.join(('Res/', filename)), res, delimiter=',')
                np.savetxt(''.join(('ResE/', filename)), resE, delimiter=',')

                del res
