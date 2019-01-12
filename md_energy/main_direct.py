# general import
import numpy as np
import scipy as sp
import scipy.constants as const
import ipdb
import copy
import pickle
import multiprocessing

# own classes
import membrane
import energy
import parameters
import plotMembrane
import monteCarlo2
import variables


var = variables.Variables()
membraneLength = 30.0  # membrane length
elementNum = 100  # numbers of elements
elementLength = membraneLength / elementNum

elementNum = np.round(membraneLength / elementLength)

nano_r = 2
nano_dist = 10

# DOPC
var.bend_mod = 9  # kbT
var.lipid_area = 72.4 * 1e-2  # angstrom2
var.stretch_mod = 23  # kbTnm-2
var.c0 = 0  # nm-1


for nano_r in np.array([1.0, 2.0, 3.0, 4.0]):
    for var.c0 in np.array([-1 / 16, 0, 1 / 16]):
        for var.stretch_mod in np.array([23]):
            for var.bend_mod in np.array([9.0]):

                # for nano_r in np.array([1.0, 1.5, 2.0, 2, 5, 3.0]):
                #     for var.c0 in np.array([-1 / 16, 0, 1 / 16]):
                #         for var.stretch_mod in np.array([2.3, 23, 230]):
                #             for var.bend_mod in np.array([5.0, 9.0, 15.0]):
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

                # np.linspace(0, 40, 41):
                for nano_dist in np.linspace(0, 20, 21):
                    var.x1 = nano_dist / 2 + var.r1
                    var.x2 = -nano_dist / 2 - var.r1

                    m = membrane.Membrane(elementPsi, elementLength,
                                          startZ, endZ)
                    e0 = energy.energy(m, var)
                    # boundaries for psi values
                    boundPsi = [(var.hc0, var.hc0 + nano_r)]
                    for eachPsi in m.psi[2:-2]:
                        boundPsi.append((-np.pi / 2, np.pi / 2))

                    x01 = sp.optimize.minimize(energy.energyPsi, np.append(m.startZ, m.psi[2:-2]), args=(
                        m, var), method='SLSQP', bounds=boundPsi, options={'disp': True, 'maxiter': 1000})

                    # x02 = sp.optimize.differential_evolution(energy.energy,
                    # m.psi, args=(m, var))

                    m.startZ = x01.x[0]
                    m.psi = np.append(np.append([0, 0], x01.x[1:]), [0, 0])
                    # print(x01)
                    filenamePlot = "R_bend_{0}_stretch_{1}_c0_{2}_nano_r_{3}_dist_{4:05.2f}.png".format(
                        var.bend_mod,  var.stretch_mod, var.c0, nano_r, nano_dist)
                    print(filenamePlot)

                    plotMembrane.plotMembraneOnly(m, var, filenamePlot)

                    energyDiff = energy.energyPerMolecule(m, var,
                                                          energy.energy(m, var) - energy.energyPlanar(m, var))

                    energy

                    result = np.append(
                        [var.bend_mod, var.stretch_mod, var.c0, nano_r,  nano_dist, x01.success, x01.fun, m.startZ], m.psi)

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

                filename = "R_bend_{0}_stretch_{1}_c0_{2}_nano_r_{3}.gz".format(
                    var.bend_mod,  var.stretch_mod, var.c0, nano_r)
                np.savetxt(''.join(('Res/', filename)), res, delimiter=',')
                np.savetxt(''.join(('ResE/', filename)), resE, delimiter=',')

                del res
