import scipy as sp
import numpy as np
import membrane
import parameters
import ipdb
import copy


def energyContactDensity(m, var):

    dst_plane = distance_from_plane(m)
    dst_sphere1 = distance_from_sphere(m, var.x1, var.r1)
    # dst_sphere2 = distance_from_sphere(m, var.x2, var.r2)

    # correction for acceptable value
    # nA_plane = m.get_Z() > 0
    # nA_sphere1 = dst_sphere1 > 0
    # nA_sphere2 = dst_sphere2 > 0

    # noAccept = nA_sphere1 * nA_sphere2 * nA_plane
    # findMinimum
    dst = np.min(np.array([dst_plane, dst_sphere1]), axis=0)
    noAccept = dst < var.hc0 / 4
    dst[noAccept] = 10 + 1000 * np.power(dst[noAccept], 4)

    return 0.5 * var.stretch_mod / np.power(var.hc0, 2) * np.power(dst - var.hc0, 2)


def distance_from_plane(m):
    return np.abs(m.get_Z())


def distance_from_sphere(m, xS, rS):
    return np.sqrt((m.get_X() - xS)**2 + (m.get_Z())**2) - rS


def energyBendingDensity(m, var):
    k_B = var.bend_mod
    c1 = m.get_c1()
    return 0.5 * k_B * np.power(c1 - var.c0, 2)


def energyPlanarDensity(m, var):
    k_B = var.bend_mod
    c1 = np.zeros_like(m.get_c1())
    return 0.5 * k_B * np.power(c1 - var.c0, 2)


def energyBending(m, var):
    return sp.integrate.trapz(energyBendingDensity(m, var), m.length)


def energyContact(m, var):
    return sp.integrate.trapz(energyContactDensity(m, var), m.length)


def energy(m, var):
    return sp.integrate.trapz(energyBendingDensity(m, var) +
                              energyContactDensity(m, var), m.length)


def energyPlanar(m, var):
    return sp.integrate.trapz(energyPlanarDensity(m, var))


def energyPerMolecule(m, var, energyT):
    return energyT / m.getLength() * var.lipid_area


def energyPsi(psiT, m, var):
    m1 = copy.deepcopy(m)
    m1.startZ = psiT[0]
    m1.psi = np.append(np.append([0, 0], psiT[1:]), [0, 0])
    # print(sp.integrate.trapz(energyBendingDensity(m1, var) +
    #                          energyContactDensity(m1, var), m1.length))
    return sp.integrate.trapz(energyBendingDensity(m1, var) +
                              energyContactDensity(m1, var), m1.length)


def energyLocal(psiR, m, idxR, var):
    m1 = copy.deepcopy(m)
    # psiz0Idx = sp.optimize.fmin_bfgs(funZ_zero, psiR, args=(m, idxR),
    #                                     full_output=0, disp=0)
    m1.startZ = psiR[0]
    m1.psi[idxR] = psiR[1:]

    return energy(m1, var)


def funZ_zero(psiR, m, idxR):
    psi = m.psi
    psi[idxR] = psiR
    dZ2 = - m.elementLength * np.sin(psi)
    dZ2 = np.append(0, dZ2)
    Z2 = m.startZ + dZ2.cumsum()
    return (m.endZ - Z2[-1])**2


def findLocalMin(m, rndIdx2, var):

    guess = np.append(m.startZ, m.psi[rndIdx2])
    # ipdb.set_trace()
    # print(guess)
    # x0 = sp.optimize.fmin(energyLocal, guess, args=(
    #    m, rndIdx2, var), full_output=0, disp=0)
    # x0 = sp.optimize.fmin_bfgs(energyLocal, guess, args=(
    #    m, rndIdx2, var), full_output=0, disp=0)
    # dicteq = {'type': 'eq', 'fun': funZ_zero,
    #          'args': (m, rndIdx2)}
    # x0 = sp.optimize.minimize(energyLocal, guess, args=(m, rndIdx2, var),
    #                          method='SLSQP',
    #                          constraints=dicteq)

    x0 = sp.optimize.minimize(energyLocal, guess, args=(
        m, rndIdx2, var), method='SLSQP')

    # print('fun', x0.fun)
    # print('x0', x0)
    # print(x0.message)
    # print('constraint', funZ_zero(x0.x, m, rndIdx2))

    m1 = copy.deepcopy(m)
    # psiz0Idx = sp.optimize.fmin_bfgs(funZ_zero, psiR, args=(m, idxR),
    #                                     full_output=0, disp=0)
    m1.startZ = x0.x[0]
    m1.psi[rndIdx2] = x0.x[1:]

    debug = 0
    if debug == 1:
        print(guess)
        print("DiffEnerg = ", energy(m, var) - energy(m1, var))
        # if (energy(m, var) - energy(m1, var)) < 0:
        #    ipdb.set_trace()
        print("estart = ", energy(m, var), " ",
              energyLocal(guess, m, rndIdx2, var))
        print("eend = ", energy(m1, var), " ",
              energyLocal(x0.x, m, rndIdx2, var))

    return x0.x
