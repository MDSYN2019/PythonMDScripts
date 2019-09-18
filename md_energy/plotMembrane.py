import matplotlib.pyplot as plt
import energy
import numpy as np
import ipdb


def plotMembrane(m, var, steps, energy_min,
                 energy_tot, filename):
    plt.ion
    plt.clf()
    # plt.close()
    #fig, (a1, a2, a3) = plt.subplots(3, 1)
    #ttl = plt.suptitle(filename)
    a1 = plt.subplot(3, 1, 1)
    angle = np.arange(0.0, 180.0, 1.0)
    angle = np.radians(angle)
    xcoor = np.append(-m.get_X()[::-1], m.get_X())
    zcoor = np.append(m.get_Z()[::-1], m.get_Z())
    # ipdb.set_trace()
    # a1.fill_between(-xcoor, zcoor, -zcoor,
    #                 facecolor='blue', alpha=0.5,
    #                 )
    a1.fill_between(xcoor, zcoor, -zcoor,
                    facecolor='blue', alpha=0.5,
                    label="phospholipids")
    a1.fill_between(var.r1 * np.cos(angle) + var.x1, var.r1 * np.sin(angle),
                    -var.r1 * np.sin(angle), facecolor='red',
                    alpha=0.8)
    a1.fill_between(var.r2 * np.cos(angle) + var.x2, var.r2 * np.sin(angle),
                    -var.r2 * np.sin(angle), facecolor='red',
                    alpha=0.8, label="nanoparticle")
    a1.axis('equal')
    lgd1 = a1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    # a1.set_xlabel("$x \mathrm{[nm]}$")
    a1.set_ylabel("$y \mathrm{[nm]}$")
    a1.set_title('Shape')
    a2 = plt.subplot(3, 1, 2)
    a2.axis([steps.min(), steps.max(), energy_tot.min(), energy_tot.max()])
    a2.plot(steps, energy_min, 'r-', label='Minimum energy')
    a2.plot(steps, energy_tot, '-b>', label='Actual energy')
    a2.set_ylabel("$U [k_B T \mathrm{nm}^{-1}]$")
    lgd2 = a2.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    #a2.legend(shadow=True, fancybox=True, loc='upper left')
    # plt.plot(step, energy_bend, '-3', label='Bending energy')
    # plt.plot(step, energy_cont, '-4', label='Contact energy')
    a3 = plt.subplot(3, 1, 3)
    a3.plot(m.get_X(), energy.energyBendingDensity(
            m, var), 'r-o', label='Bending energy')
    a3.plot(m.get_X(), energy.energyContactDensity(m, var),
            'b-o', label='Contact energy')
    a3.set_xlabel("$x  \mathrm{[nm]}$")
    a3.set_ylabel("$\lambda [k_B T \mathrm{nm}^{-2}]$")
    lgd3 = a2.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.pause(0.001)


def plotMembraneOnly(m, var, filename):
    # plt.ion
  #  plt.clf()
    plt.close()
    fig, (a1, a2) = plt.subplots(2, 1, sharex=True)
    ttl = plt.suptitle(filename)
    # a1 = plt.subplot(2, 1, 1)
    angle = np.arange(0.0, 180.0, 1.0)
    angle = np.radians(angle)
    xcoor = np.append(-m.get_X()[::-1], m.get_X())
    zcoor = np.append(m.get_Z()[::-1], m.get_Z())
    # ipdb.set_trace()
    # a1.fill_between(-xcoor, zcoor, -zcoor,
    #                 facecolor='blue', alpha=0.5,
    #                 )
    a1.fill_between(xcoor, zcoor, -zcoor,
                    facecolor='blue', alpha=0.5,
                    label="phospholipids")
    a1.fill_between(var.r1 * np.cos(angle) + var.x1, var.r1 * np.sin(angle),
                    -var.r1 * np.sin(angle), facecolor='red',
                    alpha=0.8)
    a1.fill_between(var.r2 * np.cos(angle) + var.x2, var.r2 * np.sin(angle),
                    -var.r2 * np.sin(angle), facecolor='red',
                    alpha=0.8, label="nanoparticle")
    a1.axis('equal')
    lgd1 = a1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    # a1.set_xlabel("$x \mathrm{[nm]}$")
    a1.set_ylabel("$y \mathrm{[nm]}$")
    a1.set_title('Shape')
    # a2 = plt.subplot(2, 1, 2)
    a2.plot(m.get_X(), energy.energyBendingDensity(
        m, var), 'r-o', label='Bending energy')
    a2.plot(m.get_X(), energy.energyContactDensity(m, var),
            'b-o', label='Contact energy')
    a2.set_xlabel("$x  \mathrm{[nm]}$")
    a2.set_ylabel("$\lambda [k_B T \mathrm{nm}^{-2}]$")
    a2.set_title("Energy distribution")
    lgd2 = a2.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.pause(0.001)
    plt.savefig(''.join(("Fig/", filename)), bbox_extra_artists=(
        lgd1, lgd2, ttl), bbox_inches='tight')
    plt.close()
