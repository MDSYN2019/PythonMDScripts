class Variables:
    'object containing variables'

    # bend_mod = 10        # in KT
    # hc0 = 1.47           # nm
    # stretch_mod = 0.095  # kT
    # c0 = 0
    # # parameters of nanoparticles
    # x1 = 5  # nm
    # x2 = 10  # nm
    # r1 = 1  # nm
    # r2 = 1  # nm

    def __init__(self, bend_mod=11, hc0=1.47, stretch_mod=0.095,
                 c0=0.0, x1=0.0, x2=0.0, r1=1.0, r2=1.0):
        self.bend_mod = bend_mod
        self.hc0 = hc0
        self.stretch_mod = stretch_mod
        self.c0 = c0
        self.x1 = x1
        self.x2 = x2
        self.r1 = r1
        self.r2 = r2
