from atom.model import Atom
from system.model import System
from numpy import array, pi
from plotter.plotter import Plotter
import time
import matplotlib.pyplot as plt

###########################DIAMOND##############################
a = 1.
ns = range(1, 15, 3)
ncbs = [2 * i ** 3 for i in ns]
ts = []
for N in ns:
    tcurr = []
    for i in range(1):
        tp = time.time()
        system = System([N * a / 2. * array([1., 1., 0.]),
                         N * a / 2. * array([0., 1., 1.]),
                         N * a / 2. * array([1., 0., 1.])])
        system.name = 'diamond_multi' + str(N)
        list_of_atoms = []
        for i in xrange(N):
            for j in xrange(N):
                for k in xrange(N):
                    list_of_atoms.extend([Atom('C', i * system.vectors[0] / N +
                                               j * system.vectors[1] / N +
                                               k * system.vectors[2] / N +
                                               array([0., 0., 0.])),
                                          Atom('C', i * system.vectors[0] / N +
                                               j * system.vectors[1] / N +
                                               k * system.vectors[2] / N +
                                               a / 4. * array([1., 1., 1.]))])
        system.atoms = list_of_atoms
        system.k_points = [array([pi / a, pi / a, pi / a]),
                           array([0., 0., 0.]),
                           array([0., 2 * pi / a, 0.]),
                           array([pi / 2 / a, 2 * pi / a, pi / 2 / a]),
                           array([0., 0., 0.]), ]
        system.make_k_mesh(20)
        system.parameters = {
            'C': {
                'es': 0,
                'ep': 7.4,
            },
            'CC': {
                'Vsss': -3.8,
                'Vsps': 4.44,#10.25 * sqrt(3) / 4.,
                'Vppp': -1.325,
                'Vpps': 4.9
            }
        }

        for i in xrange(len(system.atoms)):
            system.atoms[i].orbitals = ['s', 'px', 'py', 'pz']

        system.just_do_main_magic()
        plt = Plotter(system.name)
#        plt.plot_energy_bands_from_file()
        tk = time.time()
        t = tk - tp
        print 'T ' + str(N) + ' = ' + str(t)
        tcurr.append(t)
    t = sum(tcurr) / len(tcurr)
    ts.append(t)
print ts
plt.plot(ncbs, ts)
plt.show()
