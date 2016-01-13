from atom.model import Atom
from system.model import System
from numpy import array, sqrt, pi
import time
import matplotlib.pyplot as plt

#############################GRAPHENE###################################
a = 1.
ns = range(1,51,2)
ns_sq = [2 * i ** 2 for i in ns]
ts = []
for N in ns:
    tp = time.time()
    system = System([N * a / 2. * array([1., sqrt(3), 0.]),
                     N * a / 2. * array([- 1., sqrt(3), 0.])])
    system.name = 'graphene_pz_time_eigvalsh' + str(N)
    atoms_lst = []
    for i in range(N):
        for j in range(N):
            atoms_lst.extend([Atom('C', i * system.vectors[0] / N +
                                   j * system.vectors[1] / N +
                                   array([0., a / sqrt(3), 0.])),
                              Atom('C', i * system.vectors[0] / N +
                                   j * system.vectors[1] / N +
                                   array([0., 2 * a / sqrt(3), 0.]))])
    system.atoms = atoms_lst
    system.k_points = [array([0., 0., 0.]),
                       array([pi / a, -pi / sqrt(3) / a, 0]),
                       array([4 * pi / 3 / a, 0, 0]),
                       array([0., 0., 0.])]
    system.make_k_mesh(32)
    system.parameters = {
        'C': {
            'es': 0,
            'ep': 7.4,
        },
        'CC': {
            'Vsss': -3.8,
            'Vsps': 4.44,
            'Vppp': -1.325,
            'Vpps': 4.9
        }
    }

    for i in xrange(len(system.atoms)):
        system.atoms[i].orbitals = ['pz']
    system.just_do_main_magic()
    tk = time.time()
    t = (tk - tp) / 39
    ts.append(t)
print ns_sq
print ts
plt.plot(ns_sq, ts)
plt.show()

