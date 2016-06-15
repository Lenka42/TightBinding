from atom.model import Atom
from system.model import System
from numpy import array, sqrt, linspace, zeros, exp, pi
import matplotlib.pyplot as plt
import os
from itertools import izip

a = 1. # C-C bond length
n = 2
length = 2
system = System([], mode="with_vectors",
                name='zz_chain_{}_pz'.format(length))
system.atoms = []
for j in range(length):
    for i in range(n):
        dx = j * sqrt(3) * a
        dy = 3. * a * i
        dr = array([dx, dy, 0.])
        system.atoms += [
            Atom('C', dr + array([0., 0., 0.])),
            Atom('C', dr + array([a * sqrt(3) / 2., a / 2., 0.])),
            Atom('C', dr + array([a * sqrt(3) / 2., 3 * a / 2., 0.])),
            Atom('C', dr + array([0., 2 * a, 0.])),
        ]
system.spin_multiplier = 1
system.k_mesh = [array([0., 0., 0.]), ]
system.parameters = {
    'C': {
        'ep': 0.,
    },
    'CC': {
        'Vppp': -3.26,
    }
}


for i in xrange(len(system.atoms)):
    system.atoms[i].orbitals = ['pz']  # ,'dxy', 'dyz', 'dxz', 'dx2-y2', 'dz2']


system.just_do_main_magic()
with open(os.path.join(os.path.abspath('./outputs/'), system.name,
                       'energies'), 'r') as data:
    txt = data.read().strip().split(' ')
    energies = map(float, txt)

plt.vlines(energies, 0, 1)
plt.savefig(os.path.join(os.path.abspath('./outputs/'),
                             system.name, 'peaks.eps'))
plt.show()

en_min = - 9.29495037546
en_max = 9.29495037546
en_num = 200
en_mesh = linspace(en_min, en_max, en_num)
a = (en_max - en_min) / en_num
smoothed = zeros(en_num)
for en in energies:
    for i, point in enumerate(en_mesh):
        dos = 1 / a / sqrt(pi) * exp(- (point - en)**2 / a**2)
        smoothed[i] += dos
with open(os.path.join(os.path.abspath('./outputs/'), system.name,
                       'smoothed'), 'w') as f:
    f.write('\n'.join(' '.join(map(str, pair)) for pair in
                              izip(en_mesh, smoothed)))
plt.plot(en_mesh, smoothed)
plt.savefig(os.path.join(os.path.abspath('./outputs/'),
                         system.name, 'smoothed.eps'))
plt.show()

