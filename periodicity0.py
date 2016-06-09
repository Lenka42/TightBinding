from atom.model import Atom
from system.model import System
from numpy import array
import os
import matplotlib.pyplot as plt

#############################0D_Chain###################################
d = 1.
t = 1.
n = 20
system = System([], mode='standard', name='chain_0d_{}'.format(n))
system.atoms = []
for i in range(n):
    system.atoms.append(Atom('A', array([i * d, 0., 0.])))
system.k_mesh = [array([0., 0., 0.]), ]
system.parameters = {
    'A': {
        'es': 0.0,
    },
    'AA': {
        'Vsss': t,
    }
}

for i in xrange(len(system.atoms)):
    system.atoms[i].orbitals = ['s', ]

system.just_do_main_magic()
with open(os.path.join(os.path.abspath('./outputs/'), system.name,
                       'energies'), 'r') as data:
    txt = data.read().strip().split(' ')
    energies = map(float, txt)

plt.vlines(energies, 0, 1)
plt.savefig(os.path.join(os.path.abspath('./outputs/'),
                             system.name, 'peaks.eps'))
plt.show()