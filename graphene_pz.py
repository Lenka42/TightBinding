from atom.model import Atom
from system.model import System
from numpy import array, sqrt, pi, sin, cos
from plotter.plotter import Plotter
import os
import matplotlib.pyplot as plt

#############################GRAPHENE###################################
d = 1.
t = 1.
system = System([d / 2. * array([3., sqrt(3), 0.]),
                 - d / 2. * array([3, -sqrt(3), 0.])],
                mode='standard', name='graphene_pz')
system.atoms = [Atom('C', array([d, 0., 0.])),
                Atom('C', array([2 * d, 0., 0.])),
                ]

system.k_points = [array([0., 0., 0.]),
                   array([2 * pi / 3 / d, 2 * pi / 3 / sqrt(3) / d, 0]),
                   array([2 * pi / 3 / d, 0, 0]),
                   array([0., 0., 0.])]
system.make_k_mesh(100)
system.parameters = {
    'C': {
        'ep': 0.0,
    },
    'CC': {
        'Vppp': t,
    }
}


for i in xrange(len(system.atoms)):
    system.atoms[i].orbitals = ['pz', ]


system.just_do_main_magic()
plotter = Plotter(system.name)
plotter.new_plot_energy_bands_from_file()

# plot with analytical functions
with open(os.path.join(os.path.abspath('./outputs/'), system.name,
                       'analytical_energies'), 'w') as output_f:
    for k in system.k_mesh:
        energy = t * sqrt(3. + 2. * cos(sqrt(3.) * k[1] * d) +
                          4. * cos(sqrt(3.) / 2. * k[1] * d) *
                          cos(3. / 2. * k[0] * d))
        energies = [- energy, energy]
        output_f.write(' '.join(map(str, energies)) + '\n')

energies = []
first_line = True  # lena krivorukaja
n = 0
with open(os.path.join(os.path.abspath('./outputs/'), system.name,
                       'energies'), 'r') as data:
    txt = data.read().strip().split('\n')
    for line in txt:
        n += 1
        lst = line.split(' ')
        if first_line:
            for i, en in enumerate(lst):
                energies.append([en, ])
            first_line = False
        else:
            for i, en in enumerate(lst):
                energies[i].append(en)

first_line = True  # lena krivorukaja
ln = len(energies)
with open(os.path.join(os.path.abspath('./outputs/'), system.name,
                       'analytical_energies'), 'r') as data:
    txt = data.read().strip().split('\n')
    for line in txt:
        lst = line.split(' ')
        if first_line:
            for i, en in enumerate(lst):
                energies.append([en, ])
            first_line = False
        else:
            for i, en in enumerate(lst):
                energies[ln + i].append(en)

x = range(n)
for band in energies:
    print len(x), len(band)
    plt.plot(x, band)
    plt.savefig(os.path.join(os.path.abspath('./outputs/'),
                             system.name, 'band_structure.eps'))
plt.show()
