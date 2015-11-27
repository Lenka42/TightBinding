from atom.model import Atom
from system.model import System
from numpy import array, pi
from plotter.plotter import Plotter

###########################DIAMOND##############################
a = 1.
N = 2
system = System([N * a / 2. * array([1., 1., 0.]),
                 N * a / 2. * array([0., 1., 1.]),
                 N * a / 2. * array([1., 0., 1.])])
system.name = 'diamond_multi'
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
print [x.r for x in list_of_atoms]
system.k_points = [array([pi / a, pi / a, pi / a]),
                   array([0., 0., 0.]),
                   array([0., 2 * pi / a, 0.]),
                   array([pi / 2 / a, 2 * pi / a, pi / 2 / a]),
                   array([0., 0., 0.]), ]
system.make_k_mesh(100)
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
plt.plot_energy_bands_from_file()
