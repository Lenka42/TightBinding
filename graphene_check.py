from atom.model import Atom
from system.model import System
from numpy import array, sqrt, pi
from plotter.plotter import Plotter

#############################GRAPHENE###################################
a = 1.
system = System([a / 2. * array([1., sqrt(3), 0.]),
                 a / 2. * array([- 1., sqrt(3), 0.])])
system.name = 'graphene_sp'
system.atoms = [Atom('C', array([0., a / sqrt(3), 0.])),
                Atom('C', array([0., 2 * a / sqrt(3), 0.])),
                ]
system.k_points = [array([0., 0., 0.]),
                   array([pi / a, -pi / sqrt(3) / a, 0]),
                   array([4 * pi / 3 / a, 0, 0]),
                   array([0., 0., 0.])]
system.make_k_mesh(100)
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
    system.atoms[i].orbitals = ['s', 'px', 'py', 'pz']


system.just_do_main_magic()
plt = Plotter(system.name)
plt.plot_energy_bands_from_file()
