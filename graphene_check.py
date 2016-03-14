from atom.model import Atom
from system.model import System
from numpy import array, sqrt, pi
from plotter.plotter import Plotter

#############################GRAPHENE###################################
a = 1.
system = System([a / 2. * array([1., sqrt(3), 0.]),
                 a / 2. * array([- 1., sqrt(3), 0.])], mode="with_overlap")
system.name = 'graphene_sp_S'
system.atoms = [Atom('C', array([0., a / sqrt(3), 0.])),
                Atom('C', array([0., 2 * a / sqrt(3), 0.])),
                ]

system.num_of_bands = 4
system.k_points = [array([0., 0., 0.]),
                   array([pi / a, -pi / sqrt(3) / a, 0]),
                   array([4 * pi / 3 / a, 0, 0]),
                   array([0., 0., 0.])]
system.make_k_mesh(100)
system.parameters = {
    'C': {
        'es': 8.370,
        'ep': 0.0,
    },
    'CC': {
        'Vsss': -5.729,
        'Vsps': 5.618,
        'Vppp': 6.050,
        'Vpps': -3.070
    }
}
system.s_parameters = {
    'C': {
        'es': 1.,
        'ep': 1.,
        'ed': 1.,
        'lambda': 0,
    },
    'CC': {
        'Vsss': 0.102,
        'Vsps': - 0.171,
        'Vppp': -0.377,
        'Vpps': 0.070
    }
}


for i in xrange(len(system.atoms)):
    system.atoms[i].orbitals = ['s', 'px', 'py', 'pz']


system.just_do_main_magic()
plt = Plotter(system.name)
plt.plot_energy_bands_from_file()
