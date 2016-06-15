from atom.model import Atom
from system.model import System
from numpy import array, sqrt, pi
from plotter.plotter import Plotter

#############################GRAPHENE###################################
d = 1.
t = 1.
system = System([d / 2. * array([3., sqrt(3), 0.]),
                 - d / 2. * array([3, -sqrt(3), 0.])],
                mode='with_overlap', name='grapene_sp_soc')
system.atoms = [Atom('C', array([d, 0., 0.])),
                Atom('C', array([2 * d, 0., 0.])),
                ]
system.spin_multiplier = 2
system.k_points = [array([0., 0., 0.]),
                   array([2 * pi / 3 / d, 2 * pi / 3 / sqrt(3) / d, 0]),
                   array([2 * pi / 3 / d, 0, 0]),
                   array([0., 0., 0.])]
system.make_k_mesh(100)
system.parameters = {
    'C': {
        'es': -8.7,
        'ep': 0.0,
        'lambda': 0.001,
    },
    'CC': {
        'Vsss': - 6.7,
        'Vsps': 5.5,
        'Vpps': 5.1,
        'Vppp': -3.1,
    }
}

system.s_parameters = {
    'C': {
        'es': 1.,
        'ep': 1.,
        'lambda': 0.,
    },
    'CC': {
        'Vsss': 0.2,
        'Vsps': - 0.1,
        'Vpps': - 0.15,
        'Vppp': 0.12,
    }
}


for i in xrange(len(system.atoms)):
    system.atoms[i].orbitals = ['s', 'px', 'py', 'pz', ]


system.just_do_main_magic()
plotter = Plotter(system.name)
plotter.plot_energy_bands_from_file()
