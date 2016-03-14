from atom.model import Atom
from system.model import System
from numpy import array, sqrt, pi
from plotter.plotter import Plotter

##########################GRAPHENE_WITH SOC###################################
a = 1.
system = System([a / 2. * array([1., sqrt(3), 0.]),
                 a / 2. * array([- 1., sqrt(3), 0.])], mode="standard")
system.name = 'graphene_pd_soc3'
system.atoms = [Atom('C', array([0., a / sqrt(3), 0.])),
                Atom('C', array([0., 2 * a / sqrt(3), 0.])),
                ]
system.spin_multiplier = 2
system.k_points = [array([0., 0., 0.]),
                   array([pi / a, -pi / sqrt(3) / a, 0]),
                   array([4 * pi / 3 / a, 0, 0]),
                   array([0., 0., 0.])]
system.make_k_mesh(100)
system.parameters = {
    'C': {
        'ep': 1.2057,
        'ed': 24.1657,
        'lambda': 0.001
    },
    'CC': {
        'Vppp': -3.26,
        'Vpps': 0.0,
        'Vpds': 0.0,
        'Vpdp': 2.4,
        'Vdds': 0.0,
        'Vddp': 3.6,
        'Vddd': -7.4
    }
}


for i in xrange(len(system.atoms)):
    system.atoms[i].orbitals = ['pz', 'dxy', 'dyz', 'dxz']#, 'dx2-y2', 'dz2']


system.just_do_main_magic()
plt = Plotter(system.name)
plt.plot_energy_bands_from_file()
