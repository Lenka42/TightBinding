from atom.model import Atom
from system.model import System
from numpy import array, sqrt, pi
from numpy.linalg import norm
from plotter.plotter import Plotter

##########################GRAPHENE_WITH SOC###################################
a = 1. # unit vector length
system = System([a / 2. * array([1., sqrt(3), 0.]),
                 a / 2. * array([- 1., sqrt(3), 0.])],
                mode="standard", name='graphene_pd_soc_gap')
system.atoms = [Atom('C', array([0., a / sqrt(3), 0.])),
                Atom('C', array([0., 2 * a / sqrt(3), 0.])),
                ]
system.spin_multiplier = 2
distanse = 0.1
k_point = array([4 * pi / 3 / a, 0, 0])
l_point = array([pi / a, -pi / sqrt(3) / a, 0])
g_point = array([0, 0, 0])
l_side = k_point + (l_point - k_point) / norm(l_point - k_point) * distanse
g_side = k_point + (g_point - k_point) / norm(g_point - k_point) * distanse
system.k_points = [l_side, k_point, g_side]
system.make_k_mesh(200)
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
    system.atoms[i].orbitals = ['pz', 'dxy', 'dyz', 'dxz', 'dx2-y2', 'dz2']


system.just_do_main_magic()
plt = Plotter(system.name)
plt.new_plot_energy_bands_from_file()
