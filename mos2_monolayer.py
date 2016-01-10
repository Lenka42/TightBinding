from atom.model import Atom
from system.model import System
from numpy import array, sqrt, pi
from plotter.plotter import Plotter

#######################MoS2_MONOLAYER_WITH SOC#################################
a = 3.12
c = 3.11
system = System([a / 2. * array([1., sqrt(3), 0.]),
                 a / 2. * array([- 1., sqrt(3), 0.])])
system.name = 'mos2_mono_soc_2'
system.atoms = [Atom('Mo', array([0., a / sqrt(3), 0.])),
                Atom('S', array([0., 2 * a / sqrt(3), c / 2.])),
                Atom('S', array([0., 2 * a / sqrt(3), - c / 2.])),
                ]
system.spin_multiplier = 2
system.k_points = [array([0., 0., 0.]),
                   array([pi / a, -pi / sqrt(3) / a, 0]),
                   array([4 * pi / 3 / a, 0, 0]),
                   array([0., 0., 0.])]
system.make_k_mesh(100)
system.parameters = {
    'S': {
        'es': 7.6595,
        'ep': -2.1537,
        'ed': 8.7689,
        'lambda': 0.2129,
    },
    'Mo': {
        'es': 5.5994,
        'ep': 6.7128,
        'ed': 2.6429,
        'lambda': 1.0675,
    },
    'SMo': {
        'Vsss': -0.0917,
        'Vsps': 0.6656,
        'Vpps': 1.4008,
        'Vppp': -0.4812,
        'Vsds': 0.2177,
        'Vpds': -2.8732,
        'Vpdp': 0.7739,
        'Vdds': -3.1425,
        'Vddp': 2.4975,
        'Vddd': -0.3703,
    },
    'MoS': {
        'Vsss': -0.0917,
        'Vsps': -1.6515,
        'Vpps': 1.4008,
        'Vppp': -0.4812,
        'Vsds': -1.0654,
        'Vpds': 2.1898,
        'Vpdp': -1.9408,
        'Vdds': -3.1425,
        'Vddp': 2.4975,
        'Vddd': -0.3703,
    },
    'SS': {
        'Vsss': 0.3093,
        'Vsps': - 0.9210,
        'Vpps': 0.7132,
        'Vppp': - 0.1920,
        'Vsds': - 0.2016,
        'Vpds': - 0.5204,
        'Vpdp': - 0.1203,
        'Vdds': 0.8347,
        'Vddp': 0.7434,
        'Vddd': -0.1919,
    },
    'MoMo': {
        'Vsss': 0.1768,
        'Vsps': 1.0910,
        'Vpps': - 0.3842,
        'Vppp': 0.5203,
        'Vsds': - 0.5635,
        'Vpds': -0.2316,
        'Vpdp': 0.0582,
        'Vdds': 0.3602,
        'Vddp': 0.0432,
        'Vddd': 0.1108,
    }

}


for i in xrange(len(system.atoms)):
    system.atoms[i].orbitals = ['s', 'px', 'py', 'pz', 'dxy', 'dyz', 'dxz',
                                'dx2-y2', 'dz2']


system.just_do_main_magic()
plt = Plotter(system.name)
plt.plot_energy_bands_from_file()
