from atom.model import Atom
from system.model import System
from numpy import array, arange, sqrt, pi
from plotter.plotter import Plotter

a = 1.
diamond_sys = System([a / 2. * array([1., 1., 0.]),
                      a / 2. * array([0., 1., 1.]),
                      a / 2. * array([1., 0., 1.])])

diamond_sys.atoms = [Atom('C', array([0., 0., 0.])),
                     Atom('C', a / 4. * array([1., 1., 1.])),
                     ]
diamond_sys.k_points = [array([pi / a, pi / a, pi / a]),
                        array([0., 0., 0.]),
                        array([0., 2 * pi / a, 0.]),
                        array([pi / 2 / a, 2 * pi / a, pi / 2 / a]),
                        array([0., 0., 0.]),]
diamond_sys.make_k_mesh(60)
diamond_sys.parameters = {
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


for i in xrange(len(diamond_sys.atoms)):
    diamond_sys.atoms[i].orbitals = ['px', 'py', 'pz']

diamond_sys.just_do_main_magic('diamond_output')
plt = Plotter('diamond_output')
plt.plot_energy_bands_from_file()
