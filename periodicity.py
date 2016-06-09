from atom.model import Atom
from system.model import System
from numpy import array, pi
from plotter.plotter import Plotter
from dos.dos_calculator import DOSCalculator

#############################1D_Chain###################################
d = 1.
t = 1.
system = System([array([d, 0, 0])], mode='standard', name='chain2_1d')
system.atoms = [Atom('A', array([d/2, 0., 0.])),
                Atom('A', array([0., 0., 0.])),]
system.k_points = [array([- pi / d, 0., 0.]),
                   array([0., 0., 0.]),
                   array([pi / d, 0., 0.]),
                   ]
system.make_k_mesh(100)
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
plt = Plotter(system.name)
plt.new_plot_energy_bands_from_file()
doser = DOSCalculator(system.dim, system.name, 70)
doser.f()