from atom.model import Atom
from system.model import System
from numpy import array, sqrt, pi
from plotter.plotter import Plotter
from dos.dos_calculator import DOSCalculator, LDOSCalculator
from copy import deepcopy

######################ZigZag_GRAPHENE_ribbon_WITH DOS##########################
a = 1. # C-C bond length
n = 5
system = System([array([a * sqrt(3), 0., 0.])], mode="with_vectors",
                name='zz_ribbon_{}_pz_ldos'.format(n))
system.atoms = [Atom('C', array([0., 0., 0.])),
                Atom('C', array([a * sqrt(3) / 2., a / 2., 0.])),
                Atom('C', array([a * sqrt(3) / 2., 3 * a / 2., 0.])),
                Atom('C', array([0., 2 * a, 0.])),
                ]
four_atoms_cell = deepcopy(system.atoms)
shift_r = array([0., 3. * a, 0.])
for i in range(1, n):
    for atom in four_atoms_cell:
        new_atom = deepcopy(atom)
        new_atom.r = new_atom.r + i * shift_r
        system.atoms.append(new_atom)
system.spin_multiplier = 1
system.k_points = [array([- pi / sqrt(3) / a, 0., 0.]),
                   array([0., 0., 0.]),
                   array([pi / sqrt(3) / a, 0, 0]),
                   ]
system.make_k_mesh(150)
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
    system.atoms[i].orbitals = ['pz']  # ,'dxy', 'dyz', 'dxz', 'dx2-y2', 'dz2']


system.just_do_main_magic()
lst = system.find_indeces_for_ldos(atom_idx=0)
print lst
plt = Plotter(system.name)
plt.plot_energy_bands_from_file()

doser = LDOSCalculator(1, system.name, 200, 0, indeces_list=lst)
doser.f()
