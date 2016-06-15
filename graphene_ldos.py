from atom.model import Atom
from system.model import System
from numpy import array, sqrt, pi, linspace
from plotter.plotter import Plotter
from dos.dos_calculator import DOSCalculator, LDOSCalculator

#############################GRAPHENE###################################
a = 1.
system = System([a / 2. * array([1., sqrt(3), 0.]),
                 a / 2. * array([- 1., sqrt(3), 0.])],
                name="graphene_ldos_soc", mode="with_vectors")
system.atoms = [Atom('C', array([0., a / sqrt(3), 0.])),
                Atom('C', array([0., 2 * a / sqrt(3), 0.])),
                ]

system.k_points = [array([0., 0., 0.]),
                   array([pi / a, -pi / sqrt(3) / a, 0]),
                   array([4 * pi / 3 / a, 0, 0]),
                   array([0., 0., 0.])]
#system.make_k_mesh(100)
nx = 100
ny = int(2 / sqrt(3) * nx)
dx = 2 * pi / a / nx
dy = 4 * pi / sqrt(3) / a / ny
for kx in linspace(0, 2 * pi /a, nx):
    for ky in linspace(1 / sqrt(3) * kx, 4 * pi / sqrt(3) / a + 1 / sqrt(3) * kx, ny):
        system.k_mesh.append(array([kx, ky, 0.]))
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
lst = system.find_indeces_for_ldos(orbital='d')
print lst
plt = Plotter(system.name)
plt.plot_energy_bands_from_file()

doser = LDOSCalculator(system.dim, system.name, 500, 0, indeces_list=lst)
doser.f()

