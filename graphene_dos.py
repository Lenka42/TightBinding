from atom.model import Atom
from system.model import System
from numpy import array, sqrt, pi, linspace
from dos.dos_calculator import DOSCalculator

#############################GRAPHENE###################################
a = 1.
system = System([a / 2. * array([1., sqrt(3), 0.]),
                 a / 2. * array([- 1., sqrt(3), 0.])],
                name="graphene_dos_pz", mode="standard")
system.atoms = [Atom('C', array([0., a / sqrt(3), 0.])),
                Atom('C', array([0., 2 * a / sqrt(3), 0.])),
                ]

system.k_points = [array([0., 0., 0.]),
                   array([pi / a, -pi / sqrt(3) / a, 0]),
                   array([4 * pi / 3 / a, 0, 0]),
                   array([0., 0., 0.])]
#system.make_k_mesh(100)
nx = 40
ny = int(2 / sqrt(3) * nx)
dx = 2 * pi / a / nx
dy = 4 * pi / sqrt(3) / a / ny
for kx in linspace(0, 2 * pi /a, nx):
    for ky in linspace(1 / sqrt(3) * kx, 4 * pi / sqrt(3) / a + 1 / sqrt(3) * kx, ny):
        system.k_mesh.append(array([kx, ky, 0.]))
system.parameters = {
    'C': {
        'ep': 0.,
    },
    'CC': {
        'Vppp': -3.26,
    }
}


for i in xrange(len(system.atoms)):
    system.atoms[i].orbitals = ['pz']


system.just_do_main_magic()
doser = DOSCalculator(system.dim, system.name, 200)
doser.f()

