from cmath import exp
from numpy import dot
from numpy.linalg import norm
from hopping_integrals import HOPPING_INTEGRALS


class Atom(object):

    default_orbitals = ['s',
                        'px',
                        'py',
                        'pz',
                        'dxy',
                        'dyz',
                        'dxz',
                        'dx2-y2',
                        'dz2',
                        ]

    def __init__(self, name, r):
        self.name = name
        self.orbitals = None
        self.start_idx = None
        self.r = r

    def count_diagonal_matrix_elements(self, parameters, H, mult=1):
        for idx, orb in enumerate(self.orbitals):
            H[self.start_idx + idx * mult, self.start_idx + idx * mult] = \
                parameters[self.name]['e' + orb[0]]
            if mult == 2:
                H[self.start_idx + idx * mult + 1,
                  self.start_idx + idx * mult + 1] = \
                    parameters[self.name]['e' + orb[0]]

    def count_hamiltonian_matrix_elements(self, other, r, k, parameters, H, mult=1):
        params = parameters[self.name + other.name]
        params.update({'nx': r[0] / norm(r),
                       'ny': r[1] / norm(r),
                       'nz': r[2] / norm(r)})
        phase = exp(1j * dot(k, r))
        for idx1, orbital1 in enumerate(self.orbitals):
            num1 = Atom.default_orbitals.index(orbital1)
            for idx2, orbital2 in enumerate(other.orbitals):
                num2 = Atom.default_orbitals.index(orbital2)
                H[self.start_idx + idx1 * mult, other.start_idx + idx2 * mult]\
                    += HOPPING_INTEGRALS[num1][num2].subs(params) * phase
