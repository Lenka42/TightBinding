from numpy import zeros, complex_
from numpy.linalg import norm, eigvalsh
from collections import defaultdict
from itertools import izip
import os


class System(object):
    def __init__(self, vectors):
        self.name = ''
        self.spin_multiplier = 1
        self.vectors = vectors
        self.atoms = None
        self.dim = len(self.vectors)
        self.k_mesh = []
        self.nn_dict = defaultdict(list)
        self.parameters = None
        self.H_matrix_dim = None
        self.H = None
        self.k_points = None

    #TODO: cool k-d tree algorithm and second nearest neighbours
    def find_nearest_neighbours(self):
        alpha = 1.2
        atom_lst = [(i, a.r) for i, a in enumerate(self.atoms)]
        big_atom_lst = [(i, a.r) for i, a in enumerate(self.atoms)]
        for tr_idx, translation in enumerate(self.vectors):
            for at_idx, atom in enumerate(self.atoms):
                big_atom_lst.extend([(at_idx, atom.r + translation), (at_idx, atom.r - translation)])

        min_dst = 100000
        for at_idx, v in atom_lst:
            for at, vec in big_atom_lst:
                if at != at_idx:
                    min_dst = min(min_dst, norm(v - vec))
            for at, vec in big_atom_lst:
                if at != at_idx:
                    if norm(v - vec) <= alpha * min_dst:
                        self.nn_dict[at_idx].append((at, vec - v))
        print self.nn_dict

    def assign_start_indexes_to_atoms(self):
        start_idx = 0
        for atom in self.atoms:
            atom.start_idx = start_idx
            start_idx += len(atom.orbitals) * self.spin_multiplier
        self.H_matrix_dim = start_idx

    def make_k_mesh(self, n):
        k_distance = 0.
        for point in self.k_points:
            try:
                d = norm(point - previous_point)
                k_distance += d
            except NameError:
                previous_point = point
        approximate_delta = k_distance / n
        for first, second in izip(self.k_points, self.k_points[1:]):
            n_loc = int(norm(second - first) / approximate_delta)
            delta_k = (second - first) / n_loc
            loc_k_mesh = [first + i * delta_k for i in range(0, n_loc)]
            self.k_mesh += loc_k_mesh
        print len(self.k_mesh)

    def just_do_main_magic(self):
        self.find_nearest_neighbours()
        self.assign_start_indexes_to_atoms()
        with open(os.path.join(os.path.abspath('./outputs/'), self.name), 'w')\
                as output_f:
            for k in self.k_mesh:
                self.H = zeros((self.H_matrix_dim,
                                self.H_matrix_dim),
                               dtype=complex_)
                for atom_idx, atom in enumerate(self.atoms):
                    atom.count_diagonal_matrix_elements(self.parameters,
                                                        self.H,
                                                        mult=self.spin_multiplier)
                    for neighbour_atom_idx, r in self.nn_dict[atom_idx]:
                        neighbour_atom = self.atoms[neighbour_atom_idx]
                        atom.count_hamiltonian_matrix_elements(neighbour_atom,
                                                               r, k,
                                                               self.parameters,
                                                               self.H,
                                                               mult=self.spin_multiplier)
                #print self.H
                energies = eigvalsh(self.H)
                output_f.write(' '.join(map(str, k) + map(str, energies)) +
                               '\n')
