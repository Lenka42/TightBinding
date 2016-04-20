from numpy import zeros, complex_, absolute, transpose
from numpy.linalg import norm, eigh
from scipy.linalg import eigvalsh
from collections import defaultdict
from itertools import izip
import os


class System(object):
    _main_methods_dict = {
        'with_vectors': 'f_with_vectors',
        'with_overlap': 'f_with_overlap',
        'standard': 'f',
    }

    def __init__(self, vectors, mode='standard', name=''):
        self.name = name
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
        self.S = None
        self.s_parameters = None
        self.num_of_bands = None
        self.just_do_main_magic = getattr(self, self._main_methods_dict[mode])
        self.output_path = os.path.join(os.path.abspath('./outputs/'),
                                        self.name)
        print self.output_path
        if not os.path.exists(self.output_path):
            print self.output_path
            os.makedirs(self.output_path)

    # TODO: cool k-d tree algorithm and second nearest neighbours
    def find_nearest_neighbours(self):
        alpha = 1.5
        atom_lst = [(i, a.r) for i, a in enumerate(self.atoms)]
        big_atom_lst = [(i, a.r) for i, a in enumerate(self.atoms)]
        # for tr_idx, translation in enumerate(self.vectors):
        #     for at_idx, atom in enumerate(self.atoms):
        #         big_atom_lst.extend([(at_idx, atom.r + translation),
        #                              (at_idx, atom.r - translation)])
        if len(self.vectors) == 2:
            print "Finding neeighbos in 2D system"
            for i in xrange(-1, 2):
                for j in xrange(-1,2):
                    if not (i == 0 and j == 0):
                        for at_idx, atom in enumerate(self.atoms):
                            big_atom_lst.append((at_idx, atom.r + i * self.vectors[0] +
                                                 j * self.vectors[1]))
        if len(self.vectors) == 1:
            print "Finding neeighbos in 1D system"
            for i in [-1, 1]:
                for at_idx, atom in enumerate(self.atoms):
                    big_atom_lst.append((at_idx, atom.r + i * self.vectors[0]))
        min_dst = 100000
        for at_idx, v in atom_lst:
            for at, vec in big_atom_lst:
                if at != at_idx:
                    min_dst = min(min_dst, norm(v - vec))
            for at, vec in big_atom_lst:
                # if at != at_idx:
                if norm(v - vec) <= alpha * min_dst and norm(v - vec) != 0.0:
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
            print len(self.k_mesh), first
            self.k_mesh += loc_k_mesh
            print len(self.k_mesh), self.k_mesh[-1]
        self.k_mesh.append(second)
        with open(os.path.join(self.output_path, 'k_points'), 'w') as f:
            f.write('\n'.join(' '.join(map(str, k)) for k in self.k_mesh))

    def f(self):
        self.find_nearest_neighbours()
        self.assign_start_indexes_to_atoms()
        with open(os.path.join(self.output_path, 'energies'), 'w') \
                as output_f:
            for k in self.k_mesh:
                self.H = zeros((self.H_matrix_dim,
                                self.H_matrix_dim),
                               dtype=complex_)
                for atom_idx, atom in enumerate(self.atoms):
                    atom.count_diagonal_matrix_elements(
                        self.parameters,
                        self.H,
                        mult=self.spin_multiplier)
                    for neighbour_atom_idx, r in self.nn_dict[atom_idx]:
                        neighbour_atom = self.atoms[neighbour_atom_idx]
                        atom.count_hamiltonian_matrix_elements(
                            neighbour_atom,
                            r, k,
                            self.parameters,
                            self.H,
                            mult=self.spin_multiplier)
                # print self.H
                energies = eigvalsh(self.H)
                output_f.write(' '.join(map(str, energies)) + '\n')

    def f_with_overlap(self):
        self.find_nearest_neighbours()
        self.assign_start_indexes_to_atoms()
        with open(os.path.join(self.output_path, 'energies'), 'w') as output_f:
            for k in self.k_mesh:
                self.H = zeros((self.H_matrix_dim,
                                self.H_matrix_dim),
                               dtype=complex_)
                self.S = zeros((self.H_matrix_dim,
                                self.H_matrix_dim),
                               dtype=complex_)
                for atom_idx, atom in enumerate(self.atoms):
                    atom.count_diagonal_matrix_elements(
                        self.parameters,
                        self.H,
                        mult=self.spin_multiplier)
                    atom.count_diagonal_matrix_elements(
                        self.s_parameters,
                        self.S,
                        mult=self.spin_multiplier)
                    for neighbour_atom_idx, r in self.nn_dict[atom_idx]:
                        neighbour_atom = self.atoms[neighbour_atom_idx]
                        atom.count_hamiltonian_matrix_elements(
                            neighbour_atom,
                            r, k,
                            self.parameters,
                            self.H,
                            mult=self.spin_multiplier)
                        atom.count_hamiltonian_matrix_elements(
                            neighbour_atom,
                            r, k,
                            self.s_parameters,
                            self.S,
                            mult=self.spin_multiplier)
                # print self.H
                energies = eigvalsh(self.H, b=self.S)
                output_f.write(' '.join(map(str, energies)) + '\n')

    def f_with_vectors(self):
        self.find_nearest_neighbours()
        self.assign_start_indexes_to_atoms()
        states_file_path = os.path.join(self.output_path, 'states')
        energies_file_path = os.path.join(self.output_path, 'energies')
        with open(energies_file_path, 'w') as output_f, \
                open(states_file_path, 'w') as output_vector_f:
            output_vector_f.write(str(self.H_matrix_dim) + '\n')
            output_vector_f.write(' '.join(
                [str(i) + at.name + orb
                 for i, at in enumerate(self.atoms)
                 for orb in at.orbitals]) + '\n')
            for k in self.k_mesh:
                self.H = zeros((self.H_matrix_dim,
                                self.H_matrix_dim),
                               dtype=complex_)
                for atom_idx, atom in enumerate(self.atoms):
                    atom.count_diagonal_matrix_elements(
                        self.parameters,
                        self.H,
                        mult=self.spin_multiplier)
                    for neighbour_atom_idx, r in self.nn_dict[atom_idx]:
                        neighbour_atom = self.atoms[neighbour_atom_idx]
                        atom.count_hamiltonian_matrix_elements(
                            neighbour_atom,
                            r, k,
                            self.parameters,
                            self.H,
                            mult=self.spin_multiplier)
                # print self.H
                energies, vectors = eigh(self.H)
                vectors = absolute(vectors)
                output_f.write(' '.join(map(str, energies)) + '\n')
                for vec in transpose(vectors):
                    vec = absolute(vec)**2
                    vec = vec / norm(vec)
                    output_vector_f.write(' '.join(map(str, vec)) + '\n')
