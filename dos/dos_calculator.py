import matplotlib.pyplot as plt
import numpy as np
import os
from scipy import sqrt, pi, exp
from scipy.integrate import quad, nquad
from scipy.interpolate import interp1d, LinearNDInterpolator
from itertools import izip
sqrt_pi = sqrt(pi)


class DOSCalculator(object):
    _dos_methods_dict = {
        1: 'f1d',
        2: 'f2d',
    }

    def __init__(self, dim, name, en_num):
        self.dim = dim
        self.name = name
        self.en_num = en_num
        self.f = getattr(self, self._dos_methods_dict[dim])
        self.a = 0.01
        self.dos = np.zeros(self.en_num)
        self.output_path = os.path.join(os.path.abspath('./outputs/'),
                                        self.name)

    def f1d(self):
        with open(os.path.join(os.path.abspath('./outputs/'),
                               self.name, 'energies')) as f:
            lines = f.readlines()
            energies = np.array([map(float, line.strip().split())
                                 for line in lines])
        energies = np.transpose(energies)
        min_en = np.amin(energies)
        max_en = np.amax(energies)
        self.en_mesh = np.linspace(min_en, max_en, self.en_num)
        self.a = (max_en - min_en) / self.en_num
        with open(os.path.join(os.path.abspath('./outputs/'),
                               self.name, 'k_points')) as f:
            lines = f.readlines()
            self.k_points = np.array([float(line.strip().split()[0])
                                     for line in lines])
        for band in energies:
            en_on_k = interp1d(self.k_points, band)
            for i, en in enumerate(self.en_mesh):
                dos = quad(lambda k: 1 / self.a / sqrt_pi *
                           exp(- (en_on_k(k) - en)**2 / self.a**2),
                           self.k_points[0], self.k_points[-1])
                self.dos[i] += dos[0]
        with open(os.path.join(self.output_path, 'dos'), 'w') as f:
            f.write('\n'.join(' '.join(map(str, pair)) for pair in
                              izip(self.en_mesh, self.dos)))

        plt.plot(self.en_mesh, self.dos)
        plt.savefig(os.path.join(self.output_path, 'dos.eps'))
        plt.show()

    def f2d(self):
        with open(os.path.join(os.path.abspath('./outputs/'),
                               self.name, 'energies')) as f:
            lines = f.readlines()
            energies = np.array([map(float, line.strip().split())
                                 for line in lines])
        energies = np.transpose(energies)
        min_en = np.amin(energies)
        max_en = np.amax(energies)
        self.en_mesh = np.linspace(min_en, max_en, self.en_num)
        self.a = (max_en - min_en) / self.en_num

        def bounds_kx():
            return [0, 2 * pi]

        def bounds_ky(kx):
            return [kx / sqrt(3), kx / sqrt(3) + 4 * pi / sqrt(3)]

        with open(os.path.join(os.path.abspath('./outputs/'),
                               self.name, 'k_points')) as f:
            lines = f.readlines()
            self.k_points = np.array([map(float, line.strip().split()[:2])
                                      for line in lines])
        for band in energies:
            en_on_k = LinearNDInterpolator(self.k_points, band)
            for i, en in enumerate(self.en_mesh):
                def f(ky, kx):
                    return 1 / self.a / sqrt_pi * \
                           exp(- (en_on_k(kx, ky) - en) ** 2 / self.a ** 2)
                dos = nquad(f, [bounds_ky, bounds_kx])
                self.dos[i] += dos[0]
        with open(os.path.join(self.output_path, 'dos'), 'w') as f:
            f.write('\n'.join(' '.join(map(str, pair)) for pair in
                              izip(self.en_mesh, self.dos)))

        plt.plot(self.en_mesh, self.dos)
        plt.savefig(os.path.join(self.output_path, 'dos.eps'))
        plt.show()
