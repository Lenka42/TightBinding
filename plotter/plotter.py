import matplotlib.pyplot as plt
import os


class Plotter(object):
    def __init__(self, name):
        self.name = name

    def plot_energy_bands_from_file(self):
        energies = []
        first_line = True  # lena krivorukaja
        n = 0
        with open(os.path.join(os.path.abspath('./outputs/'), self.name,
                               'energies'), 'r') as data:
            txt = data.read().strip().split('\n')
            for line in txt:
                n += 1
                lst = line.split(' ')
                if first_line:
                    for i, en in enumerate(lst):
                        energies.append([en, ])
                    first_line = False
                else:
                    for i, en in enumerate(lst):
                        energies[i].append(en)

        x = range(n)

        for band in energies:
            plt.plot(x, band)
            plt.savefig(os.path.join(os.path.abspath('./outputs/'),
                                     self.name, 'band_structure.eps'))
        plt.show()
