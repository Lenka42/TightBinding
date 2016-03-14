import matplotlib.pyplot as plt
import os


class Plotter(object):
    def __init__(self, input_file):
        self.input_file = input_file

    def plot_energy_bands_from_file(self):
        k_mesh = []
        energies = []
        first_line = True  # lena krivorukaja
        n = 0
        with open(os.path.join(os.path.abspath('./outputs/'), self.input_file), 'r')\
                as data:
            txt = data.read().strip().split('\n')
            for line in txt:
                n += 1
                lst = line.split(' ')
                k_mesh.append([lst[0], lst[1], lst[2]])
                if first_line:
                    for i, en in enumerate(lst[3:]):
                        energies.append([en, ])
                    first_line = False
                else:
                    for i, en in enumerate(lst[3:]):
                        energies[i].append(en)

        x = range(n)

        for band in energies:
            plt.plot(x, band)
            plt.savefig(os.path.join(os.path.abspath('./output_images/'),
                                     self.input_file + '.png'))
        plt.show()

