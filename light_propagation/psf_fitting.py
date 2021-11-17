import sys
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt

class Fitter:
    def __init__(self, dirname):
        self.work_dir = f"./results/{dirname}/"

        self.inputx = 61
        self.inputy = 61
        self.total_depth = 28
        self.gate = 16

        with open(self.work_dir + "tssp_fit.csv", 'w') as f:
            f.write("gate,z,centroid_x,centroid_y,sigma_x,sigma_y,amplitude\n")

    def twoD_gaussian(self, coord, cen_x, cen_y, sig_x, sig_y, amp):
        x, y = coord
        function = amp * np.exp(-0.5 * (((x - cen_x) / sig_x) ** 2 +
                                ((y - cen_y) / sig_y) ** 2))

        return function.ravel()

    def find_nearest_gaussian_param(self):
        for gatenum in range(self.gate):
            tssp_map = self.work_dir + f"tssp_map/tssp(gate={gatenum}).csv"
            tssp = np.loadtxt(tssp_map, skiprows=2, delimiter=',')
            tssp = np.reshape(tssp,
                              (self.total_depth, self.inputx, self.inputy))

            for z in range(self.total_depth):
                xval, yval = np.meshgrid(
                                np.linspace(0, self.inputx - 1, self.inputx),
                                np.linspace(0, self.inputy - 1, self.inputy))

                # gaussian function parameters: cen_x, cen_y, sig_x, sig_y, amp
                initial_guess = (30, 30, 10, 10, 0.01)
                try:
                    popt, _ = opt.curve_fit(self.twoD_gaussian, (xval, yval),
                                            tssp[z].ravel(), p0=initial_guess)
                except:
                    print("No optical param")
                    continue

                with open(self.work_dir + "tssp_fit.csv", 'a') as f:
                    f.write(f"{gatenum},{z},{popt[0]},{popt[1]},")
                    f.write(f"{popt[2]},{popt[3]},{popt[4]}\n")


def main():
    fitter = Fitter(sys.argv[1])
    fitter.find_nearest_gaussian_param()


if __name__ == '__main__':
    main()
