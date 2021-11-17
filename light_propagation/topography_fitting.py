import sys
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt

class Fitter:
    def __init__(self, dirname):
        self.work_dir = f"./results/{dirname}/"
        self.dod_map = np.loadtxt(self.work_dir + f"dOD/dOD.csv",
                                  delimiter=',')

        self.size_x, self.size_y = self.dod_map.shape

    def twoD_gaussian(self, coord, cen_x, cen_y, sig_x, sig_y, amp):
        x, y = coord
        function = amp * np.exp(-0.5 * (((x - cen_x) / sig_x) ** 2 +
                                ((y - cen_y) / sig_y) ** 2))

        return function.ravel()

    def draw_nearest_gaussian(self):
        xval, yval = np.meshgrid(np.linspace(0, self.size_x - 1, self.size_x),
                                 np.linspace(0, self.size_y - 1, self.size_y))

        # gaussian function parameters: cen_x, cen_y, sig_x, sig_y, amp
        initial_guess = (30, 30, 10, 10, 0.01)
        popt, _ = opt.curve_fit(self.twoD_gaussian, (xval, yval),
                                self.dod_map.ravel(), p0=initial_guess)

        print(f"centroid_x : {popt[0]},")
        print(f"centroid_y : {popt[1]},")
        print(f"sigma_x : {popt[2]},")
        print(f"sigma_y : {popt[3]},")
        print(f"amplitude : {popt[4]}")

        best_fit = self.twoD_gaussian((xval, yval), *popt)
        best_fit = best_fit.reshape(self.size_x, self.size_y)
        if best_fit.max() != 1:
            best_fit = best_fit / best_fit.max()

        fig = plt.figure()
        plt.imshow(best_fit, cmap='hot')
        plt.colorbar()
        plt.clim(0, 1)
        fig.savefig(self.work_dir + "dOD/best_fit.png")

        plt.clf()
        plt.close()


def main():
    fitter = Fitter(sys.argv[1])
    fitter.draw_nearest_gaussian()


if __name__ == '__main__':
    main()
