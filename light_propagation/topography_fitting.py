import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt

class Fitter:
    def __init__(self, dirname):
        self.work_dir = f"./results/{dirname}/"
        self.dod_map = np.loadtxt(self.work_dir + f"dOD/dOD.csv",
                                  delimiter=',')

        self.size_x, self.size_y = self.dod_map.shape
        self.xval, self.yval = np.meshgrid(
                                np.linspace(0, self.size_x - 1, self.size_x),
                                np.linspace(0, self.size_y - 1, self.size_y))

        # gaussian function parameters: cen_x, cen_y, sig_x, sig_y, amp
        self.initial_guess = (30, 30, 10, 10, 0.01)

    def twoD_gaussian(self, coord, cen_x, cen_y, sig_x, sig_y, amp):
        x, y = coord
        function = amp * np.exp(-0.5 * (((x - cen_x) / sig_x) ** 2 +
                                ((y - cen_y) / sig_y) ** 2))

        return function.ravel()

    def find_nearest_gaussian(self):
        popt, _ = opt.curve_fit(self.twoD_gaussian,
                                (self.xval, self.yval),
                                self.dod_map.ravel(),
                                p0=self.initial_guess)

        print(popt)
        best_fit = self.twoD_gaussian((self.xval, self.yval), *popt)
        best_fit = best_fit.reshape(self.size_x, self.size_y)
        if best_fit.max() != 1:
            best_fit = best_fit / best_fit.max()

        fig = plt.figure()
        plt.imshow(best_fit, cmap='hot')
        plt.colorbar()
        plt.clim(0, 1)
        fig.savefig(self.work_dir + "dOD/fit.png")

        plt.clf()
        plt.close()

        return


def main():
    fitter = Fitter("10^9(r=2.6)")
    fitter.find_nearest_gaussian()


if __name__ == '__main__':
    main()
