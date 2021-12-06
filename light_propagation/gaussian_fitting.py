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

        self.dmua_depth_init = 14
        self.dmua_depth = 4
        self.dmua_r_min = 1
        self.dmua_r_max = 15
        self.pixel_min = 1
        self.pixel_max = 5

    def twoD_gaussian(self, coord, cen_x, cen_y, sig_x, sig_y, amp):
        x, y = coord
        function = amp * np.exp(-0.5 * (((x - cen_x) / sig_x) ** 2 +
                                ((y - cen_y) / sig_y) ** 2))

        return function.ravel()

    def find_nearest_gaussian_param(self, data):
        xval, yval = np.meshgrid(
                        np.linspace(0, self.inputx - 1, self.inputx),
                        np.linspace(0, self.inputy - 1, self.inputy))

        # gaussian function parameters: cen_x, cen_y, sig_x, sig_y, amp
        initial_guess = (30, 30, 10, 10, 0.01)
        try:
            popt, _ = opt.curve_fit(self.twoD_gaussian, (xval, yval),
                                    data.ravel(), p0=initial_guess)
        except:
            print("No optical param")
            raise Exception

        return popt

    def tssp_gparam(self):
        with open(self.work_dir + "tssp_fit.csv", 'w') as f:
            f.write("gate,time,z,sigma_x,sigma_y,")
            f.write("centroid_x,centroid_y,amplitude\n")

        for gatenum in range(self.gate):
            tssp_map = self.work_dir + f"tssp_map/tssp(gate={gatenum}).csv"
            tssp = np.loadtxt(tssp_map, skiprows=2, delimiter=',')
            tssp = np.reshape(tssp,
                              (self.total_depth, self.inputx, self.inputy))

            temp = np.zeros((self.inputx, self.inputy))
            for z in range(self.dmua_depth_init,
                           self.dmua_depth_init + self.dmua_depth):
                temp += tssp[z]

            try:
                popt = self.find_nearest_gaussian_param(temp)
            except:
                continue

            with open(self.work_dir + "tssp_fit.csv", 'a') as f:
                f.write(f"{gatenum},{gatenum*250+125},")
                f.write(f"{self.dmua_depth_init}+{self.dmua_depth - 1},")
                f.write(f"{popt[2]},{popt[3]},")
                f.write(f"{popt[0]},{popt[1]},{popt[4]}\n")

    def dod_gparam(self):
        with open(self.work_dir + "/dOD/dOD_fit.csv", 'w') as f:
            f.write("z,dmuar,pixel,sigma_x,sigma_y,")
            f.write("centroid_x,centroid_y,amplitude\n")

        for pixelsize in range(self.pixel_min, self.pixel_max + 1):
            for r in range(self.dmua_r_min, self.dmua_r_max + 1):
                dod_map = self.work_dir +\
                          f"/dOD/dOD(z={self.dmua_depth_init}+" +\
                          f"{self.dmua_depth - 1},dmuar={r}," +\
                          f"pixel={pixelsize}).csv"
                dod = np.loadtxt(dod_map, delimiter=',')

                try:
                    popt = self.find_nearest_gaussian_param(dod)
                except:
                    continue

                with open(self.work_dir + "/dOD/dOD_fit.csv", 'a') as f:
                    f.write(f"{self.dmua_depth_init}+{self.dmua_depth - 1},")
                    f.write(f"{r},{pixelsize},{popt[2]},{popt[3]},")
                    f.write(f"{popt[0]},{popt[1]},{popt[4]}\n")


def main():
    fitter = Fitter(sys.argv[1])
    fitter.tssp_gparam()
    fitter.dod_gparam()


if __name__ == '__main__':
    main()
