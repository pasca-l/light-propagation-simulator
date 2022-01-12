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
        self.total_gate = 16
        self.gate_width = 1

        self.dmua_depth_init = 14
        self.dmua_depth = 4
        self.dmua_r_min = 1
        self.dmua_r_max = 15
        self.pixel_min = 1
        self.pixel_max = 6
        self.interval_min = 1
        self.interval_max = 6

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
        tssp_topo_dir = self.work_dir +\
                        f"tssp_topography(gatewidth={self.gate_width})/"

        with open(tssp_topo_dir + "tssp_fit.csv", 'w') as f:
            f.write("gate,time,z,sigma_x,sigma_y,")
            f.write("centroid_x,centroid_y,amplitude\n")

        for i, gatenum in enumerate(range(0, self.total_gate, self.gate_width)):
            tssp_map = self.work_dir +\
                       f"tssp_topography(gatewidth={self.gate_width})/" +\
                       f"tssp(gate={gatenum}+{self.gate_width - 1}," +\
                       f"z={self.dmua_depth_init}+{self.dmua_depth}).csv"
            tssp = np.loadtxt(tssp_map, delimiter=',')

            try:
                popt = self.find_nearest_gaussian_param(tssp)
            except:
                continue

            with open(tssp_topo_dir + "tssp_fit.csv", 'a') as f:
                time = i*self.gate_width*250 + (self.gate_width*125)
                f.write(f"{gatenum},{time},")
                f.write(f"{self.dmua_depth_init}+{self.dmua_depth - 1},")
                f.write(f"{popt[2]},{popt[3]},")
                f.write(f"{popt[0]},{popt[1]},{popt[4]}\n")

    def camera_dod_gparam(self, dod_gate):
        dod_dir = self.work_dir + f"camera/dOD(gate={dod_gate})/"

        with open(dod_dir + "dOD_fit.csv", 'w') as f:
            f.write("z,dmuar,pixel,sigma_x,sigma_y,")
            f.write("centroid_x,centroid_y,amplitude\n")

        for pixelsize in range(self.pixel_min, self.pixel_max + 1):
            for r in range(self.dmua_r_min, self.dmua_r_max + 1):
                dod_map = dod_dir +\
                          f"dOD(z={self.dmua_depth_init}+" +\
                          f"{self.dmua_depth - 1},dmuar={r}," +\
                          f"pixel={pixelsize}).csv"
                dod = np.loadtxt(dod_map, delimiter=',')

                try:
                    popt = self.find_nearest_gaussian_param(dod)
                except:
                    continue

                with open(dod_dir + "dOD_fit.csv", 'a') as f:
                    f.write(f"{self.dmua_depth_init}+{self.dmua_depth - 1},")
                    f.write(f"{r},{pixelsize},{popt[2]},{popt[3]},")
                    f.write(f"{popt[0]},{popt[1]},{popt[4]}\n")

    def scan_dod_gparam(self, dod_gate):
        dod_dir = self.work_dir + f"scan/dOD(gate={dod_gate})/"

        with open(dod_dir + "dOD_fit.csv", 'w') as f:
            f.write("z,dmuar,interval,sigma_x,sigma_y,")
            f.write("centroid_x,centroid_y,amplitude\n")

        for int in range(self.interval_min, self.interval_max + 1):
            for r in range(self.dmua_r_min, self.dmua_r_max + 1):
                dod_map = dod_dir +\
                          f"dOD(z={self.dmua_depth_init}+" +\
                          f"{self.dmua_depth - 1},dmuar={r}," +\
                          f"interval={int}).csv"
                dod = np.loadtxt(dod_map, delimiter=',')

                try:
                    popt = self.find_nearest_gaussian_param(dod)
                except:
                    continue

                with open(dod_dir + "dOD_fit.csv", 'a') as f:
                    f.write(f"{self.dmua_depth_init}+{self.dmua_depth - 1},")
                    f.write(f"{r},{int},{popt[2]},{popt[3]},")
                    f.write(f"{popt[0]},{popt[1]},{popt[4]}\n")

    # # illustrates the best-fit 2d gaussian distribution
    # def draw_nearest_gaussian(self, dod_map):
    #     dod = np.loadtxt(dod_map, delimiter=',')
    #
    #     try:
    #         popt = self.find_nearest_gaussian_param(dod)
    #     except:
    #         return


def main():
    fitter = Fitter(sys.argv[1])

    # fitting tssp topography
    fitter.tssp_gparam()

    # fitting dod topography
    # gates = [6, 15]
    # for gate in gates:
    #     fitter.dod_gparam(gate)


if __name__ == '__main__':
    main()
