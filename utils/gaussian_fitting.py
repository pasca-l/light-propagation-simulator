import os
import argparse
import glob
import re
import numpy as np
import scipy.optimize as opt


class Fitter:
    def __init__(self, args):
        self.work_dir = f"../results/{args.dirname}/"
        self.dod_cam_paths = glob.glob(self.work_dir + "dOD*/camera/*.csv")
        self.dod_scan_paths = glob.glob(self.work_dir + "dOD*/scan/*.csv")
        self.ssp_paths = glob.glob(self.work_dir + "tssp_topography*/*.csv")
        self.all_paths = self.dod_cam_paths + self.dod_scan_paths +\
                         self.ssp_paths

        self.fit_dir = f"../results/{args.dirname}/fit/"
        os.makedirs(self.fit_dir, exist_ok=True)

        with open(self.fit_dir + "tssp_fit.csv", 'w') as f:
            f.write("gatewidth,gate,time,z,sigma_x,sigma_y,")
            f.write("centroid_x,centroid_y,amplitude\n")

        with open(self.fit_dir + "dod_camera_fit.csv", 'w') as f:
            f.write("gate,time,z,dmuar,pixel,sigma_x,sigma_y,")
            f.write("centroid_x,centroid_y,amplitude\n")

        with open(self.fit_dir + "dod_scan_fit.csv", 'w') as f:
            f.write("gate,time,z,dmuar,interval,sigma_x,sigma_y,")
            f.write("centroid_x,centroid_y,amplitude\n")

        self._get_model_param_from_summary()

    def _get_model_param_from_summary(self):
        with open(self.work_dir + "summary.txt") as f:
            txt = "".join(f.readlines())

        model_size_str = re.findall(r'model size\s*: \((.*)\)\n', txt).pop()
        model_size = list(map(int, re.findall(r'= (\d*)', model_size_str)))
        self.inputx, self.inputy, _ = model_size

    def twoD_gaussian(self, coord, cen_x, cen_y, sig_x, sig_y, amp):
        x, y = coord
        function = amp * np.exp(-0.5 * (((x - cen_x) / sig_x) ** 2 +
                                ((y - cen_y) / sig_y) ** 2))

        return function.ravel()

    def find_nearest_gaussian_param(self, data):
        xval, yval = np.meshgrid(np.linspace(0, self.inputx - 1, self.inputx),
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

    def find_gparam(self):
        for file in self.all_paths:
            info = re.findall(r'(?<=\().+?(?=\))', file)

            if file in self.ssp_paths:
                gate = int(info[1].split(',')[0].split('=')[1].split('+')[0])
                z = info[1].split(',')[1].split('=')[1]
                gatewidth = int(info[0].split('=')[1])
                time = gate * 250 + gatewidth * 125
            else:
                gate = int(info[0].split('=')[1])
                z = info[1].split(',')[0].split('=')[1]
                dmuar = info[1].split(',')[1].split('=')[1]
                interval = info[1].split(',')[2].split('=')[1]
                time = gate * 250 + 125

            map = np.loadtxt(file, delimiter=',')
            try:
                popt = self.find_nearest_gaussian_param(map)
            except:
                continue

            if file in self.ssp_paths:
                with open(self.fit_dir + "tssp_fit.csv", 'a') as f:
                    f.write(f"{gatewidth},{gate},{time},{z},{popt[2]},")
                    f.write(f"{popt[3]},{popt[0]},{popt[1]},{popt[4]}\n")

            else:
                if file in self.dod_cam_paths:
                    file_name = "camera"
                else:
                    file_name = "scan"

                with open(self.fit_dir + f"dod_{file_name}_fit.csv", 'a') as f:
                    f.write(f"{gate},{time},{z},{dmuar},{interval},{popt[2]},")
                    f.write(f"{popt[3]},{popt[0]},{popt[1]},{popt[4]}\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--dirname', type=str, default='data')

    fitter = Fitter(parser.parse_args())
    fitter.find_gparam()


if __name__ == '__main__':
    main()
