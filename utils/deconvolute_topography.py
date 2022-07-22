import argparse
import matplotlib.pyplot as plt
import numpy as np


class DeconvTopography:
    def __init__(self, args):
        self.work_dir = f"../results/{args.dirname}/"

        self.dod_map = np.loadtxt(args.dod_map_path, delimiter=',')
        self.psf_map = np.loadtxt(args.psf_map_path, delimiter=',')

    def single_sig_gaussian(self, center=(29,29), size=(61,61), sig=12):
        x_axis = np.linspace(0, size[0]-1, size[0]) - center[0]
        y_axis = np.linspace(0, size[1]-1, size[1]) - center[1]
        xval, yval = np.meshgrid(x_axis, y_axis)
        kernel = np.exp(-0.5 * (xval ** 2 + yval ** 2) / sig ** 2)

        return kernel

    def topography_of_deconv(self):
        f_dod_map = np.fft.fft2(self.dod_map)
        f_psf_map = np.fft.fft2(self.psf_map)

        deconv = f_dod_map / f_psf_map
        if_dod_map = np.fft.ifft2(deconv).real

        if if_dod_map.max() != 1:
            if_dod_map /= if_dod_map.max()

        plt.figure()
        plt.imshow(if_dod_map, cmap='hot')
        plt.colorbar()
        plt.clim(0, 1)

        plt.clf()
        plt.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dirname', type=str, default='data')
    parser.add_argument('--dod_map_path', type=str, required=True)
    parser.add_argument('--psf_map_path', type=str, required=True)

    deconv = DeconvTopography(parser.parse_args())
    deconv.topography_of_deconv()


if __name__ == '__main__':
    main()
