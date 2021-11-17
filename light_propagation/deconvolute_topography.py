import sys
import matplotlib.pyplot as plt
import numpy as np

class DeconvTopography:
    def __init__(self, dirname):
        self.work_dir = f"./results/{dirname}/"

        self.dod_map = np.loadtxt(self.work_dir + "dOD/dOD.csv", delimiter=',')
        self.psf_map = np.loadtxt(self.work_dir + "dOD/PSF.csv", delimiter=',')


    def single_sig_gaussian(self, center=(29,29), size=(61,61), sig=12):
        x_axis = np.linspace(0, size[0]-1, size[0]) - center[0]
        y_axis = np.linspace(0, size[1]-1, size[1]) - center[1]
        xval, yval = np.meshgrid(x_axis, y_axis)
        kernel = np.exp(-0.5 * (xval ** 2 + yyval ** 2) / sig ** 2)

        return kernel

    def topography_of_deconv(self):
        f_dod_map = np.fft.fft2(self.dod_map)
        f_psf_map = np.fft.fft2(self.psf_map)

        deconv = f_dod_map / f_psf_map
        if_dod_map = np.fft.ifft2(deconv).real

        if if_dod_map.max() != 1:
            if_dod_map = if_dod_map / if_dod_map.max()

        fig = plt.figure()
        plt.imshow(if_dod_map, cmap='hot')
        plt.colorbar()
        plt.clim(0, 1)
        fig.savefig(work_dir + "dOD/deconv.png")

        plt.clf()
        plt.close()


def main():
    deconv = DeconvTopography(sys.argv[1])
    topography_of_deconv()


if __name__ == '__main__':
    main()
