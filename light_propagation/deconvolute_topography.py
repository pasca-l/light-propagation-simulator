import matplotlib.pyplot as plt
import numpy as np

def single_sig_gaussian(center=(29,29), size=(61,61), sig=12):
    x_axis = np.linspace(0, size[0]-1, size[0]) - center[0]
    y_axis = np.linspace(0, size[1]-1, size[1]) - center[1]
    xval, yval = np.meshgrid(x_axis, y_axis)
    xxval, yyval = np.square(xval), np.square(yval)
    kernel = np.exp(-0.5 * (xxval + yyval) / np.square(sig))

    return kernel


def topography_of_deconv():
    dirname = "10^9(r=2.6)"
    work_dir = f"./results/{dirname}/"

    dod_map = np.loadtxt(work_dir + "dOD/dOD.csv", delimiter=',')
    psf_map = np.loadtxt(work_dir + "dOD/PSF.csv", delimiter=',')

    f_dod_map = np.fft.fft2(dod_map)
    f_psf_map = np.fft.fft2(psf_map)

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
    topography_of_deconv()


if __name__ == '__main__':
    main()
