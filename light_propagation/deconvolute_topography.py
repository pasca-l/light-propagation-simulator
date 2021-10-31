import matplotlib.pyplot as plt
import numpy as np

def main():
    dirname = "10^9(r=2.6)"
    work_dir = f"./results/{dirname}/"

    dod_map = np.loadtxt(work_dir + "dOD/dOD.csv", delimiter=',')
    psf_map = np.loadtxt(work_dir + "dOD/PSF.csv", delimiter=',')

    f_dod_map = np.fft.fft2(dod_map)
    f_psf_map = np.fft.fft2(psf_map)

    deconv = f_dod_map / f_psf_map
    if_dod_map = np.fft.ifft2(deconv).real

    if if_dod_map.max() != 0:
        if_dod_map = if_dod_map / if_dod_map.max()

    fig = plt.figure()
    plt.imshow(if_dod_map, cmap='hot')
    plt.colorbar()
    plt.clim(0, 1)
    fig.savefig(work_dir + "dOD/deconv.png")

    plt.clf()
    plt.close()


if __name__ == '__main__':
    main()
