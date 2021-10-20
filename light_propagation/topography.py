import matplotlib.pyplot as plt
import numpy as np

def main():
    ssp = np.loadtxt("./results/ssp.csv", skiprows=1, delimiter=',')
    # reshape to (depth, rows, columns)
    ssp = np.reshape(ssp, (38, 61, 61))

    for z in range(38):
        image = np.empty_like(ssp[z])
        if ssp[z].max() != 0:
            image = ssp[z] / ssp[z].max()
        else:
            image = ssp[z]

        fig = plt.figure()
        plt.imshow(image, cmap='hot')
        plt.colorbar()
        plt.clim(0, 1)

        fig.savefig(f"./results/ssp{z}.png")

        plt.clf()
        plt.close()


if __name__ == '__main__':
    main()
