import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import multivariate_normal


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--filepath', type=str, required=True)

    topo_path = parser.parse_args().filepath

    x, y = np.meshgrid(np.arange(-30, 31), np.arange(-30, 31))
    z = np.loadtxt(topo_path, delimiter=",")

    # show plot for data
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(x, y, z, cmap='hot', linewidth=0)
    fig.colorbar(surf, orientation='horizontal')
    ax.set_title("Surface Plot")
    plt.show()

    # show plot for gaussian distribution
    # mean = np.zeros(2)
    # sigma = np.array([[100, 0],[0, 100]])
    # w = np.c_[np.ravel(x), np.ravel(y)]
    # plot = multivariate_normal.pdf(x=w, mean=mean, cov=sigma)
    # plot = plot.reshape(x.shape) * 130

    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # surf = ax.plot_surface(x, y, plot, cmap='hot', linewidth=0)
    # fig.colorbar(surf, orientation='horizontal')
    # ax.set_title("Surface Plot")
    # plt.show()


if __name__ == '__main__':
    main()
