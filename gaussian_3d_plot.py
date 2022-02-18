import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import multivariate_normal

def main():
    x = np.arange(-30, 31)
    y = np.arange(-30, 31)

    X, Y = np.meshgrid(x, y)

    z_value = "./results/10^9(sr=0.5,dr=0.5)/" +\
              "tssp_topography(gatewidth=1)/tssp(gate=12+0,z=14+4).csv"
    z = np.loadtxt(z_value, delimiter=",")

    mean = np.zeros(2)
    sigma = np.array([[100, 0],[0, 100]])
    W = np.c_[np.ravel(X), np.ravel(Y)]
    Y_plot = multivariate_normal.pdf(x=W, mean=mean, cov=sigma)
    Y_plot = Y_plot.reshape(X.shape) * 130

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # plot 'z' for data, 'Y_plot' for function
    surf = ax.plot_surface(X, Y, Y_plot, cmap='hot', linewidth=0)
    fig.colorbar(surf, orientation='horizontal')
    ax.set_title("Surface Plot")
    plt.show()


if __name__ == '__main__':
    main()
