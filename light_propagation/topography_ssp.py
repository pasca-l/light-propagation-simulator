import os
import matplotlib.pyplot as plt
import numpy as np

class SspTopography:
    def __init__(self, dirname):
        self.work_dir = f"./results/{dirname}/"

        self.inputx = 61
        self.inputy = 61
        self.total_depth = 28
        self.gate = 16

        # make folder to save topography images
        self.ssp_topo_dir = self.work_dir + "ssp_topography/"
        if not os.path.exists(self.ssp_topo_dir):
            os.makedirs(self.ssp_topo_dir)

        self.tssp_topo_dir = self.work_dir + "tssp_topography/"
        if not os.path.exists(self.tssp_topo_dir):
            os.makedirs(self.tssp_topo_dir)

    def make_topography(self, z, ssp_map, save_name):
        image = np.zeros_like(ssp_map[z])
        if ssp_map[z].max() != 0:
            image = ssp_map[z] / ssp_map[z].max()

        fig = plt.figure()
        plt.imshow(image, cmap='hot')
        plt.colorbar()
        plt.clim(0, 1)

        fig.savefig(save_name)

        plt.clf()
        plt.close()

    def topography_of_ssp(self):
        ssp_map = self.work_dir + "ssp.csv"
        ssp = np.loadtxt(ssp_map, skiprows=1, delimiter=',')
        # reshape to (depth, rows, columns)
        ssp = np.reshape(ssp, (self.total_depth, self.inputx, self.inputy))

        for z in range(self.total_depth):
            self.make_topography(z, ssp,
                                 self.ssp_topo_dir + f"ssp(z={z}).png")

    def topography_of_tssp(self):
        for gatenum in range(self.gate):
            tssp_map = self.work_dir + f"tssp_map/tssp(gate={gatenum}).csv"
            tssp = np.loadtxt(tssp_map, skiprows=2, delimiter=',')
            tssp = np.reshape(tssp,
                              (self.total_depth, self.inputx, self.inputy))

            for z in range(self.total_depth):
                self.make_topography(z, tssp,
                                     self.tssp_topo_dir +
                                     f"tssp(gate={gatenum},z={z}).png")


def main():
    topo = SspTopography("10^9(r=0.5)")
    topo.topography_of_ssp()
    topo.topography_of_tssp()


if __name__ == '__main__':
    main()
