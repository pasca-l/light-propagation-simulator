import os
import sys
import matplotlib.pyplot as plt
import numpy as np

class SspTopography:
    def __init__(self, dirname):
        self.work_dir = f"./results/{dirname}/"

        self.inputx = 61
        self.inputy = 61
        self.total_depth = 28
        self.gate = 16
        self.depth_init = 14
        self.depth = 4

    def make_topography(self, ssp_map, save_name, csv_flag=False):
        if csv_flag:
            np.savetxt(save_name + ".csv", ssp_map, delimiter=',', fmt='%lf')

        image = np.zeros_like(ssp_map)
        if ssp_map.max() != 0:
            image = ssp_map / ssp_map.max()

        fig = plt.figure()
        plt.imshow(image, cmap='hot')
        plt.colorbar()
        plt.clim(0, 1)

        fig.savefig(save_name + ".png")

        plt.clf()
        plt.close()

    def topography_of_ssp(self):
        ssp_map = self.work_dir + "ssp.csv"
        ssp = np.loadtxt(ssp_map, skiprows=1, delimiter=',')
        # reshape to (depth, rows, columns)
        ssp = np.reshape(ssp, (self.total_depth, self.inputx, self.inputy))

        # make folder to save topography images
        ssp_topo_dir = self.work_dir + "ssp_topography/"
        if not os.path.exists(ssp_topo_dir):
            os.makedirs(ssp_topo_dir)

        for z in range(self.total_depth):
            self.make_topography(ssp[z], ssp_topo_dir + f"ssp(z={z})")

    def topography_of_tssp(self):
        tssp_topo_dir = self.work_dir + "tssp_topography/"
        if not os.path.exists(tssp_topo_dir):
            os.makedirs(tssp_topo_dir)

        for gatenum in range(self.gate):
            tssp_map = self.work_dir + f"tssp_map/tssp(gate={gatenum}).csv"
            tssp = np.loadtxt(tssp_map, skiprows=2, delimiter=',')
            tssp = np.reshape(tssp,
                              (self.total_depth, self.inputx, self.inputy))

            psf_map = np.zeros((self.inputx, self.inputy))
            for z in range(self.depth_init, self.depth_init + self.depth):
                psf_map += tssp[z]

            self.make_topography(psf_map, tssp_topo_dir +\
                                 f"tssp(gate={gatenum},z={self.depth_init}+" +\
                                 f"{self.depth})")


def main():
    topo = SspTopography(sys.argv[1])
    topo.topography_of_tssp()


if __name__ == '__main__':
    main()
