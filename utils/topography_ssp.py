import os
import argparse
import matplotlib.pyplot as plt
import numpy as np


class SspTopography:
    def __init__(self, args):
        self.work_dir = f"../results/{args.dirname}/"

        self.inputx = 61
        self.inputy = 61
        self.total_depth = 28
        self.total_gate = 16

        self.gate_width = 1
        self.depth_init = 14
        self.depth = 4

    def make_topography(self, ssp_map, save_name, csv_flag=False,
                        norm_flag=False):
        if csv_flag:
            np.savetxt(save_name + ".csv", ssp_map, delimiter=',', fmt='%lf')

        image = np.zeros_like(ssp_map)
        image = ssp_map
        if norm_flag:
            if ssp_map.max() != 0:
                image = ssp_map / ssp_map.max()

        fig = plt.figure()
        plt.imshow(image, cmap='hot')
        plt.colorbar()
        if image.max() == 0:
            plt.clim(0, 1)
        else:
            plt.clim(0, image.max())

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
        tssp_topo_dir = self.work_dir +\
                        f"tssp_topography(gatewidth={self.gate_width})/"
        if not os.path.exists(tssp_topo_dir):
            os.makedirs(tssp_topo_dir)

        for gatenum in range(0, self.total_gate, self.gate_width):
            tssp = np.zeros((self.total_depth, self.inputx, self.inputy))

            for gatewidth in range(self.gate_width):
                tssp_map = self.work_dir +\
                           f"tssp_map/tssp(gate={gatenum + gatewidth}).csv"
                temp = np.loadtxt(tssp_map, skiprows=2, delimiter=',')
                temp = np.reshape(temp,
                                  (self.total_depth, self.inputx, self.inputy))
                tssp += temp

            psf_map = np.zeros((self.inputx, self.inputy))
            for z in range(self.depth_init, self.depth_init + self.depth):
                psf_map += tssp[z]

            self.make_topography(psf_map, tssp_topo_dir +\
                                 f"tssp(gate={gatenum}+{gatewidth}," +\
                                 f"z={self.depth_init}+{self.depth})",
                                 csv_flag=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dirname', type=str, default='data')

    topo = SspTopography(parser.parse_args())
    topo.topography_of_tssp()


if __name__ == '__main__':
    main()
