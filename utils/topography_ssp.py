import os
import argparse
import re
import matplotlib.pyplot as plt
import numpy as np


class SspTopography:
    def __init__(self, args):
        self.work_dir = f"../results/{args.dirname}/"

        self.gate_width = args.gate_width
        self.attn_layer = args.layer_num
        self._get_model_param_from_summary()

    def _get_model_param_from_summary(self):
        with open(self.work_dir + "summary.txt") as f:
            txt = "".join(f.readlines())

        model_size_str = re.findall(r'model size\s*: \((.*)\)\n', txt).pop()
        model_size = list(map(int, re.findall(r'= (\d*)', model_size_str)))
        self.inputx, self.inputy, self.total_depth = model_size

        layer_w_str = re.findall(r'\d* \d* \d\.\d* \d\.\d*\n', txt)
        layer_w = [int(i.split()[1]) for i in layer_w_str]
        self.depth = layer_w[self.attn_layer]
        self.depth_init = sum(layer_w[:self.attn_layer])

        gate_num_str = re.findall(r'TSSP gate number\s*: (\d*)', txt).pop()
        self.total_gate = int(gate_num_str)

    def _make_topography(self, ssp_map, save_name, norm_flag=False):
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
        # make folder to save topography images
        ssp_topo_dir = self.work_dir + "ssp_topography/"
        if not os.path.exists(ssp_topo_dir):
            os.makedirs(ssp_topo_dir)

        ssp_map = self.work_dir + "ssp.csv"
        ssp = np.loadtxt(ssp_map, skiprows=1, delimiter=',')
        # reshape to (depth, rows, columns)
        ssp = ssp.reshape(self.total_depth, self.inputx, self.inputy)

        for z in range(self.total_depth):
            self._make_topography(ssp[z], ssp_topo_dir + f"ssp(z={z})")

    def topography_of_tssp(self, csv_flag=False):
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
                temp = temp.reshape(self.total_depth, self.inputx, self.inputy)
                tssp += temp

            psf_map = np.zeros((self.inputx, self.inputy))
            for z in range(self.depth_init, self.depth_init + self.depth):
                psf_map += tssp[z]

            filename = tssp_topo_dir +\
                        f"tssp(gate={gatenum}+{self.gate_width - 1}," +\
                        f"z={self.depth_init}+{self.depth})"

            if csv_flag:
                np.savetxt(f"{filename}.csv", psf_map, delimiter=',', fmt='%lf')
            self._make_topography(psf_map, filename)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--dirname', type=str, default='data')
    parser.add_argument('-w', '--gate_width', type=int, default=1)
    parser.add_argument('-l', '--layer_num', type=int, default=3)

    topo = SspTopography(parser.parse_args())
    topo.topography_of_tssp()


if __name__ == '__main__':
    main()
