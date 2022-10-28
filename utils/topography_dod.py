import os
import argparse
import re
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate


class DodTopography:
    def __init__(self, args):
        self.work_dir = f"../results/{args.dirname}/"
        self.gate = args.gate_at
        self.dmua_map_check = args.map_check

        self.dmua_depth_init = args.dmua_init
        self.dmua_depth = args.dmua_depth
        self.dmua = args.dmua
        self.dmua_r = args.dmua_r

        self.pixel = args.pixel
        self.interval = args.interval

        self._get_model_param_from_summary()

        ssp_map = np.loadtxt(self.work_dir +
                             f"tssp_map/tssp(gate={self.gate}).csv",
                             skiprows=2, delimiter=',')
        self.ssp = ssp_map.reshape(self.total_depth, self.inputx, self.inputy)

        self.dod_dir = self.work_dir + f"dOD(gate={self.gate})/"
        if not os.path.exists(self.dod_dir):
            os.makedirs(self.dod_dir)

        if self.dmua_map_check:
            self.dmua_dir = self.dod_dir + "dmua_map/"
            if not os.path.exists(self.dmua_dir):
                os.makedirs(self.dmua_dir)

    def _get_model_param_from_summary(self):
        with open(self.work_dir + "summary.txt") as f:
            txt = "".join(f.readlines())

        model_size_str = re.findall(r'model size\s*: \((.*)\)\n', txt).pop()
        model_size = list(map(int, re.findall(r'= (\d*)', model_size_str)))
        self.inputx, self.inputy, self.total_depth = model_size

    def add_surrounding_zeros(self, matrix):
        row, col = map(int, matrix.shape)

        zeros1 = np.zeros((int(row/2), col))
        temp = np.vstack((matrix, zeros1))
        temp = np.vstack((zeros1, temp))

        zeros2 = np.zeros((2*row-1, int(col/2)))
        temp = np.hstack((temp, zeros2))
        output = np.hstack((zeros2, temp))

        return output

    def create_dmua_map(self, ssp_map):
        r = self.dmua_r
        dmua_map = np.zeros_like(ssp_map)
        ymid, xmid = map(lambda f: int(f/2), dmua_map.shape)

        # make circle at left upper corner
        for y in range(2*r+1):
            for x in range(2*r+1):
                if (x-r)**2 + (y-r)**2 <= r*r:
                    dmua_map[y][x] = self.dmua

        # move to center
        dmua_map = np.roll(dmua_map, ymid - r, axis=0)
        dmua_map = np.roll(dmua_map, xmid - r, axis=1)
        # printing area of interest for check
        if self.dmua_map_check:
            print(dmua_map[ymid-r:ymid+r+1, xmid-r:xmid+r+1])

        return dmua_map

    def convolute_maps(self, dod_map, dmua_map, ssp_with_buffer):
        xcenter, ycenter = map(lambda f: int(f/2), ssp_with_buffer.shape)
        half_inputx, half_inputy = int(self.inputx/2), int(self.inputy/2)

        for dety in range(half_inputy, half_inputy + ycenter + 1):
            for detx in range(half_inputx, half_inputx + xcenter + 1):
                temp = np.roll(dmua_map, ycenter - dety - 1, axis=0)
                relative_dmua_map = np.roll(temp, xcenter - detx - 1, axis=1)

                i, j = dety - half_inputy, detx - half_inputx
                dod_map[i][j] = np.sum(ssp_with_buffer * relative_dmua_map)

        return dod_map

    def replace_to_representative(self, dod_map):
        offset = self.pixel // 2
        for y in range(offset, self.inputy, self.pixel):
            for x in range(offset, self.inputx, self.pixel):
                rep = dod_map[y][x]
                if self.pixel % 2 == 1:
                    dod_map[y-offset:y+offset+1, x-offset:x+offset+1] = rep
                elif self.pixel % 2 == 0:
                    dod_map[y-offset:y+offset, x-offset:x+offset] = rep

        return dod_map

    def spline_complement(self, dod_map):
        data_xaxis = np.arange(0, dod_map.shape[0], self.interval)
        data_yaxis = np.arange(0, dod_map.shape[1], self.interval)
        data_xgrid, data_ygrid = np.meshgrid(data_xaxis, data_yaxis)
        data = dod_map[data_xgrid, data_ygrid]

        f = interpolate.interp2d(data_xaxis, data_yaxis, data, kind='cubic')

        map_xaxis = np.arange(dod_map.shape[0])
        map_yaxis = np.arange(dod_map.shape[1])
        dod_map = f(map_xaxis, map_yaxis)

        return dod_map

    def make_topography(self, dod_map, save_name, norm_flag=False, raw=False):
        image = np.zeros_like(dod_map)
        image = dod_map
        if norm_flag:
            if dod_map.max() != 1:
                image = dod_map / dod_map.max()

        # save image without axis and legend
        if raw:
            plt.imsave(save_name + "_img.png", image, cmap='jet')

        # save image with axis and legend
        fig = plt.figure()
        plt.imshow(image, cmap='jet')
        plt.colorbar()
        if image.max() == 0:
            plt.clim(0, 1)
        else:
            plt.clim(0, image.max())

        fig.savefig(save_name + ".png")

        plt.clf()
        plt.close()

    def topography_of_dod(self, camera=True, spline=True, csv_flag=False):
        ssp_with_buffer = np.zeros((2 * self.inputx - 1, 2 * self.inputy - 1))
        for z in range(self.dmua_depth):
            ssp_with_buffer +=\
                self.add_surrounding_zeros(self.ssp[self.dmua_depth_init + z])

        dmua_map = self.create_dmua_map(ssp_with_buffer)
        zeros = np.zeros((self.inputx, self.inputy))
        dod_map = self.convolute_maps(zeros, dmua_map, ssp_with_buffer)

        if camera:
            camera_dir = self.dod_dir + "camera/"
            if not os.path.exists(camera_dir):
                os.makedirs(camera_dir)

            dod_topo = self.replace_to_representative(dod_map)

            save_as = camera_dir + f"dOD(z={self.dmua_depth_init}+" +\
                      f"{self.dmua_depth - 1},dmuar={self.dmua_r}," +\
                      f"pixel={self.pixel})"
            self.make_topography(dod_topo, save_as)

        if spline:
            scan_dir = self.dod_dir + "scan/"
            if not os.path.exists(scan_dir):
                os.makedirs(scan_dir)

            dod_topo = self.spline_complement(dod_map)

            save_as = scan_dir + f"dOD(z={self.dmua_depth_init}+" +\
                      f"{self.dmua_depth - 1},dmuar={self.dmua_r}," +\
                      f"interval={self.interval})"

        if csv_flag:
            np.savetxt(f"{save_as}.csv", dod_map, delimiter=',', fmt='%lf')
        self.make_topography(dod_topo, save_as)

    def topography_of_psf(self, csv_flag=False):
        psf_map = np.zeros((self.inputx, self.inputy))
        for z in range(self.dmua_depth):
            psf_map += self.ssp[self.dmua_depth_init + z]

        psf_mua_map = np.empty_like(psf_map)
        psf_mua_map.fill(self.dmua)
        psf_map = psf_map * psf_mua_map

        save_as = self.dod_dir +\
                  f"PSF(z={self.dmua_depth_init}+{self.dmua_depth - 1})"
        if csv_flag:
            np.savetxt(f"{save_as}.csv", psf_map, delimiter=',', fmt='%lf')
        self.make_topography(psf_map, save_as)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--dirname', type=str, default='data')
    parser.add_argument('-g', '--gate_at', type=int, default=12)
    parser.add_argument('-p', '--pixel', type=int, default=1)
    parser.add_argument('-i', '--interval', type=int, default=1)
    parser.add_argument('-r', '--dmua_r', type=int, default=1)
    parser.add_argument('-z', '--dmua_init', type=int, default=14)
    parser.add_argument('-w', '--dmua_depth', type=int, default=4)
    parser.add_argument('-m', '--dmua', type=float, default=0.002)
    parser.add_argument('-c', '--map_check', action='store_true')

    topo = DodTopography(parser.parse_args())
    topo.topography_of_dod()
    topo.topography_of_psf()


if __name__ == '__main__':
    main()
