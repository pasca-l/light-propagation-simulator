import os
import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate

class DodTopography:
    def __init__(self, dirname, time_gate):
        self.work_dir = f"./results/{dirname}/"
        self.gate = time_gate
        self.dmua_map_check = False

        self.inputx = 61
        self.inputy = 61
        self.total_depth = 28
        self.dmua_depth_init = 14
        self.dmua_depth = 4
        self.dmua = 0.002
        self.dmua_r_min = 1
        self.dmua_r_max = 15

        self.pixel_min = 1
        self.pixel_max = 6
        self.interval_min = 1
        self.interval_max = 6

        self.ssp_map = np.loadtxt(self.work_dir +
                                  f"tssp_map/tssp(gate={self.gate}).csv",
                                  skiprows=2, delimiter=',')
        self.ssp = np.reshape(self.ssp_map,
                              (self.total_depth, self.inputx, self.inputy))

        self.dod_dir = self.work_dir + f"dOD(gate={self.gate})/"
        if not os.path.exists(self.dod_dir):
            os.makedirs(self.dod_dir)

        if self.dmua_map_check:
            self.dmua_dir = self.dod_dir + "dmua_map/"
            if not os.path.exists(self.dmua_dir):
                os.makedirs(self.dmua_dir)

    def add_surrounding_zeros(self, matrix):
        row, col = map(int, matrix.shape)

        zeros1 = np.zeros((int(row/2), col))
        temp = np.vstack((matrix, zeros1))
        temp = np.vstack((zeros1, temp))

        zeros2 = np.zeros((2*row-1, int(col/2)))
        temp = np.hstack((temp, zeros2))
        output = np.hstack((zeros2, temp))

        return output

    def create_dmua_map(self, ssp_map, r):
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

                if self.dmua_map_check:
                    np.savetxt(self.dmua_dir + f"/dmua({i},{j}).csv",
                               relative_dmua_map[29:90, 29:90],
                               delimiter=',', fmt='%lf')

        return dod_map

    def replace_to_representative(self, dod_map, pixelsize):
        offset = pixelsize // 2
        for y in range(offset, self.inputy, pixelsize):
            for x in range(offset, self.inputx, pixelsize):
                rep = dod_map[y][x]
                if pixelsize % 2 == 1:
                    dod_map[y-offset:y+offset+1, x-offset:x+offset+1] = rep
                elif pixelsize % 2 == 0:
                    dod_map[y-offset:y+offset, x-offset:x+offset] = rep

        return dod_map

    def spline_complement(self, dod_map, intervalsize):
        data_xaxis = np.arange(0, dod_map.shape[0], intervalsize)
        data_yaxis = np.arange(0, dod_map.shape[1], intervalsize)
        data_xgrid, data_ygrid = np.meshgrid(data_xaxis, data_yaxis)
        data = dod_map[data_xgrid, data_ygrid]

        f = interpolate.interp2d(data_xaxis, data_yaxis, data, kind='cubic')

        # image of scanned measurement points
        # mask = np.zeros_like(dod_map)
        # for y in range(0, mask.shape[1], intervalsize):
        #     for x in range(0, mask.shape[0], intervalsize):
        #         data_xaxis =
        #         mask[x, y] = 1
        # mask = mask * dod_map

        map_xaxis = np.arange(dod_map.shape[0])
        map_yaxis = np.arange(dod_map.shape[1])
        map_xgrid, map_ygrid = np.meshgrid(map_xaxis, map_yaxis)
        dod_map = f(map_xaxis, map_yaxis)

        return dod_map

    def make_topography(self, dod_map, save_name, csv_flag=False,
                        norm_flag=False, raw_image=False):
        if csv_flag:
            np.savetxt(save_name + ".csv", dod_map, delimiter=',', fmt='%lf')

        image = np.zeros_like(dod_map)
        image = dod_map
        if norm_flag:
            if dod_map.max() != 1:
                image = dod_map / dod_map.max()

        # save image without axis and legend
        if raw_image:
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

    def topography_of_dod(self, camera=True, spline=True):
        ssp_with_buffer = np.zeros((2 * self.inputx - 1, 2 * self.inputy - 1))
        for z in range(self.dmua_depth_init,
                       self.dmua_depth_init + self.dmua_depth):
            ssp_with_buffer += self.add_surrounding_zeros(self.ssp[z])

        for r in range(self.dmua_r_min, self.dmua_r_max + 1):
            dmua_map = self.create_dmua_map(ssp_with_buffer, r)
            zeros = np.zeros((self.inputx, self.inputy))

            if camera:
                dod_map = self.convolute_maps(zeros, dmua_map, ssp_with_buffer)

                camera_dir = self.dod_dir + "camera/"
                if not os.path.exists(camera_dir):
                    os.makedirs(camera_dir)

                for pixel in range(self.pixel_min, self.pixel_max + 1):
                    dod_topo = self.replace_to_representative(dod_map, pixel)

                    save_as = camera_dir + f"dOD(z={self.dmua_depth_init}+" +\
                             f"{self.dmua_depth - 1},dmuar={r},pixel={pixel})"
                    self.make_topography(dod_topo, save_as)

            if spline:
                dod_map = self.convolute_maps(zeros, dmua_map, ssp_with_buffer)

                scan_dir = self.dod_dir + "scan/"
                if not os.path.exists(scan_dir):
                    os.makedirs(scan_dir)

                for int in range(self.interval_min, self.interval_max + 1):
                    dod_topo = self.spline_complement(dod_map, int)

                    save_as = scan_dir + f"dOD(z={self.dmua_depth_init}+" +\
                            f"{self.dmua_depth - 1},dmuar={r},interval={int})"
                    self.make_topography(dod_topo, save_as)

    def topography_of_psf(self):
        psf_map = np.zeros((self.inputx, self.inputy))
        for z in range(self.dmua_depth_init,
                       self.dmua_depth_init + self.dmua_depth):
            psf_map += self.ssp[z]

        psf_mua_map = np.empty_like(psf_map)
        psf_mua_map.fill(self.dmua)
        psf_map = psf_map * psf_mua_map

        save_as = self.dod_dir +\
                  f"PSF(z={self.dmua_depth_init}+{self.dmua_depth - 1})"
        np.savetxt(save_as + ".csv", psf_map, delimiter=',', fmt='%lf')
        self.make_topography(psf_map, save_as + ".png")


def main():
    gates = [15]
    for gate in gates:
        topo = DodTopography(sys.argv[1], gate)
        topo.topography_of_dod()
        topo.topography_of_psf()


if __name__ == '__main__':
    main()
