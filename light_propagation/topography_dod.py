import os
import sys
import matplotlib.pyplot as plt
import numpy as np

class DodTopography:
    def __init__(self, dirname):
        self.work_dir = f"./results/{dirname}/"
        self.dmua_map_check = False

        self.inputx = 61
        self.inputy = 61
        self.total_depth = 28
        self.gate = 15
        self.dmua_depth_init = 14
        self.dmua_depth = 4
        self.dmua = 0.002
        self.dmua_r_min = 1
        self.dmua_r_max = 15
        self.pixel_min = 1
        self.pixel_max = 5

        # self.ssp_map = np.loadtxt(self.work_dir + "ssp.csv",
        #                           skiprows=1, delimiter=',')
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

    def make_topography(self, dod_map, save_name):
        image = np.zeros_like(dod_map)
        if dod_map.max() != 1:
            image = dod_map / dod_map.max()

        # save image without axis and legend
        plt.imsave(save_name[:-4] + "_img.png", image, cmap='hot')

        # save image with axis and legend
        fig = plt.figure()
        plt.imshow(image, cmap='hot')
        plt.colorbar()
        plt.clim(0, 1)
        fig.savefig(save_name)

        plt.clf()
        plt.close()

    def topography_of_dod(self):
        ssp_with_buffer = np.zeros((2 * self.inputx - 1, 2 * self.inputy - 1))
        for z in range(self.dmua_depth_init,
                       self.dmua_depth_init + self.dmua_depth):
            ssp_with_buffer += self.add_surrounding_zeros(self.ssp[z])

        for r in range(self.dmua_r_min, self.dmua_r_max + 1):
            dmua_map = self.create_dmua_map(ssp_with_buffer, r)
            dod_map = np.zeros((self.inputx, self.inputy))

            dod_map = self.convolute_maps(dod_map, dmua_map, ssp_with_buffer)

            for pixel in range(self.pixel_min, self.pixel_max + 1):
                dod_topo = self.replace_to_representative(dod_map, pixel)

                save_as = self.dod_dir +\
                          f"dOD(z={self.dmua_depth_init}+" +\
                          f"{self.dmua_depth - 1},dmuar={r},pixel={pixel})"
                np.savetxt(save_as + ".csv", dod_topo,
                           delimiter=',', fmt='%lf')
                self.make_topography(dod_map, save_as + ".png")

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
    topo = DodTopography(sys.argv[1])
    topo.topography_of_dod()
    topo.topography_of_psf()


if __name__ == '__main__':
    main()
