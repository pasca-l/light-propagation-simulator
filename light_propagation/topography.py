import matplotlib.pyplot as plt
import numpy as np

def topography_of_ssp(total_depth=28):
    ssp = np.loadtxt("./results/ssp.csv", skiprows=1, delimiter=',')
    # reshape to (depth, rows, columns)
    ssp = np.reshape(ssp, (total_depth, 61, 61))

    for z in range(total_depth):
        image = np.zeros_like(ssp[z])
        if ssp[z].max() != 0:
            image = ssp[z] / ssp[z].max()

        fig = plt.figure()
        plt.imshow(image, cmap='hot')
        plt.colorbar()
        plt.clim(0, 1)

        fig.savefig(f"./results/ssp_topography/ssp(z={z}).png")

        plt.clf()
        plt.close()


def add_surrounding_zeros(matrix):
    row, col = map(int, matrix.shape)

    zeros1 = np.zeros((int(row/2), col))
    temp = np.vstack((matrix, zeros1))
    temp = np.vstack((zeros1, temp))

    zeros2 = np.zeros((2*row-1, int(col/2)))
    temp = np.hstack((temp, zeros2))
    output = np.hstack((zeros2, temp))

    return output


def create_dmua_map(ssp_map, dmua=0.002, r=3):
    dmua_map = np.zeros_like(ssp_map)
    ymid, xmid = map(lambda f: int(f/2), dmua_map.shape)

    # make circle at left upper corner
    for y in range(2*r+1):
        for x in range(2*r+1):
            if (x-r)**2 + (y-r)**2 < r*r:
                dmua_map[y][x] = dmua

    # move to center
    dmua_map = np.roll(dmua_map, ymid - r, axis=0)
    dmua_map = np.roll(dmua_map, xmid - r, axis=1)
    # printing area of interest for check
    # print(dmua_map[ymid-r:ymid+r+1,xmid-r:xmid+r+1])

    return dmua_map


def topography_of_dod(total_depth=28, mua_depth_init=0, mua_depth=1):
    inputx, inputy = 61, 61
    ssp = np.loadtxt("./results/ssp.csv", skiprows=1, delimiter=',')
    ssp = np.reshape(ssp, (total_depth, inputx, inputy))

    ssp_with_buffer = np.zeros((2*inputx-1, 2*inputy-1))
    xcenter, ycenter = map(lambda f: int(f/2), ssp_with_buffer.shape)
    for z in range(mua_depth_init, mua_depth):
        ssp_with_buffer += add_surrounding_zeros(ssp[z])

    dmua_map = create_dmua_map(ssp_with_buffer)
    dod_map = np.zeros((inputx, inputy))

    for dety in range(int(inputy/2), int(inputy/2) + ycenter + 1):
        for detx in range(int(inputx/2), int(inputx/2) + xcenter + 1):
            temp = np.roll(dmua_map, ycenter - dety - 1, axis=0)
            relative_dmua_map = np.roll(temp, xcenter - detx - 1, axis=1)

            i = dety - int(inputy/2)
            j = detx - int(inputx/2)
            dod_map[i][j] = np.sum(ssp_with_buffer * relative_dmua_map)
            np.savetxt(f'./results/dmua_map/dmua({i},{j}).csv',
                     relative_dmua_map[29:90,29:90], delimiter=',', fmt='%lf')

    np.savetxt('./dOD.csv', dod_map, delimiter=',', fmt='%lf')

    image = np.zeros_like(dod_map)
    if dod_map.max() != 0:
        image = dod_map / dod_map.max()

    fig = plt.figure()
    plt.imshow(image, cmap='hot')
    plt.colorbar()
    plt.clim(0, 1)

    fig.savefig(f"./results/dOD.png")

    plt.clf()
    plt.close()


def main():
    topography_of_dod()


if __name__ == '__main__':
    main()
