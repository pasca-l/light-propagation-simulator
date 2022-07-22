import sys
import glob
import numpy as np

class Profile:
    def __init__(self, dirname):
        self.work_dir = f"./results/{dirname}/"

        self.total_gate = 16
        self.gate_width = 2
        self.depth_init = 14
        self.depth = 4

        self.tssp_topo_dir = self.work_dir +\
                             f"tssp_topography(gatewidth={self.gate_width})/"

        # direction to take profile
        # x_axis (horizontal) = 0; y_axis (vertical) = 1
        self.profile_axis = 0

    def tssp_profile(self):
        profile = np.arange(-30, 31).reshape(61, 1)

        for gatenum in range(0, self.total_gate, self.gate_width):
            tssp_map = self.tssp_topo_dir +\
                       f"tssp(gate={gatenum}+{self.gate_width - 1}," +\
                       f"z={self.depth_init}+{self.depth}).csv"
            temp = np.loadtxt(tssp_map, delimiter=',')
            if self.profile_axis == 0:
                temp = temp[30, :].reshape(61, 1)
            elif self.profile_axis == 1:
                temp = temp[:, 30].reshape(61, 1)

            profile = np.concatenate([profile, temp], 1)

        save_name = self.tssp_topo_dir + f"profile(axis={self.profile_axis})"
        np.savetxt(save_name + ".csv", profile, delimiter=',', fmt='%lf')

    def folder_profile(self, folder):
        profile_folder = self.work_dir + folder
        files = glob.glob(profile_folder + "*.csv")
        profile = np.arange(-30, 31).reshape(61, 1)

        for file in files:
            temp = np.loadtxt(file, delimiter=',')
            if self.profile_axis == 0:
                temp = temp[30, :].reshape(61, 1)
            elif self.profile_axis == 1:
                temp = temp[:, 30].reshape(61, 1)

            profile = np.concatenate([profile, temp], 1)

        save_name = profile_folder + f"profile(axis={self.profile_axis})"
        np.savetxt(save_name + ".csv", profile, delimiter=',', fmt='%lf')


def main():
    genprofile = Profile(sys.argv[1])
    genprofile.tssp_profile()
    # genprofile.folder_profile("dOD(gate=12)/profile/")


if __name__ == '__main__':
    main()
