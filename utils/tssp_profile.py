import os
import argparse
import glob
import re
import numpy as np


class Profile:
    def __init__(self, args):
        self.work_dir = f"../results/{args.dirname}/"
        self.paths = glob.glob(self.work_dir + "tssp_topography*/*.csv")

        self.profile_axis = args.axis
        self.profile_position = args.position
        self.profile_dir = f"../results/{args.dirname}/profile/"
        os.makedirs(self.profile_dir, exist_ok=True)

    def tssp_profile(self):
        for file in self.paths:
            tssp_map = np.loadtxt(file, delimiter=',')
            info = re.findall(r'(?<=\().+?(?=\))', file)[-1]

            profile = np.arange(-30, 31).reshape(61, 1)

            if self.profile_axis == 0:
                temp = tssp_map[self.profile_position, :].reshape(61, 1)
            elif self.profile_axis == 1:
                temp = tssp_map[:, self.profile_position].reshape(61, 1)

            profile = np.concatenate([profile, temp], 1)

            save_name = self.profile_dir +\
                        f"profile({info},axis={self.profile_axis}," +\
                        f"position={self.profile_position}).csv"
            np.savetxt(save_name, profile, delimiter=',', fmt='%lf')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--dirname', type=str, default='data')
    parser.add_argument('-a', '--axis', type=int, default=0)
    parser.add_argument('-p', '--position', type=int, default=30)

    genprofile = Profile(parser.parse_args())
    genprofile.tssp_profile()


if __name__ == '__main__':
    main()
