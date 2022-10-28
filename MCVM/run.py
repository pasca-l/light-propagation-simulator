import os
import subprocess
import shutil
import argparse


def option_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--save_name', type=str, default='../data/')

    return parser.parse_args()


def main():
    args = option_parser()

    # prepare folder for simulation
    os.makedirs("../data/tssp_map/", exist_ok=True)
    # run simulation
    subprocess.run('gcc montecarlo.c -lm -std=c99', shell=True)
    subprocess.run('./a.out', shell=True)
    os.remove("./a.out")

    # remove unnecessary files (used to resume simulation)
    os.remove("../data/temporary_note.txt")
    os.remove("../data/data.bin")
    os.remove("../data/pd.bin")
    os.remove("../data/ssp.bin")

    # move generated data into "results" folder
    os.rename("../data/", args.save_name)
    os.makedirs("../results/", exist_ok=True)
    shutil.move(args.save_name, "../results/")


if __name__ == '__main__':
    main()