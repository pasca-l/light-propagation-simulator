import os
import subprocess
import shutil
import argparse


def option_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--save_as', type=str, default='../data/')
    parser.add_argument('-i', '--photons_in', type=int, default=10000)
    parser.add_argument('-n', '--new', action='store_true')

    return parser.parse_args()


def main():
    args = option_parser()

    # prepare folder for simulation
    os.makedirs("../data/tssp_map/", exist_ok=True)
    # run simulation
    subprocess.run('gcc montecarlo.c -lm -std=c99', shell=True)
    subprocess.run(f'./a.out {int(args.new)} {args.photons_in}', shell=True)
    os.remove("./a.out")

    # remove unnecessary files (used to resume simulation)
    os.remove("../data/temporary_note.txt")
    os.remove("../data/data.bin")
    os.remove("../data/pd.bin")
    os.remove("../data/ssp.bin")

    # move generated data into "results" folder
    os.rename("../data/", args.save_as)
    os.makedirs("../results/", exist_ok=True)
    shutil.move(args.save_as, "../results/")


if __name__ == '__main__':
    main()