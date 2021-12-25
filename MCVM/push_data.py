import os
import sys
import shutil

def main(new_name):
    data_dir = "../data/"
    renamed_dir = f"../{new_name}/"
    os.rename(data_dir, renamed_dir)
    shutil.move(renamed_dir, f"../results/{new_name}/")


if __name__ == '__main__':
    main(sys.argv[1])
