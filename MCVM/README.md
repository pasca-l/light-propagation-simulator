# MonteCarlo simulation for Voxel Model (MCVM)

## Contents
### Settings
- `model.conf`
  - model settings used in simulation

### Codes
- `montecarlo.c`
  - main simulation program
- `run.py`
  - wrapper for simulation, dumps output in "results" directory.

## Usage
1. Change settings in `model.conf` to fit assuming condition.
2. Execute `run.py` giving a directory name to save data into. (Default data directory name is `data`)
```shell
$ python run.py [-s SAVE_NAME] [-i INPUT_PHOTON_NUM] [-n]
```
> Options:
> - `-s` or `--save_name` for designating name of the resulting folder.
> - `-i` or `--photons_in` for designating total number of input photons to simulate.
> - `-n` or `--new` to start simulation from scratch.
