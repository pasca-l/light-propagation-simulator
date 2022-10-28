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
2. Execute `run.py` giving options below:
    - `-s` or `--save_name` for designating name of the resulting folder. (Default to `data`)
    - `-i` or `--photons_in` for designating total number of input photons to simulate. (Default to `10000`)
    - `-n` or `--new` to start simulation from scratch.
```shell
$ python run.py [-s SAVE_NAME] [-i INPUT_PHOTON_NUM] [-n]
```
