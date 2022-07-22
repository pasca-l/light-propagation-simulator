# MonteCarlo simulation for Voxel Model (MCVM)

## Contents
### Settings
- `settings.conf`
  - general settings for simulation
- `model.conf`
  - model settings used in simulation

### Codes
- `montecarlo.c`
  - main simulation program
- `run.py`
  - wrapper for simulation, dumps output in "results" directory giving a designated name.

## Usage
1. Change settings in `setting.conf`, and `model.conf` to fit assuming condition.
2. Execute `run.py` giving a directory name to save data into. (Default data directory name is `data`)
```shell
$ python run.py --save_name SAVE_NAME
```
