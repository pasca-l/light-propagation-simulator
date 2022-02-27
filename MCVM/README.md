# MonteCarlo simulation for Voxel Model (MCVM)

## Contents
### Settings
- settings.conf
    - general settings for simulation
- model.conf
    - model settings used in simulation

### Codes
- montecarlo.c
    - main simulation program
- exec_montecarlo.sh
    - runs: mkdir_data.py, montecarlo.c, push_data.py
    - mkdir_data.py generates a temporal directory named "data" to save simulation results
    - push_data.py moves completed simulation result to "results" directory


## How to use
1. change settings in setting.conf, and model.conf to fit assuming condition
2. execute exec_montecarlo.sh giving a directory name to save data into
