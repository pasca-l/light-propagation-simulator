# Light propagation simulator

## Project overview
This project holds codes for simulating light propagation through a turbid medium, assuming a biological tissue that scatters incident photon.

Things that can be done:
- Montecarlo simulation
    - accumulates photon weights that are used to calculate intensity and path lengths
- Topography generation
    - generates numerical map of SSP (spatial sensitivity profile) and dOD (optical density change)
- Gaussian distribution fitting
    - fits topography to gaussian distribution


## Simulators
- MCML (MonteCarlo simulation for MultiLayered model)
- MCVM (MonteCarlo simulation for Voxel Model)
