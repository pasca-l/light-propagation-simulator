# Light propagation simulator

## Project overview
This project holds codes for simulating light propagation through a turbid medium, assuming a biological tissue that scatters incident photon.

## Content
- Montecarlo simulation : accumulates photon weights that are used to calculate intensity and path lengths
- Topography generation : generates numerical map of SSP (spatial sensitivity profile) and dOD (optical density change)
- Gaussian distribution fitting : fits topography to gaussian distribution

## Simulators
- MCML (MonteCarlo simulation for MultiLayered model)
- MCVM (MonteCarlo simulation for Voxel Model)

## Codes
- topography_ssp.py
    - generates 2D topography image of (time-resolved) SSP using resulting csv file given from simulators
- topography_dod.py
    - generates topography image of dOD based on virtual absorption change area using resulting csv file given from simulators
- get_profile.py
    - cuts out values at specified axis of given image matrix
- gaussian_fitting.py
    - fits gaussian distribution to given image matrix and returns optimized function parameters
- gaussian_3d_plot.py
    - plots gaussian distribution or image in 3D
- generate_topodata.sh
    - result passed through: topography_ssp.py, topography_dod.py, gaussian_fitting.py
- deconvolute_topography.py
    - deconvolutes given dOD by SSP, by simple division in frequency domain
