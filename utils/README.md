# MonteCarlo simulation utilities

## Requirements

## Contents
- [`topography_ssp.py`](#topography_dodpy)

## Usage
### [`topography_ssp.py`](https://github.com/pasca-l/light-propagation-simulator/blob/main/utils/topography_ssp.py)
Generates 2D topography image of (time-resolved) SSP using resulting csv file given from simulators.

1. Modify 

### `topography_dod.py`
Generates topography image of dOD based on virtual absorption change area using resulting csv file given from simulators.


### `get_profile.py`
Cuts out values at specified axis of given image matrix.


### `gaussian_fitting.py`
Fits gaussian distribution to given image matrix and returns optimized function parameters.


### `gaussian_3d_plot.py`
Plots gaussian distribution or image in 3D.


### `generate_topodata.sh`
Result passed through: `topography_ssp.py`, `topography_dod.py`, `gaussian_fitting.py`.


### `deconvolute_topography.py`
Deconvolutes given dOD by SSP, by simple division in frequency domain.