# MonteCarlo simulation utilities

## Requirements


## Short description and Usage
---
[`topography_ssp.py`](https://github.com/pasca-l/light-propagation-simulator/blob/main/utils/topography_ssp.py)

Generates 2D topography image of (time-resolved) SSP using resulting csv file given from simulators.

1. Modify instance variables.
> Output of simulation gives the time-resolved SSP in the shape of ($z \times x \times y$) per gate, topography is given by summing up along certain dimensions.
>> Along the gate dimension, which is the time-resolution of measurement, is accumulated for `self.total_gate // self.gate_width` numbers of gate.
>
>> Along the $z$ dimension, which is the depth of the model, is accumulated from `self.depth_init` to `self.depth_init + self.depth`.

2. Run script, giving the name of the data directory in "results". (Default data directory name is `data`)
```
python topography_ssp.py --dirname DATA_DIRECTORY_NAME
```

---
`topography_dod.py`

Generates topography image of dOD based on virtual absorption change area using resulting csv file given from simulators.

---
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