# MonteCarlo simulation utilities

## Requirements
- numpy 1.23.1
- matplotlib 3.5.2
- scipy 1.8.1

## Description and Usage
---
[***`topography_ssp.py`***](https://github.com/pasca-l/light-propagation-simulator/blob/main/utils/topography_ssp.py)

Generates 2D topography image of (time-resolved) SSP using resulting csv file given from simulators.

1. Modify instance variables.
> Output of simulation gives the time-resolved SSP in the shape of ($z \times x \times y$) per gate, topography is given by summing up maps along certain dimensions.
>> Along the gate dimension, which is the time-resolution of measurement, is accumulated for `self.total_gate // self.gate_width` numbers of gate.
>
>> Along the $z$ dimension, which is the depth of the model, is accumulated from `self.depth_init` to `self.depth_init + self.depth`.

2. Run script, giving the name of the data directory in "results". (Default data directory name is `data`)
```
$ python topography_ssp.py --dirname DATA_DIRECTORY_NAME
```

> Output of this script is put into a directory named `tssp_topography(gatewidth=...)` with files named `tssp(gate=...,z=...)`.

---
[***`topography_dod.py`***](https://github.com/pasca-l/light-propagation-simulator/blob/main/utils/topography_dod.py)

Generates topography image of dOD based on virtual absorption change area using resulting csv file given from simulators.

1. Modify instance variables.
> dOD is calculated by convolving SSP with designated d$\mu_a$.
>> Likely with SSP topography, for computation reduction, dOD is accumulated from `self.dmua_depth_init` to `self.dmua_depth_init + self.dmua_depth` for a certain gate designated by a list of numbers given to `self.gate`.
>>
>> Unlike `topography_ssp.py`, this script is not able to sum up across the gate dimension.
>
>> d$\mu_a$ has the following properties:
>> - d$\mu_a$ is set to be in the center position, at the same depths accumulated for dOD.
>> - d$\mu_a$ is also set to be circular with radius varying from `self.dmua_r_min` to `self.dmua_r_min`, which all integer radii will be calculated within the range.
>> - The value set for d$\mu_a$ is given by `self.dmua`.
>
>> When SSP is convolved with d$\mu_a$ of a single point, which is equivalent to the product of SSP map with a constant value of d$\mu_a$, gives the point-spread function (PSF).

> dOD measurement is thought to be collected by camera or by scanning.
>> Camera measurement is given by replacing a $pixel \times pixel$ region by a representative value of that region, where $pixel$ is all integers in the range of `self.pixel_min` and `self.pixel_max`.
>
>> Scan measurement is given by spline interpolation between points with a certain $interval$, where $interval$ is all integers in the range of `self.interval_min` and `self.interval_max`.

2. Run script, giving the name of the data directory in "results". (Default data directory name is `data`)
```
$ python topography_ssp.py --dirname DATA_DIRECTORY_NAME
```

> Output of this script is put into a directory named `dOD(gate=...)`, with files named:
> 1. `dOD(z=...,dmuar=...,pixel=...)` under `camera` directory
> 2. `dOD(z=...,dmuar=...,interval=...)` under `scan` directory
> 3. `PSF(z=...)` directly below

---
[***`get_profile.py`***](https://github.com/pasca-l/light-propagation-simulator/blob/main/utils/get_profile.py)

Slice out values at specified axis of given image matrix.

0. "results" folder should contain SSP topographies, given from `topography_ssp.py`.
1. Modify instance variables.
> Profile is given by designating:
> - `self.profile_axis`, the direction of the slice
> - `self.profile_position`, the position on the axis to slice

2. Run script, giving the name of the data directory in "results". (Default data directory name is `data`)
```
$ python get_profile.py --dirname DATA_DIRECTORY_NAME
```

> The operation in the script is applied to all files under the directory named `tssp_topography(gatewidth=...)`, and the output is put into a directory named `profile` with files named `profile(gate=...,z=...,axis=...,position=...)`.

---
[***`gaussian_fitting.py`***](https://github.com/pasca-l/light-propagation-simulator/blob/main/utils/gaussian_fitting.py)

Fits gaussian distribution to given image matrix and returns optimized function parameters.

1. Run script, giving the name of the data directory in "results". (Default data directory name is `data`)
```
$ python topography_ssp.py --dirname DATA_DIRECTORY_NAME
```

> All files related to topography under the directory named `tssp_topography(gatewidth=...)`, `dOD(gate=...)` are searched, and the output is put into a directory named `fit` with files named `tssp_fit`, `dod_camera_fit`, `dod_scan_fit`.

---
[***`gaussian_3d_plot.py`***](https://github.com/pasca-l/light-propagation-simulator/blob/main/utils/gaussian_3d_plot.py)

Plots gaussian distribution or image in 3D.

1. Run script, giving the file path to show 3d plot.
```
$ python gaussian_3d_plot.py --filepath DATA_PATH
```

---
[***`deconvolute_topography.py`***](https://github.com/pasca-l/light-propagation-simulator/blob/main/utils/deconvolute_topography.py)

Deconvolutes given dOD by PSF, by simple division in frequency domain.

1. Run script, giving:
  - the name of the data directory in "results". (Default data directory name is `data`)
  - the path for dOD map to deconvolve. (required)
  - the path for PSF map as divider. (required)
```
$ python topography_ssp.py --dirname DATA_DIRECTORY_NAME --dod_map_path PATH_TO_DOD_MAP --psf_map_path PATH_TO_PSF_MAP
```