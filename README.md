
# The Quadrature Method

This repository contains a potential flow simulation and an implementation of ten dipole localization algorithms for artificial lateral lines. The code was used for our paper: "The Quadrature Method: a Novel Dipole Localisation Algorithm for Artificial Lateral Lines Compared to State of the Art" and includes the entire workflow: starting from the potential flow simulation up to generating figures for the paper.

## Getting started

### Used software
The potential flow and artificial lateral line simulation is designed for MATLAB 2018a. Especially the multi-layer perceptron implementation is quite sensitive to the version of MATLAB, as Mathworks actively developed their machine learning framework in recent versions.

The resulting predictions were visualized using Jupyter Lab with:
 - Python 3.7.8
 - Jupyter lab 2.2.9
 - Matplotlib 3.3.2
 - Seaborn 0.11.1
 - Numpy 1.19.2
 - Pandas 1.1.3

### Scripts

The workflow is based on scripts (`sxx_*.m`) that are to be run in sequence. The MATLAB scripts prepare the simulation, perform hyperparameter sweeps, use all algorithms to compute predictions using the experimental conditions of two analysis methods, and finally export all predictions to two csv files. Then, the Jupyter Lab files are used to create the figures used in the paper (as well as some others). The `try_*.m` scripts show how to use the dipole localization algorithms and may include some quick configuration comparisons.

The raw data is available in the `data` directory as `.mat` files. The `s11_export_data.m` script transforms them to `.csv` and computes the prediction errors. These csv versions are available [on zenodo](https://zenodo.org/record/4743839). The `config` directory contains all source locations, final hyper-parameter values, and hyper-parameter validation sweep outputs. The dipole localization algorithms and potential flow simulation are implemented in the `functions` directory. 

A good place to start exploring the code-base is the [`generate_parameters.m`](./functions/generate_parameters.m) function, it contains and explains all options used in the implementations.


## How to cite

### Software
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4744279.svg)](https://doi.rog/10.5281/zenodo.4744279)

```
Bot, D.M., Wolf, B.J., & van Netten, S.M. (2021). Dipole localisation algorithms for simulated artificial lateral line (Version 1.0) [Software]. Zenodo. http://doi.org/10.5281/zenodo.4744279
```

### Data
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4743839.svg)](https://doi.org/10.5281/zenodo.4743839)

```
Bot, D.M., Wolf, B.J., & van Netten, S.M. (2021). Dipole localisation predictions data set (Version 1.0) [Data set]. Zenodo. http://doi.org/10.5281/zenodo.4743839
```

### Publication
[![DOI:10.3390/s21134558](https://zenodo.org/badge/DOI/10.3390/s21134558.svg)](https://doi.org/10.3390/s21134558)

```
Bot, D.M., Wolf, B.J., van Netten, S.M. (2021). The Quadrature Method: a Novel Dipole Localisation Algorithm for Artificial Lateral Lines Compared to State of the Art. Sensors, 21(13), 4558.  https://doi.org/10.3390/s21134558
```

```
@Article{s21134558,
AUTHOR = {Bot, Daniël M. and Wolf, Ben J. and van Netten, Sietse M.},
TITLE = {The Quadrature Method: A Novel Dipole Localisation Algorithm for Artificial Lateral Lines Compared to State of the Art},
JOURNAL = {Sensors},
VOLUME = {21},
YEAR = {2021},
NUMBER = {13},
ARTICLE-NUMBER = {4558},
URL = {https://www.mdpi.com/1424-8220/21/13/4558},
ISSN = {1424-8220},
DOI = {10.3390/s21134558}
}
```

## License

[MIT license](./LICENSE)
