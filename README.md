# This repository contains python-based scripts to reproduce the calculations presented in the paper "Impact of Baseline and Cadence on AGN Variability Metrics: A Systematic Study with ZTF" (D. Martinez Collipal & S. Panda 2026).

The LC class allows the user to define a light curve object by providing arrays of times, fluxes and errors from a time series observation. The available methods for this class includes:

## Data processing methods:
### `binning(bin_size)`: This method divides the light curve into time bins of a user defined size. The representative flux for a bin corresponds to the weighted average of the fluxes within that particular time window.
    
### `sigma_clipping(sigma)`: This method rejects outliers from a light curve by filtering data points that deviate by more than $N\sigma$ from the mean flux. The threshold N can be defined by the user with the "sigma" input parameter.
    
### `decadence(target_cadence)`: This method divides the light curve into bins of the desired cadence and randomly keeps a single data point for each bin, the result is an under-sampled version of the original light curve.

## Variability methods:
### `j_index()`: This method computes the Stetson index J of a light curve, following the definition provided in Ma+2024
### `s()`: This method computes the smoothness parameter s of a light curve. It includes the corrections presented in D. Martinez Collipal & S. Panda 2026

To run the scripts it is necessary to install the package numpy, which can be installed via:
```bash
$ pip install numpy
    
