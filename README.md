# power_spectra

Creates a power spectrum (psd, power spectral density, periodogram) of X-ray 
timing data. Intended to be used with event-mode RXTE data.
Note that 'tools.py' is in the repo 'whizzy_scripts'.

Disclaimer -- There is no plan to develop and maintain this code. Please look 
into [Stingray](https://StingraySoftware.github.io)!

(I don't think the packaging for this software works.)


## Contents

### power_spectra/fit_qpo.py
Fits a QPO in a re-binned power spectrum. Uses either a Gaussian and a power law
or a Lorentzian and a power law (specified on the command line). 

### power_spectra/plot_multi_psds.py
Hacky script to plot three power spectra together on one plot.

### power_spectra/plot_powerspec.py
Linearly plots a power spectrum and saves the plot.

### power_spectra/psd_lightcurves.py
Classes for light curves and power spectra.

### power_spectra/powerspec.py
Makes a power spectrum for a light curve. Reads in a filtered event list from a
FITS or ASCII/txt/dat file, 'populates' each light curve segment, takes power 
spectrum of each segment, averages those power spectra over all segments of the
light curve, applies fractional rms^2 normalization (able to go into code to do
Leahy or absolute rms^2), computes sample frequency of averaged power spectrum, 
writes to a file.

### power_spectra/rebin_powerspec.py
Re-bins the power spectrum by frequency by a specified re-binning constant. 
Saves the re-binned power spectrum to a file, plots it logarithmically, and 
saves the plot.

### scripts/loop_powerspec.sh
Loops through a list of eventlists, making a power spectrum for each, and makes
a gif of them at the end.

### scripts/run_powerspec.sh
A bash script to run all the code: powerspec.py, rebin_powerspec.py, 
plot_powerspec.py, and fit_qpo.py.

### scripts/TimmerKoenig.ipynb
An iPython notebook that does a Timmer & Koenig simulation of a time series with
a specific power spectrum, makes it into a light curve with or without Poisson
noise, and converts back to a power spectrum. Makes many plots. Used for testing
and learning.

[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)
 
## Authors
* Abigail Stevens (UvA API)

## Copyright

All content © 2014-2017 the Authors, and is distributed under the MIT
Licence. See LICENSE for details.

