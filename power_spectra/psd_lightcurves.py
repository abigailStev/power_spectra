"""
Class for PSD Lightcurve object and normalized power spectra.
Typically imported as psd_lc. At top of program, along with other import
statements, write: "import psd_lightcurves.py as psd_lc", then you can call
psd_lc.Lightcurve and psd_lc.NormPSD!
"""
import numpy as np

__author__ = 'Abigail Stevens <A.L.Stevens at uva.nl>'
__year__ = "2015"


class Lightcurve(object):
    def __init__(self, n_bins=8192):
        self.power_array = np.zeros((n_bins, 1), dtype=np.float64)
        self.power = np.zeros(n_bins, dtype=np.float64)
        self.pos_power = np.zeros(n_bins/2+1, dtype=np.float64)
        self.mean_rate_array = 0.0
        self.mean_rate = 0.0
        self.var = 0.0   # variance of the absolute-rms-normalized power
                         # spectrum of the ref band
        self.rms = 0.0   # rms of the absolute-rms-normalized power spectrum
                         #  of the ref band
        self.rms_array = 0.0


class NormPSD(object):
    def __init__(self, n_bins=8192):
        self.noise = 0.0
        self.power = np.zeros(n_bins, dtype=np.float64)
        self.variance = 0.0
        self.rms = 0.0