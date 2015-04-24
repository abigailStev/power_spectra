#!/usr/bin/env python

from setuptools import setup

setup(
	name = 'xray_power_spectra', 
    version = '0.1',
    description = "Makes a power spectrum from an event-mode data file from "\
    	"RXTE.", 
    py_modules = ['powerspec', 'rebin_powerspec', 'plot_powerspec', 'fit_qpo', \
    	'multi_powerspec'], 
    requires = ['numpy', 'astropy', 'scipy', 'matplotlib'],
    author = 'Abigail Stevens',
    author_email = 'A.L.Stevens@uva.nl',
    classifiers = ['Programming Language :: Python :: 2.7',  
    	'Intended Audience :: Science/Research', 
    	'Topic :: Scientific/Engineering :: Astronomy']
)
