#!/usr/bin/env python

import argparse
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
from matplotlib.ticker import ScalarFormatter
import os.path
import subprocess
import rebin_powerspec as rb_psd

__author__ = 'Abigail Stevens <A.L.Stevens at uva.nl>'

plot_file = "/Users/abigailstevens/Dropbox/Research/power_spectra/out_ps/GX339-BQPO_150901_t64_64sec_adj_multi_rb.eps"

input_files = ['/Users/abigailstevens/Dropbox/Research/power_spectra/out_ps/GX339-BQPO_150901_t64_64sec_adj_3-5kev_rb.fits',
    '/Users/abigailstevens/Dropbox/Research/power_spectra/out_ps/GX339-BQPO_150901_t64_64sec_adj_5-10kev_rb.fits',
    '/Users/abigailstevens/Dropbox/Research/power_spectra/out_ps/GX339-BQPO_150901_t64_64sec_adj_10-20kev_rb.fits']

try:
    file_hdu = fits.open(input_files[0])
except IOError:
    print "\tERROR: File does not exist: %s" % input_files[0]
    exit()

table = file_hdu[1].data
rb_freq = table.field('FREQUENCY')  # frequency, in Hz
rb_rms2 = table.field('POWER')  # fractional rms^2 power
rb_err = table.field('ERROR')  # error on power
vpv_1 = rb_freq * rb_rms2
err_vpv_1 = rb_freq * rb_err

try:
    file_hdu = fits.open(input_files[1])
except IOError:
    print "\tERROR: File does not exist: %s" % input_files[1]
    exit()

table = file_hdu[1].data
rb_freq = table.field('FREQUENCY')  # frequency, in Hz
rb_rms2 = table.field('POWER')  # fractional rms^2 power
rb_err = table.field('ERROR')  # error on power
vpv_2 = rb_freq * rb_rms2
err_vpv_2 = rb_freq * rb_err

try:
    file_hdu = fits.open(input_files[2])
except IOError:
    print "\tERROR: File does not exist: %s" % input_files[2]
    exit()

table = file_hdu[1].data
rb_freq = table.field('FREQUENCY')  # frequency, in Hz
rb_rms2 = table.field('POWER')  # fractional rms^2 power
rb_err = table.field('ERROR')  # error on power
vpv_3 = rb_freq * rb_rms2
err_vpv_3 = rb_freq * rb_err


print "Power spectra: %s" % plot_file

font_prop = font_manager.FontProperties(size=18)

bf_gauss = rb_psd.make_gaussfit(rb_freq)
bf_lor = rb_psd.make_lorfit(rb_freq)

fig, ax = plt.subplots(1,1, figsize=(10,7.5), tight_layout=True, dpi=300)

ax.errorbar(rb_freq, vpv_3, yerr=err_vpv_3, lw=2, ls='-', c='green',
        elinewidth=2, capsize=2)
ax.errorbar(rb_freq, vpv_2, yerr=err_vpv_2, lw=2, ls='-', c='orange',
        elinewidth=2, capsize=2)
# ax.plot(rb_freq, vpv, lw=2, c='black')
ax.errorbar(rb_freq, vpv_1, yerr=err_vpv_1, lw=2, ls='-', c='blue',
        elinewidth=2, capsize=2)
# ax.plot(rb_freq, bf_gauss, lw=2, c='red', ls='--')
# ax.plot(rb_freq, bf_lor, lw=1, c='red', ls='-.')

ax.set_xscale('log')
ax.set_yscale('log')
# 	ax.set_xlim(rb_freq[1],np.max(rb_freq))
ax.set_xlim(rb_freq[1], 64)
ax.set_ylim(1e-5, 3e-1)
ax.xaxis.set_major_formatter(ScalarFormatter())
# 	ax.set_xlabel(r'$\nu$ [Hz]', fontproperties=font_prop)
# 	ax.set_ylabel(r'$\nu$ $\cdot$ P($\nu$) [Hz rms$^2$]', \
# 		fontproperties=font_prop)
ax.set_xlabel('Frequency (Hz)', fontproperties=font_prop)
ax.set_ylabel(r'Power$\times$ frequency (frac. rms$^{2}\times$ Hz)', \
    fontproperties=font_prop)
ax.tick_params(axis='x', labelsize=18, bottom=True, top=True, \
    labelbottom=True, labeltop=False)
ax.tick_params(axis='y', labelsize=18, left=True, right=True, \
    labelleft=True, labelright=False)
# ax.set_title(prefix, fontproperties=font_prop)

plt.savefig(plot_file)
# plt.show()
plt.close()

subprocess.call(['cp', plot_file, "/Users/abigailstevens/Dropbox/Research/CCF_paper1/"])