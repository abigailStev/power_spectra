#!/usr/bin/env python

import argparse
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
from matplotlib.ticker import ScalarFormatter
import matplotlib.ticker as ticker
from datetime import datetime
import os.path
import subprocess
from tools import type_positive_float

__author__ = "Abigail Stevens"

"""
plot_rb_powerspec.py

Plots a log power spectrum that has been re-binned in frequency, in the
frequency domain.
Type "rebin_powerspec.py -h" to see command line requirements and options.

Abigail Stevens, A.L.Stevens at uva.nl, 2013-2015

the module 'tools' is in my whizzy_scripts git repo.

"""

################################################################################
def fits_out(out_file, rb_out_file, meta_dict, mean_rate_whole, rb_freq, \
    rb_rms2, rb_err):
    """
    Writes a frequency re-binned power spectrum to a FITS file.

    """

    print "Re-binned output file: %s" % rb_out_file

    ## Updating above header for re-binned power spectrum (extension 0)
    prihdr = fits.Header()
    prihdr.set('TYPE', "Re-binned power spectrum")
    prihdr.set('DATE', str(datetime.now()), "YYYY-MM-DD localtime")
    prihdr.set('UNBINOUT', out_file, "Corresponding un-binned output.")
    prihdr.set('REBIN', meta_dict['rebin_const'], "Freqs re-binned by REBIN * prev_bin_size")
    prihdr.set('DT', meta_dict['dt'], "seconds")
    prihdr.set('N_BINS', meta_dict['n_bins'], "time bins per segment")
    prihdr.set('SEGMENTS', meta_dict['num_seg'], "segments in the whole light curve")
    prihdr.set('EXPOSURE', meta_dict['exposure'], "seconds, of light "\
        "curve")
    prihdr.set('MEANRATE', mean_rate_whole, "counts/second")
    prihdr.set('NYQUIST', meta_dict['nyquist'], "Hz")
    prihdu = fits.PrimaryHDU(header=prihdr)

    ## Making FITS table for re-binned power spectrum (extension 1)
    col1 = fits.Column(name='FREQUENCY', unit='Hz', format='E', array=rb_freq)
    col2 = fits.Column(name='POWER', unit='frac rms^2', format='E',
                       array=rb_rms2)
    col3 = fits.Column(name='ERROR', unit='frac rms^2', format='E',
                       array=rb_err)
    cols = fits.ColDefs([col1, col2, col3])
    tbhdu = fits.BinTableHDU.from_columns(cols)

    ## If the file already exists, remove it (still working on just updating it)
    assert rb_out_file[-4:].lower() == "fits", "ERROR: Re-binned output file "\
        "must have extension '.fits'."
    if os.path.isfile(rb_out_file):
# 		print "File previously existed. Removing and rewriting."
        subprocess.call(["rm", rb_out_file])

    ## Writing the re-binned power spectrum to a FITS file
    thdulist = fits.HDUList([prihdu, tbhdu])
    thdulist.writeto(rb_out_file)


################################################################################
def flx2xsp_out(rb_out, freq_min, freq_max, rb_freq, rb_rms2, rb_err):
    """
    Makes the correct output table for FLX2XSP.

    """
    dir = os.path.dirname(rb_out)
    base = os.path.basename(rb_out)[:-5]
    out_file = dir+"/"+base+"_flx2xsp.txt"
    print "Table for FLX2XSP:", out_file

    with open(out_file, 'w') as out:
        for i in range(len(rb_freq)):
            delta_nu = freq_max[i] - freq_min[i]
            pow = rb_rms2[i] * rb_freq[i] * delta_nu
            err = rb_err[i] * rb_freq[i] * delta_nu
            out.write("%f \t%f \t%.6e \t%.6e\n" % (freq_min[i], freq_max[i],
                                                   pow, err))


################################################################################
def geometric_rebinning(freq, power, err_power, rebin_const):
    """
    Re-bins the power spectrum in frequency space by some re-binning constant
    (rebin_const > 1).

    """

    ## Initializing variables
    rb_power = np.asarray([])  # List of re-binned power
    rb_freq = np.asarray([])   # List of re-binned frequencies
    rb_err = np.asarray([])	   # List of error in re-binned power
    real_index = 1.0		   # The unrounded next index in power
    int_index = 1			   # The int of real_index, added to current_m every
                               #  iteration
    current_m = 1			   # Current index in power
    prev_m = 0				   # Previous index m
    bin_power = 0.0			   # The power of the current re-binned bin
    bin_freq = 0.0			   # The frequency of the current re-binned bin
    err_bin_power2 = 0.0	   # The error squared on 'bin_power'
    bin_range = 0.0			   # The range of un-binned bins covered by this
                               #  re-binned bin
    freq_min = np.asarray([])
    freq_max = np.asarray([])

    ## Looping through the length of the array power, new bin by new bin, to
    ## compute the average power and frequency of that new geometric bin.
    ## Equations for frequency, power, and error are from A. Ingram's PhD thesis
    while current_m < len(power):
# 	while current_m < 100: # used for debugging

        ## Determining the range of indices this specific geometric bin covers
        bin_range = np.absolute(current_m - prev_m)
        ## Want mean power of data points contained within one geometric bin
        bin_power = np.mean(power[prev_m:current_m])
        ## Computing error in bin -- equation from Adam Ingram's thesis
        err_bin_power2 = np.sqrt(np.sum(err_power[prev_m:current_m] ** 2))\
            / float(bin_range)

        ## Computing the mean frequency of a geometric bin
        bin_freq = np.mean(freq[prev_m:current_m])

        ## Appending values to arrays
        rb_power = np.append(rb_power, bin_power)
        rb_freq = np.append(rb_freq, bin_freq)
        rb_err = np.append(rb_err, err_bin_power2)
        freq_min = np.append(freq_min, freq[prev_m])
        freq_max = np.append(freq_max, freq[current_m])

        ## Incrementing for the next iteration of the loop
        ## Since the for-loop goes from prev_m to current_m-1 (since that's how
        ## the range function and array slicing works) it's ok that we set
        ## prev_m = current_m here for the next round. This will not cause any
        ## double-counting bins or skipping bins.
        prev_m = current_m
        real_index *= rebin_const
        int_index = int(round(real_index))
        current_m += int_index
        bin_range = None
        bin_freq = None
        bin_power = None
        err_bin_power2 = None

    return rb_freq, rb_power, rb_err, freq_min, freq_max


################################################################################
def plot_rb(plot_file, rebin_const, prefix, rb_freq, vpv, err_vpv):
    """
    Plots the re-binned power spectrum.

    """
    print "Re-binned power spectrum: %s" % plot_file

    font_prop = font_manager.FontProperties(size=18)

    fig, ax = plt.subplots(1,1)
    ax.plot(rb_freq, vpv, lw=2)
# 	ax.errorbar(rb_freq, vpv, yerr=err_vpv, lw=2, c='orange', elinewidth=1, \
# 		capsize=1)
    ax.set_xscale('log')
    ax.set_yscale('log')
# 	ax.set_xlim(rb_freq[1],np.max(rb_freq))
    ax.set_xlim(rb_freq[1], 1e2)
    ax.set_ylim(1e-4, 1)
# 	ax.set_xlabel(r'$\nu$ [Hz]', fontproperties=font_prop)
# 	ax.set_ylabel(r'$\nu$ $\cdot$ P($\nu$) [Hz rms$^2$]', \
# 		fontproperties=font_prop)
    ax.set_xlabel('Frequency (Hz)', fontproperties=font_prop)
    ax.set_ylabel(r'Power$\times$ frequency (frac. rms$^{2}\times$ Hz)', \
        fontproperties=font_prop)
    ax.tick_params(axis='x', labelsize=16, bottom=True, top=True, \
        labelbottom=True, labeltop=False)
    ax.tick_params(axis='y', labelsize=16, left=True, right=True, \
        labelleft=True, labelright=False)
    ax.set_title(prefix, fontproperties=font_prop)

    fig.set_tight_layout(True)
# 	plt.savefig(plot_file, dpi=200)
    plt.savefig(plot_file)
# 	plt.show()
    plt.close()


################################################################################
if __name__ == "__main__":

    ###########################
    ## Parsing input arguments
    ###########################

    parser = argparse.ArgumentParser(usage="rebin_powerspec.py tab_file "\
        "rb_out_file [-o plot_file] [-p prefix] [-c rebin_constant]",
        description="Geometrically re-bins in frequency and plots a power "\
        "spectrum.", epilog="For optional arguments, default values are given "\
        "in brackets at end of description.")

    parser.add_argument('tab_file', help="The table file, in .dat or .fits "\
        "format, with frequency in column 1, fractional rms^2 power in column "\
        "2, and error on power in column 3.")

    parser.add_argument('rb_out_file', help="The FITS file to write the re-"\
        "binned power spectrum to.")

    parser.add_argument('-o', '--outfile', required=False, dest='plot_file', \
        default="./psd_rb.png", help="The output plot file name. "\
        "[./psd_rb.png]")

    parser.add_argument('-p', '--prefix', required=False, dest='prefix', \
        default="--", help="The identifying prefix for the file (proposal ID "\
        "or object nickname). [--]")

    parser.add_argument('-c', '--rebin_const', required=False, default=1.01,
        dest='rebin_const', type=type_positive_float, help="The constant by "\
        "which the data will be geometrically re-binned. [1.01]")

    args = parser.parse_args()

    assert args.rebin_const >= 1.0 , "ERROR: Re-binning constant must be >= 1."

    ##########################################
    ## Reading in power spectrum from a table
    ##########################################

    try:
        file_hdu = fits.open(args.tab_file)
    except IOError:
        print "\tERROR: File does not exist: %s" % args.tab_file
        exit()

    table = file_hdu[1].data
    freq = table.field('FREQUENCY')  # frequency, in Hz
    rms2 = table.field('POWER')  # fractional rms^2 power
    error = table.field('ERROR')  # error on power

    mean_rate_whole = file_hdu[0].header['MEANRATE']
    meta_dict = {'dt': file_hdu[0].header['DT'], \
                 'nyquist': file_hdu[0].header['NYQUIST'], \
                 'n_bins': file_hdu[0].header['N_BINS'], \
                 'detchans': file_hdu[0].header['DETCHANS'], \
                 'num_seg' : file_hdu[0].header['SEGMENTS'], \
                 'rebin_const' : args.rebin_const, \
                 'exposure' : file_hdu[0].header['EXPOSURE']}
    file_hdu.close()

    ################################################
    ## Re-binning the power spectrum by rebin_const
    ################################################

    rb_freq, rb_rms2, rb_err, freq_min, freq_max = geometric_rebinning(freq, \
        rms2, error, meta_dict['rebin_const'])

    ########################################
    ## Want to plot nu * P(nu) in log space
    ########################################

    vpv = rb_freq * rb_rms2
    err_vpv = rb_freq * rb_err

    #############
    ## Plotting!
    #############

    plot_rb(args.plot_file, meta_dict['rebin_const'], args.prefix, rb_freq, \
        vpv, err_vpv)

    ##########################################################
    ## Writing the re-binned power spectrum to an output file
    ##########################################################

    fits_out(args.tab_file, args.rb_out_file, meta_dict, mean_rate_whole, \
        rb_freq, rb_rms2, rb_err)

    flx2xsp_out(args.rb_out_file, freq_min, freq_max, rb_freq, rb_rms2, \
        rb_err)

################################################################################
