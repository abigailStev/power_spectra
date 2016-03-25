#!/usr/bin/env python
"""
Geometrically re-bin the power spectrum in frequency, plot it, and save it (as
an astropy FITS table and as a plot).

"""
from __future__ import print_function
import argparse
import numpy as np
from astropy.io import fits
from astropy.table import Table, Column
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
from matplotlib.ticker import ScalarFormatter
from datetime import datetime
import os.path
import subprocess

__author__ = "Abigail Stevens <A.L.Stevens at uva.nl>"
__year__ = "2013-2016"


################################################################################
def type_positive_float(num):
    """
	Check if an input is a positive float, as an argparse type.

	Parameters
	----------
	num : int, long, float, or double
		The number in question.

	Returns
	-------
	n : float
		The input number, if it's a positive float.

	Raises
	------
	ArgumentTypeError if n isn't a real number or a positive float.

    """
    try:
        n = float(num)
    except ValueError or TypeError:
        message = "%d is not a real number." % n
        raise argparse.ArgumentTypeError(message)

    if n >= 0:
        return n
    else:
        message = "%d is not a positive float." % n
        raise argparse.ArgumentTypeError(message)


################################################################################
def fits_out(out_file, rb_out_file, meta_dict, mean_rate_whole, rb_freq, \
    rb_rms2, rb_err):
    """
    Write a frequency re-binned power spectrum to a FITS file.

    """

    print("Re-binned output file: %s" % rb_out_file)

#     ## Updating above header for re-binned power spectrum (extension 0)
#     prihdr = fits.Header()
#     prihdr.set('TYPE', "Re-binned power spectrum")
#
#
#     ## Making FITS table for re-binned power spectrum (extension 1)
#     col1 = fits.Column(name='FREQUENCY', unit='Hz', format='E', array=rb_freq)
#     col2 = fits.Column(name='POWER', unit='frac rms^2', format='E',
#                        array=rb_rms2)
#     col3 = fits.Column(name='ERROR', unit='frac rms^2', format='E',
#                        array=rb_err)
#     cols = fits.ColDefs([col1, col2, col3])
#     tbhdu = fits.BinTableHDU.from_columns(cols)
#
#     ## If the file already exists, remove it (still working on just updating it)
#     assert rb_out_file[-4:].lower() == "fits", "ERROR: Re-binned output file "\
#         "must have extension '.fits'."
#     if os.path.isfile(rb_out_file):
# # 		print("File previously existed. Removing and rewriting.")
#         subprocess.call(["rm", rb_out_file])
#
#     ## Writing the re-binned power spectrum to a FITS file
#     thdulist = fits.HDUList([prihdu, tbhdu])
#     thdulist.writeto(rb_out_file)

    out_table = Table()
    out_table.add_column(Column(data=rb_freq, name='FREQUENCY', unit='Hz'))
    out_table.add_column(Column(data=rb_rms2, name='POWER'))
    out_table.add_column(Column(data=rb_err, name='ERROR'))

    out_table.meta['TYPE'] = "Re-binned power spectrum"
    out_table.meta['DATE'] = str(datetime.now())
    out_table.meta['UNBINOUT'] =  out_file
    out_table.meta['DT'] = np.mean(meta_dict['dt'])
    out_table.meta['N_BINS'] = meta_dict['n_bins']
    out_table.meta['SEGMENTS'] = meta_dict['n_seg']
    out_table.meta['SEC_SEG'] = meta_dict['n_seconds']
    out_table.meta['EXPOSURE'] = meta_dict['exposure']
    out_table.meta['DETCHANS'] = meta_dict['detchans']
    out_table.meta['MEANRATE'] = mean_rate_whole
    out_table.meta['RMS'] = meta_dict['rms']
    out_table.meta['NYQUIST'] = meta_dict['nyquist']
    out_table.meta['DF'] = np.mean(meta_dict['df'])
    out_table.meta['ADJUST'] = meta_dict['adjust_seg']

    out_table.write(rb_out_file, overwrite=True)


################################################################################
def flx2xsp_out(rb_out, freq_min, freq_max, rb_freq, rb_rms2, rb_err):
    """
    Make the correct output table for FLX2XSP.

    """
    dir = os.path.dirname(rb_out)
    base = os.path.basename(rb_out)[:-5]
    out_file = dir+"/"+base+"_flx2xsp.txt"
    print("Table for FLX2XSP:", out_file)

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
    Re-bin the power spectrum in frequency space by some re-binning constant
    (rebin_const > 1).

    Parameters
    ----------
    freq : np.array of floats
        1-D array of the Fourier frequencies.

    power : np.array of floats
        1-D array of the power at each Fourier frequency, with any/arbitrary
        normalization.

    err_power : np.array of floats
        1-D array of the error on the power at each Fourier frequency, with the
        same normalization as the power.

    rebin_const : float
        The constant by which the data were geometrically re-binned.

    Returns
    -------
    rb_freq : np.array of floats
        1-D array of the re-binned Fourier frequencies.

    rb_power : np.array of floats
        1-D array of the power at the re-binned Fourier frequencies, with the
        same normalization as the input power array.

    rb_err : np.array of floats
        1-D array of the error on the power at the re-binned Fourier
        frequencies, with the same normalization as the input error on power.

    freq_min : np.array of floats
        1-D array of the lower bounds of each re-binned frequency bin.

    freq_max : np.array of floats
        1-D array of the upper bounds of each re-binned frequency bin.

    """

    ## Initialize variables
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

    ## Loop through the length of the array power, new bin by new bin, to
    ## compute the average power and frequency of that new geometric bin.
    ## Equations for frequency, power, and error are from A. Ingram's PhD thesis
    while current_m < len(power):
# 	while current_m < 100: # used for debugging

        ## Determine the range of indices this specific geometric bin covers
        bin_range = np.absolute(current_m - prev_m)
        ## Want mean power of data points contained within one geometric bin
        bin_power = np.mean(power[prev_m:current_m])
        ## Compute error in bin -- equation from Adam Ingram's thesis
        err_bin_power2 = np.sqrt(np.sum(err_power[prev_m:current_m] ** 2)) / \
            float(bin_range)

        ## Compute the mean frequency of a geometric bin
        bin_freq = np.mean(freq[prev_m:current_m])

        ## Append values to arrays
        rb_power = np.append(rb_power, bin_power)
        rb_freq = np.append(rb_freq, bin_freq)
        rb_err = np.append(rb_err, err_bin_power2)
        freq_min = np.append(freq_min, freq[prev_m])
        freq_max = np.append(freq_max, freq[current_m])

        ## Increment for the next iteration of the loop
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

def make_gaussfit(rb_freq):
    p = [5.19554316,3.33440272e-01,5.55485390e-03,-1.04032202,5.39625529e-04]
    ## p[0] = mean value, p[1] = standard deviation, p[2] = scale factor
    exp_numerator = -(rb_freq - p[0])**2
    exp_denominator = 2 * p[1]**2
    G = p[2] * np.exp(exp_numerator / exp_denominator) * rb_freq
    return G

def make_lorfit(rb_freq):
    p = [5.18905778,6.32397183e-01,6.08392227e-03,-1.12134668,2.86458229e-04]
    ## p[0] = centroid frequency, p[1] = fwhm, p[2] = scale factor
    numerator = p[1] / (np.pi * 2.0)
    denominator = (rb_freq - p[0]) ** 2 + (1.0/2.0 * p[1]) ** 2
    L = (numerator / denominator) * p[2] * rb_freq
    return L

################################################################################
def plot_rb(plot_file, rebin_const, prefix, rb_freq, vpv, err_vpv):
    """
    Plot the re-binned power spectrum.

    Parameters
    ----------
    plot_file : str
        Name of the file to save the plot to.

    rebin_const : float
        The constant by which the data were geometrically re-binned.

    prefix : str
        The identifying prefix for the file (data ID or object nickname).

    rb_freq : np.array of floats
        1-D array of the frequencies of the re-binned power spectrum.

    vpv : np.array of floats
        1-D array of the power times the frequency of the re-binned power
        spectrum, i.e., nu * P(nu).

    err_vpv : np.array of floats
        1-D array of the error on vpv.

    """
    print("Re-binned power spectrum: %s" % plot_file)

    font_prop = font_manager.FontProperties(size=24)

    bf_gauss = make_gaussfit(rb_freq)
    bf_lor = make_lorfit(rb_freq)

    fig, ax = plt.subplots(1,1, figsize=(10,7.5), dpi=300)
    # fig, ax = plt.subplots(1,1)

    # ax.plot(rb_freq, vpv, lw=2, c='black')
    ax.errorbar(rb_freq, vpv, yerr=err_vpv, lw=3, ls='-', c='blue',
            elinewidth=2, capsize=2)
    # ax.plot(rb_freq, bf_gauss, lw=2, c='blue', ls='--')
    # ax.plot(rb_freq, bf_lor, lw=1, c='red', ls='-.')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(rb_freq[1],np.max(rb_freq))
    # ax.set_xlim(500, 1500)
    # ax.set_xlim(rb_freq[1], 1e2)
    # ax.set_ylim(4e-4, 1e-1)
    ax.xaxis.set_major_formatter(ScalarFormatter())
# 	ax.set_xlabel(r'$\nu$ [Hz]', fontproperties=font_prop)
# 	ax.set_ylabel(r'$\nu$ $\cdot$ P($\nu$) [Hz rms$^2$]', \
# 		fontproperties=font_prop)
    ax.set_xlabel('Frequency (Hz)', fontproperties=font_prop)
    ax.set_ylabel(r'Power$\times$ frequency (frac. rms$^{2}\times$ Hz)', \
        fontproperties=font_prop)
    ax.tick_params(axis='x', labelsize=24, bottom=True, top=True, \
        labelbottom=True, labeltop=False)
    ax.tick_params(axis='y', labelsize=24, left=True, right=True, \
        labelleft=True, labelright=False)
    # ax.set_title(prefix, fontproperties=font_prop)

    plt.savefig(plot_file)
    # plt.show()
    plt.close()


################################################################################
if __name__ == "__main__":

    #########################
    ## Parse input arguments
    #########################

    parser = argparse.ArgumentParser(usage="rebin_powerspec.py tab_file "\
            "rb_out_file [-o plot_file] [-p prefix] [-c rebin_constant]",
            description=__doc__, epilog="For optional arguments, default "\
            "values are given in brackets at end of description.")

    parser.add_argument('tab_file', help="The table file, in .fits "\
            "format, with frequency in column 1, fractional rms^2 power in "\
            "column 2, and error on power in column 3.")

    parser.add_argument('rb_out_file', help="The FITS file to write the re-"\
            "binned power spectrum to.")

    parser.add_argument('-o', '--outfile', required=False, dest='plot_file',
            default="./psd_rb.png", help="The output plot file name. "\
            "[./psd_rb.png]")

    parser.add_argument('-p', '--prefix', required=False, dest='prefix',
            default="--", help="The identifying prefix for the file (data ID "\
            "or object nickname). [--]")

    parser.add_argument('-c', '--rebin_const', required=False, default=1.01,
            dest='rebin_const', type=type_positive_float, help="The constant "\
            "by which the data will be geometrically re-binned in frequency. "\
            "[1.01]")

    args = parser.parse_args()

    assert args.rebin_const >= 1.0 , "ERROR: Re-binning constant must be >= 1."

    #######################################
    ## Read in power spectrum from a table
    #######################################

    if "_cs.fits" in args.tab_file:
        try:
            in_table = Table.read(args.tab_file)
        except IOError:
            print("\tERROR: File does not exist: %s" % args.tab_file)
            exit()

        freq = in_table['FREQUENCY']
        # power_ci = in_table['POWER_CI']
        rms2 = in_table['POWER_REF']
        error = np.zeros(np.shape(rms2))

        evt_list = in_table.meta['EVTLIST']
        meta_dict = {'dt': in_table.meta['DT'],
                     'n_bins': in_table.meta['N_BINS'],
                     'n_seg': in_table.meta['SEGMENTS'],
                     'exposure': in_table.meta['EXPOSURE'],
                     'detchans': in_table.meta['DETCHANS'],
                     'n_seconds': in_table.meta['SEC_SEG'],
                     'nyquist': in_table.meta['NYQUIST'],
                     'rebin_const' : args.rebin_const,
                     'df': in_table.meta['DF']}
        # rate_ci = np.asarray(in_table.meta['RATE_CI'].replace('[',\
        #         '').replace(']','').split(','), dtype=np.float64)
        mean_rate_whole = in_table.meta['RATE_REF']
        rms2 -= 2. / mean_rate_whole

    else:

        # try:
        #     file_hdu = fits.open(args.tab_file)
        # except IOError:
        #     print("\tERROR: File does not exist: %s" % args.tab_file)
        #     exit()
        #
        # table = file_hdu[1].data
        # freq = table.field('FREQUENCY')  # frequency, in Hz
        # rms2 = table.field('POWER')  # fractional rms^2 power
        # error = table.field('ERROR')  # error on power
        #
        # mean_rate_whole = file_hdu[0].header['MEANRATE']
        # meta_dict = {'dt': file_hdu[0].header['DT'],
        #              'nyquist': file_hdu[0].header['NYQUIST'],
        #              'n_bins': file_hdu[0].header['N_BINS'],
        #              'detchans': file_hdu[0].header['DETCHANS'],
        #              'n_seg' : file_hdu[0].header['SEGMENTS'],
        #              'rebin_const' : args.rebin_const,
        #              'exposure' : file_hdu[0].header['EXPOSURE']}
        # file_hdu.close()

        try:
            in_table = Table.read(args.tab_file)
        except IOError:
            print("\tERROR: File does not exist: %s" % args.tab_file)
            exit()

        freq = in_table['FREQUENCY']
        rms2 = in_table['POWER']
        error = in_table['ERROR']

        meta_dict = {'n_bins': in_table.meta['N_BINS'],
                     'dt': in_table.meta['DT'],
                     'nyquist': in_table.meta['NYQUIST'],
                     'detchans': in_table.meta['DETCHANS'],
                     'n_seg': in_table.meta['SEGMENTS'],
                     'exposure': in_table.meta['EXPOSURE'],
                     'rebin_const': args.rebin_const,
                     'df': in_table.meta['DF'],
                     'adjust_seg': in_table.meta['ADJUST'],
                     'rms': in_table.meta['RMS'],
                     'n_seconds': in_table.meta['SEC_SEG']}

        mean_rate_whole = in_table.meta['MEANRATE']


    ################################################
    ## Re-binning the power spectrum by rebin_const
    ################################################

    rb_freq, rb_rms2, rb_err, freq_min, freq_max = geometric_rebinning(freq,
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
