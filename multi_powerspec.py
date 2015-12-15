#!/usr/bin/env python
"""
DEPRECIATED. Absorbed into powerspec.py.

Makes an averaged power spectrum from multiple event-mode data files from RXTE.
Type "multi_powerspec.py -h" to see command line requirements and options.
"""

import argparse
import numpy as np
from astropy.io import fits
from datetime import datetime
import os.path
import subprocess
from distutils.util import strtobool
import powerspec as psd  # https://github.com/abigailStev/power_spectra
import tools  # https://github.com/abigailStev/whizzy_scripts

__author__ = "Abigail Stevens <A.L.Stevens at uva.nl>"
__year__ = "2013-2015"


################################################################################
def fits_output(out_file, data_file_list, meta_dict, freq, fracrms_power, \
        fracrms_err):
    """
    Writes power spectrum to a fits output file.

    """
    print "\nOutput sent to %s" % out_file
    detchans=64

    ## Making header for standard power spectrum
    prihdr = fits.Header()
    prihdr.set('TYPE', "Power spectrum of multiple event lists")
    prihdr.set('DATE', str(datetime.now()), "YYYY-MM-DD localtime")
    prihdr.set('IN_LIST', data_file_list)
    prihdr.set('DT', meta_dict['dt'], "seconds")
    prihdr.set('N_BINS', meta_dict['n_bins'], "Time bins per segment")
    prihdr.set('SECONDS', meta_dict['num_seconds'], "Seconds per segment")
    prihdr.set('SEGMENTS', meta_dict['num_seg'], "Segments in the whole light "\
        "curve")
    prihdr.set('EXPOSURE', meta_dict['exposure'], "seconds of data used")
    prihdr.set('DETCHANS', meta_dict['detchans'], "Number of detector energy "\
        "channels")
    prihdr.set('RMS', str(meta_dict['rms']), "Fractional rms of noise-sub PSD.")
    prihdr.set('MEANRATE', meta_dict['mean_rate'], "counts/second")
    prihdr.set('NYQUIST', meta_dict['nyquist'], "Hz")
    prihdu = fits.PrimaryHDU(header=prihdr)

    ## Making FITS table for standard power spectrum
    col1 = fits.Column(name='FREQUENCY', unit='Hz', format='1D', array=freq)
    col2 = fits.Column(name='POWER', unit='frac rms^2', format='1D', \
        array=fracrms_power)
    col3 = fits.Column(name='ERROR', unit='frac rms^2', format='1D', \
        array=fracrms_err)
    cols = fits.ColDefs([col1, col2, col3])
    tbhdu = fits.BinTableHDU.from_columns(cols)

    ## If the file already exists, remove it (still working on just updating it)
    assert out_file[-4:].lower() == "fits", "ERROR: Standard output file must "\
        "have extension '.fits'."

    ## Writing the standard power spectrum to a FITS file
    thdulist = fits.HDUList([prihdu, tbhdu])
    thdulist.writeto(out_file, clobber=True)


################################################################################
def main(infile_list, out_file, num_seconds, dt_mult, test, adjust):
    """
    Reads in one data file at a time, takes FFT of segments of light curve data,
    computes power of each segment, averages power over all segments of all data
    files, writes resulting normalized power spectrum to a file.

    """
    #####################################################
    ## Idiot checks, to ensure that our assumptions hold
    #####################################################
    
    assert num_seconds > 0, "ERROR: Number of seconds per segment must be a "\
                            "positive integer."
    assert dt_mult > 0, "ERROR: Multiple of dt must be a positive integer."

    ##################################
    ## Getting the list of data files
    ##################################
    
    data_files = [line.strip() for line in open(infile_list)]
    if not data_files:  ## If data_files is an empty list
        raise Exception("ERROR: No files in the eventlist list.")

    # adjust_segs = [494, -187, -217, 150, -305, -58, 420, 123, -691] ## for 5.16269 Hz
    adjust_segments = [932, 216, 184, 570, 93, 346, 860, 533, -324]

    ###################
    ## Initializations
    ###################

    t_res = float(tools.get_key_val(data_files[0], 0, 'TIMEDEL'))
    dt = dt_mult * t_res
    df = 1.0 / float(num_seconds)
    n_bins = num_seconds * int(1./dt)
    nyquist_freq = 1.0 / (2.0 * dt)
    total_power_sum = np.zeros(n_bins)
    sum_rate_total = 0
    total_seg = 0
    total_exposure = 0
    dt_total = np.array([])
    df_total = np.array([])
    ellsee_total = psd.Lightcurve()

    meta_dict = {'dt': dt, 't_res': t_res, 'num_seconds': num_seconds,
                 'df': df, 'nyquist': nyquist_freq, 'n_bins': n_bins,
                 'detchans': 64, 'adjust_seg': 0, 'exposure': 0}

    print "DT = %.15f seconds" % meta_dict['dt']
    print "N_bins = %d" % meta_dict['n_bins']
    print "Nyquist freq = %f" % meta_dict['nyquist']

    ellsee_total.power_array = np.zeros((meta_dict['n_bins'], 1), dtype=np.float64)
    ellsee_total.mean_rate_array = 0

    ############################
    ## THE BIG LOOP STARTS HERE
    ############################
    i = 0
    for in_file in data_files:
        if adjust:
            meta_dict['adjust_seg'] = adjust_segments[i]

        ellsee_whole, power_sum, sum_rate_whole, num_seg, exposure, dt_whole, df_whole = \
                psd.read_and_use_segments(in_file, meta_dict, test)
        
        print "Segments for this file: %d\n" % num_seg
        
        total_power_sum += power_sum
        sum_rate_total += sum_rate_whole
        total_seg += num_seg
        total_exposure += exposure
        dt_total = np.append(dt_total, dt_whole)
        df_total = np.append(df_total, df_whole)

        ellsee_total.power_array = np.hstack((ellsee_total.power_array, \
                    ellsee_whole.power_array))
        ellsee_total.mean_rate_array = np.append(ellsee_total.mean_rate_array, \
                    ellsee_whole.mean_rate_array)

        i += 1
#         print "Sum_rate_total / total_seg = ", sum_rate_total / \
#         	float(total_seg)
    
    ## End of for-loop
    print " "

    ellsee_total.power_array = ellsee_total.power_array[:,1:]
    ellsee_total.mean_rate_array = ellsee_total.mean_rate_array[1:]
    meta_dict['num_seg'] = total_seg

    ############################################################################
    ## Applying the mask from ccf to get rid of segments with negative variance
    ############################################################################

    mask = np.array([strtobool(line.strip()) for line in open("mask.txt")], \
            dtype=np.bool)
    # print mask
    # print type(mask)
    # print np.count_nonzero(mask)
    # print np.shape(ellsee_total.power_array)
    # print np.shape(ellsee_total.mean_rate_array)

    ellsee_total.power_array = ellsee_total.power_array[:,~mask]
    ellsee_total.mean_rate_array = ellsee_total.mean_rate_array[~mask]
    meta_dict['num_seg'] = total_seg - np.count_nonzero(mask)
    # print np.count_nonzero(~mask)
    for element in dt_total[mask]:
        total_exposure -= element * meta_dict['n_bins']
    meta_dict['exposure'] = total_exposure
    meta_dict['dt'] = np.mean(dt_total[~mask])
    meta_dict['df'] = np.mean(df_total[~mask])

    # print meta_dict['dt']
    # print meta_dict['df']

    print "Total exposure time =", meta_dict['exposure']

    ####################
    ## Getting averages
    ####################

    power = psd.seg_average(ellsee_total.power_array)
    mean_rate = psd.seg_average(ellsee_total.mean_rate_array)
    meta_dict['mean_rate'] = mean_rate
    print "Total segments =", meta_dict['num_seg']
    print "Mean rate total =", mean_rate

    ######################################################
    ## Normalizing the power spectrum and computing error
    ######################################################

    total_variance = np.sum(power * meta_dict['df'])
    print "Total variance:", total_variance, "(unnorm)"
    rms_total = np.sqrt(total_variance)
    print "Total RMS:", rms_total, "(unnorm)"
    
    freq, power, leahy_power, fracrms_power, fracrms_err, rms = \
        psd.normalize(power, meta_dict, mean_rate, True)

    meta_dict['rms'] = rms
    meta_dict['nyquist'] = freq[n_bins/2]
    print meta_dict['nyquist']

    ##########
    ## Output
    ##########

    fits_output(out_file, infile_list, meta_dict, freq, fracrms_power, \
            fracrms_err)


################################################################################
if __name__ == "__main__":

    ##############################################
    ## Parsing input arguments and calling 'main'
    ##############################################

    parser = argparse.ArgumentParser(usage="python multi_powerspec.py "\
        "infile_list outfile [-n NUM_SECONDS] [-m DT_MULT] [-t {0,1}]",
        description="Makes a power spectrum from multiple RXTE eventlists.",
        epilog="For optional arguments, default values are given in brackets "\
        "at end of description.")

    parser.add_argument('infile_list', help="The full path of the (ASCII/txt) "\
        "file with a list of the input files. One file per line.")

    parser.add_argument('outfile', help="The full path of the FITS "\
        "file to write the frequency and power to.")

    parser.add_argument('-n', '--num_seconds', type=tools.type_power_of_two,
        default=1, dest='num_seconds', help="Number of seconds in each Fourier"\
        " segment. Must be a power of 2, positive, integer. [1]")

    parser.add_argument('-m', '--dt_mult', type=tools.type_power_of_two,
        default=1, dest='dt_mult', help="Multiple of dt (dt is from data file)"\
        " for timestep between bins. Must be a power of 2, positive, integer. "\
        "[1]")

    parser.add_argument('-t', '--test', type=int, default=0, choices={0,1},
        dest='test', help="Int flag: 0 if computing all segments, 1 if "\
        "computing only one segment for testing. [0]")

    parser.add_argument('-a', '--adjust', default=False, action='store_true',
            dest='adjust', help="If present, artificially adjusts the "\
            "frequency of the QPO by changing the segment length. [False]")

    args = parser.parse_args()

    test = False
    if args.test == 1: 
        test = True

    main(args.infile_list, args.outfile, args.num_seconds, args.dt_mult, test, \
            args.adjust)

################################################################################
