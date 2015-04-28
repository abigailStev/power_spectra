#!/usr/bin/env python

import argparse
import numpy as np
from scipy import fftpack
from astropy.io import fits  
from datetime import datetime
import os.path
import subprocess
import powerspec as psd  # https://github.com/abigailStev/power_spectra
import tools  # https://github.com/abigailStev/whizzy_scripts

__author__ = "Abigail Stevens"

"""
        multi_powerspec.py

Makes an averaged power spectrum from multiple event-mode data files from RXTE.
Type "multi_powerspec.py -h" to see command line requirements and options.

Abigail Stevens, A.L.Stevens at uva.nl, 2013-2015

"""

################################################################################ 
def dat_output(out_file, data_file_list, meta_dict, total_exposure, mean_rate_total, freq, fracrms_power, \
    fracrms_err):
    """
            dat_output

    Writes power spectrum to a .dat output file.

    """
    print "\nOutput sent to %s" % out_file

    with open(out_file, 'w') as out:
        out.write("#\t\tPower spectrum of multiple data files")
        out.write("\n# Date(YYYY-MM-DD localtime): %s" % str(datetime.now()))
        out.write("\n# Data file list: %s" % data_file_list)
        out.write("\n# Time bin size = %.21f seconds" % meta_dict['dt'])
        out.write("\n# Number of bins per segment = %d" % meta_dict['n_bins'])
        out.write("\n# Number of seconds per segment = %d" % meta_dict['num_seconds'])
        out.write("\n# Total number of segments = %d " % meta_dict['num_seg'])
        out.write("\n# Total exposure time = %d seconds" % total_exposure)
        out.write("\n# Mean count rate = %.8f" % mean_rate_total)
        out.write("\n# Nyquist frequency = %.4f" % meta_dict['nyquist'])
        out.write("\n# ")
        out.write("\n# Column 1: Frequency [Hz]")
        out.write("\n# Column 2: Fractional rms normalized mean power")
        out.write("\n# Column 3: Fractional rms normalized error on the mean "\
                  "power")
        out.write("\n# ")
        for k in range(0, len(fracrms_power)):
            if freq[k] >= 0:
                out.write("\n%.8f\t%.8f\t%.8f" % (freq[k], fracrms_power[k],
                                                  fracrms_err[k]))
                ## End of if-statement
            ## End of for-loop
        ## End of with-block


################################################################################
def fits_output(out_file, data_file_list, meta_dict, total_exposure, \
    mean_rate_total, freq, fracrms_power, fracrms_err):
    """
                fits_output

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
    prihdr.set('SEGMENTS', meta_dict['num_seg'], "Segments in the whole light curve")
    prihdr.set('EXPOSURE', total_exposure, "seconds, of whole light "\
        "curve")
    prihdr.set('DETCHANS', meta_dict['detchans'], "Number of detector energy channels")
    prihdr.set('MEANRATE', mean_rate_total, "counts/second")
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
    if os.path.isfile(out_file):
# 		print "File previously existed. Removing and rewriting."
        subprocess.call(["rm", out_file])

    ## Writing the standard power spectrum to a FITS file
    thdulist = fits.HDUList([prihdu, tbhdu])
    thdulist.writeto(out_file)


################################################################################
def main(infile_list, out_file, num_seconds, dt_mult, test):
    """
            main

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

    meta_dict = {'dt': dt, 't_res': t_res, 'num_seconds': num_seconds, 'df': df, 'nyquist': nyquist_freq, 'n_bins': n_bins, 'detchans': 64}
    print meta_dict['dt']

    print "DT = %.15f seconds" % meta_dict['dt']
    print "N_bins = %d" % meta_dict['n_bins']
    print "Nyquist freq = %f" % meta_dict['nyquist']

    ############################
    ## THE BIG LOOP STARTS HERE
    ############################

    for in_file in data_files:

        power_sum, sum_rate_whole, num_seg = \
            psd.read_and_use_segments(in_file, meta_dict, test)
        
        print "Segments for this file: %d\n" % num_seg
        
        total_power_sum += power_sum
        sum_rate_total += sum_rate_whole
        total_seg += num_seg
#         print "Sum_rate_total / total_seg = ", sum_rate_total / \
#         	float(total_seg)
    
    ## End of for-loop
    print " "
    meta_dict['num_seg'] = total_seg

    total_exposure = meta_dict['num_seg'] * meta_dict['num_seconds']
    print "Total exposure time =", total_exposure

    ## Getting averages. Total means over all segments of all data files.
    power = total_power_sum / meta_dict['num_seg']
    mean_rate_total = sum_rate_total / meta_dict['num_seg']
    print "Total segments =", meta_dict['num_seg']
    print "Mean rate total =", mean_rate_total
    
    ######################################################
    ## Normalizing the power spectrum and computing error
    ######################################################
    total_variance = np.sum(power * df)
    print "Total variance:", total_variance, "(unnorm)"
    rms_total = np.sqrt(total_variance)
    print "Total RMS:", rms_total, "(unnorm)"
    
    freq, power, leahy_power, fracrms_power, fracrms_err = \
        psd.normalize(power, meta_dict, mean_rate_total, True)
    
    #########################
    ## Output, .dat or .fits
    #########################
    
    if out_file[-4:].lower() == "fits":
        fits_output(out_file, infile_list, meta_dict, total_exposure, \
            mean_rate_total, freq, fracrms_power, fracrms_err)
    elif out_file[-3:].lower() == "dat":	
        dat_output(out_file, infile_list, meta_dict, total_exposure, \
            mean_rate_total, freq, fracrms_power, fracrms_err)
    else:
        raise Exception("ERROR: Output file must be of type .dat or .fits.")


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

    parser.add_argument('outfile', help="The full path of the (ASCII/txt) "\
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

    args = parser.parse_args()

    test = False
    if args.test == 1: 
        test = True

    main(args.infile_list, args.outfile, args.num_seconds, args.dt_mult, test)

################################################################################
