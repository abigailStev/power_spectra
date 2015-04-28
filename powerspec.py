#!/usr/bin/env python

import argparse
import numpy as np
import scipy.fftpack as fftpack
from astropy.io import fits
from datetime import datetime
import os.path
import subprocess

import tools  # https://github.com/abigailStev/whizzy_scripts

__author__ = "Abigail Stevens"

"""
powerspec.py

Makes a power spectrum averaged over segments from an RXTE event-mode data file.
Type "powerspec.py -h" to see command line requirements and options.

Abigail Stevens, A.L.Stevens at uva.nl, 2013-2015

"""

################################################################################
def dat_out(out_file, in_file, meta_dict, mean_rate_whole, freq, fracrms_power,\
    fracrms_err, leahy_power):
    """
    Writes power spectrum to an ASCII output file.

    Parameters
    ----------
    out_file : string
        Description.

    infile : string
        Description.

    meta_dict : dict
        Description.

    mean_rate_whole : double
        Description.

    freq : np.array of doubles
        Description.

    fracrms_power : np.array of doubles
        Description.

    fracrms_err : np.array of doubles
        Description.

    leahy_power : np.array of doubles
        Description.

    Returns
    -------
    nothing

    """
    print "\nOutput file: %s" % out_file

    with open(out_file, 'w') as out:
        out.write("#\t\tPower spectrum")
        out.write("\n# Date(YYYY-MM-DD localtime): %s" % str(datetime.now()))
        out.write("\n# Data: %s" % in_file)
        out.write("\n# Time bin size = %.21f seconds" % meta_dict['dt'])
        out.write("\n# Number of bins per segment = %d" % meta_dict['n_bins'])
        out.write("\n# Number of segments per light curve = %d" % \
                  meta_dict['num_seg'])
        out.write("\n# Exposure time = %d seconds" % \
                  (meta_dict['num_seg'] * meta_dict['num_seconds']))
        out.write("\n# Mean count rate = %.8f" % mean_rate_whole)
        out.write("\n# Nyquist frequency = %.4f" % meta_dict['nyquist'])
        out.write("\n# ")
        out.write("\n# Column 1: Frequency [Hz]")
        out.write("\n# Column 2: Fractional rms^2 normalized mean power")
        out.write("\n# Column 3: Fractional rms^2 normalized error on the " \
                  "mean power")
        out.write("\n# Column 4: Leahy-normalized mean power")
        out.write("\n# ")
        for k in xrange(len(fracrms_power)):
            if freq[k] >= 0:
                out.write("\n{0:.8f}\t{1:.8f}\t{2:.8f}\t{3:.8f}".format(freq[k],
                    fracrms_power[k], fracrms_err[k], leahy_power[k]))



################################################################################
def fits_out(out_file, in_file, meta_dict, mean_rate_whole, freq, \
    fracrms_power, fracrms_err, leahy_power):
    """
    Writes power spectrum to a FITS output file.

    Parameters
    ----------
    out_file : string
        Description.

    infile : string
        Description.

    meta_dict : dict
        Description.

    mean_rate_whole : double
        Description.

    freq : np.array of doubles
        Description.

    fracrms_power : np.array of doubles
        Description.

    fracrms_err : np.array of doubles
        Description.

    leahy_power : np.array of doubles
        Description.

    Returns
    -------
    nothing
    """
    print "\nOutput file: %s" % out_file
    meta_dict['detchans']=64

    ## Making header for standard power spectrum
    prihdr = fits.Header()
    prihdr.set('TYPE', "Power spectrum")
    prihdr.set('DATE', str(datetime.now()), "YYYY-MM-DD localtime")
    prihdr.set('EVTLIST', in_file)
    prihdr.set('DT', meta_dict['dt'], "seconds")
    prihdr.set('N_BINS', meta_dict['n_bins'], "time bins per segment")
    prihdr.set('SEGMENTS', meta_dict['num_seg'], "segments in the whole light "\
        "curve")
    prihdr.set('EXPOSURE', meta_dict['num_seg'] * meta_dict['num_seconds'], \
        "seconds, of light curve")
    prihdr.set('DETCHANS', meta_dict['detchans'], "Number of detector energy "\
        "channels")
    prihdr.set('MEANRATE', mean_rate_whole, "counts/second")
    prihdr.set('NYQUIST', meta_dict['nyquist'], "Hz")
    prihdu = fits.PrimaryHDU(header=prihdr)

    ## Making FITS table for standard power spectrum
    col1 = fits.Column(name='FREQUENCY', unit='Hz', format='E', array=freq)
    col2 = fits.Column(name='POWER', unit='frac rms^2', format='E', \
        array=fracrms_power)
    col3 = fits.Column(name='ERROR', unit='frac rms^2', format='E', \
        array=fracrms_err)
    col4 = fits.Column(name='LEAHY', format='E', array=leahy_power)
    cols = fits.ColDefs([col1, col2, col3, col4])
    tbhdu = fits.BinTableHDU.from_columns(cols)

    ## If the file already exists, remove it (still working on just updating it)
    assert out_file[-4:].lower() == "fits", "ERROR: Standard output file must "\
        "have extension '.fits'."
    if os.path.isfile(out_file):
        subprocess.call(["rm", out_file])

    ## Writing the standard power spectrum to a FITS file
    thdulist = fits.HDUList([prihdu, tbhdu])
    thdulist.writeto(out_file)


################################################################################
def normalize(power, meta_dict, mean_rate, noisy):
    """
    Generates the Fourier frequencies, removes negative frequencies, normalizes
    the power by Leahy and fractional rms^2 normalizations, and computes the
    error on the fractional rms^2 power.

    Parameters
    ----------
    power : np.array of doubles
        Description.

    meta_dict : dict
        Description.

    mean_rate : double
        Description.

    noisy : boolean
        Description

    Returns
    -------
    freq : np.array of doubles
        Description.

    power : np.array of doubles
        Description.

    leahy_power : np.array of doubles
        Description.

    fracrms_power : np.array of doubles
        Description.

    fracrms_err : np.array of doubles
        Description.

    """

    ## Computing the FFT sample frequencies (in Hz)
    freq = fftpack.fftfreq(meta_dict['n_bins'], d=meta_dict['dt'])

    ## Ensuring that we're only using and saving the positive frequency
    ## values (and associated power values)
    nyq_index = meta_dict['n_bins']/2
    freq = np.abs(freq[0:nyq_index + 1])  ## because it slices at end-1, and we
        ## want to include 'nyq_index'; abs is because the nyquist freq is both
        ## pos and neg, and we want it pos here.
    power = power[0:nyq_index + 1]

    ## Computing the error on the mean power
    err_power = power / np.sqrt(float(meta_dict['num_seg']) * \
                float(meta_dict['n_bins']))

    ## Absolute rms^2 normalization
    absrms_power = 2.0 * power * meta_dict['dt'] / float(meta_dict['n_bins'])
    absrms_err = 2.0 * err_power * meta_dict['dt'] / float(meta_dict['n_bins'])

    ## Leahy normalization
    leahy_power = absrms_power / mean_rate
    leahy_err = absrms_err / mean_rate
    print "Mean value of Leahy power =", np.mean(leahy_power)  ## ~2

    ## Fractional rms^2 normalization
    fracrms_power = absrms_power / (mean_rate ** 2)
    fracrms_err = absrms_err / (mean_rate ** 2)

    ## Compute the Poisson noise level and subtract it off the power
    if noisy:
        fracrms_noise = 2.0 / mean_rate
        absrms_noise = 2.0 * mean_rate
        leahy_noise = 2.0

        if np.max(freq) > 100:
            print "Mean above 100Hz:", \
                np.mean(absrms_power[np.where(freq >= 100)])
            print "Absrms noise:", absrms_noise
        fracrms_power -= fracrms_noise
        absrms_power -= absrms_noise

    ## Find the pulsation signal frequency -- assumes that the signal dominates
    ## the power spectrum.
    signal_freq = freq[np.argmax(power*freq)]
    print "Frequency of maximum power:", signal_freq

    ## Computing and printing the rms of the signal
    # min_freq_mask = freq < signal_freq  ## we want the last 'True' element
    # max_freq_mask = freq > signal_freq  ## we want the first 'True' element
    # j_min = list(min_freq_mask).index(False)
    # j_max = list(max_freq_mask).index(True)
    # ## Extracting only the signal frequencies of the power
    # signal_pow = np.float64(fracrms_power[j_min:j_max])
    #
    # ## Computing variance and rms of the signal
    # signal_variance = np.sum(signal_pow * meta_dict['df'])
    # print "Signal variance:", signal_variance, "(frac rms2)"
    # rms_signal = np.sqrt(signal_variance)  ## should be a few % in frac rms units
    # print "Signal RMS:", rms_signal, "(frac rms2)"

    ## Variance and rms of the whole averaged power spectrum
    total_variance = np.sum(fracrms_power * meta_dict['df'])
    print "Total variance:", total_variance, "(frac rms2)"
    rms_total = np.sqrt(total_variance)
    print "Total RMS:", rms_total, "(frac rms2)"

# 	print "Mean power:", np.mean(fracrms_power[np.where(freq >= 100)])

    return freq, power, leahy_power, fracrms_power, fracrms_err


################################################################################
def make_ps(rate):
    """
            make_ps

    Computes the mean count rate, the FFT of the count rate minus the mean, and
    the power spectrum of this segment of data.

    """
    ## Computing the mean count rate of the segment
    mean_rate = np.mean(rate)

    ## Subtracting the mean rate off each value of 'rate'
    ## This eliminates the spike at 0 Hz
    rate_sub_mean = rate - mean_rate

    ## Taking the 1-dimensional FFT of the time-domain photon count rate
    ## Using the SciPy FFT algorithm, as it's faster than NumPy for large lists
    fft_data = fftpack.fft(rate_sub_mean)

    ## Computing the power
    power_segment = np.absolute(fft_data) ** 2

    return power_segment, mean_rate


################################################################################
def extracted_in(in_file, meta_dict, print_iterator, test):
    """
            extracted_powerspec

    Opens the FITS file light curve (as created in seextrct), reads the count
    rate for a segment, calls 'make_ps' to create a power spectrum, adds power
    spectra over all segments.

    """

    ## Open the fits file and load the data
    try:
        fits_hdu = fits.open(in_file)
    except IOError:
        print "\tERROR: File does not exist: %s" % in_file
        exit()

    data = fits_hdu[1].data
    fits_hdu.close()

    ## Initializations
    sum_rate_whole = 0
    power_sum = np.zeros(meta_dict['n_bins'], dtype=np.float64)
    num_seg = 0
    i = 0  # start of bin index to make segment of data for inner for-loop
    j = meta_dict['n_bins']  # end of bin index to make segment of data for inner for-loop

    print data[1].field(0) - data[0].field(0)

    assert meta_dict['dt'] == (data[1].field(0) - data[0].field(0)), "ERROR: "\
        "dt must be the same resolution as the extracted FITS data."

    ## Loop through segments of the data
    while j <= len(data.field(1)):  ## while we haven't reached the end of the
                                    ## file

        num_seg += 1

        ## Extracts the second column of 'data' and assigns it to 'rate'.
        rate = data[i:j].field(1)

        power_segment, mean_rate_segment = make_ps(rate)
        power_sum += power_segment
        sum_rate_whole += mean_rate_segment

        if num_seg % print_iterator == 0:
            print "\t", num_seg

        if (test == True) and (num_seg == 1):  # For testing
            break

        ## Clear loop variables for the next round
        rate = None
        power_segment = None
        mean_rate_segment = None

        ## Incrementing the counters and indices
        i = j
        j += meta_dict['n_bins']
        ## Since the for-loop goes from i to j-1 (since that's how the range
        ## function works) it's ok that we set i=j here for the next round.
        ## This will not cause double-counting rows or skipping rows.

    return power_sum, sum_rate_whole, num_seg


################################################################################
def fits_in(in_file, meta_dict, print_iterator, test):
    """
            fits_in

    Opens the .fits GTI'd event list, reads the count rate for a segment,
    populates the light curve, calls 'make_ps' to create a power spectrum, adds
    power spectra over all segments.

    """

    try:
        fits_hdu = fits.open(in_file)
    except IOError:
        print "\tERROR: File does not exist: %s" % in_file
        exit()

    header = fits_hdu[0].header	 ## Header is in ext 0
    data = fits_hdu[1].data  ## Data is in ext 1
    fits_hdu.close()

    sum_rate_whole = 0
    power_sum = np.zeros(meta_dict['n_bins'], dtype=np.float64)
    num_seg = 0
    lightcurve = np.asarray([])

    start_time = data.field('TIME')[0]
    final_time = data.field('TIME')[-1]
    end_time = start_time + meta_dict['num_seconds']

    ## Filter data based on pcu
# 	PCU2_mask = data.field('PCUID') == 2
# 	data = data[PCU2_mask]

    ## Filter data based on energy channel
# 	print np.shape(data)
# 	lower_bound = data.field('CHANNEL') >= 0
# 	data = data[lower_bound]
# 	upper_bound = data.field('CHANNEL') <= 27
# 	data = data[upper_bound]
# 	print np.shape(data)

    all_time = np.asarray(data.field('TIME'), dtype=np.float64)
# 	all_energy = np.asarray(data.field('CHANNEL'), dtype=np.float64)

    ################################
    ## Looping through the segments
    ################################

    while end_time <= final_time:

        time = all_time[np.where(all_time < end_time)]
# 		energy = all_energy[np.where(all_time < end_time)]

        for_next_iteration = np.where(all_time >= end_time)
        all_time = all_time[for_next_iteration]
# 		all_energy = all_energy[for_next_iteration]

        if len(time) > 0:
            num_seg += 1
            rate_1d = tools.make_1Dlightcurve(time, meta_dict['n_bins'], \
                meta_dict['dt'], start_time)
            lightcurve = np.concatenate((lightcurve, rate_1d))

            power_segment, mean_rate_segment = make_ps(rate_1d)
            assert int(len(power_segment)) == meta_dict['n_bins'], "ERROR: "\
                "Something went wrong in make_ps. Length of power spectrum "\
                "segment  != n_bins."
            power_sum += power_segment
            sum_rate_whole += mean_rate_segment

            ## Printing out which segment we're on every x segments
            if num_seg % print_iterator == 0:
                print "\t", num_seg

            ## Clearing variables from memory
            time = None
# 			energy = None
            power_segment = None
            mean_rate_segment = None
            rate_1d = None

            if test and (num_seg == 1):  # Testing
                np.savetxt('tmp_lightcurve.dat', lightcurve, fmt='%d')
                break
            start_time += meta_dict['num_seconds']
            end_time += meta_dict['num_seconds']

        elif len(time) == 0:
            print "No counts in this segment."
            start_time = all_time[0]
            end_time = start_time + meta_dict['num_seconds']
        ## End of 'if there are counts in this segment'

    return power_sum, sum_rate_whole, num_seg


################################################################################
def dat_in(in_file, meta_dict, print_iterator, test):
    """
            dat_in

    Opens the .dat GTI'd event list, reads the count rate for a segment,
    populates the light curve, calls 'make_ps' to create a power spectrum,
    adds power spectra over all segments.

    Assumes that all event times have been corrected with TIMEZERO, and
    GTI-filtered.

    """

    ## Declaring clean variables to append to for every loop iteration.
    time = np.asarray([])
    energy = np.asarray([])
    sum_rate_whole = 0
    power_sum = np.zeros(meta_dict['n_bins'], dtype=np.float64)
    num_seg = 0
    lightcurve = np.asarray([])

    ## Reading only the first line of data to get the start time of the file
    start_time = -99
    with open(in_file, 'r') as fo:
        for line in fo:
            if line[0].strip() != "#":
                line = line.strip().split()
                start_time = np.float64(line[0])
                break
    if start_time is -99:
        raise Exception("ERROR: Start time of data was not read in. Exiting.")

    end_time = start_time + meta_dict['num_seconds']
    assert end_time > start_time, "ERROR: End time must come after start time "\
        "of the segment."

    ## Open the file and read in the events into segments
    with open(in_file, 'r') as f:
        for line, next_line in tools.pairwise(f):
            if line[0].strip() != "#" and \
                float(line.strip().split()[0]) >= start_time:

                ## If the line is not a comment
                line = line.strip().split()
                next_line = next_line.strip().split()
                current_time = np.float64(line[0])
                current_chan = np.int8(line[1])
                current_pcu = np.int8(line[2])

                if current_pcu == 2:  ## Only want PCU2 here
                    time = np.append(time, current_time)
                    energy = np.append(energy, current_chan)

                next_time = float(next_line[0])
                next_end_time = end_time + meta_dict['num_seconds']

                if next_time >= end_time:  ## Triggered at the end of a segment

                    if len(time) > 0:
                        num_seg += 1
                        rate_2d, rate_1d = tools.make_lightcurve(time, energy,
                            meta_dict['n_bins'], meta_dict['dt'], start_time)
                        lightcurve = np.concatenate((lightcurve, rate_1d))

                        power_segment, mean_rate_segment = make_ps(rate_1d)
                        assert int(len(power_segment)) == meta_dict['n_bins']

                        power_sum += power_segment
                        sum_rate_whole += mean_rate_segment
                        ## Printing out which segment we're on every x segments
                        if num_seg % print_iterator == 0:
                            print "\t", num_seg

                        ## Clearing variables from memory
                        power_segment = None
                        mean_rate_segment = None
                        rate_2d = None
                        rate_1d = None
                        time = np.asarray([])
                        energy = np.asarray([])

                        if test and (num_seg == 1):  ## Testing
                            np.savetxt('lightcurve.dat', lightcurve, fmt='%d')
                            break

                        start_time += meta_dict['num_seconds']
                        end_time += meta_dict['num_seconds']

                    ## This next bit helps it handle gappy data; keep in mind
                    ## that end_time has already been incremented here
                    elif len(time) == 0:
                        start_time = next_time
                        end_time = start_time + meta_dict['num_seconds']

    return power_sum, sum_rate_whole, num_seg


################################################################################
def read_and_use_segments(in_file, meta_dict, test):
    """
            read_and_use_segments

    Opens the file, reads in the count rate, calls 'make_ps' to create a
    power spectrum. Separated from main body like this so I can easily call it
    in multi_powerspec.py. Split into 'fits_in', 'dat_in', and
    'extracted_in' for easier readability.

    """
    assert tools.power_of_two(meta_dict['n_bins']) , "ERROR: n_bins must be a "\
        "power of 2."

    print "Input file: %s" % in_file

    if meta_dict['n_bins'] == 32768:
        print_iterator = int(10)
    elif meta_dict['n_bins'] < 32768:
        print_iterator = int(20)
    elif meta_dict['n_bins'] >= 2097152:
        print_iterator = int(1)
    elif meta_dict['n_bins'] >= 1048576:
        print_iterator = int(2)
    else:
        print_iterator = int(5)

    ## Looping through length of data file, segment by segment, to compute
    ## power for each data point in the segment
    print "Segments computed:"

    ## Data from dat and fits need to be populated as a light curve before a
    ## power spectrum can be taken, whereas data from lc was made in seextrct
    ## and so it's already populated as a light curve
    if (in_file[-5:].lower() == ".fits"):
        power_sum, sum_rate_whole, num_seg = fits_in(in_file, meta_dict, \
            print_iterator, test)
    elif (in_file[-4:].lower() == ".dat"):
        power_sum, sum_rate_whole, num_seg = dat_in(in_file, meta_dict, \
            print_iterator, test)
    elif (in_file[-3:].lower() == ".lc"):
        power_sum, sum_rate_whole, num_seg = extracted_in(in_file, meta_dict,\
            print_iterator, test)
    else:
        raise Exception("ERROR: Input file type not recognized. Must be .dat, "\
            ".fits, or .lc.")

    return power_sum, sum_rate_whole, num_seg


################################################################################
def main(in_file, out_file, num_seconds, dt_mult, test):
    """
            main

    Reads in an eventlist, takes FFT of segments of the light curve, computes
    power of each segment, averages power over all segments of the light curve,
    writes resulting normalized power spectrum to a file.

    """

    ################################################
    ## Idiot checks, to ensure our assumptions hold
    ################################################

    assert num_seconds > 0, "ERROR: num_seconds must be a positive integer."
    assert dt_mult >= 1, "ERROR: dt_mult must be a positive integer."
    
    ###################
    ## Initializations
    ###################
    
    t_res = float(tools.get_key_val(in_file, 0, 'TIMEDEL'))
    dt = dt_mult * t_res
    n_bins = num_seconds * int(1.0 / dt)
    nyquist_freq = 1.0 / (2.0 * dt)
    df = 1.0 / float(num_seconds)


    meta_dict = {'dt': dt, 't_res': t_res, 'num_seconds': num_seconds, \
                 'df': df, 'nyquist': nyquist_freq, 'n_bins': n_bins, \
                 'detchans': 64}

    print "\nDT = %f seconds" % meta_dict['dt']
    print "N_bins = %d" % meta_dict['n_bins']
    print "Nyquist freq =", meta_dict['nyquist']

    ############################################################################
    ## Read in segments of a light curve or event list and make a power spectrum
    ############################################################################

    power_sum, sum_rate_whole, num_seg = read_and_use_segments(in_file, \
        meta_dict, test)
    print " "

    meta_dict['num_seg'] = num_seg

    assert isinstance(meta_dict, object)
    print "\tTotal number of segments =", meta_dict['num_seg']

    #########################################################################
    ## Dividing sums by the number of segments to get an arithmetic average.
    #########################################################################

    power = power_sum / float(meta_dict['num_seg'])
    assert int(len(power)) == meta_dict['n_bins'], "ERROR: Power should have "\
        "length n_bins."
    mean_rate_whole = sum_rate_whole / float(meta_dict['num_seg'])
    print "Mean count rate over whole lightcurve =", mean_rate_whole

    ################################
    ## Normalize the power spectrum
    ################################

    freq, power, leahy_power, fracrms_power, fracrms_err = normalize(power, \
        meta_dict, mean_rate_whole, True)
        ## True = noisy light curve (a.k.a. using real data, or simulated light
        ## curve is given a non-zero noise level)

    ###################################
    ## Output, based on file extension
    ###################################

    if out_file[-4:].lower() == "fits":
        fits_out(out_file, in_file, meta_dict, mean_rate_whole, freq, \
            fracrms_power, fracrms_err, leahy_power)
    elif out_file[-3:].lower() == "dat":
        dat_out(out_file, in_file, meta_dict, mean_rate_whole, freq, \
            fracrms_power, fracrms_err, leahy_power)
    else:
        raise Exception("ERROR: Output file must be type .dat or .fits.")


################################################################################
if __name__ == "__main__":

    ##############################################
    ## Parsing input arguments and calling 'main'
    ##############################################

    parser = argparse.ArgumentParser(usage="powerspec.py infile outfile [-n "\
        "NUM_SECONDS] [-m DT_MULT] [-t {0,1}]", description="Makes a power "\
        "spectrum from an event-mode data file from RXTE.", epilog="For "\
        "optional arguments, default values are given in brackets at end of "\
        "description.")

    parser.add_argument('infile', help="The full path of the input file with "\
        "RXTE event-mode data, with time in column 1 and rate in column 2. "\
        "FITS format must have extension .lc or .fits, otherwise assumes .dat "\
        "(ASCII/txt) format.")

    parser.add_argument('outfile', help="The full path of the .fits or .dat "\
        "file to write the frequency and power to.")

    parser.add_argument('-n', '--num_seconds', type=tools.type_power_of_two, \
        default=1, dest='num_seconds', help="Number of seconds in each Fourier"\
        " segment. Must be a power of 2, positive, integer. [1]")

    parser.add_argument('-m', '--dt_mult', type=tools.type_power_of_two, \
        default=1, dest='dt_mult', help="Multiple of dt (dt is from data file)"\
        " for timestep between bins. Must be a power of 2, positive, integer. "\
        "[1]")

# 	parser.add_argument('--test', action='store_true', dest='test', help="If "\
#       "present, only does a short test run.")

    parser.add_argument('-t', '--test', type=int, default=0, choices={0,1}, \
        dest='test', help="Int flag: 0 if computing all segments, 1 if only "\
        "computing one segment for testing. [0]")

    args = parser.parse_args()

    test = False
    if args.test == 1:
        test = True

    main(args.infile, args.outfile, args.num_seconds, args.dt_mult, test)

################################################################################
