#!/usr/bin/env python
"""
Makes a power spectrum averaged over segments from an RXTE event-mode data file.

Type "powerspec.py -h" to see command line requirements and options.

Be sure that 'tools.py' (from https://github.com/abigailStev/whizzy_scripts) is
downloaded, and it's directory is in your PYTHONPATH bash environment variable.

"""

import argparse
import numpy as np
import scipy.fftpack as fftpack
from astropy.io import fits
from datetime import datetime
import os.path
import subprocess
import psd_lightcurves as psd_lc
import tools  ## https://github.com/abigailStev/whizzy_scripts

__author__ = "Abigail Stevens <A.L.Stevens at uva.nl>"
__year__ = "2013-2015"


################################################################################
def fits_out(out_file, in_file, meta_dict, freq, fracrms_power, fracrms_err,
        leahy_power, file_desc):
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

    ## Making header for standard power spectrum
    prihdr = fits.Header()
    prihdr.set('TYPE', file_desc)
    prihdr.set('DATE', str(datetime.now()), "YYYY-MM-DD localtime")
    prihdr.set('EVTLIST', in_file)
    prihdr.set('DT', meta_dict['dt'], "seconds")
    prihdr.set('N_BINS', meta_dict['n_bins'], "time bins per segment")
    prihdr.set('SEGMENTS', meta_dict['n_seg'], "segments in the whole light "\
        "curve")
    prihdr.set('EXPOSURE', meta_dict['exposure'], "seconds of data used")
    prihdr.set('DETCHANS', meta_dict['detchans'], "Number of detector energy "\
        "channels")
    prihdr.set('RMS', meta_dict['rms'], "Fractional rms of noise-sub PSD.")
    prihdr.set('MEANRATE', meta_dict['mean_rate'], "cts/s")
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

    ## Check that the output file name has FITS file extension
    assert out_file[-4:].lower() == "fits", "ERROR: Standard output file must "\
        "have extension '.fits'."

    ## Write the power spectrum to a FITS file
    thdulist = fits.HDUList([prihdu, tbhdu])
    thdulist.writeto(out_file, clobber=True)


################################################################################
def raw_to_absrms(power, mean_rate, n_bins, dt, noisy):
    """
    Normalizes the power spectrum to absolute rms^2 normalization.

    Parameters
    ----------
    power : np.array of floats
        The raw power at each Fourier frequency.

    mean_rate : float
        The mean count rate for the light curve, in cts/s.

    n_bins : int
        Number of bins per segment of light curve.

    dt : float
        Timestep between bins in n_bins, in seconds.

    noisy : boolean
        True if there is Poisson noise in the power spectrum (i.e., from real
        data), False if there is no noise in the power spectrum (i.e.,
        simulations without Poisson noise).

    Returns
    -------
    np.array of floats
        The noise-subtracted power spectrum in absolute rms^2 units.

    """
    if noisy:
        noise = 2.0 * mean_rate
    else:
        noise = 0.0

    return power * (2.0 * dt / float(n_bins)) - noise


################################################################################
def raw_to_fracrms(power, mean_rate, n_bins, dt, noisy):
    """
    Normalizes the power spectrum to fractional rms^2 normalization.

    Parameters
    ----------
    power : np.array of floats
        The raw power at each Fourier frequency.

    mean_rate : float
        The mean count rate for the light curve, in cts/s.

    n_bins : int
        Number of bins per segment of light curve.

    dt : float
        Timestep between bins in n_bins, in seconds.

    noisy : boolean
        True if there is Poisson noise in the power spectrum (i.e., from real
        data), False if there is no noise in the power spectrum (i.e.,
        simulations without Poisson noise).

    Returns
    -------
    np.array of floats
        The noise-subtracted power spectrum in fractional rms^2 units.

    """
    if noisy:
        noise = 2.0 / mean_rate
    else:
        noise = 0.0

    return power * (2.0 * dt / float(n_bins) / (mean_rate ** 2)) - noise


################################################################################
def raw_to_leahy(power, mean_rate, n_bins, dt, noisy):
    """
    Normalizes the power spectrum to Leahy normalization.

    Parameters
    ----------
    power : np.array of floats
        The raw power at each Fourier frequency.

    mean_rate : float
        The mean count rate for the light curve, in cts/s.

    n_bins : int
        Number of bins per segment of light curve.

    dt : float
        Timestep between bins in n_bins, in seconds.

    noisy : boolean
        True if there is Poisson noise in the power spectrum (i.e., from real
        data), False if there is no noise in the power spectrum (i.e.,
        simulations without Poisson noise).

    Returns
    -------
    np.array of floats
        The noise-subtracted power spectrum in Leahy units.

    """
    if noisy:
        noise = 2.0
    else:
        noise = 0.0

    return power * (2.0 * dt / float(n_bins) / mean_rate) - noise


################################################################################
def var_and_rms(power, df):
    """
    Computes the variance and rms (root mean square) of a power spectrum.
    Assumes the negative-frequency powers have been removed.

    Parameters
    ----------
    power : np.array of floats
        The absolute rms normalized power at each of the *positive* Fourier
        frequencies.

    df : float
        The step size between Fourier frequencies.

    Returns
    -------

    float
        The variance of the power spectrum.
    float
        The RMS of the power spectrum.

    """
    variance = np.sum(power * df)
    rms = np.sqrt(variance)

    return variance, rms


################################################################################
def normalize(power, meta_dict, mean_rate, noisy):
    """
    Generates the Fourier frequencies, removes negative frequencies, normalizes
    the power by Leahy and fractional rms^2 normalizations, and computes the
    error on the fractional rms^2 power.

    Parameters
    ----------
    power : np.array of floats
        Description.

    meta_dict : dict
        Description.

    mean_rate : float
        Description.

    noisy : boolean
        Description

    Returns
    -------
    freq : np.array of floats
        Description.

    power : np.array of floats
        Description.

    leahy_power : np.array of floats
        Description.

    fracrms_power : np.array of floats
        Description.

    fracrms_err : np.array of floats
        Description.

    rms : float
        The fractional RMS of the noise-subtracted power spectrum.

    """

    ## Computing the FFT sample frequencies (in Hz)
    freq = fftpack.fftfreq(meta_dict['n_bins'], d=meta_dict['dt'])

    ## Ensuring that we're only using and saving the positive frequency
    ## values (and associated power values)
    nyq_index = meta_dict['n_bins'] / 2
    freq = np.abs(freq[0:nyq_index + 1])  ## because it slices at end-1, and we
        ## want to include 'nyq_index'; abs is because the nyquist freq is both
        ## pos and neg, and we want it pos here.
    power = power[0:nyq_index + 1]

    ## Computing the error on the mean power
    err_power = power / np.sqrt(float(meta_dict['n_seg']))

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

    ## Variance and rms of the whole averaged power spectrum
    total_variance = np.sum(fracrms_power * meta_dict['df'])
    print "Total variance:", total_variance, "(frac rms2)"
    rms_total = np.sqrt(total_variance)
    print "Total RMS:", rms_total, "(frac rms2)"

# 	print "Mean power:", np.mean(fracrms_power[np.where(freq >= 100)])

    return freq, power, leahy_power, fracrms_power, fracrms_err, rms_total


################################################################################
def make_ps(rate):
    """
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
    whole_lc = psd_lc.Lightcurve(n_bins=meta_dict['n_bins'])
    n_seg = 0
    exposure = 0
    dt_whole = np.array([])
    df_whole = np.array([])
    i = 0  # start of bin index to make segment of data for inner for-loop
    j = meta_dict['n_bins']  # end of bin index to make segment of data for inner for-loop

    # print data[1].field(0) - data[0].field(0)
    # meta_dict['dt'] = 0.015625
    # meta_dict['dt'] = 0.0078125
    # assert meta_dict['dt'] == (data[1].field(0) - data[0].field(0)), "ERROR: "\
    #     "dt must be the same resolution as the extracted FITS data."

    ## Loop through segments of the data
    while j <= len(data.field(1)):  ## while we haven't reached the end of the
                                    ## file

        n_seg += 1
        start_time = data.field(0)[i]
        end_time = data.field(0)[j-1]
        ## Extracts the second column of 'data' and assigns it to 'rate'.
        rate = data[i:j].field(1)

        power_segment, mean_rate_segment = make_ps(rate)
        assert int(len(power_segment)) == meta_dict['n_bins'], "ERROR: "\
                    "Something went wrong in make_ps. Length of power spectrum"\
                    " segment  != n_bins."

        dt_seg = (end_time - start_time) / float(meta_dict['n_bins'])
        df_seg = 1.0 / (meta_dict['n_bins'] * dt_seg)

        ## Computing variance and rms of the positive-frequency power in the
        ## reference band. Only keeping segments where the variance > 0.
        absrms_pow = raw_to_absrms(power_segment[0:meta_dict['n_bins']/2+1],
                mean_rate_segment, meta_dict['n_bins'], dt_seg, noisy=True)

        var, rms = var_and_rms(absrms_pow, df_seg)

        if var >= 0.0:
            whole_lc.power += power_segment
            whole_lc.mean_rate += mean_rate_segment

            exposure += end_time - start_time

            dt_whole = np.append(dt_whole, dt_seg)
            df_whole = np.append(df_whole, df_seg)

            if n_seg % print_iterator == 0:
                print "\t", n_seg

            if (test == True) and (n_seg == 1):  # For testing
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

    return whole_lc, n_seg, exposure, dt_whole, df_whole


################################################################################
def fits_in(in_file, meta_dict, print_iterator=int(5), test=False,
        chan_bounds=None, pcu=None):
    """
    Opens the .fits GTI'd event list, reads the count rate for a segment,
    populates the light curve, calls 'make_ps' to create a power spectrum, adds
    power spectra over all segments (average taken later).

    I take the approach: start time <= segment < end_time, to avoid double-
    counting and/or skipping events.

    Parameters
    ----------
    in_file : string
        The full path of the FITS data file being analyzed.

    meta_dict : dictionary
        Control parameters for the data analysis.

    print_iterator : int


    test : boolean
        True if only running one segment of data for testing, False if analyzing
        the whole data file. Default=False

    chan_bounds : list of ints

    pcu : int


    Returns
    -------
    whole_lc, n_seg, exposure, dt_whole, df_whole

    np.array of floats
        The sum of the power spectra across the segments.

    float
        The count rate of the light curve.

    int
        Number of segments in this data file.

    """

    try:
        fits_hdu = fits.open(in_file)
    except IOError:
        print "\tERROR: File does not exist: %s" % in_file
        exit()

    # header = fits_hdu[0].header	 ## Header is in ext 0
    data = fits_hdu[1].data  ## Data is in ext 1
    fits_hdu.close()

    whole_lc = psd_lc.Lightcurve(n_bins=meta_dict['n_bins'])
    n_seg = 0
    lightcurve = np.array([])
    exposure = 0
    dt_whole = np.array([])
    df_whole = np.array([])

    start_time = data.field('TIME')[0]
    final_time = data.field('TIME')[-1]
    end_time = start_time + meta_dict['n_seconds']

    ## Filter data based on pcu
    if pcu is not None:
        PCU2_mask = data.field('PCUID') == int(pcu)
        data = data[PCU2_mask]

    ## Filter data based on energy channel (event mode binned energy channel)
    if chan_bounds is not None:
        # print np.shape(data)
        lower_bound = data.field('CHANNEL') >= chan_bounds[0]
        data = data[lower_bound]
        upper_bound = data.field('CHANNEL') <= chan_bounds[1]
        data = data[upper_bound]
        # print np.shape(data)

    all_time = np.asarray(data.field('TIME'), dtype=np.float64)

    ################################
    ## Looping through the segments
    ################################

    while (end_time + (meta_dict['adjust_seg'] * meta_dict['dt'])) <= final_time:

        ## Adjusting segment length to artificially line up the QPOs
        end_time += (meta_dict['adjust_seg'] * meta_dict['dt'])

        time = all_time[np.where(all_time < end_time)]
        for_next_iteration = np.where(all_time >= end_time)
        all_time = all_time[for_next_iteration]

        if len(time) > 0:
            n_seg += 1
            rate_1d = tools.make_1Dlightcurve(time, meta_dict['n_bins'], \
                    start_time, end_time)
            lightcurve = np.concatenate((lightcurve, rate_1d))

            power_segment, mean_rate_segment = make_ps(rate_1d)
            assert int(len(power_segment)) == meta_dict['n_bins'], "ERROR: "\
                    "Something went wrong in make_ps. Length of power spectrum"\
                    " segment  != n_bins."

            dt_seg = (end_time - start_time) / float(meta_dict['n_bins'])
            df_seg = 1.0 / (meta_dict['n_bins'] * dt_seg)

            ## Computing variance and rms of the positive-frequency power in the
            ## reference band. Only keeping segments where the variance > 0.
            absrms_pow = raw_to_absrms(power_segment[0:meta_dict['n_bins']/2+1],
                    mean_rate_segment, meta_dict['n_bins'], dt_seg, noisy=True)

            var, rms = var_and_rms(absrms_pow, df_seg)

            if var >= 0.0:

                whole_lc.power += power_segment
                whole_lc.mean_rate += mean_rate_segment
                exposure += end_time - start_time
                dt_seg = (end_time - start_time) / float(meta_dict['n_bins'])
                dt_whole = np.append(dt_whole, dt_seg)
                df_whole = np.append(df_whole, 1.0 / (meta_dict['n_bins'] * dt_seg))

                ## Printing out which segment we're on every x segments
                if n_seg % print_iterator == 0:
                    print "\t", n_seg
                if test and (n_seg == 1):  # Testing
                    np.savetxt('tmp_lightcurve.dat', lightcurve, fmt='%d')
                    break
            ## Clearing variables from memory
            time = None
            power_segment = None
            mean_rate_segment = None
            rate_1d = None

            start_time = end_time
            end_time += meta_dict['n_seconds']

        ## This next bit deals with gappy data
        elif len(time) == 0:
            print "No counts in this segment."
            start_time = all_time[0]
            end_time = start_time + meta_dict['n_seconds']

    return whole_lc, n_seg, exposure, dt_whole, df_whole


################################################################################
def read_and_use_segments(in_file, meta_dict, test=False, chan_filt=None,
        pcu=None):
    """
    Opens the file, reads in the count rate, calls 'make_ps' to create a
    power spectrum. Separated from main body like this so I can easily call it
    in multi_powerspec.py. Split into 'fits_in' and 'extracted_in' for easier
    readability.

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

    ## Data from fits need to be populated as a light curve before a
    ## power spectrum can be taken, whereas data from lc was made in seextrct
    ## and so it's already populated as a light curve
    if ".fits" in in_file:
        whole_lc, n_seg, exposure, dt_whole, df_whole = fits_in(in_file,
                meta_dict, print_iterator=print_iterator, test=test)
                # chan_filt=chan_filt, pcu=pcu)

    elif ".lc" in in_file:
        if chan_filt is not None:
            raise Warning("Not able to filter on energy channel for an .lc"\
                    "file at this stage, since energy channel "\
                    "information is not given.")
        if pcu is not None:
            raise Warning("Not able to filter on RXTE PCU for an .lc file "\
                    "at this stage, since PCU information is not given.")

        whole_lc, n_seg, exposure, dt_whole, df_whole = extracted_in(in_file,
                meta_dict, print_iterator, test)

    else:
        raise Exception("ERROR: Input file type not recognized. Must be .fits "\
                "or .lc.")

    return whole_lc, n_seg, exposure, dt_whole, df_whole


################################################################################
def main(input_file, out_file, n_seconds, dt_mult, test=False, adjust=False,
        lo_chan=None, up_chan=None, pcu=None):
    """
    Reads in one data file at a time, takes FFT of segments of light curve data,
    computes power of each segment, averages power over all segments of all data
    files, writes resulting normalized power spectrum to a file.

    """
    #####################################################
    ## Idiot checks, to ensure that our assumptions hold
    #####################################################

    assert n_seconds > 0, "ERROR: Number of seconds per segment must be a "\
            "positive integer."
    assert dt_mult > 0, "ERROR: Multiple of dt must be a positive integer."

    ###########################
    ## Reading in data file(s)
    ###########################

    if ".txt" in input_file or ".lst" in input_file or ".dat" in input_file:
        data_files = [line.strip() for line in open(input_file)]
        if not data_files:  ## If data_files is an empty list
            raise Exception("ERROR: No files in the list of event lists.")
    else:
        data_files = [input_file]

    if adjust is True and len(data_files) == 9:
        adjust_segments = [932, 216, 184, 570, 93, 346, 860, 533, -324]
    else:
        adjust_segments = np.zeros(len(data_files))

    ###################
    ## Initializations
    ###################

    t_res = float(tools.get_key_val(data_files[0], 0, 'TIMEDEL'))

    try:
        detchans = int(tools.get_key_val(data_files[0], 0, 'DETCHANS'))
    except IOError:
        detchans = 64

    meta_dict = {'dt': dt_mult * t_res,
                 't_res': t_res,
                 'n_seconds': n_seconds,
                 'df': 1.0 / float(n_seconds),
                 'nyquist': 1.0 / (2.0 * dt_mult * t_res),
                 'n_bins': n_seconds * int(1.0 / (dt_mult * t_res)),
                 'detchans': detchans}

    print "\nDT = %f seconds" % meta_dict['dt']
    print "N_bins = %d" % meta_dict['n_bins']
    print "Nyquist freq =", meta_dict['nyquist']
    print "Testing?", test
    print "Adjusting QPO?", adjust

    # if lo_chan is not None:
    #     if up_chan is not None:
    #         assert lo_chan <= up_chan, "ERROR: lo_energy ! <= up_energy"
    #         chan_bounds = [lo_chan, up_chan]
    #
    #     else:
    #         chan_bounds = [lo_chan, detchans-1]
    #
    # else:
    #     if up_chan is not None:
    #         chan_bounds = [0, up_chan]

    total_seg = 0
    total_exposure = 0
    dt_total = np.array([])
    df_total = np.array([])
    total = psd_lc.Lightcurve(n_bins=meta_dict['n_bins'])

    #####################################
    ## Looping through files starts here
    #####################################

    for in_file, adj_seg in zip(data_files, adjust_segments):

        meta_dict['adjust_seg'] = adj_seg

        whole_lc, n_seg, exposure, dt_whole, \
                df_whole = read_and_use_segments(in_file, meta_dict, test=test)

        print "Segments for this file: %d\n" % n_seg

        total_seg += n_seg
        total_exposure += exposure
        dt_total = np.append(dt_total, dt_whole)
        df_total = np.append(df_total, df_whole)

        total.power += whole_lc.power
        total.mean_rate += whole_lc.mean_rate

    ## End of for-loop
    print " "

    meta_dict['n_seg'] = total_seg

    total.power /= float(meta_dict['n_seg'])
    total.mean_rate /= float(meta_dict['n_seg'])

    meta_dict['exposure'] = total_exposure
    meta_dict['mean_rate'] = total.mean_rate
    print "Total exposure time =", meta_dict['exposure']
    print "Total segments =", meta_dict['n_seg']
    print "Mean rate total =", meta_dict['mean_rate']

    print meta_dict['dt']
    print meta_dict['df']

    ##################################################
    ## Normalize the power spectrum and compute error
    ##################################################

    total_variance = np.sum(total.power * meta_dict['df'])
    print "Total variance:", total_variance, "(unnorm)"
    rms_total = np.sqrt(total_variance)
    print "Total RMS:", rms_total, "(unnorm)"

    freq, power, leahy_power, fracrms_power, fracrms_err, rms = \
            normalize(total.power, meta_dict, total.mean_rate, True)

    meta_dict['rms'] = rms
    print freq[meta_dict['n_bins']/2]
    meta_dict['nyquist'] = freq[meta_dict['n_bins']/2]
    print meta_dict['nyquist']

    ##########
    ## Output
    ##########

    if len(data_files) == 1:
        file_description = "Power spectrum of one observation"
    else:
        file_description = "Power spectrum of multiple observations"

    fits_out(out_file, input_file, meta_dict, freq, fracrms_power, fracrms_err,
            leahy_power, file_description)


################################################################################
if __name__ == "__main__":

    #########################################
    ## Parse input arguments and call 'main'
    #########################################

    parser = argparse.ArgumentParser(usage="powerspec.py infile outfile "\
            "[OPTIONAL ARGUMENTS]", description=__doc__, epilog="For "\
            "optional arguments, default values are given in brackets at end "\
            "of description.")

    parser.add_argument('infile', help="Could be either: 1) The full path of "\
            "the input file containing an event list: .fits has time in "\
            "column 1, energy channel in column 2, detector ID in column 3; or"\
            " .lc has time in column 1 and count rate in column 2; or 2) The "\
            "full path of the (ASCII/txt) file with a list of the input files"\
            " (as described in 1). One file per line.")

    parser.add_argument('outfile', help="The full path of the .fits file to "\
            "write the frequency and power to.")

    parser.add_argument('-n', '--n_seconds', type=tools.type_power_of_two,
            default=64, dest='n_seconds', help="Number of seconds in each "\
            "Fourier segment. Must be a power of 2, positive, integer. [64]")

    parser.add_argument('-m', '--dt_mult', type=tools.type_power_of_two,
            default=64, dest='dt_mult', help="Multiple of dt (dt is from data "\
            "file) for timestep between bins. Must be a power of 2, positive, "\
            "integer. [64]")

    parser.add_argument('-t', '--test', type=int, default=0, choices={0,1},
            dest='test', help="Int flag: 0 if computing all segments, 1 if "\
            "only computing one segment for testing. [0]")

    parser.add_argument('-a', '--adjust', default=False, action='store_true',
            dest='adjust', help="If present, artificially adjusts the "\
            "frequency of the QPO by changing the segment length. [False]")

    parser.add_argument('--le', dest='lo_energy', type=tools.type_positive_int,
            default=None, help="If present, the lower bound of the energy "\
            "channel range to compute the averaged power spectrum over. [None]")

    parser.add_argument('--ue', dest='up_energy', type=tools.type_positive_int,
            default=None, help="If present, the upper bound of the energy "\
            "channel range to compute the averaged power spectrum over. [None]")

    parser.add_argument('--pcu', dest='pcu', choices={0,1,2,3,4}, default=None,
            help="If present, the RXTE PCU to make a power spectrum for. "\
            "[None]")


    args = parser.parse_args()

    testing = False
    if args.test == 1:
        testing = True

    main(args.infile, args.outfile, args.n_seconds, args.dt_mult,
            test=testing, adjust=args.adjust, lo_chan=args.lo_energy,
            up_chan=args.up_energy, pcu=args.pcu)

################################################################################
