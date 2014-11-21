import argparse
import numpy as np
from scipy import fftpack
from astropy.io import fits
from datetime import datetime
import os

import tools

"""
		powerspec.py

Makes a power spectrum from an event-mode data file from RXTE.

Required arguments:
datafile - Name of FITS file with photon count rate data.
outfile - Name of file that the power spectrum will be written to.
rebinned_outfile - Name of file that the re-binned power spectrum will be 
	written to.

Optional arguments:
num_seconds - Number of seconds in each Fourier segment. Must be a power of 2.
rebin_const - Used to re-bin the data after the average power is computed, such 
	that bin_size[n+1] = bin_size[n] * rebin_const.
dt_mult - Multiple of 1/8192 seconds for the timestep between bins. Must be a 
	power of 2.
test - If present, only does a short test run.

Written in Python 2.7 by A.L. Stevens, A.L.Stevens@uva.nl, 2013-2014

The scientific modules imported above, as well as python 2.7, can be downloaded 
in the Anaconda package, https://store.continuum.io/cshop/anaconda/

tools is available at https://github.com/abigailStev/whizzy_scripts

"""

################################################################################
def dat_output(out_file, rebinned_out_file, in_file, dt, n_bins, nyquist_freq, 
	num_segments, mean_rate_whole, freq, rms2_power, rms2_err_power, 
	leahy_power, rebin_const, rebinned_freq, rebinned_rms2_power, 
	err_rebinned_power):
	""" 
			ascii_output
			
	Writes power spectrum and re-binned power spectrum to two ASCII output 
	files.
		
	Variables same as fits_output.
	
	"""
	print "Standard output file: %s" % out_file
	
	with open(out_file, 'w') as out:
		out.write("#\t\tPower spectrum")
		out.write("\n# Date(YYYY-MM-DD localtime): %s" % str(datetime.now()))
		out.write("\n# Data: %s" % in_file)
		out.write("\n# Time bin size = %.21f seconds" % dt)
		out.write("\n# Number of bins per segment = %d" % n_bins)
		out.write("\n# Number of segments per light curve = %d" % num_segments)
		out.write("\n# Exposure time = %d seconds" \
			% (num_segments * n_bins * dt))
		out.write("\n# Mean count rate = %.8f" % mean_rate_whole)
		out.write("\n# Nyquist frequency = %.4f" % nyquist_freq)
		out.write("\n# ")
		out.write("\n# Column 1: Frequency [Hz]")
		out.write("\n# Column 2: Fractional rms^2 normalized mean power")
		out.write("\n# Column 3: Fractional rms^2 normalized error on the mean \
			power")
		out.write("\n# Column 4: Leahy-normalized mean power")
		out.write("\n# ")
		for k in xrange(len(rms2_power)):
			if freq[k] >= 0:
# 				out.write("\n%.8f\t%.8f\t%.8f" % (freq[k], rms2_power[k], 
# 					rms2_err_power[k]))
				out.write("\n{0:.8f}\t{1:.8f}\t{2:.8f}\t{3:.8f}".format(freq[k], 
					rms2_power[k], rms2_err_power[k], leahy_power[k]))
			## End of if-statement
		## End of for-loop
	## End of with-block
	
	## Now outputting the re-binned data -- Need to do this separately since it
	## has a different number of data points from the regular power spectrum.
	
	print "Re-binned output file: %s" % rebinned_out_file
	
	with open(rebinned_out_file, 'w') as out:
		out.write("#\t\tPower spectrum")
		out.write("\n# Data: %s" % in_file)
		out.write("\n# Re-binned in frequency at (%.4f * prev bin size)" \
			% rebin_const)
		out.write("\n# Corresponding un-binned output file: %s" % out_file)
		out.write("\n# Original time bin size = %.21f seconds" % dt)
		out.write("\n# Exposure time = %d seconds" \
			% (num_segments * n_bins * dt))
		out.write("\n# Mean count rate = %.8f" % mean_rate_whole)
		out.write("\n# ")
		out.write("\n# Column 1: Frequency [Hz]")
		out.write("\n# Column 2: Fractional rms^2 normalized mean power")
		out.write("\n# Column 3: Error in fractional rms^2 normalized mean \
			power")
		out.write("\n# ")
		for k in xrange(len(rebinned_rms2_power)):
			if rebinned_freq[k] >= 0:
				out.write("\n{0:.8f}\t{1:.8f}\t{2:.8f}".format(rebinned_freq[k],
                    rebinned_rms2_power[k], err_rebinned_power[k]))
			## End of if-statement
		## End of for-loop
	## End of with-block

## End of function 'dat_output'


################################################################################
def fits_output(out_file, rebinned_out_file, in_file, dt, n_bins, nyquist_freq, 
	num_segments, mean_rate_whole, freq, rms2_power, rms2_err_power, 
	leahy_power, rebin_const, rebinned_freq, rebinned_rms2_power, 
	err_rebinned_power):
	"""
				fits_output
			
	Writes power spectrum and re-binned power spectrum to two FITS files.
		
	Passed: out_file - Name of output file for standard power spectrum.
			rebinned_out_file - Name of output file for geometrically re-binned 
				power spectrum.
			in_file - Full path of filename containing input data.
			dt - Size of time bin, in seconds (must be power of 2).
			n_bins - Number of time bins in a segment (must be power of 2).
			nyquist_freq - Nyquist frequency = 1/(2*dt)
			num_segments - Number of segments in the light curve.
			mean_rate_whole - Mean count rate of the light curve.
			freq - Frequencies (in Hz) corresponding to the power spectrum.
			leahy_power - Leahy-normalized power, averaged over all
				segments of the light curve.
			rms2_power - Fractional rms power, averaged over all 
				segments of the light curve and all light curves.
			rms2_err_power - Error on avg fractional rms power.
			rebin_const - Constant >1 by which we want to re-bin the spectrum,
				such that bin_size[n+1] = bin_size[n] * rebin_const.
			rebinned_freq - Frequencies of power spectrum re-binned
				according to rebin_const.
			rebinned_rms2_power - Power spectrum re-binned according to 
				rebin_const.
			err_rebinned_power - Error on re-binned fractional rms 
				power.

	Returns: nothing
	"""
	print "Standard output file: %s" % out_file

	## Making header for standard power spectrum
	prihdr = fits.Header()
	prihdr.set('TYPE', "Power spectrum")
	prihdr.set('DATE', str(datetime.now()), "YYYY-MM-DD localtime")
	prihdr.set('EVTLIST', in_file)
	prihdr.set('DT', dt, "seconds")
	prihdr.set('N_BINS', n_bins, "time bins per segment")
	prihdr.set('SEGMENTS', num_segments, "segments in the whole light curve")
	prihdr.set('EXPOSURE', num_segments * n_bins * dt, "seconds, of light curve")
	prihdr.set('MEANRATE', mean_rate_whole, "counts / second")
	prihdr.set('NYQUIST', nyquist_freq, "Hz")
	prihdu = fits.PrimaryHDU(header=prihdr)
	
	## Making FITS table for standard power spectrum
	col1 = fits.Column(name='FREQUENCY', unit='Hz', format='E', array=freq)
	col2 = fits.Column(name='POWER', unit='frac rms^2', format='E', \
		array=rms2_power)
	col3 = fits.Column(name='ERROR', unit='frac rms^2', format='E', \
		array=rms2_err_power)
	col4 = fits.Column(name='LEAHY', format='E', array=leahy_power)
	cols = fits.ColDefs([col1, col2, col3, col4])
	tbhdu = fits.BinTableHDU.from_columns(cols)
	
	## If the file already exists, remove it (still working on just updating it)
	assert out_file[-4:].lower() == "fits", \
		'ERROR: Standard output file must have extension ".fits".'
	if os.path.isfile(out_file):
# 		print "File previously existed. Removing and rewriting."
		os.remove(out_file)
		
	## Writing the standard power spectrum to a FITS file
	thdulist = fits.HDUList([prihdu, tbhdu])
	thdulist.writeto(out_file)	
	
	## Now outputting the re-binned data -- Need to do this separately since it
	## has a different number of data points from the regular power spectrum.
	print "Re-binned output file: %s" % rebinned_out_file

	## Updating above header for re-binned power spectrum
	prihdr.set('TYPE', "Re-binned power spectrum")
	prihdr.insert('DATE', ('UNBINOUT', out_file, "Corresponding un-binned output."))
	prihdr.insert('UNBINOUT', ('REBIN', rebin_const, "Freqs re-binned by REBIN * prev_bin_size"))
	prihdu = fits.PrimaryHDU(header=prihdr)
	
	## Making FITS table for re-binned power spectrum
	col1 = fits.Column(name='FREQUENCY', unit='Hz', format='E', \
		array=rebinned_freq)
	col2 = fits.Column(name='POWER', unit='frac rms^2', format='E', \
		array=rebinned_rms2_power)
	col3 = fits.Column(name='ERROR', unit='frac rms^2', format='E', \
		array=err_rebinned_power)
	cols = fits.ColDefs([col1, col2, col3])
	tbhdu = fits.BinTableHDU.from_columns(cols)
	
	## If the file already exists, remove it (still working on just updating it)
	assert rebinned_out_file[-4:].lower() == "fits", \
		'ERROR: Re-binned output file must have extension ".fits".'
	if os.path.isfile(rebinned_out_file):
# 		print "File previously existed. Removing and rewriting."
		os.remove(rebinned_out_file)
	
	## Writing the re-binned power spectrum to a FITS file
	thdulist = fits.HDUList([prihdu, tbhdu])
	thdulist.writeto(rebinned_out_file)	
	
	## End of function 'fits_output'
	

################################################################################
def geometric_rebinning(freq, rms2_power, rms2_err_power, rebin_const):
	"""
			geometric_rebinning
			
	Re-bins the noise-subtracted fractional rms^2 power spectrum in frequency 
	space by some re-binning constant (rebin_const>1). 

	Passed: rms2_power - Fractional rms^2 power, averaged over all segments
				in the light curve.
			rms2_err_power - Error on the fractional rms^2 power.
			freq - Frequencies (in Hz) corresponding to the power spectrum.
			rebin_const - Constant >1 by which we want to re-bin the spectrum, 
				such that bin_size[n+1] = bin_size[n] * rebin_const.
	
	Returns: rebinned_freq - Frequencies of power spectrum re-binned according 
				to rebin_const.
			 rebinned_rms2_power - Power spectrum re-binned according to 
			 	rebin_const.
			 err_rebinned_power - Error on the re-binned rms power.
	
	"""
	## Initializing variables
	rebinned_rms2_power = []	 # List of re-binned fractional rms power
	rebinned_freq = []		 # List of re-binned frequencies
	err_rebinned_power = []	 # List of error in re-binned power
	real_index = 1.0		 # The unrounded next index in power
	int_index = 1			 # The int of real_index, added to current_m every 
							 #  iteration
	current_m = 1			 # Current index in power
	prev_m = 0				 # Previous index m
	bin_power = 0.0			 # The power of the current re-binned bin
	bin_freq = 0.0			 # The frequency of the current re-binned bin
	err_bin_power2 = 0.0	 # The error squared on 'bin_power'
	bin_range = 0.0			 # The range of un-binned bins covered by this 
							 #  re-binned bin
	
	## Looping through the length of the array power, geometric bin by 
	## geometric bin, to compute the average power and frequency of that
	## geometric bin.
	## Equations for frequency, power, and error on power are from Adam Ingram's
	## PhD thesis
	while current_m < len(rms2_power):
# 	while current_m < 100: # used for debugging
		## Initializing clean variables for each iteration of the while-loop
		bin_power = 0.0  # the averaged power at each index of rebinned_power
		err_bin_power2 = 0.0  # the square of the errors on powers in this bin
		bin_range = 0.0
		bin_freq = 0.0
		
		## Looping through the data points contained within one geometric bin
		for k in xrange(prev_m, current_m):
			bin_power += rms2_power[k]
			err_bin_power2 += rms2_err_power[k] ** 2
		
		## Determining the range of indices this specific geometric bin covers
		bin_range = np.absolute(current_m - prev_m)
		## Dividing bin_power (currently just a sum of the data points) by the 
		## number of points to get an arithmetic average
		bin_power /= float(bin_range)
		
		## Computing the mean frequency of a geometric bin
		bin_freq = ((freq[current_m] - freq[prev_m]) / bin_range) + freq[prev_m]

		## If there's only one data point in the geometric bin, there's no need
		## to take an average. This also prevents it from skipping the first 
		## data point.
		if bin_range == 1:
			bin_power = rms2_power[prev_m]
			bin_freq = freq[prev_m]
		
		## Appending values to arrays
		rebinned_rms2_power.append(bin_power)
		rebinned_freq.append(bin_freq)
		## Computing error in bin -- equation from Adam Ingram's thesis
		err_rebinned_power.append(np.sqrt(err_bin_power2) / float(bin_range))
		
		## Incrementing for the next iteration of the loop
		prev_m = current_m
		real_index *= rebin_const
		int_index = int(round(real_index))
		current_m += int_index
		## Since the for-loop goes from prev_m to current_m-1 (since that's how
		## the range function works) it's ok that we set prev_m = current_m here
		## for the next round. This will not cause any double-counting bins or 
		## skipping bins.
		bin_range = None
		bin_freq = None
		bin_power = None
		err_bin_power2 = None
		
	## End of while-loop
		
	return rebinned_freq, rebinned_rms2_power, err_rebinned_power
## End of function 'geometric_rebinning'


################################################################################
def make_ps(rate):
	"""
			make_ps
	
	Computes the mean count rate, the FFT of the count rate minus the mean, and 
	the power spectrum of this segment of data.
	
	Passed: rate - Count rate for this segment of data.
	
	Returns: power_segment - Power spectra for this segment of data.
			 mean_rate - Mean count rate for this segment of data.
	
	""" 
	## Computing the mean count rate of the segment
	mean_rate = np.mean(rate)

	## Subtracting the mean rate off each value of 'rate'
	##  This eliminates the spike at 0 Hz 
	rate_sub_mean = rate - mean_rate
	
	## Taking the 1-dimensional FFT of the time-domain photon count rate
	##  Using the SciPy FFT algorithm, as it's faster than NumPy for large lists
	fft_data = fftpack.fft(rate_sub_mean)

	## Computing the power
	power_segment = np.absolute(fft_data) ** 2
	
	return power_segment, mean_rate
## End of function 'make_ps'


################################################################################
def extracted_powerspec(in_file, n_bins, dt, print_iterator, test):
	"""
			extracted_powerspec
	
	Opens the FITS file light curve (as created in seextrct), reads the count 
	rate for a segment, calls 'make_ps' to create a power spectrum, adds power 
	spectra over all segments.
	
	Passed: in_file - Name of .lc (type FITS) input file/light curve.
			n_bins - Number of integer bins per segment. Must be a power of 2.
			dt - Time step between bins, in seconds. Must be a power of 2.
			print_iterator - Prints out which segment we're on every 
				'print_iterator' segments, to keep track of progress.
			test - True if computing one segment, False if computing all.
	
	Returns: power_sum - Sum of power spectra for all segments of data from this
				input file.
			 sum_rate_whole - Sum of the mean count rate for all segments of 
			 	data from this input file.
			 num_segments - Number of segments of data from this input file.
			 
	"""	
	fits_hdu = fits.open(in_file)
	header = fits_hdu[1].header	
	data = fits_hdu[1].data
	fits_hdu.close()
	
	sum_rate_whole = 0
	power_sum = np.zeros(n_bins, dtype=np.float64)
	num_segments = 0

	i = 0  # start of bin index to make segment of data for inner for-loop
	j = n_bins  # end of bin index to make segment of data for inner for-loop

	assert dt == (data[1].field(0) - data[0].field(0)), 'ERROR: Specified dt must be same resolution as the FITS data'
# 	print "Length of FITS file:", len(data.field(1))
	while j <= len(data.field(1)):  
	
		num_segments += 1

		## Extracts the second column of 'data' and assigns it to 'rate'. 
		rate = data[i:j].field(1)
		
		power_segment, mean_rate_segment = make_ps(rate)
		power_sum += power_segment
		sum_rate_whole += mean_rate_segment
		
		if num_segments % print_iterator == 0:
			print "\t", num_segments
			
		if (test == True) and (num_segments == 1):  # For testing
			break

		rate = None
		power_segment = None
		mean_rate_segment = None
		## Incrementing the counters and indices
		i = j
		j += n_bins
		## Since the for-loop goes from i to j-1 (since that's how the range 
		## function works) it's ok that we set i=j here for the next round. 
		## This will not cause double-counting rows or skipping rows.
       
	## End of while-loop
		
	return power_sum, sum_rate_whole, num_segments
## End of function 'extracted_powerspec'


################################################################################
def fits_powerspec(in_file, n_bins, dt, print_iterator, test):
	"""
			fits_powerspec
	
	Opens the .fits GTI'd event list, reads the count rate for a segment, 
	populates the light curve, calls 'make_ps' to create a power spectrum, adds 
	power spectra over all segments.
	
	Passed: in_file - Name of .fits input file/light curve.
			n_bins - Number of integer bins per segment. Must be a power of 2.
			dt - Time step between bins, in seconds. Must be a power of 2.
			print_iterator - Prints out which segment we're on every 
				'print_iterator' segments, to keep track of progress.
			test - True if computing one segment, False if computing all.
	
	Returns: power_sum - Sum of power spectra for all segments of data from this
				input file.
			 sum_rate_whole - Sum of the mean count rate for all segments of 
			 	data from this input file.
			 num_segments - Number of segments of data from this input file.
			 
	"""	
	fits_hdu = fits.open(in_file)
	header = fits_hdu[0].header	 # Header info is in ext 0, data is in ext 1
	data = fits_hdu[1].data
	fits_hdu.close()
	
	sum_rate_whole = 0
	power_sum = np.zeros(n_bins, dtype=np.float64)
	num_segments = 0
	lightcurve = np.asarray([])

	start_time = data.field('TIME')[0]
	final_time = data.field('TIME')[-1]
	end_time = start_time + (dt * n_bins)
	
# 	print "First start time: %.21f" % start_time
# 	print "First end   time: %.21f" % end_time

	PCU2_mask = data.field('PCUID') == 2
	data = data[PCU2_mask]
	time = np.asarray(data.field('TIME'), dtype=np.float64)
	energy = np.asarray(data.field('CHANNEL'), dtype=np.float64)

	while end_time <= final_time:
			
		this_time = time[np.where(time <= end_time)]
		this_energy = energy[np.where(time <= end_time)]
		
		for_next_iteration = np.where(time > end_time)
		time = time[for_next_iteration]
		energy = energy[for_next_iteration]

# 		print "Len of time:", len(time)
# 		print "Len of this time:", len(this_time)
# 		print "\nLen time:", len(this_time)
# 		print "Start: %.21f" % start_time
# 		print "End  : %.21f" % end_time
# 		if len(this_time) == 1:
# 			print "One: %.21f" % this_time[0]# 		if len(this_time) == 0:
# 			print "Len of time = 0"
		if len(this_time) > 0:
			num_segments += 1
			rate_2d, rate_1d = tools.make_lightcurve(this_time, 
				this_energy, n_bins, dt, start_time)
			lightcurve = np.concatenate((lightcurve, rate_1d))
						
			power_segment, mean_rate_segment = make_ps(rate_1d)
			assert int(len(power_segment)) == n_bins

			power_sum += power_segment
			sum_rate_whole += mean_rate_segment
			## Printing out which segment we're on every x segments
			if num_segments % print_iterator == 0:
				print "\t", num_segments
# 				print "\t", len(time)
# 				print "\t", len(this_time)
# 				print "\t%.21f" % end_time


			## Clearing variables from memory
			this_time = None
			this_energy = None
			power_segment = None
			mean_rate_segment = None
			rate_2d = None
			rate_1d = None
			
			if test and (num_segments == 1):  # Testing
				np.savetxt('lightcurve.dat', lightcurve, fmt='%d')
				break
		## End of 'if there are counts in this segment'
		
		if num_segments == 170:
			break
		
		start_time += (n_bins * dt)
		end_time += (n_bins * dt)

	## End of while-loop
# 	print "Final end time: %.21f" % end_time
	return power_sum, sum_rate_whole, num_segments
## End of function 'fits_powerspec'


################################################################################
def dat_powerspec(in_file, n_bins, dt, print_iterator, test):
	"""
			dat_powerspec
		
	Opens the .dat GTI'd event list, reads the count rate for a segment, 
	populates the light curve, calls 'make_ps' to create a power spectrum, 
	adds power spectra over all segments.
		
	Assumes that all event times have been corrected with TIMEZERO, and 
	GTI-filtered.
	
	Passed: in_file - Name of .dat event list, still to be 'populated'.
			n_bins - Number of integer bins per segment. Must be a power of 2.
			dt - Desired time step between bins of populated light curve, in 
				seconds. Must be a power of 2.
			print_iterator - Prints out which segment we're on every 
				'print_iterator' segments, to keep track of progress.
			test - True if computing a few segments, False if computing all.
	
	Returns: power_sum - Sum of power spectra for all segments of data from this
				input file.
			 sum_rate_whole - Sum of the mean count rate for all segments of 
			 	data from this input file.
			 num_segments - Number of segments of data from this input file.
			 
	"""	
	## Declaring clean variables to append to for every loop iteration.
	time = np.asarray([])
	energy = np.asarray([])
	sum_rate_whole = 0
	power_sum = np.zeros(n_bins, dtype=np.float64)
	num_segments = 0
	
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

	end_time = start_time + (dt * n_bins)
	assert end_time > start_time, 'ERROR: End time must come after start time of the segment.'
	print "First start time: %.21f" % start_time
	print "First end   time: %.21f" % end_time

	with open(in_file, 'r') as f:
		for line, next_line in tools.pairwise(f):
			if line[0].strip() != "#" and \
				float(line.strip().split()[0]) >= start_time:  
				# If the line is not a comment
				line = line.strip().split()
				next_line = next_line.strip().split()
				current_time = np.float64(line[0])
				current_chan = np.int8(line[1])
				current_pcu = np.int8(line[2])
				
				if current_pcu == 2:  # Only want PCU2 here
					time = np.append(time, current_time)
					energy = np.append(energy, current_chan)
					
				next_time = float(next_line[0])
				next_end_time = end_time + (dt * n_bins)
				
				if next_time > end_time:  # Triggered at the end of a segment

					print "\nLen time:", len(time)
					print "Start: %.21f" % start_time
					print "End  : %.21f" % end_time
					if len(time) == 1:
						print "One: %.21f" % time[0]

					if len(time) > 0:
						num_segments += 1
						rate_2d, rate_1d = tools.make_lightcurve(time, energy, 
							n_bins, dt, start_time)
						lightcurve = np.concatenate((lightcurve, rate_1d))
					
						power_segment, mean_rate_segment = make_ps(rate_1d)
						assert int(len(power_segment)) == n_bins

						power_sum += power_segment
						sum_rate_whole += mean_rate_segment
						## Printing out which segment we're on every x segments
						if num_segments % print_iterator == 0:
							print "\t", num_segments
# 							print "\t", len(time)
							print "\t%.21f" % end_time

						## Clearing variables from memory
						power_segment = None
						mean_rate_segment = None
						rate_2d = None
						rate_1d = None
						time = np.asarray([])
						energy = np.asarray([])
					
						if test and (num_segments == 1):  # Testing
							np.savetxt('lightcurve.dat', lightcurve, fmt='%d')
							break
						if num_segments == 170:
							break
					## End of 'if there are counts in this segment'
					
					start_time += (n_bins * dt)
					end_time += (n_bins * dt)
					
					## This next bit helps it handle gappy data
					if next_time > end_time:
						while next_time > end_time:
							start_time += (n_bins * dt)
							end_time += (n_bins * dt)

				## End of 'if we're at the end of a segment'
			## End of 'if the line is not a comment'
		## End of for-loop 
	## End of with-block
	print "Final end time: %.21f" % end_time
	return power_sum, sum_rate_whole, num_segments
## End of function 'dat_powerspec'
	

################################################################################
def read_and_use_segments(in_file, n_bins, dt, test):
	"""
			read_and_use_segments
			
	Opens the file, reads in the count rate, calls 'make_ps' to create a
	power spectrum. Separated from main body like this so I can easily call it 
	in multi_powerspec.py. Split into 'fits_powerspec', 'dat_powerspec', and
	'extracted_powerspec' for easier readability.
	
	Passed: in_file - Name of input file with time in column 1 and rate in 
				column 2. Must have extension .dat, .fits, or .lc.
			n_bins - Number of integer bins per segment. Must be a power of 2.
			dt - Time step between bins, in seconds. Must be a power of 2.
			test - True if computing few segments, False if computing all.
	
	Returns: power_sum - Sum of power spectra for all segments of data from this
				input file.
			 sum_rate_whole - Sum of the mean count rate for all segments of 
			 	data from this input file.
			 num_segments - Number of segments of data from this input file.
			 
	"""
	assert tools.power_of_two(n_bins) , 'ERROR: n_bins must be a power of 2.'
	
	print "Input file: %s" % in_file
	
	if n_bins == 32768:
		print_iterator = int(10)
	elif n_bins < 32768:
		print_iterator = int(20)
	elif n_bins >= 2097152:
		print_iterator = int(1)	
	elif n_bins >= 1048576:
		print_iterator = int(2)
	else:
		print_iterator = int(5)
	
	## Looping through length of data file, segment by segment, to compute 
	## power for each data point in the segment
	print "Segments computed:"
	
	## data from dat and fits need to be populated as a light curve before a 
	## power spectrum can be taken, whereas data from lc was made in seextrct 
	## and so it's already populated as a light curve
	if (in_file[-5:].lower() == ".fits"):
		power_sum, sum_rate_whole, num_segments = fits_powerspec(in_file, 
			n_bins, dt, print_iterator, test)
	elif (in_file[-4:].lower() == ".dat"):
		power_sum, sum_rate_whole, num_segments = dat_powerspec(in_file, 
			n_bins, dt, print_iterator, test)
	elif (in_file[-3:].lower() == ".lc"):
		power_sum, sum_rate_whole, num_segments = extracted_powerspec(in_file, 
			n_bins, dt, print_iterator, test)
	else:
		raise Exception("ERROR: Input file type not recognized. Must be .dat, .fits. or .lc.")
	
	return power_sum, sum_rate_whole, num_segments
## End of function 'read_and_use_segments'


###############################################################################
def normalize(power, n_bins, dt, num_seconds, num_segments, mean_rate):
	"""
			normalize
	
	Generates the Fourier frequencies, removes negative frequencies, normalizes
	the power by Leahy and fractional rms^2 normalizations, and computes the 
	error on the fractional rms^2 power.
	
	Passed: power - The power spectrum averaged over all segments.
			n_bins - Number of bins per segment.
			dt - Timestep between bins, in seconds.
			num_seconds - Length of each Fourier segment, in seconds.
			num_segments - Number of segments the light curve is broken up into.
			mean_rate - The mean count rate over all segments of data.
	
	Returns: freq - The Fourier frequencies.
			 power - The power spectrum averaged over all segments (just the 
			 	ones with positive Fourier frequencies). 
			 leahy_power - The Leahy-normalized power.
			 rms2_power - The noise-subtracted fractional rms^2 power.
			 rms2_err_power - The error on the fractional rms^2 power.
	"""
	## Computing the FFT sample frequencies (in Hz)
	freq = fftpack.fftfreq(n_bins, d=dt)
	## Ensuring that we're only using and saving the positive frequency values 
	##  (and associated power values)
	max_index = np.argmax(freq)+1  # because in python, the scipy fft makes the 
								   # nyquist frequency negative, and we want it 
								   # to be positive! (it is actually both pos 
								   # and neg)
	freq = np.abs(freq[0:max_index + 1])  # because it slices at end-1, and we 
										  # want to include 'max_index'; abs is
										  # because the nyquist freq is both pos
										  # and neg, and we want it pos here.
	power = power[0:max_index + 1]
	
	## Computing the error on the mean power
	err_power = power / np.sqrt(float(num_segments) * float(n_bins))
	
	## Leahy normalization
	leahy_power = 2.0 * power * dt / float(n_bins) / mean_rate
	print "Mean value of Leahy power =", np.mean(leahy_power)  # ~2

	## Fractional rms^2 normalization with noise subtracted off
	rms2_power = 2.0 * power * dt / float(n_bins) / (mean_rate ** 2)
	rms2_power -= 2.0 / mean_rate
	
	## Error on fractional rms^2 power (not subtracting noise)
	rms2_err_power = 2.0 * err_power * dt / float(n_bins) / mean_rate ** 2
	
	df = 1.0 / float(num_seconds)  # in Hz
	signal_freq = freq[np.argmax(power)]  # Assumes that the signal dominates 
										  # the power spectrum.
	print "Frequency of maximum power:", signal_freq
	min_freq_mask = freq < signal_freq  # we want the last 'True' element
	max_freq_mask = freq > signal_freq  # we want the first 'True' element
	j_min = list(min_freq_mask).index(False)
	j_max = list(max_freq_mask).index(True)
	## Extracting only the signal frequencies of the power
	signal_pow = np.float64(rms2_power[j_min:j_max])
	## Computing variance and rms of the signal
	signal_variance = np.sum(signal_pow * df)
	rms = np.sqrt(signal_variance)  # should be a few % in frac rms units
	print "RMS of signal:", rms
	
	return freq, power, leahy_power, rms2_power, rms2_err_power
	## End of function 'normalize'
	
	
###############################################################################
def main(in_file, out_file, rebinned_out_file, num_seconds, rebin_const, 
	dt_mult, test):
	""" 
			main
	
	Passed: in_file - Name of input file with time in column 1 and rate in 
				column 2. FITS format must have extension .lc or .fits, 
				otherwise assumes .dat (ASCII/txt) format.
			out_file - Name of output file for standard power spectrum.
			rebinned_out_file - Name of output file for re-binned power 
				spectrum.
			num_seconds - Number of seconds in each Fourier segment. Must be a 
				power of 2.
			rebin_const - Used to re-bin the data geometrically after the 
				average power is computed, such that 
				bin_size[n+1] = bin_size[n] * rebin_const.
			dt_mult - Multiple of 1/8192 seconds for timestep between bins.
			test - True if computing 1 segment, False if computing all.
	
	Returns: nothing
	
	"""	
	assert rebin_const >= 1.0 , 'ERROR: Re-binning constant must be >= 1.' 
	
	t_res = 1.0 / 8192.0  # The time resolution of the data, in seconds
	dt = dt_mult * t_res
	n_bins = num_seconds * int(1.0 / dt)
	assert tools.power_of_two(n_bins), 'ERROR: n_bins must be a power of 2.'
	nyquist_freq = 1.0 / (2.0 * dt)
	df = 1.0 / float(num_seconds)

	print "dt = %.21f seconds" % dt
	print "n_bins = %d" % n_bins
	print "Nyquist freq =", nyquist_freq
	
	power_sum, sum_rate_whole, num_segments = read_and_use_segments(in_file, \
		n_bins, dt, test)
	
	print "\tTotal number of segments =", num_segments
	## Dividing sums by the number of segments to get an arithmetic average.
	power = power_sum / float(num_segments)
	assert int(len(power)) == n_bins, 'ERROR: Power should have length n_bins.'
	mean_rate_whole = sum_rate_whole / float(num_segments)
	print "Mean count rate over whole lightcurve =", mean_rate_whole

	freq, power, leahy_power, rms2_power, rms2_err_power = normalize(power, \
		n_bins, dt, num_seconds, num_segments, mean_rate_whole)
	
	rebinned_freq, rebinned_rms2_power, err_rebinned_power = \
		geometric_rebinning(freq, rms2_power, rms2_err_power, rebin_const)
	
	## Output, based on file extension
	if out_file[-4:].lower() == "fits":
		fits_output(out_file, rebinned_out_file, in_file, dt, n_bins, 
			nyquist_freq, num_segments, mean_rate_whole, freq, rms2_power, 
			rms2_err_power, leahy_power, rebin_const, rebinned_freq, 
			rebinned_rms2_power, err_rebinned_power)
	elif out_file[-3:].lower() == "dat":
		dat_output(out_file, rebinned_out_file, in_file, dt, n_bins, 
			nyquist_freq, num_segments, mean_rate_whole, freq, rms2_power, 
			rms2_err_power, leahy_power, rebin_const, rebinned_freq, 
			rebinned_rms2_power, err_rebinned_power)
	else:
		raise Exception("ERROR: Output file must be type .dat or .fits.")
		
## End of function 'main'
	


################################################################################
if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Makes a power spectrum from \
		an event-mode data file from RXTE.', epilog='For optional arguments,\
		default values are given in brackets at end of description.')
	parser.add_argument('infile', help='The full path of the input file with \
		RXTE event-mode data, with time in column 1 and rate in column 2. FITS \
		format must have extension .lc or .fits, otherwise assumes .dat \
		(ASCII/txt) format.')
	parser.add_argument('outfile', help='The full path of the .fits or .dat \
		file to write the frequency and power to.')
	parser.add_argument('rebinned_outfile', help='The full path of the .fits or\
		.dat file to write the re-binned frequency and power to.')
	parser.add_argument('-n', '--num_seconds', type=tools.type_power_of_two, \
		default=1, dest='num_seconds', help='Number of seconds in each Fourier \
		segment. Must be a power of 2. [1]')
	parser.add_argument('-c', '--rebin_const', type=tools.type_positive_float,\
		default=1.01, dest='rebin_const', help='Float constant by which we \
		re-bin the averaged power spectrum. [1.01]')
	parser.add_argument('-m', '--dt_mult', type=tools.type_power_of_two, \
		default=1, dest='dt_mult', help='Multiple of 1/8192 seconds for \
		timestep between bins. Must be a power of 2. [1]')
	parser.add_argument('--test', action='store_true', dest='test', help='If \
		present, only does a short test run.')
	args = parser.parse_args()
			
	main(args.infile, args.outfile, args.rebinned_outfile, args.num_seconds,
		args.rebin_const, args.dt_mult, args.test)

## End of program 'powerspec.py'
