import argparse
import math
import numpy as np
from scipy import fftpack
from astropy.io import fits
import itertools

import populate_lightcurve as lc

"""
		powerspec.py

Makes a power spectrum out of event-mode data from RXTE.

Arguments:
datafile - str - Name of FITS file with photon count rate data.
outfile - str - Name of file that the power spectrum will be written to.
rebinned_outfile - str - Name of file that the re-binned power spectrum will be written 
	to.
seconds - int - Number of seconds each segment of the light curve should be. Must be a 
	power of 2.
rebin_const - float - Used to re-bin the data geometrically after the average power is 
	computed, such that bin_size[n+1] = bin_size[n] * rebin_const.
dt - int - Multiple of 1/8192 seconds for timestep between bins.

Written in Python 2.7 by A.L. Stevens, A.L.Stevens@uva.nl, 2013-2014

All scientific modules imported above, as well as python 2.7, can be downloaded in the 
Anaconda package, https://store.continuum.io/cshop/anaconda/
I don't think argparse came with Anaconda, but I don't remember installing anything 
special to get it.

populate_lightcurve is not on GitHub yet.

"""

#################################################
## Determines if 'num' is a power of 2 (>= 1)
##  Returns 'True' if 'num' is a power of 2, 
##  Returns 'False' if 'num' is not a power of 2
#################################################
def power_of_two(num):
	"""
			power_of_two
			
	Checks if an integer is a power of two.
	
	Passed: num - int - The number in question.
	
	Returns: bool - 'True' if 'num' is a power of two, 'False' if 'num' is not.
	
	"""
	n = int(num)
	x = 2
	
	if n == 1:
		return True
	else: 
		while x < n and x < 2147483648:
			x *= 2
		return n == x
	## End of function 'power_of_two'

########################################################################
## From https://docs.python.org/2/library/itertools.html#recipes
##  Used when reading lines in the file so I can peek at the next line.
########################################################################
def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = itertools.tee(iterable)
    next(b, None)
    return itertools.izip(a, b)
    ## End of function 'pairwise'


#############################################################################
## Re-bins the power spectrum in frequency space by some re-binning constant
#############################################################################
def geometric_rebinning(rms_power_avg, rms_err_power, freq, rebin_const, orig_length_of_list):
	""" 
			geometric_rebinning
			
	Re-bins the fractional rms power spectrum in frequency space by some re-binning
	constant (rebin_const>1). 
		
	Passed: rms_power_avg - list of floats - Fractional rms power, averaged over all
				segments in the light curve.
			rms_err_power - list of floats - Error on the rms power.
			freq - list of floats - Frequencies (in Hz) corresponding to the power 
				spectrum.
			rebin_const - float - Constant >1 by which we want to re-bin the spectrum,
				such that bin_size[n+1] = bin_size[n] * rebin_const.
			orig_length_of_list - int - Length of the original power spectrum (only the 
				positive frequencies).
	
	Returns: rebinned_freq - list of floats - Frequencies of power spectrum re-binned
			 	according to rebin_const.
			 rebinned_rms_power - list of floats - Power spectrum re-binned according
			 	to rebin_const.
			 err_rebinned_power - list of floats - Error on the re-binned rms power.
	
	"""
	
	## Initializing variables
	rebinned_rms_power = []	# List of re-binned fractional rms power
	rebinned_freq = []		# List of re-binned frequencies
	err_rebinned_power = []	# List of error in re-binned power
	real_index = 1.0		# The unrounded next index in power_avg
	int_index = 1			# The int of real_index, added to current_m every iteration
	current_m = 1			# Current index in power_avg
	prev_m = 0				# Previous index m
	bin_power = 0.0			# The power of the current re-binned bin
	bin_freq = 0.0			# The frequency of the current re-binned bin
	err_bin_power2 = 0.0	# The error squared on 'bin_power'
	bin_range = 0.0			# The range of un-binned bins covered by this re-binned bin
	
	
	## Looping through the length of the array power_avg, geometric bin by geometric bin, 
	##  to compute the average power and frequency of that geometric bin
	## Equations for frequency, power, and error on power are from Adam's thesis
	while current_m < orig_length_of_list:
# 	while current_m < 400: # used for debugging
		## Initializing clean variables for each iteration of the while-loop
		bin_power = 0.0 # the averaged power at each index of rebinned_power
		err_bin_power2 = 0.0 # the square of the errors on the powers in this bin
		bin_range = 0.0
		bin_freq = 0.0
		
		## Looping through the data points contained within one geometric bin
		for k in range (prev_m, current_m):
			## Adding power data points (tiny linear bins) within a geometric bin
			##  After the while-loop, this will be divided by total number of data points
			bin_power += rms_power_avg[k]
			
# 			if rms_power_avg[k] <= 0:
# 				print "k = ", k
# 				print "rms[k] = ", rms_power_avg[k]

			## Also computing error in bin power squared, for error computation later
			err_bin_power2 += rms_err_power[k] ** 2
			## End of for-loop
		
		## Determining the range of indices this specific geometric bin covers
		bin_range = abs(current_m - prev_m)
		## Dividing bin_power (currently just a sum of the data points) by the number 
		##  of points to get an arithmetic average
		bin_power /= float(bin_range)
		
		## Computing the mean frequency of a geometric bin
		bin_freq = ((freq[current_m] - freq[prev_m]) / bin_range) + freq[prev_m]

# 		if bin_power <= 0:
# 			print "Index:", current_m
# 			print "Bin range:", bin_range
# 			print "Bin freq:", bin_freq
# 			print "Bin power:", bin_power

		## If there's only one data point in the geometric bin, there's no need to take
		##  an average. This also prevents it from skipping the first data point.
		if bin_range == 1:
			bin_power = rms_power_avg[prev_m]
			bin_freq = freq[prev_m]
		
		## Appending values to arrays
		rebinned_rms_power.append(bin_power)
		rebinned_freq.append(bin_freq)
		## Computing error in geometric bin -- equation from Adam's thesis
		err_rebinned_power.append(math.sqrt(err_bin_power2) / float(bin_range))
		
		## Incrementing for the next iteration of the loop
		prev_m = current_m
		real_index *= rebin_const
		int_index = int(round(real_index))
		current_m += int_index
		
		## Since the for-loop goes from prev_m to current_m-1 (since that's how the range 
		##  function works) it's ok that we set prev_m = current_m here for the next 
		##  round. This will not cause any double-counting of bins or missing bins.
		
		## End of while-loop
		
	return rebinned_freq, rebinned_rms_power, err_rebinned_power
	## End of function 'geometric_rebinning'
	

########################################################################################
## Writes power spectrum and geometrically re-binned power spectrum to two output files
########################################################################################
def output(out_file, rebinned_out_file, in_file, dt, n_bins, nyquist_freq, num_segments, \
		   mean_rate_whole, freq, leahy_power_avg, rms_power_avg, rms_err_power, \
		   rebin_const, rebinned_freq, rebinned_rms_power, err_rebinned_power):
	""" 
			output
			
	Writes power spectrum and re-binned power spectrum to two output files.
		
	Passed: out_file - str - Name of output file for standard power spectrum.
			rebinned_out_file - str - Name of output file for geometrically re-binned 
				power spectrum.
			in_file - str - Full path of filename containing input data.
			dt - float - Size of time bin, in seconds (must be power of 2).
			n_bins - int - Number of time bins in a segment (must be power of 2).
			nyquist_freq - float - Nyquist frequency = 1/(2*dt)
			num_segments - int - Number of segments in the light curve.
			mean_rate_whole - float - Mean count rate of the light curve.
			freq - list of floats - Frequencies (in Hz) corresponding to the power 
				spectrum.
			leahy_power_avg - list of floats - Leahy-normalized power, averaged over all
				segments of the light curve.
			rms_power_avg - list of floats - Fractional rms power, averaged over all 
				segments of the light curve and all light curves.
			rms_err_power - list of floats - Error on avg fractional rms power.
			rebin_const - float - Constant >1 by which we want to re-bin the spectrum,
				such that bin_size[n+1] = bin_size[n] * rebin_const.
			rebinned_freq - list of floats - Frequencies of power spectrum re-binned
			 	according to rebin_const.
			rebinned_rms_power - list of floats - Power spectrum re-binned according to 
				rebin_const.
			err_rebinned_power - list of floats - Error on re-binned fractional rms 
				power.

	Returns: nothing
	
	"""
	
	print "Standard output file: %s" % out_file
	
	with open(out_file, 'w') as out:
		out.write("#\t\tPower spectrum")
		out.write("\n# Data: %s" % in_file)
		out.write("\n# Time bin size = %.13f seconds" % dt)
		out.write("\n# Number of bins per segment = %d" % n_bins)
		out.write("\n# Number of segments per light curve = %d" % num_segments)
		out.write("\n# Duration of light curve used = %d seconds" % (num_segments * n_bins * dt))
		out.write("\n# Mean count rate = %.8f, over whole light curve" % mean_rate_whole)
		out.write("\n# Nyquist frequency = %.4f" % nyquist_freq)
		out.write("\n# ")
		out.write("\n# Column 1: Frequency in Hz (sample_frequency * 1.0/dt)")
		out.write("\n# Column 2: Fractional rms normalized mean power")
		out.write("\n# Column 3: Fractional rms normalized error on the mean power")
		out.write("\n# Column 4: Leahy-normalized mean power")
		out.write("\n# ")
		for k in range(0, len(freq)):
			if freq[k] >= 0:
	# 			out.write("\n%.8f\t%.8f\t%.8f" % (freq[k], rms_power_avg[k], rms_err_power[k]))
				out.write("\n%.8f\t%.8f\t%.8f\t%.8f" % (freq[k], rms_power_avg[k], rms_err_power[k], leahy_power_avg[k]))
				## End of if-statement
			## End of for-loop
		## End of with-block
	
	## Now outputting the geometric-binned data -- Need to do this separately since it 
	##  has a different number of data points from the standard un-binned power spectrum.
	
	print "Re-binned output file: %s" % rebinned_out_file
	
	with open(rebinned_out_file, 'w') as out:
		out.write("#\t\tPower spectrum")
		out.write("\n# Data: %s" % in_file)
		out.write("\n# Geometrically re-binned in frequency at (%lf * previous bin size)" % rebin_const)
		out.write("\n# Corresponding un-binned output file: %s" % out_file)
		out.write("\n# Original time bin size = %.13f seconds" % dt)
		out.write("\n# Duration of light curve used = %d seconds" % (num_segments * n_bins * dt))
		out.write("\n# Mean count rate = %.8f, over whole light curve" % mean_rate_whole)
		out.write("\n# ")
		out.write("\n# Column 1: Frequency in Hz")
		out.write("\n# Column 2: Fractional rms normalized mean power")
		out.write("\n# Column 3: Error in fractional rms normalized mean power")
		out.write("\n# ")
		for k in range(0, len(rebinned_freq)):
			if rebinned_freq[k] >= 0:
				out.write("\n%.8f\t%.8f\t%.8f" % (rebinned_freq[k], rebinned_rms_power[k], err_rebinned_power[k]))
				## End of if-statement
			## End of for-loop
		## End of with-block
	## End of function 'output'


######################################################################
## Generating the power spectrum for each segment of the light curve.
######################################################################
def each_segment(rate):
	"""
			each_segment
	
	Computes the mean count rate of this segment, the FFT of the count rate minus the 
	mean, and the power spectrum of this segment.
	
	Passed: rate - list of floats - Count rate for this segment of data.
	
	Returns: power_segment - list of floats - Power spectra for this segment of data.
			 mean_rate_segment - float - Mean count rate for this segment of data.
	
	"""
	## Computing the mean count rate of the segment
	mean_rate_segment = np.mean(rate)

	## Subtracting the mean rate off each value of 'rate'
	##  This eliminates the spike at 0 Hz 
	rate_sub_mean = [x - mean_rate_segment for x in rate]

	## Taking the 1-dimensional FFT of the time-domain photon count rate
	##  Returns 'fft_data' as a complex-valued list
	##  Using the SciPy FFT algorithm, as it is faster than NumPy for large lists
	fft_data = fftpack.fft(rate_sub_mean)

	## Computing the power
	power_segment = np.absolute(fft_data)**2
	
	return power_segment, mean_rate_segment
	## End of function 'each_segment'


######################################################################################
## Opens the FITS file light curve, reads the count rate for a segment, calls 
##  'each_segment' to create a power spectrum, adds power spectra over all segments.
######################################################################################
def fits_powerspec(in_file, n_bins, dt, print_iterator):
	"""
			fits_powerspec
	
	Opens the FITS file light curve, reads the count rate for a segment, calls 
	'each_segment' to create a power spectrum, adds power spectra over all segments.
	
	Passed: in_file - str - Name of (FITS) input file/light curve.
			n_bins - int - Number of integer bins per segment. Must be a power of 2.
			dt - float - Time step between bins, in seconds. Must be a power of 2.
			print_iterator - Prints out which segment we're on every 'print_iterator' 
				segments, to keep track of progress.
	
	Returns: power_sum - list of floats - Sum of power spectra for all segments of data 
				from this input file.
			 sum_rate_whole - float - Sum of the mean count rate for all segments of data 
			 	from this input file.
			 num_segments - int - Number of segments of data from this input file.
			 
	"""
	
	## Opens the fits file using the Astropy library 'fits.open'.
	fits_hdu = fits.open(in_file)
	## Read the header information from the FITS file into 'header'.
	header = fits_hdu[1].header	
	## Get the data from the FITS file.
	## Usually, 1 is the photon count rate, but check the header info to be safe, 
	## Or use my whizzy_scripts/fits_info.py on GitHub.
	data = fits_hdu[1].data
	fits_hdu.close()
	
	sum_rate_whole = 0
	power_sum = [0 for x in range(n_bins)]
	num_segments = 1

	i = 0 # start of bin index to make segment of data for inner for-loop
	j = n_bins # end of bin index to make segment of data for inner for-loop
	
	assert dt == (data[1].field(0) - data[0].field(0))

	while j <= len(data.field(1)): ## so 'j' doesn't overstep the length of the file

		## Initializing clean variables for each iteration of the while-loop
		rate = []
		power_segment = []
		mean_rate_segment = 0

		## Printing out which segment we're on every x segments
		if num_segments % print_iterator == 0:
			print "\t", num_segments

		## Extracts the second column of 'data' and assigns it to 'rate'. 
		## Don't be a dumbass. Don't use a for-loop.
		rate = data[i:j].field(1)

		power_segment, mean_rate_segment = each_segment(rate)
		
# 			print "Mean rate of segment =", mean_rate_segment
# 			print "Power of segment =", power_segment

		power_sum = [a+b for a,b in zip(power_segment, power_sum)]	
		sum_rate_whole += mean_rate_segment
		
# 		if num_segments == 1: # For testing purposes
# 			break
		
		## Incrementing the counters and indices
		i = j
		j += n_bins
		num_segments += 1
		## Since the for-loop goes from i to j-1 (since that's how the range function 
		## works) it's ok that we set i=j here for the next round. This will not cause 
		## any double-counting of rows or missing rows.

		## End of while-loop
		
	return power_sum, sum_rate_whole, num_segments
	## End of function 'fits_powerspec'


#########################################################################################
## Opening the (ASCII/txt/dat) event list input file, reading in the count rate, 
##  populating the light curve from the event list, calling 'each_segment' to create the 
##  power spectrum.
#########################################################################################
def ascii_powerspec(in_file, n_bins, dt, print_iterator):
	"""
			ascii_powerspec
		
	Opens the (ASCII/txt/dat) event list, reads the count rate for a segment, populates 
	the light curve, calls 'each_segment' to create a power spectrum, adds power spectra
	over all segments.
	
	Passed: in_file - str - Name of (ASCII/txt/dat) input event list, still to be 
				'populated'.
			n_bins - int - Number of integer bins per segment. Must be a power of 2.
			dt - float - Desired time step between bins of populated light curve, in 
				seconds. Must be a power of 2.
			print_iterator - Prints out which segment we're on every 'print_iterator' 
				segments, to keep track of progress.
	
	Returns: power_sum - list of floats - Sum of power spectra for all segments of data 
				from this input file.
			 sum_rate_whole - float - Sum of the mean count rate for all segments of data 
			 	from this input file.
			 num_segments - int - Number of segments of data from this input file.
			 
	"""
	
	## Declaring clean variables to append to for every loop iteration.
	time = []
	energy = []
	sum_rate_whole = 0
	power_sum = [0 for x in range(n_bins)]
	num_segments = 1
	
	## Reading only the first line of data to get the start time of the file
	with open(in_file, 'r') as fo:
		for line in fo:
			if line[0].strip() != "#":
# 				print line
				line = line.strip().split()
				start_time = float(line[0])
				break
						
	end_time = start_time + (dt * n_bins)
# 	print "Start time of file is %.15f" % start_time
# 	print "End time of first seg is %.15f" % end_time
	assert end_time > start_time
		
	with open(in_file, 'r') as f:
		for line, next_line in pairwise(f):
			if line[0].strip() != "#": ## If the line is not a comment
				line = line.strip().split()
				next_line = next_line.strip().split()
				current_time = float(line[0])
				current_chan = int(line[1])
				
				if current_chan == 2: ## Only want PCU2 here
					time.append(current_time)
					energy.append(current_chan)
# 					print line

				if (float(next_line[0]) > end_time): ## Triggered at end of a segment
# 					print "Here1"
					print len(time)
					if len(time) > 0:
# 						print "\tHere2"
						print "Start time = %.13f" % start_time
						print "End time = %.13f" % end_time
						power_segment = []
						mean_rate_segment = 0
						rate_2d, rate_1d = lc.make_lightcurve(np.asarray(time), np.asarray(energy), n_bins, dt, start_time)
						power_segment, mean_rate_segment = each_segment(list(rate_1d))
# 						print "\t", mean_rate_segment
						assert len(power_segment) == n_bins

						power_sum = [a+b for a,b in zip(power_segment, power_sum)]	
						sum_rate_whole += mean_rate_segment
				
						## Printing out which segment we're on every x segments
						if num_segments % print_iterator == 0:
							print "\t", num_segments
	# 					if num_segments == 1:  ## For testing purposes only
	# 						break
						
						## Incrementing counters and loop control variables.
						time = []
						energy = []
						num_segments += 1
						## End of 'if there are counts in this segment'
						
					start_time += (n_bins * dt)
					end_time += (n_bins * dt)
					
					## End of 'if we're at the end of a segment'
				## End of 'if the line is not a comment'
			## End of for-loop 
		## End of with-block
		
	return power_sum, sum_rate_whole, num_segments-1
	## End of function 'ascii_powerspec'
	

#########################################################################################
## Opening input file, reading in the count rate, calling 'each_segment' to create power 
##  spectrum.
##  Separated from main body like this so I can easily call it in multi_powerspec.py
## 	Split into 'fits_powerspec' and 'ascii_powerspec' for easier readability.
#########################################################################################
def make_powerspec(in_file, n_bins, dt):
	"""
			make_powerspec
			
	Opens the file, reads in the count rate, calls 'each_segment' to create power 
	spectrum. Split into 'fits_powerspec' and 'ascii_powerspec' for easier readability.
	
	Passed: in_file - str - Name of input file with time in column 1 and rate in column 2. 
				FITS format must have extension .lc or .fits, otherwise assumes .dat 
				(ASCII/txt) format.
			n_bins - int - Number of integer bins per segment. Must be a power of 2.
			dt - float - Time step between bins, in seconds. Must be a power of 2.
	
	Returns: power_sum - list of floats - Sum of power spectra for all segments of data 
				from this input file.
			 sum_rate_whole - float - Sum of the mean count rate for all segments of data 
			 	from this input file.
			 num_segments - int - Number of segments of data from this input file.
			 
	"""
	
	assert power_of_two(n_bins)
	assert power_of_two(1.0/dt)
	
	len_fname = len(in_file)
		
	if (in_file[len_fname-3:len_fname].lower() == ".lc") or \
		(in_file[len_fname-5:len_fname].lower() == ".fits"):
		using_FITS = True
	else:
		using_FITS = False
	
	print "Input file: %s" % in_file
	print "Using a FITS file:", using_FITS
	
	if n_bins == 32768:
		print_iterator = int(5)
	elif n_bins < 32768:
		print_iterator = int(10)
	else:
		print_iterator = int(1)	
	
	## Looping through length of data file, segment by segment, to compute power for each 
	##  data point in the segment
	print "Segments computed:"

	if using_FITS:
		power_sum, sum_rate_whole, num_segments = fits_powerspec(in_file, n_bins, dt, print_iterator)
	else:
		power_sum, sum_rate_whole, num_segments = ascii_powerspec(in_file, n_bins, dt, print_iterator)
	
		## End of 'if/else file is FITS' 
	
	return power_sum, sum_rate_whole, num_segments
	## End of function 'make_powerspec'
	

###################################################################################
## Reads in a FITS file, takes FFT of data, makes power spectrum, writes to a file
###################################################################################
def main(in_file, out_file, rebinned_out_file, num_seconds, rebin_const, dt):
	""" 
			make_powerspec
			
	Reads in a file, takes FFT of segments of light curve data, computes power of 
	each segment, averages power over all segments, writes data to a file. 
	
	Passed: in_file - str - Name of input file with time in column 1 and rate in column 2. 
				FITS format must have extension .lc or .fits, otherwise assumes .dat 
				(ASCII/txt) format.
			out_file - str - Name of output file for standard power spectrum.
			rebinned_out_file - str - Name of output file for re-binned power spectrum.
			num_seconds - int - Number of seconds each segment of the light curve should 
				be. Must be a power of 2.
			rebin_const - float - Used to re-bin the data geometrically after the average 
				power is computed, such that bin_size[n+1] = bin_size[n] * rebin_const.
			dt - int - Multiple of 1/8192 seconds for timestep between bins.
	
	Returns: nothing
	
	"""
	pass
	
	## Idiot checks, to ensure that our assumptions hold
	assert num_seconds > 0 # num_seconds must be a positive integer
	assert power_of_two(num_seconds)
	assert rebin_const >= 1.0 # rebin_const must be a float greater than 1
# 	print dt
	t_res = 1.0/8192.0
# 	print "%.13f" %t_res
	dt = dt * t_res
# 	print "%.13f" % dt
	n_bins = num_seconds * int(1.0 / dt)
	nyquist_freq = 1.0 /(2.0 * dt)
	
	print "dt = %.13f seconds" % dt
	print "n_bins = %d" % n_bins
	print "Nyquist freq = ", nyquist_freq
	
	power_sum, sum_rate_whole, num_segments = make_powerspec(in_file, n_bins, dt)
	
	print "Total number of segments =", num_segments

	## Dividing sums by the number of segments to get an arithmetic average.
	power_avg = [x / float(num_segments) for x in power_sum]
	mean_rate_whole = sum_rate_whole / float(num_segments)
# 	print "Power spectrum computed"

# 	print "Power_avg =", power_avg
	print "Mean rate whole =", mean_rate_whole
	
	## Making integer time bins to plot against. Giving 't' as many data points as 'rate'
	t = np.arange(n_bins)
	
	## Computing the FFT sample frequencies
	##  Returns the frequency bins in cycles/unit starting at zero, given a window length 
	##  't'
	sample_freq = np.fft.fftfreq(t.shape[-1])

	## Changing sample frequency to Hz
	freq = [x * int(1.0 / dt) for x in sample_freq]
		
	## Ensuring that we're only using and saving the positive frequency values 
	##  (and associated power values)
	max_index = freq.index(max(freq))
	freq = freq[0:max_index+1] ## So that we don't cut off the last positive data point, 
							   ## since it cuts from 0 to end-1
	power_avg = power_avg[0:max_index+1]
	
	## Computing the error on the mean power
	err_power = [x / math.sqrt(float(num_segments)*len(power_avg)) for x in power_avg]
	
	## Leahy normalization
	leahy_power_avg = [(2.0 * x * dt / ((1.0/dt) * num_seconds) / mean_rate_whole) \
		for x in power_avg]
	print "Mean value of Leahy power =", np.mean(leahy_power_avg) ## Should be ~2	
	
	## Fractional rms normalization
	rms_power_avg = [(2.0 * x * dt / ((1.0/dt) * num_seconds) / mean_rate_whole**2) - \
		(2.0 / mean_rate_whole) for x in power_avg]
# 	other_rms_power_avg = [(x / mean_rate_whole) - (2.0 / mean_rate_whole) for x in \
# 		leahy_power_avg] ## This gives the same values as the above line
	"""
	## Troubleshooting -- not a problem!
	x = 0
	while x < len(rms_power_avg):
		if abs(rms_power_avg[x] - other_rms_power_avg[x]) > 10**-14:
			print "Issue:", x, rms_power_avg[x], other_rms_power_avg[x]
		x += 1
	print ""
	"""

	## Error on fractional rms power -- don't trust this equation (yet)
	rms_err_power = [(2.0 * x * dt / ((1.0/dt) * num_seconds) / mean_rate_whole**2) \
		for x in err_power]

	orig_length_of_list = len(power_avg)
	
	## Calling the above function for geometric re-binning
	rebinned_freq, rebinned_rms_power, err_rebinned_power = geometric_rebinning(rms_power_avg, rms_err_power, freq, rebin_const, orig_length_of_list)
		
	## Calling the above function for writing to output files
	output(out_file, rebinned_out_file, in_file, dt, n_bins, nyquist_freq, num_segments, \
		mean_rate_whole, freq, leahy_power_avg, rms_power_avg, rms_err_power, \
		rebin_const, rebinned_freq, rebinned_rms_power, err_rebinned_power)

	## End of function 'main'
	

#################################################
## Parsing cmd-line arguments and calling 'main'
#################################################

if __name__ == "__main__":
	
	parser = argparse.ArgumentParser()
	parser.add_argument('infile', \
		help="The full path of the input file with RXTE event mode data, with time in column 1 and rate in column 2. FITS format must have extension .lc or .fits, otherwise assumes .dat (ASCII/txt) format.")
	parser.add_argument('outfile', \
		help="The full path of the (ASCII/txt) file to write the frequency and power to.")
	parser.add_argument('rebinned_outfile', \
		help="The full path of the (ASCII/txt) file to write the geometrically re-binned frequency and power to.")
	parser.add_argument('seconds', type=int, \
		help="Duration of segments the light curve is broken up into, in seconds. Must be an integer power of 2.")
	parser.add_argument('rebin_const', type=float, \
		help="Float constant by which we geometrically re-bin the averaged power spectrum.")
	parser.add_argument('dt', type=int, help="Multiple of 1/8192 seconds for timestep between bins.")
	args = parser.parse_args()

	main(args.infile, args.outfile, args.rebinned_outfile, args.seconds, args.rebin_const, args.dt)

## End of program 'powerspec.py'
