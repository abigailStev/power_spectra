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
datafile - Name of FITS file with photon count rate data.
outfile - Name of file that the power spectrum will be written to.
rebinned_outfile - Name of file that the re-binned power spectrum will be written to.
seconds - Number of seconds each segment of the light curve should be. Must be a power of 
	2.
rebin_const - Used to re-bin the data geometrically after the average power is computed, 
	such that bin_size[n+1] = bin_size[n] * rebin_const.
dt_mult - Multiple of 1/8192 seconds for the timestep between bins.

Written in Python 2.7 by A.L. Stevens, A.L.Stevens@uva.nl, 2013-2014

The scientific modules imported above, as well as python 2.7, can be downloaded in the 
Anaconda package, https://store.continuum.io/cshop/anaconda/

populate_lightcurve is not on GitHub yet.

"""

#########################################################################################
def power_of_two(num):
	"""
			power_of_two
			
	Checks if 'num' is a power of 2 (1 <= num < 2147483648 )
	
	Passed: num - The number in question.
	
	Returns: bool - 'True' if 'num' is a power of two, 'False' if 'num' is not.
	
	"""
	pass
	n = int(num)
	x = 2
	
	if n == 1:
		return True
	else: 
		while x < n and x < 2147483648:
			x *= 2
		return n == x
	## End of function 'power_of_two'


#########################################################################################
def pairwise(iterable):
	"""
	s -> (s0,s1), (s1,s2), (s2, s3), ...
	https://docs.python.org/2/library/itertools.html#recipes
	Used when reading lines in the file so I can peek at the next line.
	"""
	pass
	a, b = itertools.tee(iterable)
	next(b, None)
	return itertools.izip(a, b)
	## End of function 'pairwise'


#########################################################################################
def geometric_rebinning(rms_power_avg, rms_err_power, freq, rebin_const,
                        orig_length_of_list):
	"""
			geometric_rebinning
			
	Re-bins the fractional rms power spectrum in frequency space by some re-binning
	constant (rebin_const>1). 

	Passed: rms_power_avg - Fractional rms power, averaged over all segments in the light
				curve.
			rms_err_power - Error on the rms power.
			freq - Frequencies (in Hz) corresponding to the power spectrum.
			rebin_const - Constant >1 by which we want to re-bin the spectrum, such that 
				bin_size[n+1] = bin_size[n] * rebin_const.
			orig_length_of_list - Length of the original power spectrum (only the 
				positive frequencies).
	
	Returns: rebinned_freq - Frequencies of power spectrum re-binned according to 
				rebin_const.
			 rebinned_rms_power - Power spectrum re-binned according to rebin_const.
			 err_rebinned_power - Error on the re-binned rms power.
	
	"""
	pass
	## Initializing variables
	rebinned_rms_power = []	 # List of re-binned fractional rms power
	rebinned_freq = []		 # List of re-binned frequencies
	err_rebinned_power = []	 # List of error in re-binned power
	real_index = 1.0		 # The unrounded next index in power_avg
	int_index = 1			 # The int of real_index, added to current_m every iteration
	current_m = 1			 # Current index in power_avg
	prev_m = 0				 # Previous index m
	bin_power = 0.0			 # The power of the current re-binned bin
	bin_freq = 0.0			 # The frequency of the current re-binned bin
	err_bin_power2 = 0.0	 # The error squared on 'bin_power'
	bin_range = 0.0			 # The range of un-binned bins covered by this re-binned bin
	
	## Looping through the length of the array power_avg, geometric bin by geometric bin, 
	##  to compute the average power and frequency of that geometric bin
	## Equations for frequency, power, and error on power are from Adam's thesis
	while current_m < orig_length_of_list:
# 	while current_m < 400: # used for debugging
		## Initializing clean variables for each iteration of the while-loop
		bin_power = 0.0  # the averaged power at each index of rebinned_power
		err_bin_power2 = 0.0  # the square of the errors on the powers in this bin
		bin_range = 0.0
		bin_freq = 0.0
		
		## Looping through the data points contained within one geometric bin
		for k in xrange(prev_m, current_m):
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
def output(out_file, rebinned_out_file, in_file, dt, n_bins, nyquist_freq, num_segments,
           mean_rate_whole, freq, rms_power_avg, rms_err_power, leahy_power_avg, 
           rebin_const, rebinned_freq, rebinned_rms_power, err_rebinned_power):
	""" 
			output
			
	Writes power spectrum and re-binned power spectrum to two output files.
		
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
			leahy_power_avg - Leahy-normalized power, averaged over all
				segments of the light curve.
			rms_power_avg - Fractional rms power, averaged over all 
				segments of the light curve and all light curves.
			rms_err_power - Error on avg fractional rms power.
			rebin_const - Constant >1 by which we want to re-bin the spectrum,
				such that bin_size[n+1] = bin_size[n] * rebin_const.
			rebinned_freq - Frequencies of power spectrum re-binned
				according to rebin_const.
			rebinned_rms_power - Power spectrum re-binned according to 
				rebin_const.
			err_rebinned_power - Error on re-binned fractional rms 
				power.

	Returns: nothing
	
	"""
	pass
	print "Standard output file: %s" % out_file
	
	with open(out_file, 'w') as out:
		out.write("#\t\tPower spectrum")
		out.write("\n# Data: %s" % in_file)
		out.write("\n# Time bin size = %.13f seconds" % dt)
		out.write("\n# Number of bins per segment = %d" % n_bins)
		out.write("\n# Number of segments per light curve = %d" % num_segments)
		out.write("\n# Duration of light curve used = %d seconds" \
				  % (num_segments * n_bins * dt))
		out.write("\n# Mean count rate = %.8f, over whole light curve" % mean_rate_whole)
		out.write("\n# Nyquist frequency = %.4f" % nyquist_freq)
		out.write("\n# ")
		out.write("\n# Column 1: Frequency in Hz (sample_frequency * 1.0/dt)")
		out.write("\n# Column 2: Fractional rms normalized mean power")
		out.write("\n# Column 3: Fractional rms normalized error on the mean power")
		out.write("\n# Column 4: Leahy-normalized mean power")
		out.write("\n# ")
		for k in range(0, len(rms_power_avg)):
			if freq[k] >= 0:
# 				out.write("\n%.8f\t%.8f\t%.8f" % (freq[k], rms_power_avg[k], rms_err_power[k]))
				out.write("\n{0:.8f}\t{1:.8f}\t{2:.8f}\t{3:.8f}".format(freq[k], rms_power_avg[k],
                                                                        rms_err_power[k], leahy_power_avg[k]))
				## End of if-statement
			## End of for-loop
		## End of with-block
	
	## Now outputting the geometric-binned data -- Need to do this separately since it 
	##  has a different number of data points from the standard un-binned power spectrum.
	
	print "Re-binned output file: %s" % rebinned_out_file
	
	with open(rebinned_out_file, 'w') as out:
		out.write("#\t\tPower spectrum")
		out.write("\n# Data: %s" % in_file)
		out.write("\n# Geometrically re-binned in frequency at (%lf * previous bin size)\
				  " % rebin_const)
		out.write("\n# Corresponding un-binned output file: %s" % out_file)
		out.write("\n# Original time bin size = %.13f seconds" % dt)
		out.write("\n# Duration of light curve used = %d seconds" \
				  % (num_segments * n_bins * dt))
		out.write("\n# Mean count rate = %.8f, over whole light curve" % mean_rate_whole)
		out.write("\n# ")
		out.write("\n# Column 1: Frequency in Hz")
		out.write("\n# Column 2: Fractional rms normalized mean power")
		out.write("\n# Column 3: Error in fractional rms normalized mean power")
		out.write("\n# ")
		for k in range(0, len(rebinned_rms_power)):
			if rebinned_freq[k] >= 0:
				out.write("\n{0:.8f}\t{1:.8f}\t{2:.8f}".format(rebinned_freq[k],
                                                               rebinned_rms_power[k], err_rebinned_power[k]))
				## End of if-statement
			## End of for-loop
		## End of with-block
	## End of function 'output'


#########################################################################################
def each_segment(rate):
	"""
			each_segment
	
	Computes the mean count rate of this segment, the FFT of the count rate minus the 
	mean, and the power spectrum of this segment.
	
	Passed: rate - Count rate for this segment of data.
	
	Returns: power_segment - Power spectra for this segment of data.
			 mean_rate - Mean count rate for this segment of data.
	
	""" 
	pass
	print "Non-zero rate indices:", np.where(rate != 0)
	print "Num of non-zero rates:", np.shape(np.where(rate != 0))
	## Computing the mean count rate of the segment
	mean_rate = np.mean(rate)
# 	print "Shape of rate:", np.shape(rate)
	## Subtracting the mean rate off each value of 'rate'
	##  This eliminates the spike at 0 Hz 
	rate_sub_mean = rate - mean_rate
# 	print "Rate:", rate[0:4]
# 	print "Mean:", mean_rate
# 	print "Mean-subtracted rate:", rate_sub_mean[0:4]
# 	print "Shape of rate_sub_mean:", np.shape(rate_sub_mean)
	## good until here!
	
	## Taking the 1-dimensional FFT of the time-domain photon count rate
	##  Using the SciPy FFT algorithm, as it is faster than NumPy for large lists
	fft_data = fftpack.fft(rate_sub_mean)

	## Computing the power
	power_segment = np.absolute(fft_data) ** 2
	
	return power_segment, mean_rate
	## End of function 'each_segment'


#########################################################################################
def fits_powerspec(in_file, n_bins, dt, print_iterator, short_run):
	"""
			fits_powerspec
	
	Opens the FITS file light curve, reads the count rate for a segment, calls 
	'each_segment' to create a power spectrum, adds power spectra over all segments.
	
	Passed: in_file - Name of (FITS) input file/light curve.
			n_bins - Number of integer bins per segment. Must be a power of 2.
			dt - Time step between bins, in seconds. Must be a power of 2.
			print_iterator - Prints out which segment we're on every 'print_iterator' 
				segments, to keep track of progress.
			short_run - True if computing one segment, False if computing all.
	
	Returns: power_sum - Sum of power spectra for all segments of data from this input 
				file.
			 sum_rate_whole - Sum of the mean count rate for all segments of data from 
				this input file.
			 num_segments - Number of segments of data from this input file.
			 
	"""
	pass
	
	fits_hdu = fits.open(in_file)
	header = fits_hdu[1].header	
	data = fits_hdu[1].data
	fits_hdu.close()
	
	sum_rate_whole = 0
	power_sum = np.zeros(n_bins, dtype=np.float64)
	num_segments = 0

	i = 0  # start of bin index to make segment of data for inner for-loop
	j = n_bins  # end of bin index to make segment of data for inner for-loop

	assert dt == (data[1].field(0) - data[0].field(0))
	print "Length of FITS file:", len(data.field(1))
	while j <= len(data.field(1)):  # so 'j' doesn't overstep the length of the file
		
		num_segments += 1
		print "\tNew Segment; i=%d, j=%d" % (i,j)

		## Extracts the second column of 'data' and assigns it to 'rate'. 
		time = data[i:j].field(0)
		rate = data[i:j].field(1)

		power_segment, mean_rate_segment = each_segment(rate)
		print "Start time of segment: %.13f" % data[i].field(0)
		print "End time of segment: %.13f" % data[j-1].field(0)
		print "4 times: %.13f, %.13f, %.13f, %.13f" \
			% (time[50], time[51], time[52], time[53])
		print "rate:", rate[50:54]
# 		print "Mean rate of segment =", mean_rate_segment
# 		print "Power of segment =", power_segment
# 		print "Power segment: ", power_segment[0:4]
# 		print "Power sum before: ", power_sum[0:4]
		power_sum += power_segment
# 		print "Power sum after: ", power_sum[0:4]
# 		print sum_rate_whole
# 		print mean_rate_segment
		sum_rate_whole += mean_rate_segment
# 		print sum_rate_whole
		if num_segments % print_iterator == 0:
			print "\t", num_segments
			
		if (short_run == True) and (num_segments == 1):  # For testing
			break
			
		rate = None
		power_segment = None
		mean_rate_segment = None
		## Incrementing the counters and indices
		i = j
		j += n_bins
		## Since the for-loop goes from i to j-1 (since that's how the range function 
		## works) it's ok that we set i=j here for the next round. This will not cause 
		## any double-counting of rows or missing rows.
       
		## End of while-loop
		
	return power_sum, sum_rate_whole, num_segments
	## End of function 'fits_powerspec'


#########################################################################################
def ascii_powerspec(in_file, n_bins, dt, print_iterator, short_run):
	"""
			ascii_powerspec
		
	Opens the (ASCII/txt/dat) event list, reads the count rate for a segment, populates 
	the light curve, calls 'each_segment' to create a power spectrum, adds power spectra
	over all segments.
	
	Passed: in_file - Name of (ASCII/txt/dat) event list, still to be 'populated'.
			n_bins - Number of integer bins per segment. Must be a power of 2.
			dt - Desired time step between bins of populated light curve, in 
				seconds. Must be a power of 2.
			print_iterator - Prints out which segment we're on every 'print_iterator' 
				segments, to keep track of progress.
			short_run - True if computing a few segments, False if computing all.
	
	Returns: power_sum - Sum of power spectra for all segments of data from this input 
				file.
			 sum_rate_whole - Sum of the mean count rate for all segments of data from 
				this input file.
			 num_segments - Number of segments of data from this input file.
			 
	"""
	pass
	
	## Declaring clean variables to append to for every loop iteration.
	time = []
	energy = []
	sum_rate_whole = 0
	power_sum = np.zeros(n_bins, dtype=np.float64)
	num_segments = 0
	start_time = -99
	
	lightcurve = np.asarray([])
	
	## Reading only the first line of data to get the start time of the file
	with open(in_file, 'r') as fo:
		for line in fo:
			if line[0].strip() != "#":
# 				print line
				line = line.strip().split()
				start_time = float(line[0])
				break
	if start_time is -99:
		print "\tERROR: Start time of data was not read in. Exiting."
		exit()
# 	start_time = 277473713.3784303665161
	end_time = start_time + (dt * n_bins)
# 	print "Start time of file is %.15f" % start_time
# 	print "End time of first seg is %.15f" % end_time
	assert end_time > start_time
		
	with open(in_file, 'r') as f:
		for line, next_line in pairwise(f):
			if line[0].strip() != "#" and line[0].strip().split() >= start_time:  # If the line is not a comment
# 				print line.strip().split()[0]
				line = line.strip().split()
				next_line = next_line.strip().split()
				current_time = float(line[0])
				current_chan = int(line[1])
				
				if current_chan == 2:  # Only want PCU2 here
					time.append(current_time)
					energy.append(current_chan)
# 					print line

				if float(next_line[0]) > end_time:  # Triggered at end of a segment
# 					print next_line[0], end_time
# 					print "Here1"
# 					print "Length of time array at end of segment:", len(time)
					if len(time) > 0:
						num_segments += 1
						print "\tNew Segment"
						print "Start time of segment: %.13f" % start_time
						print "End time of segment: %.13f" % end_time
						rate_2d, rate_1d = lc.make_lightcurve(np.asarray(time), np.asarray(energy), n_bins, dt, start_time)
# 						print "Shape of rate_1d:", np.shape(rate_1d)
						lightcurve = np.concatenate((lightcurve, rate_1d))
# 						print "dt = ", dt
# 						print np.where(rate_1d != 0)
# 						print rate_1d[np.where(rate_1d != 0)]
# 						print len(rate_1d[np.where(rate_1d != 0)])
						power_segment, mean_rate_segment = each_segment(rate_1d)
# 						print "\t", mean_rate_segment
						assert int(len(power_segment)) == n_bins
# 						print "Shape of power segment:", np.shape(power_segment)
# 						print "Power segment: ", power_segment[0:4]
# 						print "Power sum before: ", power_sum[0:4]
						power_sum += power_segment
# 						print "Power sum after: ", power_sum[0:4]
# 						print "Shape of power sum:", np.shape(power_sum)
						sum_rate_whole += mean_rate_segment
						## Printing out which segment we're on every x segments
						if num_segments % print_iterator == 0:
							print "\t", num_segments

						## Clearing variables from memory
						time = None
						energy = None
						power_segment = None
						mean_rate_segment = None
						rate_2d = None
						rate_1d = None
						time = []
						energy = []
						
						if (short_run == True) and (num_segments == 1):  # For testing
							break
						## End of 'if there are counts in this segment'

					start_time += (n_bins * dt)
					end_time += (n_bins * dt)
					
					## End of 'if we're at the end of a segment'
				## End of 'if the line is not a comment'
			## End of for-loop 
		## End of with-block
	
	np.savetxt('lightcurve.dat', lightcurve, fmt='%d')
		
	return power_sum, sum_rate_whole, num_segments
	## End of function 'ascii_powerspec'
	

#########################################################################################
def make_powerspec(in_file, n_bins, dt, short_run):
	"""
			make_powerspec
			
	Opens the file, reads in the count rate, calls 'each_segment' to create power 
	spectrum. Separated from main body like this so I can easily call it in 
	multi_powerspec.py. Split into 'fits_powerspec' and 'ascii_powerspec' for easier 
	readability.
	
	Passed: in_file - Name of input file with time in column 1 and rate in column 2. 
				FITS format must have extension .lc or .fits, otherwise assumes .dat 
				(ASCII/txt) format.
			n_bins - Number of integer bins per segment. Must be a power of 2.
			dt - Time step between bins, in seconds. Must be a power of 2.
			short_run - True if computing few segments, False if computing all.
	
	Returns: power_sum - Sum of power spectra for all segments of data from this input 
				file.
			 sum_rate_whole - Sum of the mean count rate for all segments of data from 
				this input file.
			 num_segments - Number of segments of data from this input file.
			 
	"""
	pass
	
	assert power_of_two(n_bins)
	assert power_of_two(1.0 / dt)
	
	len_fname = len(in_file)
		
	if (in_file[len_fname - 3:len_fname].lower() == ".lc") or \
		(in_file[len_fname - 5:len_fname].lower() == ".fits"):
		using_fits = True
	else:
		using_fits = False
	
	print "Input file: %s" % in_file
	print "Using a FITS file:", using_fits
	
	if n_bins == 32768:
		print_iterator = int(10)
	elif n_bins < 32768:
		print_iterator = int(20)
	else:
		print_iterator = int(5)	
	
	## Looping through length of data file, segment by segment, to compute power for each 
	##  data point in the segment
	print "Segments computed:"

	if using_fits:
		power_sum, sum_rate_whole, num_segments = fits_powerspec(in_file, n_bins, dt,
                                                                 print_iterator, 
                                                                 short_run)
	else:
		power_sum, sum_rate_whole, num_segments = ascii_powerspec(in_file, n_bins, dt, print_iterator, short_run)
	
		## End of 'if/else file is fits format'
	
	return power_sum, sum_rate_whole, num_segments
	## End of function 'make_powerspec'
	
	
#########################################################################################
def main(in_file, out_file, rebinned_out_file, num_seconds, rebin_const, dt_mult, short_run):
	""" 
			main
	
	Passed: in_file - Name of input file with time in column 1 and rate in column 2. 
				FITS format must have extension .lc or .fits, otherwise assumes .dat 
				(ASCII/txt) format.
			out_file - Name of output file for standard power spectrum.
			rebinned_out_file - Name of output file for re-binned power spectrum.
			num_seconds - Number of seconds each segment of the light curve should be. 
				Must be a power of 2.
			rebin_const - Used to re-bin the data geometrically after the average power 
				is computed, such that bin_size[n+1] = bin_size[n] * rebin_const.
			dt_mult - Multiple of 1/8192 seconds for timestep between bins.
			short_run - True if computing 1 segment, False if computing all.
	
	Returns: nothing
	
	"""
	pass
	
	## Idiot checks, to ensure that our assumptions hold
	assert num_seconds > 0  # num_seconds must be a positive integer
	assert power_of_two(num_seconds)
	assert rebin_const >= 1.0  # rebin_const must be a float greater than 1

	t_res = 1.0 / 8192.0
 	#print "%.13f" %t_res
	dt = dt_mult * t_res
	#print dt
	#print "%.13f" % dt
	n_bins = num_seconds * int(1.0 / dt)
	nyquist_freq = 1.0 / (2.0 * dt)
	
	print "dt = %.13f seconds" % dt
	print "n_bins = %d" % n_bins
	print "Nyquist freq = ", nyquist_freq
	
	power_sum, sum_rate_whole, num_segments = make_powerspec(in_file, n_bins, dt, short_run)
	print "\tTotal number of segments =", num_segments

	## Dividing sums by the number of segments to get an arithmetic average.
	power_avg = power_sum / float(num_segments)
	assert int(len(power_avg)) == n_bins
	print "Power avg:", power_avg[0:4]
	mean_rate_whole = sum_rate_whole / float(num_segments)

	print "Mean rate whole =", mean_rate_whole

	## Computing the FFT sample frequencies (in Hz)
	t = np.arange(n_bins)
# 	print "t.shape[-1] =", t.shape[-1]
	freq = fftpack.fftfreq(n_bins, d=dt)

	## Ensuring that we're only using and saving the positive frequency values 
	##  (and associated power values)
	max_index = np.argmax(freq)
	print "Index of max freq:", max_index
# 	print "Type of 'max_index':", type(max_index)
	freq = freq[0:max_index + 1]  # because it slices at end-1, and we want to include max
	power_avg = power_avg[0:max_index + 1]
	
	## Computing the error on the mean power
	err_power = power_avg / math.sqrt(float(num_segments) * float(len(power_avg)))
	
	## Leahy normalization
	leahy_power_avg = (2.0 * power_avg * dt / ((1.0 / dt) * num_seconds) / mean_rate_whole)
	print "Mean value of Leahy power =", np.mean(leahy_power_avg)  # Should be ~2
	
	## Fractional rms normalization
	rms_power_avg = (2.0 * power_avg * dt / ((1.0 / dt) * num_seconds) / \
		mean_rate_whole ** 2) - (2.0 / mean_rate_whole)
	print "Mean value of rms power =", np.mean(rms_power_avg)
	
	## Error on fractional rms power -- don't trust this equation (yet)
	rms_err_power = (2.0 * err_power * dt / ((1.0 / dt) * num_seconds) / mean_rate_whole ** 2)
	print "rms power avg:", rms_power_avg[0:4]
	
	
	rebinned_freq, rebinned_rms_power, err_rebinned_power = \
		geometric_rebinning(rms_power_avg, rms_err_power, freq, rebin_const,
		                    int(len(power_avg)))
		
	output(out_file, rebinned_out_file, in_file, dt, n_bins, nyquist_freq, num_segments,
           mean_rate_whole, freq, rms_power_avg, rms_err_power, leahy_power_avg, 
           rebin_const, rebinned_freq, rebinned_rms_power, err_rebinned_power)
	
	freq401 = np.where(freq == 401)
	print "Bin at 401 Hz:", freq401
	print "RMS power at 401 Hz:", rms_power_avg[freq401]
	maxpow = np.argmax(rms_power_avg)
	print "Bin of max power:", maxpow
	print "Freq of max power:", freq[maxpow],"Hz"
	print "Max rms power:", rms_power_avg[maxpow]
	
	
	## End of function 'main'
	


#########################################################################################
if __name__ == "__main__":
	"""
	Parsing cmd-line arguments and calling 'main'
	"""
	parser = argparse.ArgumentParser()
	parser.add_argument('infile',
		                help="The full path of the input file with RXTE event mode data, \
		                with time in column 1 and rate in column 2. FITS format must \
		                have extension .lc or .fits, otherwise assumes .dat (ASCII/txt) \
		                format.")
	parser.add_argument('outfile',
		                help="The full path of the (ASCII/txt) file to write the \
		                frequency and power to.")
	parser.add_argument('rebinned_outfile',
		                help="The full path of the (ASCII/txt) file to write the \
		                geometrically re-binned \
		                frequency and power to.")
	parser.add_argument('seconds', type=int,
		                help="Duration of segments the light curve is broken up into, \
		                in seconds. Must be an integer power of 2.")
	parser.add_argument('rebin_const', type=float,
                        help="Float constant by which we geometrically re-bin the \
                        averaged power spectrum.")
	parser.add_argument('dt_mult', type=int,
		                help="Multiple of 1/8192 seconds for timestep between bins.")
	parser.add_argument('short_run', type=int,
		                help="1 if only computing one segment for testing, 0 if \
		                computing all segments.")
	args = parser.parse_args()
	
	assert args.short_run == 1 or args.short_run == 0
	if args.short_run == 1: 
		short_run = True
	else:
		short_run = False
		
	print "infile =", args.infile
	print "outfile =", args.outfile
	print "rebinned_outfile =", args.rebinned_outfile
	print "seconds =", args.seconds
	print "rebin_const =", args.rebin_const
	print "dt_mult =", args.dt_mult
	print "short run?", short_run
		
	main(args.infile, args.outfile, args.rebinned_outfile, args.seconds,
         args.rebin_const, args.dt_mult, short_run)

## End of program 'powerspec.py'
