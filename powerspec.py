import argparse
import numpy as np
from scipy import fftpack
from astropy.io import fits

import populate_lightcurve as lc
import tools

"""
		powerspec.py

Makes a power spectrum out of an event-mode data file from RXTE.

Arguments:
datafile - Name of FITS file with photon count rate data.
outfile - Name of file that the power spectrum will be written to.
rebinned_outfile - Name of file that the re-binned power spectrum will be 
	written to.
seconds - Number of seconds each segment of the light curve should be. Must be
	a power of 2.
rebin_const - Used to re-bin the data after the average power is computed, such 
	that bin_size[n+1] = bin_size[n] * rebin_const.
dt_mult - Multiple of 1/8192 seconds for the timestep between bins.
test - 1 if only computing one segment for testing, 0 if computing all segments.

Written in Python 2.7 by A.L. Stevens, A.L.Stevens@uva.nl, 2013-2014

The scientific modules imported above, as well as python 2.7, can be downloaded 
in the Anaconda package, https://store.continuum.io/cshop/anaconda/

populate_lightcurve is not on GitHub yet.
tools is available at https://github.com/abigailStev/whizzy_scripts

"""

################################################################################
def geometric_rebinning(freq, rms2_power_avg, rms2_err_power, rebin_const, 
						orig_length_of_list):
	"""
			geometric_rebinning
			
	Re-bins the noise-subtracted fractional rms^2 power spectrum in frequency 
	space by some re-binning constant (rebin_const>1). 

	Passed: rms2_power_avg - Fractional rms^2 power, averaged over all segments
				in the light curve.
			rms2_err_power - Error on the fractional rms^2 power.
			freq - Frequencies (in Hz) corresponding to the power spectrum.
			rebin_const - Constant >1 by which we want to re-bin the spectrum, 
				such that bin_size[n+1] = bin_size[n] * rebin_const.
			orig_length_of_list - Length of the original power spectrum (only 
				the positive frequencies).
	
	Returns: rebinned_freq - Frequencies of power spectrum re-binned according 
				to rebin_const.
			 rebinned_rms2_power - Power spectrum re-binned according to 
			 	rebin_const.
			 err_rebinned_power - Error on the re-binned rms power.
	
	"""
	pass
	## Initializing variables
	rebinned_rms2_power = []	 # List of re-binned fractional rms power
	rebinned_freq = []		 # List of re-binned frequencies
	err_rebinned_power = []	 # List of error in re-binned power
	real_index = 1.0		 # The unrounded next index in power_avg
	int_index = 1			 # The int of real_index, added to current_m every 
							 #  iteration
	current_m = 1			 # Current index in power_avg
	prev_m = 0				 # Previous index m
	bin_power = 0.0			 # The power of the current re-binned bin
	bin_freq = 0.0			 # The frequency of the current re-binned bin
	err_bin_power2 = 0.0	 # The error squared on 'bin_power'
	bin_range = 0.0			 # The range of un-binned bins covered by this 
							 #  re-binned bin
	
	## Looping through the length of the array power_avg, geometric bin by 
	## geometric bin, to compute the average power and frequency of that
	## geometric bin.
	## Equations for frequency, power, and error on power are from Adam Ingram's
	## PhD thesis
	while current_m < orig_length_of_list:
# 	while current_m < 100: # used for debugging
		## Initializing clean variables for each iteration of the while-loop
		bin_power = 0.0  # the averaged power at each index of rebinned_power
		err_bin_power2 = 0.0  # the square of the errors on powers in this bin
		bin_range = 0.0
		bin_freq = 0.0
		
		## Looping through the data points contained within one geometric bin
		for k in xrange(prev_m, current_m):
			bin_power += rms2_power_avg[k]
			err_bin_power2 += rms2_err_power[k] ** 2
			## End of for-loop
		
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
			bin_power = rms2_power_avg[prev_m]
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
def output(out_file, rebinned_out_file, in_file, dt, n_bins, nyquist_freq, 
	num_segments, mean_rate_whole, freq, rms2_power_avg, rms2_err_power, 
	leahy_power_avg, rebin_const, rebinned_freq, rebinned_rms2_power, 
	err_rebinned_power):
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
			rms2_power_avg - Fractional rms power, averaged over all 
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
	pass
	print "Standard output file: %s" % out_file
	
	with open(out_file, 'w') as out:
		out.write("#\t\tPower spectrum")
		out.write("\n# Data: %s" % in_file)
		out.write("\n# Time bin size = %.21f seconds" % dt)
		out.write("\n# Number of bins per segment = %d" % n_bins)
		out.write("\n# Number of segments per light curve = %d" % num_segments)
		out.write("\n# Duration of light curve used = %d seconds" \
			% (num_segments * n_bins * dt))
		out.write("\n# Mean count rate = %.8f, over whole light curve" % \
			mean_rate_whole)
		out.write("\n# Nyquist frequency = %.4f" % nyquist_freq)
		out.write("\n# ")
		out.write("\n# Column 1: Frequency in Hz (sample_frequency * 1.0/dt)")
		out.write("\n# Column 2: Fractional rms^2 normalized mean power")
		out.write("\n# Column 3: Fractional rms^2 normalized error on the mean \
			power")
		out.write("\n# Column 4: Leahy-normalized mean power")
		out.write("\n# ")
		for k in xrange(len(rms2_power_avg)):
			if freq[k] >= 0:
# 				out.write("\n%.8f\t%.8f\t%.8f" % (freq[k], rms2_power_avg[k], 
# 					rms2_err_power[k]))
				out.write("\n{0:.8f}\t{1:.8f}\t{2:.8f}\t{3:.8f}".format(freq[k], 
					rms2_power_avg[k], rms2_err_power[k], leahy_power_avg[k]))
				## End of if-statement
			## End of for-loop
		## End of with-block
	
	## Now outputting the re-binned data -- Need to do this separately since it
	## has a different number of data points from the regular power spectrum.
	
	print "Re-binned output file: %s" % rebinned_out_file
	
	with open(rebinned_out_file, 'w') as out:
		out.write("#\t\tPower spectrum")
		out.write("\n# Data: %s" % in_file)
		out.write("\n# Re-binned in frequency at (%.4f * prev bin size)" % \
			rebin_const)
		out.write("\n# Corresponding un-binned output file: %s" % out_file)
		out.write("\n# Original time bin size = %.21f seconds" % dt)
		out.write("\n# Duration of light curve used = %d seconds" \
			% (num_segments * n_bins * dt))
		out.write("\n# Mean count rate = %.8f, over whole light curve" % \
			mean_rate_whole)
		out.write("\n# ")
		out.write("\n# Column 1: Frequency in Hz")
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
	## End of function 'output'


################################################################################
def each_segment(rate):
	"""
			each_segment
	
	Computes the mean count rate of this segment, the FFT of the count rate 
	minus the mean, and the power spectrum of this segment.
	
	Passed: rate - Count rate for this segment of data.
	
	Returns: power_segment - Power spectra for this segment of data.
			 mean_rate - Mean count rate for this segment of data.
	
	""" 
	pass

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
	## End of function 'each_segment'


################################################################################
def fits_powerspec(in_file, n_bins, dt, print_iterator, test):
	"""
			fits_powerspec
	
	Opens the FITS file light curve, reads the count rate for a segment, calls 
	'each_segment' to create a power spectrum, adds power spectra over all 
	segments.
	
	Passed: in_file - Name of (FITS) input file/light curve.
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
	while j <= len(data.field(1)):  
	
		num_segments += 1

		## Extracts the second column of 'data' and assigns it to 'rate'. 
		rate = data[i:j].field(1)
		
		power_segment, mean_rate_segment = each_segment(rate)
		
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
	## End of function 'fits_powerspec'


################################################################################
def ascii_powerspec(in_file, n_bins, dt, print_iterator, test):
	"""
			ascii_powerspec
		
	Opens the (ASCII/txt/dat) event list, reads the count rate for a segment, 
	populates the light curve, calls 'each_segment' to create a power spectrum, 
	adds power spectra over all segments.
	
	Assumes that all event times have been corrected with TIMEZERO, and 
	GTI-filtered.
	
	Passed: in_file - Name of (ASCII/txt/dat) event list, still to be 
				'populated'.
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
	pass
	
	## Declaring clean variables to append to for every loop iteration.
	time = []
	energy = []
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
		print "\tERROR: Start time of data was not read in. Exiting."
		exit()
		
	end_time = start_time + (dt * n_bins)
	assert end_time > start_time
		
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
					time.append(current_time)
					energy.append(current_chan)

				if float(next_line[0]) > end_time:  
				# Triggered at the end of a segment
					if len(time) > 0:
						num_segments += 1
						rate_2d, rate_1d = lc.make_lightcurve(np.asarray(time), 
							np.asarray(energy), n_bins, dt, start_time)
						lightcurve = np.concatenate((lightcurve, rate_1d))
						
						power_segment, mean_rate_segment = each_segment(rate_1d)
						assert int(len(power_segment)) == n_bins
						
						power_sum += power_segment
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
						
						if (test is True) and (num_segments == 2):  # Testing
							np.savetxt('lightcurve.dat', lightcurve, fmt='%d')
							break
						## End of 'if there are counts in this segment'

					start_time += (n_bins * dt)
					end_time += (n_bins * dt)
					
					## End of 'if we're at the end of a segment'
				## End of 'if the line is not a comment'
			## End of for-loop 
		## End of with-block
	
	return power_sum, sum_rate_whole, num_segments
	## End of function 'ascii_powerspec'
	

################################################################################
def make_powerspec(in_file, n_bins, dt, test):
	"""
			make_powerspec
			
	Opens the file, reads in the count rate, calls 'each_segment' to create 
	power spectrum. Separated from main body like this so I can easily call it 
	in multi_powerspec.py. Split into 'fits_powerspec' and 'ascii_powerspec' 
	for easier readability.
	
	Passed: in_file - Name of input file with time in column 1 and rate in 
				column 2. FITS format must have extension .lc or .fits, 
				otherwise assumes .dat (ASCII/txt) format.
			n_bins - Number of integer bins per segment. Must be a power of 2.
			dt - Time step between bins, in seconds. Must be a power of 2.
			test - True if computing few segments, False if computing all.
	
	Returns: power_sum - Sum of power spectra for all segments of data from this
				input file.
			 sum_rate_whole - Sum of the mean count rate for all segments of 
			 	data from this input file.
			 num_segments - Number of segments of data from this input file.
			 
	"""
	pass
	
	assert tools.power_of_two(n_bins)
	assert tools.power_of_two(1.0 / dt)
	
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
	elif n_bins >= 2097152:
		print_iterator = int(1)	
	elif n_bins >= 1048576:
		print_iterator = int(2)
	else:
		print_iterator = int(5)
	
	## Looping through length of data file, segment by segment, to compute 
	## power for each data point in the segment
	print "Segments computed:"

	if using_fits:
		power_sum, sum_rate_whole, num_segments = fits_powerspec(in_file, 
			n_bins, dt, print_iterator, test)
	else:
		power_sum, sum_rate_whole, num_segments = ascii_powerspec(in_file, 
			n_bins, dt, print_iterator, test)
	
		## End of 'if/else file is fits format'
	
	return power_sum, sum_rate_whole, num_segments
	## End of function 'make_powerspec'
	
	
################################################################################
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
			num_seconds - Number of seconds each segment of the light curve 
				should be. Must be a power of 2.
			rebin_const - Used to re-bin the data geometrically after the 
				average power is computed, such that 
				bin_size[n+1] = bin_size[n] * rebin_const.
			dt_mult - Multiple of 1/8192 seconds for timestep between bins.
			test - True if computing 1 segment, False if computing all.
	
	Returns: nothing
	
	"""
	pass
	
	## Idiot checks, to ensure that our assumptions hold
	assert num_seconds > 0  # num_seconds must be a positive integer
	assert tools.power_of_two(num_seconds)
	assert rebin_const >= 1.0  # rebin_const must be a float greater than 1
	assert dt_mult >= 1

	t_res = 1.0 / 8192.0
 	#print "%.21f" %t_res
	dt = dt_mult * t_res
	#print dt
	#print "%.21f" % dt
	n_bins = num_seconds * int(1.0 / dt)
	nyquist_freq = 1.0 / (2.0 * dt)
	df = 1.0 / float(num_seconds)
	
	print "dt = %.21f seconds" % dt
	print "n_bins = %d" % n_bins
	print "Nyquist freq =", nyquist_freq
	
	power_sum, sum_rate_whole, num_segments = make_powerspec(in_file, n_bins, \
		dt, test)
	print "\tTotal number of segments =", num_segments

	## Dividing sums by the number of segments to get an arithmetic average.
	power_avg = power_sum / float(num_segments)
	assert int(len(power_avg)) == n_bins
	mean_rate_whole = sum_rate_whole / float(num_segments)
	print "Mean count rate over whole lightcurve =", mean_rate_whole

	## Computing the FFT sample frequencies (in Hz)
	freq = fftpack.fftfreq(n_bins, d=dt)

	## Ensuring that we're only using and saving the positive frequency values 
	##  (and associated power values)
	max_index = np.argmax(freq)
	freq = freq[0:max_index + 1]  # because it slices at end-1, and we want to 
								  # include max
	power_avg = power_avg[0:max_index + 1]
	
	## Computing the error on the mean power
	err_power = power_avg / np.sqrt(float(num_segments) * float(len(power_avg)))
	
	## Leahy normalization
	leahy_power_avg = 2.0 * power_avg * dt / float(n_bins) / mean_rate_whole
	print "Mean value of Leahy power =", np.mean(leahy_power_avg)  # ~2

	## Fractional rms^2 normalization with noise subtracted off
	rms2_power_avg = 2.0 * power_avg * dt / float(n_bins) / \
		(mean_rate_whole ** 2)
	rms2_power_avg -= (2.0 / mean_rate_whole)
	
	## Error on fractional rms^2 power (noise not subtracted off?)
	rms2_err_power = (2.0 * err_power * dt / float(n_bins) / \
		mean_rate_whole ** 2)	
	
	rebinned_freq, rebinned_rms2_power, err_rebinned_power = \
		geometric_rebinning(freq, rms2_power_avg, rms2_err_power, rebin_const,
		int(len(power_avg)))
		
	output(out_file, rebinned_out_file, in_file, dt, n_bins, nyquist_freq, 
		num_segments, mean_rate_whole, freq, rms2_power_avg, rms2_err_power, 
		leahy_power_avg, rebin_const, rebinned_freq, rebinned_rms2_power, 
		err_rebinned_power)
	
	## End of function 'main'
	


################################################################################
if __name__ == "__main__":
	"""
	Parsing cmd-line arguments and calling 'main'
	"""
	parser = argparse.ArgumentParser(description='Makes a power spectrum out \
		of an event-mode data file from RXTE.')
	parser.add_argument('-i', '--infile', required=True, dest='infile', \
		help='The full path of the input file with RXTE event-mode data, with \
		time in column 1 and rate in column 2. FITS format must have extension \
		.lc or .fits, otherwise assumes .dat (ASCII/txt) format.')
	parser.add_argument('-o', '--outfile', required=True, dest='outfile', \
		help='The full path of the (ASCII/txt) file to write the frequency and \
		power to.')
	parser.add_argument('-b', '--rebinned_outfile', required=True, 
		dest='rebinned_outfile', help='The full path of the (ASCII/txt) file \
		to write the re-binned frequency and power to.')
	parser.add_argument('-n', '--num_seconds', type=int, default=1, 
		dest='num_seconds', help='Number of seconds in each segment that the \
		light curve is broken up into. Must be a power of 2.')
	parser.add_argument('-c', '--rebin_const', type=float, default=1.01, 
		dest='rebin_const', help='Float constant by which we re-bin the \
		averaged power spectrum.')
	parser.add_argument('-m', '--dt_mult', type=int, default=1, dest='dt_mult', 
		help='Multiple of 1/8192 seconds for timestep between bins.')
	parser.add_argument('-t', '--test', type=int, default=0, choices=range(0,2), 
		dest='test', help='1 if only computing one segment for testing, 0 if \
		computing all segments.')
	args = parser.parse_args()
	
	test = False
	if args.test == 1: 
		test = True
		
	main(args.infile, args.outfile, args.rebinned_outfile, args.num_seconds,
		args.rebin_const, args.dt_mult, test)

## End of program 'powerspec.py'
