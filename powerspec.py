import argparse
import numpy as np
from scipy import fftpack
from astropy.io import fits
from datetime import datetime
import os
import tools  # https://github.com/abigailStev/whizzy_scripts

__author__ = "Abigail Stevens"
__author_email__ = "A.L.Stevens@uva.nl"
__year__ = "2013-2015"
__description__ = "Makes a power spectrum averaged over segments from an RXTE \
event-mode data file."

"""
		powerspec.py

Written in Python 2.7.

The scientific modules imported above, as well as python 2.7, can be downloaded 
in the Anaconda package, https://store.continuum.io/cshop/anaconda/

"""

################################################################################
def dat_out(out_file, in_file, dt, n_bins, nyquist_freq, num_segments, \
	mean_rate_whole, freq, rms2_power, rms2_err_power, leahy_power):
	""" 
			dat_out
			
	Writes power spectrum to an ASCII output file.
			
	"""
	print "Output file: %s" % out_file
	
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

## End of function 'dat_out'


################################################################################
def fits_out(out_file, in_file, dt, n_bins, nyquist_freq, num_segments, \
	mean_rate_whole, freq, rms2_power, rms2_err_power, leahy_power):
	"""
				fits_out
			
	Writes power spectrum to a FITS file.
	
	"""
	print "Output file: %s" % out_file

	## Making header for standard power spectrum
	prihdr = fits.Header()
	prihdr.set('TYPE', "Power spectrum")
	prihdr.set('DATE', str(datetime.now()), "YYYY-MM-DD localtime")
	prihdr.set('EVTLIST', in_file)
	prihdr.set('DT', dt, "seconds")
	prihdr.set('N_BINS', n_bins, "time bins per segment")
	prihdr.set('SEGMENTS', num_segments, "segments in the whole light curve")
	prihdr.set('EXPOSURE', num_segments * n_bins * dt, "seconds, of light curve")
	prihdr.set('MEANRATE', mean_rate_whole, "counts/second")
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

## End of function 'fits_out'


###############################################################################
def normalize(power, n_bins, dt, num_seconds, num_segments, mean_rate, noisy):
	"""
			normalize
	
	Generates the Fourier frequencies, removes negative frequencies, normalizes
	the power by Leahy and fractional rms^2 normalizations, and computes the 
	error on the fractional rms^2 power.
	
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
	if noisy:
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
	print "Signal variance:", signal_variance
	rms = np.sqrt(signal_variance)  # should be a few % in frac rms units
	print "RMS of signal:", rms
	
	return freq, power, leahy_power, rms2_power, rms2_err_power
	## End of function 'normalize'
	

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
def extracted_in(in_file, n_bins, dt, print_iterator, test):
	"""
			extracted_powerspec
	
	Opens the FITS file light curve (as created in seextrct), reads the count 
	rate for a segment, calls 'make_ps' to create a power spectrum, adds power 
	spectra over all segments.
		 
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

	assert dt == (data[1].field(0) - data[0].field(0)), \
		'ERROR: dt must be the same resolution as the extracted FITS data.'
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
## End of function 'extracted_in'


################################################################################
def fits_in(in_file, n_bins, dt, print_iterator, test):
	"""
			fits_in
	
	Opens the .fits GTI'd event list, reads the count rate for a segment, 
	populates the light curve, calls 'make_ps' to create a power spectrum, adds 
	power spectra over all segments.
			 
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
	all_time = np.asarray(data.field('TIME'), dtype=np.float64)
	all_energy = np.asarray(data.field('CHANNEL'), dtype=np.float64)

	while end_time <= final_time:
			
		time = all_time[np.where(all_time < end_time)]
		energy = all_energy[np.where(all_time < end_time)]
		
		for_next_iteration = np.where(all_time >= end_time)
		all_time = all_time[for_next_iteration]
		all_energy = all_energy[for_next_iteration]

# 		print "Len of time:", len(all_time)
# 		print "\nLen time:", len(time)
# 		print "Start: %.21f" % start_time
# 		print "End  : %.21f" % end_time
# 		if len(time) == 1:
# 			print "One: %.21f" % time[0]# 		if len(this_time) == 0:
# 			print "Len of time = 0"
		if len(time) > 0:
			num_segments += 1
			rate_2d, rate_1d = tools.make_lightcurve(time, 
				energy, n_bins, dt, start_time)
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
			time = None
			energy = None
			power_segment = None
			mean_rate_segment = None
			rate_2d = None
			rate_1d = None
			
			if test and (num_segments == 1):  # Testing
				np.savetxt('lightcurve.dat', lightcurve, fmt='%d')
				break
			start_time += (n_bins * dt)
			end_time += (n_bins * dt)
# 			if num_segments == 170:
# 				break
		elif len(time) == 0:
			start_time = all_time[0]
			end_time = start_time + (n_bins * dt)
		## End of 'if there are counts in this segment'

	## End of while-loop
# 	print "Final end time: %.21f" % end_time
	return power_sum, sum_rate_whole, num_segments
## End of function 'fits_in'


################################################################################
def dat_in(in_file, n_bins, dt, print_iterator, test):
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
# 	print "First start time: %.21f" % start_time
# 	print "First end   time: %.21f" % end_time

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
				
				if next_time >= end_time:  # Triggered at the end of a segment

# 					print "\nLen time:", len(time)
# 					print "Start: %.21f" % start_time
# 					print "End  : %.21f" % end_time
# 					if len(time) == 1:
# 						print "One: %.21f" % time[0]

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
# 							print "\t%.21f" % end_time

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

						start_time += (n_bins * dt)
						end_time += (n_bins * dt)
						
					## This next bit helps it handle gappy data; keep in mind 
					## that end_time has already been incremented here
					elif len(time) == 0:
						start_time = next_time
                		end_time = start_time + (n_bins * dt)	
					## End of 'if there are counts in this segment'
				## End of 'if it`s at the end of a segment'
			## End of 'if the line is not a comment'
		## End of for-loop 
	## End of with-block

	return power_sum, sum_rate_whole, num_segments
## End of function 'dat_in'
	

################################################################################
def read_and_use_segments(in_file, n_bins, dt, test):
	"""
			read_and_use_segments
			
	Opens the file, reads in the count rate, calls 'make_ps' to create a
	power spectrum. Separated from main body like this so I can easily call it 
	in multi_powerspec.py. Split into 'fits_in', 'dat_in', and
	'extracted_in' for easier readability.
		 
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
		power_sum, sum_rate_whole, num_segments = fits_in(in_file, 
			n_bins, dt, print_iterator, test)
	elif (in_file[-4:].lower() == ".dat"):
		power_sum, sum_rate_whole, num_segments = dat_in(in_file, 
			n_bins, dt, print_iterator, test)
	elif (in_file[-3:].lower() == ".lc"):
		power_sum, sum_rate_whole, num_segments = extracted_in(in_file, 
			n_bins, dt, print_iterator, test)
	else:
		raise Exception("ERROR: Input file type not recognized. Must be .dat, .fits, or .lc.")
	
	return power_sum, sum_rate_whole, num_segments
## End of function 'read_and_use_segments'
	
	
###############################################################################
def main(in_file, out_file, num_seconds, dt_mult, test):
	""" 
			main
	
	"""	
	
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
		n_bins, dt, num_seconds, num_segments, mean_rate_whole, True)
	
	## Output, based on file extension
	if out_file[-4:].lower() == "fits":
		fits_out(out_file, in_file, dt, n_bins, nyquist_freq, num_segments, \
			mean_rate_whole, freq, rms2_power, rms2_err_power, leahy_power)
	elif out_file[-3:].lower() == "dat":
		dat_out(out_file, in_file, dt, n_bins, nyquist_freq, num_segments, \
			mean_rate_whole, freq, rms2_power, rms2_err_power, leahy_power)
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
	parser.add_argument('-n', '--num_seconds', type=tools.type_power_of_two, \
		default=1, dest='num_seconds', help='Number of seconds in each Fourier \
		segment. Must be a power of 2. [1]')
	parser.add_argument('-m', '--dt_mult', type=tools.type_power_of_two, \
		default=1, dest='dt_mult', help='Multiple of 1/8192 seconds for \
		timestep between bins. Must be a power of 2. [1]')
	parser.add_argument('--test', action='store_true', dest='test', help='If \
		present, only does a short test run.')
	args = parser.parse_args()
			
	main(args.infile, args.outfile, args.num_seconds, args.dt_mult, args.test)

## End of program 'powerspec.py'
