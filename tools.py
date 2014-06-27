from astropy.io import fits
import numpy as np
import scipy as sp
import itertools

"""
		tools.py

This has not been rigorously tested.
There is no 'main' to this program; only helper methods to import and be called.

To assign the returned value to a variable in a bash script (in the directory containing
tools.py): 
var=$(python -c 'from tools import the_function; print the_function(passing variables)')

Written in Python 2.7 by A.L. Stevens, A.L.Stevens@uva.nl, 2013-2014

The scientific modules imported above, as well as python 2.7, can be downloaded in the 
Anaconda package, https://store.continuum.io/cshop/anaconda/

"""

#########################################################################################
def get_key_val(fits_file, ext, keyword):
	"""
			get_key_val
	
	Gets the value of a keyword from a FITS header. Keyword does not seem to be 
	case-sensitive.

	Passed: fits_file - str - The full path of the FITS file.
			ext - int - The FITS extension in which to search for the given keyword.
			keyword - str - The keyword for which you want the associated value.
			
	Returns: key_value - any - Value of the given keyword.
	
	"""
	
	ext = np.int8(ext)
	assert (ext >= 0 and ext <= 3)

	hdulist = fits.open(fits_file)
	key_value = hdulist[ext].header[keyword]
	hdulist.close()
	return key_value
	## End of function 'get_key_val'


#########################################################################################
def compute_obs_time(file_list):
	"""
			compute_obs_time
		
	Computes the total observation time of a list of observation FITS files, in seconds.
	
	Passed: file_list - str - Name of file with list of fits files of the observations.
	
	Returns: total_time - float - The total observation time.
	
	"""

	input_files = [line.strip() for line in open(file_list)]
	total_time = 0
	
	for file in input_files:
		start_time = float(get_keyword(file, 0, 'TSTART'))
		stop_time = float(get_keyword(file, 0, 'TSTOP'))
		time = stop_time - start_time
		total_time += time
		## End of for-loop

	return total_time # in seconds, but can double check units of times in FITS header
	## End of function 'compute_obs_time'

	
#########################################################################################
def read_obs_time(in_file):
	"""
		
		read_obs_time
		
	Read the total observation time from the header of a text file.
	
	Passed: in_file - str - Name of (ASCII/txt/dat) input file with observation duration 
				in line 5 (starting at 0) of the header.
		
	Returns: nothing
	
	"""
	
	with open(in_file, 'r') as f:
		for line in f:
			if line[0].strip() == "#":
				if "duration" in line.strip():
					line = line.strip().split()
# 					print line
					obs_time = float(line[len(line)-2])
					return obs_time
			else:
				return 0


#########################################################################################
def power_of_two(num):
	"""
			power_of_two
			
	Checks if 'num' is a power of 2 (1 <= num < 2147483648 )
	
	Passed: num - int - The number in question.
	
	Returns: a bool - 'True' if 'num' is a power of two, 'False' if 'num' is not.
	
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
	
	#########################################################################################
def pairwise(iterable):
	"""
			pairwise
	
	s -> (s0,s1), (s1,s2), (s2, s3), ...
	From https://docs.python.org/2/library/itertools.html#recipes
	Used when reading lines in the file so I can peek at the next line.
	
	Passed: an iterable, like a list or an open file
	
	Returns: the next two items in an iterable, like in the example a few lines above.
	
	"""

	a, b = itertools.tee(iterable)
	next(b, None)
	return itertools.izip(a, b)
	## End of function 'pairwise'
	
	
#########################################################################################
def replace_key_val(fits_file, ext, keyword, value):
	"""
			replace_key_val
			
	Replaces the value of a keyword in a FITS header with a given value.
	
	Passed: fits_file - str - Name of a FITS file.
			ext - int - The FITS extension in which you want to replace the keyword value.
			keyword - str - The keyword of the value you want to replace.
			value - any - The new value for the FITS header keyword.
	
	Returns: nothing
	
	"""
	ext = np.int8(ext)
	assert (ext >= 0 and ext <= 3)
	
	hdu = fits.open(fits_file, mode='update')
	# print "Before replacing keyword value, %s = %f" %(keyword, hdu[ext].header[keyword])
	hdu[ext].header[keyword] = value
	# print "After replacing keyword value, %s = %f" %(keyword, hdu[ext].header[keyword])
	hdu.flush()
	hdu.close()
	## End of function 'replace_key_val'
	
	
#########################################################################################
def obs_epoch_rxte(fits_file):
	"""
			obs_epoch_rxte
			
	Determines the epoch of an RXTE observation. Returns 0 if an error occurred.
	Future update: use MJD.

	WARNING:
	1. This has not been rigorously tested (barely tested, really) and is not guaranteed.
	2. It can only understand and parse 'DATE-OBS' keyword values with length 8 or 
		length 19.
	3. I'm interpreting the 'stop time' listed as the start time of the next epoch.
	
	Passed: fits_file - str - Name of an RXTE observation FITS file.
	
	Returns: epoch - int - The RXTE observation epoch of the FITS observation file.
	
	"""
	
	obs_time = get_key_val(fits_file, 0, 'DATE-OBS')
# 	print "DATE-OBS =", obs_time
# 	print "Length of DATE-OBS =", len(obs_time)
	year = -1
	month = -1
	day = -1
	hour = -1
	minute = -1
	
	if len(obs_time) == 19:
		year = int(obs_time[0:4])
		month = int(obs_time[5:7])
		day = int(obs_time[8:10])
		hour = int(obs_time[11:13])
		minute = int(obs_time[14:16])
# 		print "Year =", year
# 		print "Month =", month
# 		print "Day =", day
# 		print "Hour =", hour
# 		print "Minute =", minute
		
	elif len(obs_time) == 8:
		day = int(obs_time[0:2])
		month = int(obs_time[3:5])
		year = obs_time[6:8]
		if year[0] == '9': year = int("19"+year)
		else: year = int("20"+year)
# 		print "Day =", day
# 		print "Month =", month
# 		print "Year =", year
		hour = 0
		minute = 0
	else:
		print "\n\tERROR: Format of date is not understood."
		return 0
		
	assert (year >= 1995 and year <= 2012) ## Making sure that the date is actually when 
										   ##  RXTE was operational
	if (year is -1) or (month is -1) or (day is -1) or (hour is -1) or (minute is -1):
		print "\n\tERROR: Month, date, year, hour, or minute not properly assigned."
		return 0
	
	## Determining which epoch the date falls in
	
	if year == 1995: 
		return 1
		
	elif year == 1996:
		if month < 3:
			return 1
		elif month == 3:
			if day < 3: return 1
			elif day > 3: return 2 
			else: 
				if hour < 18: return 1
				elif hour > 18: return 2
				else:
					print "\n\tWARNING: Down to the minute in determining obs epoch. May not be correct."
					if minute < 33: return 1
					else: return 2
		elif month == 4:
			if day < 15: return 2
			elif day > 15: return 3
			else: 
				if hour < 23: return 2
				elif hour > 23: return 3
				else:
					print "\n\tWARNING: Down to the minute in determining obs epoch. May not be correct."
					if minute < 5: return 2
					else: return 3
		else:
			return 3
		
	elif year == 1997 or year == 1998:
		return 3
	
	elif year == 1999:
		if month < 3:
			return 3
		elif month == 3:
			if day < 22: return 3
			elif day > 22: return 4
			else:
				if hour < 17: return 3
				elif hour > 17: return 4
				else:
					print "\n\tWARNING: Down to the minute in determining obs epoch. May not be correct."
					if minute < 37: return 3
					else: return 4
		else: 
			return 4
		
	elif year == 2000:
		if month < 5:
			return 4
		elif month == 5:
			if day < 13: return 4
			elif day >= 13: return 5  # since it changes at 00:00
		else:
			return 5
		
	else:
		return 5
	## End of function 'obs_epoch_rxte'
	