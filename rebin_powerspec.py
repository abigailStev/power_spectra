import argparse
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
from datetime import datetime
import os
from tools import type_positive_float

__author__ = "Abigail Stevens"
__author_email__ = "A.L.Stevens@uva.nl"
__year__ = "2013-2015"
__description__ = "Plots a log power spectrum that has been re-binned in \
frequency, in the frequency domain."

"""
		plot_rb_powerspec.py

Written in Python 2.7.

All scientific modules imported above, as well as python 2.7, can be downloaded 
in the Anaconda package, https://store.continuum.io/cshop/anaconda/

tools is in the whizzy_scripts git repo.

"""

################################################################################
def fits_out(out_file, rb_out_file, dt, n_bins, nyquist_freq, num_segments, \
	mean_rate_whole, rebin_const, rb_freq, rb_rms2, rb_err):
	"""
				fits_out
			
	Writes a frequency re-binned power spectrum to a FITS file.
	
	"""
	
	print "Re-binned output file: %s" % rb_out_file

	## Updating above header for re-binned power spectrum (extension 0)
	prihdr = fits.Header()
	prihdr.set('TYPE', "Re-binned power spectrum")
	prihdr.set('DATE', str(datetime.now()), "YYYY-MM-DD localtime")
	prihdr.set('UNBINOUT', out_file, "Corresponding un-binned output.")
	prihdr.set('REBIN', rebin_const, "Freqs re-binned by REBIN * prev_bin_size")
	prihdr.set('DT', dt, "seconds")
	prihdr.set('N_BINS', n_bins, "time bins per segment")
	prihdr.set('SEGMENTS', num_segments, "segments in the whole light curve")
	prihdr.set('EXPOSURE', num_segments * n_bins * dt, "seconds, of light curve")
	prihdr.set('MEANRATE', mean_rate_whole, "counts/second")
	prihdr.set('NYQUIST', nyquist_freq, "Hz")
	prihdu = fits.PrimaryHDU(header=prihdr)
	
	## Making FITS table for re-binned power spectrum (extension 1)
	col1 = fits.Column(name='FREQUENCY', unit='Hz', format='E', \
		array=rb_freq)
	col2 = fits.Column(name='POWER', unit='frac rms^2', format='E', \
		array=rb_rms2)
	col3 = fits.Column(name='ERROR', unit='frac rms^2', format='E', \
		array=rb_err)
	cols = fits.ColDefs([col1, col2, col3])
	tbhdu = fits.BinTableHDU.from_columns(cols)
	
	## If the file already exists, remove it (still working on just updating it)
	assert rb_out_file[-4:].lower() == "fits", \
		'ERROR: Re-binned output file must have extension ".fits".'
	if os.path.isfile(rb_out_file):
# 		print "File previously existed. Removing and rewriting."
		os.remove(rb_out_file)
	
	## Writing the re-binned power spectrum to a FITS file
	thdulist = fits.HDUList([prihdu, tbhdu])
	thdulist.writeto(rb_out_file)	
	
## End of function 'fits_out'


################################################################################
def geometric_rebinning(freq, rms2_power, rms2_err_power, rebin_const):
	"""
			geometric_rebinning
			
	Re-bins the noise-subtracted fractional rms^2 power spectrum in frequency 
	space by some re-binning constant (rebin_const > 1). 

	"""
	
	## Initializing variables
	rb_rms2_power = np.asarray([])  # List of re-binned fractional rms power
	rb_freq = np.asarray([])		  # List of re-binned frequencies
	rb_err = np.asarray([])	  # List of error in re-binned power
	real_index = 1.0		  # The unrounded next index in power
	int_index = 1			  # The int of real_index, added to current_m every 
							  #  iteration
	current_m = 1			  # Current index in power
	prev_m = 0				  # Previous index m
	bin_power = 0.0			  # The power of the current re-binned bin
	bin_freq = 0.0			  # The frequency of the current re-binned bin
	err_bin_power2 = 0.0	  # The error squared on 'bin_power'
	bin_range = 0.0			  # The range of un-binned bins covered by this 
							  #  re-binned bin
	
	## Looping through the length of the array power, geometric bin by 
	## geometric bin, to compute the average power and frequency of that
	## geometric bin.
	## Equations for frequency, power, and error on power are from Adam Ingram's
	## PhD thesis
	while current_m < len(rms2_power):
# 	while current_m < 100: # used for debugging
		## Initializing clean variables for each iteration of the while-loop
		bin_power = 0.0  # the averaged power at each index of rb_rms2_power
		err_bin_power2 = 0.0  # the square of the errors on powers in this bin
		bin_range = 0.0
		bin_freq = 0.0
		
		## Determining the range of indices this specific geometric bin covers
		bin_range = np.absolute(current_m - prev_m)
		## Want mean of data points contained within one geometric bin
		bin_power = np.mean(rms2_power[prev_m:current_m])
		## Computing error in bin -- equation from Adam Ingram's thesis
		err_bin_power2 = np.sqrt(np.sum(rms2_err_power[prev_m:current_m] ** 2))\
			/ float(bin_range) 
		
		## Computing the mean frequency of a geometric bin
		bin_freq = np.mean(freq[prev_m:current_m])

		## Appending values to arrays
		rb_rms2_power = np.append(rb_rms2_power, bin_power)
		rb_freq = np.append(rb_freq, bin_freq)
		rb_err = np.append(rb_err, err_bin_power2)
		
		## Incrementing for the next iteration of the loop
		## Since the for-loop goes from prev_m to current_m-1 (since that's how
		## the range function and array slicing works) it's ok that we set 
		## prev_m = current_m here for the next round. This will not cause any
		## double-counting bins or skipping bins.
		prev_m = current_m
		real_index *= rebin_const
		int_index = int(round(real_index))
		current_m += int_index
		bin_range = None
		bin_freq = None
		bin_power = None
		err_bin_power2 = None	
	## End of while-loop
	
	return rb_freq, rb_rms2_power, rb_err
## End of function 'geometric_rebinning'


###############################################################################
if __name__ == "__main__":
	
	###########################
	## Parsing input arguments
	###########################
	
	parser = argparse.ArgumentParser(usage='rebin_powerspec.py tab_file \
rb_out_file [-o plot_file] [-p prefix] [-c rebin_constant]', description="\
Geometrically re-bins in frequency and plots a power spectrum.")
	
	parser.add_argument('tab_file', help="The table file, in .dat or .fits \
format, with frequency in column 1, fractional rms^2 power in column 2, and \
error on power in column 3.")
	
	parser.add_argument('rb_out_file', help="The FITS file to write the re-\
binned power spectrum to.")
	
	parser.add_argument('-o', '--outfile', nargs='?', dest='plot_file', \
default="rb_plot.png", help="The output plot file name.")
	
	parser.add_argument('-p', '--prefix', nargs='?', dest='prefix', \
default="X", help="The identifying prefix for the file (proposal ID or object \
nickname).")
	
	parser.add_argument('-c', '--rebin_const', nargs='?', dest='rebin_const', \
type=type_positive_float, default="x", help="The constant by which the data \
will be geometrically re-binned.")
	
	args = parser.parse_args()
	
	assert args.rebin_const >= 1.0 , 'ERROR: Re-binning constant must be >= 1.' 
	
	##########################################
	## Reading in power spectrum from a table
	##########################################
	
	if args.tab_file[-4:].lower() == ".dat":
	
		table = np.loadtxt(args.tab_file, comments='#')
		freq = np.asarray(table[:,0])  # frequency, in Hz
		rms2 = np.asarray(table[:,1])  # fractional rms^2 power
		error = np.asarray(table[:,2])  # error on power
		
	elif args.tab_file[-5:].lower() == ".fits":
	
		file_hdu = fits.open(args.tab_file)
		table = file_hdu[1].data
		file_hdu.close()
		freq = table.field('FREQUENCY')  # frequency, in Hz
		rms2 = table.field('POWER')  # fractional rms^2 power
		error = table.field('ERROR')  # error on power
		
	else:
	
		raise Exception('ERROR: File type not recognized. Must have extension \
.dat or .fits.')
	
	################################################
	## Re-binning the power spectrum by rebin_const
	################################################
	
	rb_freq, rb_rms2, rb_err = geometric_rebinning(freq, rms2, error, \
		args.rebin_const)
	
	########################################
	## Want to plot nu * P(nu) in log space 
	########################################
	
	vpv = rb_freq * rb_rms2
	err_vpv = rb_freq * rb_err
	
	#############
	## Plotting!
	#############
	
	font_prop = font_manager.FontProperties(size=16)
	print "Re-binned power spectrum: %s" % args.plot_file

	fig, ax = plt.subplots(1,1)
# 	ax.plot(rb_freq, vpv, lw=2)
	ax.errorbar(rb_freq, vpv, yerr=err_vpv, lw=2, c='blue', elinewidth=1, capsize=1)
	ax.set_xscale('log')
	ax.set_yscale('log')
# 	ax.loglog(rb_freq, vpv, lw=2, basex=10)
# 	ax.set_xlim(freq[0], 3000 )
# 	ax.set_xlim(100,3000)
	ax.set_ylim(0,)
	ax.set_xlabel(r'$\nu$ [Hz]', fontproperties=font_prop)
	ax.set_ylabel(r'$\nu$ $\cdot$ P($\nu$) [Hz rms$^2$]', \
		fontproperties=font_prop)
# 	ax.set_ylabel(r'Power, noise-subtracted fractional rms$^2$', \
# 		fontproperties=font_prop)
	ax.tick_params(axis='x', labelsize=14)
	ax.tick_params(axis='y', labelsize=14)
	ax.set_title("Power spectrum, " + args.prefix + ", Re-bin constant = " +\
		str(args.rebin_const), fontproperties=font_prop)

	## The following legend code was found on stack overflow I think
# 	legend = ax.legend(loc='lower right')
# 	## Set the fontsize
# 	for label in legend.get_texts():
# 		label.set_fontsize('small')
# 	for label in legend.get_lines():
# 		label.set_linewidth(2)  # the legend line width

	fig.set_tight_layout(True)
	plt.savefig(args.plot_file, dpi=120)
# 	plt.show()
	plt.close()
	
	##########################################################
	## Writing the re-binned power spectrum to an output file 
	##########################################################
	
	file_hdu = fits.open(args.tab_file)
	dt = file_hdu[0].header['DT']
	n_bins = file_hdu[0].header['N_BINS']
	nyquist_freq = file_hdu[0].header['NYQUIST']
	num_segments = file_hdu[0].header['SEGMENTS']
	mean_rate_whole = file_hdu[0].header['MEANRATE']
	file_hdu.close()
		
	fits_out(args.tab_file, args.rb_out_file, dt, n_bins, nyquist_freq, 
	num_segments, mean_rate_whole, args.rebin_const, rb_freq, rb_rms2, rb_err)
	
## End of program 'rebin_powerspec.py'

################################################################################
