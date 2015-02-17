import argparse
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import matplotlib.font_manager as font_manager

__author__ = "Abigail Stevens"
__author_email__ = "A.L.Stevens@uva.nl"
__year__ = "2013-2015"
__description__ = "Plots a linear power spectrum in the frequency domain."

"""
		plot_powerspec.py

Written in Python 2.7.

All modules imported above, as well as python 2.7, can be downloaded in the 
Anaconda package, https://store.continuum.io/cshop/anaconda/

"""

################################################################################
if __name__ == "__main__":
	
	###########################
	## Parsing input arguments
	###########################
	
	parser = argparse.ArgumentParser(usage='plot_powerspec.py tab_file [-o \
plot_file] [-p prefix]', description="Plots a power spectrum from a data table.\
", epilog="For optional arguments, default values are given in brackets at end \
of description.")

	parser.add_argument('tab_file', help="The table file, in .dat or .fits \
format, with frequency in column 1, fractional rms^2 power in column 2, and \
error on power in column 3.")

	parser.add_argument('-o', '--outfile', nargs='?', dest='plot_file', \
default="./psd.png", help="The output file name for the power spectrum plot. \
[./psd.png]")

	parser.add_argument('-p', '--prefix', nargs='?', dest='prefix', \
default="x", help="The identifying prefix of the data (object nickname or \
proposal ID). [x]")
		
	args = parser.parse_args()
	
	##########################################
	## Reading in power spectrum from a table
	##########################################
	
	print "Power spectrum: %s" % args.plot_file
	
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
		raise Exception('ERROR: File type not recognized. Must have extension .dat or .fits.')
		
	vpv = freq * rms2
	
	#####################################
	## Plotting the power spectrum (psd)
	#####################################
	
	font_prop = font_manager.FontProperties(size=16)
	
	fig, ax = plt.subplots(1,1)
	ax.plot(freq, rms2, linewidth=2)
# 	ax.errorbar(freq, rms2, xerr=None, yerr=error)
	ax.set_xlabel(r'$\nu$ [Hz]', fontproperties=font_prop)
	ax.set_ylabel(r'Power, noise-subtracted fractional rms$^2$', fontproperties=font_prop)
# 	ax.set_xlim(0,800)
# 	ax.set_ylim(0, )
	ax.tick_params(axis='x', labelsize=14)
	ax.tick_params(axis='y', labelsize=14)
# 	ax.set_title("Power spectrum, " + args.prefix, fontproperties=font_prop)
	
	fig.set_tight_layout(True)
	plt.savefig(args.plot_file, dpi=120)
# 	plt.show()
	plt.close()
	
## End of program 'plot_powerspec.py'

################################################################################
