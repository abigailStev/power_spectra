import argparse
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import matplotlib.font_manager as font_manager

__author__ = "Abigail Stevens"
__author_email__ = "A.L.Stevens@uva.nl"
__year__ = "2013-2014"
__description__ = "Plots a power spectrum in the frequency domain."

"""
		plot_powerspec.py

Written in Python 2.7.

All modules imported above, as well as python 2.7, can be downloaded in the 
Anaconda package, https://store.continuum.io/cshop/anaconda/

"""

###############################################################################
if __name__ == "__main__":
	
	parser = argparse.ArgumentParser(description="Plots a power spectrum from \
		a data table.")
	parser.add_argument('tab_file', help="The table file, in .dat or .fits \
		format, with frequency in column 1, fractional rms^2 power in column 2,\
		and error on power in column 3.")
	parser.add_argument('-o', '--outfile', nargs='?', dest='plot_file', \
		default="./plot.png", help="The output file name for the power spectrum\
		plot.")
	parser.add_argument('-p', '--propID', nargs='?', dest='propID', \
		default="Pxxxxx", help="The proposal ID of the data.")
	args = parser.parse_args()
	
	print "Plotting the power spectrum: %s" % args.plot_file
	
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

	font_prop = font_manager.FontProperties(size=16)

	fig, ax = plt.subplots(1,1)
	ax.plot(freq, rms2, linewidth=2)
# 	ax.errorbar(freq, rms2, xerr=None, yerr=error)
	ax.set_xlabel(r'$\nu$ [Hz]', fontproperties=font_prop)
	ax.set_ylabel(r'Power, noise-subtracted fractional rms$^2$', fontproperties=font_prop)
	ax.set_xlim(0,1000)
	ax.set_ylim(0, )
	ax.tick_params(axis='x', labelsize=14)
	ax.tick_params(axis='y', labelsize=14)
# 	ax.set_title("Power spectrum, " + args.propID, fontproperties=font_prop)
	
	plt.savefig(args.plot_file, dpi=120)
# 	plt.show()
	plt.close()
	
## End of program 'plot_powerspec.py'
