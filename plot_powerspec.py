import argparse
import numpy as np
import matplotlib.pyplot as plt

"""
		plot_powerspec.py

Plots a power spectrum.

Required arguments:
tab_file - str - Name of file holding table of data (output from powerspec.py).

Optional arguments:
plot_file - str - Name of file which the plot will be saved to in this script.
propID - str - The proposal ID of the data.

Written in Python 2.7 by A.L. Stevens, A.L.Stevens@uva.nl, 2013-2014

All modules imported above, as well as python 2.7, can be downloaded in the 
Anaconda package, https://store.continuum.io/cshop/anaconda/

"""

###############################################################################
if __name__ == "__main__":
	
	parser = argparse.ArgumentParser(description="Plots a power spectrum from \
		a data table.")
	parser.add_argument('tab_file', help="The input table file, in \
		ASCII/txt/dat format, with frequency in column 1, fractional rms^2 \
		power in column 2, and error on power in column 3.")
	parser.add_argument('-o', '--outfile', nargs='?', dest='plot_file', \
		default="./plot.png", help="The output file name for the power spectrum\
		plot.")
	parser.add_argument('-p', '--propid', nargs='?', dest='propID', \
		default="Pxxxxx", help="The proposal ID of the data.")
	args = parser.parse_args()
	
	print "Plotting the power spectrum: %s" % args.plot_file
	
	table = np.loadtxt(args.tab_file, comments='#')
	freq = np.asarray(table[:,0])  # frequency, in Hz
	rms2 = np.asarray(table[:,1])  # fractional rms^2 power
	error = np.asarray(table[:,2])  # error on avg power
	vpv = freq * rms2

	fig, ax = plt.subplots()
	ax.plot(freq, rms2, linewidth=2)
# 	plt.errorbar(freq, rms2, xerr=None, yerr=error)
	plt.xlabel(r'$\nu$ [Hz]')
	plt.ylabel(r'Noise-subtracted fractional rms$^2$ power')
	plt.xlim(0,1000)
	plt.ylim(0, )
# 	plt.xscale('symlog') # this works much better than 'log'
# 	plt.yscale('symlog')
	plt.title("Power spectrum, " + args.propID)
	plt.savefig(args.plot_file, dpi=120)
# 	plt.show()
	plt.close()
	
## End of program 'plot_powerspec.py'
