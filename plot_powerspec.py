import argparse
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import matplotlib.font_manager as font_manager
from matplotlib.ticker import MultipleLocator

__author__ = "Abigail Stevens"
__author_email__ = "A.L.Stevens at uva.nl"
__year__ = "2013-2015"
__description__ = "Linearly plots a power spectrum in the frequency domain."

"""
		plot_powerspec.py

Written in Python 2.7.

"""
################################################################################
def main(tab_file, plot_file, prefix):
	"""
			main
	
	Linearly plots a power spectrum in the frequency domain.
	
	"""
	##########################################
	## Reading in power spectrum from a table
	##########################################
		
	if args.tab_file[-4:].lower() == ".dat":
		table = np.loadtxt(tab_file, comments='#')
		freq = np.asarray(table[:,0])  # frequency, in Hz
		fracrms = np.asarray(table[:,1])  # fractional rms^2 power
		error = np.asarray(table[:,2])  # error on power
	elif args.tab_file[-5:].lower() == ".fits":
		file_hdu = fits.open(tab_file)
		table = file_hdu[1].data
		file_hdu.close()
		freq = table.field('FREQUENCY')  # frequency, in Hz
		fracrms = table.field('POWER')  # fractional rms^2 power
		error = table.field('ERROR')  # error on power
	else:
		raise Exception('ERROR: File type not recognized. Must have extension \
.dat or .fits.')
	
	
	#####################################
	## Plotting the power spectrum (psd)
	#####################################
	
	font_prop = font_manager.FontProperties(size=18)
	
	print "Power spectrum: %s" % plot_file

	fig, ax = plt.subplots(1,1)
	ax.plot(freq, fracrms, linewidth=2)
# 	ax.errorbar(freq, fracrms, xerr=None, yerr=error)

	ax.set_xlim(freq[1],np.max(freq))
	ax.set_ylim(0, 0.012)
	## Setting the axes' minor ticks. It's complicated.
	x_maj_loc = ax.get_xticks()
	y_maj_loc = ax.get_yticks()
	x_min_mult = 0.2 * (x_maj_loc[1] - x_maj_loc[0])
	y_min_mult = 0.2 * (y_maj_loc[1] - y_maj_loc[0])
	xLocator = MultipleLocator(x_min_mult)  ## location of minor ticks on the x-axis
	yLocator = MultipleLocator(y_min_mult)  ## location of minor ticks on the y-axis
	ax.xaxis.set_minor_locator(xLocator)
	ax.yaxis.set_minor_locator(yLocator)

# 	ax.set_xscale('log')
	ax.set_xlabel('Frequency (Hz)', fontproperties=font_prop)
	ax.set_ylabel(r'Power (frac. rms$^{2}$)', \
		fontproperties=font_prop)
	ax.tick_params(axis='x', labelsize=16, bottom=True, top=True, \
		labelbottom=True, labeltop=False)
	ax.tick_params(axis='y', labelsize=16, left=True, right=True, \
		labelleft=True, labelright=False)
	ax.set_title(prefix, fontproperties=font_prop)
	
	fig.set_tight_layout(True)
	plt.savefig(plot_file, dpi=100)  ## set dpi higher for better image resolution
# 	plt.show()
	plt.close()
## End of function 'main'
	

################################################################################
if __name__ == "__main__":
	
	###########################
	## Parsing input arguments
	###########################
	
	parser = argparse.ArgumentParser(usage="python plot_powerspec.py tab_file \
[-o plot_file] [-p prefix]", description="Plots a power spectrum from a data \
table.", epilog="For optional arguments, default values are given in brackets \
at end of description.")

	parser.add_argument('tab_file', help="The table file, in .dat or .fits \
format, with frequency in column 1, fractional rms^2 power in column 2, and \
error on power in column 3.")

	parser.add_argument('-o', '--outfile', nargs='?', dest='plot_file', \
default="./psd.png", help="The output file name for the power spectrum plot. \
[./psd.png]")

	parser.add_argument('-p', '--prefix', nargs='?', dest='prefix', \
default="--", help="The identifying prefix of the data (object nickname or \
proposal ID). [--]")
		
	args = parser.parse_args()
	
	main(args.tab_file, args.plot_file, args.prefix)
	

################################################################################
