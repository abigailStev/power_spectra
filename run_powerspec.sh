#!/bin/bash

################################################################################
##
## Bash script to run powerspec.py, plot_powerspec.py, and rebin_powerspec.py
##
## Don't give command line arguments. Change things in this script below.
##
## Change the directory names and specifiers before the double '#' row to best
## suit your setup.
##
## Notes: bash 3.* and Python 2.7.* (with supporting libraries) must be 
##		  installed in order to run this script. 
## 
## Written by Abigail Stevens <A.L.Stevens at uva.nl> 2013-2015
## 
################################################################################

##########################################
## Checking the number of input arguments
##########################################

if (( $# != 0 )); then
    echo -e "\tDo not give command line arguments. Usage: ./run_powerspec.sh\n"
    exit
fi

################################################################################

# Identifying prefix (object nickname or data ID)
prefix="GX339-BQPO"
# ObsID of the data
obsID="95335-01-01-01"
# Multiple of time resolution of the data for binning the light curve
dt=512
# Number of seconds per Fourier segment
numsec=64
# Whether or not to run 'testing'; 0 for no, 1 for yes
testing=0
# Constant to re-bin the power spectrum by (must be > 1)
rebin_const=1.01
adjust_seg=0
# Desired plot file extension, without the dot
plot_ext="eps"

# Your computer's home directory (gets automatically)
home_dir=$(ls -d ~)
# Today's date (gets automatically), for writing in file names
day=$(date +%y%m%d)  # make the date a string and assign it to 'day'
exe_dir="$home_dir/Dropbox/Research/power_spectra"
out_dir="$exe_dir/out_ps/${prefix}"
in_dir="$home_dir/Reduced_data/$prefix/$obsID"
in_file="$in_dir/GTId_eventlist.fits"


################################################################################
################################################################################

if [ ! -d "$out_dir" ]; then mkdir -p "$out_dir"; fi

if (( $testing == 0 )); then
	out_file="$out_dir/${obsID}_${day}_t${dt}_${numsec}sec"
	plot_file="$out_dir/${obsID}_${day}_t${dt}_${numsec}sec"
elif (( $testing == 1 )); then
	out_file="$out_dir/test_${obsID}_${day}_t${dt}_${numsec}sec"
	plot_file="$out_dir/test_${obsID}_${day}_t${dt}_${numsec}sec"
fi

rb_out_file="${out_file}_rb"
rb_plot="${plot_file}_rb"

########################
## Running powerspec.py
########################

if [ -e "$in_file" ]; then
	time python "$exe_dir"/powerspec.py "$in_file" "${out_file}.fits" \
			-n "$numsec" -m "$dt" -t "$testing" -a "$adjust_seg"
else
	echo -e "\tERROR: powerspec.py was not run. Eventlist does not exist."
fi

#########################################
## Plot the power spectrum and re-bin it
#########################################

if [ -e "${out_file}.fits" ]; then
	
	## Plotting linearly
	python "$exe_dir"/plot_powerspec.py "${out_file}.fits" \
			-o "${plot_file}.${plot_ext}" -p "$prefix"
		
 	if [ -e "${plot_file}.${plot_ext}" ]; then open "${plot_file}.${plot_ext}"
 	fi
	
	## Re-binning and plotting logarithmically
	python "$exe_dir"/rebin_powerspec.py "${out_file}.fits" \
			"${rb_out_file}.fits" -o "${rb_plot}.${plot_ext}" -p "$prefix" \
			-c "$rebin_const"
		
	if [ -e "${rb_plot}.${plot_ext}" ]; then open "${rb_plot}.${plot_ext}"; fi

	python "$exe_dir"/fit_qpo.py "${rb_out_file}.fits" --mod "G"

else
	echo -e "\tERROR: Plots were not made. Power spectrum output file does not"\
			" exist."
fi

################################################################################
## All done!
################################################################################
