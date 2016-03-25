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
## Author: Abigail Stevens <A.L.Stevens at uva.nl>, 2013-2016
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
#prefix="4U1608"
#prefix="H1743-BQPO"
#prefix="GRO1655-BQPO"

# ObsID of the data
obsID="95335-01-01-01"
# Multiple of time resolution of the data for binning the light curve
dt=64
# Number of seconds per Fourier segment
numsec=64
# Whether or not to run 'testing'; 0 for no, 1 for yes
testing=0
# Constant to re-bin the power spectrum by (must be > 1)
rebin_const=1.06
# Desired plot file extension, without the dot
plot_ext="eps"

# Your computer's home directory (gets automatically)
home_dir=$(ls -d ~)
# Today's date (gets automatically), for writing in file names
day=$(date +%y%m%d)  # make the date a string and assign it to 'day'
#day="150901"  # make the date a string and assign it to 'day'
psd_exe_dir="$home_dir/Dropbox/Research/power_spectra"
psd_out_dir="${psd_exe_dir}/out_ps/${prefix}"
#psd_out_dir="${psd_exe_dir}/out_ps"
list_dir="$home_dir/Dropbox/Lists"  ## A folder of lists; tells which files
									## we're using
in_dir="$home_dir/Reduced_data/$prefix/$obsID"
data_file="$in_dir/GTId_eventlist.fits"
out_local_name="${obsID}_${day}_t${dt}_${numsec}sec"
#data_file="$list_dir/${prefix}_eventlists_9.lst"
#out_local_name="${prefix}_${day}_t${dt}_${numsec}sec_adj"
#data_file="$list_dir/${prefix}_eventlists.lst"
#out_local_name="${prefix}_${day}_t${dt}_${numsec}sec"

################################################################################
################################################################################

if [ ! -d "${psd_out_dir}" ]; then mkdir -p "${psd_out_dir}"; fi

if (( $testing == 0 )); then
	out_file="${psd_out_dir}/${out_local_name}"
	plot_file="${psd_out_dir}/${out_local_name}"
	saved_file_list="${psd_out_dir}/${out_local_name}_filelist"

elif (( $testing == 1 )); then
	out_file="${psd_out_dir}/test_${out_local_name}"
	plot_file="${psd_out_dir}/test_${out_local_name}"
	saved_file_list="${psd_out_dir}/test_${out_local_name}_filelist"
fi

rb_out_file="${out_file}_rb"
rb_plot="${plot_file}_rb"

########################
## Running powerspec.py
########################

if [ -e "${data_file}" ]; then
	time python "${psd_exe_dir}"/powerspec.py "${data_file}" "${out_file}.fits" \
			-n "$numsec" -m "$dt" -t "$testing" --pcu 1
else
	echo -e "\tERROR: powerspec.py was not run. Data file or list of event "\
	        "lists does not exist."
fi

#########################################
## Plot the power spectrum and re-bin it
#########################################

#if [ -e "${out_file}.fits" ]; then
#
#	## Plotting linearly
#	python "${psd_exe_dir}"/plot_powerspec.py "${out_file}.fits" \
#			-o "${plot_file}.${plot_ext}" --prefix "$prefix"
#
## 	if [ -e "${plot_file}.${plot_ext}" ]; then open "${plot_file}.${plot_ext}"
## 	fi
#
	## Re-binning and plotting logarithmically
	python "${psd_exe_dir}"/rebin_powerspec.py "${out_file}.fits" \
			"${rb_out_file}.fits" -o "${rb_plot}.${plot_ext}" \
			--prefix "$prefix" -c "$rebin_const"

	if [ -e "${rb_plot}.${plot_ext}" ]; then open "${rb_plot}.${plot_ext}"; fi

#	python "${psd_exe_dir}"/fit_qpo.py "${rb_out_file}.fits" --mod "L"
#
#else
#	echo -e "\tERROR: Plots were not made. Power spectrum output file does not"\
#			" exist."
#fi

################################################################################
## All done!
################################################################################
