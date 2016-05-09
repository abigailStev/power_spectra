#!/bin/bash

################################################################################
##
## Bash script to run powerspec.py, plot_powerspec.py, rebin_powerspec.py,
## and fit_qpo.py.
##
## Runs powerspec.py for many obsIDs.
## 
## Don't give command line arguments. Change things in this script below.
##
## Change the directory names and specifiers before the double '#' row to best
## suit your setup.
##
## Notes: bash 3.* and Python 2.7.* (with supporting libraries) must be 
##		  installed in order to run this script. For the gif-making to work, 
##		  ImageMagick must be installed (open source, available on e.g. MacPorts 
##		  and HomeBrew)
## 
## Author: Abigail Stevens <A.L.Stevens at uva.nl> 2015
## 
################################################################################

##########################################
## Checking the number of input arguments
##########################################

if (( $# != 0 )); then
    echo -e "\tDo not give command line arguments. Usage: ./loop_powerspec.sh\n"
    exit
fi

################################################################################

home_dir=$(ls -d ~)

exe_dir="$home_dir/Dropbox/Research/power_spectra"
out_dir="$exe_dir/out_ps"

prefix="GX339-BQPO"

day=$(date +%y%m%d)  # make the date a string and assign it to 'day'

numsec=64
dt=64
testing=1    ## 0 for no, 1 for yes
rebin_const=1.01

obsID_list="$home_dir/Dropbox/Lists/${prefix}_obsIDs.lst"

plot_list="$out_dir/${prefix}_giflist.txt"
rb_plot_list="$out_dir/${prefix}_giflist_rb.txt"
gif_name="$out_dir/${prefix}_${day}_t${dt}_${numsec}sec.gif"
rb_gif_name="$out_dir/${prefix}_${day}_t${dt}_${numsec}sec_rb.gif"
qpofit_file="$out_dir/${prefix}_QPOfit.txt"

p_ext="png"

################################################################################
################################################################################

if [ -e "$plot_list" ]; then rm "$plot_list"; fi; touch "$plot_list"
if [ -e "$rb_plot_list" ]; then rm "$rb_plot_list"; fi; touch "$rb_plot_list"
if [ -e "$qpofit_file" ]; then rm "$qpofit_file"; fi; touch "$qpofit_file"

##########################
## Looping through obsIDs
##########################

for obsID in $( cat $obsID_list ); do

	in_dir="$home_dir/Reduced_data/$prefix/$obsID"
	in_file="$in_dir/GTId_eventlist.fits"
		
	if [ ! -d "$out_dir" ]; then mkdir -p "$out_dir"; fi

	if (( testing == 0 )); then
		out_file="$out_dir/${obsID}_${day}_t${dt}_${numsec}sec"
		plot_file="$out_dir/${obsID}_${day}_t${dt}_${numsec}sec"
	elif (( testing == 1 )); then
		out_file="$out_dir/test_${obsID}_${day}_t${dt}_${numsec}sec"
		plot_file="$out_dir/test_${obsID}_${day}_t${dt}_${numsec}sec"
	fi

	rb_out_file="${out_file}_rb"
	rb_plot="${plot_file}_rb"

	########################
	## Running powerspec.py
	########################
    echo "$in_file"

	if [ -e "$in_file" ]; then
		python "$exe_dir"/powerspec.py "$in_file" "${out_file}.fits"\
			-n "$numsec" -m "$dt" #-t "$testing"
	else
		echo -e "\tERROR: powerspec.py was not run. Eventlist does not exist."
	fi

	#################################################
	## Plotting the power spectrum and re-binning it
	#################################################

	if [ -e "${out_file}.fits" ]; then
		
		#####################
		## Plotting linearly
		#####################
		
		python "$exe_dir"/plot_powerspec.py "${out_file}.fits" \
			-o "${plot_file}.${p_ext}" -p "$prefix/$obsID"
# 		if [ -e "${plot_file}.${p_ext}" ]; then open "${plot_file}.${p_ext}"; fi
		
		echo "${plot_file}.${p_ext}" >> ${plot_list}

		###########################################
		## Re-binning and plotting logarithmically
		###########################################
		python "$exe_dir"/rebin_powerspec.py "${out_file}.fits" \
			"${rb_out_file}.fits" -o "${rb_plot}.${p_ext}" \
			-p "$prefix/$obsID" -c "$rebin_const"
# 		if [ -e "${rb_plot}.${p_ext}" ]; then open "${rb_plot}.${p_ext}"; fi
	
		echo "${rb_plot}.${p_ext}" >> ${rb_plot_list}
		
		python fit_qpo.py "${rb_out_file}.fits" --mod G --prefix "$prefix"
		
	else
		echo -e "\tERROR: Plots were not made. Power spectrum output file does"\
				"not exist."
	fi
	
	## If testing, be done after 1 obsID
	if (( testing == 1 )); then
		break
	fi

done
echo "$qpofit_file"

#open -a "TextWrangler" "$qpofit_file"
#############################################
## Make a gif of the plots using ImageMagick
#############################################

echo -e "\nGIF lists:"
echo -e "\t${plot_list}"
echo -e "\t${rb_plot_list}\n"

# convert @"${plot_list}" "${gif_name}"
# if [ -e "${gif_name}" ]; then
# 	echo "GIF made! ${gif_name}"
# 	open "${gif_name}"
# fi
# 
# convert @"${rb_plot_list}" "${rb_gif_name}"
# if [ -e "${rb_gif_name}" ]; then
# 	echo "GIF made! ${rb_gif_name}"
# 	open "${rb_gif_name}"
# fi

################################################################################
## All done!
################################################################################
