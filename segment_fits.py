from astropy.io import fits
import os.path
import subprocess
import tools

__author__ = "Abigail Stevens"


dt_mult=64
num_seconds=64
prefix="95335-01-01-01"
out_dir = "/Users/abigailstevens/Dropbox/Research/power_spectra/out_ps"
out_file = out_dir+"/"+prefix+"_t"+str(dt_mult)+str(num_seconds)+"sec"
in_file = "/Users/abigailstevens/Reduced_data/GX339-BQPO/95335-01-01-01/GTId_eventlist.fits"

t_res = float(tools.get_key_val(in_file, 0, 'TIMEDEL'))
dt = dt_mult * t_res
n_bins = num_seconds * int(1.0 / dt)
nyquist_freq = 1.0 / (2.0 * dt)
df = 1.0 / float(num_seconds)

meta_dict = {'dt': dt, 't_res': t_res, 'num_seconds': num_seconds, \
             'df': df, 'nyquist': nyquist_freq, 'n_bins': n_bins, \
             'detchans': 64}

print "\nDT = %f seconds" % meta_dict['dt']
print "N_bins = %d" % meta_dict['n_bins']
print "Nyquist freq =", meta_dict['nyquist']

try:
    fits_hdu = fits.open(in_file)
except IOError:
    print "\tERROR: File does not exist: %s" % in_file
    exit()

header = fits_hdu[0].header	 ## Header is in ext 0
data = fits_hdu[1].data  ## Data is in ext 1
fits_hdu.close()


if meta_dict['n_bins'] == 32768:
    print_iterator = int(10)
elif meta_dict['n_bins'] < 32768:
    print_iterator = int(20)
elif meta_dict['n_bins'] >= 2097152:
    print_iterator = int(1)
elif meta_dict['n_bins'] >= 1048576:
    print_iterator = int(2)
else:
    print_iterator = int(5)
sum_rate_whole = 0
power_sum = np.zeros(meta_dict['n_bins'], dtype=np.float64)
num_seg = 0
lightcurve = np.asarray([])

start_time = data.field('TIME')[0]
final_time = data.field('TIME')[-1]
end_time = start_time + meta_dict['num_seconds']

## Filter data based on pcu
# 	PCU2_mask = data.field('PCUID') == 2
# 	data = data[PCU2_mask]

## Filter data based on energy channel
# 	print np.shape(data)
# 	lower_bound = data.field('CHANNEL') >= 0
# 	data = data[lower_bound]
# 	upper_bound = data.field('CHANNEL') <= 27
# 	data = data[upper_bound]
# 	print np.shape(data)

all_time = np.asarray(data.field('TIME'), dtype=np.float64)
# 	all_energy = np.asarray(data.field('CHANNEL'), dtype=np.float64)

################################
## Looping through the segments
################################

while end_time <= final_time:

    time = all_time[np.where(all_time < end_time)]
# 		energy = all_energy[np.where(all_time < end_time)]

    for_next_iteration = np.where(all_time >= end_time)
    all_time = all_time[for_next_iteration]
# 		all_energy = all_energy[for_next_iteration]

    if len(time) > 0:
        num_seg += 1
        rate_1d = tools.make_1Dlightcurve(time, meta_dict['n_bins'], \
            meta_dict['dt'], start_time)
        lightcurve = np.concatenate((lightcurve, rate_1d))

        power_segment, mean_rate_segment = make_ps(rate_1d)
        assert int(len(power_segment)) == meta_dict['n_bins'], "ERROR: "\
            "Something went wrong in make_ps. Length of power spectrum "\
            "segment  != n_bins."
        power_sum += power_segment
        sum_rate_whole += mean_rate_segment

        ## Printing out which segment we're on every x segments
        if num_seg % print_iterator == 0:
            print "\t", num_seg

        ## Clearing variables from memory
        time = None
# 			energy = None
        power_segment = None
        mean_rate_segment = None
        rate_1d = None

        if test and (num_seg == 1):  # Testing
            np.savetxt('tmp_lightcurve.dat', lightcurve, fmt='%d')
            break
        start_time += meta_dict['num_seconds']
        end_time += meta_dict['num_seconds']

    elif len(time) == 0:
        print "No counts in this segment."
        start_time = all_time[0]
        end_time = start_time + meta_dict['num_seconds']
    ## End of 'if there are counts in this segment'

return power_sum, sum_rate_whole, num_seg