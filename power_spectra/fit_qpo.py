#!//anaconda/bin/python
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq ## Levenberg-Marquadt Algorithm
from astropy.io import fits
import argparse
import subprocess
import os.path

__author__ = "Abigail Stevens <A.L.Stevens at uva.nl>"

"""
Fits a QPO power spectrum (f*P(f)) with a power law and either a Gaussian or a
Lorentzian. Not intended to be robust, just to give an idea of how to quantify
the QPO in the power spectrum.

2015

"""

################################################################################
##
## DEFINING SPECTRAL SHAPE FUNCTIONS
## The power spectrum shape functions are *f because we're fitting f*P(f)
##
################################################################################

def powerlaw(f, p):
    ## p[3] = power law index, p[4] = scale factor
    pl = np.where(f != 0, f ** (p[3]), 0.0) * p[4] * f
    return pl

def lorentzian(f, p):
    ## p[0] = centroid frequency, p[1] = fwhm, p[2] = scale factor
    numerator = p[1] / (np.pi * 2.0)
    denominator = (f - p[0]) ** 2 + (1.0/2.0 * p[1]) ** 2
    L = (numerator / denominator) * p[2] * f
    return L

def gaussian(f, p):
    ## p[0] = mean value, p[1] = standard deviation, p[2] = scale factor
    exp_numerator = -(f - p[0])**2
    exp_denominator = 2 * p[1]**2
    G = p[2] * np.exp(exp_numerator / exp_denominator) * f
    return G

def pl_residuals(p, npn, f):
    ## Residuals from fitting a power law (this will be the background that
    ## gets subtracted)
    err = npn - powerlaw(f,p)
    return err

def L_residuals(p, npn, f):
    ## Residuals from fitting a Lorentzian (assumes background is already
    ## subtracted from npn, i.e., it's passed npn_bg_corr)
    err = npn - lorentzian(f,p)
    return err

def G_residuals(p, npn, f):
    ## Residuals from fitting a Gaussian (assumes background is already
    ## subtracted from npn, i.e., it's passed npn_bg_corr)
    err = npn - gaussian(f,p)
    return err

def L_pl_residuals(p, npn, f):
    ## Residuals from co-fitting a Lorentzian and a power law (assumes
    ## background is not previously subtracted)
    err = npn - lorentzian(f,p) - powerlaw(f,p)
    return err

def G_pl_residuals(p, npn, f):
    ## Residuals from co-fitting a Gaussian and a power law (assumes background
    ## is not previously subtracted)
    err = npn - gaussian(f,p) - powerlaw(f,p)
    return err


################################################################################
def make_plots(freq, npn_bg_corr, best_qpo, best_pl, best_resid, npn_err):

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(10,10))

    ax1.plot(freq, npn_bg_corr, 'wo')
    ax1.plot(freq, best_qpo,'r--', lw=2)
    # ax1.plot(freq, best_pl, 'b-', lw=2)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlim(np.min(freq), 100)
    ax1.set_ylim(1e-5, 1e-1)
# 	ax1.set_ylabel(r'Power $\times$ frequency (frac. rms$^{2}$ $\times$ Hz)', \
# 		fontsize=18)
    ax1.set_ylabel(r'Power $\times$ frequency', fontsize=18)

    ax2.errorbar(freq, best_resid, yerr=npn_err, linestyle=' ', ecolor='g', \
        elinewidth=2, markersize=2, marker='.')
    ax2.hlines(0, np.min(freq), np.max(freq), linestyle='dashed', lw=2)
    ax2.set_xscale('log')
    ax2.set_xlim(np.min(freq), 100)
    ax2.set_ylim(-0.003, 0.0025)
    ax2.set_xlabel('Frequency (Hz)', fontsize=18)
    ax2.set_ylabel(r'f $\cdot$ P(f) Residuals', fontsize=18)

    fig.subplots_adjust(hspace=0)
    plt.savefig('PSD_fit.png')
    # subprocess.call(["open", "PSD_fit.png"])

################################################################################
def get_fit(qpo_mod, freq, npn, npn_err):

    ##########################
    ## BACKGROUND SUBTRACTION
    ##########################

    ## Defining the 'background' part of the spectrum
    ind_bg_low = (freq > np.min(freq)) & (freq < 4.0)
    ind_bg_high = (freq > 7.0) & (freq < np.max(freq))

    freq_bg = np.concatenate((freq[ind_bg_low], freq[ind_bg_high]))
    npn_bg = np.concatenate((npn[ind_bg_low], npn[ind_bg_high]))

    ## Subtracting off mean background value (like next one but with m=0)
    print "Mean value of background", np.mean(npn_bg)
# 	npn_bg_corr = npn - np.mean(npn_bg)

    ## Fitting the background to a straight line
# 	log_freq_bg = np.log(freq_bg)
# 	log_npn_bg = np.log(npn_bg)
# 	m, c = np.polyfit(freq_bg, npn_bg, 1)
# 	print "Slope =", m,", Intercept = ", c
# 	## Removing fitted background
# 	background = m * freq + c
# 	npn_bg_corr = npn - background

    ## Fitting the background to a power law
    ## If using this, only fit a Lorentzian or Gaussian, not QPO + power law
# 	p_pl = [-1.0, 0.1]
# 	pbest = leastsq(pl_residuals, p_pl, args=(npn_bg, freq_bg), full_output=1)
# 	best_pl_p = pbest[0]
# 	print "Power law: index =", best_pl_p[0], ", scaling factor =", best_pl_p[1]
# 	best_pl = powerlaw(freq, best_pl_p)
# 	npn_bg_corr = pl_residuals(best_pl, npn, freq)

    npn_bg_corr = npn   ## If not setting these equal, will need to return
                        ## npn_bg_corr and

    ################
    ## FITTING DATA
    ################

    ## Giving initial fit parameter values
    p = [5.4651295, 0.3752125, 0.01, -0.5, 0.01]
    # p_qpo = [5.4651295, 0.3752125, 0.1]  ## Use this if subtracting off the power law first

    ## Optimizing using least squares method
    if qpo_mod == "L":
        pbest = leastsq(L_pl_residuals, p, args=(npn_bg_corr, freq), \
            full_output=1)
# 		pbest = leastsq(L_residuals, p_qpo, args=(npn_bg_corr, freq), \
# 			full_output=1)
    else:
        pbest = leastsq(G_pl_residuals, p, args=(npn_bg_corr, freq), \
            full_output=1)
# 		pbest = leastsq(G_residuals, p_qpo, args=(npn_bg_corr, freq), \
# 			full_output=1)

    ## Get the best parameters from the fit
    best_fit = pbest[0]
    print best_fit
    return best_fit


################################################################################
def main(in_file, qpo_mod, prefix):

    ################
    ## LOADING DATA
    ################

    try:
        file_hdu = fits.open(in_file)
    except IOError:
        print "\tERROR: File does not exist: %s" % in_file
        exit()

    data = file_hdu[1].data
    file_hdu.close()

    freq = data.field('FREQUENCY')
    power = data.field('POWER')
    npn = freq * power
    npn_err = data.field('ERROR') * freq

    freq = freq[1:]
    power = power[1:]
    npn = npn[1:]
    npn_err = npn_err[1:]

    #################################
    ## Get best fit of model to data
    #################################

    best_fit = get_fit(qpo_mod, freq, npn, npn_err)

    #########################
    ## Printing out best fit
    #########################

    if qpo_mod == "L":
        print "\nBest fit: Lorentzian + Power law"
        print "\tCentroid:", best_fit[0]
        fwhm = best_fit[1]
        ## Fit Lorentzian to data
        best_qpo = lorentzian(freq, best_fit)
        best_resid = L_pl_residuals(best_fit, npn, freq)

    else:
        print "\nBest fit: Gaussian + Power law"
        print "\tMean:", best_fit[0]
        print "\tStd dev:", best_fit[1]
        fwhm = 2.0 * np.sqrt(2.0 * np.log(2.0)) * best_fit[1]
        ## Fit Gaussian to data
        best_qpo = gaussian(freq, best_fit)
        best_resid = G_pl_residuals(best_fit, npn, freq)

    Q = best_fit[0] / fwhm
    print "\tFWHM:", fwhm
    print "\tQ value:", Q
    # print "\tQPO scale:", best_fit[2]
    # print "\tPL index:", best_fit[3]
    # print "\tPL scale:", best_fit[4]
    # scale_ratio = float(best_fit[4]) / float(best_fit[2])
    # print "QPO to PL scale ratio: 1.0:%.3e\n" % scale_ratio

    ## Fit power law to data
    best_pl = powerlaw(freq, best_fit)

    print "Mean of residuals:", np.mean(best_resid)

    ## Printing values to a table
    qpofit_file = os.path.dirname(in_file)+"/"+prefix+"_QPOfit.txt"
    print qpofit_file
    with open(qpofit_file, 'a') as out:
        out.write("%.4f \t %.4f \t %.4f \t %.5e\n" % (best_fit[0], fwhm, Q, \
                np.mean(best_resid)))

    ###################################################
    ## Plot power spectrum with best fit and residuals
    ###################################################

    make_plots(freq, npn, best_qpo, best_pl, best_resid, npn_err)


################################################################################
if __name__ == '__main__':

    ## Parsing input arguments and calling main

    parser = argparse.ArgumentParser(usage='infile [--mod QPO_MOD] [--prefix '\
            'PREFIX]', description="Fits a QPO power spectrum (f*P(f)) with a "\
            "power law and either a Gaussian or a Lorentzian. Not intended to "\
            "be robust, just to give an idea of how to quantify the QPO in the"\
            " power spectrum.", epilog="For optional arguments, default values"\
            " are given in brackets at end of description.")

    parser.add_argument('infile', help="Name of power spectrum to fit, in FITS"\
            " format. Usually feed it the re-binned power spectrum.")

    parser.add_argument('--mod', dest='qpo_mod', choices=['G', 'g', 'L', 'l'], \
            required=False, default='L', help="Function for QPO model: L for "\
            "Lorentzian, G for Gaussian. [L]")

    parser.add_argument('--prefix', dest='prefix', required=False, default='--',
            help="The identifying prefix of the data (object nickname or "\
            "proposal ID). [--]")

    args = parser.parse_args()

    qpo_mod = args.qpo_mod.upper()

    main(args.infile, qpo_mod, args.prefix)

################################################################################
