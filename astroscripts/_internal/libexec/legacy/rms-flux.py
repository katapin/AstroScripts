#!/usr/bin/python
#
#Calculate rms-flux relation

import argparse
import math
import numpy as np
import matplotlib.pyplot as plt
import lcrms
from myfits import *


def perform(lcurve, dt, intervlen, grouping_method, thresh_newbin, thresh_interv):
    pid=os.getpid()                 #PID for temporary files
    tmpgtifile=".%s.wi" % pid       #GTI in XRONOS format
    tmplcfile=".%s.flc" % pid       #Rebined and splitted light curve
    tmppdsfile=".%s.fps" % pid      #Power spectrum for direct inegratation

    #Create GTI
    lcrms.writegtixron(tmpgtifile,thresh_newbin,thresh_interv)

    #Prepare light curve"
    COMMAND="lcurve nser=1 cfile1=./%s window=./%s dtnb=%f nbint=%d " \
    "outfile=%s plot=no >/dev/null" % (lcurve, tmpgtifile, dt, intervlen, tmplcfile)
    if not callftools(COMMAND):
        printerr("Something is going wrong: 'lcurve' finished with"
        " error.")
        return;
    if not os.path.isfile(tmplcfile):
        print("")
        printbold("Failure: There is no valid intervals to process!\n"
        "Try to change frequency range or decrease threshold (-h to see help)\n")
        return

    ftslc = fitsopen(tmplcfile)       #Rebined light curve FITS
    lcstat = lcrms.TotalStat(ftslc)
    ftslc.close()

    allrates    = lcstat.npRates
    allrms      = lcstat.npFrms*allrates

    if np.chararray.isnumeric(grouping_method):
        Npt=int(grouping_method)
        grouping_method=int(grouping_method)
    else:
        Npt=8

    print(grouping_method)

    bin_low_lim = min(allrates) - (max(allrates)-min(allrates))/Npt/10      #Expand bin limits to avoid droping
    bin_hi_lim = max(allrates) + (max(allrates) - min(allrates))/ Npt/10    #datapoint at bin edges

    bins = np.histogram_bin_edges(allrates,grouping_method,[bin_low_lim,bin_hi_lim])   #Calculate bin edges
    Nbins = len(bins)-1

    binnumb = np.digitize(allrates,bins)   #Sift the data and get number of the bin for each interval

    X       = np.zeros(Nbins)  # Data points for plot
    X_err   = np.zeros(Nbins)
    X_bins  = np.zeros(Nbins)
    Y = np.zeros(Nbins)
    Y_err   = np.zeros(Nbins)

    printbold("%d intervals processed:" % lcstat.nInt);

    for i in range(Nbins):
        #Lists of element indices for each bin (bin number is dict key)
        elemind = [j for j, x in enumerate(binnumb) if x == i+1]

        X[i]=allrates[elemind].mean()
        X_err[i] = allrates[elemind].std()
        X_bins[i]=(bins[i]+bins[i+1])/2
        Y[i] = allrms[elemind].mean()
        Y_err[i] = allrms[elemind].std()


    plt.errorbar(X,Y,xerr=X_err,yerr=Y_err,fmt='o')
    plt.errorbar(X_bins,Y,yerr=Y_err,fmt='o')
    plt.xlabel('count rate, cnt/s')
    plt.ylabel('rms (from Frms), cnt/s')
    plt.show()


    # #Compute power spectrum
    # COMMAND="powspec cfile1=./%s window=./%s dtnb=%f nbint=%d nintfm=%d " \
    # "outfile=%s rebin=-1.1 plot=no fast=no norm=-2 >/dev/null" % (lcurve, \
    # tmpgtifile, dt, intervlen, lcstat.nInt, tmppdsfile)
    # if not callftools(COMMAND):
    #     printerr("Something is going wrong: 'powspec' finished with"
    #     " error.")
    #     return
    #
    # ftsps = fitsopen(tmppdsfile)
    # X=ftsps[1].data['FREQUENCY']
    # Y=ftsps[1].data['POWER']
    # Yerr=ftsps[1].data['ERROR']
    #
    # PdsI=np.trapz(Y,X)      #Integrated power spectrum
    # Frms_pds=sqrt(PdsI) if PdsI > 0 else 0
    #
    # print()
    # printbold("%d intervals processed:" % lcstat.nInt);
    # #lcstat.displayTable()
    # print()
    # print("Fractional rms variability in \033[01m%1.4e - %f Hz\033[0m range:" %
    # (1/(dt*intervlen), 0.5/dt))
    # print("    From\t  Frms      empirical error  theoretical error")
    # print("light curve:\t%f\t%f\t%f" % (lcstat.aFrms, lcstat.eaFrms,
    # sqrt(lcstat.asqFrmsErr/lcstat.nInt)))
    # print("power spectrum:\t%f" % (Frms_pds))
    # print("Flux (cnt/s):\t%f\t%f\t%f" % (lcstat.aRate, lcstat.eaRate, lcstat.aStd/sqrt(lcstat.nInt*lcstat.aNumGood)))

    os.remove(tmplcfile)
    os.remove(tmpgtifile)
    # os.remove(tmppdsfile)


if __name__ == '__main__':

    #Parse the arguments
    parser = argparse.ArgumentParser(description="Plot rms-flux relation "
    "for light curve")
    parser.add_argument('lcurve', nargs=1, help="light curve")
    parser.add_argument('Tmax', nargs=1, type=float, help="interval length")
    parser.add_argument('dt', nargs=1, type=float, help="time resolution")
    parser.add_argument('-g', '--grouping_method', nargs='?', type=str,
    help="Grouping method, either number of equal width bins or any another metod supported " \
    "by numpy.histogram. Default value is 8 bins")
    parser.add_argument('-n', '--newbin-threshold', nargs='?', type=int,
    help="Newbin threshold (percents), defauld value is 80")
    parser.add_argument('-i', '--interval-threshold', nargs='?', type=int,
    help="Interval threshold (percents), defauld value is 50")

    parser.add_argument('-v', '--verbose', action='store_true',
    help="Print debug information")

    argnspace=parser.parse_args(sys.argv[1:])
    lcurve=argnspace.lcurve[0]
    Tmax=argnspace.Tmax[0]
    dt_req=argnspace.dt[0]
    grouping_method=argnspace.grouping_method
    thresh_newbin = argnspace.newbin_threshold
    thresh_interv=argnspace.interval_threshold

    if Tmax<=dt_req:
        die("Wrong Tmax")

    if grouping_method is None: grouping_method='8'
    if thresh_newbin is None: thresh_newbin=80
    if thresh_interv is None: thresh_interv=50


    fts=fitsopen(lcurve)
    #if not chkislc(fts[1]):
        #die("'%s' is not a valid FITS light" % lcurve, "lcrms")

    #Check fmax to be less than Nyquist frequency
    dt_orig=float(fts[1].header['TIMEDEL'])     #Original time resolution
    if dt_orig<=dt_req:                        #If original better
        dt_new=dt_req
    else:
        dt_new=dt_orig
        printwarn("'dt' will be increased to %s" % str(dt_new), "rms-flux")


    #Check fmin
    Time=fts[1].data['TIME']
    timelen=Time[-1]-Time[0];
    if timelen<(Tmax):
        Tmax=1/timelen
        printwarn("'Tmax' will be decreased to %s" % Tmax, "lcrms")
    fts.close()


    ##################################################
    intervlen=math.ceil(Tmax/dt_new)
    if perform(lcurve, dt_new, intervlen, grouping_method, thresh_newbin, thresh_interv):
        die("Cannot process light curve","rms-flux")


