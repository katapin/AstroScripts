#!/usr/bin/python
#
#Make and combine light curves



import sys
import os
import argparse
import math
from math import sqrt
import numpy as np
from myfits import *


#Store statistics of current interval
class IntStat:
    num_good = 0                #Number of good points
    num_tot  = 0                #Total number of points 
    rate_mean   = 0             #Mean of the rate column
    rate_std    = 0             #Real standard deviation of flux
    err_mean    = 0             #Mean of the error column
    errsq_mean  = 0             #Mean square of error
    Frms        = 0             #Fractional rms variability
    Frms_err    = 0             #Theoretical arror of Frms
    
    def __init__ (self,hdu):
        Rate=hdu.data['RATE1']
        RError=hdu.data['ERROR1']
        
        ng=np.count_nonzero(~np.isnan(Rate)) 
        rm=np.nanmean(Rate)          
        rs=np.nanstd(Rate)            
        em=np.nanmean(RError)   
        e2m=np.nanmean(RError*RError)        
        rs2=rs**2
        rm2=rm**2
        
        if rs2>e2m: 
            Frms = sqrt(rs2-e2m)/rm  
            self.Frms_err = sqrt( (sqrt(2/ng)*e2m/rm2)**2 + \
            (2*sqrt(e2m/ng)*Frms/rm)**2 )/(2*Frms)
            self.Frms     = Frms           
        else:
            #Frms = sqrt(e2m-rs2)/rm  
            #self.Frms_err = sqrt( (sqrt(2/ng)*e2m/rm2)**2 + \
            #(2*sqrt(e2m/ng)*Frms/rm)**2 )/(2*Frms)
            #self.Frms     = -Frms   
            self.Frms_err = 0
            self.Frms_err = 0
        
        self.num_good       = ng
        self.num_tot        = len(Rate)               
        self.rate_mean      = rm
        self.rate_std       = rs
        self.err_mean       = em
        self.errsq_mean     = e2m
         
 
#Total statistincs of the light curve
class TotalStat:
    arIntervals     = []    #Array of intervals (IntStat objects)
    nInt            = 0     #Number of intervals (after filtering)

    npRates         = []    #Numpy array of rate
    npStd           = []    #Numpy array of standard deviations
    npErr           = []    #Array of measurement errors
    npFrms          = []    #Array of Frms

    aRate           = 0     #Average rate
    aStd            = 0     #Avarage std
    aExpStd         = 0     #Avarage expected std
    aErr            = 0     #Avegare ....
    aFrms           = 0
    aFrmsErr        = 0
    aNumGood        = 0     #Number of good points
    aNumTot         = 0     #Total nubner of points
    
    eaRate          = 0     #Error of average rate
    eaStd           = 0     #Error of average std 
    eaExpStd        = 0     #Error of ...
    eaErr           = 0
    eaFrms          = 0
    eaFrmsErr       = 0
    
    asqFexpStd      = 0     #Average square of fractional expected std: < <sig^2>/<x>^2 >
    asqFrms         = 0     #Average square of Frms
    asqFrmsErr      = 0     #Average square of Frms theoretical error
    
    def  __init__ (self, ftslc):
        nint = len(ftslc) - 1          #Primary HDU is not a light curve
        ar = []
        for i in range(1,nint+1):      #Primary HDU is not a light curve
            curi = IntStat(ftslc[i])
            ar.append(curi)
            
        #Make list of attibutes
        rate=np.array([curi.rate_mean for curi in ar])
        std=np.array([curi.rate_std for curi in ar])
        err=np.array([curi.err_mean for curi in ar])
        expstd=np.array([sqrt(curi.errsq_mean) for curi in ar])
        Frms=np.array([curi.Frms for curi in ar])
        Frms_err=np.array([curi.Frms_err for curi in ar])
        numg=np.array([curi.num_good for curi in ar])
        numt=np.array([curi.num_tot for curi in ar])

        self.npRates    = rate
        self.npStd      = std
        self.npErr      = err
        self.npFrms     = Frms

        
        self.arIntervals    = ar
        self.nInt           = nint
        self.aRate      = rate.mean()
        self.eaRate     = rate.std()/sqrt(nint)
        self.aStd       = std.mean()
        self.eaStd      = std.std()/sqrt(nint)
        self.aErr       = err.mean()
        self.eaErr      = err.std()/sqrt(nint)
        self.aExpStd    = expstd.mean()
        self.eaExpStd   = expstd.std()/sqrt(nint)
        self.aFrms      = Frms.mean()
        self.eaFrms     = Frms.std()/sqrt(nint)
        self.aFrmsErr   = Frms_err.mean()
        self.eaFrmsErr  = Frms_err.std()/sqrt(nint)
        self.aNumGood   = numg.mean()
        self.aNumTot    = numt.mean()
        
        self.asqFexpStd  = np.mean((expstd/rate)**2)
        self.asqFrms     = np.mean(Frms*Frms)
        self.asqFrmsErr  = np.mean(Frms_err*Frms_err)

        
    def displayTable(self):
        i=0
        print("+--------+--------------+----------+----------+----------+-------------+-----------------------+")
        print("| number |   points\t|   rate   |   std    |  <sig>   |sqrt(<sig^2>)|    Frms  +/- error    |")
        for curi in self.arIntervals:
            i=i+1
            print("| #%d \t | %d/%d\t| %f | %f | %f |  %f   | %f +/- %f |" % \
            (i, curi.num_good, curi.num_tot, curi.rate_mean, curi.rate_std, \
            curi.err_mean, sqrt(curi.errsq_mean), curi.Frms, curi.Frms_err ))  
        print("+--------+--------------+----------+----------+----------+-------------+-----------------------+")
        print("|average | %d/%d\t| \033[92m%f\033[0m | %f | %f |  %f   | \033[92m%f\033[0m  |  %f |" % \
        (self.aNumGood, self.aNumTot, self.aRate, self.aStd, self.aErr, 
        self.aExpStd, self.aFrms, self.aFrmsErr))
        print("| error  |    -/-   \t| %f | %f | %f |  %f   | \033[92m%f\033[0m  |  %f |" % \
        (self.eaRate, self.eaStd, self.eaErr, self.eaExpStd, self.eaFrms, self.eaFrmsErr))
        print("+--------+--------------+----------+----------+----------+-------------+-----------------------+")
    

def writegtixron(filename,thresh_newbin,thresh_interv):
    with open(filename, "w") as gtitxt:
        gtitxt.write(" 3 Windows in this < Xronos Window File >\n"
        " 0 Time Wind.: start       stop  (days)\n"
        " 0 Phase Wind.: epoch  period  (days)/ start stop (0->1) phases    max   10\n"
        " 0 Ints. Wind. for Orig. Bins in Series 1 : min  max (c/s)         max   10\n"
        " 0 Ints. Wind. for New Bins   in Series 1 : min  max (c/s)         max   10\n"
        " 0 Ints. Wind. for Intervals  in Series 1 : min  max (c/s)         max   10\n"
        " 1 Exps. Wind. for Orig. Bins in Series 1 : min  max (0->50)       max    1\n"
        "     0.00\t\t50.0     1 \n"
        " 1 Exps. Wind. for New Bins   in Series 1 : min  max (0->50)       max    1\n")
        gtitxt.write("     %s\t\t50.0     1\n" % str(thresh_newbin/100))
        gtitxt.write(" 1 Exps. Wind. for Intervals  in Series 1 : min  max (0->50)       max    1\n"
        "     %s\t\t50.0     1\n" % str(thresh_interv/100))
        
   
def calculatelcrms(lcurve, dt, intervlen, thresh_newbin, thresh_interv):
    pid=os.getpid()                 #PID for temporary files
    tmpgtifile=".%s.wi" % pid       #GTI in XRONOS format
    tmplcfile=".%s.flc" % pid       #Rebined and splitted light curve
    tmppdsfile=".%s.fps" % pid      #Power spectrum for direct inegratation
    
    #Create GTI
    writegtixron(tmpgtifile,thresh_newbin,thresh_interv)
    
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
    lcstat = TotalStat(ftslc)
    ftslc.close()

    #Compute power spectrum 
    COMMAND="powspec cfile1=./%s window=./%s dtnb=%f nbint=%d nintfm=%d " \
    "outfile=%s rebin=-1.1 plot=no fast=no norm=-2 >/dev/null" % (lcurve, \
    tmpgtifile, dt, intervlen, lcstat.nInt, tmppdsfile)
    if not callftools(COMMAND):
        printerr("Something is going wrong: 'powspec' finished with"
        " error.")
        return
        
    ftsps = fitsopen(tmppdsfile)
    X=ftsps[1].data['FREQUENCY']
    Y=ftsps[1].data['POWER']
    Yerr=ftsps[1].data['ERROR']
    
    PdsI=np.trapz(Y,X)      #Integrated power spectrum
    Frms_pds=sqrt(PdsI) if PdsI > 0 else 0
    
    print()
    printbold("%d intervals processed:" % lcstat.nInt);    
    lcstat.displayTable()
    print()
    print("Fractional rms variability in \033[01m%1.4e - %f Hz\033[0m range:" % 
    (1/(dt*intervlen), 0.5/dt))
    print("    From\t  Frms      empirical error  theoretical error")
    print("light curve:\t%f\t%f\t%f" % (lcstat.aFrms, lcstat.eaFrms, 
    sqrt(lcstat.asqFrmsErr/lcstat.nInt)))
    print("power spectrum:\t%f" % (Frms_pds))
    print("Flux (cnt/s):\t%f\t%f\t%f" % (lcstat.aRate, lcstat.eaRate, lcstat.aStd/sqrt(lcstat.nInt*lcstat.aNumGood)))
    
    os.remove(tmplcfile)
    os.remove(tmpgtifile)
    os.remove(tmppdsfile)
    

if __name__ == '__main__':    
    
    #Parse the arguments
    parser = argparse.ArgumentParser(description="Calculate freactional "
    "rms variability of light curve")
    parser.add_argument('lcurve', nargs=1, help="light curve")
    parser.add_argument('fmin', nargs=1, type=float, help="low frequency")
    parser.add_argument('fmax', nargs=1, type=float, help="high frequency")
    parser.add_argument('-n', '--newbin-threshold', nargs='?', type=int,
    help="Newbin threshold (percents), defauld value is 80")
    parser.add_argument('-i', '--interval-threshold', nargs='?', type=int,
    help="Interval threshold (percents), defauld value is 50")
    
    parser.add_argument('-v', '--verbose', action='store_true',
    help="Print debug information")
    
    argnspace=parser.parse_args(sys.argv[1:])
    lcurve=argnspace.lcurve[0]
    fmin=argnspace.fmin[0]
    fmax=argnspace.fmax[0]
    thresh_newbin = argnspace.newbin_threshold
    thresh_interv=argnspace.interval_threshold
    
    if thresh_newbin is None: thresh_newbin=80 
    if thresh_interv is None: thresh_interv=50
    
    #############################################################
    #Check input data
    
	## Match FITS filename and HDU (filename)([HDU])
	#match=re.match("(.*?)(\[(\d+)\])?", lcurve)
	#if not match:
		#die("Incorrect syntax","lcrms")
	#lcname=match.group(1)
    #lchdu = int(match.group(3)) if match.group(2) else 1
        
	##Try to open the file 
	#fts = fitsopen(lcname)
	#if len(fts)<=lchdu:
		#die("There is no HDU number '%d' in %s" % (lchdu, lcname),
		#"lcrms") 
   
    fts=fitsopen(lcurve)    
    #if not chkislc(fts[1]):
        #die("'%s' is not a valid FITS light" % lcurve, "lcrms")
      
    #Check fmax to be less than Nyquist frequency
    dt_orig=float(fts[1].header['TIMEDEL'])     #Original time resolution
    dt_fmax=1/(2*fmax)                          #Required time resolution
    if dt_orig<=dt_fmax:                        #If original better
        dt_new=dt_fmax
    else:
        dt_new=dt_orig
        printwarn("'fmax' will be decreased to %s" % str(1/(2*dt_orig)), "lcrms")
        
    
    #Check fmin
    Time=fts[1].data['TIME']
    timelen=Time[-1]-Time[0];
    if timelen<(1/fmin):
        fmin=1/timelen
        printwarn("'fmin' will be increased to %s" % fmin, "lcrms")
    fts.close()
        
    ##################################################    
    intervlen=math.ceil((1/fmin)/dt_new)
    if calculatelcrms(lcurve, dt_new, intervlen, thresh_newbin, thresh_interv):
        die("Cannot process light curve","lcrms")

                            
