#!/usr/bin/python
#
#Calculate Scargle/Lomb power specrtra for 


import sys
import os
import re
import argparse
from mypython import * 
from myfits import *
from astropy.io import fits
import matplotlib.pyplot as plt

def lspower(infile, outfile, timecol, ratecol, rebin):
    fts = fitsopen(infile)
    time=fts[1].data[timcol]
    rate=fts[1].data[timcol]
    return True


if __name__ == '__main__':    
    
    #Parse the arguments
    parser = argparse.ArgumentParser(description="Plot Good Time "
    "Intervals (GTI) as П-like functions")
    parser.add_argument('-t', '--time', nargs='?', help="Time column")
    parser.add_argument('-r', '--rate', nargs='?', help="Rate column")
    parser.add_argument('-b', '--rebin', nargs='?', help="Rebin parameter")
    parser.add_argument('infile', nargs=1, help="FITS light curve")
    parser.add_argument('outfile', nargs=1, help="FITS power spectrum")
    
    argnspace=parser.parse_args(sys.argv[1:])
    infile=argnspace.infile[0]
    outfile=argnspace.outfile[0]
    timecol=argnspace.time 
    ratecol=argnspace.rate
    rebin=argnspace.rebin
    
    if not os.path.isfile(infile):
        die("File '%s' is not found" % infile)
        
    if not os.path.isfile(outfile):
        die("File '%s' already exists" % outfile)
        
    if not timecol:
        timecol="TIME"
    if not ratecol:
        ratecol="RATE"
    if not rebin:
        rebin=1.1

    if not lspower(infile, outfile, timecol, ratecol, rebin):
        die("Cannot perform FITS-file","lspower")
    
