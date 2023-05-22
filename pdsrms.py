#!/usr/bin/python

import sys
import os
import re
import argparse
from mypython import * 
from myfits import *
from astropy.io import fits
import matplotlib.pyplot as plt

def pdsrms(infile, fmin, fmax):

    #Assumed first extension 
    fts = fitsopen(infile)
    X=fts[1].data['FREQUENCY']
    Xerr=fts[1].data['XAX_E']
    Y=fts[1].data['POWER']-noiselevel
    Yerr=fts[1].data['ERROR']
    fts.close()
    
    #COL1=X-Xerr
    #COL2=X+Xerr
    #COL3=2*Y*Xerr
    #COL4=2*Yerr*Xerr
    
    #tmpascii=".tmpxspec.dat"        #Save columns as text
    #tmprsp=".tmpresp.fts"
    #with open(tmpascii, 'w') as ftxt:
        #for i in range(0,len(COL1)):
            #ftxt.write("%e %e %e %e\n" % (COL1[i], COL2[i], COL3[i], COL4[i]))
    
    #cmd="flx2xsp '%s' '%s' '%s'" % (tmpascii, specfile, tmprsp)
    #if not callftools(cmd):
        #printerr("Can't create XSPEC fits from ascii file","fps2xsp")
        #return False
        
    ##Update header keywords
    #fxsp=fits.open(specfile,mode='update')
    #fxsp[0].header['creator']="fps2xsp.py"
    #fxsp[0].header['history']="Created by fps2xsp.py from %s" % infile
    #if noiselevel:
        #fxsp[0].header['history']="Poisson noise level %f has been subtracted" % noiselevel
    #fxsp.flush()
        
    ##Embed RMF as second and third HDU extension or just make separate FITS 
    #if rspfile:
        #os.rename(tmprsp,rspfile)
    #else:
        #frsp=fitsopen(tmprsp)
        #fxsp.append(frsp[1])
        #fxsp.append(frsp[2])
        #fxsp[1].header['RESPFILE']=specfile
        #fxsp.flush()
        #os.remove(tmprsp)
        
    #fxsp.close()
    
    #os.remove(tmpascii)
    #return True


if __name__ == '__main__': 

    #Parse the arguments
    parser = argparse.ArgumentParser(description="Calculate fractional rms variance")
    parser.add_argument('infile', nargs=1, help="FITS file with "
    "power spectrum")
    parser.add_argument('f_min', nargs='1', help="Low frequency"
    parser.add_argument('f_max', nargs='1', help="High frequency"

    
    argnspace=parser.parse_args(sys.argv[1:])
    infile=argnspace.infile[0]
    fmin=argnspace.f_min[0]
    fmax=argnspace.f_max[0]

    
    fitscheck(infile)             #Check fits power spectrum
    
    pdsrms(infile, fmin, fmax)
       
    
    
        

    
