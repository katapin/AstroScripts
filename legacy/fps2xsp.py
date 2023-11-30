#!/usr/bin/python


import sys
import os
import re
import argparse
from mypython import * 
from myfits import *
from astropy.io import fits
import matplotlib.pyplot as plt

def fps2xsp(infile, specfile, rspfile, noiselevel=0):

    #Assumed first extension 
    fts = fitsopen(infile)
    X=fts[1].data['FREQUENCY']
    Xerr=fts[1].data['XAX_E']
    Y=fts[1].data['POWER']-noiselevel
    Yerr=fts[1].data['ERROR']
    fts.close()
    
    COL1=X-Xerr
    COL2=X+Xerr
    COL3=2*Y*Xerr
    COL4=2*Yerr*Xerr
    
    tmpascii=".tmpxspec.dat"        #Save columns as text
    tmprsp=".tmpresp.fts"
    with open(tmpascii, 'w') as ftxt:
        for i in range(0,len(COL1)):
            ftxt.write("%e %e %e %e\n" % (COL1[i], COL2[i], COL3[i], COL4[i]))
    
    cmd="flx2xsp '%s' '%s' '%s'" % (tmpascii, specfile, tmprsp)
    if not callftools(cmd):
        printerr("Can't create XSPEC fits from ascii file","fps2xsp")
        return False
        
    #Update header keywords
    fxsp=fits.open(specfile,mode='update')
    fxsp[0].header['creator']="fps2xsp.py"
    fxsp[0].header['history']="Created by fps2xsp.py from %s" % infile
    if noiselevel:
        fxsp[0].header['history']="Poisson noise level %f has been subtracted" % noiselevel
    fxsp.flush()
        
    #Embed RMF as second and third HDU extension or just make separate FITS 
    if rspfile:
        os.rename(tmprsp,rspfile)
    else:
        frsp=fitsopen(tmprsp)
        fxsp.append(frsp[1])
        fxsp.append(frsp[2])
        fxsp[1].header['RESPFILE']=specfile
        fxsp.flush()
        os.remove(tmprsp)
        
    fxsp.close()
    
    os.remove(tmpascii)
    return True


if __name__ == '__main__': 

    #Parse the arguments
    parser = argparse.ArgumentParser(description="Converts power spectrum from powspec into XSPEC format.")
    parser.add_argument('infile', nargs=1, help="FITS file with "
    "power spectrum")
    parser.add_argument('outspec', nargs='?', help="Name of output"
    "xpectrum FITS suitable of XSPEC")
    parser.add_argument('outresp', nargs='?', help="Assign RMF filename, "
    "otherwise it will be builded into the spectrum FITS.")
    parser.add_argument('-s', '--subtract-noise', nargs='+', help="Subtract"
    "Poisson noise level", type=float)
    
    argnspace=parser.parse_args(sys.argv[1:])
    infile=argnspace.infile[0]
    
    #Is Poisso noise level defined?
    noiselevel=0
    namesuff=""
    if argnspace.subtract_noise:
        noiselevel=argnspace.subtract_noise[0]
        namesuff="-s%4.3f" % noiselevel     #name suffix
        
    #Is outroot defined?    
    if not argnspace.outspec:
        xspfile=os.path.splitext(infile)[0]+namesuff+".xsp"
    else:
        xspfile=argnspace.outspec
        
    rspfile=None
    if argnspace.outresp:
       rspfile=argnspace.outresp
    
    
    fitscheck(infile)             #Check fits power spectrum
    
    if not fps2xsp(infile,xspfile,rspfile,noiselevel):
       die("Cannot convert FITS-file","fps2xsp")
    
    
        

    
