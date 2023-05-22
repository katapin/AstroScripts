#!/usr/bin/python
#
#Plot GTI files

import sys
import os
import argparse
import numpy as np
from astropy.io import fits
from mypython import *
import csv


#Settings
pardef="defpar.dat"
enerdef="defener.dat"

    
if __name__ == '__main__':    
    
        
    #Parse the arguments
    parser = argparse.ArgumentParser(description="Pack XSPEC"
    "table model")
    parser.add_argument('modelname', nargs=1, help="Table model name")
    parser.add_argument('gridpath', nargs=1, help="Path to calculated model spectra")
    
    argnspace=parser.parse_args(sys.argv[1:])
    modelname=argnspace.modelname[0]           #Model name
    gridpath=argnspace.gridpath[0]             #Path
    
    
    ##
    ################  HDU with parameters  ###########################
    pardefpath=gridpath+"/"+pardef
    if not os.path.isfile(pardefpath):
        die("File '%s' is not found" % pardef)
        
    #Read defenition file
    with open(pardefpath) as csvfile:
        reader = csv.reader(csvfile)
        cols=list(zip(*reader))
        
    NVarPar=len(cols[0]);       #Number of model parameters
    arparlen=np.array([int(x) for x in cols[8]])
    parmaxlen=max(arparlen)     #Maximal length of parameter value list
    NSpec=arparlen.prod()       #Number of spectra
    col1 = fits.Column(name='NAME', format='12A', array=cols[0])
    col2 = fits.Column(name='METHOD', format='J', array=cols[1])
    col3 = fits.Column(name='INITIAL', format='E', array=cols[2])
    col4 = fits.Column(name='DELTA', format='E', array=cols[3])
    col5 = fits.Column(name='MINIMUM', format='E', array=cols[4])
    col6 = fits.Column(name='BOTTOM', format='E', array=cols[5])
    col7 = fits.Column(name='TOP', format='E', array=cols[6])
    col8 = fits.Column(name='MAXIMUM', format='E', array=cols[7])
    col9 = fits.Column(name='NUMBVALS', format='J', array=cols[8])

    
    #For each parameter file
    arparvals=[]
    for i in range(1,NVarPar+1):
        filepath=gridpath+"/defpar%ival.dat" % i
        if not os.path.isfile(filepath):
            die("File '%s' is not found" % filepath)
        with open(filepath) as parfile:
            parvals=[float(line) for line in parfile]
        parvals.extend([0.0] * (parmaxlen-len(parvals)))
        arparvals.append(parvals)
    col10 = fits.Column(name='VALUE', format="%dE" % parmaxlen, array=arparvals)

    #Create HDU
    FitsTabCols = fits.ColDefs([col1, col2, col3, col4, col5, col6, col7,
    col8, col9, col10])
    
    tbhdu1 = fits.BinTableHDU.from_columns(FitsTabCols,name="PARAMETERS")
    
    tbhdu1.header.set("NINTPARM","%d" % NVarPar)
    tbhdu1.header.set("NADDPARM","0")
    
    tbhdu1.header.set("HDUCLASS","OGIP")
    tbhdu1.header.set("HDUCLAS1","XSPEC TABLE MODEL")
    tbhdu1.header.set("HDUCLAS2","PARAMETERS")
    tbhdu1.header.set("HDUVERS","1.0.0")
    
    
    
    ################### HDU with energies ##############################
    enerdefpath=gridpath+"/"+enerdef
    if not os.path.isfile(enerdefpath):
        die("File '%s' is not found" % enerdef)

        
    #Read energies
    with open(enerdefpath) as csvfile:
        reader = csv.reader(csvfile)
        cols=list(zip(*reader))
        
    col1 = fits.Column(name='ENERG_LO', format='E', array=cols[0])
    col2 = fits.Column(name='ENERG_HI', format='E', array=cols[1])
    
    ener_lo=np.array([float(x) for x in cols[0]])
    ener_hi=np.array([float(x) for x in cols[1]])
    ener_delta=ener_hi-ener_lo
    #ener=0.5(ener_hi+ener_lo)
    
    #Create HDU
    FitsTabCols = fits.ColDefs([col1, col2])
    tbhdu2 = fits.BinTableHDU.from_columns(FitsTabCols,name="ENERGIES")
    
    tbhdu2.header.set("HDUCLASS","OGIP")
    tbhdu2.header.set("HDUCLAS1","XSPEC TABLE MODEL")
    tbhdu2.header.set("HDUCLAS2","ENERGIES")
    tbhdu2.header.set("HDUVERS","1.0.0")
    
    
    
    ################### HDU with spectra ##############################
    arparvals=[]
    arspecs=[]
    for i in range(1,NSpec+1):
        specfilepath=gridpath+"/spec%i.dat" % i
        parfilepath=gridpath+"/par%i.dat" % i
        if not os.path.isfile(specfilepath):
            die("File '%s' is not found" % specfilepath)
        if not os.path.isfile(parfilepath):
            die("File '%s' is not found" % parfilepath)
            
        with open(parfilepath) as parfile:
            parvals=[float(x) for x in parfile.readline().split()]
        with open(specfilepath) as specfile:
            #Integrate spectrum inside bins for XSPEC
            spec=np.array([float(line) for line in specfile])*ener_delta*2
        arparvals.append(parvals)
        arspecs.append(spec)
    SpecLen=len(spec)
    
    col1 = fits.Column(name='PARAMVAL', format="%dE" % NVarPar, array=arparvals)
    col2 = fits.Column(name='INTPSPEC', format='%dE' % SpecLen, unit="photons/cm^2/s", 
    array=arspecs)
    
    FitsTabCols = fits.ColDefs([col1, col2])
    tbhdu3= fits.BinTableHDU.from_columns(FitsTabCols,name="SPECTRA")
    
    tbhdu3.header.set("HDUCLASS","OGIP")
    tbhdu3.header.set("HDUCLAS1","XSPEC TABLE MODEL")
    tbhdu3.header.set("HDUCLAS2","MODEL SPECTRA")
    tbhdu3.header.set("HDUVERS","1.0.0")
        
    
    
    ################### Writhe the FITS file ###########################
    fitsname=modelname+".fts"
    hdu = fits.PrimaryHDU()       
        
    hdu.header.set("MODLNAME",modelname)
    hdu.header.set("MODLUNIT","photons/cm^2/s")
    hdu.header.set("REDSHIFT",False)
    hdu.header.set("ADDMODEL",True)
    hdu.header.set("HDUCLASS","OGIP")
    hdu.header.set("HDUCLAS1","XSPEC TABLE MODEL")
    hdu.header.set("HDUVERS","1.0.0")
    
    
    hdulist = fits.HDUList([hdu, tbhdu1, tbhdu2, tbhdu3])
    hdulist.writeto(fitsname)


