#!/usr/bin/python
#
#EVENT-file filtering script. 



import sys
import os
import argparse
from astropy.io import fits
import matplotlib.pyplot as plt
from mypython import * 
from myfits import *
import gtiplot
from xmmgeneral import *
                
def xmmflt(evtfile,NROOT):
    
    ftsevt=fitsopen(evtfile)
    evtinf=xmmgetevtinfo(ftsevt)
    if evtinf.datamode!='IMAGING':
        prinitwarn("'%s' is in '%s' mode. Unfortunately only 'IMAGING' " 
        "datamode supported." % (evtfile, evtinf.datamode))
        return False
        
    printgreen("1)Basic filtering")
    evtflt=NROOT+"_flt.evt" 
    COMMAND="evselect table='%s' withfilteredset=Y filteredset='%s' "\
    "destruct=Y keepfilteroutput=T "\
    "expression='%s'" % (evtfile, evtflt, evtinf.basic_flt_exp)
    if not callxmm(COMMAND,logfile):
        printerr("Cannot perform basic filtering: 'evselect' finished "
        "with error.")
        return False
    if not os.path.isfile(evtflt):
        printerr("Something is going wrong: %s is not "
        "created." % evtflt)
        return False
    else:
        printbold("Saved "+evtflt)
    
    #Making flash light curve
    printgreen("2)Making flash light curve")
    flshlc="lcflsh_"+NROOT+".fts"
    COMMAND="evselect table='%s' withrateset=Y rateset='%s' "\
    "maketimecolumn=Y timebinsize=100 makeratecolumn=Y "\
    "expression='%s'" % (evtflt, flshlc, evtinf.flash_flt_exp)
    if not callxmm(COMMAND,logfile):
        printerr("Cannot perform flash filtering: 'evselect' finished "
        "with error.")
        return False
    if not os.path.isfile(flshlc):
        printerr("Something is going wrong: %s is not "
        "created." % flshlc)
        return False
    else:
        printbold("Saved "+flshlc)
    
    #Making flash GTI-file
    printgreen("3)Making flash GTI-file")
    if evtinf.submode not in ['PrimeFullWindow', 'PrimeFullWindowExtended']:
        printwarn("Warning: Observation mode is not a FullFrame. "
        "The automatic count rate criterion may underestimate the significance "
        "of some faint flashes.")
    gtiflsh="gtiflsh_"+NROOT+".fts"
    COMMAND="tabgtigen table='%s' "\
    "expression='%s' gtiset='%s'" % (flshlc, evtinf.flash_gti_exp, gtiflsh)
    if not callxmm(COMMAND,logfile):
        printerr("Cannot create GTI file: 'tabgtigen' finished "
        "with error.")
        return False
    if not os.path.isfile(gtiflsh):
        printerr("Something is going wrong: %s is not "
        "created." % gtiflsh)
        return False
    else:
        printbold("Saved "+gtiflsh)
    
    #Making FITS-image
    printgreen("4)Making images")
    imgevt="img_"+NROOT+".fts"
    if not xmmimagemake(evtflt, imgevt,"X","Y"):
        printerr("Cannot create a FITS image.")
        return False
    printbold("Saved "+imgevt)
    imgevtpng="img_"+NROOT+".png"
    if fitstoimg(imgevt,imgevtpng):
        printbold("Saved "+imgevtpng)
    else:
        printwarn("Cannot convert FITS image to png.")
    
        
    if evtinf.instrkey=='EPN':
        imgevtdetxy="img_"+NROOT+"_detxy.fts"
        imgevtdetxypng="img_"+NROOT+"_detxy.png"
        if xmmimagemake(evtflt, imgevtdetxy,"DETX","DETY"):
            if fitstoimg(imgevtdetxy,imgevtdetxypng):
                printbold("Saved "+imgevtdetxypng)
                os.remove(imgevtdetxy)
            else:
                printwarn("Cannot convert a DETXY FITS image to png.")
        else:
            printerr("Cannot create DETXY FITS image.")
     
    imgflshpng="imgflsh_"+NROOT+".png"
    if gtiplot.gtiplot([gtiflsh],[1],flshlc,imgflshpng,flshlc):
        printbold("Saved "+imgflshpng)
    else:
        printwarn("Cannot save flash light curve as png image")

    #Clean unnecessary files and make filelist
    #flist=open(FILELIST,'a')
    #flist.write(evtflt+"\n")
    #flist.write(gtiflsh+"\n")
    #flist.write(imgevt+"\n")
    #flist.write(imgevtpng+"\n")
    #flist.write(imgflshpng+"\n")
    return True


if __name__ == '__main__':    
    
    #Parse the arguments
    parser = argparse.ArgumentParser(description="Perform basic "
    "filtering of XMM-Newton EVENT-files")
    parser.add_argument('evtfile', nargs=1, help="name of the "
    "raw XMM-Newton EVENT-file")
    parser.add_argument('-u', '--suffix', nargs='?', help="name suffinx "
    "for output files")
    parser.add_argument('-l', '--logfile', nargs='?', help="log file")
    argnspace=parser.parse_args(sys.argv[1:])
    evtfile=argnspace.evtfile[0]
    sufx=argnspace.suffix
    logfile=argnspace.logfile

    #Check the system variables
    if ("SAS_ODF" not in os.environ) or ("SAS_ODF" not in os.environ):
        die("'SAS_ODF' and 'SAS_CCF' variables are not defined. Please \
    define the variables and try again.")
    
    #Open EVENT-file
    ftsevt=fitsopen(evtfile)
    if not xmmchkisevt(ftsevt):
        die("'%s' is not a valid EVENT-file" % evtfile)
    
    #Check the suffix
    evtinf=xmmgetevtinfo(ftsevt)
    if sufx:
        NROOT=evtinf.instr_short_name+"-"+sufx   #Root of the filenames
    else:
        NROOT=evtinf.instr_short_name
    
    evtflt=NROOT+"_flt.evt"      #Check the filtered EVENT-file
    ftsevt.close()
    if os.path.isfile(evtflt):
        die("File %s already exist. Please change the name suffix or clear the folder." % evtflt)
 
    printcaption("%s processing..." % evtinf.instr_full_name)
    if not xmmflt(evtfile,NROOT):
        die("Cannot process '%s'" % evtfile)
    printcaption("Finished")
    print("Data mode:       %s" % evtinf.datamode)
    print("Data submode:    %s" % evtinf.submode)
    print("Filter:          %s" % evtinf.filter)
    print("Exposure ID:     %s" % evtinf.expidstr)
    print("Exposure length: %s" % evtinf.telapse)
    

