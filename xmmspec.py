#!/usr/bin/python
#
#Make spectra

import sys
import os
import argparse
from astropy.io import fits
from mypython import * 
from myfits import *
from xmmgeneral import *

def extractspectrum(evtfile, gtifile, regfile, specfile):
    regexpr=readregfile(regfile)
    if not regexpr:
        printerr("Empty region file '%s' or unsupported format. Only "
        "ds9 region format is supported." % regfile)
        return False
    
    expression="((%s) && gti(%s,TIME) && (FLAG==0))" % \
    (regexpr, gtifile)
    
    #Get evt-info
    evtinf=xmmgetevtinfo(evtfile)
    if evtinf.instrkey == 'EPN':
        #For EPIC-PN spectra
        COMMAND="evselect table='%s' energycolumn=PI " \
        "expression='%s' withspectrumset=yes spectrumset='%s' " \
        "spectralbinsize=5 withspecranges=yes specchannelmin=0 " \
        "specchannelmax=20479" % (evtfile, expression, specfile)
    else:
        #For EPIC-MOS spectra
        COMMAND="evselect table='%s' energycolumn=PI " \
        "expression='%s' withspectrumset=yes spectrumset='%s' " \
        "spectralbinsize=5 withspecranges=yes specchannelmin=0 " \
        "specchannelmax=11999" % (evtfile, expression, specfile)
    
    if not callxmm(COMMAND,logfile):
        printerr("Something is going wrong: 'evselect' finished with"
        " an error.")
        return False
    if not os.path.isfile(specfile):
        printerr("Something is going wrong: %s is not "
        "created." % specfile)
        return False
    
    #Backscale    
    COMMAND="backscale spectrumset='%s' badpixlocation='%s'" \
    % (specfile, evtfile)
    if not callxmm(COMMAND,logfile):
        printerr("Something is going wrong: 'backscale' finished with"
        " an error.")
        return False
        
    #Make images
    tmpimgfts="tmpimg_"+specfile
    tmpimgpng="tmpimg_"+specfile+".png"
    if xmmimagemake(evtfile,tmpimgfts,"X","Y",expression):
        if fitstoimg(tmpimgfts,tmpimgpng):
            printbold("Saved "+tmpimgpng)
            os.remove(tmpimgfts)
        else:
            printwarn("Cannot convert a FITS image to png.")
    else:
        printerr("Cannot create FITS image.")
    printbold("Saved "+specfile)
    return True

def makespectrum(evtfile, gtifile, objreg, bkgreg, spobj, spbkg,
rmfname, arfname):
    
    printgreen("Processing '%s' -> '%s'" % (evtfile, spobj))
    if not extractspectrum(evtfile, gtifile, objreg, spobj):
        printerr("Cannot extract object spectrum '%s'" % spobj)
        return False
        
    printgreen("Processing '%s' -> '%s'" % (evtfile, spbkg))
    if not extractspectrum(evtfile, gtifile, bkgreg, spbkg):
        printerr("Cannot extract background lightcurve '%s'" % lcbkg)
        return False
    
    #Make RMF
    printgreen("Generating RMF '%s 'for '%s'" % (spobj, rmfname))
    COMMAND="rmfgen spectrumset='%s' rmfset='%s'" % (spobj, rmfname)
    if not callxmm(COMMAND,logfile):
        printerr("Something is going wrong: 'rmfgen' finished with"
        " anerror.")
        return False
    if not os.path.isfile(rmfname):
        printerr("Something is going wrong: %s is not "
        "created." % lcnet)
        return False
    printbold("Saved "+rmfname)
        
    #Make ARF    
    printgreen("Generating RMF '%s 'for '%s'" % (spobj, arfname))
    COMMAND="arfgen spectrumset='%s' arfset='%s' withrmfset=yes " \
    "rmfset='%s' badpixlocation='%s' detmaptype=psf" % (spobj, arfname,
    rmfname, evtfile)
    if not callxmm(COMMAND,logfile):
        printerr("Something is going wrong: 'arfgen' finished with"
        " an error.")
        return False
    if not os.path.isfile(arfname):
        printerr("Something is going wrong: %s is not "
        "created." % arfname)
        return False
    printbold("Saved "+arfname)
            
    return True


if __name__ == '__main__':    
    
    #Parse the arguments
    parser = argparse.ArgumentParser(description="Accumulate spectrum, " \
    "make RMF and ARF.")
    parser.add_argument('evtfile', nargs=1, help="list of the "
    "XMM-Newton EVENT-files")
    parser.add_argument('gtifile', nargs=1, help="name of the GTI-file")
    parser.add_argument('regobj', nargs=1, help="name of the object"
    "region file")
    parser.add_argument('regbkg', nargs=1, help="name of the background"
    "region file")
    parser.add_argument('binmin', nargs=1, type=float, help="counts "
    "per grouped channel")
    parser.add_argument('-u', '--suffix', nargs='?', help="name suffinx "
    "for output files")
    parser.add_argument('-l', '--logfile', nargs='?', help="log file")
    argnspace=parser.parse_args(sys.argv[1:])
    
    evtfile=argnspace.evtfile[0]
    gtifile=argnspace.gtifile[0]
    regobj=argnspace.regobj[0]
    regbkg=argnspace.regbkg[0]
    binmin=argnspace.binmin[0]
    
    logfile=argnspace.logfile
    
    if argnspace.suffix:
        sufx="-"+argnspace.suffix
    else:
        sufx=""

    
    #############################################################
    #Check input data
        
    #Check the system variables
    if ("SAS_ODF" not in os.environ) or ("SAS_ODF" not in os.environ):
        die("'SAS_ODF' and 'SAS_CCF' variables are not defined. Please \
    define the variables and try again.","xmmflt")
    
    
    ftsevt=fitsopen(evtfile)
    if not xmmchkisevt(ftsevt):
        die("'%s' is not a valid EVENT-file" % evtfile)
        
    evinf=xmmgetevtinfo(ftsevt)
    
    spobj="spobj%s_%s.fts" % (sufx, evinf.instr_short_name)
    spbkg="spbkg%s_%s.fts" % (sufx, evinf.instr_short_name)
    rmffile="rmf%s_%s.fts" % (sufx, evinf.instr_short_name)
    arffile="arf%s_%s.fts" % (sufx, evinf.instr_short_name)
    spgrp="spgrp%s_%s_%d.fts" % (sufx, evinf.instr_short_name,int(binmin))
    imgsp="imgsp%s_%s_%d.ps" % (sufx, evinf.instr_short_name, int(binmin))
    
    if os.path.isfile(spobj):
        die("File %s already exist. Please change the name suffix "
        "or clear the folder." % spobj)
    
    if evinf.datamode!='IMAGING':
        die("'%s' is in '%s' mode. Unfortunately only 'IMAGING' " 
        "datamode supported." % (evtfile, evinf.datamode))
    if not os.path.isfile(regobj):
        die("'%s' is not found" % regobj)
    if not os.path.isfile(regbkg):
        die("'%s' is not found" % regbkg)
    ftsevt.close()
    
    #Check GTI file        
    ftsgti=fitsopen(gtifile)
    if not chkisgti(ftsgti[1]):
        die("'%s' is not a valid GTI extension" % gtifile+"[1]")
    ftsgti.close()
    
    ################################################################
    ##Perform processing
    
    printcaption("Making '%s'..." % spgrp)
    if not makespectrum(evtfile, gtifile, regobj, regbkg, spobj, spbkg, 
    rmffile, arffile):
        die("Some errors arise during EVENT-file '%s' "
        "be processed." % evtfiles[i])
    
    printgreen("Group spectrum %s: min %d counts" % (spobj, int(binmin)))
    COMMAND="specgroup spectrumset='%s' backgndset='%s' rmfset='%s' " \
    "arfset='%s' groupedset='%s' mincounts='%d'" % (spobj,
    spbkg, rmffile, arffile, spgrp, int(binmin))
    if not callxmm(COMMAND,logfile):
        die("Something is going wrong: 'specgroup' finished with"
        " an error.")
    if not os.path.isfile(spgrp):
        die("Something is going wrong: %s is not "
        "created." % lcnet)
    printbold("Saved "+spgrp)
    
    
    #TODO:
    #2)Make ps-image
    #3) pile-up
    
    #if not specplot(lcobj,lcbkg,lcnet,gtifile,imglc):
        #printerr("Can't plot light curve")
            
    printcaption("Finished")
