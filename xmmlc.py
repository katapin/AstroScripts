#!/usr/bin/python
#
#Make and combine light curves



import sys
import os
import argparse
from astropy.io import fits
from mypython import * 
from myfits import *
from xmmgeneral import *

def extractlightcurve(evtfile, gtifile, regfile, lcfile, timebin, mine, maxe, evtlc=None):
    regexpr=readregfile(regfile)
    if not regexpr:
        printerr("Empty region file '%s' or unsupported format. Only "
        "ds9 region format is supported." % regfile)
        return False
    gtilimits=getgtilimits(gtifile)
    if gtilimits:
        tstart=gtilimits[0]
        tstop=gtilimits[1]
    else:
        printerr("Corrupted GTI file")
        return False
    
    expression="((%s) && gti(%s,TIME) && PI in [%s:%s])" % (regexpr, 
    gtifile, str(int(mine*1000)), str(int(maxe*1000)))
    if evtlc:
        COMMAND="evselect table='%s' energycolumn=PI " \
        "expression='%s' withrateset=yes " \
        "rateset='%s' timebinsize=%s maketimecolumn=yes " \
        "makeratecolumn=yes timemin=%s timemax=%s " \
        "keepfilteroutput=yes withfilteredset=yes filteredset=%s" % (evtfile, expression,
        lcfile, str(timebin), str(tstart), str(tstop), evtlc) 
    else:
        COMMAND="evselect table='%s' energycolumn=PI " \
        "expression='%s' withrateset=yes " \
        "rateset='%s' timebinsize=%s maketimecolumn=yes " \
        "makeratecolumn=yes timemin=%s timemax=%s" % (evtfile, expression,
        lcfile, str(timebin), str(tstart), str(tstop)) 
    if not callxmm(COMMAND,logfile):
        printerr("Something is going wrong: 'evselect' finished with"
        " error.")
        return False
    if not os.path.isfile(lcfile):
        printerr("Something is going wrong: %s is not "
        "created." % lcfile)
        return False
        
    #Make images
    tmpimgfts="tmpimg_"+lcfile
    tmpimgpng="tmpimg_"+lcfile+".png"
    if xmmimagemake(evtfile,tmpimgfts,"X","Y",expression):
        if fitstoimg(tmpimgfts,tmpimgpng):
            printbold("Saved "+tmpimgpng)
            os.remove(tmpimgfts)
        else:
            printwarn("Cannot convert a FITS image to png.")
    else:
        printerr("Cannot create FITS image.")
    return True

def makelightcurves(evtfile, gtifile, objreg, bkgreg, timebin, mine, 
maxe,lcobj,lcbkg,lcnet,evtlc,imglc=""):
    printgreen("Processing '%s' -> '%s'" % (evtfile, lcobj))
    if not extractlightcurve(evtfile, gtifile, objreg, lcobj, timebin, 
    mine, maxe, evtlc):
        printerr("Cannot extract object lightcurve '%s'" % lcobj)
        return False
        
    printgreen("Processing '%s' -> '%s'" % (evtfile, lcbkg))
    if not extractlightcurve(evtfile, gtifile, bkgreg, lcbkg, timebin, 
    mine, maxe):
        printerr("Cannot extract background lightcurve '%s'" % lcbkg)
        return False
    
    printgreen("Processing ('%s' - '%s') -> '%s'" % (lcobj, lcbkg, lcnet))
    COMMAND="epiclccorr srctslist='%s' eventlist='%s' " \
    "outset='%s' withbkgset=yes bkgtslist='%s' " \
    "applyabsolutecorrections=yes" % (lcobj, evtfile, lcnet, lcbkg)
    if not callxmm(COMMAND,logfile):
        printerr("Something is going wrong: 'epiclccorr' finished with"
        " error.")
        return False
    if not os.path.isfile(lcnet):
        printerr("Something is going wrong: %s is not "
        "created." % lcnet)
        return False
    
    #Plot light curves    
    if imglc:
        ftslcnet=fitsopen(lcnet)
        bkgratio=ftslcnet[1].header['BKGRATIO']
        if bkgratio>1: bkgratio=1/bkgratio
        ftslcnet.close()
        if not lcplot(lcobj,imglc,lcbkg,lcnet,gtifile,bkgratio):
            printerr("Can't plot light curve")
            
    return True
                
def combinelc(lc1, lc2, reslc):
    COMMAND="lcmath %s %s %s 1 1 addsubr=yes" % (lc1, lc2, reslc)
    if not callftools(COMMAND):
        printerr("Something is going wrong: 'lcmath' finished with"
        " error.")
        return False
    if not os.path.isfile(reslc):
        printerr("Something is going wrong: %s is not "
        "created." % reslc)
        return False
    return True
    

if __name__ == '__main__':    
    
    #Parse the arguments
    parser = argparse.ArgumentParser(description="Accumulate light curve "
    "with certain time resolution and energy range. Programm can combine "
    "light curves from different detectors. Arguments 'evtfile', 'regobj' "
    "and 'regbkg'  may be lists of files (pn,mos1,mos2) "
    "but sigle GTI file must be common for all detectors.")
    parser.add_argument('evtfile', nargs=1, help="list of the "
    "XMM-Newton EVENT-files")
    parser.add_argument('gtifile', nargs=1, help="name of the GTI-file")
    parser.add_argument('regobj', nargs=1, help="name of the object"
    "region file")
    parser.add_argument('regbkg', nargs=1, help="name of the background"
    "region file")
    parser.add_argument('timebin', nargs=1, type=float, help="time "
    "resolutin of the light curve in seconds")
    parser.add_argument('mine', nargs=1, type=float, help="lower "
    "energy range limit")
    parser.add_argument('maxe', nargs=1, type=float, help="higher "
    "energy range limit")
    parser.add_argument('-u', '--suffix', nargs='?', help="name suffinx "
    "for output files")
    parser.add_argument('-l', '--logfile', nargs='?', help="log file")
    argnspace=parser.parse_args(sys.argv[1:])
    
    evtfiles=argnspace.evtfile[0].split(',')
    gtifile=argnspace.gtifile[0]
    regobjfiles=argnspace.regobj[0].split(',')
    regbkgfiles=argnspace.regbkg[0].split(',')
    timebin=argnspace.timebin[0]
    mine=argnspace.mine[0]
    maxe=argnspace.maxe[0]
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
    
    #Check length of lists
    nfiles=len(evtfiles)
    if nfiles>3:
        die("You are trying to combine more then 3 light curves. Only "
        "synchronous light curves from EPIC-PN, EPIC-MOS1 and EPIC-MOS2 "
        "can be combined.")
    if (len(regobjfiles)!=nfiles) or (len(regbkgfiles)!=nfiles):
        die("Lists of the EVENT and the region files must have the same "
        "length and order.")
        
    #Check energy range
    if mine<0.2 or mine>12:
        die("'mine' must be in the 0.2-12.0 keV range")
    if maxe<0.2 or maxe>12:
        die("'maxe' must be in the 0.2-12.0 keV range")
    if mine>=maxe:
        die("'mine' must be less then 'maxe'")
        
    #Check EVENT-files and region files
    evinf=[]
    instrlst=[]
    lcobjlst=[]
    lcbkglst=[]
    lcnetlst=[]
    evtlclst=[]
    imglclst=[]
    for i in range(nfiles):
        #Open EVENT-file
        ftsevt=fitsopen(evtfiles[i])
        if not xmmchkisevt(ftsevt):
            die("'%s' is not a valid EVENT-file" % evtfiles[i])
            
        evinf.append(xmmgetevtinfo(ftsevt))
        instrlst.append(evinf[i].instr_short_name)
        lcobjlst.append("lcobj%s_%s_%ss_%s-%s.fts" % (sufx, evinf[i].instr_short_name,
        str(timebin), str(mine), str(maxe)))
        lcbkglst.append("lcbkg%s_%s_%ss_%s-%s.fts" % (sufx, evinf[i].instr_short_name,
        str(timebin), str(mine), str(maxe)))
        lcnetlst.append("lcnet%s_%s_%ss_%s-%s.fts" % (sufx, evinf[i].instr_short_name,
        str(timebin), str(mine), str(maxe)))
        evtlclst.append("evtobj%s_%s_%ss_%s-%s.fts" % (sufx, evinf[i].instr_short_name,
        str(timebin), str(mine), str(maxe)))
        imglclst.append("imglc%s_%s_%ss_%s-%s.ps" % (sufx, evinf[i].instr_short_name,
        str(timebin), str(mine), str(maxe)))
        
        if os.path.isfile(lcobjlst[i]):
             die("File %s already exist. Please change the name suffix "
             "or clear the folder." % lcobjlst[i])
        if os.path.isfile(lcbkglst[i]):
             die("File %s already exist. Please change the name suffix "
             "or clear the folder." % lcbkglst[i])
        if os.path.isfile(lcnetlst[i]):
             die("File %s already exist. Please change the name suffix "
             "or clear the folder." % lcnetlst[i])
        
        if evinf[i].datamode!='IMAGING':
            die("'%s' is in '%s' mode. Unfortunately only 'IMAGING' " 
            "datamode supported." % (evtfile[i], evinf.datamode), "xmmlc")
        if not os.path.isfile(regobjfiles[i]):
            die("'%s' is not found" % regobjfiles[i], "xmmlc")
        if not os.path.isfile(regbkgfiles[i]):
            die("'%s' is not found" % regbkgfiles[i], "xmmlc")
        ftsevt.close()
    
    #Check GTI file        
    ftsgti=fitsopen(gtifile)
    if not chkisgti(ftsgti[1]):
        die("'%s' is not a valid GTI extension" % gtifile+"[1]",
        "xmmlc")
    ftsgti.close()
    
    #Check dublicates in instrlst
    lcres=lcnetlst[0]
    if nfiles>1:
        if nfiles!= len(dict([(item, None) for item in instrlst]).keys()):
            die("One or more light curves from the same detector. "
            "Cannot combine them.")
        lcres="lcnet%s_%s_%ss_%s-%s.fts" % (sufx, "".join(instrlst),
        str(timebin), str(mine), str(maxe))
    
    ###################################################    
    ##Perform processing
    printcaption("Making '%s'..." % lcres)
    for i in range(nfiles):
        if not makelightcurves(evtfiles[i], gtifile, regobjfiles[i],
        regbkgfiles[i], timebin, mine, maxe, lcobjlst[i], lcbkglst[i], 
        lcnetlst[i],evtlclst[i],imglclst[i]):
            die("Some errors arise during EVENT-file '%s' "
            "be processed." % evtfiles[i])
            

    if nfiles==2:
        if not combinelc(lcnetlst[0],lcnetlst[1],lcres):
            die("Cannot combine light curves '%s' and '%s'." % 
            (lcnetlst[0], lcnetlst[1]))
    elif nfiles==3:
        tmplc=".lcnet%s_tmp.fts" % sufx
        if not combinelc(lcnetlst[1],lcnetlst[2],tmplc):
            die("Cannot combine light curves '%s' and '%s'." % 
            (lcnetlst[1], lcnetlst[2]))
        if not combinelc(lcnetlst[0],tmplc,lcres):
            die("Cannot combine light curves '%s' and '%s'." % 
            (lcnetlst[0], tmplc)) 
        os.remove(tmplc)
    printbold("Saved "+lcres)
    printcaption("Finished")

            
