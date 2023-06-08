#!/usr/bin/python
"""Extract light curves from XMM-Newton observation."""



import sys, os, argparse
import myfits as my
import xmmgeneral as xmm
from astropy.io import fits
from dataclasses import dataclass, asdict

@dataclass
class ProdNames():
    """Store names for product files."""
    
    raw:str     #Raw object light curve
    bkg:str     #Background spectrum  
    net:str     #Background subtracted object light curve
    evt:str  = None   #EVT-file with only the events that used in light curve
    img:str  = None   #Quick-look image
    
def generate_names(nroot:str, timebin:float, mine:float, maxe:float):
    """Generate standard names for products.

    Parameters
    ----------
    nroot : str
        Name root.
    timebin : float
        Time resolution of the light curve.
    mine : float
        Lower limit of the energy range.
    maxe : float
        Upper limit of the energy range.

    Returns
    -------
    TYPE
        Instance of the ProdNames dataclass

    """
    tail=f'{timebin}s_{mine}-{maxe}'
    return ProdNames(
        raw=f"lcobj{nroot}_{tail}.fts",
        bkg=f"lcbkg{nroot}_{tail}.fts",
        net=f"lcnet{nroot}_{tail}.fts",
        evt=f"evtlcobj{nroot}_{tail}.fts",
        img=f"imglc{nroot}_{tail}.eps")

def xmmlc_extract_single(evtinfo:xmm.EVTinfo, gtifile:str, regfile:str, 
        lcfile:str, timebin:float, channels:tuple, evtlc=None):
    """Extract single light from the EVT-file.

    Parameters
    ----------
    evtinfo : xmm.EVTinfo
        Basic information about EVT-file being processed (incuding the type
        of the instrument). Instance of EVTinfo class.
    gtifile : str
        Path to the GTI file.
    regfile : str
        Path to the region file.
    lcfile : str
        Name of the light curve file to save the result.
    timebin : float
        Desired temporal resolution of the light curve.
    channels : tuple
        Energy channels. 
    evtlc : TYPE, optional
        Path to the filtered EVT-file with events used to create the raw light 
        curve. The default is None.

    Returns
    -------
        Returns True if no errors arose.

    """
    _no_errors=True
    _ownname=my.getownname()
    evtname=evtinfo.filepath.basename
    
    if evtinfo.datamode!='IMAGING':
        my.printerr("'{}' was taken in the '{}' mode. Only 'IMAGING' is " 
        "datamode supported yet.".format(evtname, evtinfo.datamode))
        raise NotImplementedError(f'{evtinfo.datamode} mode is not supported yet.')
        
    region_expression = xmm.xmm_read_regfile(regfile)
    if not region_expression:
        my.printerr(f"Empty region file '{regfile}' or unsupported format. Only "
        "ds9 region format is supported.")
        raise my.TaskError(_ownname, lcfile)
        
    gtilimits=my.gti_get_limits(gtifile) 
    expression='(({}) && gti({},TIME) && PI in [{:d}:{:d}])'.format(region_expression,
            gtifile, *channels) 
    
    #It's important to use timemin and timemax evselect's arguments to make
    #the combined light curves synchronous
    cmd="evselect table='{}' energycolumn=PI expression='{}' withrateset=yes "\
        "rateset='{}' timebinsize={:f} maketimecolumn=yes makeratecolumn=yes "\
        "timemin={:f} timemax={:f}".format(evtinfo.filepath, expression, lcfile,
        timebin, *gtilimits) 
    
    if evtlc:
        cmd += f" keepfilteroutput=yes withfilteredset=yes filteredset='{evtlc}'"
    
    xmm.__call_and_check_result(cmd, lcfile, 'exctraction of a light curve',
                'evselect', _ownname)
    
    #Make test images
    testimgpng=my.FilePath(lcfile).starts_with('tmpimg_').replace_extension('png')
    
    if not xmm.__make_test_images(evtinfo.filepath, testimgpng, expression, progname=_ownname):
        _no_errors=False

    return _no_errors
    
def xmmlc_make_products(evtinfo:xmm.EVTinfo, gtifile:str, objreg:str, bkgreg:str,
        prod_names:ProdNames, timebin=None, channels:tuple=None, 
        with_evtlc:bool=False, plotlc:bool=True):
    """Extract a triple of light curves.
    
    Extracts both the source and backgrounds light curves and also produces
    the net (background subtracted) light curve corrected for exposure map
    and the inctrument effective area.

    Parameters
    ----------
    evtinfo : xmm.EVTinfo
        Basic information about EVT-file being processed (incuding the type
        of the instrument). Instance of EVTinfo class..
    gtifile : str
        Path to the GTI file.
    objreg : str
        Path to the region file to extract light curve of the studied source.
    bkgreg : str
        Path to the region file to estimate the background level.
    prod_names : ProdNames
        Names of files to save specta and auxiliary files.
    timebin : TYPE, optional
        Desired temporal resolution of the light curve. If not provided then 
        it will be set equal to the frame time from the EVT-file.
    channels : tuple, optional
        Energy channels. The default is (200, 12000).
    with_evtlc : bool, optional
        Leave evtlc as well. The default is False.
    plotlc : bool, optional
        Produce quick-look image. The default is True.

    Returns
    -------
    _no_errors : TYPE
        Returns True if no errors arose.

    """
    _no_errors=True
    _ownname=my.getownname()
    evtname=evtinfo.filepath.basename
    #Determine timebin and channels
    if not timebin:
        timebin = max(evtinfo.frmtime)
        my.printwarn("'timebin' for {} is set to {}".format(evtinfo.instr_short_name, timebin))
        
    if not channels:
        channels = (200, 12000)
        my.printwarn("'channels' for {} is set to ({}, {})".format(evtinfo.instr_short_name, *channels))
        
    
    if with_evtlc:
        if not prod_names.evt:
            raise KeyError("prod_names.evt is not found. The option 'with_evtlc=True' "
               "requires the evt-file to save the result.")
        evtlc=prod_names.evt
    else:
        evtlc=None
        
        
    my.printgreen(f"Extracting the raw light curve '{evtname}' -> '{prod_names.raw}'")
    if not xmmlc_extract_single(evtinfo, gtifile, objreg, prod_names.raw, 
                                timebin, channels, evtlc):
        my.printwarn("Some minor errors arose during extraction of the raw light curve. "
             "Check the result carefully.")
        _no_errors=False
        
    my.printgreen(f"Extracting the background light curve '{evtname}' -> '{prod_names.bkg}'")
    if not xmmlc_extract_single(evtinfo, gtifile, bkgreg, prod_names.bkg, 
                                timebin, channels, None):
        my.printwarn("Some minor errors arose during extraction of the background light curve. "
             "Check the result carefully.")
        _no_errors=False
    
    my.printgreen("Producing the net light curve ('{0.raw}' - '{0.bkg}') -> '{0.net}'".format(prod_names))
    cmd="epiclccorr srctslist='{0.raw}' eventlist='{1}' outset='{0.net}' "\
    "withbkgset=yes bkgtslist='{0.bkg}' applyabsolutecorrections=yes".format(prod_names, evtinfo.filepath)
    xmm.__call_and_check_result(cmd, prod_names.net, 'background subtraction', 'evselect', _ownname)
    
    #Plot light curves    
    if plotlc:
        if not prod_names.img:
            raise KeyError("prod_names.img is not found. The option 'plotlc=True' "
               "requires the name of the eps-image file to store the result.")
        imglc=prod_names.img
        lcnet_header = fits.getheader(prod_names.net, 1)
        bkgratio=lcnet_header['BKGRATIO']
        if bkgratio>1: bkgratio=1/bkgratio
        try:
            my.fitslc_plot(prod_names.raw, gtifile, prod_names.img, prod_names.bkg, 
                           bkgratio, prod_names.net)
            my.printbold(f"Image saved in '{imglc}'")
        except:
            my.printerr("Can't plot the light curve.")
            _no_errors=False
            
    return _no_errors
                
def xmmlc_combine(lc1:str, lc2:str, reslc:str):
    """Combine two synchronous light curves.
    
    The lcmath task is used. So the resultant light curve will be a bin-to-bin
    sum of the two summands. Arguments are paths of the input ant output files.
    """
    _ownname=my.getownname()
    cmd=f"lcmath '{lc1}' '{lc2}' {reslc} 1 1 addsubr=yes"
    xmm.__call_and_check_result(cmd, reslc, 'merging of the light curves', 'lcmath', _ownname)
    
    
def _checks_for_xmmlc_make_products(prod_names:ProdNames, clobber):
    for ftype, fname in asdict(prod_names).items():
        if ftype in ['raw', 'net']:  #These are important files, warn!
            my.check_file_not_exist_or_remove(fname, overwrite=clobber, action=my.Action.DIE, 
                extra_text='Use option --clobber to override it or use --suffix.', 
                remove_warning=True)
        else:  #These are not important files, we can remove it
            my.check_file_not_exist_or_remove(fname, overwrite=True, 
                action=my.Action.WARNING, remove_warning=(not clobber))


def _main():
    
    #Parse the arguments
    parser = argparse.ArgumentParser(description="Exctact the light curve "
    "with certain time resolution and energy range. Programm can combine "
    "light curves from different detectors. Arguments 'evtfile', 'regobj' "
    "and 'regbkg' may be lists of files (pn,mos1,mos2) "
    "but a sigle GTI file must be common for all the detectors.")
    parser.add_argument('evtfile', nargs=1, help="list of the "
    "XMM-Newton EVENT-files")
    parser.add_argument('gtifile', nargs=1, help="name of the GTI-file")
    parser.add_argument('regobj', nargs=1, help="name of the object region file")
    parser.add_argument('regbkg', nargs=1, help="name of the background region file")
    parser.add_argument('timebin', nargs='?', type=float, help="time resolutin of "
        "the light curve in seconds")
    parser.add_argument('mine', nargs='?', type=float, help="lower limit of the energy range (keV)", 
                        default=0.2)
    parser.add_argument('maxe', nargs='?', type=float, help="higher limit of the energy range (keV)",
                        default=12.0)
    parser.add_argument('-u', '--suffix', nargs='?', help="name suffinx for output files")
    parser.add_argument('-l', '--logfile', nargs='?', help="log file")
    parser.add_argument('--clobber', action='store_true', help="Allow to override files")
    parser.add_argument('--with-evtlc', action='store_true', help="Leave filtered EVT-file with "
                        "events used to create the raw light curve")
    
    argnspace=parser.parse_args(sys.argv[1:])
    evtfiles=argnspace.evtfile[0].split(',')
    gtifile=argnspace.gtifile[0]
    regobjfiles=argnspace.regobj[0].split(',')
    regbkgfiles=argnspace.regbkg[0].split(',')
    timebin=argnspace.timebin
    mine=argnspace.mine
    maxe=argnspace.maxe
    logfile=argnspace.logfile
    
    if argnspace.suffix:
        sufx="-"+argnspace.suffix
    else:
        sufx=""

    
    #############################################################
    #Check input data
        
    #Check the system variables
    if ("SAS_ODF" not in os.environ) or ("SAS_ODF" not in os.environ):
        my.die("'SAS_ODF' and 'SAS_CCF' variables are not defined. Please \
    define the variables and try again.")
    
    #Check length of lists
    if not ( len(evtfiles) == len(regobjfiles)== len(regbkgfiles) ):
        my.die("Lists of the EVENT and the region files must have the same "
        "length and order.")
    nfiles=len(evtfiles)
    
    if nfiles>3:
        my.die("You are trying to combine more then 3 light curves. Only "
        "synchronous light curves from EPIC-PN, EPIC-MOS1 and EPIC-MOS2 "
        "can be combined.")
        
    evtdict={}
    regobjdict={}
    regbkgdict={}
    for i,evtfile in enumerate(evtfiles):  #Check each EVT file
        evtinfo=xmm.EVTinfo(xmm.xmm_check_file_is_evt(evtfile))
        instr=evtinfo.instr_short_name
        if instr in evtdict:        #Check for dublicates
            my.die("You are trying to extract and combine light curves from the same "
                   "detector. Can't combine them.")
        evtdict[instr] = evtinfo   #And add to a dict
        regobjdict[instr]=regobjfiles[i]
        regbkgdict[instr]=regbkgfiles[i]
        
    #Check files exist
    my.fits_check_file_is_gti(gtifile)
    for regfile in regobjfiles + regbkgfiles:
        my.check_file_exists(regfile)    
    
    #Check energy range
    if mine<0.2 or mine>12:
        my.die("'mine' argument must be in the 0.2-12.0 keV range")
    if maxe<0.2 or maxe>12:
        my.die("'maxe' argument must be in the 0.2-12.0 keV range")
    if mine>maxe:
        my.die("'mine' must be less than 'maxe'")
    
    channels=(int(mine*1000), int(maxe*1000))
    
    if logfile:
        my.logging_turn_on(logfile)
    
    _no_error=True
    if nfiles>1:  #Lcurves wil be combined
        if not timebin:
            my.die("If you want to combine light curves, you must specify the value "
                   "of the 'timebin' argument.")
            
        prodnames_dict={}
        _kwargs_generate=dict(timebin=timebin, mine=mine, maxe=maxe)
        for instr in evtdict:  #Generate names for output files and check they are already NOT exist
            prodnames_dict[instr] = generate_names(f'{sufx}_{instr}', **_kwargs_generate)
            _checks_for_xmmlc_make_products(prodnames_dict[instr], argnspace.clobber)
            
        #Trick to generate name of the combined net LC
        comb_nroot='{}_{}'.format(sufx, ''.join(sorted(evtdict)))
        _kwargs_check = dict(overwrite=argnspace.clobber, action=my.Action.DIE, 
        extra_text='Use option --clobber to override it or use --suffix.', remove_warning=True)
        lc_combined = my.check_file_not_exist_or_remove(generate_names(comb_nroot,
              **_kwargs_generate).net, **_kwargs_check)
        if nfiles==3:
            lc_mos_combined = my.check_file_not_exist_or_remove(generate_names(f'{sufx}_mos1mos2',
                  **_kwargs_generate).net, **_kwargs_check)
            
        
        #### Do the main job
        my.printcaption(f"Making '{lc_combined.basename}'...")
        for instr in evtdict:
            if not xmmlc_make_products(evtdict[instr], gtifile, regobjdict[instr], 
                           regbkgdict[instr], prodnames_dict[instr], timebin, channels, 
                           with_evtlc=argnspace.with_evtlc):
                _no_error=False
        
        if nfiles==3:
            xmmlc_combine(prodnames_dict['mos1'].net,prodnames_dict['mos2'].net, lc_mos_combined)
            xmmlc_combine(lc_mos_combined,prodnames_dict['pn'].net, lc_combined)
        else:
            lc_pair=[x.net for x in prodnames_dict.values() ]
            xmmlc_combine(*lc_pair, lc_combined)
        my.printcaption("Finished")
    else: #Extract single lcurves without combining
        prodnames = generate_names(f'{sufx}_{instr}', evtinfo.frmtime, mine, maxe)
        if not xmmlc_make_products(evtinfo, gtifile, regobjfiles[0], regbkgfiles[0], 
               prodnames, evtinfo.frmtime, channels, with_evtlc=argparse.with_evtlc):
            _no_error=False
            
    if not _no_error:
        my.printwarn("Some minor errors arose. Check the result carefully.")
            
            
if __name__ == '__main__':        
    _main()