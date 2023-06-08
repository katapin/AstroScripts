#!/usr/bin/python
"""Perform filtering of XMM-Newton EVT-files."""


import sys
import os
import argparse
import gtiplot
import myfits as my
import xmmgeneral as xmm

def xmm_perform_filtering(evtinfo, evtflt:str, nroot:str, mode:xmm.FilteringMode):
    """Perform filtering of a single EVT-file.
    
    Parameters
    ----------
    evtinfo : xmm.EVTinfo 
        Basic information about EVT-file being processed (incuding the type
        of the instrument). Instance of EVTinfo class
    evtflt : str
        Name of EVT-file to save the result
    nroot : FilePath
        Name root for auxiliary files (GTI, Flash light curve, etc.)
    mode : FilteringMode Enum
        Chooses the expression for filtering

    Raises
    ------
    NotImplementedError, TaskError, ExternalTaskError
    
    Returns
    -------
    bool
        True if no errors arose, otherwise returns False

    """
    _no_errors=True
    evtname=evtinfo.filepath.basename
    _ownname=my.getownname()
    if evtinfo.datamode!='IMAGING':
        my.printerr("'{}' was taken in the '{}' mode. Only 'IMAGING' is " 
        "datamode supported yet.".format(evtname, evtinfo.datamode))
        raise NotImplementedError(f'{evtinfo.datamode} mode is not supported yet.')
        
    my.printgreen("1) Basic filtering")
    cmd="evselect table='{}' withfilteredset=Y filteredset='{}' "\
        "destruct=Y keepfilteroutput=T expression='{}'".format(evtname, evtflt, 
        xmm.xmm_get_expression(evtinfo.instrkey, xmm.FilteringPurpose.BASIC, mode))
    xmm.__call_and_check_result(cmd, evtflt, 'basic filtering', 'evselect', _ownname)
    
    #Making flash light curve
    my.printgreen("2) Making flash light curve")
    flshlc=nroot.starts_with('lcflsh_').ends_with('.fts')
    cmd="evselect table='{}' withrateset=Y rateset='{}' maketimecolumn=Y "\
        "timebinsize=100 makeratecolumn=Y expression='{}'".format(evtname, flshlc, 
        xmm.xmm_get_expression(evtinfo.instrkey, xmm.FilteringPurpose.FLASH, mode))
    xmm.__call_and_check_result(cmd, flshlc, 'flash filtering', 'evselect', _ownname)
    
    #Making flash GTI-file
    my.printgreen("3) Making flash GTI-file")
    if evtinfo.submode not in ['PrimeFullWindow', 'PrimeFullWindowExtended']:
        my.printwarn("Warning: Observation mode is not a FullFrame. "
        "The automatic count rate criterion may underestimate the significance "
        "of some faint flashes.")
    gtiflsh=nroot.starts_with('gtiflsh_').ends_with('.fts')
    cmd="tabgtigen table='{}' expression='{}' gtiset='{}'".format(flshlc,
        xmm.xmm_get_expression(evtinfo.instrkey, xmm.FilteringPurpose.GTI,
        mode), gtiflsh)
    xmm.__call_and_check_result(cmd, gtiflsh, 'GTI-file creation', 'tabgtigen', _ownname)
        
    
    #Making FITS-image
    my.printgreen("4) Making images")
    imgevt=nroot.starts_with('img_').ends_with('.fts')
    imgevtpng=imgevt.replace_extension('png')
    if not xmm.__make_test_images(evtflt, imgevtpng, outimg_fts=imgevt, progname=_ownname):
        _no_errors=False
        
    if evtinfo.instrkey=='EPN':
        imgevtdetxy=imgevt.ends_with('_detxy')
        imgevtdetxypng=imgevtdetxy.replace_extension('png')
        if not xmm.__make_test_images(evtflt, imgevtdetxypng, colX='DETX', colY='DETY', progname=_ownname):
            _no_errors=False
     
    imgflshpng=nroot.starts_with('imgflsh_').ends_with('.png')
    try:
        gtiplot.gtiplot([(gtiflsh, 1)], imgflshpng, flshlc, flshlc.basename)
        my.printbold("Saved "+imgflshpng)
    except my.TaskError:
        my.printwarn("Cannot save flash light curve as png image.")
        _no_errors=False
        
    return _no_errors


def _main():
    
    #Parse the arguments
    parser = argparse.ArgumentParser(description="Perform basic "
    "filtering of XMM-Newton EVENT-files")
    parser.add_argument('evtfile', nargs=1, help="name of the "
    "raw XMM-Newton EVENT-file")
    parser.add_argument('-u', '--suffix', nargs='?', help="name suffinx "
    "for output files", default='')
    parser.add_argument('-l', '--logfile', nargs='?', help="log file")
    parser.add_argument('-m', '--mode', nargs='?', type=str, default='standard',
        choices=[str(x.name) for x in xmm.FilteringMode], help="Filtering mode")
    
    argnspace=parser.parse_args(sys.argv[1:])
    sufx=argnspace.suffix
    logfile=argnspace.logfile
    mode=xmm.FilteringMode[argnspace.mode]  #Convert to enum

    #Check the system variables
    if ("SAS_ODF" not in os.environ) or ("SAS_ODF" not in os.environ):
        my.die("'SAS_ODF' and 'SAS_CCF' variables are not defined. Please \
    define the variables and try again.")
    
    evtpath=xmm.xmm_check_file_is_evt(argnspace.evtfile[0])
    evtinfo=xmm.EVTinfo(evtpath)
    
    #Create name root
    nroot = my.FilePath('-'.join([x for x in [evtinfo.instr_short_name, sufx] if x]))
    
    evtflt=nroot.ends_with('_flt.evt')     #Name of the fileted file
    my.check_file_not_exist_or_remove(evtflt,overwrite=False,action=my.Action.DIE,
             extra_text='Please use name suffixes (see --help) or clean up the folder.')
    
    if logfile:
        my.logging_turn_on(logfile)
    
    my.printcaption('{} processing...'.format(evtinfo.instr_full_name))
    try:
       if not xmm_perform_filtering(evtinfo, evtflt, nroot, mode):
           my.printwarn("Some minor errors arose. Check the result carefully.")
    except Exception as ex:
        my.die(f"{ex}. Cannot process '{evtpath.basename}'")
    my.printcaption("Finished")
    
    my.printandlog("Data mode:       {}".format(evtinfo.datamode))
    my.printandlog("Data submode:    {}".format(evtinfo.submode))
    my.printandlog("Filter:          {}".format(evtinfo.filter))
    my.printandlog("Exposure ID:     {}".format(evtinfo.expidstr))
    my.printandlog("Exposure length: {:.2f}".format(evtinfo.telapse))
    my.printandlog("Used chips: {}".format(', '.join(str(x) for x in evtinfo.chips)))
    my.printandlog("Frame time (sec.): {}".format(', '.join(str(x) for x in set(evtinfo.frmtime.values()))))
    
if __name__ == '__main__':    
    _main()

