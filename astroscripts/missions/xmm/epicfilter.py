#!/usr/bin/python
"""Perform filtering of XMM-Newton EVT-files."""


import sys
import os
import argparse

import mypythonlib as mylib
from mypythonlib import FilePath, FilePathAbs
from astroscripts.gtiplot import gtiplot
from astroscripts import ExtPathAbs, TaskError

import astroscripts.missions.xmm.common as xmm
from astroscripts.missions.xmm.common import EVTinfo


def xmm_perform_filtering(evtinfo: EVTinfo, evtflt: FilePathAbs,
                          nroot: FilePath, mode: xmm.FilteringMode) -> bool:
    """Perform filtering of a single EVT-file.
    
    :param evtinfo: Basic information about EVT-file being processed
        (including the type of the instrument). Instance of EVTinfo class.
    :param evtflt: Name of EVT-file to save the result
    :param nroot: Name root for auxiliary files (GTI, Flash light curve, etc.)
    :param mode: Chooses the expression for filtering
    :returns: True if no errors arose, otherwise returns False
    :raises TaskError: if an error arises in an internal task
    :raises ExternalTaskError: if an error arises in an external task
    """
    _no_errors=True
    evtpath=evtinfo.filepath.__fspath__()
    _ownname=mylib.getownname()
    if evtinfo.datamode!='IMAGING':
        mylib.printerr("'{}' was taken in the '{}' mode. Only 'IMAGING' is " 
                    "datamode supported yet.".format(evtpath, evtinfo.datamode))
        raise NotImplementedError(f'{evtinfo.datamode} mode is not supported yet.')
        
    mylib.printgreen("1) Basic filtering")
    cmd="evselect table='{}' withfilteredset=Y filteredset='{}' "\
        "destruct=Y keepfilteroutput=T expression='{}'".format(evtpath, evtflt,
        xmm.xmm_get_expression(evtinfo.instrkey, xmm.FilteringPurpose.BASIC, mode))
    xmm._call_and_check_result(cmd, evtflt, 'basic filtering', 'evselect', _ownname)
    
    # Making flash light curve
    mylib.printgreen("2) Making flash light curve")
    flshlc=nroot.with_stem_starting('lcflsh_').with_suffix('.fts').absolute()
    cmd="evselect table='{}' withrateset=Y rateset='{}' maketimecolumn=Y "\
        "timebinsize=100 makeratecolumn=Y expression='{}'".format(evtpath, flshlc,
        xmm.xmm_get_expression(evtinfo.instrkey, xmm.FilteringPurpose.FLASH, mode))
    xmm._call_and_check_result(cmd, flshlc, 'flash filtering', 'evselect', _ownname)
    
    # Making flash GTI-file
    mylib.printgreen("3) Making flash GTI-file")
    if evtinfo.submode not in ['PrimeFullWindow', 'PrimeFullWindowExtended']:
        mylib.printwarn("Warning: Observation mode is not a FullFrame. "
                     "The automatic count rate criterion may underestimate the significance "
                     "of some faint flashes.")
    gtiflsh=nroot.with_stem_starting('gtiflsh_').with_suffix('.fts').absolute()
    cmd="tabgtigen table='{}' expression='{}' gtiset='{}'".format(flshlc,
        xmm.xmm_get_expression(evtinfo.instrkey, xmm.FilteringPurpose.GTI, mode), gtiflsh)
    xmm._call_and_check_result(cmd, gtiflsh, 'GTI-file creation', 'tabgtigen', _ownname)

    # Making FITS-image
    mylib.printgreen("4) Making images")
    imgevt=nroot.with_stem_starting('img_').with_suffix('.fts').absolute()
    imgevtpng=imgevt.with_suffix('.png').absolute()
    if not xmm._make_test_images(evtflt, imgevtpng, outimg_fts=imgevt, progname=_ownname):
        _no_errors=False
        
    if evtinfo.instrkey=='EPN':
        imgevtdetxy=imgevt.ends_with('_detxy')
        imgevtdetxypng=imgevtdetxy.replace_extension('png')
        if not xmm._make_test_images(evtflt, imgevtdetxypng, colX='DETX', colY='DETY', progname=_ownname):
            _no_errors=False
     
    imgflshpng=nroot.with_stem_starting('imgflsh_').with_stem_ending('.png').absolute()
    try:
        gtiplot(ExtPathAbs(gtiflsh, hdu=1), lcpath=flshlc, pngpath=imgflshpng, onlysave=True)
    except TaskError:
        mylib.printwarn("Cannot save flash light curve as png image.")
        _no_errors=False
        
    return _no_errors


def _main():
    
    # Parse the arguments
    parser = argparse.ArgumentParser(
        description="Perform basic filtering of XMM-Newton EVENT-files")
    parser.add_argument(
        'evtfile', nargs=1, help="name of the raw XMM-Newton EVENT-file")
    parser.add_argument(
        '-u', '--suffix', nargs='?', help="name suffinx for output files", default='')
    parser.add_argument('-l', '--logfile', nargs='?', help="log file")
    parser.add_argument('-m', '--mode', nargs='?', type=str, default='standard',
                        choices=[x.name for x in xmm.FilteringMode], help="Filtering mode")
    argnspace = parser.parse_args(sys.argv[1:])
    sufx = argnspace.suffix
    logfile = argnspace.logfile
    mode = xmm.FilteringMode[argnspace.mode]  # Convert to enum

    # Check the system variables
    if ("SAS_ODF" not in os.environ) or ("SAS_ODF" not in os.environ):
        mylib.die("'SAS_ODF' and 'SAS_CCF' variables are not defined. Please \
    define the variables and try again.")
    
    evtpath = xmm.xmm_check_file_is_evt(argnspace.evtfile[0])
    evtinfo = xmm.EVTinfo(evtpath)

    # Create name root
    nroot = mylib.FilePath('-'.join([x for x in [evtinfo.instr_short_name, sufx] if x]))

    # Name of the filted file
    evtflt = mylib.check_file_not_exist_or_remove(
        nroot.with_stem_ending('_flt.evt'), override=False, action=mylib.Actions.DIE,
        extra_text='Please use name suffixes (see --help) or clean up the folder.')

    if logfile:
        mylib.logger_turn_on(mylib.FilePath(logfile).absolute())
    
    mylib.printcaption('{} processing...'.format(evtinfo.instr_full_name))
    try:
        if not xmm_perform_filtering(evtinfo, evtflt, nroot, mode):
            mylib.printwarn("Some minor errors arose. Check the result carefully.")
    except Exception as ex:
        mylib.die(f"{ex}. Cannot process '{evtpath.name}'")
    mylib.printcaption("Finished")
    mylib.printbold(f"Info for the original EVT-file ({evtinfo.filepath.name}):", progname='')
    evtinfo.describe()

