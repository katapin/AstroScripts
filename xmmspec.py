#!/usr/bin/python
"""Extract spectra from XMM-Newton observations."""

import sys, os, argparse
import myfits as my
import xmmgeneral as xmm
from dataclasses import dataclass, asdict

@dataclass
class ProdNames():
    """Store names for product files."""
    
    src:str     #Raw object spectrum
    bkg:str     #Background spectrum  
    grp:str     #Grouped spectum
    rmf:str     #Main responce
    arf:str     #Auxiliary responce
    img:str = None   #Quick-look image
    

def generate_names(nroot, binmin):
    """Generate standard names for products.

    Parameters
    ----------
    nroot : str
        Name root
    binmin : int
        Minimal number of channels to bin

    Returns
    -------
    dict
        Instance of the ProdNames dataclass
    """
    return ProdNames(
        src=f"spobj{nroot}.fts",
        bkg=f"spbkg{nroot}.fts",
        grp=f"spgrp{nroot}_{binmin:d}.fts",
        rmf=f"rmf{nroot}.fts",
        arf=f"arf{nroot}.fts",
        img=f"imgsp{nroot}_{binmin:d}.eps")
    

def xmmspec_extract_single(evtinfo:xmm.EVTinfo, gtifile:str, regfile:str, 
        specfile:str, mode:xmm.FilteringMode):
    """Extract single spectrum from the EVT-file.
    
    Extracts a spectrum using the given GTI and region files. Also it calculates
    BACKSCAL keyword and produses test images.
    

    Parameters
    ----------
    evtinfo : xmm.EVTinfo
        Basic information about EVT-file being processed (incuding the type
        of the instrument). Instance of EVTinfo class.
    gtifile : str
        Path to the GTI file.
    regfile : str
        Path to the region file.
    specfile : str
        Name of the spectrum file to save the result.
    mode : xmm.FilteringMode
        Filtering mode (standard of strict). It determines the expression base
        that will be used to extract the spectrum.

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
        raise my.TaskError(_ownname, specfile)
    
    base_expression = xmm.xmm_get_expression(evtinfo.instrkey, xmm.FilteringPurpose.SPEC, mode)
    expression='(({}) && gti({},TIME) && ({}))'.format(region_expression, gtifile,
           base_expression)
    
    instr_type = xmm.xmm_get_instrument_type(evtinfo.instrkey)
    #Get channels and binsize of the specific instument
    if instr_type  == 'PN':
        chan_min, chan_max = xmm.SPEC_CHANNELS_PN
        binsize = xmm.SPEC_BINSIZE_PN
    elif instr_type  == 'MOS':
        chan_min, chan_max = xmm.SPEC_CHANNELS_MOS
        binsize = xmm.SPEC_BINSIZE_MOS
    else:
        raise NotImplementedError(f'Unknown instrument type {instr_type}')
        
    #Extract spectrum
    cmd=f"evselect table='{evtinfo.filepath}' energycolumn=PI expression='{expression}' "\
    f"withspectrumset=yes spectrumset='{specfile}' spectralbinsize={binsize:d} "\
    f"withspecranges=yes specchannelmin={chan_min:d} specchannelmax={chan_max:d}"
    xmm.__call_and_check_result(cmd, specfile, 'exctraction of spectrum', 'evselect',  _ownname)

    #Calculate backscale    
    cmd=f"backscale spectrumset='{specfile}' badpixlocation='{evtinfo.filepath}'"
    xmm.__call_and_check_result(cmd, specfile, 'backscale correction', 'backscale',  _ownname)
        
    #Make test images
    testimgpng=my.FilePath(specfile).starts_with('tmpimg_').replace_extension('png')
    if not xmm.__make_test_images(evtinfo.filepath, testimgpng, expression, progname=_ownname):
        _no_errors=False

    return _no_errors

def xmmspec_make_products(evtinfo:xmm.EVTinfo, gtifile:str, objreg:str, bkgreg:str, 
          prod_names:ProdNames, mode:xmm.FilteringMode):
    """Extract spectum together with all the correspoding auxiliary files.

    Parameters
    ----------
    evtinfo : xmm.EVTinfo
       Basic information about EVT-file being processed (incuding the type
       of the instrument). Instance of EVTinfo class 
    gtifile : str
        Path to the GTI file.
    objreg : str
        Path to the region file to extract spectrum of the studied source. 
    bkgreg : str
        Path to the region file to extract spectrum of the background
    prod_names : ProdNames
        Names of files to save specta and auxiliary files
    mode : FilteringMode Enum
        Filtering mode (standard of strict). It determines the expression base
        that will be used to extract the spectra.
    logfile : str, optional
        Name of the logfile

    Returns
    -------
    bool
        True if no errors arose, otherwise returns False.
    """
    _no_errors=True
    _ownname=my.getownname()
    evtname=evtinfo.filepath.basename
    my.printgreen(f"Extracting the object spectrum '{evtname}' -> '{prod_names.src}'")
    if not xmmspec_extract_single(evtinfo, gtifile, objreg, prod_names.src, mode):
        my.printwarn("Some minor errors arose during extraction of the object spectrum. "
             "Check the result carefully.")
        _no_errors=False
        
    my.printgreen(f"Extracting the background spectrum '{evtname}' -> '{prod_names.bkg}")
    if not xmmspec_extract_single(evtinfo, gtifile, bkgreg, prod_names.bkg, mode):
        my.printwarn("Some minor errors arose during extraction of the background spectrum. "
             "Check the result carefully.")
        _no_errors=False
    
    #Make RMF
    my.printgreen("Generating RMF '{0.rmf}' for '{0.src}'".format(prod_names))
    cmd="rmfgen spectrumset='{0.src}' rmfset='{0.rmf}'".format(prod_names)
    xmm.__call_and_check_result(cmd, prod_names.rmf, 'RMF production', 
            'rmfgen ',  _ownname)
        
    #Make ARF    
    my.printgreen("Generating ARF '{0.arf}' for '{0.src}'".format(prod_names))
    cmd="arfgen spectrumset='{0.src}' arfset='{0.arf}' withrmfset=yes " \
    "rmfset='{0.rmf}' badpixlocation='{evtfile}' detmaptype=psf".format(prod_names,
          evtfile=evtinfo.filepath)
    xmm.__call_and_check_result(cmd, prod_names.arf, 'ARF production', 
            'arfgen ',  _ownname)
            
    return _no_errors


def xmmspec_grouping(prod_names:ProdNames, binmin:int):
    """Perform grouping of the spcetrum.

    Parameters
    ----------
    prod_names : ProdNames
        Struct stores the filenames of the spectrum and its auxiliary filess.
    binmin : int
        Minimim counts to group.

    Returns
    -------
    None.

    """
    _ownname=my.getownname()
    my.printgreen("Perform grouping of '{0.src}: min {1:d} counts".format(prod_names, binmin))
    cmd="specgroup spectrumset='{0.src}' backgndset='{0.bkg}' rmfset='{0.rmf}' " \
    "arfset='{0.arf}' groupedset='{0.grp}' mincounts={1:d}".format(prod_names, binmin)
    xmm.__call_and_check_result(cmd, prod_names.grp, 'grouping', 'arfgen ', progname=_ownname)
    
def _checks_for_xmmspec_make_products(prod_names:ProdNames, clobber):
    for ftype, fname in asdict(prod_names).items():
        if ftype in ['src', 'grp']:  #These are important files, warn!
            my.check_file_not_exist_or_remove(fname, overwrite=clobber, action=my.Action.DIE, 
                extra_text='Use option --clobber to override it or use --suffix.', 
                remove_warning=True)
        else:  #These are not important files, we can remove it
            my.check_file_not_exist_or_remove(fname, overwrite=True, 
                action=my.Action.WARNING, remove_warning=(not clobber))
            
def _checks_for_regrouping(prod_names:ProdNames, clobber):
    for ftype, fname in asdict(prod_names).items():
        if ftype == 'grp':
            my.check_file_not_exist_or_remove(fname, overwrite=clobber, action=my.Action.DIE, 
                extra_text='Use option --clobber to override it or use --suffix.', 
                remove_warning=True)
        elif ftype == 'img':
            my.check_file_not_exist_or_remove(fname, overwrite=True, remove_warning=False)
        else: 
            my.check_file_exists(fname, action=my.Action.DIE)

def _main():
    #Parse the arguments
    parser = argparse.ArgumentParser(description="Extract the spectrum, " \
    "make corresponding RMF and ARF.")
    parser.add_argument('evtfile', nargs=1, help="list of the "
    "XMM-Newton EVENT-files")
    parser.add_argument('gtifile', nargs=1, help="name of the GTI-file")
    parser.add_argument('regobj', nargs=1, help="name of the object"
    "region file")
    parser.add_argument('regbkg', nargs=1, help="name of the background"
    "region file")
    parser.add_argument('binmin', nargs=1, type=int, help="counts "
    "per grouped channel")
    parser.add_argument('-u', '--suffix', nargs='?', help="name suffinx "
    "for output files")
    parser.add_argument('-l', '--logfile', nargs='?', help="log file")
    parser.add_argument('-m', '--mode', nargs='?', type=str, default='standard',
        choices=[str(x.name) for x in xmm.FilteringMode], help="Filtering mode")
    parser.add_argument('--clobber', action='store_true', help="Allow to override files")
    parser.add_argument('--regroup', action='store_true', help="Re-group existing spectrum")
    
    argnspace=parser.parse_args(sys.argv[1:])
    binmin  = argnspace.binmin[0]
    sufx    = '-'+argnspace.suffix if argnspace.suffix else ''
    logfile = argnspace.logfile
    mode    = xmm.FilteringMode[argnspace.mode]
    
    
    #############################################################
    #Check input data
        
    #Check the system variables
    if ("SAS_ODF" not in os.environ) or ("SAS_ODF" not in os.environ):
        my.die("'SAS_ODF' and 'SAS_CCF' variables are not defined. Please \
    define the variables and try again.")
    
    
    evtpath = xmm.xmm_check_file_is_evt(argnspace.evtfile[0])
    evtinfo=xmm.EVTinfo(evtpath)
    nroot = '{}_{}'.format(sufx, evtinfo.instr_short_name)
    prod_names = generate_names(nroot, binmin)
    
    if logfile:
        my.logging_turn_on(logfile)
    
    if not argnspace.regroup:  #Standard abalysis
        _checks_for_xmmspec_make_products(prod_names, argnspace.clobber)
        gtipath = my.fits_check_file_is_gti(argnspace.gtifile[0])
        regobj  = my.check_file_exists(argnspace.regobj[0])
        regbkg  = my.check_file_exists(argnspace.regbkg[0])
    
        _no_error=True
        my.printcaption("Making '{0.grp}'...".format(prod_names))
        if not xmmspec_make_products(evtinfo, gtipath, regobj, regbkg, prod_names, mode):
            _no_error=False
            
        try:
            xmmspec_grouping(prod_names, binmin)
        except my.TaskError:
            my.printerr("Can't perform grouping for the produced spectrum")
            _no_error=False
        
        my.printcaption("Finished")
        if not _no_error:
            my.printwarn("Some minor errors arose. Check the result carefully.")
            
    else:
        _checks_for_regrouping(prod_names, argnspace.clobber)
        xmmspec_grouping(prod_names, binmin)
 
    #TODO:
    #2)Make ps-image
    #3) pile-up
    
    #if not specplot(lcobj,lcbkg,lcnet,gtifile,imglc):
        #printerr("Can't plot light curve")
            
if __name__ == '__main__': 
    _main()