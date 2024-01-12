#!/usr/bin/python
"""Extract light curves from observation taken with XMM-Newton EPIC detectors."""

import os
import sys
from typing import Self
from dataclasses import dataclass, asdict
from astropy.io import fits

import mypythonlib as mylib
from mypythonlib import FilePath, FilePathAbs
from astroscripts.external import xronos_plot_lcurve
from astroscripts import ExtPathAbs, TaskError, gti_get_limits, fits_check_file_is_gti
import astroscripts.missions.xmm.common as xmm
from astroscripts.missions.xmm.common import EVTinfo


@dataclass
class ProdNames:
    """Store names for product files."""
    
    raw: str           # Raw object light curve
    bkg: str           # Background spectrum
    net: str           # Background subtracted object light curve
    evt: str = None   # EVT-file with only the events that used in light curve
    img: str = None   # Quick-look image

    @classmethod
    def generate(cls, nroot: str, timebin: float, mine: float, maxe: float) -> Self:
        """Generate standard names for products.

        :param nroot: Name root.
        :param timebin: Time resolution of the light curve.
        :param mine: Lower limit of the energy range.
        :param maxe: Upper limit of the energy range.
        :returns: Instance of the ProdNames dataclass
        """
        tail=f'{timebin}s_{mine}-{maxe}'
        return cls(
            raw=f"lcobj{nroot}_{tail}.fts",
            bkg=f"lcbkg{nroot}_{tail}.fts",
            net=f"lcnet{nroot}_{tail}.fts",
            evt=f"evtlcobj{nroot}_{tail}.fts",
            img=f"imglc{nroot}_{tail}.eps"
        )

    def make_abspaths(self):
        """Convert filenames to absolute paths."""
        abs={}
        for ftype, fname in asdict(self).items():
            abs[ftype] = FilePath(fname).absolute()
        self.abs = abs


def xmmlc_extract_single(
        evtinfo: xmm.EVTinfo, gtifile: ExtPathAbs, regfile: FilePathAbs,
        lcfile: FilePathAbs, timebin: float, channels: tuple,
        evtlc: FilePathAbs = None) -> bool:
    """Extract single light from the EVT-file.

    :param evtinfo: Basic information about EVT-file being processed
        (incuding the type of the instrument). Instance of EVTinfo class.
    :param gtifile: Path to the GTI file.
    :param regfile: Path to the region file.
    :param lcfile: Name of the light curve file to save the result.
    :param timebin: Desired temporal resolution of the light curve.
    :param channels: Energy channels.
    :param evtlc: Path to the filtered EVT-file with events used
        to create the raw light curve. The default is None.
    :returns: Returns True if no errors arose.
    """
    _no_errors=True
    _ownname=mylib.getownname()
    evtname, evtpath = evtinfo.filepath.name, evtinfo.filepath.fspath
    
    if evtinfo.datamode!='IMAGING':
        mylib.printerr("'{}' was taken in the '{}' mode. Only 'IMAGING' is " 
                    "datamode supported yet.".format(evtname, evtinfo.datamode))
        raise NotImplementedError(f'{evtinfo.datamode} mode is not supported yet.')
        
    region_expression = xmm.xmm_read_regfile(regfile)
    if not region_expression:
        mylib.printerr(f"Empty region file '{regfile}' or unsupported format. Only "
                    "ds9 region format is supported.")
        raise TaskError(_ownname, lcfile)
        
    gtilimits=gti_get_limits(gtifile)
    expression='(({}) && gti({},TIME) && PI in [{:d}:{:d}])'.format(region_expression,
            gtifile, *channels) 
    
    # It's important to use timemin and timemax evselect's arguments to make
    # the combined light curves synchronous
    cmd="evselect table='{}' energycolumn=PI expression='{}' withrateset=yes "\
        "rateset='{}' timebinsize={:f} maketimecolumn=yes makeratecolumn=yes "\
        "timemin={:f} timemax={:f}".format(evtpath, expression, lcfile,
        timebin, *gtilimits) 
    
    if evtlc:
        cmd += f" keepfilteroutput=yes withfilteredset=yes filteredset='{evtlc}'"
    
    xmm._call_and_check_result(cmd, lcfile, 'extracting a light curve',
                                'evselect', _ownname)
    
    # Make test images
    testimgpng=lcfile.with_stem_starting('tmpimg_').with_suffix('.png')
    
    if not xmm._make_test_images(evtinfo.filepath, testimgpng, expression, progname=_ownname):
        _no_errors=False

    return _no_errors


def xmmlc_make_products(
        evtinfo: EVTinfo, gtifile: ExtPathAbs, objreg: FilePathAbs,
        bkgreg: FilePathAbs, prod_names: ProdNames, timebin: float = None,
        channels: tuple = None, with_evtlc: bool = False, plotlc: bool = True) -> bool:
    """Extract a triple of light curves.
    
    Extracts both the source and backgrounds light curves and also produces
    the net (background subtracted) light curve corrected for exposure map
    and the instrument effective area.

    :param evtinfo: Basic information about EVT-file being processed
        (including the type of the instrument). Instance of EVTinfo class.
    :param gtifile: Path to the GTI file.
    :param objreg: Path to the region file to extract light curve of
        the studied source.
    :param bkgreg: Path to the region file to estimate the background level.
    :param prod_names: Names of files to save spectra and auxiliary files.
    :param timebin: Desired temporal resolution of the light curve.
        If not provided then it will be set equal to the frame time
        from the EVT-file.
    :param channels: Energy channels. The default is (200, 12000).
    :param with_evtlc: Save events used to build light curve.
    :param plotlc: Produce quick-look image. The default is True.
    :returns: True if no errors arose.
    """
    _no_errors=True
    _ownname=mylib.getownname()
    evtname=evtinfo.filepath.name
    if not hasattr(prod_names, 'abs'): prod_names.make_abspaths()

    # Determine timebin and channels
    if not timebin:
        timebin = max(evtinfo.frmtime.values())
        mylib.printwarn("'timebin' for {} is set to {}".format(evtinfo.instr_short_name, timebin))
        
    if not channels:
        channels = (200, 12000)
        mylib.printwarn("'channels' for {} is set to ({}, {})".format(evtinfo.instr_short_name, *channels))
        
    if with_evtlc:
        if not prod_names.evt:
            raise KeyError("prod_names.evt is not found. The option 'with_evtlc=True' "
                           "requires the evt-file to save the result.")
        evtlc=prod_names.evt
    else:
        evtlc=None
        
    mylib.printgreen(f"Extracting the raw light curve '{evtname}' -> '{prod_names.raw}'")
    if not xmmlc_extract_single(evtinfo, gtifile, objreg, prod_names.abs['raw'],
                                timebin, channels, evtlc):
        mylib.printwarn("Some minor errors arose during extraction of the raw light curve. "
                     "Check the result carefully.")
        _no_errors=False
        
    mylib.printgreen(f"Extracting the background light curve '{evtname}' -> '{prod_names.bkg}'")
    if not xmmlc_extract_single(evtinfo, gtifile, bkgreg, prod_names.abs['bkg'],
                                timebin, channels, None):
        mylib.printwarn("Some minor errors arose during extraction of the background light curve. "
                     "Check the result carefully.")
        _no_errors=False
    
    mylib.printgreen("Producing the net light curve ('{0.raw}' - '{0.bkg}') -> '{0.net}'".format(prod_names))
    cmd="epiclccorr srctslist='{raw}' eventlist='{mainevt}' outset='{net}' "\
        "withbkgset=yes bkgtslist='{bkg}' applyabsolutecorrections=yes".\
        format(mainevt=evtinfo.filepath.fspath, **prod_names.abs)
    xmm._call_and_check_result(cmd, prod_names.abs['net'], 'background subtraction',
                               'epiclccorr', _ownname)
    
    # Plot light curves
    if plotlc:
        if not prod_names.img:
            raise KeyError("prod_names.img is not found. The option 'plotlc=True' "
                           "requires the name of the eps-image file to store the result.")
        lcnet_header = fits.getheader(prod_names.net, 1)
        bkgratio=lcnet_header['BKGRATIO']
        if bkgratio>1: bkgratio=1/bkgratio   # TODO this can't be right
        mylib.printgreen(f"Making preview image of the light curves: '{prod_names.img}'")
        try:
            xronos_plot_lcurve(prod_names.abs['raw'], gtifile, prod_names.abs['img'],
                               prod_names.abs['bkg'], bkgratio, prod_names.abs['net'])
        except UnicodeError:
            mylib.printerr("Can't plot the light curve.")
            _no_errors=False
            
    return _no_errors


def xmmlc_combine(lc1: FilePathAbs, lc2: FilePathAbs, reslc: FilePathAbs):
    """Combine two synchronous light curves.
    
    The lcmath task is used. So the resultant light curve will be a bin-to-bin
    sum of the two commands. Arguments are paths of the input ant output files.
    """
    _ownname=mylib.getownname()
    cmd=f"lcmath '{lc1}' '{lc2}' {reslc} 1 1 addsubr=yes"
    xmm._call_and_check_result(cmd, reslc, 'merging of the light curves', 'lcmath', _ownname)
    
    
def _checks_for_xmmlc_make_products(prod_names:ProdNames, clobber):
    prod_names.make_abspaths()
    for ftype, fname in prod_names.abs.items():
        if ftype in ['raw', 'net']:  # These are important files, warn!
            mylib.check_file_not_exist_or_remove(
                fname, override=clobber, action=mylib.Actions.DIE,
                extra_text='Use option --clobber to override it or use --suffix.', 
                remove_warning=True)
        else:  # These are not important files, we can remove it
            mylib.check_file_not_exist_or_remove(
                fname, override=True,
                action=mylib.Actions.WARNING, remove_warning=(not clobber))


def _main():
    import argparse
    
    # Parse the arguments
    parser = argparse.ArgumentParser(
        description="Extract the light curve "
        "with certain time resolution and energy range. Program can combine "
        "light curves from different detectors. Arguments 'evtfile', 'regobj' "
        "and 'regbkg' may be lists of files (pn,mos1,mos2) "
        "but a single GTI file must be common for all the detectors.")
    parser.add_argument('evtfile', nargs=1, help="list of the "
                        "XMM-Newton EVENT-files")
    parser.add_argument('gtifile', nargs=1, help="name of the GTI-file")
    parser.add_argument('regobj', nargs=1, help="name of the object region file")
    parser.add_argument('regbkg', nargs=1, help="name of the background region file")
    parser.add_argument(
        'timebin', nargs='?', type=float, help="time resolution of "
        "the light curve in seconds")
    parser.add_argument('mine', nargs='?', type=float, help="lower limit of the energy range (keV)", 
                        default=0.2)
    parser.add_argument('maxe', nargs='?', type=float, help="higher limit of the energy range (keV)",
                        default=12.0)
    parser.add_argument('-u', '--suffix', nargs='?', help="name suffinx for output files")
    parser.add_argument('-l', '--logfile', nargs='?', help="log file")
    parser.add_argument('--clobber', action='store_true', help="Allow to override files")
    parser.add_argument(
        '--with-evtlc', action='store_true', help="Leave filtered EVT-file with "
        "events used to create the raw light curve")
    
    argnspace=parser.parse_args(sys.argv[1:])
    arg_evtfiles=argnspace.evtfile[0].split(',')
    arg_gtifile=argnspace.gtifile[0]
    arg_regobjfiles=argnspace.regobj[0].split(',')
    arg_regbkgfiles=argnspace.regbkg[0].split(',')
    timebin=argnspace.timebin
    mine=argnspace.mine
    maxe=argnspace.maxe
    logfile=argnspace.logfile
    
    if argnspace.suffix:
        sufx="-"+argnspace.suffix
    else:
        sufx=""

    #############################################################
    #### Check input data
        
    # Check the system variables
    if ("SAS_ODF" not in os.environ) or ("SAS_ODF" not in os.environ):
        mylib.die("'SAS_ODF' and 'SAS_CCF' variables are not defined. Please "
               "define the variables and try again.")
    
    # Check length of lists
    if not ( len(arg_evtfiles) == len(arg_regobjfiles)== len(arg_regbkgfiles) ):
        mylib.die("Lists of the EVENT and the region files must have the same "
               "length and order.")
    nfiles=len(arg_evtfiles)
    
    if nfiles>3:
        mylib.die("You are trying to combine more then 3 light curves. Only "
               "synchronous light curves from EPIC-PN, EPIC-MOS1 and EPIC-MOS2 "
               "can be combined.")
        
    evtdict={}
    regobjdict={}
    regbkgdict={}
    for i, strevtfile in enumerate(arg_evtfiles):  # Check each EVT file
        evtinfo=EVTinfo(xmm.xmm_check_file_is_evt(strevtfile))
        instr=evtinfo.instr_short_name
        if instr in evtdict:        # Check for duplicates
            mylib.die("You are trying to extract and combine light curves from the same "
                   "detector. Can't combine them.")
        evtdict[instr] = evtinfo   # And add to a dict
        regobjdict[instr] = mylib.check_file_exists(arg_regobjfiles[i])
        regbkgdict[instr] = mylib.check_file_exists(arg_regbkgfiles[i])
        
    # Check GTI file exists
    gtifile=fits_check_file_is_gti(arg_gtifile)

    # Check energy range
    if mine<0.2 or mine>12:
        mylib.die("'mine' argument must be in the 0.2-12.0 keV range")
    if maxe<0.2 or maxe>12:
        mylib.die("'maxe' argument must be in the 0.2-12.0 keV range")
    if mine>maxe:
        mylib.die("'mine' must be less than 'maxe'")
    
    channels=(int(mine*1000), int(maxe*1000))
    
    if logfile:
        mylib.logger_turn_on(FilePath(logfile).absolute())
    
    _no_error=True
    if nfiles>1:  # Lcurves wil be combined
        if not timebin:
            mylib.die("If you want to combine the light curves, you must "
                   "explicitly specify the 'timebin' argument.")
            
        prodnames_dict={}
        _kwargs_generate=dict(timebin=timebin, mine=mine, maxe=maxe)
        for instr in evtdict:  # Generate names for output files and check they are already NOT exist
            prodnames_dict[instr] = ProdNames.generate(f'{sufx}_{instr}', **_kwargs_generate)
            _checks_for_xmmlc_make_products(prodnames_dict[instr], argnspace.clobber)
            
        # Trick to generate name of the combined net LC
        comb_nroot='{}_{}'.format(sufx, ''.join(sorted(evtdict)))  # Nameroot for combined lc
        _kwargs_check = dict(override=argnspace.clobber, action=mylib.Actions.DIE,
        extra_text='Use option --clobber to override it or use --suffix.', remove_warning=True)
        lc_combined = mylib.check_file_not_exist_or_remove(
            ProdNames.generate(comb_nroot, **_kwargs_generate).net,
            **_kwargs_check)  # this is pnmos1mos2 if nfiles==3
        if nfiles == 3:
            lc_mos_combined = mylib.check_file_not_exist_or_remove(
                ProdNames.generate(f'{sufx}_mos1mos2', **_kwargs_generate).net,
                **_kwargs_check)   # this is mos1mos2

        #### Do the main job
        mylib.printcaption(f"Making '{lc_combined.name}'...")
        for instr in evtdict:  # Make light curves for each instrument
            if not xmmlc_make_products(
                    evtdict[instr], gtifile, regobjdict[instr],
                    regbkgdict[instr], prodnames_dict[instr], timebin,
                    channels, with_evtlc=argnspace.with_evtlc):
                _no_error=False
        
        if nfiles == 3:
            xmmlc_combine(*[prodnames_dict[x].abs['net'] for x in ('mos1', 'mos2')], reslc=lc_mos_combined)
            xmmlc_combine(lc_mos_combined, prodnames_dict['pn'].abs['net'], lc_combined)
        else:
            lc_pair=[x.net for x in prodnames_dict.values() ]
            xmmlc_combine(*lc_pair, lc_combined)
    else: # Extract single lcurve without combining
        instr, evtinfo = evtdict.popitem()
        if not timebin:
            timebin = max(evtinfo.frmtime.values())
            mylib.printwarn("'timebin' for {} is set to {}".format(instr, timebin))
        prodnames = ProdNames.generate(f'{sufx}_{instr}', timebin, mine, maxe)
        _checks_for_xmmlc_make_products(prodnames, argnspace.clobber)

        mylib.printcaption(f"Making '{prodnames.net}'...")
        if not xmmlc_make_products(
                evtinfo, gtifile, mylib.check_file_exists(arg_regobjfiles[0]),
                mylib.check_file_exists(arg_regbkgfiles[0]), prodnames, timebin,
                channels, with_evtlc=argnspace.with_evtlc):
            _no_error=False
            
    mylib.printcaption("Finished")
    if not _no_error:
        mylib.printwarn("Some minor errors arose. Check the result carefully.")
            
