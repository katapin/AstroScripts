"""The module providing general functions and basic settings for XMM-Newton data processing."""


################################################################
#### Configuration

FILTERING_EXPRESSIONS = dict(
BASIC_FLT_STD_PN     = "#XMMEA_EP && (PATTERN<=4) && (PI in [200:12000])",
BASIC_FLT_STD_MOS    = "#XMMEA_EM && (PATTERN<=12) && (PI in [200:12000])",
BASIC_FLT_STRICT_PN  = "(FLAG==0) && (PATTERN<=4) && (PI in [200:12000])",
BASIC_FLT_STRICT_MOS = "(FLAG==0) && (PATTERN<=12) && (PI in [200:12000])",

FLASH_FLT_STD_PN     = "#XMMEA_EP && (PATTERN==0) && (PI in [10000:12000])",
FLASH_FLT_STD_MOS    = "#XMMEA_EM && (PATTERN==0) && (PI>10000)",
FLASH_FLT_STRICT_PN  = "(FLAG==0) && (PATTERN==0) && (PI in [10000:12000])",
FLASH_FLT_STRICT_MOS = "(FLAG==0) && (PATTERN==0) && (PI>10000)",

FLASH_GTI_PN  = "RATE<=0.4",
FLASH_GTI_MOS = "RATE<=0.35",

SPEC_STD_PN     = "(FLAG==0) && (PATTERN<=4)",
SPEC_STD_MOS    = "#XMMEA_EM && (PATTERN<=12)",
SPEC_STRICT_PN  = "(FLAG==0) && (PATTERN==0)",
SPEC_STRICT_MOS = "(FLAG==0) && (PATTERN==0)"
)

## Spectra extraction hardcoded setting. Don't change these numbers (search  
# the XMM-Newton documentation with a keyword 'specchannelmax'.)

SPEC_CHANNELS_PN  = (0, 20479)
SPEC_CHANNELS_MOS = (0, 11999)
SPEC_BINSIZE_PN   = 5
SPEC_BINSIZE_MOS  = 5

# PATH_WORK_BIN="/home/Kjsdja/work/bin/"
# XMMCALL_SCRIPT=PATH_WORK_BIN+"xmmcall.sh"
# FILELIST="Filelist.lst"

################################################################

import os, re
from enum import StrEnum
from typing import Union, Type

from astropy.io import fits
from astropy.time import Time
from astropy import units as u

from mypythonlib import (
    FilePath, FilePathAbs, TempfilePath, Actions, printerr, printwarn, callandlog,
    printandlog, getownname, logger_turn_on
)

import astroscripts
from astroscripts import ExtPath, ExtPathAbs, TaskError
from astroscripts.external import ExternalTaskError, check_result_file_appeared, fitsimg_to_png


## Enumerations helping to choose the appropriate filtering expression
class FilteringPurpose(StrEnum):
    """Enum helping to choose filtering expression depending on specific goal."""
    
    BASIC  = 'BASIC_FLT'
    FLASH  = 'FLASH_FLT'
    GTI    = 'FLASH_GTI'
    SPEC   = 'SPEC'
    
class FilteringMode(StrEnum):
    """Enum helping to choose filtering expression depending on processing mode."""
    
    standard = 'STD'
    strict   = 'STRICT'
    

INSTRKEYS=['EPN','EMOS1','EMOS2']
INSTRFULLNAMES={'EPN':'EPIC-PN', 'EMOS1':'EPIC-MOS1','EMOS2':'EPIC-MOS2'}
INSTRSHORTNAMES={'EPN':'pn','EMOS1':'mos1','EMOS2':'mos2'}


##############################################
#### Functions to use in child modules #######
def xmm_get_expression(instrkey: str, purpose: FilteringPurpose, mode: FilteringMode):
    """Return expression for specific task."""
    if type(purpose) != FilteringPurpose:
        raise TypeError("'purpose' argument must be FilteringPurpose type")
    if type(mode) != FilteringMode:
        raise TypeError("'mode' argument must be FilteringPurpose type")
    instr = xmm_get_instrument_type(instrkey)
    
    if purpose is FilteringPurpose.GTI:
        return FILTERING_EXPRESSIONS['{}_{}'.format(purpose.value, instr)]
    else:
        return FILTERING_EXPRESSIONS['{}_{}_{}'.format(purpose.value, mode.value, instr)]


def xmm_get_instrument_type(instrkey:str):
    """Convert instrument keyword from FITS header to more convenient format.""" 
    if instrkey == 'EPN':
        instr='PN'
    elif instrkey in ('EMOS1', 'EMOS2'):
        instr='MOS'
    else:
        raise ValueError(f"Unknown instrument type '{instrkey}'")
    return instr


def callxmmsas(cmd, stdin='', separate_logfile='', return_code=False):
    """Call a program withing XMM-SAS, pass input and log result."""
    return callandlog(cmd, stdin=stdin, separate_logfile=separate_logfile,
          extra_start='source ~/.bashrc \nsasinit >/dev/null\n',
          return_code=return_code, progname='callxmmsas')


def _call_and_check_result(cmd: str, targetfile: FilePathAbs, operation: str,
                           xmmtaskname: str, progname: str, separate_logfile: str = '') -> bool:
    if not callxmmsas(cmd, separate_logfile=separate_logfile):
        printerr(f"Cannot perform {operation}: '{xmmtaskname}' finished with errors.", progname)
        raise ExternalTaskError(xmmtaskname, filename=targetfile, caller=progname)
    check_result_file_appeared(targetfile, progname)
    return True


def _make_test_images(
        evtpath: FilePathAbs, outimg_png: FilePathAbs, expression: str = '',
        colX: str = 'X', colY: str = 'Y', outimg_fts: FilePathAbs = None,
        progname: str = None) -> object:
    """Create images filtered with the same GTI and aperture as the main product."""
    imgfts = outimg_fts or TempfilePath.generate_from(outimg_png, new_suffix='.fts')
    try:
        xmm_make_image(evtpath, imgfts, colX, colY, expression)
        fitsimg_to_png(imgfts, outimg_png)
    except TaskError:
        printwarn("Can't create the testimage.", progname)
        return False
    return True
            

def xmm_read_regfile(regfile: FilePathAbs) -> str:
    """Convert content of ds9 region files to filtering expression.

    Only ds9 region format is supported. Regions can be a circle, annulus or box.
    
    :param regfile: Path to ds9 region file.
    :returns: Expression to use with the 'evselect' test.
    """
    arealist=[]
    _ownname=getownname()
    with open(regfile, "r") as freg:
        if freg.readline().strip() != "# Region file format: DS9 version 4.1":
            printerr("Only DS9 version 4.1' region format is supported.", _ownname)
            raise TaskError(_ownname, filename=regfile)
            
        # Skip second line with region file specification
        freg.readline()
        
        # Coordinate system
        if freg.readline().strip()!='physical':
            printerr("Wrong coordinate system. Coorinate system must be 'physical'.",_ownname)
            raise TaskError(_ownname, filename=regfile)
        for line in freg.readlines():
            areastr=line.strip()    # String with area definition
            match=re.match("^(\w+)\((.*)\)$", areastr)
            if match:
                shape=match.group(1)        # Area name
                if shape in ['circle','annulus']:
                    arealist.append("(X,Y) IN %s" % areastr)
                elif shape=='box':
                    bxp=match.group(2).split(',')  # Box parameters
                    arealist.append("(X,Y) IN box(%s,%s,%s,%s,%s)" % (bxp[0],
                    bxp[1],float(bxp[2])/2,float(bxp[3])/2,bxp[4]))
                else:
                    printerr(f"Unknown region shape '{shape}'. Can't read the region file", _ownname)
                    raise TaskError(_ownname, filename=regfile)
            else:
                printerr(f"Cannot parse region file: wrong line '{areastr}'",_ownname)
                raise TaskError(_ownname, filename=regfile)
    return "||".join(arealist)    # Expression to select region


def xmm_make_image(evtfile: FilePathAbs, imgfile: FilePathAbs,
                   colxname: str, colyname: str, expression: str = '',
                   xbin: int = 600, ybin: int = 600, logfile = ''):
    """Create a FITS image from the EVT-file."""
    expression_statement=f"expression='{expression}'" if expression else ''
    cmd="evselect table='{}' imagebinning=imageSize withimageset=yes "\
        "imageset='{}' {} xcolumn='{}' ycolumn='{}' ximagesize={:d} "\
        "yimagesize={:d}".format(evtfile.fspath, imgfile.fspath, expression_statement,
        colxname, colyname, xbin, ybin)

    if _call_and_check_result(cmd, imgfile, 'image creation', 'evselect', getownname(), logfile):
        return True
    return False


def xmm_get_mjdobs(path: os.PathLike | str, time_object: bool = False) -> tuple:
    """Get analogues of the MJD-OBS and MJD-END keywords which is normally absent in the XMM data."""
    with fits.open(path) as fts:
        if astroscripts.fits_check_hdu_is_spectrum(fts[1], action=Actions.NOTHING) or xmm_check_hdu_is_evt(fts[1]):
            ti1 = Time(fts[0].header['DATE-OBS'], format='fits')  # Zero-extension!
            ti2 = Time(fts[0].header['DATE-END'], format='fits')  # Zero-extension!
        elif astroscripts.fits_check_hdu_is_lc(fts[1], silent=True):
            mjdref = fts[1].header['MJDREF']
            tstart = fts[1].header['TSTART']
            tstop = fts[1].header['TSTOP']
            ti1 = Time(mjdref, format="mjd") + tstart * u.s
            ti2 = Time(mjdref, format="mjd") + tstop * u.s
        else:
            raise NotImplementedError(f"Unsupported file type of '{path}'")
    return (ti1, ti2) if time_object else (ti1.mjd, ti2.mjd)
        
        
def xmm_check_file_is_evt(filepath: Union[ExtPath, FilePath, str],
                          action: Actions = Actions.DIE, progname=None) -> ExtPathAbs | None:
    """Check whether the first extension is an XMM event file and return abspath."""
    return astroscripts.main._helper_check_file_is(xmm_check_hdu_is_evt, 'XMM evt-file', filepath,
                                             default_hdu='EVENTS', action=action, progname=progname)


def xmm_check_hdu_is_evt(HDU, action: Union[Actions, str] = Actions.WARNING,
                         exception_class : Type[Exception] = TypeError, **kwargs):
    """Check if the HDU is an EVT-extension."""
    extra_text='Probably this is not a valid XMM-Newton EVENT-file.'
    if 'EXTNAME' not in HDU.header:
        Actions(action).do("There is no 'EXTNAME' keyword in the HDU header. "
                          "Probably this FITS file is wrong.", exception_class=exception_class, **kwargs)
        return False
    if HDU.header['EXTNAME'] != 'EVENTS':
        Actions(action).do("The extenion name is not 'EVENTS'. "+extra_text,
                          exception_class=exception_class, **kwargs)
        return False
    if 'INSTRUME' not in HDU.header:
        Actions(action).do("There is no 'INSTRUME' keyword in the 'EVENTS' header. "+extra_text,
                          exception_class=exception_class, **kwargs)
        return False        
    if 'DATAMODE' not in HDU.header:
        Actions(action).do("There is no 'DATAMODE' keyword in the 'EVENTS' header. "+extra_text,
                          exception_class=exception_class, **kwargs)
        return False
    if 'SUBMODE' not in HDU.header:
        Actions(action).do("There is no 'SUBMODE' keyword in the 'EVENTS' header. "+extra_text,
                          exception_class=exception_class, **kwargs)
        return False
    if HDU.header['INSTRUME'] not in INSTRKEYS:
        Actions(action).do("Unknown instrument '{}'.".format(HDU.header['INSTRUME']),
                          exception_class=exception_class, **kwargs)
        return False
    return True


class EVTinfo:
    """Struct to store most important information from XMM evt-file."""
    
    def __init__(self, evtpath: ExtPathAbs):
        with fits.open(evtpath) as fts:
            header = fts['EVENTS'].header
            instrkey=header['INSTRUME']
            self.chips = set(fts['EVENTS'].data['CCDNR'])
            self.frmtime={}
            for chip in self.chips:
                extname=f'EXPOSU{chip:02d}'
                self.frmtime[chip]=fts[extname].header['FRMTIME']/1000
        self.instrkey=instrkey
        self.instr_full_name=INSTRFULLNAMES[instrkey]
        self.instr_short_name=INSTRSHORTNAMES[instrkey]
        self.datamode = header['DATAMODE']
        self.submode  = header['SUBMODE']
        self.filter   = header['FILTER']
        self.expidstr = header['EXPIDSTR']
        self.telapse  = header['TELAPSE']
        self.times    = xmm_get_mjdobs(evtpath, time_object=True)
        self.filepath = evtpath

    def describe(self):
        """Print summary."""
        _ownname='evtinfo'
        printandlog("Data mode:       {}".format(self.datamode),_ownname)
        printandlog("Data submode:    {}".format(self.submode),_ownname)
        printandlog("Filter:          {}".format(self.filter),_ownname)
        printandlog("Exposure ID:     {}".format(self.expidstr),_ownname)
        printandlog("Exposure length: {:.2f}".format(self.telapse),_ownname)
        printandlog("Start time:      {0.fits} (MJD{0.mjd:.5f})".format(self.times[0]),_ownname)
        printandlog("Stops time:      {0.fits} (MJD{0.mjd:.5f})".format(self.times[1]),_ownname)
        printandlog("Used chips: {}".format(', '.join(str(x) for x in self.chips)),_ownname)
        printandlog("Frame time (sec.): {}".format(', '.join(str(x) for x
                                            in set(self.frmtime.values()))), _ownname)
      
    


