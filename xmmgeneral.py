"""The module providing general functions and basic settings for XMM-Newton data processing."""


#### Configuration
################################################################

## Filteing expression

BASIC_FLT_STD_PN="#XMMEA_EP && (PATTERN<=4) && (PI in [200:12000])"
BASIC_FLT_STD_MOS="#XMMEA_EM && (PATTERN<=12) && (PI in [200:12000])"
BASIC_FLT_STRICT_PN="(FLAG==0) && (PATTERN<=4) && (PI in [200:12000])"
BASIC_FLT_STRICT_MOS="(FLAG==0) && (PATTERN<=12) && (PI in [200:12000])"

FLASH_FLT_STD_PN="#XMMEA_EP && (PATTERN==0) && (PI in [10000:12000])"
FLASH_FLT_STD_MOS="#XMMEA_EM && (PATTERN==0) && (PI>10000)"
FLASH_FLT_STRICT_PN="(FLAG==0) && (PATTERN==0) && (PI in [10000:12000])"
FLASH_FLT_STRICT_MOS="(FLAG==0) && (PATTERN==0) && (PI>10000)"

FLASH_GTI_PN="RATE<=0.4"
FLASH_GTI_MOS="RATE<=0.35"

SPEC_STD_PN="(FLAG==0) && (PATTERN<=4)"
SPEC_STD_MOS="#XMMEA_EM && (PATTERN<=12)"
SPEC_STRICT_PN="(FLAG==0) && (PATTERN==0)"
SPEC_STRICT_MOS="(FLAG==0) && (PATTERN==0)"

## Spectra extraction hardcoded setting. Don't change these numbers (search  
#the XMM-Newton documentation with a keyword 'specchannelmax'.) 

SPEC_CHANNELS_PN  = (0, 20479)
SPEC_CHANNELS_MOS = (0, 11999)
SPEC_BINSIZE_PN   = 5
SPEC_BINSIZE_MOS  = 5

# PATH_WORK_BIN="/home/Kjsdja/work/bin/"
# XMMCALL_SCRIPT=PATH_WORK_BIN+"xmmcall.sh"
# FILELIST="Filelist.lst"

################################################################

import os, re
import myfits as my
from astropy.io import fits
from enum import StrEnum


## Enumerations for helping to choose the appropriate filtering expression
class FilteringPurpose(StrEnum):
    """Enum helping to choose filtering expression depending on specific goal."""
    
    BASIC  = 'BASIC_FLT'
    FLASH  = 'FLASH_FLT'
    GTI    = 'FLASH_GTI'
    SPEC   = 'SPEC'
    
class FilteringMode(StrEnum):
    """Enum helping to choose filtering expression depedning on processing mode."""
    
    standard  = 'STD'
    strict = 'STRICT'
    

INSTRKEYS=['EPN','EMOS1','EMOS2']
INSTRFULLNAMES={'EPN':'EPIC-PN', 'EMOS1':'EPIC-MOS1','EMOS2':'EPIC-MOS2'}
INSTRSHORTNAMES={'EPN':'pn','EMOS1':'mos1','EMOS2':'mos2'}

def xmm_get_expression(instrkey, purpose:FilteringPurpose, mode:FilteringMode):
    """Return exporession for specific task."""
    if type(purpose) != FilteringPurpose:
        raise TypeError("'purpose' argument must be FilteringPurpose type")
    if type(mode) != FilteringMode:
        raise TypeError("'mode' argument must be FilteringPurpose type")
    instr = xmm_get_instrument_type(instrkey)
    
    if purpose != FilteringPurpose.GTI:
        return globals()['{}_{}_{}'.format(purpose.value, mode.value, instr)]
    else:
        return globals()['{}_{}'.format(purpose.value, instr)]
        
def xmm_get_instrument_type(instrkey):
    """Convert instrument keyword from FITS header to more convenient format.""" 
    if instrkey == 'EPN':
        instr='PN'
    elif instrkey in ('EMOS1', 'EMOS2'):
        instr='MOS'
    else:
        raise ValueError(f"Unknown instrument type '{instrkey}'")
    return instr
 

def callxmm(cmd, stdin='', separate_logfile='', return_code=False):
    """Call Chandra ciao programm, pass input and log result."""
    return my.callandlog(cmd, stdin=stdin, separate_logfile=separate_logfile, 
          extra_start='source ~/.bashrc \nsasinit >/dev/null\n',
          return_code=return_code, progname='callxmm')   


def __call_and_check_result(cmd, targetfile, operation, xmmtaskname, progname, separate_logfile=''):
    if not callxmm(cmd, separate_logfile=separate_logfile):
        my.printerr(f"Cannot perform {operation}: '{xmmtaskname}' finished with errors.", progname)
        raise my.ExternalTaskError(xmmtaskname, filename=targetfile, caller=progname)
    if not os.path.exists(targetfile):
        my.printerr(f"Something is going wrong: '{targetfile}' has not been created.", progname)
        raise my.TaskError(progname, filename=targetfile)
    my.printbold("Saved "+targetfile, progname)
    return True    

def __make_test_images(evtpath, outimg_png, expression='', colX='X', colY='Y' , outimg_fts='', progname=None):
    imgfts = outimg_fts  or my.TempFile.generate_from(outimg_png, new_extension='fts')
    try:
        xmm_make_image(evtpath, imgfts, colX, colY, expression)
        if my.fitsimg_to_png(imgfts, outimg_png):
            my.printbold(f"Image saved in '{outimg_png}'", progname)
    except my.TaskError:
        my.printwarn("Can't create the testimage.", progname)
    return True
            

def xmm_read_regfile(regfile: str) -> str:
    """Convert content of ds9 region files to filtering expression.

    Only ds9 region format is supported. Regions can be circle, annulus or box.
    
    Parameters
    ----------
    regfile : str
        Path to ds9 region file.

    Returns
    -------
    str
        Expression to use with the 'evselect' test.
    """
    arealist=[]
    _ownname=my.getownname()
    with open(regfile, "r") as freg:
        if freg.readline().strip() != "# Region file format: DS9 version 4.1":
            my.printerr("Only DS9 version 4.1' region format is supported.", _ownname)
            raise my.TaskError(_ownname, filename=regfile)
            
        #Skip second line with region file specification
        freg.readline()
        
        #Coordinate system
        if freg.readline().strip()!='physical':
            my.printerr("Wrong coordinate system. Coorinate system must be 'physical'.",_ownname)
            raise my.TaskError(_ownname, filename=regfile)
        for line in freg.readlines():
            areastr=line.strip()    #String with area definition
            match=re.match("^(\w+)\((.*)\)$", areastr)
            if match:
                shape=match.group(1)        #Area name
                if shape in ['circle','annulus']:
                    arealist.append("(X,Y) IN %s" % areastr)
                elif shape=='box':
                    bxp=match.group(2).split(',') #Box parameters
                    arealist.append("(X,Y) IN box(%s,%s,%s,%s,%s)" % (bxp[0],
                    bxp[1],float(bxp[2])/2,float(bxp[3])/2,bxp[4]))
                else:
                    my.printerr(f"Unknown region shape '{shape}'. Can't read the region file", _ownname)
                    raise my.TaskError(_ownname, filename=regfile)
            else:
                my.printerr(f"Cannot parse region file: wrong line '{areastr}'",_ownname)
                raise my.TaskError(_ownname, filename=regfile)
    return "||".join(arealist)    #Expression to select region    
    
def xmm_make_image(evtfile, imgfile, colxname:str, colyname:str, expression:str='', xbin:int=600, 
                 ybin:int=600, logfile=''):
    """Create a FITS image from the EVT-file."""
    expression_statement=f"expression='{expression}'" if expression else ''
    cmd="evselect table='{}' imagebinning=imageSize withimageset=yes "\
        "imageset='{}' {} xcolumn='{}' ycolumn='{}' ximagesize={:d} "\
        "yimagesize={:d}".format(evtfile, imgfile, expression_statement, 
        colxname, colyname, xbin, ybin)
        
    if __call_and_check_result(cmd, imgfile, 'image creation', 'evselect', my.getownname(), logfile):
        return True
    
    return False
        
        
def xmm_check_file_is_evt(filepath, hdu=1, action=my.Action.DIE, progname=None):
    """Check whether the first extension is an XMM event file and return abspath."""
    return my._helper_check_file_is(xmm_check_hdu_is_evt, 'XMM evt-file', filepath, 
                hdu=hdu, action=action, progname=progname)
    
def xmm_check_hdu_is_evt(HDU):
    """Check if the HDU is an EVT-extetnsion."""
    if 'INSTRUME' not in HDU.header:
        my.printwarn("There is no 'INSTRUME' keyword in the 'EVENTS' header. "
        "Probably this is not a valid XMM-Newton EVENT-file.")
        return False        
    if 'DATAMODE' not in HDU.header:
        my.printwarn("There is no 'DATAMODE' keyword in the 'EVENTS' header. "
        "Probably this is not a valid XMM-Newton EVENT-file.")
        return False
    if 'SUBMODE' not in HDU.header:
        my.printwarn("There is no 'SUBMODE' keyword in the 'EVENTS' header. "
        "Probably this is not a valid XMM-Newton EVENT-file.")
        return False
    if HDU.header['INSTRUME'] not in INSTRKEYS:
        my.printwarn("Unknown instrument '{}'.".format(HDU.header['INSTRUME']),
        "chkisevt")
        return False
    return True

class EVTinfo:
    """Struct to store most important information from XMM evt-file."""
    
    def __init__(self, evtpath):
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
        self.filepath = evtpath
      
    


