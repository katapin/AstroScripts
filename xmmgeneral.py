#This is module providing general functions and basic settings for 
#automatical XMM-Newton data processing

#Configuration
################################################################

BASIC_FLT_EXP_PN="#XMMEA_EP && (PATTERN<=4) && (PI in [200:12000])"
BASIC_FLT_EXP_MOS="#XMMEA_EM && (PATTERN<=12) && (PI in [200:12000])"
FLASH_FLT_EXP_PN="#XMMEA_EP && (PATTERN==0) && (PI in [10000:12000])"
FLASH_FLT_EXP_MOS="#XMMEA_EM && (PATTERN==0) && (PI>10000)"

#BASIC_FLT_EXP_PN="(FLAG==0) && (PATTERN<=4) && (PI in [200:12000])"
#BASIC_FLT_EXP_MOS="(FLAG==0) && (PATTERN<=12) && (PI in [200:12000])"
#FLASH_FLT_EXP_PN="(FLAG==0) && (PATTERN==0) && (PI in [10000:12000])"
#FLASH_FLT_EXP_MOS="(FLAG==0) && (PATTERN==0) && (PI>10000)"

FLASH_GTI_EXP_PN="RATE<=0.4"
FLASH_GTI_EXP_MOS="RATE<=0.35"

PATH_WORK_BIN="/home/Kjsdja/work/bin/"
XMMCALL_SCRIPT=PATH_WORK_BIN+"xmmcall.sh"
FILELIST="Filelist.lst"

################################################################

import sys
import re
from astropy.io import fits
from mypython import *
from myfits import *
import subprocess

INSTRKEYS=['EPN','EMOS1','EMOS2']
INSTRFULLNAMES={'EPN':'EPIC-PN', 'EMOS1':'EPIC-MOS1','EMOS2':'EPIC-MOS2'}
INSTRSHORTNAMES={'EPN':'pn','EMOS1':'mos1','EMOS2':'mos2'}

BASIC_FLT_EXP={'EPN':BASIC_FLT_EXP_PN,'EMOS1':BASIC_FLT_EXP_MOS,
               'EMOS2':BASIC_FLT_EXP_MOS}
FLASH_FLT_EXP={'EPN':FLASH_FLT_EXP_PN,'EMOS1':FLASH_FLT_EXP_MOS,
               'EMOS2':FLASH_FLT_EXP_MOS}
FLASH_GTI_EXP={'EPN':FLASH_GTI_EXP_PN,'EMOS1':FLASH_GTI_EXP_MOS,
               'EMOS2':FLASH_GTI_EXP_MOS}

def callxmm(cmd,logfile=""):

    printbold(cmd,"callxmm")
    if logfile:
        res=subprocess.call("source ~/.bashrc \nsasinit " \
        ">/dev/null\n"+cmd+" | tee -a "+logfile,shell=True)
    else:
        res=subprocess.call("source ~/.bashrc \nsasinit >/dev/null\n"+
        cmd,shell=True)
    return not res

def readregfile(regfile):
    ##Only ds9 region format supported!
    
    arealist=[]
    freg=open(regfile, "r")
    if freg.readline().strip()!= \
    "# Region file format: DS9 version 4.1":
        printerr("Only DS9 version 4.1' region format is supported.", \
        "readregfile")
        return None
    #Skip second line with region file specification
    freg.readline()
    #Coordinate system
    if freg.readline().strip()!='physical':
        printerr("Wrong coordinate system. Coorinate system must be "
        "'physical'.", "readregfile")
        return None
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
                printerr("Unknown region shape '%s'" % shape, 
                "readregfile")
                return None
        else:
            printerr("Cannot parse region file: wrong line '%s'" \
             % areastr, "readregfile")
            return None
    freg.close()
    return "||".join(arealist)    #Expression to select region    
    
def xmmimagemake(evtfile,imgfile,colxname,colyname,expression="", \
xbin=600,ybin=600):
    
    if expression:
        COMMAND="evselect table='%s' imagebinning=imageSize " \
        "withimageset=yes imageset='%s' expression='%s' xcolumn=%s " \
        "ycolumn=%s ximagesize=%s yimagesize=%s" % (evtfile, imgfile, 
        expression, colxname, colyname, xbin, ybin)
    else:
        COMMAND="evselect table='%s' imagebinning=imageSize " \
        "imageset='%s' withimageset=yes xcolumn=%s ycolumn=%s " \
        "ximagesize=%s yimagesize=%s" % (evtfile, imgfile, colxname, 
        colyname, xbin, ybin)
    if not callxmm(COMMAND):
        printerr("Cannot create an image: 'evselect' finished "
        "with error.", "makeimage")
        return False
    if not os.path.isfile(imgfile):
        printerr("Something is going wrong: %s is not "
        "created." % imgfile,"makeimage")
        return False
    return True

        
def xmmchkisevt(fts):
    try:
        fts['EVENTS']
    except KeyError:
        printwarn("There is no 'EVENTS' extension in this FITS file. "
        "Probably this is not a valid XMM-Newton EVENT-file.","chkisevt")
        return False
    HDU=fts['EVENTS']
    if not chkkey(HDU,'INSTRUME'):
        printwarn("There is no 'INSTRUME' keyword in the 'EVENTS' header. "
        "Probably this is not a valid XMM-Newton EVENT-file.","chkisevt")
        return False        
    if not chkkey(HDU,'DATAMODE'):
        printwarn("There is no 'DATAMODE' keyword in the 'EVENTS' header. "
        "Probably this is not a valid XMM-Newton EVENT-file.","chkisevt")
        return False
    if not chkkey(HDU,'SUBMODE'):
        printwarn("There is no 'SUBMODE' keyword in the 'EVENTS' header. "
        "Probably this is not a valid XMM-Newton EVENT-file.","chkisevt")
        return False
    if HDU.header['INSTRUME'] not in INSTRKEYS:
        printwarn("Unknown instrument '%s'." % HDU.header['INSTRUME'],
        "chkisevt")
        return False
    return True

class evtinfo:
    def __init__(self, HDU):
        instrkey=HDU.header['INSTRUME']
        self.instrkey=instrkey
        self.instr_full_name=INSTRFULLNAMES[instrkey]
        self.instr_short_name=INSTRSHORTNAMES[instrkey]
        self.basic_flt_exp=BASIC_FLT_EXP[instrkey]
        self.flash_flt_exp=FLASH_FLT_EXP[instrkey]
        self.flash_gti_exp=FLASH_GTI_EXP[instrkey]
        self.datamode=HDU.header['DATAMODE']
        self.submode=HDU.header['SUBMODE']
        self.filter=HDU.header['FILTER']
        self.expidstr=HDU.header['EXPIDSTR']
        self.telapse=HDU.header['TELAPSE']

def xmmgetevtinfo(fts):
    
    #Is fts a filename or a file descriptor?
    if type(fts) == str:
        f=fits.open(fts)
        info=evtinfo(f['EVENTS'])
        f.close()
    else:
        info=evtinfo(fts['EVENTS'])
    return info
    

      
    


