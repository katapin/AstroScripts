#!/usr/bin/python
#
#preview chandra data. 

##TODO Добавить возможность анализам мульэксопзиционный файлов

import sys
import os
import argparse
import fileinput
import glob
import uuid
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord
from regions import read_ds9
import matplotlib
import matplotlib.pyplot as plt
import datetime
from mypython import * 
from myfits import *
import gtiplot
from xmmgeneral import *    

class DataIsDamagedError(Exception):
    def __init__(self,obsid,filekey,path,msg=''):
        if not msg:
            msg="ObsID=%s is damaged: %s-file is not found (%s)" % (obsid, filekey, path)
        self.msg=msg
        self.obsid=obsid
        self.filekey=filekey
        self.path=path
    def __str__(self):
        return self.msg
        
class DataIsAbsentError(Exception):
    def __init__(self,obsid,path,msg=''):
        if not msg:
            msg="ObsID=%s is not found in the storage (%s)" % (obsid,path)
        self.msg=msg
        self.obsid=obsid
        self.path=path
    def __str__(self):
        return self.msg
        
class SourceOutsideFOV(Exception):
    def  __init__(self):
        pass
        
class Chandra_observation:
    def __init__(self,obsid,datapath):
        self.obsid=obsid
        self.__datapath=datapath
        self.__repropath=self.__datapath+"/"+self.obsid+"/repro"
        self.__auxilpath=self.__datapath+"/"+self.obsid+"/broad"
        self.files={}
        
        self.__search_for_data()
        self.__load_metadata()
        
    #Load data from evt2 file
    def __search_for_data(self):
        if not os.path.isdir(self.__datapath+"/"+self.obsid):
            raise DataIsAbsentError(self.obsid,self.__datapath)
        
        __nameroot="%05d" % int(self.obsid)
        self.files={'evt2': self.__repropath+"/acisf"+__nameroot+"_repro_evt2.fits.gz",
            'fluximg': self.__auxilpath+"/"+__nameroot+"_broad_flux.img.gz",
            'expmap': self.__auxilpath+"/"+__nameroot+"_broad_thresh.expmap.gz",
            'psfmap': self.__auxilpath+"/"+__nameroot+"_broad_thresh.psfmap.gz"}
            
        for key,val in self.files.items():
            if not os.path.isfile(val):
                raise DataIsDamagedError(self.obsid,key,val)
                
    #Check EVT2 file    
    def __load_metadata(self):
        try:
            with fits.open(self.files['evt2']) as fts:
                self.main_header=fts[0].header
                self.evt_header=fts[1].header
                self.totalevents=len(fts[1].data)
                self.mjd_obs=fts[0].header['MJD-OBS']
                self.date_obs=fts[0].header['DATE-OBS']
                self.exposure=fts[1].header['EXPOSURE']
                self.livetime=fts[1].header['LIVETIME']
                self.datamode=fts[1].header['DATAMODE']
                self.detnam=fts[1].header['DETNAM']
        except Exception:
            raise DataIsDamagedError(self.obsid, 'evt2', self.files['evt2'], "ObsID=%s is "
            "damaged: can't read evt2-file" % self.obsid)
            
    def extract_events(self,expression,outfile):
        command="dmcopy infile=\"%s[%s]\" outfile=\"%s\"" %(self.files['evt2'], expression, outfile)
        if (not callciao(command)) or (not os.path.isfile(outfile)):
            raise Exception("Can't extract events")
        return True
        
    def extract_spatial(self,filekey,outfile,src_reg,bkg_reg=None):
        bkg_str=""      #Extra command to extract background
        if bkg_reg:
            bkg_str="bkg=\"%s[bin sky=region(%s)]\"" % (self.files[filekey], bkg_reg)
        command="dmextract infile=\"%s[bin sky=region(%s)]\" %s outfile=\"%s\" opt=generic" %(self.files[filekey], src_reg, bkg_str, outfile)
        if (not callciao(command)) or (not os.path.isfile(outfile)):
            raise Exception("Can't make products")
        return True
        
    def save_png_image(self,filekey,outfile,region_lst=[]):
        extra_commands=""
        for regfile in region_lst:
            extra_commands+="-regions \"%s\" " % regfile
        fitstoimg(self.files[filekey],outfile,ds9_extra_commands=extra_commands)
            
    def extract_spectrum(self):
        pass
    def extract_curve(self):
        pass

            
class Source_in_obs:
    def __init__(self,obs,src_reg,bkg_reg=""):
        self.src=src_reg
        self.bkg=bkg_reg
        self.obs=obs
        self.obsid=obs.obsid
        
        self.__check_hit_theFOV()
            
    def __check_hit_theFOV(self):    #Does source hit the FOV
        tmpfile=self.__make_tmpfile()
        self.obs.extract_spatial('expmap',tmpfile,self.src)
        with fits.open(tmpfile) as fts:
            expmap_sur_bri=fts[1].data['SUR_BRI']
        os.remove(tmpfile)
        if (expmap_sur_bri>0):
            self.src_effective_exposure=expmap_sur_bri
        else:
            raise SourceOutsideFOV()
            
    def __make_tmpfile(self):
            return "/tmp/"+uuid.uuid4().hex+".fts"
        
    def measure_rates(self):
        tmpfile=self.__make_tmpfile()
        self.obs.extract_spatial('evt2',tmpfile,self.src,self.bkg)
        with fits.open(tmpfile) as fts:
            result_dict={'src_area' : fts[1].data['AREA'][0],
            'raw_counts': fts[1].data['COUNTS'][0],
            'raw_counts_err': fts[1].data['ERR_COUNTS'][0],
            'raw_rate': fts[1].data['COUNT_RATE'][0],
            'raw_rate_err': fts[1].data['COUNT_RATE_ERR'][0],
            'exposure': fts[1].data['EXPOSURE'][0]}
            if self.bkg:
                result_dict.update( {'bkg_area': fts[1].data['BG_AREA'][0],
                'bkg_counts': fts[1].data['BG_COUNTS'][0],
                'bkg_counts_err': fts[1].data['BG_ERR'][0],
                'bkg_rate': fts[1].data['BG_RATE'][0],
                'net_counts': fts[1].data['NET_COUNTS'][0],
                'net_counts_err': fts[1].data['NET_ERR'][0],
                'net_rate': fts[1].data['NET_RATE'][0],
                'net_rate_err': fts[1].data['ERR_RATE'][0]} )
        os.remove(tmpfile)
        self.measured_rates=result_dict
        return result_dict
            
            
        

def do_chandra_preview_object(obsidlist,storage_path,report_path,srcreg,bkgreg):
    obslst=[]
    sor_in_obs=[]
    print("Reading data:")
    for obsid in obsidlist:
        try:
            curobs=Chandra_observation(obsid,storage_path)
            obslst.append(curobs)
            print("ObsID=%05d, Date=%s (MJD=%s) LiveTime=%s" % (int(obsid), curobs.date_obs, 
             curobs.mjd_obs, curobs.livetime))
        except DataIsAbsentError:
            printwarn("ObsID=%05d - Not found" % int(obsid))
        except DataIsDamagedError as ex:
            printwarn("ObsID=%05d - Data is damaged. File %s is not found at path %s" % (int(obsid), ex.filekey, ex.path))
            
            
    obslst.sort(key=lambda x: x.mjd_obs)  #Sort by date      
    
    reportfile=open(report_path+'/report.txt', 'w')
    
    print("#%s" % ' '.join(sys.argv[0:]),file=reportfile)
    print("#ObsID\tDETNAM\t\t\tDate\t\t\tMJD Exposure\t\tRAW_CNT ERR BKG_CNT ERR\t\tBKG_RATE\tNET_RATE NET_RATE_ERR",file=reportfile)        
    for obs in obslst:
        try:
            soinob=Source_in_obs(obs,srcreg,bkgreg)
            sor_in_obs.append(soinob)
            rates=soinob.measure_rates()
            area_ratio=rates['src_area']/rates['bkg_area']
            bkg_cnt=area_ratio*rates['bkg_counts']
            bkg_cnt_err=area_ratio*rates['bkg_counts_err']
            bkg_rate=area_ratio*rates['bkg_rate']
            print("%05d %-12s %s %5.1f %5.1f\t%7d\t%5d\t%1.3f\t%1.3f\t%1.3e\t%1.3e\t%1.3e" % ( int(obs.obsid), obs.detnam, obs.date_obs, float(obs.mjd_obs), float(obs.exposure),
            rates['raw_counts'],rates['raw_counts_err'],bkg_cnt,bkg_cnt_err,bkg_rate,
            rates['net_rate'],rates['net_rate_err']),file=reportfile)
        except SourceOutsideFOV:
            print("#%05d %-12s %s %5.1f %5.1f ======================================== Source outside FOV" % ( int(obs.obsid), obs.detnam,
            obs.date_obs, float(obs.mjd_obs), float(obs.exposure)),file=reportfile)
        png_name=report_path+"/%05d.png" % int(obs.obsid)
        obs.save_png_image('fluximg',png_name,[srcreg,bkgreg])
        
    if len(sor_in_obs)==0: return 1
        
    mjd = [x.obs.mjd_obs for x in sor_in_obs ]
    netrate= [x.measured_rates['net_rate'] for x in sor_in_obs ]   
    netrate_err= [x.measured_rates['net_rate_err'] for x in sor_in_obs ]   
    dates=[x.obs.date_obs for x in sor_in_obs ]
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.xaxis_date()
    ax.errorbar(Time(dates).plot_date,netrate,yerr=netrate_err,fmt='o')
    ax.set_ylabel('Net count rate')
    ax.set_xlabel('Date')
    fig.autofmt_xdate()
    plt.draw()
    plt.savefig(report_path+"/chandra_lcuvrve.png")
       
    reportfile.close()
    
    
def create_ds9reg(filename,sky_coords,radii,key):
    co=sky_coords
    center_str="%d:%d:%2.4f,%d:%d:%2.3f" %(co.ra.hms.h,co.ra.hms.m,co.ra.hms.s,
        co.dec.dms.d,abs(co.dec.dms.m),abs(co.dec.dms.s))
    with open(filename,'w') as freg:
        print("# Region file format: DS9 version 4.1",file=freg)
        print("global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1",file=freg)
        print("fk5",file=freg)
        if (key=='circle'):
            print("circle(%s,%1.4f\")" % (center_str, radii[0]),file=freg)
        else:
            print("annulus(%s,%1.4f\",%1.4f\")" % (center_str, radii[0],radii[1]),file=freg)



if __name__ == '__main__':    
    
    #Parse the arguments
    parser = argparse.ArgumentParser(description="Quick check does "
    "the source flass inside the FOVs of listed ObsIDs ")
    parser.add_argument('obsidlist', nargs=1, help="list of ObsIDs")
    parser.add_argument('-s', '--src-reg', nargs='?', help="source region file")
    parser.add_argument('-b', '--bkg-reg', nargs='?', help="bkg region file")
    parser.add_argument('-c', '--coord', nargs='?', help="source coordinates")
    parser.add_argument('-r', '--rad', nargs='?', help="radii of source and background", default="3,8,24")
    parser.add_argument('--storage-path', nargs='?', default=os.getenv('HOME')+"operative_data/Xray/Chandra/")
    parser.add_argument('--report-dir', nargs='?', help="direcrory to save report")
    
    argnspace=parser.parse_args(sys.argv[1:])
    listfile=argnspace.obsidlist[0]
    srcregfile=argnspace.src_reg
    bkgregfile=argnspace.bkg_reg
    arg_coord_str=argnspace.coord
    radii=[int(x) for x in argnspace.rad.split(',')]
    storage_path=argnspace.storage_path
    report_dir=argnspace.report_dir
    
    if not os.path.isdir(storage_path): die("The storage is not found. "
        "Use '--storage-path' to replace default settings.")
    if not os.listdir(storage_path): die("The storage is empty. "
        "Use '--storage-path' to set valid storage path.")
    
    pid = os.getpid()
    if not srcregfile:
        if not arg_coord_str:
            die("Please set source coordinates or pass region file")
        srcregfile="src_%d.reg" % pid
        create_ds9reg(srcregfile,SkyCoord(arg_coord_str),[radii[0]],'circle')
        
    if not bkgregfile:
        bkgregfile="bkg_%d.reg" % pid
        src_regions_obj = read_ds9(srcregfile)
        create_ds9reg(bkgregfile,src_regions_obj[0].center,radii[1:],'annulus')
        
        
    #Prepare ObsID list    
    obsidtextlist=[];
    try:
        for line in fileinput.input(listfile):
            if not line.split(): continue   #Ignore empty lines
            if line[0] == '#': continue     #Ignore commented lines
            obsidtextlist+=line.strip().split(',')
    except FileNotFoundError:
        die("File %s is not found" % listfile)
    if not obsidtextlist:  die("List is empty, nothing to do")
    
    if not report_dir:
        report_dir="chandra_preview_object-%d" % pid
    os.mkdir(report_dir)

    if do_chandra_preview_object(obsidtextlist,storage_path,report_dir,srcregfile,bkgregfile):
        print("Finished with errors. Report is stored in %s" % report_dir)
    else:
        print("Finished! Report is stored in %s" % report_dir)
    
    if not argnspace.src_reg: os.remove(srcregfile)
    if not argnspace.bkg_reg: os.remove(bkgregfile)
