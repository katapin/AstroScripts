#!/usr/bin/python
#
#Plot GTI files

import sys
import os
import re
import argparse
from mypython import * 
from myfits import *
from astropy.io import fits
import matplotlib.pyplot as plt

def gtiplot(gtilst,hdulst,lcname,pngname,title=""):
    color=['m', 'g', 'c', 'r']         #Color list
    gtilev=1                           #Default GTI level
    
    #Plot light curve
    if lcname:
        ftslc=fitsopen(lcname)
        time=ftslc[1].data['TIME']
        rate=ftslc[1].data['RATE']
        gtilev=max(rate)*1.05
        plt.plot(time, rate,'b.-')
        plt.ylabel("RATE")
        plt.title(lcname)
        
    #Create rectagular porfiles of GTI
    for i in range(len(gtilst)):
        ftsgti = fitsopen(gtilst[i])
        startcol=ftsgti[hdulst[i]].data['START']
        stopcol=ftsgti[hdulst[i]].data['STOP']
        ftsgti.close()
        X=[]
        Y=[]
        for k in range(0, len(startcol)):
            X.append(startcol[k])
            Y.append(gtilev*0.05*i)
            X.append(startcol[k])
            Y.append(gtilev+gtilev*0.05*i)
            X.append(stopcol[k])
            Y.append(gtilev+gtilev*0.05*i)
            X.append(stopcol[k])
            Y.append(gtilev*0.05*i)
           
        #plt.plot(X, Y,'.-',c='m',label=str(i+1))
        plt.plot(X, Y,'.-',c=color[i],label=str(i+1))
        i=i+1
        
    plt.ylim(-0.05*gtilev,gtilev+gtilev*0.05*i)
    plt.xlabel("TIME")
    plt.title(title)
    
    if pngname:
        plt.draw()
        plt.savefig(pngname)
        if not os.path.isfile(pngname):
            printwarn("Cannot save the figurm re into '%s'." % pngname,
            "gtiplot")
            return False
    else:
        plt.show()
    return True
   
    
if __name__ == '__main__':    
    
    #Parse the arguments
    parser = argparse.ArgumentParser(description="Plot Good Time "
    "Intervals (GTI) as П-like functions")
    parser.add_argument('-s', '--save', nargs='?', help="save plot as  "
    "png image")
    parser.add_argument('-r', '--rate', nargs='?', help="plot GTI over "
    "the background light curve")
    parser.add_argument('gtifiles', nargs='+', help="GTI extension in "
    "a format of 'filname[HDU]'")
    
    argnspace=parser.parse_args(sys.argv[1:])
    pngname=argnspace.save          #Image filename
    lcname=argnspace.rate           #Light curve

    gtilst=[]
    hdulst=[]
    #Parse each input arguments and check them
    for arg in argnspace.gtifiles:
        # (filename)([HDU])
        match=re.match("(.*?)(\[(\d+)\])?$", arg)
        if not match:
            print("Error: Incorrect syntax")
        
        filename=match.group(1)
        if not os.path.isfile(filename):
            die("File '%s' is not found" % filename)
        HDU = int(match.group(3)) if match.group(2) else 1
            
        #Open the file 
        ftsgti = fitsopen(filename)
        if len(ftsgti)<=HDU:
            die("There is no HDU number '%d' in %s" % (HDU, filename))
        
        if not chkisgti(ftsgti[HDU]):
            die("'%s' is not a valid GTI extension" % 
            (filename+"["+str(HDU)+"]"),"gtiplot")
        ftsgti.close()
        gtilst.append(filename)
        hdulst.append(HDU)
    
    if lcname:
        ftslc=fitsopen(lcname)
        if not chkislc(ftslc[1]):
            printwarn("'%s' is not a valid light curve extension" % 
            (lcname+"["+str(1)+"]"),"gtiplot")
        ftslc.close()
    
    if not gtiplot(gtilst,hdulst,lcname,pngname):
        die("Cannot plot the figure","gtiplot")
    

    
    
