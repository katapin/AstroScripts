#!/usr/bin/python


import sys
import os
import re
import argparse
from mypython import * 
from myfits import *
from astropy.io import fits
import matplotlib.pyplot as plt

def fitsplot(ftslst,hdulst,colxlst,colylst,pngname,xtitle="",ytitle="",
title=""):

    fig = plt.figure()
    ax = fig.add_subplot(111)
    cl=['b', 'r', 'g', 'c', 'm']         #Color list
    
    for i in range(len(ftslst)):
        fts = fitsopen(ftslst[i])
        X=fts[hdulst[i]].data[colxlst[i]]
        Y=fts[hdulst[i]].data[colylst[i]]
        plt.plot(X, Y,'.-',c=cl[i],label=str(i+1))
    
    plt.legend(loc='upper right');        
    plt.xlabel(xtitle)
    plt.ylabel(ytitle)
    plt.title(title)
    if pngname:
        plt.draw()
        plt.savefig(pngname)
        if not os.path.isfile(pngname):
            printwarn("Cannot save the figurm re into '%s'." % pngname,
            "fitsplot")
            return False
    else:
        plt.show()
    return True


if __name__ == '__main__': 

    #Parse the arguments
    parser = argparse.ArgumentParser(description="Plot FITS tables")
    parser.add_argument('-s', '--save', nargs='?', help="save plot as  "
    "png image")
    parser.add_argument('files', nargs='+', help="FITS file in "
    "a format of 'filname[HDU][ColX,ColY]'")
    argnspace=parser.parse_args(sys.argv[1:])
    pngname=argnspace.save                      #Image filename

    ftslst=[]
    hdulst=[]
    ColXlst=[]
    ColYlst=[]
    #Parse each input arguments
    for arg in argnspace.files:
        # (filename)([HDU])([X,Y])
        match=re.match("(.*?)(\[(\d+)\])?\[(\w+)\,(\w+)\]", arg)
        if not match:
            die("Incorrect syntax","fitsplot")
            
        filename=match.group(1)
        HDU = int(match.group(3)) if match.group(2) else 1
        ColX=match.group(4)
        ColY=match.group(5)
        
        #Open the file 
        fts = fitsopen(filename)
        if len(fts)<=HDU:
            die("There is no HDU number '%d' in %s" % (HDU, filename),
            "fisplot")
                            
        #Check the columns do exist
        cols=fts[HDU].columns
        if not ColX in cols.names:
            die("There is no column '%s' in '%s[%s]'" % (ColX, filename,
            str(HDU)))
        if not ColY in cols.names:
            die("There is no column '%s' in '%s[%s]'" % (ColY, filename,
            str(HDU)))    
        ftslst.append(filename) 
        hdulst.append(HDU)
        ColXlst.append(ColX)
        ColYlst.append(ColY)
        
    if not fitsplot(ftslst,hdulst,ColXlst,ColYlst,pngname,ColX,ColY):
       die("Cannot plot the figure","fitsplot")
    
