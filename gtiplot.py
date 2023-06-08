#!/usr/bin/python
#
#Plot GTI files

import sys, os, re
import argparse
import myfits as my
from astropy.io import fits
import matplotlib.pyplot as plt

__progname=os.path.basename(__file__)

def gtiplot(file_hdu_pairs:list, pngname:str=None, lcname:str=None, title:str=""):
    """Plot GTI intervals together with the light curve.

    Parameters
    ----------
    file_hdu_pairs : list
        List of tuples (filename, hdu_number) to precess.
    pngname : str, optional
        Name of a png-file to store the result
    lcname : str, optional
        Path to the light curve
    title : str, optional
        Title of the figure.
    """
    color=['m', 'g', 'c', 'r']         #Color list
    gtilev=1                           #Default GTI level
    
    #Plot light curve
    if lcname:
        with fits.open(lcname) as ftslc:
            time=ftslc[1].data['TIME']
            rate=ftslc[1].data['RATE']
        gtilev=max(rate)*1.05
        plt.plot(time, rate,'b.-')
        plt.ylabel("RATE")
        plt.title(lcname)
        
    #Create rectagular porfiles of GTI
    for i,pair in enumerate(file_hdu_pairs):
        gtifile, hdu = pair
        with fits.open(gtifile) as ftsgti:
            startcol=ftsgti[hdu].data['START']
            stopcol=ftsgti[hdu].data['STOP']
            
        X=[]
        Y=[]
        for k in range(len(startcol)):
            X += [startcol[k], startcol[k], stopcol[k], stopcol[k]]
            Y += [gtilev*0.05*i, gtilev+gtilev*0.05*i, gtilev+gtilev*0.05*i, gtilev*0.05*i]

        plt.plot(X, Y,'.-',c=color[i],label=str(i+1))
        
    plt.ylim(-0.05*gtilev,gtilev+gtilev*0.05*(i+1))
    plt.xlabel("TIME")
    plt.title(title)
    
    if pngname:
        plt.draw()
        plt.savefig(pngname)
        if not os.path.isfile(pngname):
            my.printwarn("Cannot save the figure into '{pngname}'.", __progname)
            raise my.TaskError(__progname, filename=pngname)
    else:
        plt.show()
    return True
   
    
    
def _main():
    #Parse the arguments
    parser = argparse.ArgumentParser(description="Plot Good Time "
    "Intervals (GTI) as П-like functions")
    parser.add_argument('-s', '--save', nargs='?', help="save plot as  "
    "png image")
    parser.add_argument('-r', '--rate', nargs='?', help="plot GTI over "
    "the background light curve")
    parser.add_argument('gtifiles', nargs='+', help="GTI extension in "
    "a format 'filname[HDU]'")
    
    argnspace=parser.parse_args(sys.argv[1:])
    pngname=argnspace.save          #Image filename
    lcname=argnspace.rate           #Light curve

    pairs=[]
    #Parse each input arguments and check them
    for curentry in argnspace.gtifiles:
        curpair = my.fits_parse_expression(curentry)
        if curpair == None:
            my.die("Can't parse input expression. Incorrect syntax")
            
        hdu = int(curpair[1]) if curpair[1]!=None else 1
        filepath=my.fits_check_file_is_gti(curpair[0], hdu)
        pairs.append((filepath, hdu))
        
    if pngname:
        my.check_file_not_exist_or_remove(pngname, overwrite=True)
    
    if lcname:
        my.fits_check_file_is_lc(lcname)
    
    try:
        gtiplot(pairs, pngname, lcname)
    except Exception as ex:
        my.die(f'{ex}. Cannot plot the figure',__progname)
    
if __name__ == '__main__':    
    _main()

    
    
