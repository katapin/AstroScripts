#!/usr/bin/python3

"""View multiple FITS images in ds9."""

import os, sys, argparse
import base, ds9 



#### Defaults ####
ARG_DEFAULTS_CMAP='bb'
ARG_DEFAULTS_SCALE='log'
ARG_DEFAULTS_TITLE='NBI'



def show_iteratively(ds9obj, filelist:list, *, preprocessor=None, asker=None,
        args_ds9_load_region:list=None):
    for file in filelist:
        try:
            ds9.load_image(ds9obj, file)
        except ds9.DS9Error:
            base.printerr(f"Can't load the image '{file}'")
        

def _main():
    parser = argparse.ArgumentParser(description="Iteratively shows FITS images in ds9")
    parser.add_argument('input', nargs='+', help="One or multiple FITS images or ASCII list")
    parser.add_argument('--cmap', nargs='?', help="Name of the colormap in ds9", default=ARG_DEFAULTS_CMAP)
    parser.add_argument('--scale', nargs='?', help="Scale parameter of ds9", default=ARG_DEFAULTS_SCALE)
    parser.add_argument('--region', nargs='*', help="Load region file(s) for every image")
    parser.add_argument('--region-another', nargs='*', help="Plot region(s) with a different color")
    parser.add_argument('--title', nargs='?', help="Title of the ds9 window. It's "
        "possible to connect to the existing session. The title will be random if "
        "this argument is empty.", default=ARG_DEFAULTS_TITLE)
    parser.add_argument('--ds9-once', nargs='?', help="Extra arguments for ds9, run at start")
    parser.add_argument('--ds9-per-file', nargs='?', help="Extra arguments for ds9, run for each file")
    
    argnspace=parser.parse_args(sys.argv[1:])
    
    if len(argnspace.input)==1: #try to guess type of the only infile 
        infile = base.check_file_exists_or_die(argnspace.input[0])
        if base.check_file_is_FITS(infile):  #The file is a FITS
            filelist=[infile]
        else:                                    #The file is an ASCII list
            try: 
                filelist=base.read_ascii(infile)
            except UnicodeError:  #It's not an ASCII list...
                base.die(f"Can't determine type of the input file '{infile}'")
            if len(filelist) == 0:
                print(f"The list '{infile}' is empty. Noting to do. Exiting...")
                sys.exit(0)
            
    
    ds9obj=ds9.start(title=argnspace.title, cmap=argnspace.cmap, scale=argnspace.scale)
    show_iteratively(ds9obj, filelist)

if __name__ == '__main__':    
    _main()
