#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 10 11:54:52 2023

@author: kjsdja
"""

import sys, argparse
import xmmgeneral as xmm
import myfits as my


def _main():
    parser = argparse.ArgumentParser(description="Print basic properties of "
    "XMM-Newton EVENT-files.")
    parser.add_argument('evtfile', nargs='+', help="names of the "
    "XMM-Newton EVENT-files")
    
    argnspace=parser.parse_args(sys.argv[1:])
    evtfiles=argnspace.evtfile
    
    for file in evtfiles:
        filepath=xmm.xmm_check_file_is_evt(file)
        my.printgreen("{}:".format(filepath.basename))
        evtinfo = xmm.EVTinfo(file)
        evtinfo.show()
        
if __name__ == '__main__':    
    _main()
