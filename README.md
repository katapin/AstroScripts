# AstroScripts
My collection of scripts for the analysis of astronomical data. The package provides 
a common python library and a number of executable files divided into several sets:
### xmm
Scripts for processing XMM-Newton/EPIC observations:

``xmminfo`` shows short info about an EVT-file\
``xmmflt``  perform basic filtering of the EVT-file\
``xmmlc``   extracts a light curve from the EVT-file\
``xmmspec`` extracts a spectrum

### swift
Scripts for Swift/XRT observations:

``swift_wget`` downloads from the NASA archive a preprocessed EVT2-file together 
with auxiliary files needed for the expomap generation\
``swift_wget_uk`` downloads the same files from the UK Swift datacentre's archive
(processed with more fresh software)\
``swift_archive_clone`` clone archive files for the specific ObsID\
``swift_preprocess`` preprocess the downloaded archive data on the local machine\
``swift_expmap`` produces the combined expomap

### chandra 
Scripts for Chandra observations:

``chandra_proper_download`` wrapper for the download_chandra_obsid task 
which allows to download observations in conditions of a poor
internet connection. Also, it can use proxy via the proxychains program.

### fits
``gtiplot.py`` plots Good Time Interval files as П-profiles together with
a light curve. It's useful for checks of the GTI-files produced by the filtering
procedure aimed at rejection of the moments with high particle background in 
XMM-Newton observations.

``regview.py`` iteratively shows FITS images together with the region file 
passed. Also, one can pass common scale and zoom options. 

## Requirements
 - numpy
 - astropy
 - matplotlib
 - scipy 
 - pyds9 

## Installation and usage
Preferred installation method is the use of [pipx](https://pipx.pypa.io/stable/).
To activate a specific collection together with the package virtual environment
execute ``astroscripts`` with the collection name argument:
```
$ pipx intstall https://github.com/katapin/AstroScripts.git
$ astroscipts <collection-name>
```


