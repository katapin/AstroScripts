#!/usr/bin/env python
"""Plot Good Time Intervals (GTI) as П-like functions."""

import sys
from astropy.io import fits

import mypythonlib as mylib

from .main import ExtPath, ExtPathAbs, FilePathAbs, \
    fits_check_file_is_lc, fits_check_file_is_gti
from .plot import PlotPair, plot_and_save


def gtiplot(files: ExtPathAbs | tuple[ExtPathAbs] | list[ExtPathAbs],
            lcpath: FilePathAbs = None, title: str = '', pngpath: FilePathAbs = None,
            onlysave: bool = False) -> bool:
    """Plot GTI intervals together with the light curve.

    Parameters
    ----------

    files: TBL_extended_path | tuple
        A dataclass instance describing path to the table, column names and table format, a tuple
        of these structs.
    lcpath : str, optional
        Path to the light curve
    title : str, optional
        Title of the figure.
    pngpath : str, optional
        Name of a png-file to save the result. The default is None.
    onlysave
    """
    _ownname = mylib.getownname()
    color = ['m', 'g', 'c', 'r']  # Color list
    gtilev = 1  # Default GTI level

    if not isinstance(files, (tuple, list)):
        files = files,  # Make a tuple

    data, pltopts = [], {}
    _title = None

    # Plot light curve
    if lcpath:
        with fits.open(lcpath) as ftslc:
            rate = ftslc[1].data['RATE']
            data.append(PlotPair(ftslc[1].data['TIME'], ftslc[1].data['Rate'],
                                 pltopts=dict(fmt='.-', color='b', label=None)))
        gtilev = max(rate) * 1.05
        pltopts['ylabel'] = "RATE"
        _title = lcpath.name

    # Create rectangular profiles of GTI
    for i, item in enumerate(files):
        with fits.open(item) as ftsgti:
            startcol = ftsgti[item.hdu].data['START']
            stopcol = ftsgti[item.hdu].data['STOP']

        X, Y = [], []
        for curstart, curstop in zip(startcol, stopcol):
            X += [curstart, curstart, curstop, curstop]
            Y += [gtilev * 0.05 * i, gtilev + gtilev * 0.05 * i, gtilev + gtilev * 0.05 * i, gtilev * 0.05 * i]
        data.append(PlotPair(X, Y, pltopts=dict(fmt='.-', color=color[i], label=str(i + 1))))
        # pp.plot()

        # data.append(dict(X=X, Y=Y, Xerr=None, Yerr=None, fmt='.-', col=color[i], lab=str(i + 1)))

    pltopts['xlabel'] = 'Time'
    pltopts['title'] = title or _title
    if len(files) > 1: pltopts['legend'] = None
    plot_and_save(data, pngpath, onlysave=onlysave,
                  progname=_ownname, **pltopts)
    return True

    
if __name__ == '__main__':
    import argparse
    # Parse the arguments
    parser = argparse.ArgumentParser(
        description="Plot Good Time Intervals (GTI) as П-like functions")
    parser.add_argument(
        '-s', '--save', nargs='?', help="save plot as png image")
    parser.add_argument(
        '-r', '--rate', nargs='?', help="plot GTI over the background light curve")
    parser.add_argument(
        'gtifiles', nargs='+', help="GTI extension in a format 'filename[HDU]'")
    
    argnspace = parser.parse_args(sys.argv[1:])
    pngname = argnspace.save          # Image filename
    lcname = argnspace.rate           # Light curve

    files_to_plot = []
    # Parse each input arguments and check them
    for curpath in argnspace.gtifiles:
        files_to_plot.append(fits_check_file_is_gti(ExtPath(curpath)))

    if pngname:
        pngname = mylib.check_file_not_exist_or_remove(pngname, overwrite=True)
    
    if lcname:
        lcname = fits_check_file_is_lc(lcname)  # Change from None to FilePathAbs
    
    try:
        gtiplot(files_to_plot, lcpath=lcname, pngpath=pngname)
    except Exception as ex:
        mylib.die(f'{ex}. Cannot plot the figure.')
    
    
