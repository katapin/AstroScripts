#!/usr/bin/python
"""Plot FITS and ASCII tables."""


import sys, argparse
import myfits as my
from astropy.io import fits, ascii
from textwrap import dedent


# def _helper_plot(data:list, pngname=None, progname=None, **kwargs):
#     import matplotlib.pyplot as plt
#     for item in data:
#         plt.errorbar(item['X'], item['Y'], xerr=item['Xerr'], yerr=item['Yerr'],
#              fmt=item['fmt'], color=item['col'], label=item['lab'])
#
#     # plt.legend(loc='upper right');
#     plt.xlabel(kwargs['xlabel'])
#     plt.ylabel(kwargs['ylabel'])
#     plt.title(kwargs['title'])
#     if kwargs['legend']: plt.legend()
#     if pngname:
#         plt.draw()
#         plt.savefig(pngname)
#         _helper_check_result_file_appeared(pngname, progname)
#     else:
#         plt.show()


def tbl_plot(files_and_columns: TBL_extended_path | tuple, xlabel: str = '', ylabel: str = '',
             title: str = '', pngname: str = None):
    """Plot columns from a FITS table or from an ASCII table.

    Parameters
    ----------
    files_and_columns: TBL_extended_path | tuple
        A dataclass instance describing path to the table, column names and table format, a tuple
        of these structs.
    xlabel : str, optional
        Custom X-axis label. If empty, column names will be used.
    ylabel : str, optional
        Custom Y-axis label. If empty, column names will be used.
    title : TYPE, optional
        Title of the plot. The default is ''.
    pngname : str, optional
        Name of a png-file to save the figure. The default is None.
    """
    _ownname = my.getownname()
    col = ['b', 'r', 'g', 'c', 'm']  # Color list

    if type(files_and_columns) != tuple and type(files_and_columns) != list:
        files_and_columns = files_and_columns,

    # Can't make one-column and multiple-columns plots in the same figure
    col_count = {len(x.columns) for x in files_and_columns}
    if min(col_count) == 1 and len(col_count) > 1:
        raise ValueError('Inconsistent number of columns to plot.')

    data = [];
    xdict, ydict = {}, {}
    for i, item in enumerate(files_and_columns):
        # Open  the table file
        if item.format == 'fits':
            table = fits.getdata(item.path, item.hdu)
        else:
            table = ascii.read(item.path, format=item.format)
        # Prepare columns
        Xerr, Yerr = None, None
        match len(item.columns):
            case 0:
                raise IndexError('Empty array: No columns provided.')
            case 1:
                Y = table[item.columns[0]];
                X = range(len(Y));
                xdict[0] = 'Point number'
                ydict[i] = item.columns[0]
            case 2:
                X = table[item.columns[0]];
                Y = table[item.columns[1]];
                xdict[i] = item.columns[0]
                ydict[i] = item.columns[1]
            case 3:
                X = table[item.columns[0]];
                Y = table[item.columns[1]];
                Yerr = table[item.columns[2]]
                xdict[i] = item.columns[0]
                ydict[i] = '{}[{}]'.format(item.columns[1], item.columns[2])
            case 4:
                X = table[item.columns[0]];
                Xerr = table[item.columns[1]];
                Y = table[item.columns[2]];
                Yerr = table[item.columns[3]]
                xdict[i] = '{}[{}]'.format(item.columns[0], item.columns[1])
                ydict[i] = '{}[{}]'.format(item.columns[2], item.columns[3])
            case _:
                raise ValueError('Too many columns to plot.')

        data.append(dict(X=X, Y=Y, Xerr=Xerr, Yerr=Yerr))
        data[-1]['fmt'] = item.fmt if hasattr(item, 'fmt') else '.'
        data[-1]['col'] = item.fmt if hasattr(item, 'col') else col[i]
        data[-1]['lab'] = item.fmt if hasattr(item, 'lab') else str(i + 1)

    # Prepare auto labels
    auto_xlabel = xdict[0] if len(set(xdict.values())) == 1 else str(xdict).replace('{', '').replace('}', '').replace(
        "'", '')
    auto_ylabel = ydict[0] if len(set(ydict.values())) == 1 else str(ydict).replace('{', '').replace('}', '').replace(
        "'", '')
    _helper_plot(data, pngname, _ownname, xlabel=(xlabel or auto_xlabel),
                 ylabel=(ylabel or auto_ylabel), title=title, legend=True if len(data) > 1 else False)

    return True


if __name__ == '__main__': 

    #Parse the arguments
    parser = argparse.ArgumentParser(description="Plot FITS and ASCII tables.",
             formatter_class=argparse.RawDescriptionHelpFormatter,
             epilog=dedent("""\
    Examples:
        file.fts[1]['Col1']                  plot 'Col1' from the HDU #1 versus row number
        file.fts['DATA']['X','Y']          - plot 'Y' versus 'X' from HDU named 'DATA'               
        txtfile.dat['X','Y','Yer']         - plot 'Y' with errorbars versus 'X' from an ASCII file
        file.fts[1]['X','Xerr,'Y','Yer']   - plot 'Y' versus 'X' both with errorbars
        
    Most useful formats of ascii tables: fixed_width, no_header,
    fixed_width_no_header, commented_header, csv, tab, daophot.
    """))
    parser.add_argument('files', nargs='+', help="FITS file in "
    "a format of 'filname[HDU][ColX,ColY]'")
    parser.add_argument('-s', '--save', nargs='?', help="save plot as  "
    "a png image")
    parser.add_argument('-f', '--format', nargs='?', help="file format: 'fits' or "
    "any format supported by astropy.io.ascii. The defualt is 'fits'.", default='fits')
    parser.add_argument('-x', '--xlabel', nargs='?', help="custom X-axis label.", metavar='str')
    parser.add_argument('-y', '--ylabel', nargs='?', help="custom Y-axis label.", metavar='str')
    parser.add_argument('-t', '--title', nargs='?', help="title of the figure.", metavar='str')
    parser.add_argument('-l', '--with-lines', action='store_true', 
        help="connect datapoints with lines")
    
    argnspace=parser.parse_args(sys.argv[1:])
    pngname=argnspace.save                      #Image filename

    files_to_plot=[]
    #Parse each input arguments
    for curentry in argnspace.files:
        result = my.parse_extended_path(curentry, format=argnspace.format)
        if result == None:
            my.die(f"Can't parse the input expression '{curentry}'. Incorrect syntax")
            
        if not result.columns:
            my.die(f"No column names passed in the expression '{curentry}'.")
            
        if len(result.columns)>4:
            my.die('Too many columns to plot.')
            
        #Check the file exists
        if argnspace.format == 'fits':
            try:
                result.path=my.fits_check_file_is_fits(result.path, result.hdu, action=my.Action.EXCEPTION)
            except TypeError:  #File is not a FITS, may be it is an ascii table
                my.die(f"It seems the file '{result.path}' you're trying to read is not a FITS file. Use the '--format' "\
                       "argument to open an ASCII table.")
            all_columns = fits.getdata(result.path, result.hdu).columns
        else:
            result.path=my.check_file_exists(result.path)
            try:   #Wrong format, may be it's a typo
                tbl=ascii.read(result.path,format=argnspace.format)
            except ValueError as ex:  #print message from astropy and exit
                my.die(str(ex))
            all_columns = tbl.columns
        
        #Check if columns really exist
        for col in result.columns:
            if col not in all_columns:
                my.die("Column '{}' is not found in {}. Might be the chosen "\
                       "table format was incorrect".format(col, 
                           result.path))
        result.fmt='.-' if argnspace.with_lines else '.'        
        files_to_plot.append(result)  #All if OK, add the file
        
    if pngname:
        pngname=my.check_file_not_exist_or_remove(pngname, overwrite=True)
        
    try:
        my.tbl_plot(files_to_plot, xlabel=argnspace.xlabel, ylabel=argnspace.ylabel,
                    title=argnspace.title, pngname=pngname)
    except Exception as ex:
        my.die(f'{ex}. Cannot plot the figure.')
                

