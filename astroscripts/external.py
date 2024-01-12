"""Provide facilities for calling external programs."""

import mypythonlib as mylib
from mypythonlib import FilePathAbs, check_file_not_exist_or_remove
from .extpath import ExtPathAbs
from .exceptions import TaskError, ExternalTaskError
from .other import gti_get_limits


def callshell(cmd: str, stdin: str = '', separate_logfile: FilePathAbs = '',
              return_code: bool = False) -> bool | tuple[bool, int]:
    """Call shell program, pass the input and log the result."""
    return mylib.callandlog(cmd, stdin=stdin, separate_logfile=separate_logfile,
                            return_code=return_code, progname='callshell')


def callftools(cmd: str, stdin: str = '', separate_logfile: FilePathAbs = '',
               return_code: bool = False) -> bool | tuple[bool, int]:
    """Call ftools program, pass the input and log the result."""
    return mylib.callandlog(cmd, stdin=stdin, separate_logfile=separate_logfile,
                            extra_start='source ~/.bashrc \nheainit\n',
                            return_code=return_code, progname='callftools')


def callciao(cmd: str, stdin: str = '', separate_logfile: FilePathAbs = '',
             return_code: bool = False) -> bool | tuple[bool, int]:
    """Call Chandra's ciao program, pass the input and log the result."""
    return mylib.callandlog(cmd, stdin=stdin, separate_logfile=separate_logfile,
                            extra_start='source ~/.bashrc \nciaoinit >/dev/null\n',
                            return_code=return_code, progname='callciao')


def _check_result_file_appeared(filepath: FilePathAbs, progname: str) -> None:
    if not filepath.exists():
        mylib.printerr(f"Something is going wrong: '{filepath}' has not been created.", progname)
        raise TaskError(progname, filename=filepath)
    mylib.printbold(f"Saved '{filepath.name}'", progname)


#### Some task associated with external programs


def fitsimg_to_png(ftspath: FilePathAbs, imgpath: FilePathAbs = None, ds9_extra_commands: str = '',
                   convert_extra_commands: str = '') -> bool:
    """Convert FITS image to png with ds9."""
    _ownname = mylib.getownname()
    tmpimg = mylib.TempfilePath.generate('.ps')
    if imgpath is None: imgpath = ftspath.with_suffix('png')

    if not callshell(
            "ds9 {} -nopanner -nomagnifier -noinfo -colorbar no -view buttons no -geometry 1360x768 " 
            "-zoom to fit -cmap bb -scale log exp 10000 -scale log {} -print destination file " 
            "-print filename {} -print -exit".format(ftspath, ds9_extra_commands, tmpimg)):
        mylib.printerr("Cannot call 'ds9'.", _ownname),
        raise ExternalTaskError('ds9', filename=ftspath, caller=_ownname)

    if not callshell(
            "convert -density 300 {} {} -background white " 
            "-flatten -fuzz 5%% -trim {}".format(tmpimg, convert_extra_commands, imgpath)):
        mylib.printerr(f"Cannot call 'convert', '{imgpath.name}' is not created.", _ownname)
        raise ExternalTaskError('convert', filename=imgpath.name, caller=_ownname)

    _check_result_file_appeared(imgpath, _ownname)
    return True


def xronos_plot_lcurve(
        lcraw: FilePathAbs, gtitable: ExtPathAbs, outepsimg: FilePathAbs,
        lcbkg: FilePathAbs = None, bkgratio: float = None, lcnet: FilePathAbs = None,
        binsize: int = 500):
    """Make an eps image of the light curve using the 'lcurve' task from XRONOS.

    Parameters
    ----------
    lcraw : str
        Path to the raw object light curve (extracted from the source aperture
        but not background subtracted).
    gtifile_with_hdu : str
        GTI file to obtain X-axis limits for the plot. Maight in a form of gtifile.fts[#hdu]
    outepsimg : str
        Path to eps file to store the result
    lcbkg : str, optional
        Path to the light curve extracted from the background aperture. The default is None.
    bkgratio : float, optional
        Area correction factor for the background light curve. It's mandatory if
        the 'lcbkg' argument is provided. The default is None.
    lcnet : str, optional
        Path to the net (background subtracted) light curve. The default is None.
    binsize : int, optional
        Desired temporal resolution in seconds. The default is 500 sec.

    Returns
    -------
    bool
        Return True.

    """
    _ownname = mylib.getownname()
    add_quotes = lambda s: "'" + s + "'"

    gtilimits = gti_get_limits(gtitable)
    nbin = (gtilimits[1] - gtilimits[0]) / binsize + 10  # Bins per interval

    lcurves = [add_quotes(lcraw.fspath)]  # list which will help to generate the argument
    # string for the lcurve task below

    # If lcbkg is needed, get bkgratio and calc new background light curve
    if lcbkg:
        if not bkgratio:
            raise TypeError("Argument 'lcbkg' requires 'bkgratio' but None is given.")
        mylib.printbold(f'BKGRATIO={bkgratio}')
        tmplcbkg = mylib.TempfilePath.generate_from('lcbkg.fts')
        if not (callftools("ftcalc '{}' '{}' RATE RATE*{:f}".format(lcbkg, tmplcbkg, bkgratio)) or
                callftools("ftcalc '{}' '{}' ERROR ERROR*{:f} clobber=yes ".format(tmplcbkg, tmplcbkg,
                                                                                   bkgratio))):
            mylib.printerr("Cannot apply correction for the BACKSCALE")
            raise ExternalTaskError('ftcalc', filename=lcbkg, caller=_ownname)
        lcurves.append(add_quotes(tmplcbkg.fspath))

    # Convert GTI to XRONOS format
    tmpgtixron = mylib.TempfilePath.generate_from('tmpgti.wi')
    if not callftools(f"gti2xronwin -i '{gtitable.name}' -o '{tmpgtixron}'"):
        mylib.printerr("Can't convert GTI to XRONOS format")
        raise ExternalTaskError('gti2xronwin', filename=lcbkg, caller=_ownname)

    # Prepare PCO file
    tmppco = mylib.TempfilePath.generate_from('tmppco.pco')
    with open(tmppco, 'w') as fpco:
        fpco.write("plot off\n")
        fpco.write("Col 1 on 2\n")
        if lcbkg: fpco.write("Col 2 on 3\n")
        if lcnet:
            fpco.write("Col 4 on 4\n")
            lcurves.append(add_quotes(lcnet.fspath))
        fpco.write("Lab T {}\n".format(outepsimg.name))
        fpco.write("Lab Y Rate, counts/sec\n")
        fpco.write("R\n")
        fpco.write(f"H {outepsimg}/CPS\n")
    fpco.close()

    cmd = "lcurve {:d} {} window='{}' dtnb={} nbint={} outfile=' ' plot=yes " \
          "plotdev='/dev/null/PS' plotdnum=3 plotfile='{}' </dev/null " \
          ">/dev/null".format(len(lcurves), ' '.join(lcurves), tmpgtixron, binsize,
                              int(nbin), tmppco)
    if not callftools(cmd):
        mylib.printerr("Can't plot the light curve")
        raise ExternalTaskError('lcurve', filename=lcraw, caller=_ownname)

    _check_result_file_appeared(outepsimg, _ownname)
    return True


def grppha_group_spectrum(infile: FilePathAbs, outfile='', group_min_counts=1, clobber=False):
    """Call grppha to group the spectrum."""
    _ownname = mylib.getownname()

    if not outfile:
        outfile = infile.with_stem_ending('_grp{}'.format(group_min_counts))

    check_file_not_exist_or_remove(outfile, override=clobber, action='exception')
    if not callftools('grppha {} {}'.format(infile, outfile),
                      stdin='group min {}\nexit\n'.format(group_min_counts)):
        mylib.printerr("Can't plot the light curve")
        raise ExternalTaskError('lcurve', filename=infile, caller=_ownname)

    _check_result_file_appeared(outfile, _ownname)
    return True
