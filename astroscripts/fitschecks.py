"""Provide check whether the FITS-file follows the standard."""

from typing import Type, Union, Callable
from astropy.io import fits
from mypythonlib import FilePath, Actions
from .extpath import ExtPath, ExtPathAbs

__all__ = [
    "file_is_fits",
    "file_is_lc",
    "file_is_gti",
    "hdu_is_lc",
    "hdu_is_gti",
    "hdu_is_spectrum"
]


def _helper_check_file_is(func_deep_check: Callable, filetype: str,
                          filepath: Union[str, FilePath, ExtPath],
                          default_hdu: Union[str, int], action: Actions, progname) -> ExtPathAbs | None:
    extpath = ExtPath(filepath).absolute()
    if extpath.extname.hdu is None: extpath.extname.hdu = default_hdu
    hdu = extpath.extname.hdu
    try:
        with fits.open(extpath) as fts:
            if isinstance(hdu, int):    # hdu argument is a number (index)
                if hdu >= len(fts):
                    action.do(f"'{extpath.name}' doesn't have an extension #{hdu}.",
                                      exception_class=IndexError, progname=progname)
                    return None
            elif isinstance(hdu, str):  # hdu argument is a dict keyword
                if hdu not in fts:
                    action.do(f"'{extpath.name}' doesn't have an extension "
                              f"named '{hdu}'.", exception_class=KeyError, progname=progname)
                    return None
            else:
                raise TypeError("Wrong type of the 'hdu' argument")

            hduobj=fts[hdu]
            _deepcheck_action = Actions.WARNING if action is not Actions.NOTHING else Actions.NOTHING
            if func_deep_check(hduobj, action=_deepcheck_action):
                return extpath
            else:
                hdustr = f'#{hdu}' if isinstance(hdu, int) else f"'{hdu}'"
                action.do(f"'{extpath.name}' doesn't seem to be a valid {filetype} "
                          f"(extension {hdustr}).", exception_class=TypeError, progname=progname)
    except FileNotFoundError:
        action.do(f"File '{extpath.name}' is not found.", exception_class=FileNotFoundError,
                  progname=progname)
    except OSError:
        action.do(f"'{extpath.name}' is not a FITS file.", exception_class=TypeError,
                  progname=progname)
    return None


def file_is_fits(filepath: Union[ExtPath, FilePath, str],
                 action: Union[Actions, str] = Actions.DIE, progname=None) -> ExtPathAbs | None:
    """Check if the file is readable FITS and return its absolute path."""
    return _helper_check_file_is(lambda x: True, 'FITS', filepath,
                                 action=Actions(action), default_hdu=0, progname=progname)


def file_is_lc(filepath: Union[ExtPath, FilePath, str],
               action: Union[Actions, str] = Actions.DIE, progname=None) -> ExtPathAbs | None:
    """Check whether the FITS-file is a light curve."""
    return _helper_check_file_is(hdu_is_lc, 'light curve', filepath,
                                 default_hdu=1, action=Actions(action), progname=progname)


def file_is_gti(filepath: Union[ExtPath, FilePath, str],
                action: Union[Actions, str] = Actions.DIE, progname=None) -> ExtPathAbs | None:
    """Check whether the FITS-file is a light curve."""
    return _helper_check_file_is(hdu_is_gti, 'GTI file', filepath,
                                 default_hdu=1, action=Actions(action), progname=progname)


def _helper_check_hdu_is(HDU, class1_key:str, type_string:str, classmax:int, action:Actions, **kwargs):
    if 'HDUCLASS' not in HDU.header:
        action.do("There is no 'HDUCLASS' keyword in the HDU header. "
                 "Probably this FITS file does not conform OGIP standard.", **kwargs)
        return False, {}
    if HDU.header['HDUCLASS'].upper() != 'OGIP':
        action.do("'HDUCLASS' keyword is not 'OGIP'.", **kwargs)
        return False, {}
    classkeys={}
    for classkey in ['HDUCLAS'+str(i) for i in range(1,classmax+1)]:
        if classkey not in HDU.header:
            action.do(f"There is no '{classkey}' keyword in the HDU header. "
                     "Probably this FITS file does not conform OGIP standard.", **kwargs)
            return False, {}
        classkeys[classkey] = HDU.header[classkey].upper()
    if HDU.header['HDUCLAS1'].upper() != class1_key:
        action.do("The 'HDUCLAS1' keyword is not '{}'. "\
                 "Probably this is not a {} extension.".format(class1_key, type_string), **kwargs)
        return False, classkeys
    return True, classkeys


def hdu_is_spectrum(HDU, withclasskeys: bool = False,
                    action: Union[Actions, str] = Actions.WARNING,
                    exception_class : Type[Exception] = TypeError, **kwargs):
    """Check if the HDU is a spectrum following the OGIP standard."""
    checkres, classkeys = _helper_check_hdu_is(HDU, 'SPECTRUM', 'spectral', 3, Actions(action),
                                               exception_class = exception_class, **kwargs)
    if not checkres:
        return (False, classkeys) if withclasskeys else False
    #TODO: other checks
    return (True, classkeys) if withclasskeys else True


def hdu_is_lc(HDU, withclasskeys=False,
              action: Union[Actions, str] = Actions.WARNING,
              exception_class: Type[Exception] = TypeError, **kwargs):
    """Check if the HDU is a light curve following the OGIP standard."""
    checkres, classkeys = _helper_check_hdu_is(HDU, 'LIGHTCURVE', 'light curve', 3, Actions(action),
                                               exception_class = exception_class, **kwargs)
    if not checkres:
        return (False, classkeys) if withclasskeys else False
    if 'TIME' not in HDU.columns.names:
        Actions(action).do("There is no column 'TIME'.", exception_class = exception_class, **kwargs)
        return (False, classkeys) if withclasskeys else False
    if 'RATE' not in HDU.columns.names:
        Actions(action).do("There is no column 'RATE'.", exception_class = exception_class, **kwargs)
        return (False, classkeys) if withclasskeys else False
    return (True, classkeys) if withclasskeys else True


def hdu_is_gti(HDU, withclasskeys=False,
               action: Union[Actions, str] = Actions.WARNING,
               exception_class: Type[Exception] = TypeError, **kwargs):
    """Check if the HDU is a GTI-file following the OGIP standard."""
    checkres, classkeys = _helper_check_hdu_is(HDU, 'GTI', 'GTI', 2, Actions(action),
                                               exception_class = exception_class, **kwargs)
    if not checkres:
        return (False, classkeys) if withclasskeys else False
    if 'START' not in HDU.columns.names:
        Actions(action).do("There is no column 'START'.", exception_class = exception_class, **kwargs)
        return (False, classkeys) if withclasskeys else False
    if 'STOP' not in HDU.columns.names:
        Actions(action).do("There is no column 'STOP'.", exception_class = exception_class, **kwargs)
        return (False, classkeys) if withclasskeys else False
    if len(HDU.data) == 0:
        Actions(action).do("File seems to have a valid GTI extension, but it doesn't "
                          "contain any data.", exception_class = exception_class, **kwargs)
        return (False, classkeys) if withclasskeys else False
    return (True, classkeys) if withclasskeys else True