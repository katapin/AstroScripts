"""Load light curve from FITS files."""

import numpy as np
from mypythonlib import Actions, printwarn
from ..other import fits_keywords_getmany
# from ._classes import LCrateBinnedEven
from astropy.io import ascii
from os import PathLike
from .core import TUnits, LCType, LCurve
from astropy.time import Time as apyTime

#TODO list
# 1) Steal the astropy.io.ascii.load() implemetation of the loading factory
# 2) Switch from custom warnings to Python stdlib warnings

_registered_loaders={}


def register(name: str):
    """Register class as a lc-loader."""
    def decorator(cls: type):
        _registered_loaders[name] = cls
        return cls
    return decorator

class BaseLoader:
    def me(self):
        print('Im', self.__class__)

@register('XMM')
class XmmLoader(BaseLoader):
    pass

@register('SWIFT')
class SwiftLoader(BaseLoader):
    pass



class _LC_read_helper():
    """Provide telescope dependent column mapping for LCurveBinnedEven.from_fits()."""

    keyword_to_read = ['OBS_ID', 'EXPOSURE', 'LIVETIME', 'DATE-OBS', 'DATE-END',
                       'OBJECT', 'MJD-OBS', 'TOTCTS', 'CHANMIN', 'CHANMAX', 'BACKSCAL']

    def __init__(self, HDU, absent_keywords={}):
        telescope = HDU.header.get('TELESCOP', 'generic')
        self.telesope = telescope
        self.HDU = HDU
        self.absent_keywords = absent_keywords
        processor = self.__getattribute__(telescope.capitalize())  # Get specific functions
        processor()

    def _do_job(self, colmap: dict, keyword_to_read: list = []):
        # colmap = dict(time='TIME', rate='COUNT_RATE', etc.), it will be used to
        # unpack arguments for LCurveBinnedEven constructor as arg=vector
        self.required_columns = {arg: self.HDU.data[colname] for arg, colname in colmap.items()}
        self.extra_columns = {colname: self.HDU.data[colname] for colname in self.HDU.columns.names \
                              if colname not in colmap.values()}

        # Read time keywords (mandatory) from FITS header
        # self.timekeys = (fits_keywords_getmany(self.HDU, LCrateBinnedEven._timekeys_default,
        #                                             self.absent_keywords, action=Actions.EXCEPTION))
        self.TIMEDEL = self.timekeys['TIMEDEL']

        # Read optiocal keywords from FITS header
        self.keywords = {'TELESCOP': self.telesope}
        self.keywords.update(fits_keywords_getmany(self.HDU, keyword_to_read,
                                                        action=Actions.NOTHING))

    def Generic(self):
        colmap = dict(time='TIME', rate='RATE', rate_err='ERROR', fracexp='FRACEXP')
        self._do_job(colmap)

    def Swift(self):
        colmap = dict(time='TIME', rate='RATE', rate_err='ERROR', fracexp='FRACEXP')
        self._do_job(colmap, self.keyword_to_read)

    def Chandra(self):
        colmap = dict(time='TIME', rate='COUNT_RATE', rate_err='COUNT_RATE_ERR')
        self.absent_keywords['TIMEZERO'] = 0
        self._do_job(colmap, self.keyword_to_read)
        # Chandra doesn't have the FRACEXP column but has EXPOSURE instead
        self.required_columns['fracexp'] = self.HDU.data['EXPOSURE'] / self.TIMEDEL
        # TODO make time_err from TIME_END - TIME_START

    def Xmm(self):
        self.absent_keywords.update({'TIMEZERO': 0, 'TIMEPIXR': 0.5})
        colmap = dict(time='TIME', rate='RATE', rate_err='ERROR')
        if 'FRACEXP' not in self.HDU.columns.names:
            printwarn("You are trying to use not correct XMM-Newton light curve "
                         "which is not recomended. 'FRACEXP' column will be fake.")
            self._do_job(colmap, self.keyword_to_read)
            self.required_columns['fracexp'] = np.ones_like(self.required_columns['time'])
        else:
            self._do_job(colmap | {'fracexp': 'FRACEXP'}, self.keyword_to_read)


def _detect_columns(colnames: list[str]):
    # Find out whether the names are automatic (col1, col2, ...)
    res = {}
    _s = set(name[:3] for name in colnames)
    if len(_s) > 1:   # Yes, the names are not automatic
        colmap = {
            'time': ('TIME', 'MJD', 'JD'),
            'data': ('RATE', 'FLUX', 'COUNTS', 'MAG')
        }
        for field in colmap:
            # Build the list of candidate columns
            cand = [col for col in colnames if col.upper() in colmap[field]]
            if len(cand) > 1:
                raise TypeError(
                    "Automatic column detection yielded an ambiguous result for the "
                    f"field '{field}'. Please provide manual column mapping."
                )
            res[field] = cand[0]

            # Try to detect the column with errors
            l = len(cand[0])   # Length of the found column name
            err = [
                col for col in colnames if
                col.replace('_','').upper()[l:l+2] == 'ER'
            ]  # For example, RATE_ERR, RATE_ER or RATEER
            if len(err) == 1:
                res[field+'_err'] = err[0]
    else:   # Automatic names
        res['time'] = 'col1'
        if len(colnames) == 2:
            res['data'] = 'col2'
        elif len(colnames) == 3:
            res['data'] = 'col2'
            res['data_err'] = 'col3'
        elif len(colnames) > 3:
            res['time_err'] = 'col2'
            res['data'] = 'col3'
            res['data_err'] = 'col4'
    # len(res) == 1 is obviously wrong result, but we'll raise
    # an exception latter
    return res


def read_ascii(filepath: PathLike | str, *, column_mapping: dict[str, str] = None,
               timezero: apyTime = None, tunits: TUnits = None,
               lctype: LCType = None) -> LCurve:
    tbl = ascii.read(filepath)
    colmap = column_mapping or _detect_columns(tbl.colnames)
    if len(colmap) < 2:
        raise TypeError(
            f"Can't identify columns for '{filepath}'. Please provide "
            "manual column mapping."
        )

    coldict = {colname: tbl[colmap[colname]] for colname in colmap}

    # Treat time column
    if colmap['time'] in ('JD', 'MJD'):
        tunits = tunits or TUnits.day
        if colmap['time'] == 'MJD' and timezero is None:    # Time is absolute MJDs
            timezero = apyTime(tbl['MJD'][0], format='mjd')   # Take the timezero
            coldict['time'] -= tbl['MJD'][0]                  # And move the vector

    if tunits is None:
        raise TypeError("Can't guess 'tunits'. Please provide it manually.")

    # Treat data column
    try:
        lctype = lctype or LCType(colmap['data'].lower())
    except ValueError:
        raise TypeError("Can't guess 'lctype'. Please provide it manually.")

    return LCurve.from_columns(coldict, timezero=timezero, tunits=tunits,
                               lctype=lctype)

