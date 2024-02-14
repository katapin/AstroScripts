"""Provide facilities to manipulate light curves."""

import numpy as np
from numpy import ndarray
from numpy.typing import ArrayLike, NDArray
from astropy.io import fits
from astropy.time import (Time as apyTime,
                          TimeDelta as apyTimeDelta)
# from astropy import units
from typing import List, Union, Self
from enum import StrEnum

from functools import singledispatchmethod
from mypythonlib import printwarn
from mypythonlib.collections import ROdict
from .. import fitschecks as checks
from ..plot import PlotPair, PlotVector
from ..extpath import ExtPath
from ._common import LCError
from .._internal.columnbased import (
    _make_vector,
    ColumnBasedMinimal,
    ColumnBased
)

# from ._loader import _LC_read_helper

#TODO list
# 1) make time column an object of its own class:
#   - make a special object to store Timezero
#   - __eq__ must return True if both time columns data and timezeros coincides
#   - delegate it check_synchronicity()'s job

__all__ = [
    "LCrate",
    "LCrateBinned",
    "LCrateBinnedEven",
    "LCcntBinned",
    "LCcntBinnedEven"
]


# def _check_tunits(txt: str) -> str:
#     """Check whether the string represent valid type units."""
#     if not isinstance(txt, str):
#         raise TypeError('Time units must of the string type.')
#     if txt not in ('d', 's'):
#         raise ValueError(f"Unknown units: '{txt}'.")
#     return txt

class TUnits(StrEnum):
    """Allowed units for TimeVector."""
    day = 'd'
    sec = 's'


class LCType(StrEnum):
    """Types of light curves' data column."""
    rate = 'rate'
    counts = 'counts'
    mag = 'mag'


class TimeVector(ColumnBasedMinimal):
    """Class to store time data of time series.

    :param ticks:  numpy array or other iterable with time stamps
        representing time since 'timezero'. Must be sorted in
        ascending order.
    :param units:  units of time stamps. Either second (s) or days (d).
    :param timezero:  astropy.time.Time object or None.

    Since timezero stored in the astropy object, this class
    naturally treats times in different time systems.
    """

    def __init__(self, ticks: ArrayLike, units: TUnits,
                 timezero: apyTime = None):

        vec: ndarray = _make_vector(
            ticks, vecname='ticks', none_is_allowed=False, from_numbers=False
        )  # make a numpy arrray

        # Check the ticks are monotonically increase
        if not np.all(np.diff(vec) >= 0):
            raise ValueError("Time ticks must be only in ascending order.")

        # Check whether time units are days or seconds
        if TUnits(units) is TUnits.day:
            vec = vec*86400     # ticks are always stored in seconds

        super().__init__({'ticks': vec})
        self.timezero = timezero   # It will also call self._rebuild_hook()

    def _rebuild_hook(self):
        """To be called after slicing to normalize the 'ticks' vector."""
        # Shift ticks to zero
        ticks = self._columns['ticks']
        if self._timezero and ticks.size > 0:
            if ticks[0] != 0:
                dt = ticks[0]
                self._columns['ticks'] = ticks - dt
                self._timezero = self._timezero + apyTimeDelta(dt, format='sec')
        self._ticks = self._columns['ticks']     # make a shortcut

    def __eq__(self, other: Self) -> bool:
        if (self._timezero is not None) and (other._timezero is not None):
            tzcheck = self._timezero.isclose(other._timezero)
        elif self._timezero is other._timezero is None:
            tzcheck = True   # both None
        else:  # One is None but the other is not
            return False
        return super().__eq__(other) and tzcheck

    def __repr__(self):
        if self._timezero is not None:
            tz = f"'MJD{self._timezero.mjd}'"
        else:
            tz = "None"
        clsname = self.__class__.__name__
        return f"{clsname}(timezero={tz}, ticks={self._ticks!r}"

    def ticks_from(self, time: apyTime, units: TUnits):
        """Return timestamps referenced to 'time' in desired time units."""
        units = TUnits(units)  # Raise an error if the string is incorrect
        if time is self._timezero:  # 'None is None' is also here
            dt = 0
        else:
            if self._timezero is None:
                raise TypeError("The light curve doesn't the 'timezero' reference.")
            dt = (self._timezero - time).sec
        ticks = self._ticks + dt     # It's a new ndarray even if dt is 0
        if units is TUnits.sec:
            return ticks
        else:  # days
            return ticks/86400

    def seconds_from(self, time: apyTime):
        """Return timestamps in seconds referenced to 'time'."""
        return self.ticks_from(time, TUnits.sec)

    def days_from(self, time: apyTime):
        """Return timestamps in days referenced to 'time'."""
        return self.ticks_from(time, TUnits.day)

    def get_ticks(self, units: TUnits):
        """Return ticks in the desired units"""
        return self.ticks_from(self._timezero, units)

    @property
    def timezero(self):
        """Return timezero object."""
        # TODO make sure that astropy.Time objects are immutable
        return self._timezero

    @timezero.setter
    def timezero(self, val: apyTime):
        if val is not None and not isinstance(val, apyTime):
            raise TypeError(f"The 'timezero' argument must be a "
                            f"{type(apyTime).__name__} object or None.")
        self._timezero = val
        self._rebuild_hook()  # Shift ticks


    @property
    def tickss(self):
        """Return time stamps in seconds (copy)."""
        return self.get_ticks(TUnits.sec)

    @property
    def ticksd(self):
        """Return time stamps in days (copy)."""
        return self.get_ticks(TUnits.day)


class LCurve(ColumnBased):

    def __init__(self, timevec: TimeVector, data: ArrayLike, *, lctype: LCType,
                 tunits: TUnits, time_err: ArrayLike = None,
                 data_err: ArrayLike = None, binwidth: ArrayLike = None,
                 extra_columns: dict[str, ArrayLike] = None,
                 meta: dict = None):
        _loc = locals()

        self._lctype = LCType(lctype)
        self._tunits = TUnits(tunits)
        self._timevec = timevec

        reflen = len(timevec)
        _pd = {k: _make_vector(_loc[k], vecname=k, reflen=reflen) for
               k in ('time_err', 'data', 'data_err', 'binwidth')
               if _loc[k] is not None}  # Dict for protected columns

        if extra_columns is not None:
            _pd |= {k: _make_vector(extra_columns[k], vecname=k, reflen=reflen)
                    for k in extra_columns}

        super().__init__(columns=None, pcolumns=_pd)

    @classmethod
    def from_columns(cls, columns: dict[str, ArrayLike], *, timezero: apyTime,
                     lctype: LCType, tunits: TUnits) -> Self:
        """Create light curve object from dict of columns."""
        time = columns.pop('time')
        timevec = TimeVector(time, units=tunits, timezero=timezero)
        return cls(timevec, lctype=lctype, tunits=tunits, **columns)

    @property
    def lctype(self):
        return self._lctype

    @property
    def timevec(self):
        """Return time vector object."""
        return self._timevec

    def time_from(self, time: apyTime):
        """Return timestamps since 'time'."""
        return self._timevec.ticks_from(time, self._tunits)

    @property
    def time(self):
        """Return timestamps."""
        return self.time_from(self._timevec.timezero)

    @property
    def tunits(self):
        return self._tunits

    @tunits.setter
    def tunits(self, val: TUnits):
        if (enval := TUnits(val)) is self._tunits:
            return
        fact = {
            TUnits.sec: 86400,
            TUnits.day: 1/86400
        }  # factor
        self._tunits = enval
        for colname in ('time_err', 'binwidth'):
            if colname in self._columns:
                self._columns[colname] *= fact[enval]





class LCrate:
    """Base class for light curves.
    
    Light curve represented as a pair of columns: (TIME, RATE), optionally 
    errors for both columns can also be included for plotting purposes. Time 
    spacing doesn't matter. Default time properties are: TIMEUNITS='s', TIMESYS='UTC',
    TIMEZERO=0, MJDREF=0.
    """
    
    _colnames_main = ('time', 'time_err', 'rate', 'rate_err')  # essentual  columns
    _timekeys_default={'MJDREF':0, 'TIMEZERO':0, 'TIMESYS':'UTC', 'TIMEUNIT':'s'}
    _timekeys_allowed_values={'TIMEUNIT':['s','d'],
                              'TIMESYS':['UTC','TT']}
    
    _precision = 1e-5  # Accuracy for float comparison
    
    def __init__(self, time:ndarray, rate:ndarray, *, time_keywords:dict, 
                 keywords:dict={}, time_err:ndarray=None, rate_err:ndarray=None, 
                 extra_columns:dict[str, ndarray]={}):
        # if values_type not in ['RATE', 'COUNTS']:
        #     raise ValueError("Unknown type of the light curve. It can be "\
        #             "either 'RATE' or 'COUNTS'")

        _columns = {}
        refshape = (len(time),)  # reference shape
        
        # Copy columns. We don't want to update() the dictionary, but we want to compare
        # possible duplicates. It's for compatibility with subclasses.
        for colname,val in tuple(dict(__time=time, rate=rate, rate_err=rate_err, 
                    time_err=time_err).items()) + tuple(extra_columns.items()):
            if val is not None:
                if not isinstance(val, ndarray):  #This check is neaded for futher
                    #advanced sampling. Simple python lists or tuples are not appropreate here. 
                    raise TypeError("Data columns must be of the numpy array type")
                if val.shape != refshape: 
                    raise ValueError("All the arrays must be 1D of the the same length.")
                if colname in _columns:  #Compare dublicates
                    if np.any(_columns[colname] != val):
                        raise ValueError("The 'extra_columns' argument contains "\
                                 f"another '{colname}' array with different data.")
                    continue
                _columns[colname] = val
        self._columns = ColumnStorage(**_columns)
        
        # Check presence of time keywords
        if diff:= set(self._timekeys_default).difference(time_keywords):
            raise KeyError(f"Time keywords {diff} not found in 'time_keywords' dict")
        # if diff:=set(time_keywords).difference(self._timekeys_default):
        #     raise KeyError(f"Time keywords {diff} are unknown")
            
        _timekeys={} # Check validity of time keywords
        for key,val in time_keywords.items():
            if key in self._timekeys_allowed_values:
                if val not in self._timekeys_allowed_values[key]:
                    raise ValueError(f"{key} must be in {self._timekeys_allowed_values[key]}.")
            _timekeys['__'+key] = val  # Make all the timekeys readonly
            
        self.keywords = keywords.copy()   # Ordinary  python dict
        self._timekeys = ROdict(**_timekeys)
        _tk = self._timekeys.dict()  # Timekeys with normal names (without '__')
        self._mjdref = TimeVector(_tk['MJDREF'], format='mjd', scale=_tk['TIMESYS'].lower())
        self._timeunit = getattr(units, _tk['TIMEUNIT'])
        self._timezero = self._mjdref + _tk['TIMEZERO']*self._timeunit
   
    def __len__(self):
        """Return number of points in the light curve."""
        return len(self.columns['time'])
    
    def __str__(self):
        """Return string representation."""
        text=self._make_title() + '\n'
        text+='Time zero: MJD{} ({}) \n'.format(self._timezero.mjd, self._timezero.iso)
        text+='Columns: {}\n'.format(str(list(self._columns.keys())))
        text+=f'Length: {len(self)} bins\n'
        text+='==============================\n'
        cols = [x for x in self._colnames_main if x in self._columns]
        text += ' | '.join(x.capitalize() for x in cols) # Caption of the table
        indices = list(range(len(self)))
        prev = 0
        for i in set(indices[:10] + indices[-10:]):
            if (i-prev)>1:
                text+='\n'+'\t'.join(['---']*len(cols))
            text+='\n'+'\t'.join(str(self._columns._dict[x][i]) for x in cols)
            prev = i
        return text
            

    def get_column(self, colname: str):
        """Return column by its name."""
        if colname in ('_time_err', '_rate_err'):
            return self._columns.get(colname[1:])
        elif colname[:9] == '_timeabs_':
            timeobjlst = [self._timezero+x*self._timeunit for x in self._columns['time']]
            if colname == '_timeabs_astropy':
                return np.array(timeobjlst)
            elif colname == '_timeabs_mjd':
                return np.array([x.mjd for x in timeobjlst])
            elif colname == '_timeabs_datetime':
                return np.array([x.datetime for x in timeobjlst])
        elif colname == '_timeabserr_mjd':
            if self._timekeys['TIMEUNIT'] == 'd':
                return self._columns.get('time_err')
            else:  #time is in seconds
                return self._columns.get('time_err')/86400
        return self._columns[colname]
    
    #### Properties
    @property
    def timeabs_astropy(self):
        """Return timestamps as astropy.time.Time objects."""
        return self.get_column('_timeabs_astropy')
    
    @property
    def timeabs_datetime(self):
        """Return timestamps as python datetime objects."""
        return self.get_column('_timeabs_datetime')
    
    @property
    def timeabs_mjd(self):
        """Return timestamps as MJDs."""
        return self.get_column('_timeabs_mjd')
    
    @property 
    def time(self):
        """Return TIME column of the light curve."""
        return self._columns['time']
    
    @property 
    def rate(self):
        """Return RATE column of the light curve."""
        return self._columns['rate']
    
    @property 
    def rate_err(self):
        """Return ERROR column of the light curve if it exists or zero vector."""
        return self.get_column('_rate_err')
    
    @property 
    def time_err(self):
        """Return error of the time column or zero vector."""
        return self.get_column('_time_err')
    
    @property 
    def columns(self):
        """Return read-only dict with columns passed to constructor."""
        return self._columns
    
    @property 
    def colnames(self):
        """Return list of names of existing columns."""
        return list(self._columns.keys())
    
    @property 
    def time_keywords(self):
        """Return read-only dict with time keywords."""
        return self._timekeys
    
    @property 
    def timezero(self):
        """Return astropy.time.Time object corresponding to the first timestamp."""
        return self._timezero
    
    @property 
    def timeunits(self):
        """Return astropy.units object corresponind to TIMEUNIT keyword."""
        return self._timeunit
    
    @property 
    def extra_columns(self) -> dict:
        """Return dict with 'extra_columns' argument passed to the constructor."""
        maincols, extracols = self._split_columns(self._columns.dict())
        return extracols
        
    @classmethod 
    def get_time_defaults(cls):
        """Return dict with mandatory time keywords and their defaults."""
        return cls._timekeys_default.copy()
    
    #### Copying and sampling
    @classmethod 
    def _split_columns(cls, allcols:dict)-> tuple[dict, dict] :
        """Split column to 'main' and 'extra' for specific subclass.""" 
        maincols={}  #Main columns
        allcols = allcols.copy()
        for key in tuple(allcols):  #Without tuple dict.pop() is not allowed
            if key in cls._colnames_main:   
                maincols[key] = allcols.pop(key)
        return maincols, allcols
    
    @singledispatchmethod
    def __getitem__(self, arg):
        """Get subset of the light curve."""
        raise TypeError("Unsupported type of key: {}.".format(type(arg)))
    
    @__getitem__.register    
    def _(self, index: int):  # Get one point
        #TODO
        raise NotImplementedError('__getitem__ for int')
    
    @__getitem__.register    
    def _(self, arg: Union[slice, list, ndarray]):
        obj = self.copy()
        for key in self._columns:
            obj._columns['_'+key] = self._columns['_'+key][arg]
        return obj
        
    def copy(self):
        """Clone the object."""
        cls = self.__class__
        obj = cls.__new__(cls)
        for key, val in self.__dict__.items():
            if hasattr(val, 'copy'):
                obj.__dict__[key] = val.copy()
            else:
                obj.__dict__[key] = val
        return obj
        
    #### Make PlotPair
    def _make_title(self):
        keys=self.keywords
        empty = 'Unknown'
        text = "Light curve: " + keys.get('FILENAME','').name + '\n'
        text += 'Object: {}, '.format(keys.get('OBJECT',empty)) 
        text += '{}-'.format(keys.get('TELESCOP',empty)) 
        text += '{}'.format(keys.get('OBS_ID',empty)) 
        return text
    
    def _prepare_plot_meta(self, abstime:bool, dates:bool):
        """Return labels and real columns for time and rate."""
        units={'s':'seconds','d':'days'}
        tzero=self._timezero
        colmap = {'rate':'rate'}
        colmap['rate_err']='rate_err' if 'rate_err' in self._columns else None
        text='Time'
        if abstime:
            if not dates:   #MJD
                text+=', MJD'
                colmap['time'] = '_timeabs_mjd'
                colmap['time_err'] = '_timeabserr_mjd' if 'time_err' in \
                    self._columns else None
            else:
                colmap['time'] = '_timeabs_datetime'
                colmap['time_err'] = None
        else:
            text+=', {} since '.format(units[self._timekeys['TIMEUNIT']])
            text+= tzero.iso if dates else f'MJD{tzero.mjd:.5f}' 
            colmap['time'] = 'time'
            colmap['time_err'] = 'time_err' if 'time_err' in self._columns else None
        labels={'time': text, 'rate':'Rate, counts/s'}
        return colmap, labels 
    
    def make_PlotPair(self, *, X:str=None, Xerr:str=None, Y:str=None, Yerr:str=None,
                    abstime:bool=True, dates:bool=False) ->PlotPair:
        """Return PlotPair object containing desired columns.
        
        If no column names passed, it generate a 'Rate vs Time' plot. Any 
        combinations of avalible columns can be passed. By default,
        Time is just the time vector stored in the light curve's table. One can
        convert it to absolute time (MJD or dates) using special options.

        Parameters
        ----------
        X : str, optional
            Column for X axis. If None the defaults will be used.
        Xerr : str, optional
            Column for error of X axis. If None the defaults will be used
        Y : str, optional
            Column for Y axis. If None the defaults will be used
        Yerr : str, optional
            Column for error of Y axis. If None the defaults will be used
        abstime : Bool, optional
            Use absolute time for time columns. The default is True.
        dates : TYPE, optional
            Use dates insted of MJDs. The default is False.

        Returns PlotPair object
        """
        loc = locals()
        #Column defaults
        defcols = dict(X='time', Xerr='time_err', Y='rate', Yerr='rate_err') 
        realcols, labels = self._prepare_plot_meta(abstime, dates)
        
        #Plotting options for PlotPair constructor
        pltopts={'title':self._make_title(),'fmt':'+'} 
        ppargs={'pltopts':pltopts}  #Arguments for PlotPair constructor
        for key in ('X', 'Y'):
            keyerr = key+'err'
            colname, colerrname = loc[key], loc[keyerr]  
            if colname is None:
                colname = defcols[key]  #Take default if None,  e.g. 'time'
                if colerrname is None:  #Take default error if both col end colerr are None
                    colerrname = defcols[keyerr]
                    
            pltopts[key.lower()+'label']= labels.get(colname,colname)
            
            colname = realcols.get(colname, colname) #Replace, for example, 'time' -> '_timeabs_mjd'
            colerrname = realcols.get(colerrname, colerrname)
            if loc[keyerr] is not None and colerrname is None:
                printwarn("The option '{keyerr}={loc[keyerr]}' will be ignored.")

            if colerrname:  #Empty '' may be used to omit existing errors
                ppargs[key] = PlotVector.create_with_errors_symmetric(self.get_column(colname), 
                                    self.get_column(colerrname))
            else:
                ppargs[key] = PlotVector(self.get_column(colname))
            
        colors={'TOTAL':'k', 'NET':'b', 'BKG':'r'}
        if hduclass:=self.keywords.get('HDUCLAS2',None):
            pltopts['color']=colors[hduclass]

        return PlotPair(**ppargs)
    
    ####  BKG curve
    
    # def _check_curve_type(self, other):
    #     if not isinstance(other, self.__class__):
    #         raise TypeError("The second light curve must have the same type as "\
    #                         "the first one.")
                
    def get_bkgratio(self, bkgcurve):
        """Return ratio of BKG to SRC collection areas."""
        # self._check_curve_type(bkgcurve)
        if 'BACKSCAL' in self._keywords and 'BACKSCAL' in bkgcurve._keywords:
            return self._bkgcurve['BACKSCAL']/self._keywords['BACKSCAL']
        else:
            None
            
    def check_synchronicity(self, other, silent=False):
        """Check wether two light curves are synchronous."""
        #TODO check time+TIMZERO first and than check timevectors
        # self._check_curve_type(other)
        printfun =  printwarn if not silent else lambda x: None
        if len(self) != len(other):
            printfun("Light curves have different lengths")
            return False
        if np.any(np.abs(self.time - other.time) > self._precision):
            printfun("Light curves have different time vectors.")
            return False
        if self._timekeys != other._timekeys:
            printfun("Light curves seems to be synchronicity but have different "\
                     "time metadata.")
            return False
        return True
            
    
    def make_gti(self):
        pass
                    
                
class LCrateBinned(LCrate):
    """Light curve with known widths of time bins.
    
    The timestamps don't have to be evenly spaced, and bin widths can be
    any. It's assumed that all the timestamps refer to the middle of their time
    bins. The TIMEPIXR keyword is ignored if present.
    """
    
    _colnames_main = LCrate._colnames_main + tuple(['binwidth'])
    
    
    class Groups(list):
        """Allow to save metadata in a list of groups."""
        
        def __init__(self, lst, *, lcurve):
            super().__init__(lst)
            self.lcurve = lcurve
    
    def __init__(self, time:ndarray, rate:ndarray, binwidth:ndarray, *,
                   time_keywords:dict, keywords:dict={}, time_err:ndarray=None, 
                   rate_err:ndarray=None, extra_columns:dict[str, ndarray]={}):
        super().__init__(time, rate, time_err=time_err, rate_err=rate_err, 
                 time_keywords=time_keywords, keywords=keywords,
                 extra_columns={'binwidth':binwidth} | extra_columns)
        
    @property 
    def binwidth(self):
        """Return vector with real exposure time of each bin."""
        return self._columns['binwidth']
    
    def reallen(self):
        """Return number of points with nonzero exposure."""
        return np.count_nonzero(self._columns['binwidth'])
    
    def remove_empty(self):
        """Return a new light curve with bins with only nonzero exposure."""
        mask = (self._columns['binwidth'] > 0.0)
        return self[mask]

    def toCounts(self, precision:float=None):
        """Return LCurveCounts of rates can be converted in ints with desired precision."""
        precounts = self._columns._dict['rate']*self._columns._dict['binwidth']
        counts = np.round(precounts)
        _precision = precision or self._precision
        if np.any(np.abs(counts-precounts) >= _precision):
            badbins = np.nonzero(np.abs(counts-precounts) >= _precision)[0]
            raise LCError("Can't convert light curve from rates to counts. At least "\
                      f"{len(badbins)} of {len(self)} bins are far from integers.")
                
        cols = self._columns.dict()
        #Remove old rate columns because they may contradict to counts/TIMEDEL
        #because of rounding errors.
        del cols['rate']  
        if 'rate_err' in cols: del cols['rate_err']
        
        maincols, extracols = self._class_cnt._split_columns(cols)
        return self._class_cnt(counts=counts.astype(int), time_keywords=self.time_keywords, 
                    keywords=self.keywords, extra_columns=extracols, **maincols)
    
class LCrateBinnedEven(LCrateBinned):
    """Evenly binned light curve with all reqired OGIP keywords."""
    
    _colnames_main = LCrate._colnames_main + tuple(['fracexp'])
    _timekeys_default = LCrateBinned._timekeys_default | {'TIMEPIXR':0.5,
             'TIMEDEL':1}
    
    def __init__(self, time:ndarray, rate:ndarray, fracexp:ndarray, *,
                   time_keywords:dict, keywords:dict={}, time_err:ndarray=None, 
                   rate_err:ndarray=None, extra_columns:dict[str, ndarray]={}):
        
        #Check step of the timestamps
        if 'TIMEDEL' not in time_keywords:
            raise KeyError("Required time keyword 'TIMEDEL' is not found.")
        timedel = time_keywords['TIMEDEL']
        if np.abs(timedel - np.min(time[1:]-time[0:-1])) > self._precision:
            raise ValueError("It seems the time vectror is not evenly spaced or "\
                 "passed 'TIMEDEL' values is wrong.")
        
        binwidth = fracexp*timedel
        super().__init__(time, rate, binwidth, time_keywords=time_keywords, 
                 keywords=keywords, time_err=time_err, rate_err=rate_err, 
                 extra_columns={'fracexp':fracexp} | extra_columns)
        
        # if time_keywords['TIMEPIXR']  != 0.5:  #TODO Move to rebin?
        #     raise NotImplementedError('Only TIMEPIXR=0.5 is supported yet.')
            
    @classmethod
    def from_fits(cls, path: ExtPath, absent_keywords={}):
        """Load an evenly binned light curve from a FITS file."""
        with fits.open(path.fspath) as fts:
            HDU=fts[path.hdu]
            res, classkeys = checks.hdu_is_lc(HDU, withclasskeys=True) ##TODO action=exceptions
            if classkeys['HDUCLAS3'] == 'COUNTS':
                raise NotImplementedError('Only light curves in rates are '
                          'supported (HDUCLAS3=RATE).')
            data = _LC_read_helper(HDU)  #Parse data from the specific telescope
        data.keywords.update(classkeys)
        data.keywords['FILENAME'] = path
        return cls(**data.required_columns, time_keywords=data.timekeys, 
                       keywords=data.keywords, extra_columns=data.extra_columns)
                            
    @property 
    def fracexp(self):
        """Return FRACEXP column."""
        return self._columns['fracexp'] 
    
    @property
    def TIMEDEL(self):
        """Return temporal resolution of the light curve."""
        return self._timekeys['TIMEDEL']
    
    def getgroups_continuous(self, maxgap:int, minlen:int) -> List[List[int]]:
        """Return groups of points corresponding to continuous intervals.

        Parameters
        ----------
        maxgap : float
            Maximum of absent bins between to subgroups to consider them as still belonging 
            to the single continues group.
        minlen : float
            Minimum number of bins in the group. Shorter groups will be dropped.

        Returns
        -------
        groups : TYPE
            DESCRIPTION.

        """
        #This grouping method assumed constant time stap and requires TIMEDEL keyword 
        lc = self.remove_empty()
        groups=[]
        time = lc._columns._dict['time']
        diff = (time[1:] - time[:-1])/lc.TIMEDEL
        indices = np.nonzero(diff > maxgap+1)[0] #indexed of gaps wider than allowed
        prevind = -1 #previous index
        for curind in list(indices) + [len(time)-1]:
            if (curind - prevind) >= minlen:
                groups.append(list(range(prevind+1,curind+1)))
            prevind = curind
        
        return self.Groups(groups, lcurve=self)
    
    def correct_timepix(self):
        # timekeys['TIMEZERO'] += (0.5 - timekeys['TIMEPIXR'])*self.TIMEDEL
        pass
       
    def rebin(self, groups): 
        """Rebin the light curve using list of bin groups."""
        if self.time_keywords['TIMEPIXR'] != 0.5:
            raise TypeError("Light curves with TIMEPIXR != 0.5 are not "
                    "supported. Use correct_timepix method first.")
            
        lc = super().rebin(groups)
        del lc.time_keywords['_TIMEDEL'], lc.time_keywords['_TIMEPIXR']
        return lc
    
        
class _LCcntMixin():
    """Mixin class to add counts columns and associated facilities."""
        
    @property 
    def counts(self):
        """Return vector with real exposure time of each bin."""
        return self._columns['counts'] 
    
    def rebin(self, groups):
        """Rebin the light curve using list of bin groups."""
        return super().rebin(groups).toCounts()
    
    def getgroups_mincounts(self, mincounts:int, maxbinwidth:float):
        
        lc = self.remove_empty()  #remove bins with zero exposure
        time=lc._columns._dict['time']
        counts=lc._columns._dict['counts']
        bw=lc._columns._dict['binwidth']
        
        #Add fake [hoton] at the end of the light curve
        time = np.append(time, time[-1]+10*maxbinwidth)
        counts = np.append(counts, 1)       #put one photon
        bw = np.append(bw, bw[-1])     #just copy the last binwodth
        
        
        groups, collected = [], [] #lists for groups and acculumated photons
        buff=0   #photon buffer
        gs=0     #first element of current group
        prev=-1  #Index of previous found photon
        for k in np.nonzero(counts>0)[0]:   #loop over bins with photons
            if buff >= mincounts:    #Buffer is already full
                ge=prev    #Set end of the group to positon of previoues photon
                groups.append(list(range(gs,ge+1)))  #grab bins (icludnig ge itself)
                collected.append(buff)   #Number of acollected photons
                buff=0     #Clean the buffer 
                gs=ge+1    #Next group will start from bin next to last photon 
            while (time[k]-time[gs]+0.5*(bw[k]+bw[gs])) > maxbinwidth: #Here we assumed 
                #that timestamps refer to bins' centers
                for j in range(k-1, gs-1, -1):  #reverse loop to find ge
                    if (time[j]-time[gs]+0.5*(bw[j]+bw[gs])) <= maxbinwidth:
                        ge=j
                        groups.append(list(range(gs,ge+1)))
                        collected.append(buff)
                        buff=0
                        gs=ge+1
                        break
            buff+= counts[k]
            prev = k  #Save index of the added photon
            
        grpobj = self.Groups(groups, lcurve=self)
        grpobj.counts  = collected
        return grpobj
        
    
class LCcntBinned(_LCcntMixin, LCrateBinned): 
    
    _colnames_main = ('time', 'time_err', 'counts', 'binwidth')
    
    def __init__(self, time:ndarray, counts:ndarray, binwidth:ndarray, *,
                   time_keywords:dict, keywords:dict={}, time_err:ndarray=None, 
                   extra_columns:dict[str, ndarray]={}):
        
        super().__init__(time, counts/binwidth, binwidth, time_keywords=time_keywords, 
                keywords=keywords, time_err=time_err, rate_err=np.sqrt(counts)/binwidth,
                extra_columns={'counts':counts} | extra_columns)
        
        
class LCcntBinnedEven(_LCcntMixin, LCrateBinnedEven):       
    
    _colnames_main = ('time', 'time_err', 'counts', 'fracexp')
    
    def __init__(self, time:ndarray, counts:ndarray, fracexp:ndarray, *,
                   time_keywords:dict, keywords:dict={}, time_err:ndarray=None, 
                   extra_columns:dict[str, ndarray]={}):
        
        #We don't have binwidth column now, but parent init will produce it.
        #Let's call it with zero rate and rate_err
        super().__init__(time, np.zeros_like(time), fracexp, time_keywords=time_keywords, 
                keywords=keywords, time_err=time_err, rate_err=np.zeros_like(time),
                extra_columns={'counts':counts} | extra_columns)
        
        #Now let's add appropriate rate and rate_err columns
        self._columns['rate'] = counts/self.binwidth
        self._columns['rate_err'] = np.sqrt(counts)/self.binwidth
    

# Set corresponding LCcnt-classes for LCrates
LCrateBinned._class_cnt = LCcntBinned    
LCrateBinnedEven._class_cnt = LCcntBinnedEven








