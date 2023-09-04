"""Provide facilities to manipulate light curves."""

import numpy as np
from numpy import ndarray
import scipy
from astropy.io import fits
from astropy.time import Time
from astropy import units
from typing import List, Union
from functools import singledispatchmethod
import mypython as my
from mypython import Actions, _to_abspath
from .plot import PlotPair, PlotVector, PlotPairUplimits
from . import main 
from .main import ExtPath


#TODO list
#1)Class from xmm light curve that can excrete the back ground LC from
# BACKV and BACKE columns in NET lc-FITS.


#### Some thoughts about rebining algorithms
#
#About TIMEPIXR
#Let's suppose timestep is 2s. If TIMEPIXR=0.5 the timestamps will's as follows
#  | bin-1 | bin-2 | bin-3 | bin-4 | bin-5 | bin-6 | bin-7 | bin-8 |....   
#  |   ^   |   ^   |   ^   |   ^   |   ^   |   ^   |   ^   |   ^   |
#      0s      2s      4s      6s  |   8s  |  10s  |  12s  |  14s  |  
#The TIMEZERO keyword in the FITS  header refers to the timstamp 0s (and not to
#the star time of the exposure!). If we, for example, everage 4 bins, we'll obtain:
#  |           long bin            |           long bin            |.... 
#  |   ^           ^               |               ^               |           
#    timezero      3s              |              11s              |
#Now the TIMEZERO no longer refers to the bin center, but timestamps are 
#correct and still refer the bin center. Thus, only TIMEPIXR=0.5 remain valid
#after the rebinning
#
#
############################
#About getgroups_continous()
#Let's suppose maxgap = 1 (in units of TIMEDEL), one bin can be missed 
# 
# detector in operation: |+++|+++|+++|   |+++|+++|   |   |   |+++|+++|   |   |+++|+++|+++|   |   |+++|| END
#    time              : | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10| 11| 12| 13| 14| 15| 16| 17| 18|| END
#    indixes           : | 0 | 1 | 2 |   | 3 | 4 |           | 5 | 6 |       | 7 | 8 | 9 |       | 10||
#
#Data (after removing empty bins)
# array index  : 0 1 2 3 4 5  6  7  8  9  10
# time column  : 0 1 2 4 5 9  10 13 14 15 18
# time[1:]     : 1 2 4 5 9 10 13 14 15 18
# time[:-1]    : 0 1 2 4 5 9  10 13 14 15
# time diff    : 1 1 2 1 4 1  3  1  1  3
#Indices of diff>(maxgap+1)  : 4 6 9
#group0 = range(0, 4+1)  - > 0 1 2 3 4
#group1 = range(4+1, 6+1) -> 5 6
#group2 = range(6+1, 9+1) -> 7 8 9
#group3 = range(9+1, len) -> 10
#
#
#############################
#About getgroups_mincounts()
#
# detector in operation: |+++|+++|+++|+++|+++|+++|+++|+++|+++|   |+++|+++|+++|+++|+++|   |+++|+++|   |+++|+++|+++|+++|+++||END
# photons arrveved     : |   | 1 |   | 2 | 2 |   | 2 |   |   |   |   |   |   |   |   |   |   | 1 |   | 1 | 1 |   |   |   ||END 
#    time              : | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10| 11| 12| 13| 14| 15| 16| 17| 18| 19| 20| 21| 22| 23||END
#    indixes           : | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |   | 9 | 10| 11| 12| 13|   | 14| 15|   | 16| 17| 18| 19| 20||END
# Let's add a photon far-far away (e.g. at time=1000) for algorithmic purposes
# Cells with arrived counts > 0 
# indices of photons : 1 3 4 6 15 16 17 21
# timestamps         : 1 3 4 6 17 19 20 1000
# Case 1) mincounts = 4, maxbinlength=5s  
# Starting...
# Set gs=0 -- index of group start, prev = -1 --- index of previous photon
# Begin k-loop over list with photon indixes 
# prev=-1, k=1, put photon to 'buff', because 'buff'<4 and time[k]-time[gs]+bw < maxbinlength
# prev=1, k=3, put photons 
# prev=3, k=4, put photons bacause before if 'buff' was still <4
# prev=4, k=6, buf=5. Set ge=prev, group.append(range(gs,ge+1)), set buff=0,
#   set gs=ge+1, put ptohotos to buf 
# prev=6, k=15, buf=2, time[k]-time[gs]+bw> maxbinlength
#   Begin loop while (time[gs] - time[k]+bw) > maxbinlength
#     Begin reverse j-loop from k-1(=14) to gs(=5)
#     j=14, time[14]=16, time[gs]=5, (14-5+bw)>5
#     ....
#     j=8, time[8]=8,  time[gs]=5, (10-8+bw)=4
#         set ge=j, group2.append(range(gs,ge+1)), set buf=0, set gs=ge+1=9,
#         j-loop ended
#     Begin reverse j-loop from k-1=(14) to gs(=9)
#     j=14, time[14]=16, time[gs]=10, (16-10+bw)>5
#     j=13, time[14]=13, time[gs]=10, (13-10+bw)=5
#         set ge=j, group2.append(range(gs,ge+1)), set buf=0, set gs=ge+1=14,
#         j-loop ended
#     while loop ended
#   put photon at k=15 to buf, buf=1 
# prev=15, k=16, put photon
# prev=16, k=17, put photon
# prev=17, k=21, its the fake photon with time[21]=1000
#   Begin loop while (time[gs] - time[k]+bw) > maxbinlength
#     Begin reverse j-loop from k-1(=20) to gs(=14)
#     j=20, time[23]=24, time[gs]=16
#     ....
#     j=17

class LCError(Exception):
    """Class for errors arising during manipulation with light curves.""" 
    
    pass


class ColumnStorage(my.ROdict):
    """Read-only dict to store ndarray vectors."""
    
    def _item_change_existing_as_normal(self, key, val):
        self._dict[key][:] =  val
        
    def _item_get_missing_as_normal(self, key):
        raise KeyError(f"Column '{key}' doesn't exist.")

    def _item_del_missing_as_normal(self, key):
        raise KeyError(f"Column '{key}' doesn't exist.")

    def __str__(self):
        """Return string representation."""
        text = self.__class__.__name__ + ': {\n'
        lst=[]
        for k,v in self._dict.items():
            lst.append( (f"'{k}'*:" if k in self._readonly else f"'{k}':") +\
                      repr(v))
        return text + '\n'.join(lst) +' }'


class _LC_read_helper():
    """Provide telescope dependent column mapping for LCurveBinnedEven.from_fits()."""
    
    keyword_to_read = ['OBS_ID', 'EXPOSURE', 'LIVETIME', 'DATE-OBS', 'DATE-END',
                       'OBJECT', 'MJD-OBS', 'TOTCTS', 'CHANMIN', 'CHANMAX', 'BACKSCAL']
    
    def __init__(self, HDU, absent_keywords={}):
        telescope = HDU.header.get('TELESCOP','generic')
        self.telesope = telescope
        self.HDU = HDU
        self.absent_keywords = absent_keywords
        processor = self.__getattribute__(telescope.capitalize()) #Get specific functions
        processor()
        
    def _do_job(self, colmap:dict, keyword_to_read:list=[]):
        #colmap = dict(time='TIME', rate='COUNT_RATE', etc.), it will be used to
        #unpack arguments for LCurveBinnedEven constructor as arg=vector
        self.required_columns = {arg:self.HDU.data[colname] for arg,colname in colmap.items()}
        self.extra_columns = {colname:self.HDU.data[colname] for colname in self.HDU.columns.names \
             if colname not in colmap.values()}
            
        #Read time keywords (mandatory) from FITS header 
        self.timekeys=(main.fits_keywords_getmany(self.HDU, LCrateBinnedEven._timekeys_default, 
               self.absent_keywords, action=Actions.EXCEPTION))
        self.TIMEDEL = self.timekeys['TIMEDEL']
        
        #Read optiocal keywords from FITS header 
        self.keywords = {'TELESCOP':self.telesope}
        self.keywords.update(main.fits_keywords_getmany(self.HDU, keyword_to_read, 
               action=Actions.NOTHING))
        
    def Generic(self):
        colmap = dict(time='TIME', rate='RATE', rate_err='ERROR', fracexp='FRACEXP')
        self._do_job(colmap)
        
    def Swift(self):
        colmap = dict(time='TIME', rate='RATE', rate_err='ERROR', fracexp='FRACEXP')
        self._do_job(colmap, self.keyword_to_read)
        
    def Chandra(self):
        colmap = dict(time='TIME', rate='COUNT_RATE', rate_err='COUNT_RATE_ERR')
        self.absent_keywords['TIMEZERO']=0
        self._do_job(colmap, self.keyword_to_read)
        #Chandra doesn't have the FRACEXP column but has EXPOSURE instead
        self.required_columns['fracexp'] = self.HDU.data['EXPOSURE']/self.TIMEDEL 
        #TODO make time_err from TIME_END - TIME_START
        
    def Xmm(self):
        self.absent_keywords.update({'TIMEZERO':0, 'TIMEPIXR':0.5})
        colmap = dict(time='TIME', rate='RATE', rate_err='ERROR')
        if 'FRACEXP' not in self.HDU.columns.names:
            my.printwarn("You are trying to use not correct XMM-Newton light curve "\
                 "which is not recomended. 'FRACEXP' column will be fake.")
            self._do_job(colmap, self.keyword_to_read)
            self.required_columns['fracexp'] = np.ones_like(self.required_columns['time'])
        else:
            self._do_job(colmap|{'fracexp':'FRACEXP'}, self.keyword_to_read)
        

class LCrate():
    """Base class for light curves.
    
    Light curve represented as a pair of columns: (TIME, RATE), optionally 
    errors for both columns can also be included for plotting purposes. Time 
    spacing doesn't matter. Default time properties are: TIMEUNITS='s', TIMESYS='UTC',
    TIMEZERO=0, MJDREF=0.
    """
    
    _colnames_main = ('time', 'time_err', 'rate', 'rate_err') #important columns
    _timekeys_default={'MJDREF':0, 'TIMEZERO':0, 'TIMESYS':'UTC', 'TIMEUNIT':'s'}
    _timekeys_allowed_values={'TIMEUNIT':['s','d'],
                              'TIMESYS':['UTC','TT']}
    
    _precision = 1e-5  #Accuracy for float comparison
    
    def __init__(self, time:ndarray, rate:ndarray, *, time_keywords:dict, 
                 keywords:dict={}, time_err:ndarray=None, rate_err:ndarray=None, 
                 extra_columns:dict[str, ndarray]={}):
        # if values_type not in ['RATE', 'COUNTS']:
        #     raise ValueError("Unknown type of the light curve. It can be "\
        #             "either 'RATE' or 'COUNTS'")

        _columns={}
        refshape = (len(time),)  #reference shape
        
        #Copy columns. We don't want to update() the dictionary, but we want to compare
        #possible duplicates. It's for compatibility with subclasses.
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
        
        #Check presence of time keywords
        if diff:= set(self._timekeys_default).difference(time_keywords):
            raise KeyError(f"Time keywords {diff} not found in 'time_keywords' dict")
        # if diff:=set(time_keywords).difference(self._timekeys_default):
        #     raise KeyError(f"Time keywords {diff} are unknown")
            
        _timekeys={} #Check validity of time keywords
        for key,val in time_keywords.items():
            if key in self._timekeys_allowed_values:
                if val not in self._timekeys_allowed_values[key]:
                    raise ValueError(f"{key} must be in {self._timekeys_allowed_values[key]}.")
            _timekeys['__'+key] = val  #Make all the timekeys readonly
            
        self.keywords = keywords.copy()   #Ordinary  python dict
        self._timekeys = my.ROdict(**_timekeys)
        _tk = self._timekeys.dict()  #Timekeys with normal names (without '__')
        self._mjdref = Time(_tk['MJDREF'], format='mjd', scale=_tk['TIMESYS'].lower())
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
        text += ' | '.join(x.capitalize() for x in cols) #Caption of the table
        indices = list(range(len(self)))
        prev = 0
        for i in set(indices[:10] + indices[-10:]):
            if (i-prev)>1:
                text+='\n'+'\t'.join(['---']*len(cols))
            text+='\n'+'\t'.join(str(self._columns._dict[x][i]) for x in cols)
            prev = i
        return text
            

    def get_column(self, colname:str):
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
    def _(self, index: int):  #Get one point
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
                my.printwarn("The option '{keyerr}={loc[keyerr]}' will be ignored.")

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
        """Return ratio  of BKG to SRC collection areas."""
        # self._check_curve_type(bkgcurve)
        if 'BACKSCAL' in self._keywords and 'BACKSCAL' in bkgcurve._keywords:
            return self._bkgcurve['BACKSCAL']/self._keywords['BACKSCAL']
        else:
            None
            
    def check_synchronicity(self, other, silent=False):
        """Check wether two light curves are synchronous."""
        #TODO check time+TIMZERO first and than check timevectors
        # self._check_curve_type(other)
        printfun =  my.printwarn if not silent else lambda x: None
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
    
    def rebin(self, groups):
        """Rebin the light curve using list of bin groups."""
        if not groups:
            raise ValueError('List of groups cannot be empty')
        if hasattr(groups, 'lcurve'):
            if np.any(self.time != groups.lcurve.time):
                my.printwarn("It seems that this groups was produced for "
                        "a different light curve.")
        
        lc=self.remove_empty()  
        time = lc._columns['time']
        rate = lc._columns['rate']
        width = lc._columns['binwidth']
        
        if 'rate_err' not in lc._columns:
            raise NotImplementedError("Light curves without error column are "\
                      "not supported yet.")
        error = lc._columns['rate_err']
        
        #new arrays
        nlen = len(groups)
        colnames= ['time', 'time_err', 'rate', 'rate_err', 'binwidth']
        cols = {k:np.zeros(nlen) for k in colnames}
        #TODO account for existing time_err
        for i, gr in enumerate(groups):
            cols['time'][i] = 0.5*(time[gr[0]]+time[gr[-1]])
            cols['time_err'][i] = time[gr[-1]]-time[gr[0]] + \
                                       0.5*(width[gr[0]]+width[gr[-1]])
            bw = width[gr].sum()
            # rate[i] = (rt1*bw1 + rt2*bw2 + .. rt_n*bw_n)/sum(bw1...bw_n)
            cols['rate'][i] = np.sum(rate[gr]*width[gr])/bw
            #sqrt ( (s1*b1/smn(bw))^2 + (s2*b2/smn(bw))^2 + ... ) =
            # = 1/sum(bw) * sqrt( (s1*b1)^2 + (s2*b2)^2)
            cols['rate_err'][i] = np.sqrt(np.sum(np.power(error[gr]*width[gr],2)))/bw
            cols['binwidth'][i]=bw
            
        timekeys = self._timekeys.dict()
        return LCrateBinned(time_keywords=timekeys, keywords=self.keywords, **cols)
                
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
            res, classkeys = main.fits_check_hdu_is_lc(HDU, withclasskeys=True) ##TODO action=exceptions
            if classkeys['HDUCLAS3'] == 'COUNTS':
                raise NotImplementedError('Only light curves in rates are '\
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
    

#Set corresponding LCcnt-classes for LCrates     
LCrateBinned._class_cnt = LCcntBinned    
LCrateBinnedEven._class_cnt = LCcntBinnedEven
    
#### General functions

def read(filepath):
    """Factory function to read a light curve from FITS."""
    pass

def rate_uplimit(rawcnt:int, bkgcnt:float, exp:float, prob:float=0.9) ->float:
    """Calculate upper limit to count rate.
    
    Calculate quantile of the count rate corresponding to the given probability.
    Poisson distribution is assumed.

    Parameters
    ----------
    rawcnt : int
        Number of raw counts accumulated in the source aperture.
    bkgcnt : float
        Number of background counts rescaled to the source aperture.
    exp : float
        Exposure time.
    prob : float, optional
        Probability. The default is 0.9 (90% upper limit).

    Returns float value.
    """
    return (scipy.optimize.fsolve(lambda x: scipy.stats.poisson.cdf(rawcnt, \
                bkgcnt+x)-(1-prob), rawcnt-bkgcnt)[0])/exp
                        

def _rebin_prechecks(src, bkg, bkgratio, attrlist:list, collist:list)-> float : 
    """Make checks before rebinning."""
    
    for attr in attrlist:
        if not hasattr(src, attr):
            raise TypeError("This light curve cannot be rebined. (It doesn't have "\
                            f"the '{attr}()' method.)")
    for curve in (src, bkg):
        if curve is not None:
            for col in collist:
                if col not in curve.colnames:
                    raise TypeError("This light curve cannot be rebined. (It doesn't have "\
                                    f"required column '{col}'.)")
    if bkg:
        if not src.check_synchronicity(bkg, silent=False):
            raise TypeError("Can't subtract a non-synchronous light curve.")
        if bkgratio is None:
            if bkgratio:= src.get_bkgratio(bkg) is None:
                raise TypeError("Information about collection areas has not "\
                    "been found in the lcurves' metadata. So the 'bkgratio' argument "\
                    "is mandatory.")
    else:
        bkgratio = 1

    if bkgratio < 1:
        my.printwarn("The provided bkgratio=BKG/SRC < 1. Please, make sure that "\
                     "everything is correct.")
    return bkgratio


def rebin(lcurve, groups, bkgcurve=None, bkgratio=None):
    pass

        

def rebin_mincounts(lcurve, mincounts:int, maxbinwidth:float,
               bkgcurve=None, bkgratio=None, prob=0.9):
    
    bkgratio = _rebin_prechecks(lcurve, bkgcurve, bkgratio, 
                ['getgroups_mincounts'],['counts','binwidth'])
    
    groups = lcurve.getgroups_mincounts(mincounts, maxbinwidth)
    srclc2 = lcurve.rebin(groups)
    rawcnt = srclc2.counts   #Count in 'raw' light curve
    bw = srclc2.binwidth
    
    extracols={}    #dict for constructor
    keys = srclc2.keywords
    
    if bkgcurve:
        bkglc2 = lcurve.rebin(groups)
        bkgcnt = bkglc2.counts
        counts = rawcnt - bkgcnt/bkgratio
        extracols['counts_raw'] = rawcnt
        keys["HDUCLAS2"] = 'NET'
    else:
        bkgcnt=np.zeros_like(rawcnt)
        counts = rawcnt
    uplim = np.zeros_like(bw)      #Must be np.float!
    for i in range(len(rawcnt)):
        uplim[i] = rate_uplimit(rawcnt[i], bkgcnt[i]/bkgratio, bw[i], prob)
    extracols['rate_uplimits'] = uplim
        
    return LCcntBinned(srclc2.time, counts, bw, time_err=srclc2.time_err, 
                    time_keywords=srclc2.time_keywords, keywords=keys,
                    extra_columns=extracols)
    
        


def plot_uplimits(lcurve, *, drop_shorter:float=0, drop_empty_shorter=0.0,
                  uplim_below=4.0, plotpair_timeopt=None, returnobj=None, 
                  print_report=True):
    
    _rebin_prechecks(lcurve, None, None, [], ['counts','rate_uplimits'])
    cnt = lcurve.counts
    bw = lcurve.binwidth
    uplim = lcurve.columns['rate_uplimits']
    nbins = len(cnt)
    
    #Two masks to extract normal rates and uplimits
    mask1, mask2 = [False]*nbins, [False]*nbins
    counters={k:0 for k in ['empty', 'drop_empty', 'drop_short','uplim', 'normal']}
    for i in range(nbins):
        if  cnt[i] == 0:   #Empty bin
            counters['empty'] +=1
            if bw[i] <= drop_empty_shorter:
                counters['drop_empty'] +=1
                continue
        if cnt[i] < uplim_below:
            if bw[i] <= drop_shorter:
                counters['drop_short'] +=1
                continue
            mask2[i] = True
            counters['uplim'] +=1
        else:
            counters['normal'] +=1
            mask1[i] = True
            
    lcrt = lcurve[mask1]
    lcul = lcurve[mask2]
    
    if print_report == True:
        print("Statistics:")
        print("Counts total: {:d}".format(cnt.sum()))
        print("Max counts/bin: {:d}".format(cnt.max()))
        print("Initial length: {:d} bins".format(nbins))
        print("Normal bins: {:d}".format(counters['normal']))
        print("Only uplimins: {:d}".format(counters['uplim']))
        print("Dropped: {:d}".format(counters['drop_empty']+counters['drop_short']))
        print("  Empty: {:d}".format(counters['drop_empty']))
        print("  Short: {:d}".format(counters['drop_short']))
        
    
    if returnobj == 'lcurves':
        return {'rates':lcrt, 'uplimits':lcul}
    
    if not plotpair_timeopt:
        plotpair_timeopt=dict(abstime=True, dates=False)
    
    pprt = lcrt.make_PlotPair(**plotpair_timeopt)
    _ppul = lcul.make_PlotPair(Y='rate_uplimits', **plotpair_timeopt)
    ppul = PlotPairUplimits(_ppul.X, _ppul.Y, arrow_lengths='15%', pltopts=_ppul.pltopts)
    
    if returnobj == 'plotpairs':
        return {'rates':pprt, 'uplimits':ppul}
    else:
        import matplotlib.pyplot as plt
        ppul.oplot(plt, color='green')
        pprt.plot(plt)
        
        


#### legacy #################### 


def _basic_checks(hdu, *, action=Actions.EXCEPTION, refkeys=None):
    """Check the presence of the mandatory keywords for evenly binned light curve."""
    res={}
    for key in ('TIMEPIXR', 'TIMEDEL', 'TIMEZERO'): 
        if key not in hdu.header:
            my._do_action(action, f"Mandatory keyword f'{key}' is not found.",
                exception_class=KeyError)
        res[key] = hdu.header[key]
        
    return res




def rebin_continues_intervals(objfile:str, bkgfile:str=None, *, bkgscale:float=None):
    """Rebin light curve to have one point per each countinuos interval.

    Parameters
    ----------
    objfile : str
        Path to the source light curve FITS file.
    bkgfile : str, optional
        Path to the background light curve FITS file. The default is None
    bkgscale : float, optional
        Background correction factor, i.e. the ratio of areas of the source
        and background apertures. The default is None.

    Returns
    -------
    None.

    """

    
    with fits.open(objfile) as objfts:
        keywords=_make_checks(objfile)
        if keywords['TIMEPIXR']  != 0.5:  ##TODOMove to make_checks
            raise NotImplementedError('Only TIMEPIXR=0.5 is supported yet.')
        step = keywords['TIMEDEL']
        time=objfts[1].data['Time']
        rate=objfts[1].data['Rate'] 
        err=objfts[1].data['Error']  
        
    if bkgfile:
        if bkgscale is None:
            raise TypeError("'bkgfile' argument requires the 'bkgscale' but it is None.")
        with fits.open(bkgfile) as bkgfts:
            b_keywords=_make_checks(objfile)
            for key in keywords:  ## TODO Move this loop to _make_checks refkeys
                if keywords[key] != b_keywords[key]:
                    raise ValueError("The background light curve keyword has "
                        f"different value of the '{key}' keyword.")
            b_rate=bkgfts[1].data['Rate'] 
            b_err=bkgfts[1].data['Error']  
      
    groups_raw = _make_groups(time, step)
    groups=[]
    for gr in groups_raw:
      if len(gr)>1:
        groups.append(gr)   
    
    l=len(groups)
    x=np.zeros(l)
    ex=np.zeros(l)
    yraw=np.zeros(l)
    eyraw=np.zeros(l)
    ynet=np.zeros(l)
    eynet=np.zeros(l)
    ybkg=np.zeros(l)
    eybkg=np.zeros(l)
    i=0
    for i in range(l):
        gr=groups[i]
      
        # print('%d:' % i)
        # print(gr)
      
        # if i==0: print(time[gr])
        x1=time[gr[0]]
        x2=time[gr[-1]] 
        x[i]=0.5*(x1+x2)
        ex[i]=0.5*(x2-x1+step)       #  0.5* ( x2+0.5*STEP - (x1-0.5*STEP) )
        yraw[i]=np.mean(rate[gr])
        ybkg[i]=np.mean(b_rate[gr])*bkgscale
        eyraw[i]=np.sqrt(np.sum( np.power(err[gr],2) ))/len(gr)   # 1/n * sqrt ( s1^2 + s2^2 + ... sn^2 )
        eybkg[i]=np.sqrt(np.sum( np.power(b_err[gr],2) ))/len(gr)*bkgscale
        i+=1
        
    ynet=yraw-ybkg
    eynet=np.sqrt(eyraw**2+eybkg**2) 
    
    i_min=np.argmin(ynet)
    i_max=np.argmax(ynet)
    print(ynet[i_min],ynet[i_max])
    print('Min={:2.3f} at {:e} ({:d} points), Max={:2.3f} at {:e} ({:d} points)'.format(np.amin(ynet),x[i_min],len(groups[i_min]),np.amax(ynet),x[i_max],len(groups[i_max])))
    
    res={}
    for var in ['x','ex','yraw','ybkg','ynet','eyraw','eybkg','eynet']:
      res[var]=locals()[var]
    return res
      
def rebin_smart(objfile:str, mincounts:int, maxbinlength:float):
    from math import sqrt
    with fits.open(objfile) as objfts:
        keywords=_basic_checks(objfts[1])
        step = keywords['TIMEDEL']
        time  = objfts[1].data['Time']
        rate  = objfts[1].data['Rate'] 
        frexp = objfts[1].data['Fracexp'] 

    newlc=[]
    binwidth = buf = t1 = t2 = 0

        
        
    res={}
    res['x']= np.array([pt['time'] for pt in newlc])
    res['yraw']= np.array([pt['rate'] for pt in newlc])
    res['eyraw']= np.array([pt['error'] for pt in newlc])
    return res
        
        