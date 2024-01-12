#/usr/bin/env python3

"""Plotting facilities."""

from typing import TypeVar
import numpy as np
from numpy import ndarray
import matplotlib.pyplot as plt
from typing import Union
from collections.abc import Iterable
from functools import singledispatch, singledispatchmethod
from mypythonlib import printwarn
from .main import FilePathAbs
from .external import _check_result_file_appeared


def _vec_checks(vec, vecname=None, refshape=None, *, none_is_allowed=False, 
                from_numbers=False) -> ndarray | None:
    """Check if ndarary can be created from the input argument."""
    if vec is None:
        if none_is_allowed == True:
            return None
        else:
            raise ValueError("The vector {} cannot be None.".format(
                   "'"+vecname+"' " if vecname else ''))
            
    if from_numbers == True and refshape:  # Stretch a number to vector
        if not isinstance(vec, (float, int)):
            vec = np.ones(refshape)*vec
            return vec
    if not isinstance(vec, (list, tuple, ndarray)):
        raise TypeError("Data columns must be of a sequence (list, tuple, etc.) "
                        "or the numpy array type.")
    if len(vec) == 0:
        raise ValueError("The vector {}can be None but cannot be empty.".format(
                         "'"+vecname+"' " if vecname else ''))
    if isinstance(vec, (list, tuple)):  #Convert to ndarray
        vec = np.array(vec)
    if refshape and (vec.shape != refshape) or (vec.ndim != 1): 
        raise ValueError("All the arrays must be 1D of the the same length.")
    return vec


class PlotPoint(float):
    """Data point with errors."""
    
    def __new__(cls, *args):
        """Create object inherit from float."""
        return super().__new__(cls, args[0])
    
    def __init__(self, value: float, error_lower: float, error_upper: float):
        self.erlo = error_lower
        self.erup = error_upper
        
    #TODO arithmetic operations taking into account errors


class PlotVector(): 
    """Class to store errorbars together with the main data.
    
    The errorbars are stored as the lower and upper boundaries os the main values. 
    By default, if some points don't satisfy values_lower<=values or
    values_upper>=values, the ValueError exception will be raised. Such a bad 
    points can be replaced by the values 0.5*(lower+upper) if the option 
    safe=True is set.
    
    Parameters
    ----------
    values : ndarray
        Vector with main values.
    values_lower : ndarray, optional
        Vector with lower values. The default is None.
    values_upper : ndarray, optional
        Vector with lower values. The default is None.
    safe : bool, optional
        Correct wrong values instead of raising exception. The default is False.
    """
    
    def __init__(self, values: ndarray, values_lower: ndarray = None, 
                 values_upper: ndarray = None, *, safe: bool = False):
        _d = dict(main=_vec_checks(values, 'values', none_is_allowed=False).copy())
        refshape  = (len(values),)
        _bl, _bu = values_lower is None,  values_upper is None  #convert to bool
        if (_bl or _bu ) != (_bl and _bu) :
            raise TypeError("Both lower and upper values required or neither.")
        for name,vec in {'lower':values_lower, 'upper':values_upper}.items():
            if (vec:= _vec_checks(vec, name, refshape, none_is_allowed=True)) is not None:
                _d[name] = vec.copy()
            else:
                _d[name] = _d['main']
                
        #Treat the problem values
        mask = np.logical_or(_d['lower'] > _d['main'], _d['upper'] < _d['main'])
        if np.any(mask):
            if safe == False:
                raise ValueError("Incorrect values. Errorbars can't have negative lengths.")
            _d['main'][mask] = 0.5*(_d['lower'][mask]+_d['upper'][mask])
        self._values = _d
     
    @classmethod
    def create_with_errors(cls, values: ndarray, errors_lower: Union[ndarray, float], 
                 errors_upper: Union[ndarray, float]):
        """Create vector object with errorbars."""
        main = _vec_checks(values, 'values', none_is_allowed=False).copy()
        refshape  = (len(main),)
        lower = main - _vec_checks(errors_lower, 'errors_lower', refshape,
                       none_is_allowed=False, from_numbers=True)
        
        upper = main + _vec_checks(errors_lower, 'errors_lower', refshape,
                       none_is_allowed=False, from_numbers=True)
        return cls(main, lower, upper)
    
    @classmethod
    def create_with_errors_symmetric(cls, values: ndarray, errors: Union[ndarray, float]):
        """Create vector object with symmetric errorbars."""
        return cls.create_with_errors(values, errors, errors)
            
    def __len__(self):
        """Return vector length."""
        return len(self._values['main'])
    
    @singledispatchmethod 
    def __getitem__(self, arg):
        """Return self[key]."""
        raise TypeError("Unsupported type of key: {}.".format(type(arg)))
        
    @__getitem__.register     
    def _(self, index: int):
        val = self._values['main'][index]
        erlo = val - self._values['lower'][index] 
        erup = self._values['upper'][index]  - val
        return PlotPoint(val, erlo, erup)
    
    @__getitem__.register    
    def __getitem__(self, arg: Union[slice, list, ndarray]):
        obj = self.copy()
        for vecname in ['main', 'lower', 'upper']:
            obj._values[vecname] = obj._values[vecname][arg]
        return obj

    def __add__(self, val):
        """Add constant to each vector value."""
        obj = self.copy()
        for vecname in ['main', 'lower', 'upper']:
            obj._values[vecname] = obj._values[vecname]+val
        return obj
            
    def __mul__(self, val):
        """Multiply each vector value by float."""
        obj = self.copy()
        for vecname in ['main', 'lower', 'upper']:
            obj._values[vecname] = obj._values[vecname]*val
        return obj
				
    def isgood(self):  #TODO
        return not np.all(self.data.mask)  #All values masked, it's an empty array 
    
    @property 
    def has_errors(self):
        """Return True if errors are nonzero."""
        _v = self._values
        return not np.allclose(_v['lower'],_v['upper'], rtol=1e-5, atol=0.0)
            
    
    @property
    def values(self) -> ndarray :
        """Return array with main data values."""
        return self._values['main'].copy()
    
    @property 
    def errors_lower(self) -> ndarray :
        """Return array with lower errors or zero array."""
        return self._values['main'] - self._values['lower']
        
    @property 
    def errors_upper(self) -> ndarray :
        """Return array with upper errors or zero array."""
        return self._values['upper'] - self._values['main'] 
    
    @property 
    def errors_mean(self) -> ndarray :
        """Return array with upper errors or zero array."""
        return 0.5*(self._values['upper'] - self._values['lower'])
        
    @property 
    def values_lower(self) -> ndarray :
        """Return array with lower values."""
        return self._values['lower'].copy()

    @property 
    def values_upper(self) -> ndarray : 
        """Return array with upper values."""
        return self._values['upper'].copy()

    def copy(self):
        """Clone the object."""
        cls = self.__class__
        obj = cls.__new__(cls)
        obj._values = {'main': self.values, 'lower': self.values_lower,
                       'upper': self.values_upper}
        return obj
    

# Type for vectors in PlotPair
TVec = TypeVar('TVec', PlotVector, ndarray, list, tuple)


class _PlotPairBase():
    """Base class for PlotPair and PlotPairUplimits."""
    def __init__(self, X: TVec, Y: TVec, *, pltopts: dict = None):
        self._X = PlotVector(X) if not isinstance(X, PlotVector) else X
        self._Y = PlotVector(Y) if not isinstance(Y, PlotVector) else Y
        if len(X) != len(Y):
            raise ValueError('Vectors must have the same length')
        self.pltopts = pltopts or {}

    def __len__(self):
        """Return light of the vectors."""
        return len(self._X)
    
    def copy(self):
        """Clone the object."""
        cls = self.__class__
        obj = cls.__new__(cls)
        obj._X = self._X.copy()
        obj._Y = self._Y.copy()
        obj.pltopts = self.pltopts.copy()
        return obj
    
    @singledispatchmethod
    def __getitem__(self, arg):
        """Return self[key]."""
        raise TypeError("Unsupported type of key: {}.".format(type(arg)))
        
    @__getitem__.register
    def _(self, index: int):
        raise NotImplementedError()
    
    @__getitem__.register
    def _(self, arg: Union[slice, list, ndarray]):
        obj = self.copy()
        obj._X = obj._X[arg]
        obj._Y = obj._Y[arg]
        return obj
    
    @property 
    def X(self):
        """Return vector for the X axis."""
        return self._X
    
    @property 
    def Y(self):
        """Return vector for the Y axis."""
        return self._Y
    
    def rescaleX(self, factor):
        """Multiply the X vector by 'factor'."""
        self._X = self._X*factor
        
    def rescaleY(self, factor):
        """Multiply the Y vector by 'factor'."""
        self._Y = self._Y*factor
        
    def moveX(self, shift):
        """Add 'shift' to each element of the X vector."""
        self._X = self._X + shift
    
    def moveY(self, shift):
        """Add 'shift' to each element of the Y vector."""
        self._Y = self._Y + shift
    
    def _plot(self, pltopts:dict, plthandle, drawnow:bool):
        
        # _non_errorbars_options_mapping = ('xlabel', 'ylabel', 'title', 'xscale', 'yscale')
        opts = self.pltopts.copy()          # Copy to not change the original dict
        if pltopts: opts.update(pltopts)    # Override only for current plot
        
        if plthandle is None: 
            plthandle=plt
            drawnow = True
            
        ebopts={}       # Options for pyplot.errorbar
        popts={}        # Options for pyplot itself: title, labels, etc.
        for key,val in opts.items():  # Split plot options
            if hasattr(plt, key):
                popts[key]=val                   
            else:
                ebopts[key]=val
                
        # Plot!
        plthandle.errorbar(self.X.values, self.Y.values,
                xerr=(self.X.errors_lower, self.X.errors_upper),
                yerr=(self.Y.errors_lower, self.Y.errors_upper), **ebopts)
        
        for key, val in popts.items():
            if val is not None:
                getattr(plt, key)(val)  # Apply pyplot options
            else:
                getattr(plt, key)()

        if drawnow is True:
            plt.show()
    
    def plot(self, plthandle=None, **pltopts):
        """Create figure and immediately show it."""
        self._plot(pltopts, plthandle=plthandle, drawnow=True)
        
    def oplot(self, plthandle, **pltopts):
        """Overplot the existing figure."""
        self._plot(pltopts, plthandle=plthandle, drawnow=False)


class PlotPair(_PlotPairBase):
    """Pair of vectors to create a plot.
    
    It stores values for both X and Y axes as well as the metadata for plotting.
    Besides the main values each vector can also store lengths of errorbars that
    will be automatically passed to pyplot.errorbars() task.
    
    Parameters
    ----------
    X : PlotVector
        Values for X axis
    Y : PlotVector
        Values for X axis
    pltopts : dict, optional
        Dict of options to pass pyplot and pyplot.errorbars.
    """

    pass
    # def __init__(self, X: PlotVector, Y: PlotVector, *, pltopts: dict = None):
    #     super().__init__(X, Y, pltopts=pltopts)

    # TODO make stepped line plot


class PlotPairUplimits(_PlotPairBase):
    """PlotPair version for plotting upper limits.
    
    Upper limits will be plotted as arrows. Arrow lengths can be tuned with the
    parameter 'arrow_lengths'. It can be expressed either as an absolute (constant)
    value (useful for plots with linear scale) or as a fraction (percentages) of 
    the main upper limit values (for log-scale plots).
    """
    
    def __init__(self, X: TVec, Y_uplim: TVec, *,
                 arrow_lengths: Union[float, str], pltopts: dict = None):
        
        super().__init__(X, Y_uplim, pltopts=pltopts)
        if self._Y.has_errors is True:
            printwarn("Errors of 'Y_uplim' will be ignored.")
        self.arrow_lengths = arrow_lengths  # This calls setter
        self.pltopts.update(dict(uplims=True, fmt='_', capsize=0, elinewidth=1))
        
    def copy(self):
        """Clone the object."""
        obj = super().copy()
        obj._arrow_lengths = self. _arrow_lengths
        return obj
        
    @property
    def arrow_lengths(self):
        """Length of arrows for uplimits."""
        return self._arrow_lengths
    
    @arrow_lengths.setter
    def arrow_lengths(self, val):
        vec = self._Y.values
        if isinstance(val, str):  # '15%'
            if len(val) < 2 or val[-1] != '%':
                raise ValueError(f"Wrong format of 'arrow_lengths': '{val}'")
            factor = float(val[:-1])/100   # 0.15
            self._Y = PlotVector.create_with_errors_symmetric(vec, vec*factor)
        elif isinstance(val, (float, int)):
            self._Y = PlotVector.create_with_errors_symmetric(vec, val)
        else:
            raise TypeError("Wrong time of 'arrow_lengths'")
        self._arrow_lengths = val


@singledispatch
def plot_and_save(data, imgpath: FilePathAbs | None = None, *, onlysave: bool= False,
                  progname: str = None, **pltopts):
    """Show the figure and save image to file."""
    raise TypeError("Unsupported type of the 'data' argument.")


@plot_and_save.register(list | tuple)
def _plot_and_save_list(data: list[PlotPair] | tuple[PlotPair],
                       imgpath: FilePathAbs | None = None, *,
                       onlysave: bool= False, progname: str = None, **pltopts):
    data[-1].pltopts.update(pltopts)
    for pp in data:
        pp.oplot(plt)

    if imgpath:
        fig = plt.gcf()
        fig.savefig(imgpath.fspath)
        _check_result_file_appeared(imgpath, progname)

    if onlysave:
        if not imgpath:
            raise TypeError("Argument 'onlysave=True' requires 'imgpath' but "
                            "None is given.")
        plt.clf()
    else:
        plt.show()


@plot_and_save.register
def _plot_and_save_single(data: PlotPair | None, imgpath: FilePathAbs = None, *,
                         onlysave: bool= False, progname: str = None, **pltopts):
    _plot_and_save_list([data] if data else [], imgpath, onlysave=onlysave,
                       progname=progname, **pltopts)

