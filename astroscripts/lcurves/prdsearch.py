
import sys
from math import atan2, pi
from enum import StrEnum
import numpy as np
from numpy.typing import ArrayLike, NDArray
from .._internal.columnbased import ColumnBased, _make_vector
from .core import TUnits, LCurve
from . import _c_callers as clib


class Backend(StrEnum):
    CLIB = 'clib'
    NUMPY = 'numpy'


class Periodogram(ColumnBased):

    def __init__(self, freq: ArrayLike, data: ArrayLike = None,
                 tunits: TUnits = TUnits.sec):
        self._tunits = tunits
        freqvec = _make_vector(freq, 'freq')
        datavec = _make_vector(data or 0.0, 'data', from_numbers=True,
                               reflen=len(freqvec))
        super().__init__({'data': datavec}, {'freq':freqvec})

    @classmethod
    def from_lcurve(cls, lcurve: LCurve):
        tunits = lcurve.tunits
        T = lcurve.T

        # TODO Nyquist frequency from evenly binned light curve
        dt = np.median(np.diff(lcurve.time))
        freq = np.arange(1/T, 0.5/dt, 1/T)
        return cls(freq, tunits=tunits)

    @property
    def tunits(self):
        return self._tunits

    @property
    def periods(self):
        return 1/self.freq


def lspertest(lcurve: LCurve, freq: NDArray | Periodogram = None,
              progress: bool = True, backend: Backend = Backend.CLIB):
    """Calculate LS-spectrum in Horne & Baliunas (1986) normalization."""
    # Prepare a new Periodogram object with correct tunits
    if freq is None:
        prd = Periodogram.from_lcurve(lcurve)
    elif isinstance(freq, Periodogram):
        prd = freq.copy()
    elif isinstance(freq, np.ndarray):
        prd = Periodogram(freq, tunits=lcurve.tunits)
    else:
        raise TypeError(f"Unknown type of the 'freq' argument ({type(freq).__name__}).")

    time = lcurve.timevec.get_ticks(prd.tunits)
    prd.columns['data'] = {
        Backend.CLIB:  clib.call_lspower,
        Backend.NUMPY: _npbacked_lspower
    }[backend](time, lcurve.data, prd.freq, progress)
    return prd


def _npbacked_lspower(time: NDArray, data: NDArray, freq: NDArray,
                      progress: bool = True) -> NDArray:
    """Calculate Lomb-Scargle spectrum with numpy."""
    ls = np.zeros_like(freq)
    time = time - time[0]
    print(f"Progress: 00%", end='')
    sys.stdout.flush()
    for i, fr in enumerate(freq):
        arg0 = 4*np.pi*fr*time
        omtau = 0.5*atan2(np.sum(np.sin(arg0)), np.sum(np.cos(arg0)))
        arg=2*fr*pi*time-omtau
        COS = np.cos(arg)
        SIN = np.sin(arg)
        yCOS = data @ COS
        ySIN = data @ SIN
        ls[i] = yCOS**2/(COS @ COS) + ySIN**2/(SIN @ SIN)
        if i % int(freq.size/100) == 0:
            print("\b\b\b\b{:3.0f}%".format(100*i/freq.size), end='')
            sys.stdout.flush()
    print("\b\b\b\b{:3.0f}%".format(100))
    sys.stdout.flush()
    return ls

