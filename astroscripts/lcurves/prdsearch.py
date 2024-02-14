
import numpy as np
from numpy.typing import ArrayLike, NDArray
from .._internal.columnbased import ColumnBased, _make_vector
from .core import TUnits, LCurve


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
        time = lcurve.time
        time -= time[0]
        T = time[-1]

        # TODO Nyquist frequency from evenly binned light curve
        dt = np.median(np.diff(time))
        freq = np.arange(1/T, 0.5/dt, 1/T)
        return cls(freq, tunits=tunits)

    @property
    def tunits(self):
        return self._tunits

def lspertest(lcurve: LCurve, freq: NDArray | Periodogram = None):

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


