from ctypes import POINTER, CDLL, c_float, c_uint32, c_ubyte, c_bool
import numpy as np
from numpy.typing import NDArray
from mypythonlib.packaging import getrootpath
import astroscripts


libpath = getrootpath(__package__) / '_internal/libprdsearch.so'
clib = CDLL(libpath)
c_float_p = POINTER(c_float)

def _caller(clibfunc, inX: NDArray, inY: NDArray, outX: NDArray,
                    progress: bool = True) -> NDArray:
    """Call functions from prdsearch.c"""
    if inX.size != inY.size:
        raise ValueError("'time' and 'data' vectors must have the same length.")
    inX = inX - inX[0]
    f_inX = inX.astype(c_float)   # Convert to float (the original array could be int)
    p_inX = f_inX.ctypes.data_as(c_float_p)   # Pointers

    f_inY = inY.astype(c_float)
    p_inY = f_inY.ctypes.data_as(c_float_p)

    f_outX = outX.astype(c_float)
    p_outX = f_outX.ctypes.data_as(c_float_p)

    f_outY = np.zeros_like(outX, dtype=c_float)
    p_outY = f_outY.ctypes.data_as(c_float_p)

    inlen = c_uint32(len(inX))
    outlen =  c_uint32(len(outX))
    nthreads = c_ubyte(astroscripts.CPUCORES)

    res = clibfunc(p_inX, p_inY, p_outX, p_outY, inlen, outlen,
                       c_bool(progress), nthreads)
    if res != len(outX):
        raise RuntimeError("External library function return inconsistent result.")
    return f_outY


def call_lspower(time: NDArray, data: NDArray, freq: NDArray,
                 progress: bool = True) -> NDArray:
    """Calculate Lomb-Scargle spectrum."""
    return _caller(clib.lspower, time, data, freq, progress)


