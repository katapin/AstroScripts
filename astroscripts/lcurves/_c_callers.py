import ctypes
from ctypes import c_float, c_uint32 
import numpy as np
from numpy.typing import NDArray
from mypythonlib.packaging import getrootpath


libpath = getrootpath(__package__) / '_internal/libprdsearch.so'
clib = ctypes.CDLL(libpath)
c_float_p = ctypes.POINTER(ctypes.c_float)

def call_lspower(time: NDArray, data: NDArray, freq: NDArray):
    """Calculate Lomb-Scargle spectrum."""
    f_time = time.astype(c_float)   # Convert to float (the original array could be int)
    p_time = f_time.ctypes.data_as(c_float_p)   # Pointers 

    f_data = data.astype(c_float)  
    p_data = f_data.ctypes.data_as(c_float_p) 

    f_freq = freq.astype(c_float)   
    p_freq = f_freq.ctypes.data_as(c_float_p) 

    f_power = np.zeros_like(freq, dtype=c_float)
    p_power = f_power.ctypes.data_as(c_float_p) 

    inlen = c_uint32(len(time))
    outlen =  c_uint32(len(freq)) 

    clib.lspower(p_time, p_data, p_freq, p_power, inlen, outlen)
    print(f_power)




    