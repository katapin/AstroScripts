import ctypes
from mypythonlib.packaging import getrootpath


libpath = getrootpath(__package__) / '_internal/libprdsearch.so'
clib = ctypes.CDLL(libpath)

def call_lspower():
    """Calculate Lomb-Scargle spectrum."""