#!/usr/bin/python
#
#Z-stat

#import scipy.stats, general_random
from scipy.interpolate import interp1d
import numpy as np
from numpy import *
from myfits import *
import matplotlib.pyplot as plt

import sys


    
        
def zfun(phases, nhar=1):
    """calculates Z statistics for phases distributed
    in [0-1] interval. 
    """
    N = len(phases)
    return 2./N*sum([np.sum(cos(2*pi*x*phases))**2+np.sum(sin(2*pi*x*phases))**2 for x in range(1,nhar+1)])

def vz(toa, dt=False, **kwargs):
    """docstring for vz"""
    duration = toa.max()-toa.min()
    if not dt:
        dt = duration/len(toa)
    fnq = 0.5/dt    
    df = 1./duration
    trial_freq = arange(df,fnq,df)
    trial_periods = (1./trial_freq)
    phases = (toa-toa[0])/trial_periods[:,None]%1.
    periods, z =  trial_periods, apply_along_axis(lambda x: zfun(x,**kwargs),1,phases)
    #returns trial periods, z value and full probability that given peak is due to the presence of periodic sigal
    return periods, z, (1.-exp(-z/2.))**len(z)
    
if __name__ == '__main__':   
    lc=sys.argv[1]
    ftslc=fitsopen(lc)
    time=ftslc[1].data['TIME']
    (period, z , ex) = vz(time)
    plt.plot(period,ex)
    plt.show()
