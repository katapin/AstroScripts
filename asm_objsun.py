#!/usr/bin/python
#
#A few string about this program


import sys
import os
import re
from astropy.io import fits
import astropy
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import get_sun
from astropy.time import Time

import matplotlib.pyplot as plt
import numpy as np


Css433 = SkyCoord('19 11 49.57 +04 58 57.9', unit=(u.hourangle, u.deg))

x=[]
y=[]



for MJD in np.arange(50090,50095,0.01):
    time=Time(MJD, format='mjd', scale='utc')
    Csun=get_sun(time)
    sep=Css433.separation(Csun)
    (d,m,c)=sep.dms
    x.append(MJD)
    y.append(d+m/60+c/3600)
    print("%s %s" %(MJD,y[-1]))
    
plt.plot(x, y,'b.-')
plt.show()

ind=y.index(min(y))
print("%s %s" % (x[ind], y[ind]))

