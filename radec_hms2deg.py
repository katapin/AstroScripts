#!/usr/bin/python

import sys
from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy import units as u

data=ascii.read(sys.argv[1],format='no_header')
for r in data:
	coo=SkyCoord(r[0],r[1],unit=(u.hour,u.deg))
	print("%s %s" % (coo.ra.deg, coo.dec.deg))
