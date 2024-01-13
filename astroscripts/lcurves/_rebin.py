
import numpy as np
from mypythonlib import printwarn
from ._classes import LCrateBinned, LCcntBinned
from ._common import _bkgration_msg
from .statistics import rate_uplimit
from ..plot import PlotPairUplimits

# TODO list
# 1)Class from xmm light curve that can extract the background LC from
# the BACKV and BACKE columns in NET lc-FITS.

"""
Some thoughts about rebinning algorithms.

About TIMEPIXR.
Let's suppose timestep is 2s. If TIMEPIXR=0.5 the timestamps will's as follows
 | bin-1 | bin-2 | bin-3 | bin-4 | bin-5 | bin-6 | bin-7 | bin-8 |....   
 |   ^   |   ^   |   ^   |   ^   |   ^   |   ^   |   ^   |   ^   |
     0s      2s      4s      6s  |   8s  |  10s  |  12s  |  14s  |  
The TIMEZERO keyword in the FITS  header refers to the timstamp 0s (and not to
the star time of the exposure!). If we, for example, everage 4 bins, we'll obtain:
 |           long bin            |           long bin            |.... 
 |   ^           ^               |               ^               |           
   timezero      3s              |              11s              |
Now the TIMEZERO no longer refers to the bin center, but timestamps are 
correct and still refer the bin center. Thus, only TIMEPIXR=0.5 remain valid
after the rebinning


###########################
About getgroups_continous()
Let's suppose maxgap = 1 (in units of TIMEDEL), one bin can be missed 

detector in operation: |+++|+++|+++|   |+++|+++|   |   |   |+++|+++|   |   |+++|+++|+++|   |   |+++|| END
   time              : | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10| 11| 12| 13| 14| 15| 16| 17| 18|| END
   indixes           : | 0 | 1 | 2 |   | 3 | 4 |           | 5 | 6 |       | 7 | 8 | 9 |       | 10||

Data (after removing empty bins)
array index  : 0 1 2 3 4 5  6  7  8  9  10
time column  : 0 1 2 4 5 9  10 13 14 15 18
time[1:]     : 1 2 4 5 9 10 13 14 15 18
time[:-1]    : 0 1 2 4 5 9  10 13 14 15
time diff    : 1 1 2 1 4 1  3  1  1  3
Indices of diff>(maxgap+1)  : 4 6 9
group0 = range(0, 4+1)  - > 0 1 2 3 4
group1 = range(4+1, 6+1) -> 5 6
group2 = range(6+1, 9+1) -> 7 8 9
group3 = range(9+1, len) -> 10


############################
About getgroups_mincounts()

detector in operation: |+++|+++|+++|+++|+++|+++|+++|+++|+++|   |+++|+++|+++|+++|+++|   |+++|+++|   |+++|+++|+++|+++|+++||END
photons arrived      : |   | 1 |   | 2 | 2 |   | 2 |   |   |   |   |   |   |   |   |   |   | 1 |   | 1 | 1 |   |   |   ||END 
   time              : | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10| 11| 12| 13| 14| 15| 16| 17| 18| 19| 20| 21| 22| 23||END
   indexes           : | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |   | 9 | 10| 11| 12| 13|   | 14| 15|   | 16| 17| 18| 19| 20||END
Let's add a photon far-far away (e.g. at time=1000) for algorithmic purposes
Cells with arrived counts > 0 
indices of photons : 1 3 4 6 15 16 17 21
timestamps         : 1 3 4 6 17 19 20 1000
Case 1) mincounts = 4, maxbinlength=5s  
Starting...
Set gs=0 -- index of group start, prev = -1 --- index of previous photon
Begin k-loop over list with photon indexes 
prev=-1, k=1, put photon to 'buff', because 'buff'<4 and time[k]-time[gs]+bw < maxbinlength
prev=1, k=3, put photons 
prev=3, k=4, put photons because before if 'buff' was still <4
prev=4, k=6, buf=5. Set ge=prev, group.append(range(gs,ge+1)), set buff=0,
  set gs=ge+1, put photons to buf 
prev=6, k=15, buf=2, time[k]-time[gs]+bw> maxbinlength
  Begin loop while (time[gs] - time[k]+bw) > maxbinlength
    Begin reverse j-loop from k-1(=14) to gs(=5)
    j=14, time[14]=16, time[gs]=5, (14-5+bw)>5
    ....
    j=8, time[8]=8,  time[gs]=5, (10-8+bw)=4
        set ge=j, group2.append(range(gs,ge+1)), set buf=0, set gs=ge+1=9,
        j-loop ended
    Begin reverse j-loop from k-1=(14) to gs(=9)
    j=14, time[14]=16, time[gs]=10, (16-10+bw)>5
    j=13, time[14]=13, time[gs]=10, (13-10+bw)=5
        set ge=j, group2.append(range(gs,ge+1)), set buf=0, set gs=ge+1=14,
        j-loop ended
    while loop ended
  put photon at k=15 to buf, buf=1 
prev=15, k=16, put photon
prev=16, k=17, put photon
prev=17, k=21, its the fake photon with time[21]=1000
  Begin loop while (time[gs] - time[k]+bw) > maxbinlength
    Begin reverse j-loop from k-1(=20) to gs(=14)
    j=20, time[23]=24, time[gs]=16
    ....
    j=17

"""


def _rebin_prechecks(src, bkg, bkgratio: float,
                     attrlist: list, collist: list) -> float:
    """Make checks before rebinning."""

    for attr in attrlist:
        if not hasattr(src, attr):
            raise TypeError("This light curve cannot be rebined. (It doesn't have the "
                            f"the '{attr}()' method.)")
    for curve in (src, bkg):
        if curve is not None:
            for col in collist:
                if col not in curve.colnames:
                    raise TypeError("This light curve cannot be rebined. (It doesn't have the "
                                    f"required column '{col}'.)")
    if bkg:
        if not src.check_synchronicity(bkg, silent=False):
            raise TypeError("Can't subtract a non-synchronous light curve.")
        if bkgratio is None:
            if bkgratio := src.get_bkgratio(bkg) is None:
                raise TypeError(
                    "Information about collection areas has not "
                    "been found in the lcurves' metadata. So the 'bkgratio' argument "
                    "is mandatory.")
    else:
        bkgratio = 1

    _bkgration_msg(bkgratio)

    return bkgratio


def rebin(self, groups):
    """Rebin the light curve using list of bin groups."""
    if not groups:
        raise ValueError('List of groups cannot be empty')
    if hasattr(groups, 'lcurve'):
        if np.any(self.time != groups.lcurve.time):
            printwarn("It seems that this groups was produced for "
                         "a different light curve.")

    lc = self.remove_empty()
    time = lc._columns['time']
    rate = lc._columns['rate']
    width = lc._columns['binwidth']

    if 'rate_err' not in lc._columns:
        raise NotImplementedError("Light curves without error column are " \
                                  "not supported yet.")
    error = lc._columns['rate_err']

    # new arrays
    nlen = len(groups)
    colnames = ['time', 'time_err', 'rate', 'rate_err', 'binwidth']
    cols = {k: np.zeros(nlen) for k in colnames}
    # TODO account for existing time_err
    for i, gr in enumerate(groups):
        cols['time'][i] = 0.5 * (time[gr[0]] + time[gr[-1]])
        cols['time_err'][i] = time[gr[-1]] - time[gr[0]] + \
                              0.5 * (width[gr[0]] + width[gr[-1]])
        bw = width[gr].sum()
        # rate[i] = (rt1*bw1 + rt2*bw2 + .. rt_n*bw_n)/sum(bw1...bw_n)
        cols['rate'][i] = np.sum(rate[gr] * width[gr]) / bw
        # sqrt ( (s1*b1/smn(bw))^2 + (s2*b2/smn(bw))^2 + ... ) =
        # = 1/sum(bw) * sqrt( (s1*b1)^2 + (s2*b2)^2)
        cols['rate_err'][i] = np.sqrt(np.sum(np.power(error[gr] * width[gr], 2))) / bw
        cols['binwidth'][i] = bw

    timekeys = self._timekeys.dict()
    return LCrateBinned(time_keywords=timekeys, keywords=self.keywords, **cols)


def rebin_mincounts(lcurve, mincounts: int, maxbinwidth: float,
                    bkgcurve=None, bkgratio=None, prob=0.9):
    bkgratio = _rebin_prechecks(lcurve, bkgcurve, bkgratio,
                                ['getgroups_mincounts'], ['counts', 'binwidth'])

    groups = lcurve.getgroups_mincounts(mincounts, maxbinwidth)
    srclc2 = lcurve.rebin(groups)
    rawcnt = srclc2.counts  # Count in 'raw' light curve
    bw = srclc2.binwidth

    extracols = {}  # dict for constructor
    keys = srclc2.keywords

    if bkgcurve:
        bkglc2 = lcurve.rebin(groups)
        bkgcnt = bkglc2.counts
        counts = rawcnt - bkgcnt / bkgratio
        extracols['counts_raw'] = rawcnt
        keys["HDUCLAS2"] = 'NET'
    else:
        bkgcnt = np.zeros_like(rawcnt)
        counts = rawcnt
    uplim = np.zeros_like(bw)  # Must be np.float!
    for i in range(len(rawcnt)):
        uplim[i] = rate_uplimit(rawcnt[i], bkgcnt[i] / bkgratio, bw[i], prob)
    extracols['rate_uplimits'] = uplim

    return LCcntBinned(srclc2.time, counts, bw, time_err=srclc2.time_err,
                       time_keywords=srclc2.time_keywords, keywords=keys,
                       extra_columns=extracols)


def plot_uplimits(lcurve, *, drop_shorter: float = 0, drop_empty_shorter=0.0,
                  uplim_below=4.0, plotpair_timeopt=None, returnobj=None,
                  print_report=True):
    _rebin_prechecks(lcurve, None, None, [], ['counts', 'rate_uplimits'])
    cnt = lcurve.counts
    bw = lcurve.binwidth
    uplim = lcurve.columns['rate_uplimits']
    nbins = len(cnt)

    # Two masks to extract normal rates and uplimits
    mask1, mask2 = [False] * nbins, [False] * nbins
    counters = {k: 0 for k in ['empty', 'drop_empty', 'drop_short', 'uplim', 'normal']}
    for i in range(nbins):
        if cnt[i] == 0:  # Empty bin
            counters['empty'] += 1
            if bw[i] <= drop_empty_shorter:
                counters['drop_empty'] += 1
                continue
        if cnt[i] < uplim_below:
            if bw[i] <= drop_shorter:
                counters['drop_short'] += 1
                continue
            mask2[i] = True
            counters['uplim'] += 1
        else:
            counters['normal'] += 1
            mask1[i] = True

    lcrt = lcurve[mask1]
    lcul = lcurve[mask2]

    if print_report is True:
        print("Statistics:")
        print("Counts total: {:d}".format(cnt.sum()))
        print("Max counts/bin: {:d}".format(cnt.max()))
        print("Initial length: {:d} bins".format(nbins))
        print("Normal bins: {:d}".format(counters['normal']))
        print("Only uplimins: {:d}".format(counters['uplim']))
        print("Dropped: {:d}".format(counters['drop_empty'] + counters['drop_short']))
        print("  Empty: {:d}".format(counters['drop_empty']))
        print("  Short: {:d}".format(counters['drop_short']))

    if returnobj == 'lcurves':
        return {'rates': lcrt, 'uplimits': lcul}

    if not plotpair_timeopt:
        plotpair_timeopt = dict(abstime=True, dates=False)

    pprt = lcrt.make_PlotPair(**plotpair_timeopt)
    _ppul = lcul.make_PlotPair(Y='rate_uplimits', **plotpair_timeopt)
    ppul = PlotPairUplimits(_ppul.X, _ppul.Y, arrow_lengths='15%', pltopts=_ppul.pltopts)

    if returnobj == 'plotpairs':
        return {'rates': pprt, 'uplimits': ppul}
    else:
        import matplotlib.pyplot as plt
        ppul.oplot(plt, color='green')
        pprt.plot(plt)

