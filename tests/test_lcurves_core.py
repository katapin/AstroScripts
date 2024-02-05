import pytest
import numpy as np
from numpy import array
from astropy.time import Time
from astroscripts.lcurves import _classes as lc


tz = Time(60000, format='mjd')  # default timezero


class Test_Time:
    @pytest.mark.parametrize(
        "ticks, tunits, timezero, output",
        [
            (array([0, 1, 2]), 's', tz, dict(ticks=array([0, 1, 2]), timezero=tz)),
            (array([0, 1, 2]), 's', None, dict(ticks=array([0, 1, 2]), timezero=None)),
            (array([2, 3, 4]), 's', None, dict(ticks=array([2, 3, 4]), timezero=None)),
            (array([0, 1, 2]), 'd', tz, dict(ticks=array([0, 86400, 2*86400]), timezero=tz)),
            (array([86400, 86401, 86402]), 's', tz,
             dict(ticks=array([0, 1, 2]), timezero=Time(60001, format='mjd')))
        ]
    )
    def test_init_valid(self, ticks, tunits, timezero, output):
        tv = lc.TimeVector(ticks, units=tunits, timezero=timezero)
        assert np.all(tv.tickss == output['ticks'])
        if timezero is not None:
            assert abs((tv.timezero - output['timezero']).sec) <= 1e-9  # compare float values
        else:
            assert tv.timezero is output['timezero']

    @pytest.mark.parametrize(
        "ticks, tunits, timezero, excls",
        [
            (array([0, 1, 2]), 's', '60000', TypeError),
            (array([0, 1, 2]), 'min', tz, ValueError),
            (array([1, 0, 2]), 's', tz, ValueError),
        ]
    )
    def test_init_invalid(self, ticks, tunits, timezero, excls):
        with pytest.raises(excls) as exception:
            lc.TimeVector(ticks, units=tunits, timezero=timezero)

    @pytest.mark.parametrize(
        "ticks1, timezero1, ticks2, timezero2, res",
        [
            ([0, 1, 2], tz, [0, 1, 2], tz, True),
            ([0, 1, 2], tz, [0, 2, 3], tz, False),
            ([0, 1, 2], tz, [2, 3, 4], tz, False),
            ([0, 1, 2], tz, [0, 1, 2], None, False),
            ([0, 1, 2], None, [0, 1, 2], None, True),
            ([0, 1, 2], Time(60002, format='mjd'), [2, 3, 4], tz, True),
        ]
    )
    def test_eq(self, ticks1, timezero1, ticks2, timezero2, res):
        tv1 = lc.TimeVector(ticks1, 'd', timezero1)
        tv2 = lc.TimeVector(ticks2, 'd', timezero2)
        assert ((tv1 == tv2) == res)

    @pytest.mark.parametrize(
        "ticks, timezero, text",
        [
            (array([0, 1, 2]), tz,
             "TimeVector(timezero='MJD60000.0', ticks=array([0, 1, 2])"),
            (array([0, 1, 2]), None,
             "TimeVector(timezero=None, ticks=array([0, 1, 2])"),
        ]
    )
    def test_repr(self, ticks, timezero, text):
        tv = lc.TimeVector(ticks, 's', timezero)
        assert repr(tv) == text

    def test_tickss(self):
        inar = array([0, 1, 2])
        tv=lc.TimeVector(inar, 's', tz)
        ticks = tv.tickss
        assert np.all(ticks == inar)
        ticks[1] = 55
        assert np.all(tv.tickss == inar)  # Check internal array didn't change

    def test_ticksd(self):
        inar = array([0, 1, 2])
        tv=lc.TimeVector(inar, 'd', tz)
        ticks = tv.ticksd
        assert np.all(ticks == inar)
        ticks[1] = 55
        assert np.all(tv.ticksd == inar)  # Check internal array didn't change

