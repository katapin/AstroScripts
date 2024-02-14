import pytest
import numpy as np
from numpy import array
from astroscripts._internal.columnbased import (
    _make_vector,
    ColumnBased,
    ColumnBasedMinimal,
)

class Test_make_vector():

    @pytest.mark.parametrize(
        "input,expected,noneallowed,from_numbers,reflen",
        [
            ([1, 2, 3], array([1, 2, 3]), False, False, None),
            ([1, 2, 3], array([1, 2, 3]), True, False, None),
            ([1, 2, 3], array([1, 2, 3]), False, True, None),
            ([1, 2, 3], array([1, 2, 3]), False, True, 3),
            ((1, 2, 3), array([1, 2, 3]), False, False, None),
            (array([1, 2, 3]), array([1, 2, 3]), False, False, None),
            (iter([1, 2, 3]), array([1, 2, 3]), None, False, None),
            (range(5), array([0, 1, 2, 3, 4]), None, False, None),
            (None, None, True, False, None),
            (2, array([2, 2, 2]), False, True, 3),
        ]
    )
    def test_valid_input(self, input, expected, noneallowed, from_numbers, reflen):
        output = _make_vector(input, none_is_allowed=noneallowed,
                              from_numbers=from_numbers, reflen=reflen)
        print(type(output))
        assert type(output) == type(expected)
        assert np.all(output == expected)

    @pytest.mark.parametrize(
        "input, noneallowed, from_numbers, excls, reflen",
        [
            (None, False, False, TypeError, None),
            (None, False, True, TypeError, None),
            ('qwe', True, False, TypeError, None),
            ('qwe', True, True, TypeError, None),
            # ([], False, False, ValueError, None),     # cannot be empty
            # ([], True, False, ValueError, None),
            ([1,2,3], False, False, ValueError, 4),   # wrong length
            (array([[1,2,3],[4,5,6]]), False, False, ValueError, None),  # 2D is not allowed
            ([[1,2,3],[4,5,6]], False, False, ValueError, None),         # as well as nested lists
            (5, False, True, TypeError, None),        # from numbers but reflen is not passed
        ]
    )
    def test_invalid_input(self, input, noneallowed, from_numbers, excls, reflen):
        with pytest.raises(excls) as exception:
            _make_vector(input, vecname='values', none_is_allowed=noneallowed,
                         from_numbers=from_numbers, reflen=reflen)


def mkcols(names: list[str], data):
    """Prepare columns for tests."""
    return {n: (_make_vector(a if not isinstance(a, int) else range(a)))
            for n,a in zip(names, data)}


class TestColumnBasedMinimal():
    ColumnBasedMinimal.__abstractmethods__ = set()

    @pytest.mark.parametrize(
        "cols, res",
        [
            (mkcols(['rate'], [5]), 5),
            (None, 0)  # Empty object
        ]
    )
    def test_len(self, cols, res):
        cbm = ColumnBasedMinimal(cols)
        assert len(cbm) == res

    @pytest.mark.parametrize(
        "cols, index, res",
        [
            (mkcols(['X', 'Y'], [5, 5]), 2, {'X': 2, 'Y': 2}),
        ]
    )
    def test_getitem_point(self, cols, index, res):
        cbm = ColumnBasedMinimal(cols)
        assert cbm[index] == res

    @pytest.mark.parametrize(
        "cols, mask, res",
        [
            (  # boolean mask
                mkcols(['X', 'Y'], [5, 5]),
                [False, True, True, False, False],
                mkcols(['X', 'Y'], [[1, 2], [1, 2]])
            ),
            (  # indexing array
                    mkcols(['X', 'Y'], [5, 5]),
                    [1, 2],
                    mkcols(['X', 'Y'], [[1, 2], [1, 2]])
            ),
            (  # indexing array
                    mkcols(['X', 'Y'], [5, 5]),
                    slice(1, 3),
                    mkcols(['X', 'Y'], [[1, 2], [1, 2]])
            ),
        ]
    )
    def test_getitem_subsample(self, cols, mask, res):
        cbm = ColumnBasedMinimal(cols)
        cmb2 = ColumnBasedMinimal(res)
        new = cbm[mask]
        assert all(list(new._columns['X']) == res['X'])
        assert new == cmb2

    def test_copy(self):
        # TODO
        pass



class TestColumnBased():
    ColumnBased.__abstractmethods__ = set()
    @pytest.mark.parametrize(
        "cols, pcols",
        [
            (mkcols(['rate'], [5]), None),    # no protected
            (None, mkcols(['time'], [5])),    # no normal
            (mkcols(['rate'], [5]), mkcols(['time'], [5])),    # expose
        ]
    )
    def test_init(self, cols, pcols):
        cb = ColumnBased(cols, pcols)
        cols, pcols = cols or {}, pcols or {}
        allcols = list(cols.keys()) + list(pcols.keys())
        assert list(cb._columns) == allcols   # Normal and protected are stored together
        assert list(cb.pcolnames) == list(pcols.keys())

    @pytest.mark.parametrize(
        "cols, pcols, cols2, pcols2, res",
        [
            (
                mkcols(['rate'], [5]),
                None,
                mkcols(['rate'], [5]),
                None,
                True
            ), # OK
            (
                mkcols(['rate'], [5]),
                None,
                mkcols(['rate'],[[1, 1, 1, 1, 1]]),
                None,
                False
            ),  # Different data
            (
                    mkcols(['rate'], [5]),
                    None,
                    mkcols(['rate'], [4]),
                    None,
                    False
            ),  # Different length
            (
                mkcols(['rate'], [5]),
                None,
                mkcols(['rate2'], [5]),
                None,
                False
            ),  # Different names
            (
                    mkcols(['rate'], [5]),
                    mkcols(['time'], [5]),
                    mkcols(['rate'], [5]),
                    None,
                    False
            ),  # Different protected columns
            (
                    mkcols(['rate'], [5]),
                    mkcols(['time'], [5]),
                    mkcols(['rate'], [5]),
                    mkcols(['time'], [5]),
                    True
            ),  # OK
            (
                    mkcols(['rate', 'time'], [5, 5]),
                    None,
                    mkcols(['rate'], [5]),
                    mkcols(['time'], [5]),
                    False
            ),  # Names are the same but not protected
        ]
    )
    def test_eq(self, cols, pcols, cols2, pcols2, res):
        cb = ColumnBased(cols, pcols)
        cb2 = ColumnBased(cols2, pcols2)
        assert (cb == cb2) is res

    @pytest.mark.parametrize(
        "cols, pcols",
        [
            (mkcols(['rate'], [5]), None),    # no protected
            (None, mkcols(['time'], [5])),    # no normal
            (mkcols(['rate'], [5]), mkcols(['time'], [5])),    # expose
        ]
    )
    def test_getattr(self, cols, pcols):
        """Check whether columns accessible as dynamic attributes."""
        cb = ColumnBased(cols, pcols)
        for coldict in (cols, pcols):
            if coldict:
                colname = list(coldict.keys())[0]
                assert hasattr(cb, colname)

                colvals = getattr(cb, colname)
                assert np.all(colvals == coldict[colname])
                assert colvals is not coldict[colname]

    @pytest.mark.parametrize(
        "cols, pcols",
        [
            (mkcols(['rate'], [5]), None),    # no protected
            (None, mkcols(['time'], [5])),    # no normal
            (mkcols(['rate'], [5]), mkcols(['time'], [5])),    # expose
        ]
    )
    def test_columns_get(self, cols, pcols):
        """Check whether columns accessible as dynamic attributes."""
        cb = ColumnBased(cols, pcols)
        for coldict in (cols, pcols):
            if coldict:
                colname = list(coldict.keys())[0]
                assert (colname in cb.columns)

                colvals = cb.columns[colname]
                assert np.all(colvals == coldict[colname])
                assert (colvals is not coldict[colname])


    @pytest.mark.parametrize(
        "cols, pcols, excls",
        [
            (mkcols(['rate'], [5]), None, None),          # no protected
            (None, mkcols(['time'], [5]), ValueError),    # no normal
        ]
    )
    def test_columns_modify(self, cols, pcols, excls):
        """Check whether columns accessible as dynamic attributes."""
        cb = ColumnBased(cols, pcols)
        for coldict in (cols, pcols):
            if coldict:
                colname = list(coldict.keys())[0]
                if colname in cb.pcolnames:  # The column is protected
                    with pytest.raises(excls) as exception:
                        cb.columns[colname][0] = 55
                else:  # The column is normal
                    cb.columns[colname][0] = 55
                    res = coldict[colname].copy()  # Modify initial input value
                    res[0] = 55
                    assert np.all(cb.columns[colname] == res)


    @pytest.mark.parametrize(
        "cols, pcols, new, excls",
        [
            (
                mkcols(['rate'], [5]),
                None,
                array([1, 1, 1, 1, 1]),
                None
            ),  # everything is OK
            (
                mkcols(['rate'], [5]),
                None,
                array([1, 1, 1, 1]),
                ValueError
            ),  # wrong length
            (
                None,
                mkcols(['time'], [5]),
                array([1, 1, 1, 1, 1]),
                ValueError
            ),  # trying to modify protected
        ]
    )
    def test_columns_set(self, cols, pcols, new, excls):
        """Check whether columns accessible as dynamic attributes."""
        cb = ColumnBased(cols, pcols)
        for coldict in (cols, pcols):
            if coldict:
                colname = list(coldict.keys())[0]
                if excls:  # We're expecting an error
                    with pytest.raises(excls) as exception:
                        cb.columns[colname] = new
                else:
                    cb.columns[colname] = new
                    assert np.all(cb.columns[colname] == new)

