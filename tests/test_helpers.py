import pytest
import numpy as np
from numpy import array
from astroscripts._internal.helpers import _make_vector


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
            # (iter([1, 2, 3]), array([1, 2, 3]), None, False),
            (None, None, True, False, None),
            (2, array([2, 2, 2]), False, True, 3),
        ]
    )
    def test_valid_input(self, input, expected, noneallowed, from_numbers, reflen):
        output = _make_vector(input, none_is_allowed=noneallowed,
                              from_numbers=from_numbers, reflen=reflen)
        assert type(output) == type(expected)
        assert np.all(output == expected)

    @pytest.mark.parametrize(
        "input,noneallowed,excls, reflen",
        [
            (None, False, TypeError, None),
            ('qwe', True, TypeError, None),
            ([], False, ValueError, None),
            ([], True, ValueError, None),
            ([1,2,3], False, ValueError, 4),
            (array([[1,2,3],[4,5,6]]), False, ValueError, None),
            ([[1,2,3],[4,5,6]], False, ValueError, None),
        ]
    )
    def test_invalid_input(self, input, noneallowed, excls, reflen):
        with pytest.raises(excls) as exception:
            _make_vector(input, vecname='values', none_is_allowed=noneallowed,
                         reflen=reflen)

