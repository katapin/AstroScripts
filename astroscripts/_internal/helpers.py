
import numpy as np
from numpy import ndarray
from typing import Iterable, Sized


def _make_vector(vec: Iterable | float | None, vecname: str = None, *,
                 none_is_allowed: bool = False, from_numbers: bool = False,
                 reflen: int | None = None) -> ndarray | None:
    """Create new ndarray vector.

    Takes an iterable object or a number, returns a ndarray with own data.
    """
    _vecname = f"'{vecname}' " if vecname else ''
    _typertext = f"Data column {_vecname}must be a sequence (list, tuple, etc.) of numbers "\
                 "or a numpy array"
    if vec is None:
        if none_is_allowed is True:
            return None
        else:
            raise TypeError(f"The vector {_vecname}cannot be None.")

    if isinstance(vec, Iterable):
        vec2 = vec
        if not isinstance(vec, Sized):   # numpy can't create arrays from iterator
            vec2 = tuple(vec)

        if len(vec2) == 0:
            raise ValueError(
                "The vector {}{}cannot be empty.".format(
                    _vecname,
                    'can be None but ' if none_is_allowed else ''
                )
            )

        res = np.array(vec2)  # create the array to return

        if (res.ndim > 1 or    # the original vec was multidimensional (nested)
                (reflen and reflen != len(res))):   # different length
            raise ValueError("All the arrays must be 1D of the the same length.")
        if res.ndim == 1:
            return res

    else:  # probably, this is a number
        if from_numbers and reflen:
            try:
                flt = float(vec)   # test for a number
                res = np.ones((1, reflen)) * flt  # Create 1D array
                return res
            except Exception:
                raise TypeError(_typertext + ' or a single number.')

    # vec was a string or something not digestible by numpy
    raise TypeError(_typertext + '.')
