
import numpy as np
from numpy.typing import NDArray, ArrayLike
from typing import Iterable, Sized
from abc import ABC, abstractmethod
from functools import singledispatchmethod
from mypythonlib.collections import ROdict


def _make_vector(vec: ArrayLike | None, vecname: str = None, *,
                 none_is_allowed: bool = False, from_numbers: bool = False,
                 reflen: int | None = None) -> NDArray[float] | None:
    """Create new ndarray vector.

    Takes an iterable object or a number, returns an ndarray with own data.
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


class ColumnStorage(ROdict):
    """Read-only dict to store ndarray vectors."""

    def _item_change_existing_as_normal(self, key, val):
        self._dict[key][:] = val

    def _item_get_missing_as_normal(self, key):
        raise KeyError(f"Column '{key}' doesn't exist.")

    def _item_del_missing_as_normal(self, key):
        raise KeyError(f"Column '{key}' doesn't exist.")

    def __str__(self):
        """Return string representation."""
        text = self.__class__.__name__ + ': {\n'
        lst = []
        for k, v in self._dict.items():
            lst.append((f"'{k}'*:" if k in self._readonly else f"'{k}':") + \
                       repr(v))
        return text + '\n'.join(lst) + ' }'


class ColumnBased(ABC):
    """Abstract class for column-based objects like light curves,
     spectra, plots, etc. Provide sampling via masks and
     slices. Assumes that the data is stored in _columns and
     _meta attributes."""

    @abstractmethod
    def __init__(self, columns: ColumnStorage, meta: ROdict):
        self._columns = columns
        self._meta = meta

    @singledispatchmethod
    def __getitem__(self, key):
        raise TypeError("Unsupported type of key: {}.".format(type(key)))

    @__getitem__.register
    def _(self, index: int):
        """Must be implemented by the subclass."""
        return NotImplemented

    @__getitem__.register
    def __getitem__(self, arg: slice, list, ndarray):
        obj = object.__new__(self.__class__)
        obj._meta = self._meta
        # for col in self._columns
        for vecname in ['main', 'lower', 'upper']:
            obj._values[vecname] = obj._values[vecname][arg]
        return obj

    import copy

