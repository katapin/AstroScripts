
from copy import deepcopy
from abc import ABC, abstractmethod
from itertools import chain
from collections.abc import KeysView
from functools import singledispatchmethod
from typing import Iterable, Sized, Self
# from mypythonlib.collections import ROdict

import numpy as np
from numpy.typing import NDArray, ArrayLike


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

        #C an the vector by empty?
        # if len(vec2) == 0:
        #     raise ValueError(
        #         "The vector {}{}cannot be empty.".format(
        #             _vecname,
        #             'can be None but ' if none_is_allowed else ''
        #         )
        #     )

        res = np.array(vec2)  # create the array to return

        if (res.ndim > 1 or    # the original vec was multidimensional (nested)
                (reflen and reflen != len(res))):   # different length
            raise ValueError("All the arrays must be 1D of the the same length but "
                             f"{_vecname}doesn't.")
        if res.ndim == 1:
            return res

    elif from_numbers:  # probably, this is a number
        if not isinstance(reflen, int):
            raise TypeError(f"Integer 'reflen' is expected but {reflen} is received.")
        try:
            flt = float(vec)   # test for a number
            res = np.ones((1, reflen)) * flt  # Create 1D array
            return res
        except Exception:
            raise TypeError(_typertext + ' or a single number.')

    # vec was a string or something not digestible by numpy
    raise TypeError(_typertext + f' (not {type(vec)}).')


# class ColumnStorage(ROdict):
#     """Read-only dict to store ndarray vectors."""
#
#     def _item_change_existing_as_normal(self, key, val):
#         self._dict[key][:] = val
#
#     def _item_get_missing_as_normal(self, key):
#         raise KeyError(f"Column '{key}' doesn't exist.")
#
#     def _item_del_missing_as_normal(self, key):
#         raise KeyError(f"Column '{key}' doesn't exist.")
#
#     def __str__(self):
#         """Return string representation."""
#         text = self.__class__.__name__ + ': {\n'
#         lst = []
#         for k, v in self._dict.items():
#             lst.append((f"'{k}'*:" if k in self._readonly else f"'{k}':") + \
#                        repr(v))
#         return text + '\n'.join(lst) + ' }'


# class _ColDesc():
#     """Descriptor"""
#
#
#
#     def __set_name__(self, owner, name):
#         self.name = name
#
#     def __get__(self, obj, objtype=None):
#         print(obj)
#         if not obj:
#             return
#         if obj._exposecols is False:
#             raise AttributeError(f"There is no attribute 'self.name'")
#
#         return self.Proxy(obj)


class ColumnBasedMinimal(ABC):
    """Abstract class for column-based objects like light curves,
    spectra, plots, etc. Provide sampling via masks and
    slices. It assumes that all the passed columns are
    valid ndarrays of the same length (it doesn't check this).
    """

    @abstractmethod
    def __init__(self, columns: dict[str, NDArray] | None = None):
        self._columns = columns.copy() if columns else {}

    def __len__(self):
        """Return length of the columns."""
        if len(self._columns):
            fk = tuple(self._columns.keys())[0]  # first key
            return len(self._columns[fk])
        else:  # There are no columns
            return 0

    def __eq__(self, other: Self) -> bool:
        if type(self) is not type(other):
            return NotImplemented
        if self._columns.keys() != other._columns.keys():
            return False
        for colname in self._columns:
            c1, c2 = self._columns[colname], other._columns[colname]
            if len(c1) != len(c2):
                return False
            if not np.all(np.isclose(c1, c2)):
                return False
        return True

    @singledispatchmethod
    def __getitem__(self, arg):
        raise TypeError("Unsupported type of key: {}.".format(type(arg)))

    @__getitem__.register(int)
    def _(self, arg: int):
        return self._getitem_point(arg)

    @__getitem__.register
    def _(self, arg: slice | list | np.ndarray):
        return self._getitem_subsample(arg)

    def _getitem_point(self, index: int) -> dict[str, float]:
        """Get a single point (in multidimensional space)."""
        return {n: float(ar[index]) for n, ar in self._columns.items()}

    def _getitem_subsample(self, arg: slice | ArrayLike) -> Self:
        """Get a subsample via slice or mask."""
        obj = self.copy()
        for colname in obj._columns:
            obj._columns[colname] = obj._columns[colname][arg]
        if hasattr(obj, '_rebuild_hook'):  # The object may require additional tuning
            obj._rebuild_hook()
        return obj

    def copy(self) -> Self:
        """Return a deep copy of the object."""
        # It's fantastic but self.__proxy.owner = self is treated correctly!
        return deepcopy(self)

class ColumnBased(ColumnBasedMinimal):
    """Extends ColumnBasedMinimal. Provide access to columns
    through dynamic attributes of the same name and through
    dedicated self.columns object. The latter method returns
    a read-only view of the ndarray for the protected columns
    and a normal view for the remaining. The attribute access
    always returns a read-only view.
    """

    class __ColumnProxy:
        """The proxy object to be returned via self.columns,
        provides dict-like access to the columns."""
        def __init__(self, owner):
            self.owner = owner

        def __contains__(self, key: str) -> bool:
            return key in self.owner._columns

        def __get_column_raw(self, key: str) -> NDArray:
            if key not in self.owner._columns:
                raise KeyError(f"Column '{key}' does not exist")
            return self.owner._columns[key]

        def __getitem__(self, key: str) -> NDArray:
            view = self.__get_column_raw(key).view()
            if key in self.owner._pcolnames:        # If protected
                view.flags.writeable = False        # make the view read-only
            return view

        def __setitem__(self, key: str, value: NDArray):
            if key in self.owner._pcolnames:
                # We use ValueError instead of TypeError to match
                # the numpy behaviour when writing to the writeable=False
                # destination array
                raise ValueError(f"assignment ('{key}') destination is read-only")
            self.__get_column_raw(key)[:] = value   # numpy broadcasting guaranties the same shape

        def __repr__(self) -> str:
            return 'Columns([{}])'.format(', '.join(
                f"'{col}'" for col in self.owner._columns.keys())
            )

    @abstractmethod
    def __init__(self, columns: dict[str, NDArray] | None = None,
                 pcolumns: dict[str, NDArray] | None = None):
        """
        :param columns: dict of ordinary columns
        :param pcolumns: dict of protected columns
        """
        super().__init__(columns)
        self._pcolnames = set()        # Names of the protected columns
        if pcolumns:
            self._columns.update(pcolumns)
            self._pcolnames.update(pcolumns.keys())

        self.__proxy = self.__ColumnProxy(self)    # Dict-like access

    def __eq__(self, other: Self):
        return (super().__eq__(other) and
            self._pcolnames == other._pcolnames)

    def __dir__(self):
        """Provide autocomplete for dynamic attributes in ipython."""
        return chain(object.__dir__(self), self._columns.keys())

    def __getattr__(self, item):
        try:
            col = self.columns[item]
            # Make every column read-only independently on its protection status
            col.flags.writeable = False
        except KeyError:
            raise AttributeError(f"There is no attribute '{item}'")
        return col

    @property
    def colnames(self) -> list:
        """Return names of all the columns."""
        return list(self._columns.keys())

    @property
    def pcolnames(self) -> list:
        """Return names of the protected columns."""
        return list(self._pcolnames)

    @property
    def columns(self):
        return self.__proxy
