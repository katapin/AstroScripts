"""Extended paths: provide paths with HDU and table columns specifications."""

import os
import re
import copy
from pathlib import PurePath
from typing import Union, Self
from mypythonlib import FilePath, FilePathAbs, printwarn

__all__ = [
    "ExtFileName",
    "ExtPath",
    "ExtPathAbs"
]


class ExtFileName:
    """Extended filename to store also hdu number and columns of FITS/ASCII tables.

    Notes. The PurePath's constructor calls str(os.fspath(obj)) for each
    of the input arguments and then uses the resultant string with '/'.join().
    The function os.fspath() do nothing for subclasses of str and calls
    __fspath__() otherwise. Therefore, if ExtFileName is a subclass of str,
    it must return the pure filename (i.e. without hdu or filter) via
    its __str__ method. If it is of any other type, the filepath must be
    returned by __fspath__, and the magic __str__ may have any output.
    """

    def __init__(self, name: str, hdu: int | str = None, filter=None):
        self._name = self._name_check(name)
        self.hdu = hdu
        self.filter = filter

    @classmethod
    def _name_check(cls, name):
        if not isinstance(name, str) or \
                FilePath(name).name != name:
            raise ValueError(f"Incorrect name '{repr(name)}'. Only pure "
                             "string names are allowed (not paths or other objects).")
        if '[' in name or ']' in name:
            printwarn(f"you are using square brackets in name: {name}.",
                         cls.__name__)
        return name

    def __str__(self):
        lst = [self._name]
        if self._hdu:
            lst.append("[{!r}]".format(self._hdu))
        if self.filter:
            lst.append("[{!r}]".format(self.filter))
        return ''.join(lst)

    def __repr__(self):
        return '<{}({!r}, hdu={!r}, filter={!r})>'.format(
            self.__class__.__name__, self._name, self._hdu, self.filter)

    def __fspath__(self):
        """Implement os.PathLike protocol."""
        return self._name

    def __rtruediv__(self, other):
        if isinstance(other, (str, os.PathLike)):
            return ExtPath(other, self)
        return NotImplemented

    def get_opts(self):
        """Return hdu and filter arguments as a dict."""
        return dict(hdu=self._hdu, filter=self.filter)

    @property
    def name(self):
        return self._name

    @property
    def hdu(self):
        return self._hdu

    @hdu.setter
    def hdu(self, val):
        if not isinstance(val, (int, str, type(None))):
            raise TypeError("HUD key must be either a integer or a string")
        self._hdu = val

    # @property
    # def filter(self):
    #     return self._filter

    # @filter.setter
    # def filter(self, val):
    #     self._filter = val

    @classmethod
    def from_string(cls, expression: str, without_hdu=False):
        """Create the object from the parse_string() result."""
        return cls(*cls.parse_string(expression, without_hdu))

    @staticmethod
    def _process_filter(match1, without_hdu):
        if without_hdu is True:
            filter = match1.group(2)
            if match1.group(3):
                raise ValueError(f"Inadmissible expression {match1.group(3)} "
                                 "in the second pair of brackets.")
        else:
            filter = match1.group(3)

        if filter and len(filter) > 2 and filter[0] == filter[-1] and \
                (filter[0] == '"' or filter[0] == "'"):
            filter = filter[1:-1]

        return filter

    @classmethod
    def parse_string(cls, expression: str, without_hdu: bool = False,
                     only_path: bool = False):
        """Parse if expression has a format "file.fts[hdu]['filter']".

        Here hdu may be an integer representing an ordering number of the HDU in
        the FITS file or a string representing the HDU's name. The existence
        of the HDU won't be checked. Path and filter can be any strings.

        :param expression: String expression to analyse.
        :param without_hdu: If True it allows to omit the hdu key and
            consider the string in the first square brackets as
            the filtering exporession.
        :param only_path: Ignore the hdu and filter keys and return
            only the filepath. It may be useful for isolation of
            the pure filepath from the total expression without
            printing warnings.
        :return: The tuple of three elements: tuple(path, hdu, filter),
            where the path van be empty string while hdu and filter
            can be None.
        :raises ValueError:when the hdu string cannot be parsed.
        """
        hdu, filter = None, None
        match1 = re.match("(.*?)(?:\[(.*?)])?(?:\[(.*?)])?$", expression)  # It always gives some result,
        # at least match1.group(1) always exists but may be not a filename
        if only_path:
            return match1.group(1), hdu, filter

        if without_hdu is False:
            if match1.group(2):  # Parse first brackets
                hdustr = match1.group(2).replace('"', '').replace("'", '')  # remove all quotes
                if match2 := re.match("""^(\d+)|(\w+)$""", hdustr):  # int or string
                    if match2.group(1):
                        hdu = int(match2.group(1))
                    elif match2.group(2):
                        hdu = str(match2.group(2))
                else:  # Can't parse hdu key
                    without_hdu = True
                    printwarn(f"the string {match1.group(2)} will be "
                                 "interpreted as the filtering expression.", cls.__name__)

        filter = cls._process_filter(match1, without_hdu)
        return match1.group(1), hdu, filter


class ExtPath(FilePath):
    """Extend FilePath to store FITS's hdu and filtering expression

    The object can be constructed from path, path + ExtName, string or
    other ExtPath, where path is any os.Pathlike object; HDU and
     filter can be passed as keyword arguments."""

    def __new__(cls, *args: Union[str, PurePath, ExtFileName], **kwargs):
        superargs = list(args) if len(args) > 0 else ['']
        last = superargs[-1]  # the last argument
        if isinstance(last, ExtFileName):
            extname = ExtFileName(last.name, **kwargs) if kwargs else copy.copy(last)
        elif isinstance(last, PurePath):
            extname = copy.copy(last._extname) if hasattr(last, '_extname') else \
                ExtFileName(last.name, **kwargs)
        elif isinstance(last, str):  # Try to parse string
            # Split to the pure string path and its extended part
            strpath, hdu, filter = ExtFileName.parse_string(last, only_path=bool(kwargs))
            extname_extra_args = kwargs or {'hdu': hdu, 'filter': filter}
            extname = ExtFileName(FilePath(strpath).name, **extname_extra_args)
            superargs[-1] = strpath  # save the string without the tail [hdu][filter]
        else:
            raise TypeError(f"Unsupported type of argument: {type(last)}.")
        obj = super().__new__(cls, *superargs)
        obj._extname = extname
        return obj

    def _from_parsed_parts(self, drv, root, parts):
        """Construct object from parts.

        This function overrides the PurePath's classmethod(!) used to
        create derivative subpaths from the original path (for self.with_name,
        self.parent, etc.). Now, it will be an ordinary method which gets extname
        from 'self' and puts it into the child object.
        """
        ancestor = {ExtPath: FilePath, ExtPathAbs: FilePathAbs}
        cls = self.__class__
        name = os.fspath(parts[-1]) if len(parts) > 0 else ''
        #   for self.parent                   for self.with_parent
        if len(parts) != len(self._parts) and name != self.name:
            cls = ancestor[cls]
        obj = FilePath._from_parsed_parts.__func__(cls, drv, root, parts)
        obj._extname = self._extname
        return obj

    def __fspath__(self):
        """Return the normal path for os.fspath()."""
        return super().__str__()

    def __str__(self):
        """Return string representation."""
        tail = str(self._extname)[len(self._extname.name):]
        return super().__str__() + tail

    def __repr__(self):
        return "{}({!r}/{})".format(self.__class__.__name__,
                                    str(self.parent), repr(self._extname))

    def __truediv__(self, other):
        """Return self / other."""
        raise TypeError(f"Unsupported operation for {self.__class__}")

    def __rtruediv__(self, other):
        """Return other / self."""
        if isinstance(other, (str, os.PathLike)):
            return self.__class__(other, self)

    @property
    def extname(self):
        """Return name as a ExtFileName object."""
        return self._extname

    @property
    def hdu(self):
        return self._extname.hdu

    def with_extname(self: Self, new_extname: ExtFileName) -> Self:
        return self.__class__(self.parent, new_extname)


class ExtPathAbs(ExtPath, FilePathAbs):
    pass


ExtPath._class_abspath = ExtPathAbs
ExtPathAbs._class_abspath = ExtPathAbs

