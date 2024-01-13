from mypythonlib import printwarn


class LCError(Exception):
    """Class for errors arising during manipulation with light curves."""
    pass


def _bkgration_msg(bkgratio) -> None:
    """Print warning if bkgratio < 1."""
    if bkgratio < 1:
        printwarn("The provided bkgratio=S(bkg)/S(src) < 1. Please, make sure that "
                     "everything is correct.")