"""General purpose functions."""

from astropy.io import fits
from mypythonlib import Actions
from .extpath import ExtPath


def fits_keywords_getmany(HDU, keynames: list[str], defaults: dict = None, *,
                          action=Actions.EXCEPTION) -> dict[str]:
    """Read multiple keyword from the FITS header.

    The aim of this function to provide a functional analogues to dict.get()
    for headers' keywords.

    :param HDU: One object from the list obtained via astropy.io.fits.open().
    :param keynames: List of the keywords' names to read from FITS header
    :param defaults: If the keyword is missing use value from this dict.
        The default is None.
    :param action: Perform 'action' if the keyword is missing both in the
        header and in the defaults. The default is EXCEPTION
    :return: dict with keywords and their values
    """
    result={}
    for key in keynames:
        if key in HDU.header:  # Search for KEY or a KEYI+KEYF pair
            result[key] = HDU.header[key]
        elif (key+'I' in HDU.header) and (key+'F' in HDU.header):
            result[key] = HDU.header[key+'I'] + HDU.header[key+'F']
        else:   # Try to get it from user dict
            if key not in (defaults or {}):
                Actions(action).do(f"Required keyword '{key}' is not found",
                                   exception_class=KeyError)
            else:
                result[key] = defaults[key]
    return result


def gti_get_limits(gtifile: ExtPath):
    """Return boundaries of the GTI interval.

    :param gtifile: Path to the GTI FITS file
    :returns: a pair of float values
    """
    with fits.open(gtifile) as ftsgti:
        start=ftsgti[gtifile.hdu].data['START']
        stop=ftsgti[gtifile.hdu].data['STOP']
    return [start[0], stop[-1]]     # Get starttime from the first row and stoptime from the last row


