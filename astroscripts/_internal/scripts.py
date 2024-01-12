"""Provide executables for the package (described in pyproject.toml's [project.scripts] section)."""

import os
import sys
import argparse
import importlib.metadata
from tempfile import NamedTemporaryFile

from mypythonlib import die
from mypythonlib.packaging import getmainname

SCRIPTROOTDIR = 'libexec'      # folder with scripts inside this package
AUXILINIT = {                  # additional init command for the specific collection
    'xmm': 'sasinit',
    'chandra': 'ciaoinit',
}

def _get_activate_path():
    """Return path to the venv's activate script."""
    bindir = os.path.dirname(sys.executable)  # Dir with python interpreter
    actpath = bindir + '/activate'
    return actpath if os.path.exists(actpath) else None


def init():
    """Active venv and set $PATH variable.

    This functions starts a child bash process with activated python
    virtual environment and modified $PATH variable. To do this, it
    creates own temporary init bash script.
    """
    parser = argparse.ArgumentParser(
        description="Initialize specific collection of scripts")
    parser.add_argument(
        'name', nargs='*', help="name of collection")
    parser.add_argument(
        '-v', '--version', action='store_true', help="Print package version")
    parser.add_argument(
        '-x', '--auxil', action='store_true', help="Additionally call the collection "
        "specific init script.")
    argnspace = parser.parse_args(sys.argv[1:])

    pkgname = getmainname(__package__)   # Name of the main package
    pathlines = []  # Store paths to the claimed collections to append to $PATH

    # Print version info
    if argnspace.version:
        print('{}, version {}'.format(pkgname, importlib.metadata.version(pkgname)))
        exit(0)

    # Probe the folder that we intend to make accessible
    try:
        rootpath = os.path.dirname(__file__) + '/' + SCRIPTROOTDIR
        allnames = os.listdir(rootpath)
    except:
        die("Can't find folder with scripts. The package might be corrupted.", progname=pkgname)

    # Process 'all' option
    if argnspace.name == 'all' in argnspace.name:
        names = allnames
    else:
        names = argnspace.name

    initshlines = [             # Content of the init script
        'source ~/.bashrc\n',
        '_PS1_bak=${PS1}\n']
    if actpath := _get_activate_path():
        initshlines.append(f'source {actpath}\n')  # add 'activate' scripts if it exists
    initshlines += [
        'PS1="\033[1;34m${_PS1_bak}\033[0m"\n',
        'heainit\n'  # Init HEASOFT
    ]

    for folder in names:
        if folder not in allnames:
            die(f"Script collection '{folder}' is not found.", progname=pkgname)
        pathlines.append(rootpath + '/' + folder)
        if argnspace.auxil:  # Add collection-specific init script
            initshlines.append(AUXILINIT.get(folder, '')+'\n')


    with NamedTemporaryFile(mode='a') as initsh:
        initshlines.append('export PATH={}:$PATH\n'.format(
            ':'.join(pathlines)))
        initsh.writelines(initshlines)
        initsh.flush()
        os.system(f'bash --init-file {initsh.name}')
