"""Provide executables for the package."""

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


def init():
    parser = argparse.ArgumentParser(
        description="Initialize collections of scripts")
    parser.add_argument(
        'name', nargs='*', help="name of collection")
    parser.add_argument(
        '-v', '--version', action='store_true', help="Print package version")
    argnspace = parser.parse_args(sys.argv[1:])

    initlines = [
        'source ~/.bashrc\n',
        'PS1="\033[1;34m${PS1}\033[0m"\n'
        'heainit\n'
    ]

    pkgname = getmainname(__package__)   # Name of the main package
    pathlines = []

    # Print version info
    if argnspace.version:
        print('{}, version {}'.format(pkgname, importlib.metadata.version(pkgname)))
        exit(0)

    try:
        rootpath = os.path.dirname(__file__) + '/' + SCRIPTROOTDIR
        allnames = os.listdir(rootpath)
    except:
        die("Can't find folder with scripts. The package might be corrupted.", progname=pkgname)

    # Process 'all' option
    if len(argnspace.name) == 0 or 'all' in argnspace.name:
        names = allnames
    else:
        names = argnspace.name

    for folder in names:
        if folder not in allnames:
            die(f"Script collection '{folder}' is not found.", progname=pkgname)
        pathlines.append(rootpath + '/' + folder)
        initlines.append(AUXILINIT.get(folder,'')+'\n')

    with NamedTemporaryFile(mode='a') as initsh:
        initlines.append('export PATH={}:$PATH\n'.format(
            ':'.join(pathlines)))
        initsh.writelines(initlines)
        initsh.flush()
        os.system(f'bash --init-file {initsh.name}')
