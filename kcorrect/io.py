# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Input and output utility.
"""
from __future__ import print_function, division
import os
import numpy as np
from .globaldefs import FTYPE


def path_to_file(fname, dirs):
    """Get the path to a file.

    The file is searched for in directories in order given in `dirs`.

    Parameters
    ----------
    fname : str
        File name to search for.
    dirs : list of directories
        List of directories to search in.

    Returns
    -------
    path : str
        The path to the file found or None if file is not found.
    """
    for path in ('/'.join([d, fname]) for d in dirs):
        if os.path.exists(path):
            return path
    return None


def read_basel(fname):
    """Read a spectrum from a Basel spectrum file.

    Parameters
    ----------
    fname : str
        Path to the Basel spectrum file.

    Notes
    -----
    See IDL procedure seds/k_read_basel.pro (kcorrect v4_2).
    """
    npoint = 1221  # hard-coded!!!
    
    f = open(fname)
    elems = f.read().split()
    lamb = np.empty(npoint, dtype=FTYPE)
    for i in range(npoint):
        lamb[i] = float(elems[i])

    flux, modelno, teff, logg, mh, vturb, xh = [], [], [], [], [], [], []
    ite = iter(elems[npoint:])
    try:
        while True:
            modelno.append( float(ite.next()) )
            teff.append( float(ite.next()) )
            logg.append( float(ite.next()) )
            mh.append( float(ite.next()) )
            vturb.append( float(ite.next()) )
            xh.append( float(ite.next()) )
            tmpflux = np.zeros(npoint, dtype=FTYPE)
            for i in range(npoint):
                tmpflux[i] = float(ite.next())
            flux.append(tmpflux)
    except StopIteration:
        pass
    return np.asarray(lamb, dtype=FTYPE), np.asarray(flux, dtype=FTYPE)
