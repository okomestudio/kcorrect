# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Global definitions.
"""
from __future__ import print_function, division
import os
from numpy import float32


KCORRECT_DIR = os.path.relpath(os.path.dirname(__file__))
if KCORRECT_DIR is None:
    raise ValueError('The environment variable KCORRECT_DIR not defined.')
elif not os.path.exists(KCORRECT_DIR):
    raise IOError('Problem accessing {0}.'.format(KCORRECT_DIR))


# default directories for some data sets
DEFAULT_FILTER_DIR = KCORRECT_DIR + '/data/filters'
DEFAULT_TEMPLATE_DIR = KCORRECT_DIR + '/data/templates'

# default cosmology (Omega_matter, Omega_lambda, h_100)
COSMO_DEFAULT = (0.3, 0.7, 0.7)

# default redshift range (zmin, zmax, # of points)
ZRANGE_DEFAULT = (0.0, 2.0, 1000)

# template name
VNAME = 'default'

# within the kcorrect C library, float C type is consistently used.
FTYPE = float32

# filter definition file extension
FILTEREXT = '.par'
