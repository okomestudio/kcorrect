# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import print_function, division
import numpy as np
from .globaldefs import DEFAULT_TEMPLATE_DIR, FTYPE, ZRANGE_DEFAULT, VNAME
from .io import get_file_path
from .filter import FilterList
from .ext import k_read_ascii_table, k_projection_table


def load_vmatrix(vname=VNAME, file_paths=None):
    """Load spectral template information.

    Parameters
    ----------
    vname : str
        Name of templates.

    Returns
    -------
    A tuple of wavelengths and vmatrix.

    Notes
    -----
    This function loads so-called vmatrix in the IDL version.
    
    See IDL procedure utils/k_load_vmatrix.pro (kcorrect v4_2).
    """
    vfile = 'vmatrix.{:s}.dat'.format(vname)
    lfile = 'lambda.{:s}.dat'.format(vname)
    file_paths = (([] if file_paths is None else file_paths)
                  + [DEFAULT_TEMPLATE_DIR])
    vmatrix = k_read_ascii_table(get_file_path(vfile, file_paths))
    lamb = k_read_ascii_table(get_file_path(lfile, file_paths))
    return lamb, vmatrix


def make_projection_table(templates, filter_list, zrange=ZRANGE_DEFAULT):
    """Make a table of redshift-dependent projection.

    Returns a tuple of redshift values and corresponding projection
    table.

    Parameters
    ----------
    templates : (lamb_edges, flux)
        A spectrum to be convolved with filters.
    filter_list : a FilterList object
        Filters to convolve with.
    zrange : (z_min, z_max, nz)
        Tuple of z_min, z_max, and a number of redshifts between them.

    Notes
    -----
    This function generates so-called rmatrix in the IDL version.

    See IDL procedure fit/k_projection_table.pro (kcorrect v4_2).
    """
    zmin, zmax, nz = zrange
    zvals = zmin + (zmax - zmin) * (np.arange(nz) + 0.5) / nz
    filter_curves = filter_list.filters
    band_shift = 0.

    # k_projection_table.pro returns (zvals, rmatrix), so just
    # following the convention.
    return (zvals, k_projection_table(filter_curves,
                                      templates[0],
                                      templates[1],
                                      zvals,
                                      band_shift))


class ProjectionTable(tuple):
    """SED projection table.
    
    An instance of this class keeps the projection table for a set of
    filters and redshift specifications.
    """

    def __new__(cls, filter_list, zrange=ZRANGE_DEFAULT, vname=VNAME,
                file_paths=None):
        templates = load_vmatrix(vname, file_paths)
        pt = make_projection_table(templates, filter_list, zrange)
        return super(ProjectionTable, cls).__new__(cls, pt)

    def __init__(self, filter_list, zrange=ZRANGE_DEFAULT, vname=VNAME,
                 file_paths=None):
        self.filter_list = filter_list
        self.zrange = zrange
        self.vname = vname
        self.file_paths = file_paths
        self.templates = load_vmatrix(vname, file_paths)


class ProjectionTableDB(dict):
    """A collection of ProjectionTable instances."""
    
    def __init__(self, zrange=ZRANGE_DEFAULT, vname=VNAME, file_paths=None):
        self.zrange = zrange
        self.vname = vname
        self.file_paths = file_paths
        super(ProjectionTableDB, self).__init__({})

    def __getitem__(self, filter_list):
        if not isinstance(filter_list, FilterList):
            raise KeyError('Key must be a FilterList instance.')

        if filter_list in self.keys():
            return super(ProjectionTableDB, self).__getitem__(filter_list)

        # projection table does not exist yet, so make a new one
        self[filter_list] = ProjectionTable(filter_list,
                                            self.zrange,
                                            self.vname,
                                            self.file_paths)

        return super(ProjectionTableDB, self).__getitem__(filter_list)


# global master projection table database
PTABLE_MASTER = ProjectionTableDB()
