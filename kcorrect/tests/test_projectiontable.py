# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import print_function, division
from kcorrect.filter import FilterList
from kcorrect.projectiontable import ProjectionTable, ProjectionTableDB
import pytest


def test_general():
    filter_list1 = FilterList(['sdss_u0', 'sdss_g0', 'sdss_r0', 'sdss_i0',
                              'sdss_z0'])
    filter_list2 = FilterList(['sdss_u0', 'sdss_g0', 'sdss_r0', 'sdss_i0'])

    ptdb = ProjectionTableDB()

    a = ptdb[filter_list1]
    b = ptdb[filter_list1]
    c = ptdb[filter_list2]

    assert a is b
    assert a is not c
