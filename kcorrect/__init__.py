# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
kcorrect -- Python port

This is an astropy affiliated package containing a Python port of 
the K-correction program written by Michael Blanton (NYU).  
Currently this only implements the most basic aspects of 
K-correction.
"""

try:
    from .version import version as __version__
except ImportError:
    # TODO: Issue a warning...
    __version__ = ''
# The version number can be found in the "version" variable of version.py

# set up the test command
from astropy.tests.helper import TestRunner
_test_runner = TestRunner(__path__[0])
del TestRunner
test = _test_runner.run_tests
