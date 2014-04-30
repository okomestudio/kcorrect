# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals

from distutils.core import Extension
import os.path
from os.path import join
import sys

from astropy import setup_helpers
from astropy.extern import six


KCROOT = os.path.relpath(os.path.dirname(__file__))


if six.PY3:
    def string_escape(s):
        s = s.decode('ascii').encode('ascii', 'backslashreplace')
        s = s.replace(b'\n', b'\\n')
        s = s.replace(b'\0', b'\\0')
        return s.decode('ascii')
elif six.PY2:
    def string_escape(s):
        # string_escape has subtle differences with the escaping done in Python
        # 3 so correct for those too
        s = s.encode('string_escape')
        s = s.replace(r'\x00', r'\0')
        return s.replace(r"\'", "'")


def get_extensions():
    from astropy.version import debug

    ######################################################################
    # DISTUTILS SETUP
    cfg = setup_helpers.DistutilsExtensionArgs()

    cfg['include_dirs'].extend(['numpy', join(KCROOT, "include")])
    #cfg['define_macros'].extend([
    #    ('ECHO', None),
    #    ('WCSTRIG_MACRO', None),
    #    ('ASTROPY_WCS_BUILD', None),
    #    ('_GNU_SOURCE', None),
    #    ('WCSVERSION', WCSVERSION)])

    ## if (not setup_helpers.use_system_library('wcslib') or
    ##     sys.platform == 'win32'):
    ##     write_wcsconfig_h()

    ##     wcslib_path = join("cextern", "wcslib")  # Path to wcslib
    ##     wcslib_cpath = join(wcslib_path, "C")  # Path to wcslib source files
    ##     wcslib_files = [  # List of wcslib files to compile
    ##         'flexed/wcsbth.c',
    ##         'flexed/wcspih.c',
    ##         'flexed/wcsulex.c',
    ##         'flexed/wcsutrn.c',
    ##         'cel.c',
    ##         'lin.c',
    ##         'log.c',
    ##         'prj.c',
    ##         'spc.c',
    ##         'sph.c',
    ##         'spx.c',
    ##         'tab.c',
    ##         'wcs.c',
    ##         'wcserr.c',
    ##         'wcsfix.c',
    ##         'wcshdr.c',
    ##         'wcsprintf.c',
    ##         'wcsunits.c',
    ##         'wcsutil.c']
    ##     cfg['sources'].extend(join(wcslib_cpath, x) for x in wcslib_files)
    ##     cfg['include_dirs'].append(wcslib_cpath)
    ## else:
    ##     cfg.update(setup_helpers.pkg_config(['wcslib'], ['wcs']))

    #cfg.update(setup_helpers.pkg_config(['kcorrect'], ['kcorrect']))

    astropy_kc_files = ['gaussj.c',
                        'iterate_lf.c',
                        'k_binspec.c',
                        'k_brent.c',
                        'k_choldc.c',
                        'k_cholsl.c',
                        'k_evolve.c',
                        'k_fileopen.c',
                        'k_filter_struct.c',
                        'k_fit_nonneg.c',
                        'k_fit_photoz.c',
                        'k_fit_spec.c',
                        'k_fit_spec_linear.c',
                        'k_interpolate.c',
                        'k_load_filters.c',
                        'k_locate.c',
                        'k_midpnt.c',
                        'k_nonneg_solve.c',
                        'k_polint.c',
                        'k_projection_table.c',
                        'k_qromo.c',
                        'k_read_ascii_table.c',
                        'k_reconstruct_maggies.c',
                        'k_strparse.c',
                        'k_utils.c',
                        'k_yanny_readone.c',
                        'k_zbrent.c',
                        'lf_WH_interp.c',
                        'lf_calc_vmax.c',
                        'lf_eep.c',
                        'lf_eepfit.c',
                        'lf_select_eep.c',
                        'lf_set_AB.c',
                        'lf_sum_AB.c',
                        'phierrors_lf.c',
                        'philike.c',
                        'ztransform.c']                        

    cfg['sources'].extend(join(KCROOT, 'src', x) for x in astropy_kc_files)

    ## if debug:
    ##     cfg['define_macros'].append(('DEBUG', None))
    ##     cfg['undef_macros'].append('NDEBUG')
    ##     if (not sys.platform.startswith('sun') and
    ##         not sys.platform == 'win32'):
    ##         cfg['extra_compile_args'].extend(["-fno-inline", "-O0", "-g"])
    ## else:
    ##     # Define ECHO as nothing to prevent spurious newlines from
    ##     # printing within the libwcs parser
    ##     cfg['define_macros'].append(('NDEBUG', None))
    ##     cfg['undef_macros'].append('DEBUG')

    ## if sys.platform == 'win32':
    ##     # These are written into wcsconfig.h, but that file is not
    ##     # used by all parts of wcslib.
    ##     cfg['define_macros'].extend([
    ##         ('YY_NO_UNISTD_H', None),
    ##         ('_CRT_SECURE_NO_WARNINGS', None),
    ##         ('_NO_OLDNAMES', None),  # for mingw32
    ##         ('NO_OLDNAMES', None),  # for mingw64
    ##         ('__STDC__', None)  # for MSVC
    ##     ])

    ## if sys.platform.startswith('linux'):
    ##     cfg['define_macros'].append(('HAVE_SINCOS', None))

    cfg['sources'] = [str(x) for x in cfg['sources']]

    cfg = dict((str(key), val) for key, val in six.iteritems(cfg))

    return [Extension(str('kcorrect._clib'), **cfg)]


def get_package_data():
    # Installs the testing data files
    api_files = ['kcorrect.h']
    api_files = [join('include', x) for x in api_files]

    return {str('kcorrect'): api_files}
