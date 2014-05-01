# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import print_function, division
import numpy as np
from .filter import FilterList
from .globaldefs import COSMO_DEFAULT, FTYPE, KCORRECT_DIR, ZRANGE_DEFAULT
from .io import read_basel

from .projectiontable import PTABLE_MASTER
from .utils.cosmology import ztodm
#from .clib import k_fit_nonneg, k_fit_photoz, k_reconstruct_maggies
from ext import k_reconstruct_maggies
from kcorrect.utils.spec import project_filters, lamb_to_edges


## def add_absigma_to_maggie_ivar(maggie, maggie_ivar, absigma):
##     """Add AB magnitude uncertainty to inverse variance in maggies.

##     This function adds in quadrature AB magnitude uncertainty to
##     inverse variances and returns the inverse variances in maggies.

##     Parameters
##     ----------
##     maggie, maggie_ivar : array_like
##         AB maggies and their inverse variance.
##     absigma : array_like
##         1-sigma uncertainties in AB magnitudes.

##     Notes
##     -----
##     See utils/k_minerror.pro of kcorrect v4_2.
##     """
##     min_err = np.ma.array(absigma, copy=False, dtype=FTYPE)
##     maggie = np.ma.array(maggie, copy=False, dtype=FTYPE)
##     maggie_ivar = np.ma.array(maggie_ivar, copy=False, dtype=FTYPE)
##     maggie_ivar = np.ma.masked_equal(maggie_ivar, 0.)
##     factor = 2.5 / np.ma.log(10.)
##     err = factor / np.ma.sqrt(maggie_ivar) / maggie
##     err2 = err**2 + min_err**2
##     maggie_ivar = factor * factor / (maggie * maggie * err2)
##     return maggie_ivar.filled(0.)


## def vega2ab(filter_list, ref='kurucz'):
##     """Compute the magnitude offsets between Vega and AB magnitudes.

##     This function computes conversion offset from Vega to AB
##     magnitudes such that

##     m(AB) = m(Vega) + vega2ab

##     where vega2ab is the offset provided by this function.

##     Parameters
##     ----------
##     filter_list : a FilterList object
##         List of filters for photometry.
##     ref : {'kurucz', 'hayes'}, optional
##         Stellar atmosphere model of reference source.  The default is
##         'kurucz'.

##     Notes
##     -----
##     See IDL procedure utils/k_vega2ab.pro (kcorrect v4_2).
##     """
##     if ref == 'kurucz':
##         veganame = 'lcbvega.ori'
##         fname = ''.join([KCORRECT_DIR, '/data/basel/', veganame])
##         lamb, flux = read_basel(fname)
##         lamb = lamb * 10.
##         cspeed = 2.99792e+18
##         flux = np.pi * 4. * flux * cspeed / lamb**2
##         flux = flux * 6.5043898e-17
##     elif ref == 'hayes':
##         fname = ''.join([KCORRECT_DIR, '/data/filters/',
##                          'hoggraw/hayes/hayes.txt'])
##         with open(fname) as f:
##             lamb, lflux = [], []
##             for each in f:
##                 if each.startswith('#'):
##                     continue
##                 x, y = each.split()
##                 lamb.append(float(x))
##                 lflux.append(float(y))
##         lamb = np.asarray(lamb, dtype=FTYPE)
##         flux = np.asarray([10**(-0.4 * np.asarray(lflux)) * 4.65e-9])
##     else:
##         raise ValueError('Invalid ref: must be "kurucz" (default) or "hayes".')
##     lamb_e = lamb_to_edges(lamb)
##     maggies = project_filters(lamb_e, flux, filter_list)
##     maggies = maggies[0]
##     mag = -2.5 * np.log10(maggies)
##     return mag


## def ab2maggie(mag, dmag=None):
##     """Convert AB magnitude into AB maggies.

##     Parameters
##     ----------
##     mag : array_like
##         Array of AB magnitudes.
##     dmag : array_like, optional
##         Array of uncertainties in AB magnitudes.

##     Returns
##     -------
##     The function returns photometry in AB maggies.  If `dmag` is also
##     given, it is interpreted as magnitude uncertainty and the
##     corresponding uncertainty in AB maggies is returned.
##     """
##     mag = np.asarray(mag, dtype=FTYPE)
##     maggie = 10**(-0.4 * mag)
##     if dmag is None:
##         return maggie

##     dmag = np.asarray(dmag, dtype=FTYPE)
##     if dmag.shape == mag.shape:
##         dmaggie = 0.4 * np.log(10.) * np.fabs(maggie * dmag)
##         return maggie, dmaggie
##     else:
##         raise ValueError('Invalid dmag array shape.')


def maggie2ab(maggie, dmaggie=None):
    """Convert AB maggie into AB magnitude.

    If the uncertainties are given as the second argument, the
    uncertainties in magnitude are also returned.

    Parameters
    ----------
    maggie : array_like
        Array of AB maggies.
    dmaggie : array_like, optional
        Array of uncertainties in AB maggies.
    """
    maggie = np.asarray(maggie).astype(FTYPE)
    mag = -2.5 * np.log10(maggie)

    if dmaggie is None:
        return mag

    dmaggie = np.asarray(dmaggie, dtype=FTYPE)

    if dmaggie.shape == maggie.shape:
        dmag = 2.5 * np.fabs(dmaggie / maggie) / np.log(10.)
        return mag, dmag
    else:
        raise ValueError('Invalid dmaggie array shape.')


## def fit_nonneg(maggies, maggies_ivar, redshift, ptable,
##                maxiter=50000, tolerance=1e-6):
##     """
##     Notes
##     -----
##     See IDL procedure fit/k_fit_nonneg.pro in kcorrect v4_2.
##     """
##     return k_fit_nonneg(maggies, maggies_ivar, redshift,
##                         ptable[0], ptable[1], maxiter, tolerance, 0)


## def fit_photoz(maggies, maggies_ivar, ptable, zpriors=None, lpriors=None, 
##                maxiter=50000, tolerance=1e-6):
##     lpriors = np.asarray(np.zeros(2) if lpriors is None else lpriors,
##                          dtype=FTYPE)
##     zpriors = np.asarray(ZRANGE_DEFAULT[:2] if zpriors is None else zpriors,
##                          dtype=FTYPE) 
##     return k_fit_photoz(maggies,
##                         maggies_ivar,
##                         zpriors,
##                         lpriors,
##                         ptable[0], ptable[1],
##                         maxiter, tolerance, 0)


def reconstruct_maggie(coeffs, redshift, ptable):
    """Compute maggies given the fit coefficients and redshift.

    Parameters
    ----------
    coeffs -- 
    redshift --
    ptable --

    Notes
    -----
    The definition of an AB maggie is

    m = -2.5 log10( maggie )     (1)

    where m is in AB magnitudes.  From the definition of AB magnitude,
    this leads to

              \int dlamb lamb f_lamb(lamb) R(lamb)
    maggie = --------------------------------------
             \int dlamb lamb g_lamb(R,lamb) R(lamb)

    where lamb is the wavelength, R the filter transmission curve, and
    g is the reference spectrum, i.e., for the AB scale, the reference
    spectrum is flat in fnu space.

    What this function returns, however, is

                      \int dlamb (1+z) lamb L_lamb[lamb/(1+z)] R(lamb)
    reconst. maggie = --------------------------------------------------
                           \int dlamb lamb g_lamb(R,lamb) R(lamb)

    (see Eq.(13) of Hogg et al. 2002, astro-ph/0210394v1) where z is
    the redshift, and L_lamb is the rest-frame spectrum whose absolute
    scale is determined to match the observed flux at the redshift for
    which the given coeffs were obtained via non-negative fitting.
    This definition makes the expression for the k correction (not
    band-shifted) sensible, i.e.,

                               /  reconst. maggie(z = z) \
    k correction = -2.5 log10 | ------------------------- | .
                               \  reconst. maggie(z = 0) /

    This means, however, the reconstructed maggies cannot be converted
    into *observed* AB magnitude simply via (1), as the normalization
    of the term (1+z) lamb L_lamb[lamb/(1+z)] is only correct for the
    redshift at which coeffs are obtained via the nonnegative fitting.

    See IDL procedure fit/k_reconstruct_maggies.pro in kcorrect v4_2.
    """
    return k_reconstruct_maggies(coeffs, redshift, ptable[0], ptable[1])


## def abmag2flam(mag, lam, magerr=None):
##     """
##     Convert AB magnitude into flux in erg / s / cm^2 / angstrom.

##     Parameters
##     ----------
##     mag : array_like
##         Photometry in AB magnitude.
##     lam : array_like
##         Wavelength in angstrom at which the flux is measured.
##     magerr : array_like, optional
##         AB magnitude uncertainty.
##     """
##     flam = 10**( (mag + 48.60) / -2.5 ) * (+3e18 / lam**2)
##     if magerr is None:
##         return flam
##     dflam = flam * (1. - 10**(-0.4 * magerr))
##     return flam, dflam






class Photo(object):
    """Base class for basic photometry.

    Internally, kcorrect photometry uses AB maggie, which is not a
    very familiar system for astronomers.  One wants to subclass Photo
    to implement photometric conversion from the conventional
    photometry system to AB maggie by overriding the _set_input method.

    This class also implements other essential photometry routines.
    """

    def __init__(self, *args, **kwargs):
        """
        Parameters
        ----------
        cosmo : tuple of float, optional
            The cosmological density parameters and the Hubble
            constant normalized to 100 km/s/Mpc, as in (Omega_matter,
            Omega_Lambda, h100).  The default is COSMO_DEFAULT.
        ptable : ProjectionTableDB object, optional
            Template projection table.  If a non-standard projection
            table (i.e., with a custom redshift intervals and/or
            template spectra, etc.) is to be used in reconstructing AB
            maggies, then a custom ProjectionTableDB object can be
            specified here.  The default is PTABLE_MASTER.
        """
        self._cosmo = kwargs.get('cosmo', COSMO_DEFAULT)
        self._ptable = kwargs.get('ptable', PTABLE_MASTER)
        self._redshift = None
        self._set_input(*args, **kwargs)
        self._find_coeffs()

    def __len__(self):
        return self.nobj

    @property
    def nobj(self):
        """Number of objects."""
        return self._input_maggie.shape[0]

    @property
    def nfilter(self):
        """Number of filters."""
        return self._input_maggie.shape[1]

    def _set_input(self, maggie, maggie_ivar, filter_list, *args, **kwargs):
        """Set input data.

        The method ensures the input data integrity as well as the
        unit of input photometry and inverse variance, which must be
        in AB maggies.  The method should be overridden in subclass
        and implements a photometric conversion if the input are not
        in AB maggies.

        Parameters
        ----------
        filter_list : FilterList object
        maggie, maggie_ivar : array_like
            List of photometry in maggies and inverse variances.
        """
        # check input data integrity
        if not isinstance(filter_list, FilterList):
            raise TypeError('Expecting a FilterList object for filter_list.')
        maggie = np.atleast_2d(np.asarray(maggie, dtype=FTYPE))
        maggie_ivar = np.atleast_2d(np.asarray(maggie_ivar, dtype=FTYPE))
        if maggie.shape != maggie_ivar.shape:
            s = ('maggie and maggie_ivar must have the same shape.')
            raise ValueError(s)
        if maggie.shape[1] != len(filter_list):
            s = ('maggie must have a shape (n, m), for n objects observed in'
                 ' m bandpasses.')
            raise ValueError(s)

        # these attributes are always expected to exist
        self._input_maggie = maggie
        self._input_maggie_ivar = maggie_ivar
        self._filter_list = filter_list

    def _find_coeffs(self):
        """Do fitting to find best-fit coefficients."""
        # these attributes are always expected to exist
        self._coeffs = None
        self._chi2 = None
        self._iter = None

    def appmag(self, filter_list=None, redshift=None):
        """Compute apparent AB magnitude from best fit SED.

        The apparent AB magnitudes of best-fit SEDs are observed
        through the given filters.  The best-fit SEDs can obtionally
        be redshifted to arbitrary redshift(s).

        Parameters
        ----------
        filter_list : FilterList object, optional
            A set of filters through which best-fit SEDs are observed.
            The default is the input filters.
        redshift : int, float, or array_like, optional
            When a redshift or an array of redshifts (with the number
            of elements equal to the number of objects) is given, the
            best-fit SEDs are redshifted to the given redshifts before
            being observed.
        """
        filter_list = self._filter_list if filter_list is None else filter_list
        
        if redshift is None:
            # simply observe the input objects with given filters
            maggie = reconstruct_maggie(self._coeffs, self._redshift,
                                        self._ptable[filter_list])
            return maggie2ab(maggie)

        # when redshifts are given explicitly, observe objects at
        # those redshifts.
        if isinstance(redshift, (int, float)):
            z = redshift
            redshift = np.empty(self.nobj, dtype=FTYPE)
            redshift[:] = z
        else:
            redshift = np.reshape(redshift, (self.nobj,))

        maggie_S = reconstruct_maggie(self._coeffs, redshift,
                                      self._ptable[filter_list])
        m_S = maggie2ab(maggie_S)
        dm1 = ztodm(self._redshift, self._cosmo)
        dm2 = ztodm(redshift, self._cosmo)

        # because of the return value of reconstruct_maggie, the flux
        # needs to be rescaled by the ratio of distances
        return m_S - dm1 + dm2

    def model_spectrum(self, index=None, lmin=None, lmax=None):
        """Synthesize the best-fit SED(s) in the rest frame.

        Parameters
        ----------
        index : array of indices
            The index(ices) of the object for which spectrum is obtained.
        lmin : float
            Minimum wavelength in angstroms.
        lmax : float
            Maximum wavelength in angstroms.

        Returns
        -------
        The wavelengths and spectrum are returned.  The output
        spectrum is in units of erg s^-1 cm^-2 A^-1.
        """
        lamb, vmatrix = self._ptable[self._filter_list].templates
        lamb = 0.5 * (lamb[:-1] + lamb[1:])
        vmatrix = np.transpose(vmatrix)

        if (lmin is not None) or (lmax is not None):
            mask = np.ones(lamb.shape, dtype=bool)
            mask = mask if lmin is None else mask * np.greater(lamb, lmin)
            mask = mask if lmax is None else mask * np.less(lamb, lmax)
            lamb = lamb[mask]
            vmatrix = np.compress(mask, vmatrix, axis=0)

        index = (np.arange(0, len(self)) if index is None
                 else np.atleast_1d(np.asarray(index, dtype=int)))

        data = []
        for coeffs in self._coeffs[index]:
            flux = np.add.reduce(np.transpose(vmatrix * coeffs))
            data.append(flux)

        return (lamb, np.asarray(data) if len(data) > 1 else data[0])
