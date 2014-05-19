cimport cext
import numpy as np
cimport numpy as np
from globaldefs import FTYPE, ITYPE
import ctypes


class FilterCurves(object):
    """Interface for filter information.

    Given a set of filter info, transform it into the form that the C
    extension can understand.

    """

    def __init__(self, input):
        # the number of filters
        nk = len(input)
        # number of data points in each curve
        filter_n = np.empty(nk, dtype=ITYPE)
        # max number of curves
        maxn = 0

        for i in xrange(nk):
            lamb = input[i][0]
            passbd = input[i][1]
            filter_n[i] = lamb.size()
            maxn = max(maxn, lamb.size())

        # array of wavelengths
        filter_lambda = np.empty(maxn * nk, dtype=FTYPE)
        # array of transmission curve
        filter_pass = np.empty(maxn * nk, dtype=FTYPE)

        for i in xrange(nk):
            filter_lambda[:maxn * i] = input[i][0]
            filter_pass[:maxn * i] = input[i][1]
            filter_lambda[maxn*i:] = 0.0
            filter_pass[maxn*i:] = 0.0

        self.nk = nk
        self.filter_lambda = filter_lambda
        self.filter_pass = filter_pass
        self.filter_n = filter_n
        self.maxn = maxn


def k_projection_table(filter_curves,
                       np.ndarray lamb,
                       np.ndarray vmatrix,
                       np.ndarray zvals,
                       float band_shift):
    """k_projection_table(fcrv,lambda,vmatrix,zvals,band_shift)

    DESCRIPTION

    Create the rmatrix, a lookup table which speeds analysis

    INPUT

    TO BE WRITTEN
    """
    filter_curves = FilterCurves(filter_curves)

    cdef int nk = filter_curves.nk
    cdef int nl = lamb.shape[0]
    cdef int nv = vmatrix.shape[0]
    cdef int nz = zvals.shape[0]

    cdef np.ndarray rmatrix = np.zeros((nk, nv, nz), dtype=FTYPE)
    cdef np.ndarray filter_n = filter_curves.filter_n
    cdef np.ndarray filter_lambda = filter_curves.filter_lambda
    cdef np.ndarray filter_pass = filter_curves.filter_pass

    cext.k_projection_table(<float*>rmatrix.data, nk, nv,
                            <float*>vmatrix.data,
                            <float*>lamb.data,
                            nl - 1,
                            <float*>zvals.data, nz,
                            <cext.IDL_LONG*>filter_n.data,
                            <float*>filter_lambda.data,
                            <float*>filter_pass.data,
                            band_shift,
                            filter_curves.maxn);
    return rmatrix


def k_reconstruct_maggies(np.ndarray coeffs,
                          np.ndarray redshift,
                          np.ndarray zvals,
                          np.ndarray rmatrix):
    """k_reconstruct_maggies(coeffs,redshift,zvals,rmatrix)

    DESCRIPTION

      Interpolate the rmatrix

    INPUT
    
      TO BE WRITTEN
    """
    cdef int nz = zvals.shape[0]
    cdef int nk = rmatrix.shape[0]
    cdef int nv = rmatrix.shape[1]
    cdef int ngalaxy = redshift.shape[0]

    cdef np.ndarray reconstruct_maggies = np.empty((ngalaxy, nk), dtype=FTYPE)

    cext.k_reconstruct_maggies(<float*>zvals.data, nz,
                               <float*>rmatrix.data, nk, nv,
                               <float*>coeffs.data,
                               <float*>redshift.data,
                               <float*>reconstruct_maggies.data,
                               ngalaxy)
    return reconstruct_maggies


def k_fit_nonneg(np.ndarray maggies,
                 np.ndarray maggies_ivar,
                 np.ndarray redshift,
                 np.ndarray zvals,
                 np.ndarray rmatrix,
                 int maxiter,
                 float tolerance,
                 int verbose):
    """k_fit_nonneg(maggies,maggies_ivar,redshift,zvals,rmatrix,maxiter,tolerange,verbose)

    DESCRIPTION

    Nonnegative solutions

    INPUT

    TO BE WRITTEN
    """
    cdef int nk = rmatrix.shape[0]
    cdef int nv = rmatrix.shape[1]
    cdef int nz = zvals.shape[0]
    cdef int ngalaxy = redshift.shape[0]
    cdef int niter = 0

    cdef np.ndarray coeffs = np.empty((ngalaxy, nv), dtype=FTYPE)
    cdef np.ndarray chi2 = np.empty(ngalaxy, dtype=FTYPE)

    cext.k_fit_nonneg(<float*>coeffs.data,
                      <float*>rmatrix.data, nk, nv,
                      <float*>zvals.data, nz,
                      <float*>maggies.data, <float*>maggies_ivar.data,
                      <float*>redshift.data, ngalaxy,
                      tolerance,
                      maxiter,
                      &niter,
                      <float*>chi2.data,
                      verbose,
                      0)

    return coeffs, chi2, niter


def k_fit_photoz(np.ndarray maggies,
                 np.ndarray maggies_ivar,
                 np.ndarray zpriors,
                 np.ndarray lpriors,
                 np.ndarray zvals,
                 np.ndarray rmatrix,
                 int maxiter,
                 float tolerance,
                 int verbose):
    """k_fit_photoz(maggies,maggies_ivar,redshift,zvals,rmatrix,maxiter,tolerange,verbose)
    
    DESCRIPTION
    
    Get photometric redshifts
    
    INPUT
    
    TO BE WRITTEN

    """
    cdef int nk = rmatrix.shape[0]
    cdef int nv = rmatrix.shape[1]
    cdef int nz = zvals.shape[0]
    cdef int ngalaxy = maggies.shape[0]
    cdef int npriors = lpriors.shape[0]
    cdef int niter = 0

    cdef np.ndarray photoz = np.empty(ngalaxy, dtype=FTYPE)
    cdef np.ndarray coeffs = np.empty((ngalaxy, nv), dtype=FTYPE)
    cdef np.ndarray chi2 = np.empty(ngalaxy, dtype=FTYPE)
    
    cext.k_fit_photoz(<float*>photoz.data,
                      <float*>coeffs.data,
                      <float*>rmatrix.data, nk, nv,
                      <float*>lpriors.data,
                      <float*>zpriors.data,
                      npriors,
                      <float*>zvals.data, nz,
                      <float*>maggies.data,
                      <float*>maggies_ivar.data, ngalaxy,
                      tolerance,
                      maxiter, 
                      &niter,
                      <float*>chi2.data,
                      verbose)

    return photoz, coeffs, chi2, niter


def k_read_ascii_table(char* fname):
    """k_read_ascii_table(fname)

    Read an ASCII file in the following format, which is:
 
    <ndim> <size_{0}> ... <size_{ndim-1}>
    <entry_0>
    <entry_1>
    ...
    <entry_n>
    ...
    <entry_{size_0*size_1*..*size_{ndim-1}-1>
 
    where the table element [k,j,i] (for ndim==3) would be the entry
    element n, where n=i*size_2*size1+j*size2+k.

    INPUT

    fname -- Path to an ASCII file

    """
    cdef float* table = NULL
    cdef int nd = 0
    cdef int* sizes = NULL

    cext.k_read_ascii_table(&table,
                            &nd,
                            &sizes,
                            fname)

    n = 1
    shape = []
    for i in xrange(nd):
        x = sizes[i]
        n *= x
        shape.append(x)

    cdef np.ndarray t = np.empty(n, dtype=FTYPE)
    for i in xrange(n):
        t[i] = table[i]

    return t.reshape(shape)


def ztor(float z, float o, float p):
    return cext.ztor(z, o, p)


def ztodV(float z, float omega0, float omegal0):
    return ztodV(z, omega0, omegal0)


def ztoV(float z, float omega0, float omegal0):
    return ztoV(z, omega0, omegal0)


def Vtoz(float V, float omega0, float omegal0):
    return Vtoz(V, omega0, omegal0)


def rtoz(float r, float omega0, float omegal0):
    return rtoz(r, omega0, omegal0)


def z2dm(float z, float omega0, float omegal0):
    return z2dm(z, omega0, omegal0)


def dm2z(float dm, float omega0, float omegal0):
    return dm2z(dm, omega0, omegal0)


def z2add(float z, float omega0, float omegal0):
    return z2add(z, omega0, omegal0)


def z2t(float z, float omega0, float omegal0):
    return z2t(z, omega0, omegal0)


def t2z(float z, float omega0, float omeagl0):
    return t2z(z, omega0, omeagl0)
