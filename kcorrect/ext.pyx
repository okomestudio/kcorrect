cimport cext
import numpy as np


def k_projection_table(filter_curves,
                       lamb, vmatrix,
                       zvals, band_shift):
    """k_projection_table(fcrv,lambda,vmatrix,zvals,band_shift)

    DESCRIPTION

    Create the rmatrix, a lookup table which speeds analysis

    INPUT

    TO BE WRITTEN
    """
    nk = len(filter_curves)
    nl = filter_curves[0].size()
    nv = lamb.shape[0]
    nz = zvals.size()

    rmat = np.zeros((nk, nv, nz))

    ret = cext.k_projection_table(bandshift,)
    return

    ## IDL_LONG k_projection_table(float *rmatrix, IDL_LONG nk, IDL_LONG nv,
    ##                             float *vmatrix, float *lambda, IDL_LONG nl,
    ##                             float *zvals, IDL_LONG nz, IDL_LONG *filter_n,
    ##                             float *filter_lambda, float *filter_pass,
    ##                             float band_shift, IDL_LONG maxn);


## def k_reconstruct_maggies(coeffs, redshift, ptable[0], ptable[1]):
##     """k_reconstruct_maggies(coeffs,redshift,zvals,rmatrix)

##     DESCRIPTION

##       Interpolate the rmatrix

##     INPUT
    
##       TO BE WRITTEN
##     """
##     d = []
##     nd = 0
##     ret = cext.k_reconstruct_maggies(ptable[0],
##                                      ptable[0].shape[0],
##                                      ptable[1],
##                                      ptable[1].shape[0],
##                                      ptable[1].shape[1],
##                                      coeffs,
##                                      redshift,
##                                      d,
##                                      nd)
##     return d

    ## IDL_LONG k_reconstruct_maggies(float *zvals, IDL_LONG nz, float *rmatrix,    
    ##                                IDL_LONG nk, IDL_LONG nv, float *coeffs,
    ##                                float *galaxy_z, float *reconstruct_maggies,
    ##                                IDL_LONG ngalaxy)

    ## IDL_LONG k_read_ascii_table(float **table, IDL_LONG *ndim, IDL_LONG **sizes,



def ztor(z, o, p):
    return cext.ztor(z, o, p)


