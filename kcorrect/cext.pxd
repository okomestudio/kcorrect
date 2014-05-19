"""Expose C library."""


cdef extern from "include/kcorrect.h":

    ctypedef int IDL_LONG

    IDL_LONG k_projection_table(float *rmatrix, IDL_LONG nk, IDL_LONG nv,
                                float *vmatrix, float *lamb, IDL_LONG nl,
                                float *zvals, IDL_LONG nz, IDL_LONG *filter_n,
                                float *filter_lambda, float *filter_pass,
                                float band_shift, IDL_LONG maxn)

    IDL_LONG k_reconstruct_maggies(float *zvals, IDL_LONG nz, float *rmatrix,
                                   IDL_LONG nk, IDL_LONG nv, float *coeffs,
                                   float *galaxy_z, float *reconstruct_maggies,
                                   IDL_LONG ngalaxy)

    IDL_LONG k_fit_nonneg(float *coeffs, float *rmatrix, IDL_LONG nk,
                          IDL_LONG nv, float *zvals, IDL_LONG nz,
                          float *maggies, float *maggies_ivar,
                          float *redshift, IDL_LONG ngalaxy,
                          float tolerance, IDL_LONG maxiter, IDL_LONG *niter, 
                          float *chi2, IDL_LONG verbose, IDL_LONG dontinit)

    IDL_LONG k_fit_photoz(float *photoz, float *coeffs, float *rmatrix, 
                          IDL_LONG nk, IDL_LONG nv,
                          float *lprior, float *zprior,
                          IDL_LONG nprior, float *zvals, 
                          IDL_LONG nz, float *maggies, float *maggies_ivar, 
                          IDL_LONG ngalaxy, float tolerance, IDL_LONG maxiter, 
                          IDL_LONG *niter, float *chi2, IDL_LONG verbose)
    
    IDL_LONG k_read_ascii_table(float **table,
                                IDL_LONG *ndim,
                                IDL_LONG **sizes,
                                char *fname)

    float ztor(float z, float omega0, float omegal0)
    float ztodV(float z, float omega0, float omegal0)
    float ztoV(float z, float omega0, float omegal0)
    float Vtoz(float V, float omega0, float omegal0)
    float rtoz(float r, float omega0, float omegal0)
    float z2dm(float z, float omega0, float omegal0)
    float dm2z(float dm, float omega0, float omegal0)
    float z2add(float z, float omega0, float omegal0)
    float z2t(float z, float omega0, float omegal0)
    float t2z(float z, float omega0, float omeagl0)
