cdef extern from "include/kcorrect.h":

    ctypedef int IDL_LONG

    IDL_LONG k_projection_table(float *rmatrix, IDL_LONG nk, IDL_LONG nv,
                                float *vmatrix, float *lambda, IDL_LONG nl,
                                float *zvals, IDL_LONG nz, IDL_LONG *filter_n,
                                float *filter_lambda, float *filter_pass,
                                float band_shift, IDL_LONG maxn);

    ## IDL_LONG k_reconstruct_maggies(float *zvals, IDL_LONG nz, float *rmatrix,    
    ##                                IDL_LONG nk, IDL_LONG nv, float *coeffs,
    ##                                float *galaxy_z, float *reconstruct_maggies,
    ##                                IDL_LONG ngalaxy)

    IDL_LONG k_read_ascii_table(float **table, IDL_LONG *ndim, IDL_LONG **sizes,

    float ztor(float z, float omega0, float omegal0)
    float ztodV(float z, float omega0, float omegal0)
