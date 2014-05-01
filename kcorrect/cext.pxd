cdef extern from "include/kcorrect.h":
    float ztor(float z, float omega0, float omegal0)
    float ztodV(float z, float omega0, float omegal0)
