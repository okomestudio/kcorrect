cimport _ext



cdef float ohztor(float z, float o, float p):
    return _ext.ztor(z, o, p)


def ztor(z, o, p):
    return ohztor(z, o, p)
