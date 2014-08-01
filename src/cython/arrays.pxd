"""
Cython header for arrays module
"""
cimport cython

###################################################################################################

cdef class IntArray:
    cdef int* array
    cdef int __size
    cdef int __capacity
    cdef int getSize(self)
    cdef int getLast(self)
    cdef int at(self, int index)
    cdef void addValue(self, int index, int value)
    cdef void append(self, int value)

###################################################################################################

cdef class DoubleArray:
    cdef double* array
    cdef int __size
    cdef int __capacity
    cdef int getSize(self)
    cdef double getLast(self)
    cdef double at(self, int index)
    cdef void addValue(self, int index, double value)
    cdef void append(self, double value)

###################################################################################################
