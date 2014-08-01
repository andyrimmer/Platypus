cimport arrays
from arrays cimport IntArray
from arrays cimport DoubleArray

###################################################################################################

cdef class Histogrammer:
    cdef IntArray bins
    cdef DoubleArray counts
    cdef int granularity
    cdef int _histbin(self, int value)
    cdef void fast_enter(self, int value, int number)
    cdef void enter(self, int value, double number=*)

###################################################################################################
