import samtoolsWrapper
cimport samtoolsWrapper
from samtoolsWrapper cimport AlignedRead

###################################################################################################

cdef class BamFileIterator(object):
    cdef public object iterator
    cdef public AlignedRead currentValue
    cpdef next(self)

###################################################################################################

cdef class MultiBamFileIterator(object):
    cdef public list files
    cdef list queue
    cdef int nReadsInQueue
    cdef AlignedRead lastRead

###################################################################################################

cdef class MultiBamFileReader(object):
    cdef public object files
    cpdef getRName(self, int refId)
    cpdef MultiBamFileIterator fetch(self, char* region, int start, int end)

###################################################################################################
