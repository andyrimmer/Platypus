import cython

cimport htslibWrapper
cimport fastafile

from htslibWrapper cimport Samfile
from htslibWrapper cimport ReadIterator
from htslibWrapper cimport cAlignedRead
from fastafile cimport FastaFile

###################################################################################################

cdef class ReadArray:
    cdef cAlignedRead** array
    cdef cAlignedRead** windowStart
    cdef cAlignedRead** windowEnd
    cdef int __size
    cdef int __capacity
    cdef int __longestRead
    cdef int getSize(self)
    cdef void append(self, cAlignedRead* value)
    cdef void setWindowPointers(self, int start, int end)
    cdef void setWindowPointersBasedOnMatePos(self, int start, int end)
    cdef int getLengthOfLongestRead(self)
    cdef int countReadsCoveringRegion(self, int start, int end)

###################################################################################################

cdef class bamReadBuffer:
    cdef char* chrom
    cdef int chromID
    cdef int* filteredReadCountsByType
    cdef int isSorted
    cdef int startBase
    cdef int endBase
    cdef int windowStartBase
    cdef int windowEndBase
    cdef int maxReads
    cdef int minMapQual
    cdef int minBaseQual
    cdef int minFlank
    cdef int trimReadFlank
    cdef int verbosity
    cdef int minGoodBases
    cdef int trimOverlapping
    cdef int trimAdapter
    cdef int trimSoftClipped
    cdef cAlignedRead* lastRead
    cdef bytes sample
    cdef ReadArray reads
    cdef ReadArray badReads
    cdef ReadArray brokenMates
    cdef void setWindowPointers(self, int start, int end, int refStart, int refEnd, char* refSeq, int qualBinSize)
    cdef void recompressReadsInCurrentWindow(self, int refStart, int refEnd, char* refSeq, int qualBinSize, int compressReads)
    cdef void addReadToBuffer(self, cAlignedRead* theRead)
    cdef int countImproperPairs(self)
    cdef int countAlignmentGaps(self)
    cdef void sortReads(self)
    cdef void sortBrokenMates(self)
    cdef void logFilterSummary(self)
    cdef int countReadsCoveringRegion(self, int start, int end)

###################################################################################################
