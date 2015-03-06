cimport cython
cimport htslibWrapper
from htslibWrapper cimport cAlignedRead

cdef int hash_nucs
cdef int hash_size
cdef int max_sequence_length

cdef short* hash_sequence(char* sequence, int sequenceLength) nogil
cdef void hash_sequence_multihit(char* sequence, int sequenceLength, short** hash_table, short** next_array) nogil
cdef int mapAndAlignReadToHaplotype(char* read, char* quals, int readStart, int hapStart, int readLen, int hapLen, short* haplotypeHash, short* haplotypeNextArray, short* readHash, char* haplotype, int gapExtend, int nucprior, char* localGapOpen, int* mapCounts, int mapCountsLen, int hapFlank, int doCalculateFlankScore) nogil
cdef void hashReadForMapping(cAlignedRead* read) nogil
