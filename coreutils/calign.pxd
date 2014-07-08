cimport cython
cimport samtoolsWrapper
from samtoolsWrapper cimport cAlignedRead

cdef int hash_nucs
cdef int hash_size
cdef int max_sequence_length

cdef short* hash_sequence(char* sequence, int sequenceLength)
cdef void hash_sequence_multihit(char* sequence, int sequenceLength, short** hash_table, short** next_array)
cdef int mapAndAlignReadToHaplotype(char* read, char* quals, int readStart, int hapStart, int readLen, int hapLen, short* haplotypeHash, short* haplotypeNext, short* readHash, char* haplotype, int gapExtend, int nucprior, char* errorModel, int useHomopolQ, char* aln1, char* aln2, int* bestMappingPosition)
cdef void hashReadForMapping(cAlignedRead* read)
