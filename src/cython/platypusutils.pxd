
cimport variant
cimport fastafile

from variant cimport Variant
from fastafile cimport FastaFile

cdef double logFactorial(int x)
cdef double binomial(int x, int size, double prob)
cdef double betaBinomialCDF(int k, int n, int alpha, int beta)
#cdef int permutations(int nObjects, int nChosen)
cdef list loadBAMData(list bamFiles, bytes chrom, int start, int end, options, list samples, dict samplesByID, dict samplesByBAM, char* refSeq)
cdef list getAlignmentErrorsBetweenReadsAndBestHaplotypes(list readBuffers, list haplotypes)
cpdef tuple pruned_read_start_end(read, int minq, int minAnchor)
cpdef tuple pruned_ref_start_end(read, int minq, int minAnchor)
cdef int isHaplotypeValid(tuple variants)
cdef Variant leftNormaliseIndel(Variant variant, FastaFile refFile, int maxReadLength)
