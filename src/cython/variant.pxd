
###################################################################################################

import cython
from fastafile cimport FastaFile
from samtoolsWrapper cimport cAlignedRead

###################################################################################################

cdef public int PLATYPUS_VAR
cdef public int FILE_VAR
cdef public int ASSEMBLER_VAR

###################################################################################################

@cython.final
@cython.freelist(1000)
cdef class Variant:
    cdef:
        public bytes refName
        public bytes added
        public bytes removed
        public bytes bamAdded   # Insertion as reported in BAM file
        public bytes bamRemoved # Deletion as reported in BAM file
        public int refPos
        public int bamMinPos # Min Ref position as reported in BAM file
        public int bamMaxPos # Max Ref position as reported in BAM file
        public int minRefPos
        public int maxRefPos
        public int nSupportingReads
        public int varSource
        public int hashValue
        public int nAdded
        public int nRemoved
        public int varType
        double indelPrior(self, FastaFile refFile, int type)
        double calculatePrior(self, FastaFile refFile)
        void addVariant(self, Variant other)
        int overlaps(self, Variant other)

###################################################################################################
###################################################################################################

@cython.final
cdef class VariantCandidateGenerator:
    cdef int CIGAR_M
    cdef int CIGAR_I
    cdef int CIGAR_D
    cdef int CIGAR_N
    cdef int CIGAR_S
    cdef int CIGAR_H
    cdef int CIGAR_P
    cdef int minMapQual
    cdef int minBaseQual
    cdef int minFlank
    cdef int refId
    cdef int maxCoverage
    cdef int verbosity
    cdef int genSNPs
    cdef int genIndels
    cdef int maxReadLength
    cdef long int rStart
    cdef long int rEnd
    cdef long int refSeqStart
    cdef long int refSeqEnd
    cdef char* refSeq
    cdef dict variantHeap
    cdef bytes pyRefSeq
    cdef bytes rname
    cdef FastaFile refFile
    cdef object options
    cdef int qualBinSize

    cdef void addVariantToList(self, Variant var)
    cdef void getSnpCandidatesFromReadSegment(self, cAlignedRead* read, char* readSeq, char* readQual, int readStart, int readOffset, int refOffset, int lenSeqToCheck, int minFlank)
    cdef void getVariantCandidatesFromSingleRead(self, cAlignedRead* read)
    cdef void addCandidatesFromReads(self, cAlignedRead** readStart, cAlignedRead** readEnd)
    cdef list getCandidates(self, int minReads)

###################################################################################################
