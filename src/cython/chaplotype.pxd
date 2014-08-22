import cython

cimport fastafile
cimport samtoolsWrapper
cimport variant

from samtoolsWrapper cimport cAlignedRead
from fastafile cimport FastaFile
from variant cimport Variant

cdef int computeOverlapOfReadAndHaplotype(int hapStart, int hapEnd, cAlignedRead* theRead)

@cython.final
@cython.freelist(500)
cdef class Haplotype:
    cdef:
        bytes refName
        int startPos
        int endPos
        int maxReadLength
        int minVarPos
        int maxVarPos
        tuple variants
        bytes haplotypeSequence
        char* cHaplotypeSequence
        FastaFile refFile
        object options
        double* likelihoodCache
        int lenCache
        int lastIndividualIndex
        bytes referenceSequence
        int hash
        int hapLen
        int verbosity
        int endBufferSize
        char* cHomopolQ
        short* localGapOpen
        short* hapSequenceHash # For local re-mapping of reads to haplotype
        short* hapSequenceNextArray # same
        int* mapCounts # Store counts when mapping read to this haplotype
        int mapCountsLen
        double* alignReads(self, int individualIndex, cAlignedRead** start, cAlignedRead** end, cAlignedRead** badReadStart, cAlignedRead** badReadsEnd, cAlignedRead** brokenReadsStart, cAlignedRead** brokenReadsEnd, int useMapQualCap)
        double alignSingleRead(self, cAlignedRead* theRead, int useMapQualCap)
        char* getReferenceSequence(self, prefix=*)
        char* getMutatedSequence(self)
        dict vcfINFO(self)
        bytes getSequenceContext(self, Variant variant)
        int homopolymerLengthForOneVariant(self, Variant variant)
        list homopolymerLengths(self)
        void annotateWithGapOpen(self)
