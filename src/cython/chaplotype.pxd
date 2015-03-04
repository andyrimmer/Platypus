import cython

cimport fastafile
cimport htslibWrapper
cimport variant

from htslibWrapper cimport cAlignedRead
from fastafile cimport FastaFile
from variant cimport Variant

cdef int computeOverlapOfReadAndHaplotype(int hapStart, int hapEnd, cAlignedRead* theRead)

@cython.final
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
        bytes shortHaplotypeSequence
        char* cHaplotypeSequence
        FastaFile refFile
        object options
        double* likelihoodCache
        int lenCache
        int lastIndividualIndex
        bytes referenceSequence
        bytes shortReferenceSequence
        Variant longVar
        int hash
        int hapLen
        int verbosity
        int endBufferSize
        char* cHomopolQ
        char* localGapOpen
        short* hapSequenceHash # For local re-mapping of reads to haplotype
        short* hapSequenceNextArray # same
        int* mapCounts # Store counts when mapping read to this haplotype
        int mapCountsLen
        double* alignReads(self, int individualIndex, cAlignedRead** start, cAlignedRead** end, cAlignedRead** badReadStart, cAlignedRead** badReadsEnd, cAlignedRead** brokenReadsStart, cAlignedRead** brokenReadsEnd, int useMapQualCap)
        double alignSingleRead(self, cAlignedRead* theRead, int useMapQualCap)
        char* getReferenceSequence(self, prefix=*)
        char* getMutatedSequence(self)
        char* getShortHaplotypeSequence(self)
        dict vcfINFO(self)
        bytes getSequenceContext(self, Variant variant)
        int homopolymerLengthForOneVariant(self, Variant variant)
        list homopolymerLengths(self)
        void annotateWithGapOpen(self)
