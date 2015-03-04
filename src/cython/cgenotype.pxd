import cython
import chaplotype
import htslibWrapper

from chaplotype cimport Haplotype
from htslibWrapper cimport cAlignedRead

@cython.final
cdef class DiploidGenotype:
    cdef Haplotype hap1
    cdef Haplotype hap2
    cdef double hap1Like
    cdef double hap2Like
    cdef double calculateDataLikelihood(DiploidGenotype self, cAlignedRead** start, cAlignedRead** end, cAlignedRead** badReadsStart, cAlignedRead** badReadsEnd, cAlignedRead** brokenMatesStart, cAlignedRead** brokenMatesEnd, int individualIndex, int nIndividuals, double* gof, int useMapQualCap)

@cython.final
cdef class HaploidGenotype:
    cdef Haplotype hap1
    cdef double calculateDataLikelihood(self, cAlignedRead** start, cAlignedRead** end, int individualIndex, int nIndividuals)

cdef list generateAllGenotypesFromHaplotypeList(list haplotypes)
