
import chaplotype
import samtoolsWrapper

from chaplotype cimport Haplotype
from samtoolsWrapper cimport cAlignedRead

cdef class DiploidGenotype:
    cdef Haplotype hap1
    cdef Haplotype hap2
    cdef double hap1Like
    cdef double hap2Like
    cdef double calculateDataLikelihood(DiploidGenotype self, cAlignedRead** start, cAlignedRead** end, cAlignedRead** badReadsStart, cAlignedRead** badReadsEnd, cAlignedRead** brokenMatesStart, cAlignedRead** brokenMatesEnd, int individualIndex, int nIndividuals, double* gof, int printAlignments)

cdef class HaploidGenotype:
    cdef Haplotype hap1
    cdef double calculateDataLikelihood(self, cAlignedRead** start, cAlignedRead** end, int individualIndex, int nIndividuals)

cdef list generateAllGenotypesFromHaplotypeList(list haplotypes)
