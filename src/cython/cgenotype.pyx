"""
Fast implementation of the Genotype class.
"""

from __future__ import division
import logging
import cython
cimport cython

cimport chaplotype
cimport samtoolsWrapper
cimport variant

from chaplotype cimport Haplotype
from variant cimport Variant

###################################################################################################

logger = logging.getLogger("Log")

###################################################################################################

cdef double mLTOT = -0.23025850929940459    # Minus log ten over ten
cdef double log10E = 0.43429448190325182
cdef double log10Half = -0.3010299956639812
cdef double log10Two = 0.3010299956639812
cdef double logTwo = 0.69314718055994529
cdef double logHalf = -0.69314718055994529

###################################################################################################

cdef extern from "math.h":
    double exp(double)
    double log(double)
    double log10(double)

###################################################################################################

@cython.final
@cython.freelist(2000)
cdef class HaploidGenotype(object):
    """
    This class represents a haploid genotype. It stores one
    haplotypes.
    """
    def __init__(self, Haplotype hap1):
        """
        Constructor. Takes haplotypes.
        """
        self.hap1 = hap1

    def __str__(self):
        """
        Generate a string representation of the genotype class
        """
        theString = "{ ["

        for v in self.hap1.variants:
            theString += v.shortRepr()

        theString += "] }"
        return theString.ljust(30)

    def __repr__(self):
        """
        Generate a string representation of the genotype class
        """
        return self.__str__()

    cdef double calculateDataLikelihood(self, cAlignedRead** start, cAlignedRead** end, int individualIndex, int nIndividuals):
        """
        """
        cdef double likelihood = 0.0
        cdef Py_ssize_t readIndex = 0
        cdef double like1 = 0.0
        #cdef double* arr1 = self.hap1.alignReads(individualIndex, start, end, nReads, nIndividuals)

        #for readIndex from 0 <= readIndex < nReads:
        #    like1 = arr1[readIndex]
        #    likelihood += log(like1)

        return likelihood

###################################################################################################

@cython.final
@cython.freelist(2000)
cdef class DiploidGenotype(object):
    """
    This class represents a diploid genotype. It stores two
    haplotypes.
    """
    def __init__(self, Haplotype hap1, Haplotype hap2):
        """
        Constructor. Takes haplotypes.
        """
        self.hap1 = hap1
        self.hap2 = hap2

    def __contains__(self, Variant v):
        """
        Check if a variant exists in this genotype.
        """
        if v in self.hap1.variants or v in self.hap2.variants:
            return True
        else:
            return False

    def __str__(self):
        """
        Generate a string representation of the genotype class
        """
        theString = "{ ["

        for v in self.hap1.variants:
            theString += v.shortRepr()

        theString += "] , ["

        for v in self.hap2.variants:
            theString += v.shortRepr()

        theString += "] }"

        return theString.ljust(75)

    def __repr__(self):
        """
        Generate a string representation of the genotype class
        """
        return self.__str__()

    cdef double calculateDataLikelihood(DiploidGenotype self, cAlignedRead** start, cAlignedRead** end, cAlignedRead** badReadsStart, cAlignedRead** badReadsEnd, cAlignedRead** brokenReadsStart, cAlignedRead** brokenReadsEnd, int individualIndex, int nIndividuals, double* gof):
        """
        """
        cdef double likelihood = 0.0
        cdef Py_ssize_t readIndex = 0
        cdef double like1 = 0.0
        cdef double like2 = 0.0
        cdef double* arr1 = self.hap1.alignReads(individualIndex, start, end, badReadsStart, badReadsEnd, brokenReadsStart, brokenReadsEnd, False)
        cdef double* arr2 = self.hap2.alignReads(individualIndex, start, end, badReadsStart, badReadsEnd, brokenReadsStart, brokenReadsEnd, False)
        cdef double goodnessOfFitValue = 0.0
        cdef double logLike1 = 0.0
        cdef double logLike2 = 0.0
        cdef int nReads = end - start
        cdef int nBadReads = badReadsEnd - badReadsStart
        cdef int nBrokenPairReads = brokenReadsEnd - brokenReadsStart
        cdef int totalReads = nReads + nBadReads + nBrokenPairReads

        self.hap1Like = 0
        self.hap2Like = 0

        for readIndex from 0 <= readIndex <= totalReads:
            like1 = arr1[readIndex] # These are log values
            like2 = arr2[readIndex] # These are log values

            if like1 == 999 and like2 == 999:
                break

            logLike1 = log10E*(like1) # Convert from log to log10
            logLike2 = log10E*(like2) # Convert from log to log10
            self.hap1Like += logLike1
            self.hap2Like += logLike2
            goodnessOfFitValue += max(logLike1, logLike2)

            # Special-case optimisation for homozygous genotype.
            if arr1 == arr2:
                likelihood += like1

            # Approximation good to 0.1%: take the highest log value
            elif abs(like1 - like2) >= 3:
                likelihood += (logHalf + max(like1, like2))

            # Haplotypes give the same likelihood. This happens when either they are both rubbish fits, and
            # hit the likelihood cap, or are both perfect fits, or the read lies in a position such that it cannot
            # be used to distinguish between the two haplotypes.
            elif abs(like1 - like2) <= 1e-3:
                likelihood += like1

            # Currently does full log and exp computation.
            else:
                likelihood += log(0.5*(exp(like1) + exp(like2)))

        if gof and nReads > 0:
            gof[individualIndex] = (-10*goodnessOfFitValue) / nReads
        elif gof:
            gof[individualIndex] = 0
        else:
            pass

        return likelihood

###################################################################################################

cdef list generateAllGenotypesFromHaplotypeList(list haplotypes):
    """
    Return a list of genotypes, corresponding to all allowed combinations of
    haplotypes in the specified region.

    e.g. if nHaplotypes = 2, and n = 2, then the list returned is
    [(0,0), (0,1), (1,1)]

    i.e. we are allowed 2 sets of haplotype 0, 2 sets of haplotype 1, or one of
    each. This would be appropriate for a single diploid individual if we were
    considering 2 haplotypes.
    """
    cdef list genotypes = []
    cdef int nHaplotypes = len(haplotypes)
    cdef int indexOne = 0
    cdef int indexTwo = 0
    cdef Haplotype haplotypeOne
    cdef Haplotype haplotypeTwo

    for indexOne from 0 <= indexOne < nHaplotypes:
        for indexTwo from indexOne <= indexTwo < nHaplotypes:
            haplotypeOne = haplotypes[indexOne]
            haplotypeTwo = haplotypes[indexTwo]
            genotypes.append(DiploidGenotype(haplotypeOne, haplotypeTwo))

    return genotypes

###################################################################################################
