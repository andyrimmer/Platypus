import variant
import htslibWrapper
import cwindow

from variant cimport Variant
from htslibWrapper cimport cAlignedRead
from cwindow cimport bamReadBuffer
from fastafile cimport FastaFile

###################################################################################################

cdef class Population:

    # Class variables
    cdef:
        dict vcfInfo
        dict vcfFilter
        dict varsByPos
        dict variantPosteriors
        list genotypeCalls
        list variants
        list haplotypes
        list genotypes
        double** genotypeLikelihoods
        double** goodnessOfFitValues
        double** EMLikelihoods
        double* frequencies
        double* newFrequencies
        int** haplotypeIndexes
        int* nReads
        double* maxLogLikelihoods
        double* freqsPrimeByHapIndex
        int nVariants
        int nHaplotypes
        int nIndividuals
        int maxHaplotypes
        int maxGenotypes
        int nGenotypes
        int verbosity
        int useEMLikelihoods
        list readBuffers
        FastaFile refFile
        object options

    # Class functions
    cdef:
        double calculatePosterior(self, Variant var, int flatPrior=*)
        void computeVariantPosteriors(self)
        void computeVariantINFO(self)
        void reset(self)
        void computeVariantFILTER(self)
        void callGenotypes(self)
        void call(self, int maxIters, int computeVCFFields)
        double EMiteration(self, double* freq, double* newFreqs)
        void setup(self, list variants, list haplotypes, list genotypes, int nInd, int verbosity, list readBuffers)

###################################################################################################
