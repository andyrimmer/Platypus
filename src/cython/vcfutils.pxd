cimport fastafile
cimport cpopulation

from cpopulation cimport Population
from fastafile cimport FastaFile
from chaplotype cimport Haplotype
from variant cimport Variant

cdef void outputCallToVCF(dict varsByPos, dict vcfInfo, dict vcfFilter, list haplotypes, list genotypes, double* haplotypeFrequencies, double** genotypeLikelihoods, double** gofValues, int** haplotypeIndexes, list readBuffers, int nIndividuals, vcfFile, FastaFile refFile, outputFile, options, list allVariants, int windowStart, int windowEnd)
cdef tuple refAndAlt(char* chrom, int POS, list variants, FastaFile refFile)
cdef dict vcfINFO(double* haplotypeFrequencies, dict variantPosteriors, list genotypeCalls, list genotypes, list haplotypes, list readBuffers, int nHaplotypes, options, FastaFile refFile)
cdef dict vcfFILTER(list genotypeCalls, list haplotypes, dict vcfInfo, dict varsByPos, options)
cdef void outputHLACallToVCF(list haplotypes, list readBuffers, int nIndividuals, FastaFile refFile, outputFile, options, int windowStart, int windowEnd)