"""
This module implements the core EM algorithm for inferring population
haplotype frequencies.
"""
from __future__ import division

import cython
cimport cython

import logging
import math
import cwindow
import chaplotype
import random

cimport variant
cimport chaplotype
cimport cgenotype
cimport vcfutils
cimport cerrormodel
cimport samtoolsWrapper

from variant cimport Variant
from chaplotype cimport Haplotype
from cgenotype cimport DiploidGenotype
from cwindow cimport bamReadBuffer
from samtoolsWrapper cimport cAlignedRead

logger = logging.getLogger("Log")

###################################################################################################

ctypedef long long size_t

###################################################################################################

cdef double PI = math.pi
cdef double mLTOT = -0.23025850929940459    # Minus log ten over ten

###################################################################################################

cdef extern from "stdlib.h":
    void free(void *)
    void* malloc(size_t)
    void* calloc(size_t,size_t)
    void* realloc(void *,size_t)
    void* memset(void* ptr, int value, size_t num)

###################################################################################################

cdef extern from "math.h":
    double exp(double)
    double log(double)
    double log10(double)
    double fabs(double)
    double round(double)

###################################################################################################

cdef void my_free(void* thePointer):
    """
    Cython wrapper. Used for profiling.
    """
    free(thePointer)

###################################################################################################

cdef void* my_malloc(size_t theSize):
    """
    Cython wrapper. Used for profiling.
    """
    return malloc(theSize)

###################################################################################################

cdef void* my_calloc(size_t theSize1, size_t theSize2):
    """
    Cython wrapper. Used for profiling.
    """
    return calloc(theSize1, theSize2)

###################################################################################################

cdef class Population:
    """
    The Population Caller takes a list of haplotypes and a list of reads and
    computes a MLE for the population frequency of each of the variants.
    """
    def __init__(self, options):
        """
        Constructor allocates storage.
        """
        self.options = options
        self.verbosity = options.verbosity
        self.nIndividuals = options.nInd
        self.maxHaplotypes = options.maxHaplotypes
        self.maxGenotypes = options.maxGenotypes
        self.genotypeCalls = []
        self.vcfInfo = {}
        self.vcfFilter = {}
        self.variantPosteriors = {}
        self.varsByPos = {}
        self.useEMLikelihoods = options.useEMLikelihoods
        self.refFile = None

        # Allocate storage once
        self.freqsPrimeByHapIndex = <double*>(my_calloc(self.options.maxHaplotypes, sizeof(double)))
        self.genotypeLikelihoods = <double**>(my_calloc(self.nIndividuals, sizeof(double*)))
        self.EMLikelihoods = <double**>(my_calloc(self.nIndividuals, sizeof(double*)))
        self.haplotypeIndexes = <int**>(my_calloc(self.options.maxGenotypes, sizeof(int*)))
        self.nReads = <int*>(my_calloc(self.nIndividuals, sizeof(int)))
        self.maxLogLikelihoods = <double*>(my_calloc(self.nIndividuals, sizeof(double)))
        self.goodnessOfFitValues = <double**>(my_calloc(self.options.maxGenotypes, sizeof(double*)))
        self.frequencies = <double*>(my_calloc(self.options.maxHaplotypes, sizeof(double)))
        self.newFrequencies = <double*>(my_calloc(self.options.maxHaplotypes, sizeof(double)))

        cdef int individualIndex = 0
        cdef int genotypeIndex = 0

        for individualIndex from 0 <= individualIndex < self.nIndividuals:
            self.genotypeLikelihoods[individualIndex] = <double*>(my_calloc(self.options.maxGenotypes, sizeof(double)))
            self.EMLikelihoods[individualIndex] = <double*>(my_calloc(self.options.maxGenotypes, sizeof(double)))

        for genotypeIndex from 0 <= genotypeIndex < self.options.maxGenotypes:
            self.haplotypeIndexes[genotypeIndex] = <int*>(my_calloc(2, sizeof(int)))
            self.goodnessOfFitValues[genotypeIndex] = <double*>(my_calloc(self.nIndividuals, sizeof(double)))

    def __dealloc__(self):
        """
        De-allocate memory. Make sure to do this correctly for 2-d
        arrays!
        """
        cdef int index = 0
        cdef int genIndex = 0
        cdef int hapIndex = 0

        for index from 0 <= index < self.nIndividuals:
            my_free(self.EMLikelihoods[index])
            my_free(self.genotypeLikelihoods[index])

        for genIndex from 0 <= genIndex < self.options.maxGenotypes:
            my_free(self.haplotypeIndexes[genIndex])
            my_free(self.goodnessOfFitValues[genIndex])

        my_free(self.haplotypeIndexes)
        my_free(self.genotypeLikelihoods)
        my_free(self.goodnessOfFitValues)
        my_free(self.frequencies)
        my_free(self.newFrequencies)
        my_free(self.nReads)
        my_free(self.maxLogLikelihoods)
        my_free(self.EMLikelihoods)
        my_free(self.freqsPrimeByHapIndex)

    cdef void computeVariantINFO(self):
        """
        """
        self.vcfInfo = vcfutils.vcfINFO(self.frequencies, self.variantPosteriors, self.genotypeCalls, self.genotypes, self.haplotypes, self.readBuffers, self.nHaplotypes, self.options, self.refFile)

    cdef void computeVariantFILTER(self):
        """
        FILTER field for vcf output, two level dictionary of variant - info field - value
        """
        self.vcfFilter = vcfutils.vcfFILTER(self.genotypeCalls, self.haplotypes, self.vcfInfo, self.varsByPos, self.options)

    cdef void reset(self):
        """
        Reset all relevant member variables of the population instance, so it is ready for re-use in the
        next window.
        """
        self.vcfInfo = {}
        self.vcfFilter = {}
        self.varsByPos = {}
        self.variantPosteriors = {}
        self.genotypeCalls = []
        self.variants = []
        self.haplotypes = []
        self.genotypes = []

        memset(self.freqsPrimeByHapIndex, 0, sizeof(double)*self.maxHaplotypes)
        memset(self.nReads, 0, sizeof(int)*self.nIndividuals)
        memset(self.maxLogLikelihoods, 0, sizeof(double)*self.nIndividuals)
        memset(self.frequencies, 0, sizeof(double)*self.maxHaplotypes)
        memset(self.newFrequencies, 0,  sizeof(double)*self.maxHaplotypes)

        cdef int individualIndex = 0
        cdef int genotypeIndex = 0

        for individualIndex from 0 <= individualIndex < self.nIndividuals:
            memset(self.genotypeLikelihoods[individualIndex], 0, sizeof(double)*self.maxGenotypes)
            memset(self.EMLikelihoods[individualIndex], 0, sizeof(double)*self.maxGenotypes)

        for genotypeIndex from 0 <= genotypeIndex < self.options.maxGenotypes:
            memset(self.haplotypeIndexes[genotypeIndex], 0, sizeof(int)*2)
            memset(self.goodnessOfFitValues[genotypeIndex], 0, sizeof(double)*self.nIndividuals)

    cdef void setup(self, list variants, list haplotypes, list genotypes, int nInd, int verbosity, list readBuffers):
        """
        Constructor for the population caller. The initial loop over reads, calculating
        likelihoods for each read, given a particular haplotype, is done here.

        We cache posterior values, the indexes in 'haplotypes' of haplotypes for each genotype,
        which are used in each step of the EM algorithm but do not change.
        """
        self.nGenotypes = len(genotypes)
        self.nVariants = len(variants)
        self.nHaplotypes = len(haplotypes)

        if not nInd == len(readBuffers):
            logger.error("nInd != len(readBuffers). Something is very wrong here. Quitting now.")
            logger.info("Variants = %s" %(variants))
            logger.info("nVariants = %s. nHaplotypes = %s. nGenotypes = %s" %(len(variants), self.nHaplotypes, self.nGenotypes))
            raise StandardError, "Error in cPopulation.setup"

        if not self.nGenotypes <= self.options.maxGenotypes:
            logger.error("nGenotypes > options.maxGenotypes. (%s > %s). Something is very wrong here. Quitting now." %(self.nGenotypes, self.options.maxGenotypes))
            logger.info("Variants = %s" %(variants))
            logger.info("nVariants = %s. nHaplotypes = %s. nGenotypes = %s" %(len(variants), self.nHaplotypes, self.nGenotypes))
            raise StandardError, "Error in cPopulation.setup"

        if not self.nHaplotypes <= self.options.maxHaplotypes:
            logger.error("nHaplotypes > options.maxHaplotypes. (%s > %s). Something is very wrong here. Quitting now." %(self.nHaplotypes, self.options.maxHaplotypes))
            logger.info("Variants = %s" %(variants))
            logger.info("nVariants = %s. nHaplotypes = %s. nGenotypes = %s" %(len(variants), self.nHaplotypes, self.nGenotypes))
            raise StandardError, "Error in cPopulation.setup"

        self.variants = variants
        self.haplotypes = haplotypes
        self.genotypes = genotypes
        self.readBuffers = readBuffers
        cdef Haplotype h = haplotypes[0]
        self.refFile = h.refFile

        if verbosity >=3:
            logger.debug("Population constructed with %s haplotypes, %s genotypes and %s samples" %(self.nHaplotypes, self.nGenotypes, len(self.readBuffers)))

        cdef int genotypeIndex = 0
        cdef int haplotypeIndex = 0
        cdef int individualIndex = 0
        cdef int index = 0
        cdef int hap1Index = 0
        cdef int hap2Index = 0
        cdef double logLikelihood = 0.0
        cdef double likelihood = 0.0
        cdef Haplotype hap1
        cdef Haplotype hap2
        cdef Haplotype hap
        cdef DiploidGenotype genotype

        # For debugging
        cdef Py_ssize_t debugReadIndex = 0
        cdef double debugLike1 = 0.0
        cdef double debugLike2 = 0.0
        cdef double* debugArr1 = NULL
        cdef double* debugArr2 = NULL
        cdef cAlignedRead* debugRead
        cdef Haplotype debugHap1
        cdef Haplotype debugHap2
        cdef double debugLikelihood = 0.0

        cdef dict hapIndDict = {}

        # Cache haplotype indexes by haplotype
        for hapIndex,hap in enumerate(self.haplotypes):
            hapIndDict[hap] = hapIndex

        # Store the haplotype indices in an array
        for genotypeIndex from 0 <= genotypeIndex < self.nGenotypes:
            genotype = genotypes[genotypeIndex]

            if individualIndex == 0:
                hap1 = genotype.hap1
                hap2 = genotype.hap2
                hap1Index = <int>(hapIndDict[hap1])
                hap2Index = <int>(hapIndDict[hap2])
                self.haplotypeIndexes[genotypeIndex][0] = hap1Index
                self.haplotypeIndexes[genotypeIndex][1] = hap2Index

        # Now we compute the genotype likelihoods
        cdef int nReadsThisInd = 0
        cdef bamReadBuffer theBuffer

        for individualIndex from 0 <= individualIndex < self.nIndividuals:

            theBuffer = self.readBuffers[individualIndex]
            nReadsThisInd = theBuffer.reads.windowEnd - theBuffer.reads.windowStart
            self.nReads[individualIndex] = nReadsThisInd
            self.maxLogLikelihoods[individualIndex] = -1e7

            for genotypeIndex from 0 <= genotypeIndex < self.nGenotypes:
                genotype = <DiploidGenotype>genotypes[genotypeIndex]

                if nReadsThisInd == 0:
                    self.genotypeLikelihoods[individualIndex][genotypeIndex] = 1.0
                else:
                    logLikelihood = genotype.calculateDataLikelihood(theBuffer.reads.windowStart, theBuffer.reads.windowEnd, theBuffer.badReads.windowStart, theBuffer.badReads.windowEnd, theBuffer.brokenMates.windowStart, theBuffer.brokenMates.windowEnd, individualIndex, self.nIndividuals, self.goodnessOfFitValues[genotypeIndex])

                    if logLikelihood > self.maxLogLikelihoods[individualIndex]:
                        self.maxLogLikelihoods[individualIndex] = logLikelihood

                    self.genotypeLikelihoods[individualIndex][genotypeIndex] = logLikelihood

        # Re-scale all genotype likelihoods using the maximum value
        for individualIndex from 0 <= individualIndex < self.nIndividuals:
            for genotypeIndex from 0 <= genotypeIndex < self.nGenotypes:
                if self.nReads[individualIndex] != 0:
                    self.genotypeLikelihoods[individualIndex][genotypeIndex] = max(1e-300, exp(self.genotypeLikelihoods[individualIndex][genotypeIndex] - self.maxLogLikelihoods[individualIndex]))
                else:
                    self.genotypeLikelihoods[individualIndex][genotypeIndex] = 1.0

        # DEBUG OUTPUT
        if verbosity >=3:
            logger.debug("Printing all haplotypes in population and their sequences...")

            for index,hap1 in enumerate(self.haplotypes):
                logger.debug("%s\t%s" %(index, hap1))

            for index,hap1 in enumerate(self.haplotypes):
                logger.debug("%s\t%s" %(index, hap1.haplotypeSequence))

            logger.debug("Done printing all haplotypes in population...")

        if verbosity >= 4:

            logger.debug("")
            logger.debug("####################################################################")
            logger.debug("Read alignment Likelihood debug information (top 10 genotypes)")
            logger.debug("####################################################################")
            logger.debug("")
            logger.debug("Sample\tPhred-likelihood\tNumber of reads\tGenotype")

            for individualIndex from 0 <= individualIndex < self.nIndividuals:

                topGenotypes = []
                theBuffer = self.readBuffers[individualIndex]
                nReadsThisInd = theBuffer.reads.windowEnd - theBuffer.reads.windowStart

                for genotypeIndex from 0 <= genotypeIndex < self.nGenotypes:
                    genotype = genotypes[genotypeIndex]

                    if self.genotypeLikelihoods[individualIndex][genotypeIndex] > 1e-300:
                        logLikelihood = log(self.genotypeLikelihoods[individualIndex][genotypeIndex])
                    else:
                        logLikelihood = 1e7*mLTOT

                    topGenotypes.append( (int(0.5+(logLikelihood/mLTOT)), genotype) )

                topGenotypes = sorted(topGenotypes)[0:10]

                for phred,genotype in topGenotypes:
                    logger.debug("%s\t%s\t%s\t%s" %(theBuffer.sample, phred, nReadsThisInd, genotype))

                    # Extremely verbose output. Print liklelihoods, for each genotype, haplotype and read
                    if verbosity >= 5:

                        debugHap1 = genotype.hap1
                        debugHap2 = genotype.hap2
                        debugArr1 = debugHap1.alignReads(individualIndex, theBuffer.reads.windowStart, theBuffer.reads.windowEnd, theBuffer.badReads.windowStart, theBuffer.badReads.windowEnd, theBuffer.brokenMates.windowStart, theBuffer.brokenMates.windowEnd, False)
                        debugArr2 = debugHap2.alignReads(individualIndex, theBuffer.reads.windowStart, theBuffer.reads.windowEnd, theBuffer.badReads.windowStart, theBuffer.badReads.windowEnd, theBuffer.brokenMates.windowStart, theBuffer.brokenMates.windowEnd, False)

                        logger.debug("Verbose output: logging likelihoods for each read...")
                        logger.debug("")
                        logger.debug("For Genotype %s" %(genotype))
                        logger.debug("Hap1 start = %s. end = %s. Hap1 start = %s. end = %s." %(genotype.hap1.startPos,genotype.hap1.endPos,genotype.hap2.startPos,genotype.hap2.endPos))
                        logger.debug("Logging haplotype sequences...")
                        logger.debug(str(debugHap1.getMutatedSequence())[50:-50])
                        logger.debug(str(debugHap2.getMutatedSequence())[50:-50])
                        logger.debug("")
                        logger.debug("Sample\tRead\tLL Hap1\tLL Hap2\tGL So Far\tRead MapQ\tRead Start\tRead End")

                        debugLikelihood = 0.0

                        for debugReadIndex from 0 <= debugReadIndex < nReadsThisInd:
                            debugLike1 = debugArr1[debugReadIndex]
                            debugLike2 = debugArr2[debugReadIndex]
                            debugRead = theBuffer.reads.windowStart[debugReadIndex]
                            debugLikelihood += log(0.5*(exp(debugLike1) + exp(debugLike2)))
                            logger.debug("%s\t%s\t%1.2f\t%1.2f\t%1.2f\t%s\t%s\t%s" %(theBuffer.sample, debugReadIndex, -10*debugLike1, -10*debugLike2, debugLikelihood, debugRead.mapq, debugRead.pos, debugRead.end))

                        logger.debug("")
                        logger.debug("#####################################################################################################")
                        logger.debug("")

    cdef double EMiteration(self, double* freq, double* newFreqs):
        """
        Perform one EM update, given initial frequency freq. Returns the updated frequency
        """
        cdef double* genotypeLikelihoodsThisIndividual
        cdef double* csr
        cdef double maxChange = 0.0
        cdef double csrSum = 0.0
        cdef double tmp = 0.0
        cdef double thisCSR = 0.0
        cdef double oldFreq = 0.0
        cdef double freqChange = 0.0
        cdef double newFreq = 0.0
        cdef double likelihoodThisGenotype = 0.0
        cdef int i = 0
        cdef int j = 0
        cdef int k = 0
        cdef int s = 0
        cdef int r = 0
        cdef int nIndWithData = 0

        for i from 0 <= i < self.nIndividuals:

            if self.nReads[i] == 0:
                continue

            genotypeLikelihoodsThisIndividual = self.genotypeLikelihoods[i]
            csr = self.EMLikelihoods[i]
            csrSum = 0.0
            nIndWithData += 1

            for j from 0 <= j < self.nGenotypes:

                likelihoodThisGenotype = genotypeLikelihoodsThisIndividual[j]
                s = self.haplotypeIndexes[j][0]
                r = self.haplotypeIndexes[j][1]
                thisCSR = likelihoodThisGenotype*freq[s]*freq[r]*(1 + (r != s))
                csr[j] = thisCSR
                csrSum += thisCSR

            # Normalisation
            if csrSum > 0.0:
                for j from 0 <= j < self.nGenotypes:
                    csr[j] /= csrSum

        # New, better implementation. Linear in nIndividuals and nGenotypes
        # Set all temp frequencies to 0
        for k from 0 <= k < self.nHaplotypes:
            newFreqs[k] = 0.0

        for i from 0 <= i < self.nIndividuals:

            if self.nReads[i] == 0:
                continue

            csr = self.EMLikelihoods[i]

            for j from 0 <= j < self.nGenotypes:
                thisCSR = csr[j]
                s = self.haplotypeIndexes[j][0]
                r = self.haplotypeIndexes[j][1]
                newFreqs[s] += thisCSR
                newFreqs[r] += thisCSR

        for k from 0 <= k < self.nHaplotypes:
            newFreqs[k] = newFreqs[k]/(2*nIndWithData)
            freqChange = fabs(freq[k] - newFreqs[k])

            if freqChange > maxChange:
                maxChange = freqChange

            freq[k] = newFreqs[k]

        return maxChange

    cdef double calculatePosterior(self, Variant var, int flatPrior=0):
        """
        Calculate for ecah variant the posterior probability that that variant (in standard format)
        is segregating in the population
        """
        cdef Haplotype hap
        cdef tuple vsf
        cdef double* genotypeLikelihoodsThisIndividual
        cdef double* freqsPrimeByHapIndex = self.freqsPrimeByHapIndex
        cdef double prior = 0.0

        if flatPrior == 1:
            prior = 0.5
        else:
            prior = var.calculatePrior( self.refFile )

        cdef double totalProb = 0.0
        cdef double sumFreqs = 0.0
        cdef double thisFreq = 0.0
        cdef double ratio = 0.0
        cdef double likelihoodThisGenotype = 0.0
        cdef double variantPosterior = 0.0
        cdef double phredPosterior = 0.0
        cdef double sumProbVariantThisIndividual = 0.0
        cdef double sumProbNoVariantThisIndividual = 0.0
        cdef double sumLogProbVariant = 0.0
        cdef double sumLogProbNoVariant = 0.0
        cdef double r_prime = 0.0
        cdef double s_prime = 0.0
        cdef double factor = 0.0
        cdef int i = 0
        cdef int j = 0
        cdef int r = 0
        cdef int s = 0
        cdef int nReadsThisInd = 0
        cdef int logOfMinFloat = -708
        cdef bamReadBuffer theBuffer

        if self.verbosity >= 3:
            logger.debug("")
            logger.debug("#########################################################################")
            logger.debug("Posterior calculation debug information")
            logger.debug("#########################################################################")
            logger.debug("")
            logger.debug("Computing posterior for variants %s. N haplotypes = %s. n Ind = %s" %(var, self.nHaplotypes, self.nIndividuals))

        # Re-scale haplotype frequencies for those not containing this variant, and
        # store them.

        for i from 0 <= i < self.nHaplotypes:

            hap = self.haplotypes[i]
            vsf = hap.variants

            if var not in vsf:
                freqsPrimeByHapIndex[i] = self.frequencies[i]
                sumFreqs += self.frequencies[i]
            else:
                freqsPrimeByHapIndex[i] = 0.0

        if self.verbosity >= 3:
            logger.debug("Sum of frequencies of haplotypes containing variant %s = %s" %(var, sumFreqs))

        if self.verbosity >= 3:
            logger.debug("Haplotype\tUn-scaled freq\tScaled freq")

        # Avoid divide by zero errors
        if sumFreqs > 0:

            for i from 0 <= i < self.nHaplotypes:

                if self.verbosity >= 3:
                    logger.debug("%s\t%s\t%s" %(self.haplotypes[i], self.frequencies[i], freqsPrimeByHapIndex[i]))

                freqsPrimeByHapIndex[i] /= sumFreqs

        # Compute the variant posterior from the frequencies, and the per-sample
        # likelihood values for all haplotypes.

        if self.verbosity >= 4:
            logger.debug("")
            logger.debug("Sample\tLikelihood\tfreqHap1\tfreqHap2\tscaledFreq1\tscaledFreq2\tsumVar\tsumNoVar\tGenotype")

        for i from 0 <= i < self.nIndividuals:

            nReadsThisInd = self.nReads[i]

            if nReadsThisInd == 0:
                continue

            genotypeLikelihoodsThisIndividual = self.genotypeLikelihoods[i]
            sumProbVariantThisIndividual = 0.0
            sumProbNoVariantThisIndividual = 0.0

            for j from 0 <= j < self.nGenotypes:

                likelihoodThisGenotype = genotypeLikelihoodsThisIndividual[j]

                r = self.haplotypeIndexes[j][0]
                s = self.haplotypeIndexes[j][1]

                factor = 1.0

                if r != s:
                    factor = 2.0

                sumProbVariantThisIndividual += (factor * self.frequencies[r] * self.frequencies[s] * likelihoodThisGenotype)

                # r_prime and s_prime will be 0.0 for all haplotypes containing the variant
                r_prime = freqsPrimeByHapIndex[r]
                s_prime = freqsPrimeByHapIndex[s]
                sumProbNoVariantThisIndividual += (factor * r_prime * s_prime * likelihoodThisGenotype)

                if self.verbosity >= 4:
                    theBuffer = self.readBuffers[i]
                    logger.debug("%s\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%s" %(theBuffer.sample,likelihoodThisGenotype,self.frequencies[r],self.frequencies[s],r_prime,s_prime,sumProbVariantThisIndividual,sumProbNoVariantThisIndividual,self.genotypes[j]))

            if sumProbVariantThisIndividual > 0:
                sumLogProbVariant += log(sumProbVariantThisIndividual)
            else:
                sumLogProbVariant += logOfMinFloat

            if sumProbNoVariantThisIndividual > 0:
                sumLogProbNoVariant += log(sumProbNoVariantThisIndividual)
            else:
                sumLogProbNoVariant += logOfMinFloat

        ratio = max(1e-300, exp(sumLogProbNoVariant - sumLogProbVariant))

        if self.verbosity >= 3:
            variantPosterior = prior / (prior + ratio * (1.0 - prior))
            logger.debug("For variant %s, post = %.20f. 1-post = %s. max(1e-20, 1-post) = %s" %(var, variantPosterior, 1.0-variantPosterior, max(1e-20, 1.0-variantPosterior)))
            logger.debug("For variant %s, Prior = %s. Ratio = %s. Post = %s. SumLogProbNoVar = %s. SumLogProbVar = %s" %(var,prior,ratio,variantPosterior,sumLogProbNoVariant,sumLogProbVariant))
            logger.debug("Variant phred score = %s" %(min(200, round(-10.0 * (log10(ratio*(1.0-prior)) - log10(prior + ratio*(1.0 - prior)))))))

        return round(-10.0 * (log10(ratio*(1.0-prior)) - log10(prior + ratio*(1.0 - prior))))

    cdef void computeVariantPosteriors(self):
        """
        """
        cdef Haplotype haplotype
        cdef Variant v
        cdef set done = set()
        cdef double posterior = 0.0

        for haplotype in self.haplotypes:

            for v in haplotype.variants:

                if v in done:
                    continue

                posterior = self.calculatePosterior(v)

                if posterior >= self.options.minPosterior:
                    self.variantPosteriors[v] = posterior

                    try:
                        self.varsByPos[v.refPos].append(v)
                    except:
                        self.varsByPos[v.refPos] = [v]

                done.add(v)

    cdef void callGenotypes(self):
        """
        Identify best genotype calls
        """
        cdef:
            int index = 0
            int indexOfBestGenotype = -1
            double maxEMLikelihood = 0.0
            double likelihoodThisGenotype = 0.0
            DiploidGenotype gt
            bamReadBuffer theBuffer

        for index from 0 <= index < self.nIndividuals:

            if self.nReads[index] == 0:
                self.genotypeCalls.append(None)
                continue

            indexOfBestGenotype = -1
            maxEMLikelihood = 0.0
            likelihoodThisGenotype = 0.0

            if self.verbosity >= 4:
                logger.debug("")
                logger.debug("####################################################################")
                logger.debug("EM Likelihood debug information")
                logger.debug("####################################################################")
                logger.debug("")
                logger.debug("Sample\tEM Likelihood\tGenotype")
                logger.debug("useEM = %s" %(self.useEMLikelihoods))

            for genotypeIndex from 0 <= genotypeIndex < self.nGenotypes:

                if self.useEMLikelihoods == 1:
                    likelihoodThisGenotype = self.EMLikelihoods[index][genotypeIndex]
                else:
                    likelihoodThisGenotype = self.genotypeLikelihoods[index][genotypeIndex]

                # Debug information
                if self.verbosity >= 4:
                    theBuffer = self.readBuffers[index]
                    logger.debug("%s\t%s\t%s" %(theBuffer.sample,likelihoodThisGenotype,self.genotypes[genotypeIndex]))

                if indexOfBestGenotype == -1 or likelihoodThisGenotype > maxEMLikelihood:
                    maxEMLikelihood = likelihoodThisGenotype
                    indexOfBestGenotype = genotypeIndex

            gt = self.genotypes[indexOfBestGenotype]

            if self.verbosity >= 3:
                theBuffer = self.readBuffers[index]
                logger.debug("Called genotype %s for indiviual %s" %(gt, theBuffer.sample))

            self.genotypeCalls.append(gt)

    cdef void call(self, int maxIters, int computeVCFFields):
        """
        Estimate the haplotype frequencies and genotypes for this population and locus.
        Terminate the algorithm once the maximum change across all components of the
        frequency is < eps, or after maxIters iterations.
        """
        cdef double eps = min(1e-3, 1.0 / (self.nIndividuals*2*2))
        cdef double maxChange = eps + 1
        cdef double uniformFreq = 1.0/self.nHaplotypes
        cdef int iters = 0
        cdef int index = 0

        # Uniform initial frequency
        for index from 0 <= index < self.nHaplotypes:
            self.frequencies[index] = uniformFreq

        if self.verbosity >= 4:
            logger.debug("")
            logger.debug("####################################################################")
            logger.debug("EM Iteration debug information")
            logger.debug("####################################################################")
            logger.debug("")

        while maxChange > eps and iters < maxIters:
            maxChange = self.EMiteration(self.frequencies, self.newFrequencies)
            iters += 1

            # More debug information
            if self.verbosity >= 4:
                logger.debug("Done %s EM iterations. Max change = %s" %(iters, maxChange))
                logger.debug("Haplotype\tFrequency")
                for k from 0 <= k < self.nHaplotypes:
                    logger.debug("%s\t%s" %(self.haplotypes[k], self.frequencies[k]))

        if self.verbosity >= 3:
            logger.debug("EM Stats: Final Max Frequency Change = %s. nIterations = %s" %(maxChange, iters))

        self.callGenotypes()
        self.computeVariantPosteriors()

        if computeVCFFields != 0 and len(self.variantPosteriors.keys()) > 0:
            self.computeVariantINFO()
            self.computeVariantFILTER()

        if self.options.alignScoreFile == "":
            return

        # Print alignment scores of reads against haplotypes to file
        cdef bamReadBuffer theBuffer
        cdef double* arr
        cdef Haplotype hap
        cdef cAlignedRead** pStartRead
        cdef cAlignedRead** pEndRead

        fo = open(self.options.alignScoreFile, "a")

        for individualIndex from 0 <= individualIndex < self.nIndividuals:
            theBuffer = self.readBuffers[individualIndex]
            nReadsThisInd = theBuffer.reads.windowEnd - theBuffer.reads.windowStart
            windowStart =  theBuffer.startBase
            windowEnd   =  theBuffer.endBase
            fo.write("Individual\t%d\t%d\t%d\n" %(individualIndex, len(self.haplotypes), nReadsThisInd))

            for hapIdx, hap in enumerate(self.haplotypes):
                thisStr = hap.getMutatedSequence()
                thisLen = len(thisStr)
                fo.write("%d %d %s %f\n" %(hap.startPos, hap.endPos, thisStr[hap.endBufferSize:thisLen-hap.endBufferSize+1], self.frequencies[hapIdx]))

            for hap in self.haplotypes:
                arr = hap.alignReads(individualIndex, theBuffer.reads.windowStart, theBuffer.reads.windowEnd, theBuffer.badReads.windowStart, theBuffer.badReads.windowEnd, theBuffer.brokenMates.windowStart, theBuffer.brokenMates.windowEnd, self.verbosity)
                arrP = []
                for readIndex from 0 <= readIndex < nReadsThisInd:
                    arrP.append("%1.3E" %(-10*arr[readIndex]))
                fo.write("%s\n" %("\t".join(arrP)))
            mappingQual = []
            pStartRead = theBuffer.reads.windowStart
            pEndRead = theBuffer.reads.windowEnd
            while pStartRead != pEndRead:
                pRead = pStartRead[0]
                mappingQual.append(pRead.mapq)
                pStartRead +=1
            fo.write("%s\n" %("\t".join(map(str,mappingQual))))
        fo.close()

###################################################################################################
