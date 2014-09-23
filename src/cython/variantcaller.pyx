"""
This module contains the top-level variant-calling
code, including a VariantCaller class.
"""
from __future__ import division

import logging
import multiprocessing
import window
import vcf
import variantutils
import time
import datetime
import platypusutils
import vcfutils
import chaplotype
import sys

cimport platypusutils
cimport vcfutils
cimport fastafile
cimport samtoolsWrapper
cimport cwindow
cimport cpopulation
cimport variant
cimport variantFilter
cimport chaplotype

from chaplotype cimport Haplotype
from cgenotype cimport generateAllGenotypesFromHaplotypeList,DiploidGenotype
from cwindow cimport bamReadBuffer
from samtoolsWrapper cimport Samfile,cAlignedRead
from cpopulation cimport Population
from variant cimport Variant
from variant cimport VariantCandidateGenerator
from fastafile cimport FastaFile
from variantFilter cimport filterVariants
from variantFilter cimport filterVariantsInWindow
from variantFilter cimport filterVariantsByCoverage
from variantFilter cimport padVariants
from variantFilter cimport computeVariantReadSupportFrac
from vcfutils cimport outputCallToVCF,refAndAlt
from platypusutils cimport loadBAMData
from assembler cimport assembleReadsAndDetectVariants
from bisect import bisect
from platypusutils cimport leftNormaliseIndel
from platypusutils cimport betaBinomialCDF
from platypusutils import open

logger = logging.getLogger("Log")
vcfHeader = [('fileDate',datetime.date.fromtimestamp(time.time())), ('source','Platypus_Version_%s' %(platypusutils.PLATYPUS_VERSION))]

###################################################################################################

cdef extern from "math.h":
    double exp(double)
    double log(double)
    double log10(double)
    double fabs(double)
    double round(double)

###################################################################################################

cdef extern from "string.h":
  ctypedef int size_t
  void *memcpy(void *dst,void *src,size_t len)
  int strncmp(char *s1,char *s2,size_t len)
  char *strncpy(char *dest,char *src, size_t len)
  size_t strlen(char *s)
  int memcmp( void * s1, void *s2, size_t len )

###################################################################################################

cdef void callVariantsInWindow(dict window, options, FastaFile refFile, list readBuffers, Population pop, int start, int end, char* refSeq) except *:
    """
    """
    cdef Haplotype hap
    cdef Variant v
    cdef bamReadBuffer theReadBuffer
    cdef bytes chrom = window['chromosome']
    cdef list variants = window["variants"]
    cdef int windowStart = window['startPos']
    cdef int windowEnd = window['endPos']
    cdef int nReadsThisWindow = 0
    cdef int qualBinSize = options.qualBinSize

    cdef Haplotype refHaplotype = Haplotype(chrom, windowStart, windowEnd, (), refFile, options.rlen, options)

    pop.reset()

    for theReadBuffer in readBuffers:
        theReadBuffer.setWindowPointers(windowStart, windowEnd, start, end, refSeq, qualBinSize)
        nReadsThisWindow += (theReadBuffer.reads.windowEnd - theReadBuffer.reads.windowStart)

    if nReadsThisWindow == 0 and not options.outputRefCalls:
        if options.verbosity >= 3:
            logger.debug("No coverage in window %s:%s-%s. Skipping" %(chrom, windowStart, windowEnd))
        return

    if nReadsThisWindow > options.maxReads:
        logger.debug("Skipping pathalogical window %s:%s-%s with %s reads" %(chrom, windowStart, windowEnd, nReadsThisWindow))
        return

    if len(variants) > options.maxVariants:

        if options.skipDifficultWindows:
            logger.debug("Skipping window with %s variants" %(len(variants)))
            return

        elif options.filterVarsByCoverage:
            #logger.debug("There are %s variants before filtering" %(len(window['variants'])))
            filterVariantsByCoverage(window, chrom, windowStart, windowEnd, refFile, options, variants, refHaplotype, readBuffers)
            #logger.debug("There are %s variants after filtering" %(len(window['variants'])))

    # Create haplotype list using all data from all samples. Always consider the reference haplotype.
    cdef list allVarHaplotypes = variantFilter.getHaplotypesInWindow(window, nReadsThisWindow, refFile, options.maxReads, options.minMapQual, options.minBaseQual, options.maxHaplotypes, options.maxVariants, options.rlen, options.verbosity, readBuffers, options)
    #cdef list allUniqueHaplotypes = list(set([refHaplotype] + allVarHaplotypes))
    cdef list allUniqueHaplotypes = mergeHaplotypes([refHaplotype] + allVarHaplotypes, refFile)
    cdef list allGenotypes = generateAllGenotypesFromHaplotypeList(allUniqueHaplotypes)

    cdef int nUniqueHaplotypes = len(allUniqueHaplotypes)

    if nUniqueHaplotypes <= 1:

        # This occurs if we found mutation candidates, but the sequences of the mutated haplotypes are all the same as the
        # reference.
        if not options.outputRefCalls:

            #logger.debug("No non-ref haplotypes to check in window %s. All var haps = %s. all Unique haps = %s" %(window, allVarHaplotypes, allUniqueHaplotypes))
            #logger.debug("REF Then VAR")
            #logger.debug(refHaplotype.haplotypeSequence)

            #for hap in allVarHaplotypes:
            #    logger.debug(hap.haplotypeSequence)

            return

    pop.setup(variants, allUniqueHaplotypes, allGenotypes, options.nInd, options.verbosity, readBuffers)

    cdef int maxEMIterations = 100
    pop.call(maxEMIterations, 1)

###################################################################################################

def countTotalReadsInRegion(list readBuffers):
    """
    Count and return the total number of reads (good, bad, broken)
    loaded into the read buffers.
    """
    cdef bamReadBuffer theReadBuffer
    totalReads = 0
    cigarSize = 0
    seqQualSize = 0
    readSize = 0
    pointerSize = 0
    totalBufferSize = 0
    totalSeqQualSize = 0

    for theReadBuffer in readBuffers:
        thisBufferSize = theReadBuffer.reads.getSize()
        thisBadBufferSize = theReadBuffer.badReads.getSize()
        totalReads += (theReadBuffer.reads.getSize())
        totalReads += (theReadBuffer.badReads.getSize())
        totalReads += (theReadBuffer.brokenMates.getSize())
        totalBufferSize += (theReadBuffer.reads.__capacity)
        totalBufferSize += (theReadBuffer.badReads.__capacity)
        totalBufferSize += (theReadBuffer.brokenMates.__capacity)

        for i in range(thisBufferSize):
            totalSeqQualSize += (strlen(theReadBuffer.reads.array[i].seq) + 1)
            totalSeqQualSize += (strlen(theReadBuffer.reads.array[i].qual) + 1)

            totalSeqQualSize += 2*theReadBuffer.reads.array[i].cigarLen

            seqQualSize += (strlen(theReadBuffer.reads.array[i].seq) + 1)
            seqQualSize += (strlen(theReadBuffer.reads.array[i].qual) + 1)

            cigarSize += 2*theReadBuffer.reads.array[i].cigarLen*sizeof(theReadBuffer.reads.array[0].cigarOps[0])

        for i in range(thisBadBufferSize):
            totalSeqQualSize += (strlen(theReadBuffer.badReads.array[i].seq) + 1)
            totalSeqQualSize += (strlen(theReadBuffer.badReads.array[i].qual) + 1)

            totalSeqQualSize += 2*theReadBuffer.badReads.array[i].cigarLen

            seqQualSize += (strlen(theReadBuffer.badReads.array[i].seq) + 1)
            seqQualSize += (strlen(theReadBuffer.badReads.array[i].qual) + 1)

            cigarSize += 2*theReadBuffer.badReads.array[i].cigarLen*sizeof(theReadBuffer.reads.array[0].cigarOps[0])

        # Size of reads
        totalSeqQualSize += sizeof(theReadBuffer.reads.array[0][0])*thisBufferSize
        totalSeqQualSize += sizeof(theReadBuffer.badReads.array[0][0])*thisBadBufferSize

        readSize += sizeof(theReadBuffer.reads.array[0][0])*thisBufferSize
        readSize += sizeof(theReadBuffer.badReads.array[0][0])*thisBadBufferSize

        # Size of pointers to reads
        totalSeqQualSize += sizeof(theReadBuffer.reads.array[0])*theReadBuffer.reads.__capacity
        totalSeqQualSize += sizeof(theReadBuffer.reads.array[0])*theReadBuffer.badReads.__capacity

        pointerSize += sizeof(theReadBuffer.reads.array[0])*theReadBuffer.reads.__capacity
        pointerSize += sizeof(theReadBuffer.reads.array[0])*theReadBuffer.badReads.__capacity

    totalSize = seqQualSize + cigarSize + readSize + pointerSize
    logger.debug("Sizes: SeqQual = %s. cigar = %s. read = %s. pointers = %s. Total = %s" %(seqQualSize, cigarSize, readSize, pointerSize, totalSize))

    return totalReads,totalBufferSize,totalSeqQualSize

###################################################################################################

cdef int doWeNeedToAssembleThisRegion(list readBuffers, bytes chrom, int start, int end, options, char* refSeq):
    """
    Scan through the specified region, and report whether or not we need to run the assember to
    assemble it. This decision will be based on a number of metrics: insert size distribution,
    mapping qualities, numbers of mis-alignments/mismatches etc.
    """
    cdef double totalReads = 0.0
    cdef double readsThisSample = 0.0
    cdef double badReadsThisSample = 0.0
    cdef double totalMismatches = 0.0
    cdef double totalAlignmentGaps = 0.0
    cdef double nImproperReads = 0.0
    cdef int qualBinSize = options.qualBinSize

    cdef bamReadBuffer theReadBuffer

    for theReadBuffer in readBuffers:
        theReadBuffer.setWindowPointers(start, end, start, end, refSeq, qualBinSize)

    if options.assembleAll:
        return True

    # Check all samples, but as soon as we see one which triggers the assembly we
    # return, so we assemble everything is even one of the samples triggers.
    for theReadBuffer in readBuffers:

        totalReads += (theReadBuffer.reads.windowEnd - theReadBuffer.reads.windowStart)
        readsThisSample = (theReadBuffer.reads.windowEnd - theReadBuffer.reads.windowStart)
        badReadsThisSample = (theReadBuffer.badReads.windowEnd - theReadBuffer.badReads.windowStart)

        if readsThisSample == 0:
            continue

        totalAlignmentGaps = theReadBuffer.countAlignmentGaps()
        nImproperReads = theReadBuffer.countImproperPairs()

        if totalAlignmentGaps / readsThisSample > 2:
            logger.info("Need to assemble %s:%s-%s. nReads = %s. nGaps = %s. n per read = %s" %(chrom, start, end, readsThisSample, totalAlignmentGaps, totalAlignmentGaps/readsThisSample))
            return 1

        if nImproperReads / (readsThisSample + badReadsThisSample) > 0.1:
            logger.info("Need to assemble %s:%s-%s. nReads = %s. nImproper = %s. fraction = %s" %(chrom, start, end, readsThisSample, nImproperReads, nImproperReads/(readsThisSample+badReadsThisSample)))
            return 1

    # No signal for assembly
    return 0

###################################################################################################

cdef list mergeHaplotypes(list haplotypes, FastaFile refFile):
    """
    Merge haplotypes by sequence. If we see 2 identical haplotypes, then take the one with
    the fewest variants, i.e. we want to report complex long haplotypes as single variants, if
    possible.
    """
    cdef list sortedHaps = sorted(haplotypes)
    cdef list mergedHaplotypes = []

    cdef Haplotype hap = None
    cdef Haplotype lastHap = None

    cdef double priorOne = 1.0
    cdef double priorTwo = 1.0
    cdef Variant theVar

    #logger.debug("Printing un-merged haplotypes")

    for hap in sortedHaps:
        #logger.debug(hap)

        # Start.
        if lastHap is None:
            lastHap = hap

        # Same as last hap: take hap with least variants
        elif hap == lastHap:

            logger.debug("Merging haplotypes %s and %s" %(hap, lastHap))

            # Take the haplotype with the set of variants
            # which have the lowest combined prior.
            priorOne = 1.0
            priorTwo = 1.0

            for theVar in lastHap.variants:
                priorOne *= theVar.calculatePrior(refFile)

            for theVar in hap.variants:
                priorTwo *= theVar.calculatePrior(refFile)

            # lastHap has better prior. Keep that one
            if priorOne > priorTwo:
                pass

            # hap has better prior. Keep that one
            elif priorTwo > priorOne:
                lastHap = hap

            # Priors are close or the same. Doesn't really matter
            else:
                pass

        # Add whichever hap we kept to merged list
        else:
            mergedHaplotypes.append(lastHap)
            lastHap = hap

    if lastHap is not None:
        mergedHaplotypes.append(lastHap)

    #logger.debug("")
    #logger.debug("Printing merged haplotypes...")
    #for hap in mergedHaplotypes:
    #    logger.debug(hap)
    #logger.debug("")

    return mergedHaplotypes

###################################################################################################

#cdef int computeCoverageThreshold(list readBuffers, options):
#    """
#    Compute and return an estimate
#    """
#    cdef bamReadBuffer theReadBuffer
#
#    for theReadBuffer in readBuffers:
#        if theReadBuffer.reads.getSize() > 0:
#            if minSampleCoverage == -1:
#                minSampleCoverage = (theReadBuffer.reads.getSize() * theReadBuffer.reads.array[0].rlen) / (end-start)
#            else:
#                minSampleCoverage = min(minSampleCoverage, (theReadBuffer.reads.getSize() * theReadBuffer.reads.array[0].rlen) / (end-start))
#
#    cdef int coverageThreshold = max(options.minReads, options.minVarFreq*minSampleCoverage)
#    logger.debug("The minimum mean coverage in this region is %s. Small variants with coverage < %s will be removed" %(minSampleCoverage, coverageThreshold))

###################################################################################################

cdef list generateVariantsInRegion(bytes chrom, int start, int end, bamFiles, FastaFile refFile, options, windowGenerator, outputFile, vcf, list readBuffers):
    """
    Generate a list of variants, based on BAM file data and any input VCFs, in the specified region,
    and return these.
    """
    cdef VariantCandidateGenerator varCandGen
    cdef bamReadBuffer theReadBuffer
    cdef list rawBamVariants = []
    cdef list vcfFileVariants = []
    cdef list assemblerVariants = []
    cdef bytes refSequenceBytes = refFile.getSequence(chrom,start,end)
    cdef char* refSequence = refSequenceBytes
    cdef Variant v
    cdef int maxReadLength = options.rlen
    cdef int longestRead = 0
    cdef int minSampleCoverage = -1
    cdef double minVarFreq = options.minVarFreq

    if options.verbosity >= 3:
        totalReadsInRegion,totalBufferSizeInRegion,totalSeqQualSize = countTotalReadsInRegion(readBuffers)
        logger.debug("There are %s reads (buffer size = %s. Total reads size = %s bytes) in the region %s:%s-%s" %(totalReadsInRegion, totalBufferSizeInRegion, totalSeqQualSize, chrom, start, end))


    # Generate initial variant candidate list from BAM files. Do this always, so we
    # get read-counts for variants in the BAMs
    if options.getVariantsFromBAMs:

        allSampleVarCandGen = VariantCandidateGenerator((chrom,start,end), refFile, options.minMapQual, options.minFlank, options.minBaseQual, options.maxReads, maxReadLength, options, options.verbosity, options.genSNPs, options.genIndels)

        # We want to check all BAMs and accumulate data before filtering
        for theReadBuffer in readBuffers:

            # Make generator per sample
            varCandGen = VariantCandidateGenerator((chrom,start,end), refFile, options.minMapQual, options.minFlank, options.minBaseQual, options.maxReads, maxReadLength, options, options.verbosity, options.genSNPs, options.genIndels)

            if theReadBuffer.reads.getLengthOfLongestRead() > longestRead:
                longestRead = theReadBuffer.reads.getLengthOfLongestRead()

            varCandGen.addCandidatesFromReads(theReadBuffer.reads.array, theReadBuffer.reads.array + theReadBuffer.reads.getSize())
            #varCandGen.addCandidatesFromReads(theReadBuffer.badReads.array, theReadBuffer.badReads.array + theReadBuffer.badReads.getSize())

            if options.verbosity >= 3:
                logger.debug("Processed sample %s. Detected %s unfiltered variant candidates so far" %(theReadBuffer.sample, len(varCandGen.variantHeap)))

            ## Filter by per-sample coverage on each file. No variant should be kept unless it has >= coverageThreshold supporting reads
            ## in at least one sample
            for key,v in varCandGen.variantHeap.iteritems():
                varReadFrac = computeVariantReadSupportFrac(v, theReadBuffer)
                #logger.info("Var frac for variant %s is %s" %(v, varReadFrac))

                # Add to total list, merging if necessary
                if varReadFrac >= minVarFreq:
                    allSampleVarCandGen.addVariantToList(v)
                # Always add indels
                elif v.nAdded != v.nRemoved:
                    allSampleVarCandGen.addVariantToList(v)
                else:
                    continue

            if options.verbosity >= 3:
                logger.debug("This reduces to %s candidates after coverage-filtering" %(len(varCandGen.variantHeap)))

        # Add candidates from all samples to list
        rawBamVariants.extend(allSampleVarCandGen.getCandidates(0))

        # options.rlen over-rides the data, if set, otherwise we use the longest read
        # in the data
        if longestRead > 0:
            if longestRead >= options.maxSize:
                logger.warning("Found very long read (%s bases). Capping max read length at --maxSize (%s)." %(longestRead, options.maxSize))
                logger.warning("Longer reads will be used for alignments but not to determine window boundaries")
                logger.warning("Increase --maxSize if you want longer reads to be used for determining widow boundaries")
                logger.warning("But keep --maxSize below 5000 otherwise problems will occur downstream and may lead to crashes")
                options.rlen = options.maxSize
            else:
                options.rlen = longestRead

    maxReadLength = options.rlen

    # Get variants from external VCF source file
    if options.sourceFile:
        variantSource = variantutils.VariantCandidateReader(options.sourceFile, options)
        vcfFileVariants.extend(variantSource.Variants(chrom,start,end))

    # Get variants from assembler
    if options.assemble:

        assemRegionShift = max(100, min(1000, options.assemblyRegionSize//2))
        #assemRegionShift = max(100, min(1000, options.assemblyRegionSize//4))
        logger.debug("Assembly region windows will be tiled. Each region will start %s bases after the last one" %(assemRegionShift))

        for assemIndex,assemStart in enumerate(xrange(start, end, assemRegionShift)):
        #for assemIndex,assemStart in enumerate(xrange(start, end, options.assemblyRegionSize//2)):
        #for assemIndex,assemStart in enumerate(xrange(start, end, options.assemblyRegionSize)):

            assemEnd = min(assemStart + options.assemblyRegionSize, end)
            refStart = max(0, assemStart-options.assemblyRegionSize)
            refEnd = assemEnd+options.assemblyRegionSize
            refSequenceBytes = refFile.getSequence(chrom, refStart, refEnd)
            refSequence = refSequenceBytes

            needToAssemble = doWeNeedToAssembleThisRegion(readBuffers, chrom, assemStart, assemEnd, options, refSequence)

            if not needToAssemble:
                continue

            #logger.debug("Assembling region %s:%s-%s (ref seq %s:%s-%s)" %(chrom, assemStart, assemEnd, chrom, refStart, refEnd))
            assemblerVariants.extend(assembleReadsAndDetectVariants(chrom, assemStart, assemEnd, refStart, refEnd, readBuffers, refSequence, options))
            #logger.debug("Variants Found by assembler in region are %s" %(assemblerVariants))

    cdef list allVarCands = rawBamVariants + vcfFileVariants + assemblerVariants
    cdef list leftNormVars = sorted( [leftNormaliseIndel(v, refFile, maxReadLength) for v in allVarCands] )
    cdef list filteredVariants = filterVariants(leftNormVars, refFile, maxReadLength, options.minReads, options.maxSize, options.verbosity, options)

    # Need to work out how to consistently report these.
    logger.debug("There are %s filtered variant candidates in reads which overlap the region %s:%s-%s" %(len(filteredVariants), chrom, start, end))

    cdef list paddedVariants = padVariants(filteredVariants, refFile, chrom)
    return paddedVariants

###################################################################################################

cdef void callVariantsInRegion(bytes chrom, int start, int end, bamFiles, FastaFile refFile, options, windowGenerator, outputFile, vcf, list samples, dict samplesByID, dict samplesByBAM, Population pop):
    """
    Given a set of BAM files, and a sensibly-sized genomic region (typically ~1MB), call variants in the specified region.
    """
    refFile.setCacheSequence(chrom, start-(10*options.rlen), end+(10*options.rlen))

    cdef long long int contigLength = refFile.refs[chrom].SeqLength
    cdef long long int maxContigPos = contigLength - 1
    cdef bytes refSequenceBytes = refFile.getSequence(chrom,start, min(end+5*options.rlen, refFile.refs[chrom].SeqLength-1))
    cdef char* refSequence = refSequenceBytes
    cdef dict window
    cdef list readBuffers
    cdef bamReadBuffer readBuffer

    try:
        readBuffers = loadBAMData(bamFiles, chrom, start, end, options, samples, samplesByID, samplesByBAM, refSequence)

    except Exception, e:
        logger.error('Exception in region %s:%s-%s. Error was %s' %(chrom, start, end, e))
        logger.warning("Region %s:%s-%s will be skipped" %(chrom, start, end))
        return

    if readBuffers is None:
        logger.info("Skipping region %s:%s-%s as data could not be loaded" %(chrom, start, end))
        return

    cdef list allSortedVariants = generateVariantsInRegion(chrom, start, end, bamFiles, refFile, options, windowGenerator, outputFile, vcf, readBuffers)
    cdef Variant var,nextVar
    cdef list sortedVarList

    for windowIndex,window in enumerate(windowGenerator.WindowsAndVariants(chrom, start, end, maxContigPos, allSortedVariants, options)):

        try:

            chrom = window['chromosome']
            windowStart = window['startPos']
            windowEnd = window['endPos']

            if windowEnd - windowStart > options.maxSize and len(window['variants']) > 0:
                logger.info("Skipping very large window %s:%s-%s of size %s. Max window size is %s (set in option --maxSize)" %(chrom, windowStart, windowEnd, windowEnd - windowStart, options.maxSize))
                continue

            if len(window['variants']) > 0:
                callVariantsInWindow(window, options, refFile, readBuffers, pop, start, end, refSequence)

            if len(window['variants']) > 0 and len(pop.variantPosteriors.keys()) > 0:
                outputCallToVCF(pop.varsByPos, pop.vcfInfo, pop.vcfFilter, pop.haplotypes, pop.genotypes, pop.frequencies, pop.genotypeLikelihoods, pop.goodnessOfFitValues, pop.haplotypeIndexes, pop.readBuffers, pop.nIndividuals, vcf, refFile, outputFile, options, pop.variants, window['startPos'], window['endPos'])

                # Need to make ref calls between variants and at start and end of windows
                if options.outputRefCalls and len(pop.varsByPos.keys()) > 1:
                    lastPos = None
                    lastVars = None

                    for index,(pos,theseVars) in enumerate(pop.varsByPos.iteritems()):
                        if index > 0:
                            lastVarPos = max([v.maxRefPos for v in lastVars])
                            nextVarPos = min([v.minRefPos for v in theseVars]) + 1 # To account for VCF 1-indexing vs Variant 0-indexing

                            if nextVarPos - lastVarPos > 1:
                                for refBlockStart in xrange(lastVarPos+1, nextVarPos, options.refCallBlockSize):
                                    refBlockEnd = min(refBlockStart + options.refCallBlockSize, nextVarPos-1)

                                    if refBlockStart == refBlockEnd:
                                        continue

                                    thisWindow = dict(chromosome=chrom,startPos=refBlockStart, endPos=refBlockEnd,variants=[],nVar=0)
                                    outputRefCall(chrom, pop, vcf, refFile, outputFile, windowIndex, thisWindow, options, readBuffers)
                        lastPos = pos
                        lastVars = theseVars
            else:
                if options.outputRefCalls:
                    outputRefCall(chrom, pop, vcf, refFile, outputFile, windowIndex, window, options, readBuffers)

            if options.compressReads:
                for readBuffer in readBuffers:
                    readBuffer.recompressReadsInCurrentWindow(start, end, refSequence, options.qualBinSize, options.compressReads)

        except Exception, e:
            logger.exception('Exception in window %d-%d. Error was %s' %(window['startPos'],window['endPos'], e))
            logger.warning("Window %s:%s-%s will be skipped" %(chrom, window['startPos'],window['endPos']))

###################################################################################################

cdef int sumMapQualsTimesReadLength(list readBuffers):
    """
    Compute and return the sum of mapping quality phred scores of all good-quality reads
    mapping to this region multuplied by the read lengths.
    """
    cdef int result = 0
    cdef bamReadBuffer theBuffer
    cdef cAlignedRead** rStart = NULL
    cdef cAlignedRead** rEnd = NULL

    for i from 0 <= i < len(readBuffers):
        theBuffer = readBuffers[i]
        rStart = theBuffer.reads.windowStart
        rEnd = theBuffer.reads.windowEnd

        while rStart != rEnd:
            result += (rStart[0].mapq*rStart[0].rlen)
            rStart += 1

    return result

###################################################################################################

def outputRefCall(bytes chrom, Population pop, vcfFile, FastaFile refFile, outputFile, int windowIndex, dict window, options, list readBuffers):
    """
    """
    cdef int windowStart = window['startPos']
    cdef int windowEnd = window['endPos']
    cdef int windowSize = windowEnd - windowStart
    cdef list variants = window['variants']
    cdef int nIndividuals = len(readBuffers)
    cdef bamReadBuffer theBuffer

    cdef int minCov = -1
    cdef int theStart = 0
    cdef int granularity = 1

    for index,theBuffer in enumerate(readBuffers):
        
        for theStart in range(windowStart, windowEnd, granularity):
            if minCov == -1:
                minCov = theBuffer.countReadsCoveringRegion(theStart, theStart+1)
            else:
                minCov = min(minCov, theBuffer.countReadsCoveringRegion(theStart, theStart+1))

    # What should the qual value be for this?
    #
    # 1) If there is no data, then qual == 0.
    #
    # 2) If there are reads but no candidates, then qual is given by a beta-binomial p-value based on the coverage.
    #
    # 3) If there are reads and candidates, then the posterior must be related to the posterior of the best variant
    #    candidate. Try qual == phred( min( (1.0 - max(varPosteriors), sum(mapQuals)/sizeof(window) ) ) )

    cdef int phredPValue = int(-10*log10(betaBinomialCDF(0, minCov, 20, 20)))

    if minCov == 0:
        qual = 0
    else:
        if len(variants) == 0:
            qual = phredPValue
        else:
            maxPost = max( [pop.calculatePosterior(v, 1) for v in variants] )
            maxProbVar = 1.0 - 10**(-0.1*maxPost)
            probRef = 1.0 - maxProbVar
            qual = min(int(round(-10.0 * log10(1.0-probRef))), phredPValue)


    # Sort out VCF output
    ref = refFile.getSequence(chrom, windowStart, windowStart+1)
    alt = None

    # Kludge. To keep within the VCF standard, REF and ALT must be different. At least one of
    # the fields will always be 'N', which should help to avoid confusion with normal variant calls
    if ref == "N":
        alt = ["T"]
    else:
        alt = ["N"]

    id = "."
    linefilter = ["REFCALL"]
    lineinfo = {}
    lineinfo['END'] = [windowEnd]
    lineinfo['Size'] = [windowSize]

    # I'll fill in some of these later.
    lineinfo['FR'] = ["."]
    lineinfo['MMLQ'] = ["."]
    lineinfo['HP'] = ["."]
    lineinfo['TCR'] = ["."]
    lineinfo['WE'] = ["."]
    lineinfo['WS'] = ["."]
    lineinfo['Source'] = ["."]
    lineinfo['FS'] = ["."]
    lineinfo['START'] = ["."]
    lineinfo['PP'] = ["."]
    lineinfo['TR'] = ["."]
    lineinfo['NF'] = ["."]
    lineinfo['TCF'] = ["."]
    lineinfo['NR'] = ["."]
    lineinfo['TC'] = ["."]
    lineinfo['MGOF'] = ["."]
    lineinfo['SbPval'] = ["."]
    lineinfo['ReadPosRankSum'] = ["."]
    lineinfo['MQ'] = ["."]
    lineinfo['QD'] = ["."]
    lineinfo['SC'] = ["."]
    lineinfo['BRF'] = ["."]
    lineinfo['HapScore'] = ["."]

    # Most of these are empty
    theFormat = ['GT:GL:GOF:GQ:NR:NV']

    vcfDataLine = {'chrom':chrom,'pos':windowStart,'ref':ref,'alt':alt,'id':id,'info':lineinfo,'filter':linefilter,'qual':qual,'format':theFormat}

    for i from 0 <= i < nIndividuals:

        theBuffer = readBuffers[i]
        thisSample = theBuffer.sample

        # Missing genotype call. Probably due to zero coverage for this individual.
        if theBuffer.reads.windowEnd - theBuffer.reads.windowStart == 0:
            vcfDataLine[thisSample] = dict(GT=[[".", "/", "."]], GL=[-1,-1,-1], GQ=[-1],GOF=[-1],NR=[0], NV=[0])
        else:
            vcfDataLine[thisSample] = dict(GT=[[".", "/", "."]], GL=[-1,-1,-1], GQ=[-1],GOF=[-1],NR=[theBuffer.reads.windowEnd - theBuffer.reads.windowStart], NV=[0])

    vcfFile.write_data(outputFile, vcfDataLine)

###################################################################################################

class PlatypusSingleProcess(object):
    """
    Simple class to represent a single Platypus process.
    """
    samples = None
    samplesByID = None
    samplesByBAM = None
    bamFiles = None

    def __init__(self, fileName, options, regions, continuing=False, samples=None, bamFileNames=None, samplesByID=None, samplesByBAM=None, bamFiles=None, theLocks=None):
        """
        Constructor. Store options and output file name
        """
        self.options = options
        self.fileName = fileName
        self.regions = regions
        self.continuing = continuing

        cdef Samfile theBamFile

        if bamFiles is not None and theLocks is not None:
            for theBamFile,theLock in zip(bamFiles,theLocks):
                theBamFile.lock = theLock

        if self.options.verbosity >= 3:
            logger.info("Searching for variants in the following regions: %s" %(self.regions))

        # Need to load data etc.
        self.bamFileNames = options.bamFiles

        if len(self.bamFileNames) == 1 and not self.bamFileNames[0].lower().endswith(".bam"):
            logger.debug("Treating --bamFiles argument, %s as a text file containing a list of BAM file names" %(self.bamFileNames))
            self.bamFileNames = platypusutils.getBAMFileNamesFromTextFile(self.bamFileNames[0])

        self.samples,self.samplesByID,self.samplesByBAM,self.bamFiles = platypusutils.getSampleNamesAndLoadIterators(self.bamFileNames, self.regions, options)

        if options.refFile.endswith(".gz") or options.refFile.endswith(".bz2") or options.refFile.endswith(".bgz"):
            logger.error("Reference file-name (%s) looks like a compressed file-name. Please un-compress the reference FASTA file before running Platypus" %(options.refFile))
            raise StandardError, "Invalid reference FASTA file supplied"

        self.refFile = fastafile.FastaFile(options.refFile, options.refFile + ".fai")

        # Set these to sensible values, based on population size
        if self.options.maxHaplotypes == -1:
            self.options.maxHaplotypes = 257

        # maxHaplotypes has to be capped
        self.options.originalMaxHaplotypes = self.options.maxHaplotypes
        self.options.maxHaplotypes = min(257, self.options.maxHaplotypes)
        self.options.maxGenotypes = min(33153, platypusutils.nCombinationsWithReplacement(self.options.maxHaplotypes, 2))
        #self.options.maxGenotypes = 33153
        #self.options.maxHaplotypes = min(257, self.options.maxHaplotypes)

        logger.debug("Max haplotypes used for initial haplotype filtering = %s" %(self.options.originalMaxHaplotypes))
        logger.debug("Max haplotypes used for genotype generation = %s" %(self.options.maxHaplotypes))
        logger.debug("Max genotypes = %s" %(self.options.maxGenotypes))

        options.nInd = len(set(self.samples))
        
        if options.alignScoreFile != "":
            logger.info("Alignment scores of reads with haplotypes are written to %s" %(options.alignScoreFile))
            fo = open(options.alignScoreFile, "w")
            fo.write("#Alignment scores of reads against haplotypes within a each window for each sample\n")
            fo.close()



    def run(self):
        """
        Run a single Platypus process
        """
        # Cache reference sequence for this region. This should result in a substantial speed-up for
        # the fastafile.getSequence function.
        self.vcf = vcf.VCF()
        self.vcf.setheader(vcfHeader + [('platypusOptions', str(self.options))])
        self.vcf.setsamples(sorted(set(self.samples)))
        self.vcf.setinfo(vcfutils.vcfInfoSignature)
        self.vcf.setfilter(vcfutils.vcfFilterSignature)
        self.vcf.setformat(vcfutils.vcfFormatSignature)
        self.windowGenerator = window.WindowGenerator()

        if self.continuing:
            logger.info("Opening file %s for appending" %(self.fileName))
            self.outputFile = open(self.fileName, 'a')
        else:
            if self.fileName == "-":
                self.outputFile = sys.stdout
            else:
                self.outputFile = open(self.fileName, 'w')
            self.vcf.writeheader(self.outputFile)

        cdef Population pop = Population(self.options)

        for index,(chrom,start,end) in enumerate(self.regions):

            if index % 10 == 0:
                logger.info("Processing region %s:%s-%s. (Only printing this message every 10 regions of size %s)" %(chrom, start, end, self.options.bufferSize))
            elif self.options.verbosity >= 3:
                logger.info("Processing region %s:%s-%s" %(chrom, start, end))
            else:
                pass

            callVariantsInRegion(chrom, start, end, self.bamFiles, self.refFile, self.options, self.windowGenerator, self.outputFile, self.vcf, self.samples, self.samplesByID, self.samplesByBAM, pop)

        if self.fileName != "-":
            self.outputFile.close()

###################################################################################################

class PlatypusMultiProcess(multiprocessing.Process):
    """
    Simple class to represent a single Platypus process, which is run as part of a
    multi-process job.
    """
    def __init__(self, fileName, options, regions, continuing=False, samples=None, bamFileNames=None, samplesByID=None, samplesByBAM=None, bamFiles=None, theLocks=None):
        """
        Constructor. Store options and output file name
        """
        multiprocessing.Process.__init__(self)
        self.singleProcess = PlatypusSingleProcess(fileName, options, regions, continuing, samples, bamFileNames, samplesByID, samplesByBAM, bamFiles, theLocks)

    def run(self):
        """
        Run a single Platypus process
        """
        self.singleProcess.run()

##################################################################################################
