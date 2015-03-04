"""
Various classes and functions for use in generating and processing
variant candidates.
"""
from __future__ import division

import logging
import fastafile
import cython
import math

cimport chaplotype
cimport variant
cimport fastafile
cimport cpopulation

from variant cimport Variant,PLATYPUS_VAR,ASSEMBLER_VAR,FILE_VAR
from chaplotype cimport Haplotype,computeOverlapOfReadAndHaplotype
from fastafile cimport FastaFile
from operator import attrgetter
from cpopulation cimport Population
from cgenotype cimport generateAllGenotypesFromHaplotypeList
from cgenotype cimport HaploidGenotype,DiploidGenotype
from htslibWrapper cimport cAlignedRead
from htslibWrapper cimport cAlignedRead
from htslibWrapper cimport Read_IsReverse
from htslibWrapper cimport Read_IsPaired
from htslibWrapper cimport Read_IsProperPair
from htslibWrapper cimport Read_IsDuplicate
from htslibWrapper cimport Read_IsUnmapped
from htslibWrapper cimport Read_MateIsUnmapped
from htslibWrapper cimport Read_MateIsReverse
from htslibWrapper cimport Read_IsQCFail
from htslibWrapper cimport Read_IsReadOne
from htslibWrapper cimport Read_IsSecondaryAlignment
from htslibWrapper cimport Read_SetQCFail
from cwindow cimport bamReadBuffer
from platypusutils cimport leftNormaliseIndel
from platypusutils cimport isHaplotypeValid
from itertools import combinations
from heapq import heappush,heappop,heappushpop

nSupportingReadsGetter = attrgetter("nSupportingReads")
logger = logging.getLogger("Log")

###################################################################################################

cdef extern from "math.h":
    double exp(double)
    double sqrt(double)
    double log(double)
    double log2(double)
    double log10(double)
    double fabs(double)

###################################################################################################

cdef list padVariants(list sortedVariants, FastaFile refFile, bytes chrom):
    """
    All variants which overlap will be padded so that they start at the same place. This
    is just for later convenience when calling.
    """
    cdef Variant thisVar = None
    cdef Variant nextVar = None
    cdef bytes padding   = None
    cdef list paddedVars = []
    
    for nextVar in sortedVariants:
        if thisVar is None:
            thisVar = nextVar
            paddedVars.append(thisVar)
        else:
            # Overlapping variants. These need to be padded
            if thisVar.maxRefPos >= nextVar.minRefPos and thisVar.refPos < nextVar.refPos:
                padding           = refFile.getSequence(chrom, thisVar.minRefPos+1, nextVar.minRefPos+1)
                nextVar.minRefPos = thisVar.minRefPos
                nextVar.refPos    = thisVar.refPos
                nextVar.removed   = padding + nextVar.removed
                nextVar.added     = padding + nextVar.added
                nextVar.nAdded    = len(nextVar.added)
                nextVar.nRemoved  = len(nextVar.removed)
                nextVar.hashValue = -1
                paddedVars.append(nextVar)
                #logger.info("Done padding. Now vars are %s and %s" %(nextVar, thisVar))
            # No overlap
            else:
                paddedVars.append(nextVar)
            
            # Always keep 'thisVar' as the variant which extends furthest. Variants are sorted by
            # minRefPos, not maxRefPos.
            if nextVar.maxRefPos > thisVar.maxRefPos:
                thisVar = nextVar
    
    return paddedVars

###################################################################################################

cdef list filterVariants(list varList, FastaFile refFile, int maxReadLength, int minSupport, int maxDiff, int verbosity, options):
    """
    Generator function, this calls the Candidates generator function, and gets a list of
    sorted variant candidates. This list is then merged such that only unique variant candidates
    are returned, and supporting reads are accumulated in the returned candidates.

    No additional filtering is done at this step.
    """
    cdef Variant v = None
    cdef Variant lastVariant = None
    cdef list filteredVariants = []
    cdef int minReads = options.minReads
    cdef int maxSize = options.maxSize

    for v in varList:

        # Start.
        if lastVariant is None:
            lastVariant = v

        # Same as last variant: add supporting reads.
        elif v == lastVariant:
            lastVariant.addVariant(v)

        # Not the same: add to filtered list, if it passes filters.
        else:
            support = lastVariant.nSupportingReads
            source = lastVariant.varSource
            varSize = max(lastVariant.nAdded, lastVariant.nRemoved)

            # Skip Vars with < minSupport supporting reads, if they only come from Platypus, and they are small (< 20 bases)
            if support < minSupport and varSize < 15 and source & PLATYPUS_VAR and not (source & ASSEMBLER_VAR) and not(source & FILE_VAR):
                pass

            # Skip Vars with < minReads supporting reads, if they only come from Platypus, and they are large (< 20 bases)
            elif support < minReads and varSize >= 15 and source & PLATYPUS_VAR and not (source & ASSEMBLER_VAR) and not(source & FILE_VAR):
                pass

            # Can't currently deal with very large variants due to hashing issues.
            elif varSize > maxSize:
                pass

            else:

                #nSupportingSamples = len(lastVariant.supportingSamples)

                #if support/nSupportingSamples < minSupport:
                #    if verbosity >= 3:
                #        logger.debug("Variant %s failed only multi-sample filter. nReads = %s. nSamples = %s. Ratio = %s." %(lastVariant, support, nSupportingSamples, support/nSupportingSamples))
                #    pass

                #else:
                if verbosity >= 3:
                    logger.debug("Adding variant %s to filtered list" %(lastVariant))

                filteredVariants.append(lastVariant)

            lastVariant = v

    if lastVariant is not None:

        support = lastVariant.nSupportingReads
        source = lastVariant.varSource

        # Skip Vars with < minSupport supporting reads, if they only come from Platypus
        if support < minSupport and source & PLATYPUS_VAR and not (source & ASSEMBLER_VAR) and not(source & FILE_VAR):
            pass
        else:
            if verbosity >= 3:
                logger.debug("Adding variant %s to filtered list" %(lastVariant))

            filteredVariants.append(lastVariant)

    return sorted(filteredVariants)

###################################################################################################

cdef void filterVariantsInWindow(dict thisWindow, bytes chrom, int windowStart, int windowEnd, FastaFile refFile, options, list variants, Haplotype refHaplotype, list readBuffers):
    """
    In some windows there are too many variants to crunch through the complete Platypus model, so we need to
    filter the variant candidates in those windows, to get rid of noise.

    This filtering is currently a 2-step process:

    1) For pathological cases (>50 variants in a window), first sort the variants by coverage and take the top 50
       variants. This will not work for assembler variants, or candidates from VCF input files, or difficult indels.

    2) Filter using the population model, but only on a sub-set of haplotypes (the set of haplotypes containing one
       variant each). Take the variants with the highest posterior values.
    """
    cdef int nVar = len(variants)

    if nVar > 50:
        logger.debug("Taking top 50 variant by coverage")
        variants = sorted(sorted(variants, key=nSupportingReadsGetter, reverse=True)[0:50])
        nVar = len(variants)

    logger.debug("There are %s variants. Filtering variants with EM method" %(nVar))
    cdef list varHaps = [refHaplotype] + [Haplotype(chrom, windowStart, windowEnd, (v,), refFile, options.rlen, options) for v in variants]
    cdef list varGens = generateAllGenotypesFromHaplotypeList(varHaps)
    cdef Population pop = Population(options)
    pop.setup(variants, varHaps, varGens, options.nInd, options.verbosity, readBuffers)
    pop.call(options.maxEMIterations, 0)
    cdef list varsByPost = []

    for v in variants:
        post = pop.calculatePosterior(v)
        varsByPost.append( (post,v) )

    varsByPost.sort(reverse=True)
    variants = sorted([theTuple[1] for theTuple in varsByPost[:options.maxVariants]])
    thisWindow['variants'] = variants

###################################################################################################

cdef double computeBestScoreForHaplotype(list readBuffers, Haplotype hap):
    """
    """
    cdef bamReadBuffer readBuff
    cdef cAlignedRead** readBegin = NULL
    cdef cAlignedRead** readEnd = NULL
    cdef double bestScoreThisHap = -1e20
    cdef double scoreThisHapAndSample = 0.0

    for readBuff in readBuffers:
        scoreThisHapAndSample = 0.0
        readBegin = readBuff.reads.windowStart
        readEnd = readBuff.reads.windowEnd

        while readBegin != readEnd:
            scoreThisHapAndSample += hap.alignSingleRead(readBegin[0], False)
            readBegin += 1

        bestScoreThisHap = max(bestScoreThisHap, scoreThisHapAndSample)

    return bestScoreThisHap

###################################################################################################

cdef double computeBestScoreForGenotype(list readBuffers, DiploidGenotype gt, int windowSize, int targetCoverage):
    """
    """
    cdef bamReadBuffer readBuff
    cdef cAlignedRead** readBegin = NULL
    cdef cAlignedRead** readEnd = NULL
    cdef cAlignedRead** badReadsBegin = NULL
    cdef cAlignedRead** badReadsEnd = NULL
    cdef cAlignedRead** brokenReadsBegin = NULL
    cdef cAlignedRead** brokenReadsEnd = NULL
    cdef double bestScoreThisHap = -1e20
    cdef double scoreThisHapAndSample = 0.0
    cdef double scoreThisHapAndSample2 = 0.0
    cdef int nIndividuals = len(readBuffers)
    cdef int meanCoverage = 0
    cdef int sampleRate = 0

    assert targetCoverage > 0

    for individualIndex,readBuff in enumerate(readBuffers):

        scoreThisHapAndSample = 0.0
        readBegin = readBuff.reads.windowStart
        readEnd = readBuff.reads.windowEnd
        badReadsBegin = readBuff.badReads.windowStart
        badReadsEnd = readBuff.badReads.windowEnd
        brokenReadsBegin = readBuff.brokenMates.windowStart
        brokenReadsEnd = readBuff.brokenMates.windowEnd

        if readBegin == readEnd:
            continue

        meanCoverage = readBegin[0].rlen * (readEnd - readBegin) // windowSize
        sampleRate = max(1, meanCoverage // targetCoverage)
        #logger.info("nReads = %s. Mean coverage for sample %s in window is %s. Window size = %s. Target cov = %s. Sampling one read per %s" %(readEnd - readBegin, individualIndex, meanCoverage, windowSize, targetCoverage, sampleRate))

        while readBegin < readEnd:
            score1 = gt.hap1.alignSingleRead(readBegin[0], False)
            score2 = gt.hap2.alignSingleRead(readBegin[0], False)
            scoreThisHapAndSample += log(0.5*(exp(score1) + exp(score2)))
            #readBegin += 1
            readBegin += sampleRate

        #scoreThisHapAndSample = gt.calculateDataLikelihood(readBegin, readEnd, badReadsBegin, badReadsEnd, brokenReadsBegin, brokenReadsEnd, individualIndex, nIndividuals, NULL)
        bestScoreThisHap = max(bestScoreThisHap, scoreThisHapAndSample)

    return bestScoreThisHap

###################################################################################################

cdef int doesAtLeastOneReadSupportHaplotype(list readBuffers, Haplotype hap):
    """
    Return True if at least one read out of the lot gives an alignment score of 0.
    """
    cdef bamReadBuffer readBuff
    cdef cAlignedRead** readBegin = NULL
    cdef cAlignedRead** readEnd = NULL
    cdef double scoreThisHapAndSample = 0.0
    cdef int found = 0
    cdef int hapVarLen = hap.maxVarPos - hap.minVarPos
    cdef int readLen = 0

    for readBuff in readBuffers:
        scoreThisHapAndSample = 0.0
        readBegin = readBuff.reads.windowStart
        readEnd = readBuff.reads.windowEnd

        while readBegin != readEnd:

            readLen = readBegin[0].end - readBegin[0].pos
            readOverlap = computeOverlapOfReadAndHaplotype(hap.minVarPos, hap.maxVarPos, readBegin[0])

            if Read_IsQCFail(readBegin[0]) or readOverlap < min(readLen//2, hapVarLen):
                pass
            else:
                scoreThisHapAndSample = hap.alignSingleRead(readBegin[0], False)

                if scoreThisHapAndSample > -1e-2:
                    return True

            readBegin += 1

    return False

###################################################################################################

cdef int countReadsSupportingHaplotype(list readBuffers, Haplotype hap):
    """
    Return the number of reads supporting a given haplotype
    """
    cdef bamReadBuffer readBuff
    cdef cAlignedRead** readBegin = NULL
    cdef cAlignedRead** readEnd = NULL
    cdef double scoreThisHapAndSample = 0.0
    cdef int nSupporting = 0
    cdef int hapVarLen = hap.maxVarPos - hap.minVarPos
    cdef int readLen = 0

    for readBuff in readBuffers:
        scoreThisHapAndSample = 0.0
        readBegin = readBuff.reads.windowStart
        readEnd = readBuff.reads.windowEnd

        while readBegin != readEnd:

            readLen = readBegin[0].end - readBegin[0].pos
            readOverlap = computeOverlapOfReadAndHaplotype(hap.minVarPos, hap.maxVarPos, readBegin[0])

            if Read_IsQCFail(readBegin[0]) or readOverlap < min(readLen//2, hapVarLen):
                pass
            else:
                scoreThisHapAndSample = hap.alignSingleRead(readBegin[0], False)

                if scoreThisHapAndSample > -1e-2:
                    nSupporting += 1

            readBegin += 1

    return nSupporting

###################################################################################################

cdef double computeVariantReadSupportFrac(Variant variant, bamReadBuffer readBuffer):
    """
    Compute and return the fraction of reads which support the specified variant. For now,
    I am simply using the number of supporting reads specified in the variant object over the
    """
    cdef int varReadPos = variant.refPos
    cdef int nVarReads = variant.nSupportingReads
    cdef int nTotalReads = readBuffer.countReadsCoveringRegion(varReadPos, varReadPos+1)

    if nTotalReads == 0:
        return 0.0
        #logger.error("Variant %s supported by %s reads, but nTotalReads = %s" %(variant, nVarReads, nTotalReads))
    cdef double varFrac = <double>(nVarReads)/nTotalReads
    #logger.info("Var read frac for var %s is %s (%s/%s)" %(variant, varFrac, nVarReads, nTotalReads))
    return varFrac

###################################################################################################

cdef list getFilteredHaplotypes(dict thisWindow, bytes chrom, int windowStart, int windowEnd, FastaFile refFile, options, list variants, Haplotype refHaplotype, list readBuffers):
    """
    Here we attempt to filter variants in windows with many (>maxVariants) candidates. This is done by constructing
    all possible haplotypes, and computing likelihoods for each haplotype. Haplotypes (and therefore their variants)
    are kept if they
    """
    cdef Haplotype hap
    cdef Variant tempVar
    cdef tuple varThisHap
    cdef tuple varsThisHap2
    cdef tuple varsFromBothSets
    cdef list allHaps = []
    cdef list hapsByBestScore = []
    cdef list tempOldHaps = []
    cdef list varsSortedByCoverage = []
    cdef double bestScoreThisHap = 0.0
    cdef double scoreThisHapAndSample = 0.0
    cdef int originalMaxHaplotypes = options.originalMaxHaplotypes - 1 # Ref will be added later
    cdef int maxHaplotypes = options.maxHaplotypes - 1 # Ref will be added later
    cdef int nVarsInHap = 0
    cdef int maxReadLength = options.rlen
    cdef int nVarsProcessed = 0
    cdef int nVars = len(variants)
    cdef int varChunkSize = 1 # How many vars to process in one go?
    cdef int nTempOldHaps = 0
    cdef int nHapsTried = 0
    cdef int nHapsDone = 0
    cdef int nValidHaps = 0
    cdef int verbosity = options.verbosity
    cdef int windowSize = windowEnd - windowStart
    cdef int targetCoverage = options.coverageSamplingLevel
    cdef DiploidGenotype gt = DiploidGenotype(refHaplotype, refHaplotype)

    # If nVar is small, or we've already checked, and the number of valid haplotypes is less than the max, return all combinations
    if nVars <= log2(maxHaplotypes) or (options.filterVarsByCoverage and options.maxVariants <= log2(maxHaplotypes)):

        for nVarsInHap from 1 <= nVarsInHap <= nVars:
            for varsThisHap in combinations(variants, nVarsInHap):
                if isHaplotypeValid(varsThisHap):
                    nValidHaps += 1
                    hap = Haplotype(chrom, windowStart, windowEnd, varsThisHap, refFile, maxReadLength, options)
                    allHaps.append(hap)

                    ## Check for complex haps
                    #if nVars > 3:

                    #    # Accept
                    #    if doesAtLeastOneReadSupportHaplotype(readBuffers, hap):
                    #        allHaps.append(hap)

                    #    # Reject
                    #    else:
                    #        pass

                    ## Accept all haps when nVars is small
                    #else:
                    #    allHaps.append(hap)

        #if nVars > 3:
        #    logger.info("%s vars, %s valid haps. %s accepted haps" %(nVars, nValidHaps, len(allHaps)))

        return allHaps

    # Otherwise, we need to do some filtering. First split the variant list into 2 lists, with all the well-supported
    # variants in ones list, and the poorly-supported (small nSupportingReads) variants in the other.
    else:
        varsSortedByCoverage = sorted(variants, key=nSupportingReadsGetter, reverse=True)
        #logger.debug("Vars by coverage are ...")
        #for v in varsSortedByCoverage:
        #    logger.debug(v)

        #logger.info("Taking top %s haps in window %s:%s-%s out of (%s vars, %s haps)." %(maxHaplotypes, chrom, windowStart, windowEnd, nVars, math.pow(2,nVars) -1))

        while nVarsProcessed < nVars:

            tempVar = varsSortedByCoverage[nVarsProcessed]
            tempOldHaps = sorted(hapsByBestScore)
            nTempOldHaps = len(tempOldHaps)
            varThisHap = (tempVar,)

            #logger.info("temp var = %s. n old haps = %s. nHaps tried = %s. nHaps done = %s." %(tempVar, nTempOldHaps, nHapsTried, nHapsDone))

            gt.hap2 = Haplotype(chrom, windowStart, windowEnd, varThisHap, refFile, maxReadLength, options)
            bestScoreThisHap = computeBestScoreForGenotype(readBuffers, gt, windowSize, targetCoverage)
            nHapsDone += 1

            if len(hapsByBestScore) < originalMaxHaplotypes:
                heappush(hapsByBestScore, (bestScoreThisHap, varThisHap))
            else:
                heappushpop(hapsByBestScore, (bestScoreThisHap, varThisHap))

            for score,varsThisHap2 in tempOldHaps:

                # Create a set of variants which is the combination of this set with one of the top haplotypes
                varsFromBothSets = tuple(sorted(varThisHap + varsThisHap2))
                nHapsTried += 1

                if isHaplotypeValid(varsFromBothSets):

                    gt.hap2 = Haplotype(chrom, windowStart, windowEnd, varsFromBothSets, refFile, maxReadLength, options)
                    bestScoreThisHap = computeBestScoreForGenotype(readBuffers, gt, windowSize, targetCoverage)

                    if len(hapsByBestScore) < originalMaxHaplotypes:
                        heappush(hapsByBestScore, (bestScoreThisHap, varsFromBothSets))
                    else:
                        heappushpop(hapsByBestScore, (bestScoreThisHap, varsFromBothSets))

                    nHapsDone += 1

            nVarsProcessed += 1

    #logger.debug("Summarising top 10 haplotypes...")
    #for index,(score,varsThisHap) in enumerate(sorted(hapsByBestScore, reverse=True)):
    #    if index < 10:
    #        logger.debug("Score = %s. Vars = %s" %(score, varsThisHap))
    #    else:
    #        break
    #logger.debug("Done.")

    #for score,varsThisHap in hapsByBestScore:

    for index,(score,varsThisHap) in enumerate(sorted(hapsByBestScore, reverse=True)):
        if index < maxHaplotypes:
            hap = Haplotype(chrom, windowStart, windowEnd, varsThisHap, refFile, maxReadLength, options)
            allHaps.append(hap)
        else:
            break

    #logger.info("Done.")
    return allHaps

###################################################################################################

cdef list calcNPossibleHaplotypes(list variants):
    """
    Check all possible haplotypes for consistency, and return
    the number of consistent haplotypes.
    """
    cdef int nVars = len(variants)
    cdef list haps = []

    for nVarsInHap from 1 <= nVarsInHap <= nVars:
        for varsThisHap in combinations(variants, nVarsInHap):
            if isHaplotypeValid(varsThisHap):
                haps.append(varsThisHap)

    return haps

###################################################################################################

cdef int recurse(list variants, dict nonOverlaps, int maxHaps):
    """
    """
    cdef int haps = 0
    cdef int newHaps = 0
    cdef Variant variant

    for variant in variants:
        haps += 1
        newHaps = recurse(nonOverlaps[variant], nonOverlaps, maxHaps)

        if newHaps == -1 or (haps + newHaps > maxHaps):
            return -1
        else:
            haps += newHaps

    return haps

###################################################################################################

cdef int countHapsByNVarsThatDontOverlap(list variants, int maxHaps):
    """
    Loop through variants, and count those that don't overlap with others, according to the
    usual rules (variants are assumed to be sorted). Return the count of non-overlapping variants.
    These are the ones which can be considered together in haplotypes.

    This is quadratic in nVars
    """
    nonOverlaps = {}
    cdef Variant var_i
    cdef Variant var_j

    for i in range(len(variants)):
        var_i = variants[i]
        nonOverlaps[var_i] = []
        for j in range(i+1, len(variants)):
            var_j = variants[j]
            if not var_i.overlaps(var_j):
                nonOverlaps[var_i].append(var_j)

    return recurse(nonOverlaps.keys(), nonOverlaps, maxHaps)

###################################################################################################

cdef void filterVariantsByCoverage(dict thisWindow, bytes chrom, int windowStart, int windowEnd, FastaFile refFile, options, list variants, Haplotype refHaplotype, list readBuffers):
    """
    Simple filtering: pick top maxVariants vars by coverage, prioritising indels.
    """
    cdef int maxVar = options.maxVariants
    cdef Variant thisVar
    cdef list filteredVars = []
    cdef list temp = []

    cdef int nVars = len(variants)
    cdef int maxSupport = max( [thisVar.nSupportingReads for thisVar in variants] )
    cdef int nPossHaps = 0
    cdef int nPossHapsFromVars = 0

    logger.debug("Taking top %s variants (out of %s) by coverage in window %s:%s-%s" %(maxVar, len(variants), chrom, windowStart, windowEnd))

    #cdef int nHaps = countHapsByNVarsThatDontOverlap(variants, options.maxHaplotypes)
    #logger.debug("There are %s valid haps in this window" %(nHaps))

    # This actually means that the total number of valid haplotypes is <= maxHaplotypes
    #if nHaps > 0:
    #    return

    cdef Haplotype oneVarHap
    cdef double score = 0.0

    for thisVar in variants:
        #oneVarHap = Haplotype(thisWindow['chromosome'], thisWindow['startPos'], thisWindow['endPos'], (thisVar,), refFile, options.rlen, options)
        #score = computeBestScoreForGenotype(readBuffers, refHaplotype, oneVarHap)
        #temp.append( (score, thisVar) )

        ##logger.info("Score for var %s is %s" %(thisVar, score))

        # Prioritise assembler variants
        if thisVar.varSource == ASSEMBLER_VAR:
            temp.append( (maxSupport+1, thisVar) )
        else:
            temp.append( (thisVar.nSupportingReads, thisVar) )

    # Sort by coverage, descending
    temp.sort(reverse=True)

    # Take top candidates
    filteredVars = sorted( [x[1] for x in temp[0:maxVar]] )

    if options.verbosity >= 3:
        logger.debug("")
        logger.debug("Window variants before filtering = %s" %(variants))
        logger.debug("Window variants after filtering = %s" %(filteredVars))
        logger.debug("")

    thisWindow['variants'] = filteredVars

###################################################################################################

cdef list getHaplotypesInWindow(dict window, int nReads, FastaFile refFile, int maxCoverage, int minMapQual, int minBaseQual, int maxHaplotypes, int maxVariants, int maxReadLength, int verbosity, list readBuffers, options):
    """
    Takes a window with some potential variants. Using all data in this region, calculate what the best possible
    haplotyes are and return them. If there are zero reads covering this window, then simply return the reference
    haplotype; if the mean coverage is > maxCoverage, then raise an exception.

    Returns:

    haplotypes: a list of haplotypes
    """
    cdef bytes windowChr = window['chromosome']
    cdef int windowStart = window['startPos']
    cdef int windowEnd = window['endPos']
    cdef list variants = window['variants']
    cdef int nVar = len(variants)
    cdef int readLength = 0
    cdef Haplotype refHaplotype = Haplotype(windowChr, windowStart, windowEnd, (), refFile, maxReadLength, options)

    # Make sure that we log regions of zero coverage. This should not really happen.
    if nReads == 0:
        if not options.outputRefCalls:
            logger.warning("Coverage is zero in window %s:%s-%s. Variants considered in this window are %s" %(windowChr, windowStart, windowEnd, variants))
        return [refHaplotype]

    return getFilteredHaplotypes(window, windowChr, windowStart, windowEnd, refFile, options, variants, refHaplotype, readBuffers)

###################################################################################################

cdef list getAllHLAHaplotypesInRegion(bytes chrom, int windowStart, int windowEnd, FastaFile refFile, options, list variants, Haplotype refHaplotype, list readBuffers):

    """
    get all haplotypes in a window
    """
    cdef Haplotype hap
    cdef Variant tempVar, thisVar
    cdef tuple varThisHap
    cdef tuple varsThisHap2
    cdef tuple varsFromBothSets
    cdef list allHaps = []
    cdef list hapsByBestScore = []
    cdef list tempOldHaps = []
    cdef double bestScoreThisHap = 0.0
    cdef double scoreThisHapAndSample = 0.0
    cdef int originalMaxHaplotypes = options.originalMaxHaplotypes - 1 # Ref will be added later
    cdef int maxHaplotypes = options.maxHaplotypes - 1 # Ref will be added later
    cdef int nVarsInHap = 0
    cdef int maxReadLength = options.rlen
    cdef int nVars = len(variants)
    cdef int varChunkSize = 1 # How many vars to process in one go?
    cdef int nTempOldHaps = 0
    cdef int nValidHaps = 0
    cdef int verbosity = options.verbosity
    cdef int windowSize = windowEnd - windowStart
    cdef int targetCoverage = options.coverageSamplingLevel
    cdef DiploidGenotype gt = DiploidGenotype(refHaplotype, refHaplotype)

    cdef list assemblerVars = []
    cdef list assemblerHaps = []
    cdef list outputHaps = []
    cdef int nAssemblerVars = 0   

    cdef int nHaps = 0
    cdef Haplotype thisHap
    cdef Haplotype bestHap
    cdef int nHapsProcessed = 0

    for tempVar in variants:
        if tempVar.varSource == FILE_VAR :
             hap = Haplotype(chrom, windowStart, windowEnd, (tempVar, ), refFile, maxReadLength, options)
             allHaps.append(hap)
    nHaps = len(allHaps)

    # If nVar is small, or we've already checked, and the number of valid haplotypes is less than the max, return all combinations
    # Otherwise, we need to do some filtering. 
    maxHaplotypes = 150
    if nHaps <= maxHaplotypes :
        return allHaps
    else:
        while nHapsProcessed < nHaps:
            thisHap = allHaps[nHapsProcessed]
            tempOldHaps = sorted(hapsByBestScore)
            nTempOldHaps = len(tempOldHaps)
            
            gt.hap2 = thisHap
            bestScoreThisHap = computeBestScoreForHaplotype(readBuffers, thisHap)#computeBestScoreForGenotype(readBuffers, gt, windowSize, targetCoverage)
            if len(hapsByBestScore) < originalMaxHaplotypes:
                heappush(hapsByBestScore, (bestScoreThisHap, thisHap))
            else:
                heappushpop(hapsByBestScore, (bestScoreThisHap, thisHap))
            nHapsProcessed += 1

        for index,(score, thisHap) in enumerate(sorted(hapsByBestScore, reverse=True)):
            if index < maxHaplotypes/2:
                outputHaps.append(thisHap)
            else:
                break

        (score, bestHap) = sorted(hapsByBestScore, reverse=True)[0]
        gt.hap1 = bestHap
        nHapsProcessed = 0

        while nHapsProcessed < nHaps:
            thisHap = allHaps[nHapsProcessed]
            gt.hap2 = thisHap
            bestScoreThisHap = computeBestScoreForGenotype(readBuffers, gt, windowSize, targetCoverage)
            
            if len(hapsByBestScore) < originalMaxHaplotypes:
                heappush(hapsByBestScore, (bestScoreThisHap, thisHap))
            else:
                heappushpop(hapsByBestScore, (bestScoreThisHap, thisHap))
            nHapsProcessed += 1
  
    for index,(score, thisHap) in enumerate(sorted(hapsByBestScore, reverse=True)):
        if index < maxHaplotypes/2:
            outputHaps.append(thisHap)
        else:
            break

    return outputHaps

###################################################################################################
cdef  Variant normaliseVar(Variant thisVar):
        """
        Trimming leading and trailing bases of the variants
        """
         
        if thisVar.nRemoved ==1:
            return thisVar
        #Trimming leading bases
        cdef char* added = thisVar.added
        cdef char* removed = thisVar.removed
        cdef int refPos = thisVar.refPos
        while  len(added) >0 and len(removed)>0 and removed[0] == added[0]:
            added +=1
            removed +=1
            refPos +=1
           
        #Trimming trailing bases
        while  len(added) >0 and len(removed)>0 and removed[len(removed)-1] == added[len(added)-1]:
            added[len(added)-1] = '\0'
            removed[len(removed)-1] = '\0'

        return Variant(thisVar.refName, refPos, removed, added, thisVar.nSupportingReads, thisVar.varSource)
###################################################################################################
cdef  Variant trimLongVar(Variant thisVar, int windowStart, int windowEnd):
        """
        Trimming leading and trailing bases of the variants
        """
         
        if thisVar.nRemoved ==1:
            return thisVar
        #Trimming leading bases
        cdef char* added = thisVar.added
        cdef char* removed = thisVar.removed
        cdef int refPos = thisVar.refPos
        cdef int diff = 0
        #trim at the beginning
        if len(added) == len(removed):
            if  refPos + len(removed) > windowEnd:
                diff =  refPos + len(removed) - windowEnd
                added[len(added) - diff] = '\0'
                removed[len(removed) - diff] = '\0'
            if refPos < windowStart:
                diff = windowStart - refPos
                added += diff
                removed += diff  
        while  len(added) >0 and len(removed)>0 and removed[0] == added[0]:
            added +=1
            removed +=1
            refPos +=1
           
        #Trimming trailing bases
        while  len(added) >0 and len(removed)>0 and removed[len(removed)-1] == added[len(added)-1]:
            added[len(added)-1] = '\0'
            removed[len(removed)-1] = '\0'

        return Variant(thisVar.refName, refPos, removed, added, thisVar.nSupportingReads, thisVar.varSource)
###################################################################################################
cdef list getAllAssemblerHaplotypesInRegion(bytes chrom, int windowStart, int windowEnd, FastaFile refFile, options, list variants, Haplotype refHaplotype, list readBuffers):

    """
    get all haplotypes in a window
    """
    cdef Haplotype hap
    cdef Variant tempVar, thisVar
    cdef tuple varThisHap
    cdef tuple varsThisHap2
    cdef tuple varsFromBothSets
    cdef list allHaps = []
    cdef list hapsByBestScore = []
    cdef list tempOldHaps = []
    cdef double bestScoreThisHap = 0.0
    cdef double scoreThisHapAndSample = 0.0
    cdef int originalMaxHaplotypes = options.originalMaxHaplotypes - 1 # Ref will be added later
    cdef int maxHaplotypes = options.maxHaplotypes - 1 # Ref will be added later
    cdef int nVarsInHap = 0
    cdef int maxReadLength = options.rlen
    cdef int nVars = len(variants)
    cdef int varChunkSize = 1 # How many vars to process in one go?
    cdef int nTempOldHaps = 0
    cdef int nValidHaps = 0
    cdef int verbosity = options.verbosity
    cdef int windowSize = windowEnd - windowStart
    cdef int targetCoverage = options.coverageSamplingLevel
    cdef DiploidGenotype gt = DiploidGenotype(refHaplotype, refHaplotype)

    cdef list assemblerVars = []
    cdef list assemblerHaps = []
    cdef list outputHaps = []
    cdef int nAssemblerVars = 0   

    cdef int nHaps = 0
    cdef Haplotype thisHap
    cdef Haplotype bestHap
    cdef int nHapsProcessed = 0

    for tempVar in variants:
        if tempVar.varSource == ASSEMBLER_VAR:
            thisVar = trimLongVar(tempVar, windowStart, windowEnd) 
            assemblerVars.append(thisVar)

    nAssemblerVars = len(assemblerVars)
    #Generate all haplotypes in the region from two sources 
    for nVarsInHap from 1<= nVarsInHap <=nAssemblerVars:
        for varsThisHap in combinations(assemblerVars, nVarsInHap):
            if isHaplotypeValid(varsThisHap):       
                hap = Haplotype(chrom, windowStart, windowEnd, varsThisHap, refFile, maxReadLength, options)
                assemblerHaps.append(hap)
                
    allHaps.extend(assemblerHaps)
    nHaps = len(allHaps)
    # If nVar is small, or we've already checked, and the number of valid haplotypes is less than the max, return all combinations
    # Otherwise, we need to do some filtering. 
    if nHaps <= maxHaplotypes:
        return allHaps
    else:
        while nHapsProcessed < nHaps:
            thisHap = allHaps[nHapsProcessed]
            tempOldHaps = sorted(hapsByBestScore)
            nTempOldHaps = len(tempOldHaps)
            gt.hap2 = thisHap
            bestScoreThisHap = computeBestScoreForHaplotype(readBuffers, thisHap)#computeBestScoreForGenotype(readBuffers, gt, windowSize, targetCoverage)
            if len(hapsByBestScore) < originalMaxHaplotypes:
                heappush(hapsByBestScore, (bestScoreThisHap, thisHap))
            else:
                heappushpop(hapsByBestScore, (bestScoreThisHap, thisHap))
            nHapsProcessed += 1

        for index,(score, thisHap) in enumerate(sorted(hapsByBestScore, reverse=True)):
            if index < maxHaplotypes/2:
                outputHaps.append(thisHap)
            else:
                break

        (score, bestHap) = sorted(hapsByBestScore, reverse=True)[0]
        gt.hap1 = bestHap
        nHapsProcessed = 0

        while nHapsProcessed < nHaps:
            thisHap = allHaps[nHapsProcessed]
            gt.hap2 = thisHap
            bestScoreThisHap = computeBestScoreForGenotype(readBuffers, gt, windowSize, targetCoverage)
            
            if len(hapsByBestScore) < originalMaxHaplotypes:
                heappush(hapsByBestScore, (bestScoreThisHap, thisHap))
            else:
                heappushpop(hapsByBestScore, (bestScoreThisHap, thisHap))
            nHapsProcessed += 1
  
    for index,(score, thisHap) in enumerate(sorted(hapsByBestScore, reverse=True)):
        if index < maxHaplotypes/2:
            outputHaps.append(thisHap)
        else:
            break

    return outputHaps

###################################################################################################
