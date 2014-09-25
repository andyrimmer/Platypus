"""
Some utility functions and constants that don't belong anywhere
else.
"""

import math
import logging
import gzip
import os

cimport samtoolsWrapper
cimport cwindow
cimport variant

from samtoolsWrapper cimport Samfile
from samtoolsWrapper cimport cAlignedRead
from samtoolsWrapper cimport IteratorRow
from samtoolsWrapper cimport createRead
from samtoolsWrapper cimport Read_IsReverse
from samtoolsWrapper cimport Read_IsPaired
from samtoolsWrapper cimport Read_IsProperPair
from samtoolsWrapper cimport Read_IsDuplicate
from samtoolsWrapper cimport Read_IsUnmapped
from samtoolsWrapper cimport Read_MateIsUnmapped
from samtoolsWrapper cimport Read_MateIsReverse
from samtoolsWrapper cimport Read_IsQCFail
from samtoolsWrapper cimport Read_IsReadOne
from samtoolsWrapper cimport Read_IsSecondaryAlignment
from samtoolsWrapper cimport compressRead
from cwindow cimport bamReadBuffer
from chaplotype cimport Haplotype
from variant cimport Variant
from fastafile cimport FastaFile

###################################################################################################

PLATYPUS_VERSION = "0.7.9.3"

###################################################################################################

logger = logging.getLogger("Log")

###################################################################################################

cdef extern from "stdlib.h":
  void free(void *)
  void *alloca(size_t)
  void *calloc(size_t,size_t)
  int c_abs "abs" (int)

cdef extern from "math.h":
    double exp(double)
    double log(double)
    double log10(double)

cdef double PI = math.pi
cdef double mLTOT = -0.23025850929940459    # Minus log ten over ten

cdef int CIGAR_M = 0
cdef int CIGAR_I = 1
cdef int CIGAR_D = 2
cdef int CIGAR_N = 3
cdef int CIGAR_S = 4

###################################################################################################

def open(fileName, mode, compressLevel=9):
    """
    Function that allows transparent usage of dictzip, gzip and ordinary files
    """
    if fileName.endswith(".gz") or fileName.endswith(".GZ"):
        return gzip.GzipFile(fileName, mode, compressLevel)
    else:
        return file(fileName, mode)

###################################################################################################

def getSampleNamesAndLoadIterators(bamFileNames, regions, options):
    """
    Return a list of all the sample labels in the BAM files.
    """
    samples = []
    samplesByID = {}
    samplesByBAM = {}
    bamFiles = []
    nFiles = len(bamFileNames)
    cdef Samfile bamFile

    for index,fileName in enumerate(bamFileNames):

        if fileName[-4:] != ".bam":
            logger.error("Input file %s is not a BAM file" %(fileName))
            raise StandardError, "Input file %s is not a BAM file" %(fileName)

        bamFile = Samfile(fileName)
        bamFile._open('rbh', True)
        bamFiles.append(bamFile)

        try:
            theHeader = bamFile.header
            readGroupTags = theHeader['RG']

            if len(readGroupTags) > 1:
                logger.debug("Found multiple read group tags in file %s" %(fileName))

            samplesInBAM = set([x['SM'] for x in readGroupTags])
            samplesByBAM[bamFile] = list(samplesInBAM)

            for tag in readGroupTags:
                samplesByID[tag['ID']] = tag['SM']

            if len(samplesInBAM) > 1:
                logger.info("Found multiple samples in BAM file %s" %(fileName))
                logger.info("Samples are %s" %(samplesInBAM))
                samples.extend(samplesInBAM)

            elif len(samplesInBAM) == 0:
                raise StandardError, "No sample information in RG tag in file %s" %(fileName)

            else:
                sampleName = theHeader['RG'][0]['SM']
                samples.append(sampleName)
                logger.debug("Adding sample name %s, using BAM RG:SM tag in file %s" %(sampleName, fileName))

            del(theHeader)

        except StandardError, e:
            logger.debug("Error in BAM header sample parsing. Error was \n%s\n" %(e))
            sampleName = fileName.split("/")[-1][0:-4]
            samples.append(sampleName)
            samplesByBAM[bamFile] = sampleName
            logger.debug("Adding sample name %s, from BAM file %s" %(sampleName, fileName))

        # Always close the BAM file
        finally:
            bamFile.close()

    return samples,samplesByID,samplesByBAM,bamFiles

###################################################################################################

def getBAMFileNamesFromTextFile(fileName):
    """
    Reads a list of BAM file names from a text file
    """
    fileNames = []
    theTextFile = open(fileName, 'r')

    for line in theTextFile:
        line = line.strip()

        if line[-4:] == ".bam":
            fileNames.append(line)

    theTextFile.close()
    return fileNames

###################################################################################################

cdef double logFactorial(int x):
    """
    Return the logarithm of the factorial of x. Uses Stirling's approximation for large
    values.
    """
    cdef double ans = 0
    cdef int i = 0
    cdef double y = x

    if x < 15:
        for i from 1 <= i <= x:
            ans += log(i)
        return ans
    else:
        ans = y*log(y) + log(2.0*PI*y)/2 - y + (pow(y, -1))/12 - (pow(y,-3))/360 + (pow(y,-5))/1260 - (pow(y,-7))/1680 + (pow(y,-9))/1188
        return ans

###################################################################################################

def betaFunction(x,y):
    """
    For the special case of positive integers, the beta function is
    simply a ratio of factorials:

    B(x,y) = ( (x-1)!(y-1)! ) / ( (x+y-1)! )
    """
    #numerator = factorial(x-1)*factorial(y-1)
    #denominator = factorial(x+y-1)
    #return numerator/denominator
    logNumerator = logFactorial(x-1) + logFactorial(y-1)
    logDenominator = logFactorial(x+y-1)
    return exp(logNumerator - logDenominator)

###################################################################################################

cdef double logBetaFunction(int x, int y):
    """
    Return the natural logarithm of the beta function.
    """
    cdef double logNumerator = logFactorial(x-1) + logFactorial(y-1)
    cdef double logDenominator = logFactorial(x+y-1)
    return logNumerator - logDenominator

###################################################################################################

def betaPDF(alpha, beta, x):
    """
    The probability density function of a beta distribution.
    """
    return (x**(alpha-1)*(1-x)**(beta-1)) / betaFunction(alpha,beta)

###################################################################################################

def poch(x, n):
    """
    The pochhammer symbol.
    """
    if n == 0:
        return 1

    answer = x

    for i in xrange(1,n):
        answer *= (x+i)

    return answer

###################################################################################################

def logPoch(x,n):
    """
    Natural logarithm of the pochhammer symbol
    """
    if n == 0:
        return log(1)

    answer = log(x)

    for i in xrange(1,n):

        if x + i == 0:
            continue

        answer += log(abs(x+i))

    return answer

###################################################################################################

cdef double threeFTwo(int k, int n, int alpha, int beta):
    """
    The generalised hypergeometric function 3F_2(a; b; k).

    a_1 = 1
    a_2 = alpha + k + 1
    a_3 = k - n + 1
    b_1 = -beta - n + k + 2
    b_2 = k + 2
    z = 1
    """
    cdef double a_1 = 1.0
    cdef double a_2 = alpha + k + 1.0
    cdef double a_3 = k - n + 1.0
    cdef double b_1 = k + 2.0
    cdef double b_2 = -beta - n + k + 2.0
    cdef double z = 1.0
    cdef double theSum = 0.0
    cdef double lastTerm = 1.0

    theSum = lastTerm

    for i in xrange(1, abs(k-n+1)+1):
        newTerm = lastTerm * (a_2 + i - 1) * (a_3 + i - 1) / ( (b_1 + i - 1)*(b_2 + i - 1)  )
        theSum += newTerm
        lastTerm = newTerm

    return theSum

###################################################################################################

def py_betaBinomialCDF(k, n, alpha, beta):
    """
    Thin Python wrapper.
    """
    return betaBinomialCDF(k, n, alpha, beta)

###################################################################################################

cdef double betaBinomialCDF(int k, int n, int alpha, int beta):
    """
    Return the cumulative probability of the beta-binomial distribution.
    """
    if k == n:
        return 1.0

    cdef double numerator = logBetaFunction(beta+n-k-1, alpha+k+1) + log(threeFTwo(k, n, alpha, beta))
    cdef double denominator = logBetaFunction(alpha,beta) + logBetaFunction(n-k,k+2) + log(n+1)
    return max(1e-30, 1.0-exp(numerator - denominator))

###################################################################################################

cdef double binomial(int x, int size, double prob):
    """
    Optimised binomial probability function.
    """
    if x == size and prob == 1:
        return 1.0
    elif x != size and prob == 1:
        return 0.0
    elif x == 0 and prob == 0:
        return 1.0
    elif x == 0 and prob == 1:
        return 0.0
    elif x == 0 and size == 0:
        return 1.0

    cdef double logBinomCoefficient = logFactorial(size) - (logFactorial(x) + logFactorial(size-x))
    cdef double logBinomProb = x*log(prob) + (size-x)*log(1.0-prob)

    return exp(logBinomCoefficient + logBinomProb)

###################################################################################################

def nPermutations(nObjects, nChosen):
    """
    Return the number of combinations of nChosen objects selected from
    a pool of nObjects, where the order is significant, and objects are
    not replaced.
    """
    return int(round(exp( logFactorial(nObjects) - (logFactorial(nObjects - nChosen))), 2))

###################################################################################################

def nCombinations(nObjects, nChosen):
    """
    Return the number of combinations of nChosen objects selected from
    a pool of nObjects, where the order is not significant, and objects are
    not replaced.
    """
    return int(round(exp( logFactorial(nObjects) - (logFactorial(nChosen) + logFactorial(nObjects - nChosen))), 2))

###################################################################################################

def nPermutationsWithReplacement(nObjects, nChosen):
    """
    Return the number of combinations of nChosen objects selected from
    a pool of nObjects, where the order is significant, and objects are
    replaced.
    """
    return nObjects**nChosen

###################################################################################################

def nCombinationsWithReplacement(nObjects, nChosen):
    """
    Return the number of combinations of nChosen objects selected from
    a pool of nObjects, where the order is not significant, and objects are
    replaced.
    """
    return int(round((exp( logFactorial(nChosen + nObjects - 1) - (  logFactorial(nChosen) + logFactorial(nObjects - 1) ) )), 2))

###################################################################################################

cpdef tuple pruned_read_start_end(read, int minq, int minAnchor):
    """
    Calculates the start and end of the read, pruning low-quality bases at either end
    """
    cdef bytes qual = read.qual
    cdef char* cqual = qual
    cdef int readlen = len(qual)
    cdef int BASEQ = 33
    # find pruned start and end on read
    cdef int readstart = 0
    cdef int readend = readlen
    while readend > 0 and cqual[readend-1] - BASEQ < minq:
        readend -= 1
    readend -= minAnchor
    while readstart < readend and cqual[readstart] - BASEQ < minq:
        readstart += 1
    readstart += minAnchor
    return readstart, readend

###################################################################################################

cpdef tuple pruned_ref_start_end(read, int minq, int minAnchor):
    """
    Calculates the reference start and end of the read, pruning low-quality bases at either end,
    and removing an anchor sequence at either end
    """
    cdef list cigar = read.cigar
    cdef int refpos = read.pos
    cdef int readstart
    cdef int readend
    readstart, readend = pruned_read_start_end(read, minq, minAnchor)
    # deal with very low quality reads
    if readstart >= readend or not cigar: 
        return refpos, refpos
    # now convert start/end to refstart/refend
    cdef int readpos = 0
    cdef int refstart = -1
    cdef int refend = -1
    cdef int op
    cdef int le
    for op, le in cigar:
        if op == CIGAR_M or op == CIGAR_I or op == CIGAR_S:
            # does start lie in this segment?
            if readpos <= readstart < readpos+le:
                # compute ref position and store
                if op == CIGAR_M:
                    refstart = refpos + (readstart - readpos)
                else:
                    refstart = refpos
            # does end lie in this segment?
            if readpos < readend <= readpos+le:
                # compute ref position and store
                if op == CIGAR_M:
                    refend = refpos + (readend - readpos)
                else:
                    refend = refpos
            readpos += le
        if op == CIGAR_M or op == CIGAR_D or op == CIGAR_N:
            refpos += le
            # deal with the case of a deletion at the right boundary
            if readpos == readend:
                refend = refpos
    # problem
    assert refend != -1 and refstart != -1
    return refstart, refend

###################################################################################################

cdef list loadBAMData(list bamFiles, bytes chrom, int start, int end, options, list samples, dict samplesByID, dict samplesByBAM, char* refSeq):
    """
    Take a list of BAM files, and a genomic region, and reuturn a list of buffers, containing the
    reads for each BAM file in that region.

    If a single file is spcecified, and its extension is not .bam, then this is assumed to be a text
    file containing a list of input BAMs.
    """
    cdef bamReadBuffer theReadBuffer
    cdef cAlignedRead* theRead = NULL
    cdef IteratorRow readIterator
    cdef Samfile reader
    cdef list readBuffers = []
    cdef list uniqueSamples = sorted(set(samples))
    cdef list uniqueBAMs = sorted(set(bamFiles))
    cdef dict bamsBySample = {}
    cdef dict buffersBySample = {}
    cdef dict brokenMateCoordsByBAM = {}
    cdef str sample
    cdef str sampleThisRead
    cdef int totalReads = 0
    cdef int maxReads = options.maxReads
    cdef int fetchBrokenMates = options.assembleBrokenPairs
    cdef int chromID = -1
    cdef int verbosity = options.verbosity
    cdef int compressReads = options.compressReads
    cdef int qualBinSize = options.qualBinSize
    cdef char* rgID = NULL

    if len(set(samples)) == len(set(bamFiles)):
        logger.debug("There is one sample in each BAM file. No merging is required")

        for sample in uniqueSamples:
            bamsThisSample = sorted([x for x in uniqueBAMs if sample in samplesByBAM[x]])
            assert len(bamsThisSample) == 1, "Something is screwy here"
            reader = bamsThisSample[0]

            # Need to lock here when sharing BAM files
            if reader.lock is not None:
                reader.lock.acquire()

            theReadBuffer = bamReadBuffer(chrom, start, end, options)
            theReadBuffer.sample = bytes(sample)
            readIterator = reader.fetch(chrom, start, end)
            brokenMateCoords = []

            while readIterator.cnext():

                theRead = createRead(readIterator.b, 0, NULL)

                if chromID == -1:
                    chromID = theRead.chromID

                theReadBuffer.addReadToBuffer(theRead)

                if compressReads:
                    compressRead(theRead, refSeq, start, end, qualBinSize, 0)

                totalReads += 1

                if fetchBrokenMates:
                    if (not Read_IsProperPair(theRead)) or Read_IsUnmapped(theRead) or Read_MateIsUnmapped(theRead):

                        if theRead.mateChromID == theRead.chromID:
                            brokenMateCoords.append( (reader.getrname(theRead.mateChromID), theRead.mateChromID, theRead.matePos) )
                        else:
                            # Can't get these, as we don't know where they map. This happens if BAM is split and unmapped
                            # reads or broken mates are in another file.
                            if theRead.mateChromID == -1:
                                pass
                            else:
                                brokenMateCoords.append( (reader.getrname(theRead.mateChromID), theRead.mateChromID, theRead.matePos) )

                if totalReads % 250000 == 0:
                    logger.debug("Loaded %s reads in region %s:%s-%s" %(totalReads, chrom, start, end))

                if totalReads >= maxReads:
                    # Explicitly clear up memory, as Cython doesn't seem to do this
                    logger.warning("Too many reads (%s) in region %s:%s-%s. Quitting now. Either reduce --bufferSize or increase --maxReads." %(totalReads, chrom, start, end))
                    return None

            if fetchBrokenMates:
                logger.info("There are %s broken pairs in BAM %s in region %s:%s-%s" %(len(brokenMateCoords), reader.filename, chrom, start, end))
                brokenMateCoords.sort()
                queries = mergeQueries(brokenMateCoords)

                for qChrom,qStart,qEnd in queries:
                    if verbosity >= 3:
                        logger.debug("Querying broken mates %s:%s-%s" %(qChrom, qStart, qEnd))
                    readIterator = reader.fetch(qChrom, qStart, qEnd)
                    theRead = NULL

                    while readIterator.cnext():
                        if readIterator.b.core.mtid == chromID and start <= readIterator.b.core.mpos <= end:
                            theRead = createRead(readIterator.b, 0, NULL)
                            assert theRead != NULL
                            theReadBuffer.brokenMates.append(theRead)

            readBuffers.append(theReadBuffer)
            reader.close()

            # Need to release lock here when sharing BAM files
            if reader.lock is not None:
                reader.lock.release()

    # We need to merge data from multiple BAM files, or split BAM files by sample. Either way, we check the sample for each read and
    # add it to the relevant buffer.
    else:
        for sample in uniqueSamples:
            theReadBuffer = bamReadBuffer(chrom, start, end, options)
            theReadBuffer.sample = bytes(sample)
            buffersBySample[sample] = theReadBuffer

        for reader in uniqueBAMs:

            # Need to lock here when sharing BAM files
            if reader.lock is not None:
                reader.lock.acquire()

            readIterator = reader.fetch(chrom, start, end)
            brokenMateCoords = []

            while readIterator.cnext():
                theRead = createRead(readIterator.b, 1, &rgID)

                if chromID == -1:
                    chromID = theRead.chromID

                sampleThisRead = samplesByID[rgID]
                theReadBuffer = buffersBySample[sampleThisRead]
                theReadBuffer.addReadToBuffer(theRead)
                free(rgID)
                rgID = NULL

                if compressReads:
                    compressRead(theRead, refSeq, start, end, qualBinSize, 0)

                totalReads += 1

                if fetchBrokenMates:
                    if (not Read_IsProperPair(theRead)) or Read_IsUnmapped(theRead) or Read_MateIsUnmapped(theRead):
                        if theRead.mateChromID == theRead.chromID:
                            brokenMateCoords.append( (reader.getrname(theRead.mateChromID), theRead.mateChromID, theRead.matePos) )
                        else:
                            # Can't get these, as we don't know where they map. This happens if BAM is split and unmapped
                            # reads or broken mates are in another file.
                            if theRead.mateChromID == -1:
                                pass
                            else:
                                brokenMateCoords.append( (reader.getrname(theRead.mateChromID), theRead.mateChromID, theRead.matePos) )

                if totalReads >= maxReads:
                    logger.warning("Too many reads (%s) in region %s:%s-%s. Quitting now. Either reduce --bufferSize or increase --maxReads." %(totalReads, chrom, start, end))
                    return None

            if fetchBrokenMates:
                brokenMateCoords.sort()
                queries = mergeQueries(brokenMateCoords)
                logger.info("There are %s broken pairs in BAM %s in region %s:%s-%s" %(len(brokenMateCoords), reader.filename, chrom, start, end))

                for qChrom,qStart,qEnd in queries:

                    if verbosity >= 3:
                        logger.debug("Querying broken mates %s:%s-%s" %(chrom, qStart, qEnd))
                    readIterator = reader.fetch(qChrom, qStart, qEnd)
                    theRead = NULL

                    while readIterator.cnext():
                        if readIterator.b.core.mtid == chromID and start <= readIterator.b.core.mpos <= end:
                            theRead = createRead(readIterator.b, 1, &rgID)
                            assert theRead != NULL
                            sampleThisRead = samplesByID[rgID]
                            free(rgID)
                            rgID = NULL
                            theReadBuffer = buffersBySample[sampleThisRead]
                            theReadBuffer.brokenMates.append(theRead)

            reader.close()

            # Need to release lock here when sharing BAM files
            if reader.lock is not None:
                reader.lock.release()

        for theReadBuffer in buffersBySample.values():
            readBuffers.append(theReadBuffer)

    cdef list sortedBuffers = []

    for theReadBuffer in readBuffers:
        if theReadBuffer.reads.getSize() > 0:
            theReadBuffer.chromID = theReadBuffer.reads.array[0].chromID

        if theReadBuffer.brokenMates.getSize() > 0:
            theReadBuffer.sortBrokenMates()

        if not theReadBuffer.isSorted:
            theReadBuffer.sortReads()

        theReadBuffer.logFilterSummary()
        sortedBuffers.append( (theReadBuffer.sample, theReadBuffer) )

    sortedBuffers.sort()

    # Return buffers sorted by sample name
    return [x[1] for x in sortedBuffers]

###################################################################################################

cdef list mergeQueries(list coords):
    """
    Merge individual read queries into regions.
    """
    queries = []

    for mateChrom,mateChromID,matePos in coords:
        if len(queries) == 0:
            queries.append( [mateChrom, matePos, matePos + 1] )
        elif mateChrom == queries[-1][0]:
            if matePos - queries[-1][2] < 1e4 and matePos - queries[-1][1] < 1e5:
                queries[-1][2] = matePos + 1
            else:
                queries.append( [mateChrom, matePos, matePos + 1] )
        else:
            queries.append( [mateChrom, matePos, matePos + 1] )

    return queries

###################################################################################################

cdef list getAlignmentErrorsBetweenReadsAndBestHaplotypes(list readBuffers, list haplotypes):
    """
    We could allow for iterative calling of variants, by doing one pass, and then taking
    the remaining best variant candidates, and re-aligning the reads to these,
    """
    cdef bamReadBuffer theReadBuffer
    cdef list newVariants = []
    cdef tuple alignment
    cdef cAlignedRead** bufferStart
    cdef cAlignedRead** bufferEnd
    cdef cAlignedRead* theRead
    cdef Haplotype refHap = haplotypes[0]

    for theReadBuffer in readBuffers:
        bufferStart = theReadBuffer.reads.windowStart
        bufferEnd = theReadBuffer.reads.windowEnd

        while bufferStart != bufferEnd:
            theRead = bufferStart[0]
            #alignment = alignWithTraceback(theRead, refHap.cHaplotypeSequence, refHap.hapSequenceHash, refHap.startPos, refHap.varSeqLen, refHap.cLocalGapOpenQ, 3, 2, refHap.hapLen, False)
            #print alignment

###################################################################################################

cdef int isHaplotypeValid(tuple variants):
    """
    Check if this is a valid haplotype. If the variants overlap then we
    can't make a haplotype. Also optionally check if the standard form variant lies
    in the region provided - need to check at this level, rather than the variant level
    since this may depend on the reference sequence, and on the other variants in
    the haplotype
    """
    cdef int nVariants = len(variants)

    # Reference haplotype is always valid. And a single variant can't really overlap
    # with anything.
    if nVariants <= 1:
        return True

    cdef int index = 0
    cdef Variant thisVar
    cdef Variant nextVar

    # Check the normal variant starts and ends. The variants are sorted by co-ordinate, so if
    # the end of variant[i-1] occurs after the beginning of variant[i], then something is wrong.
    # Note, that the end of variant v is actually v.refPos + len(v.removed) - len(v.added).
    for index from 0 <= index < nVariants:

        thisVar = variants[index]

        # The end of the last variant can't overlap with anything.
        if index + 1 == nVariants:
            break

        nextVar = variants[index+1]

        # This should never happen
        if thisVar.minRefPos > nextVar.minRefPos:
            logger.error("Variants %s and %s are out of order. This should never happen." %(thisVar, nextVar))
            raise StandardError, "Variants out of order in haplotype!"

        # If this occurs then the haplotype is invalid. This will only happen if a deletion deletes the ref pos
        # of the next variant.
        if thisVar.maxRefPos > nextVar.minRefPos:
            return False

        #if thisVar.maxRefPos == nextVar.minRefPos:
        #    logger.warning("Funky variant combo. %s followed by %s. mins = %s,%s. maxs = %s,%s. nAdd = %s,%s. nRem = %s,%s" %(thisVar, nextVar, thisVar.minRefPos, nextVar.minRefPos, thisVar.maxRefPos, nextVar.maxRefPos, thisVar.nAdded,nextVar.nAdded,thisVar.nRemoved,nextVar.nRemoved))

        # This could be, for example, a SNP adjacent to an insertion/deletion
        elif thisVar.maxRefPos == nextVar.minRefPos:

            # Allow a SNP at the same labelled base as a deletion, as the deletion really deletes the base after
            # the reference position.
            if (thisVar.nAdded == thisVar.nRemoved) and (nextVar.nAdded < nextVar.nRemoved):
                continue

            # Allow a SNP at the same labelled base as an insertion, as the insertion really inserts after
            # the reference position.
            elif (thisVar.nAdded == thisVar.nRemoved) and (nextVar.nAdded > nextVar.nRemoved):
                continue

            # 2 SNPs/Insertions/Deletions at same position is not ok.
            else:
                return False

        # All is well with this pair of variants
        else:
            continue

    # If we get here, then all is well.
    return True

###################################################################################################

cdef Variant leftNormaliseIndel(Variant variant, FastaFile refFile, int maxReadLength):
    """
    Shift all indels as far to the left as possible. When this fails due to running out
    of sequence, then report the original position in the BAM.
    """
    cdef int nAdded = variant.nAdded
    cdef int nRemoved = variant.nRemoved

    # For SNPs, and multi-nucleotide substitutions, we don't need to do anything
    if nAdded == 1 and nRemoved == 1:
        return variant

    # Don't try to left-normalise replacements
    elif nAdded > 0 and nRemoved > 0 and nAdded != nRemoved:
        return variant

    # Don't try to left-normalise MNPs
    elif nAdded > 1 and nRemoved > 1 and nAdded == nRemoved:
        return variant

    # Hack. Leave variants alone if they occur at the extreme left end of the contigs, otherwise
    # we run out of sequence.
    if variant.refPos < 100:
        return variant

    # We look at the reference for this many bases either side of the variant, so that we can adjust its position
    # accordingly, if there are homopolymers etc.
    cdef int window = max(nAdded, nRemoved) + maxReadLength
    cdef int seqMax = refFile.refs[variant.refName].SeqLength - 1
    cdef int windowMin = max(1, variant.refPos - window)
    cdef int windowMax = min(variant.refPos + window, seqMax)

    cdef bytes bytesRefSeq = refFile.getSequence(variant.refName, windowMin, windowMax)
    cdef bytes bytesHapSeq = bytesRefSeq[0 : (variant.refPos - windowMin) + 1]

    bytesHapSeq += variant.added
    bytesHapSeq += bytesRefSeq[(variant.refPos - windowMin + nRemoved) + 1:]

    # This will invariably happen at the end of chromosomes, e.g. MT
    if bytesHapSeq[-1] != bytesRefSeq[-1] and windowMax != seqMax:
        raise StandardError, "Variant %s not correctly normalised. \nRef = %s\nHap = %s" %(variant, bytesRefSeq, bytesHapSeq)

    cdef char* refSequence = bytesRefSeq
    cdef char* hapSequence = bytesHapSeq

    cdef int lenRef = len(refSequence)
    cdef int lenHap = len(hapSequence)

    cdef int index = 0
    cdef int minLenRefHap = 0

    if lenRef < lenHap:
        minLenRefHap = lenRef
    else:
        minLenRefHap = lenHap

    # First loop forwards to find out how far to the right we can push it
    for index from 0 <= index < minLenRefHap:
        if hapSequence[index] != refSequence[index]:
            break;

    cdef int maxPos = windowMin + index + nRemoved
    cdef int newPos = -1

    cdef char refChar
    cdef char hapChar

    cdef int hapIndex = 0
    cdef int refIndex = 0

    cdef bytes newAdded = bytes("")
    cdef bytes newRemoved = bytes("")

    cdef int effecctSize = 0
    cdef int insStart = 0
    cdef int delStart = 0
    cdef int lenNewAdded = 0
    cdef int lenNewRemoved = 0

    cdef Variant newVar

    # Loop backwards through the mutated and reference sequences
    for index from 0 <= index < minLenRefHap:

        hapIndex = (lenHap - index) - 1
        refIndex = (lenRef - index) - 1

        refChar = refSequence[refIndex]
        hapChar = hapSequence[hapIndex]

        # If the sequences are the same at this point, keep going.
        # Once we hit a difference, use that position.
        if hapChar != refChar:

            if nAdded > 0:
                newPos = windowMin + (lenRef - 1) - index
                insStart = newPos - windowMin + 1 # Position of first inserted base in haplotype sequence string

                #if insStart > len(hapSequence) or insStart + nAdded > len(hapSequence):
                #    logger.info("Shit. len hap seq = %s. insStart = %s. nAdded = %s. nRemoved = %s. newAdded = %s" %(len(hapSequence), insStart, nAdded, nRemoved, hapSequence[insStart: insStart + nAdded]))

                newAdded = hapSequence[insStart: insStart + nAdded]

            if nRemoved > 0:
                newPos = windowMin + (lenRef - 1) - index - (nRemoved)
                delStart = newPos - windowMin + 1 # Position of first deleted base in reference sequence
                newRemoved = refSequence[delStart: delStart + nRemoved]

            if newPos > variant.refPos:
                logger.error("Old pos = %s new pos = %s" %(variant.refPos, newPos))
                logger.error(variant)
                logger.error(refSequence)
                logger.error(hapSequence)

            newVar = Variant(variant.refName, newPos, newRemoved, newAdded, variant.nSupportingReads, variant.varSource)
            newVar.bamMinPos = newPos
            newVar.bamMaxPos = maxPos
            newVar.bamAdded = variant.bamAdded
            newVar.bamRemoved = variant.bamRemoved

            lenNewAdded = len(newAdded)
            lenNewRemoved = len(newRemoved)

            if lenNewAdded != nAdded or lenNewRemoved != nRemoved:
                logger.error("New variant in standard format %s is broken" %(newVar))
                logger.error("Original variant was %s" %(variant))
                logger.error(refSequence)
                logger.error(hapSequence)
                raise StandardError, "Error in variant conversion to standard format"

            return newVar

    # If we get to here, and we haven't returned, then it's almost certainly an error
    logger.warning("Could not left-normalise variant %s. Using position as reported in BAM" %(variant))
    logger.info("\n")
    logger.info(variant)
    logger.info(bytesRefSeq)
    logger.info(bytesHapSeq)
    logger.info("\n")
    return variant

###################################################################################################

def getRegions(options):
    """
    Extract the regions to act on. Either it's a list of regions in the format "chrX:start-end". If start-end
    isn't specified, we do the whole chr. Or it's None, in which case we read the file header and do everything.
    """
    if options.refFile.endswith(".gz") or options.refFile.endswith(".bz2") or options.refFile.endswith(".bgz"):
        logger.error("Reference file-name (%s) looks like a compressed file-name. Please un-compress the reference FASTA file before running Platypus" %(options.refFile))
        raise StandardError, "Invalid reference FASTA file supplied"

    cdef FastaFile refFile = FastaFile(options.refFile, options.refFile + ".fai", parseNCBI = options.parseNCBI)
    fileName = None

    if len(options.bamFiles) == 1 and options.bamFiles[0][-4:] != ".bam":

        theTextFile = open(options.bamFiles[0], 'r')

        for line in theTextFile:
            line = line.strip()
            if line[-4:] == ".bam":
                fileName = line
                break
    else:
        fileName = options.bamFiles[0]

    if fileName[-4:] != ".bam":
        logger.error("Input file %s is not a BAM file" %(fileName))
        raise StandardError, "Input file %s is not a BAM file" %(fileName)

    file = samtoolsWrapper.Samfile(fileName)
    file._open('rbh', loadIndex=True)
    finalRegions = []
    regions = []

    if options.regions is not None and os.path.exists(options.regions[0]):

        # Text file with regions in format chr:start-end
        if options.regions[0].endswith(".txt"):

            logger.info("Interpreting --regions argument (%s) as a text file with regions in the format chr:start-end" %(options.regions[0]))

            with open(options.regions[0], 'r') as theFile:
                for line in theFile:
                    chrom,region = line.split(":")
                    start = int(region.split("-")[0])
                    end = int(region.split("-")[1])

                    regions.append( (chrom,start,end) )

        # BED file with regions in format chr\start\tend
        elif options.regions[0].endswith(".bed"):

            logger.info("Interpreting --regions argument (%s) as a BED file with regions in the format chr\tstart\tend" %(options.regions[0]))

            with open(options.regions[0], 'r') as theFile:
                for line in theFile:
                    try:
                        cols = line.split("\t")
                        chrom = cols[0]
                        start = int(cols[1])
                        end = int(cols[2])
                        regions.append( (chrom,start,end) )
                    except:
                        logger.debug("Could not parse line in regions file (%s). Skipping..." %(options.regions[0]))
                        logger.debug("Line was %s" %(line))
                        continue

    elif options.regions == None:

        try:
            header = file.header
            regions = [ (d['SN'], 1, d['LN']) for d in header['SQ'] ]
        except:

            for region,regionTuple in refFile.refs.iteritems():
                regions.append((region, 1, regionTuple.SeqLength))

    else:
        for region in options.regions:

            split = region.split(":")
            chrom = bytes(split[0])

            if len( split ) == 2 :
                [ start, end ] = split[1].split("-")
                regions.append((chrom,int(start),int(end)))

                if regions[-1][2] - regions[-1][1] > 1e9:
                    logger.error("Input region (%s) is too long. Try again" %(region))
                    raise StandardError, "Invalid input region: %s" (region)

            elif len(split) == 1:
                start = 1
                try:
                    header = file.header
                    pivot = dict(zip([d['SN'] for d in header['SQ']], [d['LN'] for d in header['SQ']]))
                    end = pivot[chrom]
                    regions.append((chrom,int(start),int(end)))
                except:
                    regions = []
                    for region,regionTuple in refFile.refs.iteritems():
                        if region == chrom:
                            regions.append((region, 1, regionTuple.SeqLength))
            else:
                regions.append((chrom,None,None))

    if len(regions) == 0:
        logger.error("Platypus found no regions to search. Check that you are using the correct reference FASTA file or have specified the 'regions' argument correctly")
    elif len(regions) < 100:
        logger.debug("The following regions will be searched: %s" %(regions))
    else:
        logger.debug("%s regions will be searched" % len(regions))

    file.close()
    finalRegions = []

    # Break-down large regions to save on memory, and to prevent windows
    # becoming ridiculous.
    for region in regions:

        regionLen = None

        try:
            regionLen = refFile.refs[region[0]].SeqLength
        except KeyError:
            logger.debug("Reference sequence %s is not represented in reference fasta file. Skipping this sequence" %(region[0]))
            continue

        # An invalid region.
        if region[1] is not None and region[1] > regionLen:
            logger.warning("Skipping region %s, as start co-ord (%s) is > region length (%s)" %(region,region[1],regionLen))
            continue

        if region[1] is None or region[2] is None:

            for i in range(1, region[2], options.bufferSize):
                start = i
                end = min(i+options.bufferSize, region[2])
                finalRegions.append( (region[0], start, end) )

        elif region[2] - region[1] > options.bufferSize:
            for i in range(region[1], region[2], options.bufferSize):
                start = i
                end = min(i+options.bufferSize, region[2])
                finalRegions.append( (region[0], start, end) )
        else:
            finalRegions.append(region)

    if options.verbosity >= 3:
        logger.debug("The following genomic regions will be searched: %s" %(finalRegions))

    return finalRegions

###################################################################################################
