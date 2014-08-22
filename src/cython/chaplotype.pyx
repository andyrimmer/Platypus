"""
module containing various classes and functions for use in generating
and processing haplotypes.
"""

from __future__ import division

import cython
import logging
import math
import random

cimport cython
cimport variant
cimport calign
cimport fastafile
cimport samtoolsWrapper
cimport cerrormodel

from calign cimport hash_sequence, hash_sequence_multihit, hashReadForMapping
from calign cimport mapAndAlignReadToHaplotype
from fastafile cimport FastaFile
from samtoolsWrapper cimport cAlignedRead
from samtoolsWrapper cimport Read_IsReverse
from samtoolsWrapper cimport Read_IsPaired
from samtoolsWrapper cimport Read_IsProperPair
from samtoolsWrapper cimport Read_IsDuplicate
from samtoolsWrapper cimport Read_IsUnmapped
from samtoolsWrapper cimport Read_MateIsUnmapped
from samtoolsWrapper cimport Read_MateIsUnmapped
from samtoolsWrapper cimport Read_IsQCFail
from samtoolsWrapper cimport Read_IsReadOne
from samtoolsWrapper cimport Read_IsSecondaryAlignment
from variant cimport Variant
from calign cimport hash_nucs,hash_size

###################################################################################################

logger = logging.getLogger("Log")

###################################################################################################

cdef double PI = math.pi
cdef double mLTOT = -0.23025850929940459    # Minus log ten over ten

cdef extern from "math.h":
    double exp(double)
    double log(double)
    double log10(double)
    double fabs(double)
    double sqrt(double)
    double pow(double, double)

cdef extern from "stdlib.h":
    void free(void *)
    void *malloc(size_t)
    void *calloc(size_t,size_t)
    void *realloc(void *,size_t)
    void *memset(void *buffer, int ch, size_t count )

###################################################################################################

# Some nasty global variables
cdef list per_base_indel_errors = [2.9e-5, 2.9e-5, 2.9e-5, 2.9e-5, 4.3e-5, 1.1e-4, 2.4e-4, 5.7e-4, 1.0e-3, 1.4e-3] + [ 1.4e-3 + 4.3e-4*(n-10) for n in range(11,50) ]

# homopolymer indel error model
cdef bytes homopolq = bytes(''.join([chr(int(33.5 + 10*log( (idx+1)*q )/log(0.1) )) for idx,q in enumerate(per_base_indel_errors)]))

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

cdef void* my_realloc(void* thePointer, size_t newSize):
    """
    Cython wrapper. Used for profiling.
    """
    return realloc(thePointer, newSize)

###################################################################################################

cdef int computeOverlapOfReadAndHaplotype(int hapStart, int hapEnd, cAlignedRead* theRead):
    """
    Compute and return the number of bases by which a read overlaps the haplotype of interest.
    """
    cdef int readStart = theRead[0].pos
    cdef int readEnd = theRead[0].end
    cdef int overlapStart = max(hapStart, readStart)
    cdef int overlapEnd = min(hapEnd, readEnd)

    if overlapEnd > overlapStart:
        return overlapEnd - overlapStart
    else:
        return -1

###################################################################################################

@cython.final
@cython.freelist(500)
cdef class Haplotype:
    """
    Class to encapsulate a single haplotype. This will store all the
    variants pertaining to this haplotype, as well as the reference
    sequence, all supporting reads, and the start and end positions of
    this haplotype in the reference.
    """
    def __init__(self, bytes refName, int startPos, int endPos, tuple variants, FastaFile refFile, int maxReadLength, options):
        """
        Constructor. Takes a tuple of variants and a
        fasta file of the referene sequence.
        Variants are sorted on instantiation
        """
        self.refName = refName
        self.refFile = refFile
        self.variants = variants
        self.hash = -1
        self.localGapOpen = NULL
        self.haplotypeSequence = None
        self.startPos = max(0, startPos)
        self.endPos = min(endPos, self.refFile.refs[self.refName].SeqLength-1)
        self.maxReadLength = maxReadLength
        self.endBufferSize = min(2*maxReadLength, 200) # Cap the buffer size at a reasonable length
        self.verbosity = options.verbosity
        self.options = options
        self.lastIndividualIndex = -1

        cdef Variant v

        if len(variants) > 0:
            self.minVarPos = min([v.minRefPos for v in variants])
            self.maxVarPos = max([v.maxRefPos for v in variants])

            if self.minVarPos == self.maxVarPos:
                self.maxVarPos += 1
        else:
            self.minVarPos = self.startPos
            self.maxVarPos = self.endPos

        self.referenceSequence = self.refFile.getSequence(self.refName, self.startPos - self.endBufferSize, self.endPos + self.endBufferSize)

        if len(self.variants) == 0:
            self.haplotypeSequence = self.referenceSequence
        else:
            leftBuffer = self.refFile.getSequence(self.refName, self.startPos - self.endBufferSize, self.startPos)
            rightBuffer = self.refFile.getSequence(self.refName, self.endPos, self.endPos + self.endBufferSize)
            self.haplotypeSequence = leftBuffer + self.getMutatedSequence() + rightBuffer

        self.cHaplotypeSequence = self.haplotypeSequence
        self.hapLen = len(self.cHaplotypeSequence)

        #if self.referenceSequence == self.haplotypeSequence and len(self.variants) != 0:
        #    logger.error("Haplotype is broken. Var seq and ref seq are the same, and variants are %s" %(list(self.variants)))

        if self.hapLen > hash_size:
            logger.error("Haplotype with vars %s has len %s. Start is %s. End is %s. maxReadLen = %s" %(self.variants, self.hapLen, self.startPos, self.endPos, maxReadLength))
            logger.debug(self.haplotypeSequence)
            raise StandardError, "Haplotype is too long. Max allowed length is %s" %(hash_size)

        self.cHomopolQ = homopolq
        self.hapSequenceHash = NULL
        self.hapSequenceNextArray = NULL
        self.likelihoodCache = NULL
        self.lenCache = 0
        self.mapCounts = <int*>malloc( (2*(self.hapLen+maxReadLength))*sizeof(int) )
        self.mapCountsLen = 2*(self.hapLen + maxReadLength)

    def __dealloc__(self):
        """
        Clean up cache.
        """
        if self.likelihoodCache != NULL:
            my_free(self.likelihoodCache)

        if self.hapSequenceHash != NULL:
            my_free(self.hapSequenceHash)

        if self.hapSequenceNextArray != NULL:
            my_free(self.hapSequenceNextArray)

        if self.localGapOpen != NULL:
            my_free(self.localGapOpen)

        if self.mapCounts != NULL:
            my_free(self.mapCounts)

    def __copy__(self):
        """
        Make sure this never gets called for haplotypes.
        """
        raise StandardError, "Oh no! The bridge is gone!"

    def __richcmp__(Haplotype self, Haplotype other, int opCode):
        """
        Comparison function:

        Are two haplotypes equal? Only return true if the mutated
        sequences are exactly equal.
        """
        # <
        if opCode == 0:
            if self.refName < other.refName:
                return True
            elif self.refName == other.refName and self.startPos < other.startPos:
                return True
            elif self.refName == other.refName and self.startPos == other.startPos and self.haplotypeSequence < other.haplotypeSequence:
                return True
            else:
                return False
        # <=
        elif opCode == 1:
            if self.refName > other.refName:
                return False
            elif self.refName == other.refName and self.startPos > other.startPos:
                return False
            elif self.refName == other.refName and self.startPos == other.startPos and self.haplotypeSequence > other.haplotypeSequence:
                return False
            else:
                return True
        # >
        elif opCode == 4:
            if self.refName > other.refName:
                return True
            elif self.refName == other.refName and self.startPos > other.startPos:
                return True
            elif self.refName == other.refName and self.startPos == other.startPos and self.haplotypeSequence > other.haplotypeSequence:
                return True
            else:
                return False
        # >=
        elif opCode == 5:
            if self.refName < other.refName:
                return False
            elif self.refName == other.refName and self.startPos < other.startPos:
                return False
            elif self.refName == other.refName and self.startPos == other.startPos and self.haplotypeSequence < other.haplotypeSequence:
                return False
            else:
                return True
        # ==
        if opCode == 2:
            if self.refName != other.refName:
                return False
            elif self.startPos != other.startPos:
                return False
            else:
                thisSeq = self.haplotypeSequence
                otherSeq = other.haplotypeSequence
                return thisSeq == otherSeq

        # !=
        elif opCode == 3:
            if self.refName != other.refName:
                return True
            elif self.startPos != other.startPos:
                return True
            else:
                thisSeq = self.haplotypeSequence
                otherSeq = other.haplotypeSequence
                return thisSeq != otherSeq
        else:
            raise StandardError, "Op code %s not implemented in haplotype__richcmp__()" %(opCode)

    def __hash__(self):
        """
        Implementing this function allows haplotypes to be hashed, and so stored in
        a set or dictionary. The supporting reads are not included in the hashing, as
        we want two haplotypes to give the same hash id if they have the same positions
        and sequences.
        """
        if self.hash == -1:
            self.hash = hash((self.refName, self.startPos, self.endPos, self.haplotypeSequence))

        return self.hash

    cdef double* alignReads(self, int individualIndex, cAlignedRead** start, cAlignedRead** end, cAlignedRead** badReadsStart, cAlignedRead** badReadsEnd, cAlignedRead** brokenReadsStart, cAlignedRead** brokenReadsEnd, int useMapQualCap):
        """
        """
        cdef int readIndex = 0
        cdef double score = 0.0
        cdef int nReads = end - start
        cdef int nBadReads = badReadsEnd - badReadsStart
        cdef int nBrokenReads = brokenReadsEnd - brokenReadsStart
        cdef int totalReads = nReads + nBadReads + nBrokenReads
        cdef int readOverlap = 0
        cdef int readLen = 0
        cdef double* temp = NULL

        # Either first time, or new individual
        if individualIndex != self.lastIndividualIndex:

            if self.likelihoodCache == NULL:
                self.likelihoodCache = <double*>(my_malloc((totalReads+1)*sizeof(double)))
                self.lenCache = totalReads

                if self.likelihoodCache == NULL:
                    logger.error("Could not allocate haplotype cache")
                    raise StandardError, "Out of memory in cHaplotype.alignReads"
            else:
                if totalReads >= self.lenCache:
                    temp = <double*>realloc(self.likelihoodCache, 2*totalReads*sizeof(double))

                    if temp == NULL:
                        logger.error("Could not reallocate haplotype cache")
                        raise StandardError, "Out of memory in cHaplotype.alignReads"

                    self.likelihoodCache = temp
                    self.lenCache = 2*totalReads

            self.lastIndividualIndex = individualIndex

            while start != end:

                readOverlap = computeOverlapOfReadAndHaplotype(self.startPos, self.endPos, start[0])

                if Read_IsQCFail(start[0]) or readOverlap < hash_nucs:
                    self.likelihoodCache[readIndex] = 0
                else:
                    score = alignReadToHaplotype(start[0], self, useMapQualCap)
                    self.likelihoodCache[readIndex] = score

                start += 1
                readIndex += 1

            while badReadsStart != badReadsEnd:

                readOverlap = computeOverlapOfReadAndHaplotype(self.startPos, self.endPos, badReadsStart[0])

                if Read_IsQCFail(badReadsStart[0]) or readOverlap < hash_nucs:
                    self.likelihoodCache[readIndex] = 0
                else:
                    score = alignReadToHaplotype(badReadsStart[0], self, useMapQualCap)
                    self.likelihoodCache[readIndex] = score

                badReadsStart += 1
                readIndex += 1

            # It doesn't make sense to check overlap for the broken mates, as their mapping positions don't make
            # sense in this context.
            while brokenReadsStart != brokenReadsEnd:
                score = alignReadToHaplotype(brokenReadsStart[0], self, useMapQualCap)
                self.likelihoodCache[readIndex] = score
                brokenReadsStart += 1
                readIndex += 1

            self.likelihoodCache[readIndex] = 999 # End marker

        return self.likelihoodCache

    cdef inline double alignSingleRead(self, cAlignedRead* theRead, int useMapQualCap):
        """
        Returns the alignment score for a single read. If 'useMapQualCap' is True, then read likelihood
        is capped using the mapping quality of the read. Otherwise it is capped at 1e-300.
        """
        return alignReadToHaplotype(theRead, self, useMapQualCap)

    cdef char* getReferenceSequence(self, prefix = 0):
        """
        Return the refernece sequence for the region covered by this haplotype. pretty shows where the variants are.
        """
        if prefix == 0 and self.referenceSequence != None:
            return self.referenceSequence

        seqMax = self.refFile.refs[self.refName].SeqLength - 1
        self.referenceSequence = self.refFile.getSequence(self.refName, max(0, self.startPos - prefix), min(self.endPos + prefix, seqMax))
        return self.referenceSequence

    cdef char* getMutatedSequence(self):
        """
        Return the reference sequence mutated with all the variants being
        considered in this haplotype.

        Just to remind ourselves: SNPs are reported at the index of the changed
        base (0-indexed internally, and 1-indexed in VCF). Insertions are reported
        such that the index is that of the last reference base before the insertion.
        Deletions are reported such that the index is that of the first deleted base.
        """
        cdef Variant v
        cdef Variant firstVar

        if self.haplotypeSequence is None:

            currentPos = self.startPos

            # Get sequence up to one base before the first variant
            firstVar = self.variants[0]
            bitsOfMutatedSeq = [self.refFile.getSequence(self.refName, currentPos, firstVar.refPos)]
            currentPos = firstVar.refPos

            for v in self.variants:

                # Move up to one base before the next variant, if we're not already there.
                if v.refPos > currentPos:
                    bitsOfMutatedSeq.append(self.refFile.getSequence(self.refName, currentPos, v.refPos))
                    currentPos = v.refPos

                # SNP/Mult-SNP/Complex
                if v.nAdded == v.nRemoved:
                    bitsOfMutatedSeq.append(v.added)
                    currentPos += v.nRemoved

                # Arbitrary length-changing sequence replacement
                else:
                    if v.refPos == currentPos:
                        bitsOfMutatedSeq.append(self.refFile.getCharacter(self.refName, v.refPos))
                        currentPos += 1

                    currentPos += v.nRemoved
                    bitsOfMutatedSeq.append(v.added)

            # Is this ok when currentPos == endPos?
            if currentPos > self.endPos:
                logger.error("cpos = %s end pos = %s. Variants are %s" %(currentPos, self.endPos, self.variants))

            if currentPos < self.endPos:
                bitsOfMutatedSeq.append(self.refFile.getSequence(self.refName, currentPos, self.endPos))

            self.haplotypeSequence = bytes(''.join(bitsOfMutatedSeq))

        return self.haplotypeSequence

    cdef list homopolymerLengths( self ):
        """
        Return a list of homopolymer lengths for the sequence surrounding
        each variant in this haplotype.
        """
        cdef tuple vsf = self.variants
        if len( vsf ) == 0:
            return []
        else:
            return [ self.homopolymerLengthForOneVariant( v ) for v in vsf ]

    cdef int homopolymerLengthForOneVariant(self, Variant variant):
        """
        Calculate and return the length of the largest homopolymer
        touching this variant. Compute homopolymer lengths on the
        left and right, and return the largest.
        """
        varChrom = variant.refName
        varPos = variant.refPos

        leftRefSeq = self.refFile.getSequence(varChrom, varPos-20, varPos)
        rightRefSeq = self.refFile.getSequence(varChrom, varPos+1, varPos + 21)

        if len(leftRefSeq) == 0 or len(rightRefSeq) == 0:
            return 0

        leftHpSize = 0
        rightHpSize = 0

        firstLeftChar = leftRefSeq[-1]
        firstRightChar = rightRefSeq[0]

        for char in reversed(leftRefSeq):
            if char == firstLeftChar:
                leftHpSize += 1
            else:
                break

        for char in rightRefSeq:
            if char == firstRightChar:
                rightHpSize += 1
            else:
                break

        if firstLeftChar != firstRightChar:
            return max(leftHpSize, rightHpSize)
        else:
            return leftHpSize + rightHpSize

    cdef bytes getSequenceContext(self, Variant variant):
        """
        Return the sequence surrounding this variant's position.
        """
        varChrom = variant.refName
        varPos = variant.refPos
        return self.refFile.getSequence(varChrom, varPos-10, varPos + 11)

    cdef dict vcfINFO(self):
        """
        Information to go in the vcf file INFO field - a two level dictionary of variant - field - value
        This can be augmented at the individual/population level with more information, or some of it can be
        removed before printing the output

        Include
        HP = Honopolymer tract length
        NR = Number of supporting reads
        CD = Coverage depth
        """
        cdef Variant variant
        cdef dict INFO = {}

        homopolymerLengths = self.homopolymerLengths()
        vsf = self.variants

        for varIndex, variant in enumerate( self.variants ):
            HP = homopolymerLengths[varIndex]
            SC = self.getSequenceContext(variant)

            INFO[vsf[varIndex]] = {'HP':[HP], 'SC':[SC]}

        return INFO

    def __str__(self):
        """
        Generate a string representation of the haplotype. This is useful for
        debugging.
        """
        if len(self.variants) == 0:
            return '  Haplotype(*Reference*) %s:%s-%s' %(self.refName, self.startPos, self.endPos)

        vars = [str(v) for v in self.variants]
        string = "  Haplotype(" + ",".join(vars) + ") %s:%s-%s" %(self.refName, self.startPos, self.endPos)

        return string

    def __repr__(self):
        """
        The representation function. Called when printing the screen.
        """
        return self.__str__()

    cdef void annotateWithGapOpen(self):
        """
        Annotate this haplotype with a context-specific gap open penalty, using the
        homopolymer model
        """
        # Only do this one per haplotype
        if self.localGapOpen != NULL:
            return

        cdef int index = self.hapLen
        cdef char* seq = self.cHaplotypeSequence
        cdef char* errorModel = self.cHomopolQ
        cdef int homopol = -1
        cdef int homopollen = 0

        self.localGapOpen = <short*>(malloc(self.hapLen*sizeof(short)))

        homopol = -1
        homopollen = 0

        while index > 0:
            index -= 1

            if seq[index] == homopol:
                homopollen += <int>(not(not(errorModel[homopollen+1])))
            else:
                homopollen = 0

            self.localGapOpen[index] = 4*(<int>(errorModel[homopollen]) - (<int>'!'))
            homopol = seq[index];

            if homopol == 'N':
                homopol = 0

###################################################################################################

cdef double alignReadToHaplotype(cAlignedRead* read, Haplotype hap, int useMapQualCap):
    """
    This is the basic, banded-alignment routine that forms the heart of Platypus. This function decides where to anchor the
    read sequence to the specified haplotype, and calls the fastAlignmentRoutine function, which performs a banded alignment.

    If we don't anchor the read sequence to the correct part of the haplotype, then all the results, particularly for indels,
    will be rubbish.
    """
    cdef char* hapSeq = hap.cHaplotypeSequence
    cdef char* aln1 = NULL
    cdef char* aln2 = NULL
    cdef char* errorModel = NULL
    cdef int hapStart = hap.startPos - hap.endBufferSize
    cdef int gapExtend = 3
    cdef int nucprior = 2
    cdef int hapLen = hap.hapLen

    cdef char* readSeq = read.seq
    cdef char* readQuals = read.qual

    cdef int readStart = read.pos
    cdef int readLen = read.rlen
    cdef int mapQual = read.mapq

    cdef int lenOfHapSeqToTest = readLen + 15
    cdef int alignScore = 0

    cdef double probMapWrong = mLTOT*read.mapq  # A log value
    cdef double probMapRight = log(1.0 - exp(mLTOT*read.mapq)) # A log value

    # Arbitrary cap
    cdef double likelihoodCap = 0.0

    if useMapQualCap == True:
        likelihoodCap = probMapWrong
    else:
        likelihoodCap = -300 # Arbitrary cap close to min float value

    # Make sure read hash exists when needed
    if read.hash == NULL:
        hashReadForMapping(read)

    # Make sure haplotype hash exists when needed
    if hap.hapSequenceHash == NULL:
        hash_sequence_multihit(hap.cHaplotypeSequence, hap.hapLen, &hap.hapSequenceHash, &hap.hapSequenceNextArray)

    if hap.localGapOpen == NULL:
        hap.annotateWithGapOpen()

    alignScore = mapAndAlignReadToHaplotype(readSeq, readQuals, readStart, hapStart, readLen, hapLen, hap.hapSequenceHash, hap.hapSequenceNextArray, read.hash, hapSeq, gapExtend, nucprior, hap.localGapOpen, hap.mapCounts, hap.mapCountsLen)

    return max(mLTOT*alignScore + probMapRight, likelihoodCap)

###################################################################################################

cdef void logAlignmentOfReadToHaplotype(Haplotype hap, char* readSeq, char* readQuals, char* aln1, char* aln2, int mapPos, int alignScore, options):
    """
    Does exactly what it says on the tin: uses the logger to output a string representation of
    the final alignment of this read to this haplotype.
    """
    cdef list alignmentChars = [" "] * mapPos
    cdef str hapSeq = str(hap.haplotypeSequence)
    cdef str aln1Str = str(aln1)
    cdef str aln2Str = str(aln2)
    cdef int minQual = options.minBaseQual
    cdef int qualIndex = 0
    cdef int hapIndex = mapPos
    cdef int hadInsertion = False

    for theChar1,theChar2 in zip(aln1Str,aln2Str):

        if hapIndex >= len(hapSeq):
            break

        # Gap in read: deletion of ref bases.
        if theChar2 == "-":
            alignmentChars.append("-")
            hapIndex += 1

        # Gap in hap: deletion of ref bases.
        elif theChar1 == "-":
            #alignmentChars.append("-")
            #hapOffset -= 1
            qualIndex += 1
            hadInsertion = True

        elif theChar2 == hapSeq[hapIndex]:

            if hadInsertion:
                alignmentChars.append("!")
            else:
                alignmentChars.append(".")

            hadInsertion = False

            qualIndex += 1
            hapIndex += 1

        else:
            newChar = None

            if readQuals[qualIndex] > minQual:
                newChar = theChar2.lower()
            else:
                newChar = "n"

            if hadInsertion:
                newChar = newChar.upper()

            alignmentChars.append( newChar )
            hadInsertion = False

            qualIndex += 1
            hapIndex += 1

    logger.debug("".join(alignmentChars) + " Score = %s" %(alignScore))

###################################################################################################
