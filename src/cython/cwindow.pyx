"""
Fast cython implementation of some windowing functions.
"""

# This seems to affect c-division as well now...
#from __future__ import division

import cython

cimport htslibWrapper

import logging
import math

from htslibWrapper cimport destroyRead
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
from htslibWrapper cimport Read_SetCompressed
from htslibWrapper cimport Read_IsCompressed
from htslibWrapper cimport compressRead
from htslibWrapper cimport uncompressRead
from bisect import bisect_left
from heapq import heappush,heappop

###################################################################################################

logger = logging.getLogger("Log")

###################################################################################################

cdef int LOW_QUAL_BASES=0
cdef int UNMAPPED_READ=1
cdef int MATE_UNMAPPED=2
cdef int MATE_DISTANT=3
cdef int SMALL_INSERT=4
cdef int DUPLICATE=5
cdef int LOW_MAP_QUAL=6

###################################################################################################

cdef extern from "stdlib.h":
    void free(void *)
    void *malloc(size_t)
    void *calloc(size_t,size_t)
    void *realloc(void *,size_t)
    void qsort(void*, size_t, size_t, int(*)(void*, void*))

###################################################################################################

cdef extern from "math.h":
    double exp(double)
    double log(double)
    double log10(double)
    double fabs(double)
    int abs(int)

###################################################################################################

cdef extern from "string.h":
  ctypedef int size_t
  void *memcpy(void *dst,void *src,size_t len)
  int strncmp(char *s1,char *s2,size_t len)
  char *strncpy(char *dest,char *src, size_t len)
  size_t strlen(char *s)
  int memcmp( void * s1, void *s2, size_t len )

###################################################################################################

@cython.profile(False)
cdef inline int readPosComp(const void* x, const void* y) nogil:
    """
    Comparison function for use in qsort, to sort reads by their start positions and
    then end positions.
    """
    cdef cAlignedRead** readOne = <cAlignedRead**>(x)
    cdef cAlignedRead** readTwo = <cAlignedRead**>(y)

    #if readOne[0].pos != readTwo[0].pos
    return readOne[0].pos - readTwo[0].pos
    #else:
    #    return readOne[0].end - readTwo[0].end

###################################################################################################

@cython.profile(False)
cdef int readMatePosComp(const void* x, const void* y):
    """
    Comparison function for use in qsort, to sort reads by their mate positions, as long as
    the mates are on the same chromosome.
    """
    cdef cAlignedRead** readOne = <cAlignedRead**>(x)
    cdef cAlignedRead** readTwo = <cAlignedRead**>(y)

    # Sorting is broken for reads with mates on different chromosomes
    assert readOne[0].mateChromID == readTwo[0].mateChromID
    return  readOne[0].matePos - readTwo[0].matePos

###################################################################################################

cdef class ReadArray:
    """
    Simple structure to wrap a raw C array, with some bounds
    checking.
    """
    def __init__(self, int size):
        """
        Allocate an array of size 'size', with initial values
        'init'.
        """
        self.array = <cAlignedRead**>(malloc(size*sizeof(cAlignedRead*)))
        assert self.array != NULL, "Could not allocate memory for ReadArray"
        
        self.__size = 0 # We don't put anything in here yet, just allocate memory
        self.__capacity = size
        self.__longestRead = 0
        self.windowStart = NULL
        self.windowEnd = NULL
        
        # Always initialise to NULL
        for index from 0 <= index < size:
            self.array[index] = NULL
    
    def __dealloc__(self):
        """
        Free memory
        """
        cdef int index = 0
        if self.array != NULL:
            for index from 0 <= index < self.__size:
                if self.array[index] != NULL:
                    destroyRead(self.array[index])
                    self.array[index] = NULL
            
            free(self.array)
    
    cdef int getSize(self):
        """
        Return the size of the array
        """
        return self.__size
    
    cdef void append(self, cAlignedRead* value):
        """
        Append a new value to the array, re-allocating if necessary.
        """
        cdef cAlignedRead** temp = NULL

        if self.__size == self.__capacity:
            temp = <cAlignedRead**>(realloc(self.array, 2*sizeof(cAlignedRead*)*self.__capacity))

            if temp == NULL:
                raise StandardError, "Could not re-allocate ReadArray"
            else:
                self.array = temp
                self.__capacity *= 2

        self.array[self.__size] = value
        self.__size += 1

        cdef int readStart = value.pos
        cdef int readEnd = value.end
        cdef int readLength = readEnd - readStart

        if readLength > self.__longestRead:
            self.__longestRead = readLength

    cdef int countReadsCoveringRegion(self, int start, int end):
        """
        Return the number of reads which overlap this region, where 'end' is not
        included, and 'start' is.
        """
        cdef cAlignedRead** rStart
        cdef cAlignedRead** rEnd
        cdef int firstOverlapStart = -1
        cdef int startPosOfReads = -1
        cdef int endPosOfReads = -1

        # Set pointers for good reads
        if self.__size == 0:
            return 0
        else:
            firstOverlapStart = max(1, start - self.__longestRead)
            startPosOfReads = bisectReadsLeft(self.array, firstOverlapStart, self.__size)
            endPosOfReads = bisectReadsLeft(self.array, end, self.__size)

            while startPosOfReads < self.__size and self.array[startPosOfReads].end <= start:
                startPosOfReads += 1

            rStart = self.array + startPosOfReads
            rEnd = min(self.array + endPosOfReads, self.array + self.__size)

            if startPosOfReads > endPosOfReads:
                logger.error("Start pos = %s. End pos = %s. Read start pos = %s. end pos = %s" %(start, end, startPosOfReads, endPosOfReads))
                logger.error("There are %s reads here." %(self.__size))
                raise StandardError, "This should never happen. Read start pointer > read end pointer!!"

            return (rEnd - rStart)

    cdef void setWindowPointers(self, int start, int end):
        """
        Set the pointers 'windowStart' and 'windowEnd' to
        point to the relevant first and last +1 reads as specified by
        the co-ordinates.
        """
        cdef int firstOverlapStart = -1
        cdef int startPosOfReads = -1
        cdef int endPosOfReads = -1

        # Set pointers for good reads
        if self.__size == 0:
            self.windowStart = self.array
            self.windowEnd = self.array
        else:
            firstOverlapStart = max(1, start - self.__longestRead)
            startPosOfReads = bisectReadsLeft(self.array, firstOverlapStart, self.__size)
            endPosOfReads = bisectReadsLeft(self.array, end, self.__size)

            while startPosOfReads < self.__size and self.array[startPosOfReads].end <= start:
                startPosOfReads += 1

            self.windowStart = self.array + startPosOfReads
            self.windowEnd = min(self.array + endPosOfReads, self.array + self.__size)

            if startPosOfReads > endPosOfReads:
                logger.info("Start pos = %s. End pos = %s. Read start pos = %s. end pos = %s" %(start, end, startPosOfReads, endPosOfReads))
                logger.info("There are %s reads here." %(self.__size))
                raise StandardError, "This should never happen. Read start pointer > read end pointer!!"

    cdef void setWindowPointersBasedOnMatePos(self, int start, int end):
        """
        Set the pointers 'windowStart' and 'windowEnd' to
        point to the relevant first and last +1 reads as specified by
        the co-ordinates of the mates of the reads in this array.
        """
        cdef int firstOverlapStart = -1
        cdef int startPosOfReads = -1
        cdef int endPosOfReads = -1

        # Set pointers for good reads
        if self.__size == 0:
            self.windowStart = self.array
            self.windowEnd = self.array
        else:
            firstOverlapStart = max(1, start - self.__longestRead)
            startPosOfReads = bisectReadsLeft(self.array, firstOverlapStart, self.__size, 1)
            endPosOfReads = bisectReadsLeft(self.array, end, self.__size, 1)

            #while startPosOfReads < self.__size and self.array[startPosOfReads].end <= start:
            #    startPosOfReads += 1

            self.windowStart = self.array + startPosOfReads
            self.windowEnd = min(self.array + endPosOfReads, self.array + self.__size)

            if startPosOfReads > endPosOfReads:
                logger.info("Start pos = %s. End pos = %s. Read start pos = %s. end pos = %s" %(start, end, startPosOfReads, endPosOfReads))
                logger.info("There are %s reads here." %(self.__size))
                raise StandardError, "This should never happen. Read start pointer > read end pointer!!"

    cdef int getLengthOfLongestRead(self):
        """
        Simple wrapper function.
        """
        return self.__longestRead

###################################################################################################

cdef int bisectReadsLeft(cAlignedRead** reads, int testPos, int nReads, int testMatePos=0):
    """
    Specialisation of bisection algorithm for array of
    read pointers.
    """
    cdef int low = 0
    cdef int high = nReads
    cdef int mid = 0

    while low < high:

        mid = (low + high) / 2

        if not testMatePos:
            if reads[mid].pos < testPos:
                low = mid + 1
            else:
                high = mid
        else:
            if reads[mid].matePos < testPos:
                low = mid + 1
            else:
                high = mid

    return low

###################################################################################################

cdef int bisectReadsRight(cAlignedRead** reads, int testPos, int nReads, int testMatePos=0):
    """
    Specialisation of bisection algorithm for array of
    read pointers.
    """
    cdef int low = 0
    cdef int high = nReads
    cdef int mid = 0

    while low < high:

        mid = (low + high) / 2

        if not testMatePos:
            if testPos < reads[mid].pos:
                high = mid
            else:
                low = mid + 1
        else:
            if testPos < reads[mid].matePos:
                high = mid
            else:
                low = mid + 1

    return low

###################################################################################################

cdef int checkAndTrimRead(cAlignedRead* theRead, cAlignedRead* theLastRead, int minGoodQualBases, int* filteredReadCountsByType, int minMapQual, int minBaseQual, int minFlank, int trimOverlapping, int trimAdapter, int trimReadFlank, int trimSoftClipped):
    """
    Performs various quality checks on the read, and trims read (i.e. set q-scores to zero). Returns
    true if read is ok, and false otherwise.
    """
    if Read_IsSecondaryAlignment(theRead):
        Read_SetQCFail(theRead)
        return False
    
    if theRead.mapq < minMapQual:
        filteredReadCountsByType[LOW_MAP_QUAL] += 1
        Read_SetQCFail(theRead)
        return False

    # Filter out low quality reads, i.e. those with < minGoodQualBases with Q scores >= 20
    cdef int nQualsBelowMin = 0
    cdef int index = 0

    for index from 0 <= index < theRead.rlen:
        if theRead.qual[index] < minBaseQual:
            nQualsBelowMin += 1

    if theRead.rlen - nQualsBelowMin < minGoodQualBases:
        filteredReadCountsByType[LOW_QUAL_BASES] += 1
        Read_SetQCFail(theRead)
        return False

    # Remove unmapped reads
    if Read_IsUnmapped(theRead):
        filteredReadCountsByType[UNMAPPED_READ] += 1
        Read_SetQCFail(theRead)
        return False

    # Remove broken pairs, i.e. pairs where the mate is mapped to a different chromosome or the
    # mate is unmapped
    if filteredReadCountsByType[MATE_UNMAPPED] != -1:

        if Read_IsPaired(theRead) and Read_MateIsUnmapped(theRead):
            filteredReadCountsByType[MATE_UNMAPPED] += 1
            return False

    if filteredReadCountsByType[MATE_DISTANT] != -1:

        if Read_IsPaired(theRead) and (theRead.chromID != theRead.mateChromID or (not Read_IsProperPair(theRead))):
            filteredReadCountsByType[MATE_DISTANT] += 1
            return False

    # If the insert size is < read length, then we almost certainly have adapter contamination, in which case we should
    # skip these reads, as they may be mapped to the wrong location in the genome.
    if filteredReadCountsByType[SMALL_INSERT] != -1:

        if Read_IsPaired(theRead) and (theRead.insertSize != 0 and abs(theRead.insertSize) < theRead.rlen):
            filteredReadCountsByType[SMALL_INSERT] += 1
            Read_SetQCFail(theRead)
            return False

    # Check if this read is actually a duplicate. TODO: store library tag and check.
    if filteredReadCountsByType[DUPLICATE] != -1:
        
        if Read_IsDuplicate(theRead):
            filteredReadCountsByType[DUPLICATE] += 1
            Read_SetQCFail(theRead)
            return False
        elif theLastRead != NULL:
            if theRead.pos == theLastRead.pos and theRead.rlen == theLastRead.rlen:

                # For paired reads, check mate's position
                if Read_IsPaired(theRead):

                    if theLastRead.matePos == theRead.matePos:
                        filteredReadCountsByType[DUPLICATE] += 1
                        Read_SetQCFail(theRead)
                        return False

                # For single reads, just check pos and length of reads
                else:
                    filteredReadCountsByType[DUPLICATE] += 1
                    Read_SetQCFail(theRead)
                    return False

    ## Any read that gets passed here will be used, but low quality tails will be trimmed, and
    ## any overlapping pairs, from small fragments, will have the overlapping bits trimmed.
     
    # Trim low-quality tails for forward reads
    if not Read_IsReverse(theRead):

        for index from 1 <= index <= theRead.rlen:
            if index < trimReadFlank or theRead.qual[theRead.rlen - index] < 5:
                theRead.qual[theRead.rlen - index] = 0
            else:
                break
    # Trim low-quality tails for reverse reads
    else:
        for index from 0 <= index < theRead.rlen:

            if index < trimReadFlank or theRead.qual[index] < 5:
                theRead.qual[index] = 0
            else:
                break
    
    cdef int absIns = abs(theRead.insertSize)

    # Trim overlapping part of forward read, in pairs where the read length is greater than the insert size
    # N.B Insert size is from start of forward read to end of reverse read, i.e. fragment size. This is done to
    # remove duplicate information, which gives systematic errors when pcr errors have occured in library prep.

    if trimOverlapping == 1 and (Read_IsPaired(theRead) and absIns > 0 and (not Read_IsReverse(theRead)) and Read_MateIsReverse(theRead) and absIns < 2*theRead.rlen):
        for index from 1 <= index <= min(theRead.rlen, (2*theRead.rlen - theRead.insertSize) + 1):
            theRead.qual[theRead.rlen - index] = 0
    
    # Trim the end of any read where the insert size is < read length. If these have not been
    # already filtered out then they need trimming, as adapter contamination will cause a
    # high FP rate otherwise.
    if trimAdapter == 1 and (Read_IsPaired(theRead) and absIns > 0 and absIns < theRead.rlen):
        
        if Read_IsReverse(theRead):
            for index from 1 <= index < theRead.rlen - absIns + 1:
                theRead.qual[theRead.rlen - index] = 0
        else:
            for index from absIns <= index < theRead.rlen:
                theRead.qual[index] = 0
    
    
    # Check for soft-clipping (present in BWA reads, for example, but not Stampy). Soft-clipped
    # sequences should be set to QUAL = 0, as they may include contamination by adapters etc.
    cdef int cigarIndex = 0
    cdef int j = 0
    cdef int cigarOp = -1
    cdef int cigarLen = -1
    
    if trimSoftClipped == 1:
        index = 0
        
        for cigarIndex in range(theRead.cigarLen):
            
            cigarOp = theRead.cigarOps[2*cigarIndex]
            cigarLen = theRead.cigarOps[2*cigarIndex + 1]
            
            # Skip good sequence. 0 is match. 1 is insertion.
            if cigarOp == 0 or cigarOp == 1:
                index += theRead.cigarOps[2*cigarIndex + 1]
            
            # 4 is soft-clipping
            elif cigarOp == 4:
                # Set quals to zero across all sequence flagged as soft-clipped
                for j in range(theRead.cigarOps[2*cigarIndex + 1]):
                    theRead.qual[index] = 0
                    index += 1
    
    return True

###################################################################################################

cdef class bamReadBuffer(object):
    """
    Utility class for bufffering reads from a single BAM file, so we only make a single pass
    through the data in each BAM in the loop through windows.
    """
    def __init__(self, char* chrom, int start, int end, options):
        """
        Constructor.
        """
        cdef int initialSize = max(100, ( (end - start) / options.rlen))
        self.isSorted = True
        self.reads = ReadArray(initialSize)
        self.badReads = ReadArray(initialSize)
        self.brokenMates = ReadArray(initialSize)
        self.filteredReadCountsByType = <int*>(calloc(7, sizeof(int)))
        self.chrom = chrom
        self.startBase = start
        self.endBase = end
        self.maxReads = options.maxReads
        self.minBaseQual = options.minBaseQual
        self.minFlank = options.minFlank
        self.trimReadFlank = options.trimReadFlank
        self.minMapQual = options.minMapQual
        self.minGoodBases = options.minGoodQualBases
        self.verbosity = options.verbosity
        self.trimOverlapping = options.trimOverlapping
        self.trimAdapter = options.trimAdapter
        self.trimSoftClipped = options.trimSoftClipped
        self.lastRead = NULL

        if options.filterDuplicates == 0:
            self.filteredReadCountsByType[DUPLICATE] = -1

        if options.filterReadsWithUnmappedMates == 0:
            self.filteredReadCountsByType[MATE_UNMAPPED] = -1

        if options.filterReadsWithDistantMates == 0:
            self.filteredReadCountsByType[MATE_DISTANT] = -1

        if options.filterReadPairsWithSmallInserts == 0:
            self.filteredReadCountsByType[SMALL_INSERT] = -1

    cdef void logFilterSummary(self):
        """
        Useful debug information about which reads have been filtered
        out.
        """
        if self.verbosity >= 3:
            logger.debug("Sample %s has %s good reads" %(self.sample, self.reads.getSize()))
            logger.debug("Sample %s has %s bad reads" %(self.sample, self.badReads.getSize()))
            logger.debug("Sample %s has %s broken mates" %(self.sample, self.brokenMates.getSize()))
            logger.debug("N low map quality reads = %s" %(self.filteredReadCountsByType[LOW_MAP_QUAL]))
            logger.debug("N low qual reads = %s" %(self.filteredReadCountsByType[LOW_QUAL_BASES]))
            logger.debug("N un-mapped reads = %s" %(self.filteredReadCountsByType[UNMAPPED_READ]))
            logger.debug("N reads with unmapped mates = %s" %(self.filteredReadCountsByType[MATE_UNMAPPED]))
            logger.debug("N reads with distant mates = %s" %(self.filteredReadCountsByType[MATE_DISTANT]))
            logger.debug("N reads pairs with small inserts = %s" %(self.filteredReadCountsByType[SMALL_INSERT]))
            logger.debug("N duplicate reads = %s" %(self.filteredReadCountsByType[DUPLICATE]))

            if self.trimOverlapping == 1:
                logger.debug("Overlapping segments of read pairs were clipped")
            else:
                logger.debug("Overlapping segments of read pairs were not clipped")

            if self.trimAdapter == 1:
                logger.debug("Overhanging bits of reads (poss. adapter sequences) were clipped")
            else:
                logger.debug("Overhanging bits of reads (poss. adapter sequences) were not clipped")

    def __dealloc__(self):
        """
        Clean up memory
        """
        free(self.filteredReadCountsByType)

    cdef void addReadToBuffer(self, cAlignedRead* theRead):
        """
        Add a new read to the buffer, making sure to re-allocate
        memory when necessary.
        """
        # Temp variable for checking that re-alloc works
        cdef int readOk = 0
        cdef int minGoodBasesThisRead = 0
        cdef int readStart = -1
        cdef int readEnd = -1
        cdef int readLength = -1

        if theRead == NULL:
            if self.verbosity >= 3:
                logger.debug("Found null read")
            return
        else:
            
            # TODO: Check that this works for duplicates when first read goes into bad reads pile...
            if self.lastRead != NULL:
                readOk = checkAndTrimRead(theRead, self.lastRead, self.minGoodBases, self.filteredReadCountsByType, self.minMapQual, self.minBaseQual, self.minFlank, self.trimOverlapping, self.trimAdapter, self.trimReadFlank, self.trimSoftClipped)

                if self.lastRead.pos > theRead.pos:
                    self.isSorted = False
            else:
                readOk = checkAndTrimRead(theRead, NULL, self.minGoodBases, self.filteredReadCountsByType, self.minMapQual, self.minBaseQual, self.minFlank, self.trimOverlapping, self.trimAdapter, self.trimReadFlank, self.trimSoftClipped)
            
            self.lastRead = theRead

            # Put read into bad array, and remove from good array
            if not readOk:
                self.badReads.append(theRead)

            # Put read into good array
            else:
                self.reads.append(theRead)

    cdef int countAlignmentGaps(self):
        """
        Count and return the number of indels seen
        by the mapper in all good and bad reads.
        """
        cdef cAlignedRead** start = self.reads.windowStart
        cdef cAlignedRead** end = self.reads.windowEnd
        cdef cAlignedRead** bStart = self.badReads.windowStart
        cdef cAlignedRead** bEnd = self.badReads.windowEnd

        cdef int nGaps = 0
        cdef int i = 0

        while start != end:
            for i from 0 <= i < start[0].cigarLen:
                if 1 <= start[0].cigarOps[2*i] <= 4:
                    nGaps += 1

            start += 1

        while bStart != bEnd:
            for i from 0 <= i < bStart[0].cigarLen:
                if 1 <= bStart[0].cigarOps[2*i] <= 4:
                    nGaps += 1
            bStart += 1

        return nGaps

    cdef int countImproperPairs(self):
        """
        Count and return the number of reads (Good and bad) that
        are members of improper pairs.
        """
        cdef cAlignedRead** start = self.reads.windowStart
        cdef cAlignedRead** end = self.reads.windowEnd
        cdef cAlignedRead** bStart = self.badReads.windowStart
        cdef cAlignedRead** bEnd = self.badReads.windowEnd

        cdef int nImproper = 0

        while start != end:
            if not Read_IsProperPair(start[0]):
                nImproper += 1
            start += 1

        while bStart != bEnd:
            if not Read_IsProperPair(bStart[0]):
                nImproper += 1
            bStart += 1

        return nImproper

    cdef int countReadsCoveringRegion(self, int start, int end):
        """
        Return the number of 'good' reads covering this region.
        """
        return self.reads.countReadsCoveringRegion(start, end)

    cdef void setWindowPointers(self, int start, int end, int refStart, int refEnd, char* refSeq, int qualBinSize):
        """
        Set the windowStart and windowEnd pointers to point to the first
        and last+1 reads covering this window.
        """
        self.reads.setWindowPointers(start, end)
        self.badReads.setWindowPointers(start, end)
        self.brokenMates.setWindowPointersBasedOnMatePos(start, end)

        cdef cAlignedRead** theStart = NULL
        cdef cAlignedRead** theEnd = NULL

        # If the reads in this window are compressed, then uncompress them
        theStart = self.reads.windowStart
        theEnd = self.reads.windowEnd

        if theStart != theEnd and Read_IsCompressed(theStart[0]):

            while theStart != theEnd:
                uncompressRead(theStart[0], refSeq, refStart, refEnd, qualBinSize)
                theStart += 1

            theStart = self.badReads.windowStart
            theEnd = self.badReads.windowEnd

            while theStart != theEnd:
                uncompressRead(theStart[0], refSeq, refStart, refEnd, qualBinSize)
                theStart += 1

            theStart = self.brokenMates.windowStart
            theEnd = self.brokenMates.windowEnd

            while theStart != theEnd:
                uncompressRead(theStart[0], refSeq, refStart, refEnd, qualBinSize)
                theStart += 1

    cdef void recompressReadsInCurrentWindow(self, int refStart, int refEnd, char* refSeq, int qualBinSize, int compressReads):
        """
        Set the windowStart and windowEnd pointers to point to the first
        and last+1 reads covering this window.
        """
        cdef cAlignedRead** theStart = NULL
        cdef cAlignedRead** theEnd = NULL

        # If the reads in this window are compressed, then uncompress them
        theStart = self.reads.windowStart
        theEnd = self.reads.windowEnd

        while theStart != theEnd:

            # Compresss reads is required
            if compressReads:
                compressRead(theStart[0], refSeq, refStart, refEnd, qualBinSize, 1)

            # Otherwise, just clear the hash
            if theStart[0].hash != NULL:
                free(theStart[0].hash)
                theStart[0].hash = NULL

            theStart += 1

        theStart = self.badReads.windowStart
        theEnd = self.badReads.windowEnd

        while theStart != theEnd:

            # Compresss reads is required
            if compressReads:
                compressRead(theStart[0], refSeq, refStart, refEnd, qualBinSize, 1)

            # Otherwise, just clear the hash
            if theStart[0].hash != NULL:
                free(theStart[0].hash)
                theStart[0].hash = NULL

            theStart += 1

        theStart = self.brokenMates.windowStart
        theEnd = self.brokenMates.windowEnd

        while theStart != theEnd:

            # Compresss reads is required
            if compressReads:
                compressRead(theStart[0], refSeq, refStart, refEnd, qualBinSize, 1)

            # Otherwise, just clear the hash
            if theStart[0].hash != NULL:
                free(theStart[0].hash)
                theStart[0].hash = NULL

            theStart += 1

    cdef void sortReads(self):
        """
        Sort the contents of the reads array by position
        """
        if not self.isSorted:
            if self.reads.getSize() > 0:
                qsort(self.reads.array, self.reads.getSize(), sizeof(cAlignedRead**), readPosComp)

            if self.badReads.getSize() > 0:
                qsort(self.badReads.array, self.badReads.getSize(), sizeof(cAlignedRead**), readPosComp)

            self.isSorted = True

    cdef void sortBrokenMates(self):
        """
        Sort the contents of the brokenMates array by the co-ordinates of their
        mates.
        """
        qsort(self.brokenMates.array, self.brokenMates.getSize(), sizeof(cAlignedRead**), readMatePosComp)

###################################################################################################

