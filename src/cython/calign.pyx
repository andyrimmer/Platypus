"""
Cython module containing glue code for the alignment routines used in
Platypus.
"""

from __future__ import division

cimport cython

import logging
import samtoolsWrapper
cimport cerrormodel
cimport samtoolsWrapper

from samtoolsWrapper cimport cAlignedRead

logger = logging.getLogger("Log")

###################################################################################################

# set the size of the hash.  Ensure that hash_size == math.pow(4,hash_nucs)
cdef int hash_nucs = 7
cdef int hash_size = 4**hash_nucs
cdef int max_sequence_length = hash_size
cdef int mask = (1 << 2*hash_nucs) - 1

ctypedef long long size_t
cdef double mLTOT = -0.23025850929940459    # Minus log ten over ten

###################################################################################################

cdef extern from "align.h":
    int fastAlignmentRoutine(char* seq1, char* seq2, char* qual2, int len1, int len2, int gapextend, int nucprior, short* localgapopen)

###################################################################################################

cdef extern from "stdlib.h":
    void free(void *)
    void *malloc(size_t)
    void *calloc(size_t,size_t)
    void *memset(void *buffer, int ch, size_t count )

###################################################################################################

cdef extern from "string.h":
    int strncmp (char*, char*, size_t)

cdef extern from "math.h":
    double exp(double)
    double log(double)

###################################################################################################

cdef unsigned int my_hash(char* seq):
    """
    Encodes first hash_nucs nucleotides (A,C,G,T or a,c,g,t) into an integer
    """
    cdef int i, c
    cdef unsigned int h = 0

    for i in range(hash_nucs):
        # Just a simple hash function.
        c = seq[i] & 7   # a,A->1  c,C->3  g,G->7  t,T->4

        if c == 7:
            c = 2

        h = (h << 2) + <unsigned int>(c & 3)

    return h

###################################################################################################

cdef void hash_sequence_multihit(char* sequence, int sequenceLength, short** hash_table, short** next_array):
    """
    Builds a hash table of the byte sequence, allowing multiple hits
    """
    hash_table[0] = <short*>(calloc(hash_size, sizeof(short)))
    next_array[0] = <short*>(calloc(hash_size, sizeof(short)))

    if sequenceLength < hash_nucs:
        return

    if sequenceLength >= max_sequence_length:
        logger.warning("Trying to hash sequence that is too long. Length = %s. Limit is %s. Something is wrong here." %(sequenceLength, max_sequence_length))

    cdef int next_empty = 1
    cdef int hidx
    cdef int i
    cdef int j

    for i in range(sequenceLength - hash_nucs):
        hidx = my_hash(sequence+i)

        # add position i into hash table under hash hidx
        #*position[next_empty] = i
        # invariant: next_empty == i+1; therefore the position array is not necessary
        if hash_table[0][hidx] == 0:
            hash_table[0][hidx] = next_empty
        else:
            # find the end of the list
            j = hash_table[0][hidx]
            while next_array[0][j] != 0:
                j = next_array[0][j]
            next_array[0][j] = next_empty
        # point to next empty location
        next_empty += 1

###################################################################################################

cdef short* hash_sequence(char* sequence, int sequenceLength):
    """
    Returns a hash of the byte sequence
    """
    cdef char* seq = sequence
    cdef short* h = <short*>(calloc(hash_size, sizeof(short)))
    cdef int i
    cdef unsigned int hidx

    if sequenceLength < hash_nucs:
        return h

    if sequenceLength >= max_sequence_length:
        logger.warning("Trying to hash sequence that is too long. Length = %s. Limit is %s. Something is wrong here." %(sequenceLength, max_sequence_length))

    for i in range(sequenceLength - hash_nucs):

        #hidx = my_hash(seq+i)

        if i == 0:
            hidx = my_hash(seq + i)
        else:
            hidx = my_hash_update(seq[i+hash_nucs-1], hidx)

        if h[hidx] != 0:
            h[hidx] = -1
        else:
            h[hidx] = i+1

    return h

###################################################################################################

cdef void hashReadForMapping(cAlignedRead* read):
    """
    Hash the read, and store the hash for use in the mapReadToHaplotype
    function.
    """
    cdef int i
    read.hash = <short*>(calloc(read.rlen, sizeof(short)))

    for i in range(read.rlen - hash_nucs):

        if i == 0:
            read.hash[i] = my_hash(read.seq + i)
        else:
            read.hash[i] = my_hash_update(read.seq[i+hash_nucs-1], read.hash[i-1])

###################################################################################################

cdef int my_hash_update(char character, int oldVal):
    """
    Incrementally update hash value with one new character
    """
    global mask
    cdef int c = character & 7   # a,A->1  c,C->3  g,G->7  t,T->4

    if c == 7:
        c = 2

    return ( (oldVal << 2) & mask) + <unsigned int>(c & 3)

###################################################################################################

cdef int mapAndAlignReadToHaplotype(char* read, char* quals, int readStart, int hapStart, int readLen, int hapLen, short* haplotypeHash, short* haplotypeNextArray, short* readHash, char* haplotype, int gapExtend, int nucprior, short* localGapOpen, int* mapCounts, int mapCountsLen):
    """
    Map a single read to a small-ish (~1kb) sequence. This is used to find the best anchor point
    in a haplotype sequence for the specified read.
    Returns the read's index into the sequence.  May be negative.  Assumes a forward read direction
    """
    assert hapLen < max_sequence_length
    assert haplotypeHash != NULL
    assert haplotypeNextArray != NULL
    assert readHash != NULL
    assert readLen >= hash_nucs

    cdef int lenOfHapSeqToTest = readLen + 15
    cdef int maxcount = 0
    cdef int maxpos = 0
    cdef int i = -1
    cdef int pos = -1
    cdef int hapIdx = 0
    cdef int bestMappingPosition = -1
    cdef int count = -1
    cdef int hapLenForAlignment = 0
    cdef int indexOfReadIntoHap = -1
    cdef int readStartInHap = -1
    cdef int bestScore = 1000000
    cdef int alignScore = 1000000

    # If we have an exact match in the original position then simply return 0 now, as we can never
    # have a better match than that.
    indexOfReadIntoHap = readStart - hapStart

    if strncmp(read, haplotype + indexOfReadIntoHap, readLen) == 0:
        return 0

    # Reset counts
    memset(mapCounts, 0, sizeof(int)*mapCountsLen)

    # Do the mapping.
    for i in range(readLen - hash_nucs):
        hapIdx = haplotypeHash[ readHash[i] ]

        while hapIdx != 0:

            pos = hapIdx-1
            pos -= i     # account for index into read

            if pos >= hapLen:
                logger.warning("pos = %s. sl = %s. readHash[i] = %s" %(pos, hapLen, readHash[i]))

            count = mapCounts[pos + readLen] + 1
            mapCounts[pos + readLen] = count

            if count > maxcount:
                maxcount = count
                maxpos = pos

            # move to next hash hit
            hapIdx = haplotypeNextArray[ hapIdx ]

    if maxcount > 0:
        for i in range(hapLen + readLen):

            # Found match with best score (there may be many of these).
            if mapCounts[i] == maxcount:
                indexOfReadIntoHap = i - readLen

                # Align read around this position. Make sure we can't go off the end.
                if indexOfReadIntoHap >= -1*readLen and (indexOfReadIntoHap + readLen + 15 < hapLen):
                    readStartInHap = max(0,indexOfReadIntoHap-8)
                    hapLenForAlignment = readLen + 15 # This is fixed by the alignment algorithm
                    alignScore = fastAlignmentRoutine(haplotype + readStartInHap, read, quals, hapLenForAlignment, readLen, gapExtend, nucprior, localGapOpen + readStartInHap)

                    if alignScore < bestScore:
                        bestScore = alignScore
                        bestMappingPosition = indexOfReadIntoHap

                        # Short-circuit this loop if we find an exact match
                        if bestScore == 0:
                            return bestScore

    # Now try original mapping position. If the read is past the end of the haplotype then
    # don't allow the algorithm to align off the end of the haplotype. This will only happen in the case
    # of very large deletions.
    indexOfReadIntoHap = min(readStart - hapStart, hapLen - readLen - 15)

    # Only try this if the mapper hasn't already found this position
    if indexOfReadIntoHap != bestMappingPosition:
        readStartInHap = max(0,indexOfReadIntoHap-8)
        alignScore = fastAlignmentRoutine(haplotype + readStartInHap, read, quals, readLen+15, readLen, gapExtend, nucprior, localGapOpen + readStartInHap)

        if alignScore < bestScore:
            bestScore = alignScore
            bestMappingPosition = indexOfReadIntoHap

    if bestScore < 0:
        logger.error("Alignment score is %s" %(bestScore))

    return bestScore

###################################################################################################
