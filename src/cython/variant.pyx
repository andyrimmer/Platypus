"""
Utility module, containing classes to handle reading and writing
files of variants.
"""

import math
import datetime
import logging
import cython
import samtoolsWrapper
import fastafile

cimport cython
cimport fastafile
cimport samtoolsWrapper
cimport cerrormodel

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
from samtoolsWrapper cimport Read_IsCompressed
from samtoolsWrapper cimport compressRead
from samtoolsWrapper cimport uncompressRead

from fastafile cimport FastaFile

logger = logging.getLogger("Log")

###################################################################################################

cdef double mLTOT = -0.23025850929940459    # Minus log ten over ten

###################################################################################################

cdef int PLATYPUS_VAR = 1
cdef int FILE_VAR = 2
cdef int ASSEMBLER_VAR = 4

###################################################################################################

cdef int SNP = 0
cdef int MNP = 1
cdef int INS = 2
cdef int DEL = 3
cdef int REP = 4

cdef list varTypes = ["SNP", "MNP", "INS", "DEL", "REP"]

###################################################################################################

cdef extern from "math.h":
    double exp(double)
    double log(double)
    double log10(double)
    double fabs(double)
    double sqrt(double)

###################################################################################################

cdef dict indel_prior_model = {1: "LIGC@:62/-*'&%$", 
                               2: "LIGDB@><9630.,+**)(''&&%%%$$$", 
                               3: "LIGA@B@><;8763220/.-,+++)*))(((''''&&&&&&%%%%%%%%$$$$$$$", 
                               4: "LIGA@???=<886533210/.--,+**))))((('''''&&&&&&&&%%%%%%%%%%%$$$$$$$$", 
                               5: 'LIGA@??>=>=;966543210///-,,++*', 
                               6: 'LIGA@??>>=<=;:764532210/----,++', 
                               7: 'LIGA@??>>==<;;987543210/....-,,,++++', 
                               8: 'LIGA@??>>==<<;9876432200/..--,,,+++', 
                               9: 'LIGA@??>>==<<;;9966432100//../..----,,,,,++++++', 
                               10: 'LIGA@??>>==<<;;:986432110//..----,,,,++++', 
                               11: 'LIGA@??>>==<<<;;:87642210////..--,,,,,+++', 
                               12: 'LIGA@??>>==<<<;;;:986532110000/...-----,,,,,+++++', 
                               13: 'LIGA@??>>==<<<;;;::987543111000/////.......--------,,,,,,,,,,,,,+++++++++', 
                               14: 'LIGA@??>>==<<<;;;::987642210/0/.....-------,,,,,,,,+++++++', 
                               15: 'LIGA@??>>==<<<;;;;::988754322110000////////.......------------,,,,,,,,,,,,,,,,,++++++++++', 
                               16: 'LIGA@??>>==<<<;;;;:::98765321110////........-------,,,,,,,,,,,,,,+++++++++', 
                               17: 'LIGA@??>>==<<<;;;;::::988764433211110000000///////.............-----------------,,,,,,,,,,,,,,,,,,,', 
                               18: 'LIGA@??>>==<<<;;;:::::998875433221111000000///////.............-----------------,,,,,,,,,,,,,,,,,,,', 
                               19: 'LIGA@??>>==<<<;;;;::::999887654433222221111111100000000//////////////..................------------', 
                               20: 'LIGA@??>>==<<<;;;;::::9999876543322111000000///////............-----------------,,,,,,,,,,,,,,,,,,,', 
                               21: 'LIGA@??>>==<<<;;;;::::9999988765544433322222221111111100000000000000//////////////////.............', 
                               22: 'LIGA@??>>==<<<;;;;::::9999987765432221000000////////...........-----------------,,,,,,,,,,,,,,,,,,,', 
                               23: 'LIGA@??>>==<<<;;;;::::9999998776543322111100000000////////................-------------------,,,,,,', 
                               24: 'LIGA@??>>==<<<;;;;::::9999998887654433322111111100000000/////////////...................-----------'}

# Currently hard-coded insertion/deletion priors
cdef double complex_deletion_prior = 5e-5
cdef double complex_insertion_prior = 5e-6

###################################################################################################

@cython.final
@cython.freelist(1000)
cdef class Variant(object):
    """
    Class to encapsulate information for all common variant types. The basic
    idea is to tread all variants as replacements, i.e. a certain sequence of
    bases is removed from the reference, and a certain sequence is added. A SNP is
    a replacement of 1 base with 1 different base, a deletion is a replacement of one
    sequence with a smaller sequence (or none at all), and an insertion is a replacement
    of one sequence with a longer sequence.
    """
    def __init__(self, bytes refName, int refPos, char* removed, char* added, int nSupportingReads, int varSource):
        """
        Constructor. Takes the name of the reference sequence (e.g. 'chr1'),
        the position of the variant (given by the left-most co-ordinate of the
        changed sequence (so a SNP co-ordinate is the exact location of the SNP,
        an indel co-ordinate is the location of the first inserted base, and a deletion
        co-ordinate is the location of the first deleted base)), the removed and added
        sequences, and the supporting read.
        """
        # Make sure we don't go negative here.
        refPos = max(0, refPos)
        self.refName = refName
        self.nAdded = len(added)
        self.nRemoved = len(removed)
        self.varSource = varSource

        # SNP
        if self.nRemoved == 1 and self.nAdded == 1:
            self.minRefPos = refPos
            self.maxRefPos = refPos
            self.varType = SNP

        # Multi-nucleotide substitution
        elif self.nRemoved == self.nAdded:
            self.minRefPos = refPos
            self.maxRefPos = refPos + self.nAdded - 1
            self.varType = MNP

        # Arbitrary sequence replacement
        else:
            self.minRefPos = refPos
            self.maxRefPos = refPos + self.nRemoved

            if self.nRemoved == 0:
                self.varType = INS
            elif self.nAdded == 0:
                self.varType = DEL
            else:
                self.varType = REP

        self.refPos = refPos
        self.bamMinPos = refPos
        self.bamMaxPos = refPos
        self.removed = removed
        self.added = added
        self.bamAdded = added
        self.bamRemoved = removed
        self.nSupportingReads = nSupportingReads
        self.hashValue = -1

    cdef double indelPrior(self, FastaFile refFile, int indel_length_and_type):
        """
        Calculate indel prior, based on sequence context
        """
        context = 100
        leftPos = max(0, self.refPos - context)
        rightPos = self.refPos + context          # no range check here

        cdef int relRefPos = self.refPos - leftPos
        cdef bytes sequence
        cdef bytes sizes
        cdef bytes displacements
        cdef char qbase = 33
        cdef char prior = (<char*>indel_prior_model[1][0])[0] - qbase
        cdef char newprior
        cdef unsigned char size
        cdef unsigned char prior_tractlength = -1
        cdef object disp
        cdef double dprior

        try:
            sequence = refFile.getSequence( self.refName, leftPos+1, rightPos+1 )
        except IndexError:
            sequence = <bytes>""

        (sizes, displacements) = cerrormodel.calculate_size_and_displacement(sequence, True)
        #logger.info("Sequence = %s" %(sequence))
        #logger.info("Sizes = %s. Displacements = %s" %(str(sizes), str(displacements)))

        for i in range(relRefPos-1, relRefPos+1):
            disp = (<char*>displacements)[i]
            #logger.info("Disp = %s" %(disp))

            if disp in indel_prior_model:
                size = (<unsigned char*>sizes)[i]
                #logger.info("Size = %s" %(size))

                if size > len(indel_prior_model[disp]):
                    size =  len(indel_prior_model[disp])
                    #logger.info("Size2 = %s" %(size))

                # get prior -- the -1 is to account for the model starting with size=1, not 0
                newprior = (<char*>(indel_prior_model[ disp ]))[ size-1 ] - qbase
                #logger.info("new prior = %s" %(newprior))

                if newprior < prior:
                    prior = newprior
                    prior_tractlength = size

        # convert into a double
        dprior = math.pow(0.1, prior/10.0)
        #logger.info("dprior = %s" %(dprior))

        # Take into account length for indels in non-repetitive contexts (tract length = 1,2,3)
        # use geometric distribution with exponent 0.75; normalize by factor 1.0-0.75 to sum to 1
        # this ensures that long insertions/deletions in complex sequence get a low prior probability
        if prior_tractlength <= 3:

            # Old indel prior
            #dprior = math.pow(0.75, abs(indel_length_and_type)-1) * (1.0 - 0.75)

            # New indel prior
            if indel_length_and_type < 0:
                dprior = complex_deletion_prior * math.pow(0.75, (-indel_length_and_type)-1) * (1.0 - 0.75)
                #logger.info("Here 1. dprior = %s" %(dprior))
            else:
                dprior = complex_insertion_prior * math.pow(0.75, indel_length_and_type-1) * (1.0 - 0.75) * math.pow(0.33, indel_length_and_type)
                #logger.info("Here 2. dprior = %s. length and type = %s." %(dprior,indel_length_and_type))
                #logger.info("cip = %s. pow1 = %s. pow2 = %s" %(complex_insertion_prior, math.pow(0.75, indel_length_and_type-1), math.pow(0.33, indel_length_and_type)))

        #logger.debug("Prior for Indel %s is %s. Context is %s." %(self, dprior, refFile.getSequence(self.refName, self.refPos-10, self.refPos+10)))
        return dprior

    cdef double calculatePrior(self, FastaFile refFile):
        """
        Calculate and return the prior probability for this
        variant.
        """
        cdef double prior = 0.0

        # Basic prior for SNPs is 1e-3, and for indels is 1e-4
        if self.nAdded == 1 and self.nRemoved == 1:
            # Use uniform prior for all 3 possibilities
            prior = 1e-3 / 3

        # Multi-nucleotide substitution
        elif self.nAdded == self.nRemoved:
            # Made-up prior
            nDiffs = len([(x,y) for (x,y) in zip(self.added,self.removed) if x != y])
            prior = 5e-5*(0.1**(nDiffs-1)) * (1.0 - 0.1)

        # Insertion
        elif self.nAdded > 0 and self.nRemoved == 0:
            # most indels are slippage events, so approximate using a dirac prior
            # on the sequence (and don't check correctness of the sequence!)
            prior = self.indelPrior( refFile, self.nAdded )
            #logger.info("Prior for variant %s is %s. Old prior = %s." %(self, prior, 1e-4*0.25**(self.nAdded)))
            #prior = 1e-4*0.25**self.nRemoved

        # Deletion
        elif self.nAdded == 0 and self.nRemoved > 0:
            # most indels are slippage events, so approximate using a dirac prior
            # on the sequence (and don't check correctness of the sequence!)
            prior = self.indelPrior( refFile, -self.nRemoved )
            #logger.info("Prior for variant %s is %s. Old prior = %s." %(self, prior, 1e-4*0.6**self.nRemoved))
            #prior = 1e-4*0.6**self.nRemoved

        # Replacement -- made-up prior for now
        else:
            prior = 5e-6

        # Cap the prior at some sensible value, otherwise large variants get ridiculous
        # priors. We need better priors for ALU events and similar large indels.
        return max(prior, 1e-10)

    cdef void addVariant(self, Variant other):
        """
        Add supporting data from another variant instance.
        """
        self.nSupportingReads += other.nSupportingReads
        self.varSource |= other.varSource
        self.bamMinPos = min(self.bamMinPos, other.bamMinPos)
        self.bamMaxPos = max(self.bamMaxPos, other.bamMaxPos)

    def __hash__(self):
        """
        Implementing this function allows variants to be hashed, and so stored in
        a set or dictionary. The supporting reads are not included in the hashing, as
        we want two variants to give the same hash id if they have the same position
        and added/removed sequences.
        """
        if self.hashValue == -1:
            self.hashValue = hash( (self.refName, self.refPos, self.removed, self.added) )

        return self.hashValue

    def __richcmp__(Variant self, Variant other, int opCode):
        """
        Comparison function:

        Are two variants equal? Only return true if the positions and added and
        removed sequences are exactly equal.
        """
        cdef int thisRefPos = self.refPos
        cdef int otherRefPos = other.refPos
        cdef int thisType = self.varType
        cdef int otherType = other.varType
        cdef int thisNAdded = self.nAdded
        cdef int otherNAdded = other.nAdded
        cdef int thisNRemoved = self.nRemoved
        cdef int otherNRemoved = other.nRemoved
        cdef bytes thisRefName = self.refName
        cdef bytes otherRefName = other.refName
        cdef bytes thisAdded = self.added
        cdef bytes otherAdded = other.added
        cdef bytes thisRemoved = self.removed
        cdef bytes otherRemoved = other.removed

        # <
        if opCode == 0:
            if thisRefName < otherRefName:
                return True
            elif thisRefName == otherRefName and thisRefPos < otherRefPos:
                return True
            elif thisRefName == otherRefName and thisRefPos == otherRefPos and thisType < otherType:
                return True
            elif thisRefName == otherRefName and thisRefPos == otherRefPos and thisType == otherType and thisNRemoved < otherNRemoved:
                return True
            else:
                return False
        # <=
        elif opCode == 1:
            if thisRefName > otherRefName:
                return False
            elif thisRefName == otherRefName and thisRefPos > otherRefPos:
                return False
            elif thisRefName == otherRefName and thisRefPos == otherRefPos and thisType > otherType:
                return False
            elif thisRefName == otherRefName and thisRefPos == otherRefPos and thisType == otherType and thisNRemoved > otherNRemoved:
                return False
            else:
                return True
        # ==
        elif opCode == 2:
            if thisRefName == otherRefName and thisRefPos == otherRefPos and thisAdded == otherAdded and thisRemoved == otherRemoved:
                return True
            else:
                return False
        # !=
        elif opCode == 3:
            if thisRefName != otherRefName or thisRefPos != otherRefPos or thisAdded != otherAdded or thisRemoved != otherRemoved:
                return True
            else:
                return False
        # >
        elif opCode == 4:
            if thisRefName > otherRefName:
                return True
            elif thisRefName == otherRefName and thisRefPos > otherRefPos:
                return True
            elif thisRefName == otherRefName and thisRefPos == otherRefPos and thisType > otherType:
                return True
            elif thisRefName == otherRefName and thisRefPos == otherRefPos and thisType == otherType and thisNRemoved > otherNRemoved:
                return True
            else:
                return False
        # >=
        elif opCode == 5:
            if thisRefName < otherRefName:
                return False
            elif thisRefName == otherRefName and thisRefPos < otherRefPos:
                return False
            elif thisRefName == otherRefName and thisRefPos == otherRefPos and thisType < otherType:
                return False
            elif thisRefName == otherRefName and thisRefPos == otherRefPos and thisType == otherType and thisNRemoved < otherNRemoved:
                return False
            else:
                return True

    def __str__(self):
        """
        Another way of printing the variant.
        """
        string = "%s(%s:%s-%s" %(varTypes[self.varType], self.refName,self.minRefPos,self.maxRefPos)

        if self.nRemoved > 0:
            string += (' -' + self.removed)
        if self.nAdded > 0:
            string += (' +' + self.added)

        string += " nReads = %s, Source= %s" %(self.nSupportingReads, self.varSource)
        string += ")"
        return string

    def __repr__(self):
        """
        __repr__ is called when you do "print variant", and will return a short string
        describing the variant.
        """
        return str(self)

    def shortRepr(self):
        """
        A more succinct way of printing the variant.
        """
        string = "%s(%s:%s-%s" %(varTypes[self.varType], self.refName,self.minRefPos,self.maxRefPos)

        if self.nRemoved > 0:
            string += (' -' + self.removed)
        if self.nAdded > 0:
            string += (' +' + self.added)

        string += ")"
        return string

    cdef int overlaps(self, Variant other):
        """
        """
        if other.minRefPos < self.minRefPos < other.maxRefPos:
            return True

        elif self.minRefPos < other.minRefPos < self.maxRefPos:
            return True

        elif self.minRefPos == other.minRefPos:

            # 2 SNPS overlap in this case. MNPs covered by previous check.
            if self.nAdded == self.nRemoved and other.nAdded == other.nRemoved:
                return True

            # 1 SNP and 1 indel, at the same position. No overlap.
            elif (self.varType == SNP and other.nAdded != other.nRemoved) or (other.varType == SNP and self.nAdded != self.nRemoved):
                return False

            # 2 Indels or one MNP and one Indel. Overlap
            else:
                return True

        elif self.minRefPos == other.maxRefPos:

            # 2 SNPS/MNPs overlap in this case
            if self.nAdded == self.nRemoved and other.nAdded == other.nRemoved:
                return True

            # 1 SNP/MNP and 1 indel, at the same position. No overlap if SNP/MNP is first
            elif (self.nAdded != self.nRemoved and other.nAdded == other.nRemoved):
                return False

            # 2 Indels or Indel followed by MNP/SNP. Overlap
            else:
                return True

        elif self.maxRefPos == other.minRefPos:

            # 2 SNPS/MNPs overlap in this case
            if self.nAdded == self.nRemoved and other.nAdded == other.nRemoved:
                return True

            # 1 SNP/MNP and 1 indel, at the same position. No overlap if SNP/MNP is first
            elif (other.nAdded != other.nRemoved and self.nAdded == self.nRemoved):
                return False

            # 2 Indels or Indel followed by MNP/SNP. Overlap
            else:
                return True

        # No overlap.
        else:
            return False

###################################################################################################

@cython.final
cdef class VariantCandidateGenerator(object):
    """
    A class to generate variant candidates from a bunch of reads.
    """
    def __init__(self, tuple region, FastaFile referenceFile, int minMapQual, int minFlank, int minBaseQual, int maxCoverage, int maxReadLength, options,  int verbosity=2, int genSNPs=1, int genIndels=1):
        """
        Constructor. Takes read buffer, reader, reference, and a set of quality arguments for the variant
        candidates. Create a storage splace for variant candidates, and store the values of some flags which
        are used in the pysam CIGAR information (these should really be module level variables).
        """
        self.CIGAR_M = 0 # Match
        self.CIGAR_I = 1 # Insertion
        self.CIGAR_D = 2 # Deletion
        self.CIGAR_N = 3 # Skipped region from reference
        self.CIGAR_S = 4 # Soft clipping. Sequence is present in read
        self.CIGAR_H = 5 # Hard clipping. Sequence is not present in read
        self.CIGAR_P = 6 # Padding. Used for padded alignment

        self.minMapQual = minMapQual
        self.minBaseQual = minBaseQual
        self.minFlank = minFlank
        self.variantHeap = {} # List of variants
        self.refFile = referenceFile
        self.rname = region[0]
        self.refSeqStart = max(1, region[1]-2000) # Don't try to fetch anything < 0
        self.refSeqEnd = min(region[2]+2000, self.refFile.refs[region[0]].SeqLength-1) # Don't try to fetch anything > seqLength
        self.pyRefSeq = self.refFile.getSequence(self.rname, self.refSeqStart, self.refSeqEnd) # Cache this
        self.refSeq = self.pyRefSeq
        self.rStart = region[1]
        self.rEnd = region[2]
        self.maxCoverage = maxCoverage
        self.maxReadLength = maxReadLength
        self.verbosity = verbosity
        self.genSNPs = genSNPs
        self.genIndels = genIndels
        self.options = options
        self.qualBinSize = options.qualBinSize

    cdef void addVariantToList(self, Variant var):
        """
        Check if the variant is already on the heap; if so, increment the number of
        supporting reads; if not, add it.
        """
        cdef Variant theVar = self.variantHeap.get(var)

        if theVar is None:
            self.variantHeap[var] = var

            if self.verbosity >= 4:
                logger.debug("Adding new variant %s to candidate list" %(var))

        else:
            theVar.addVariant(var)

            if self.verbosity >= 4:
                logger.debug("Adding variant %s to existing variant in candidate list" %(var))

    cdef void getSnpCandidatesFromReadSegment(self, cAlignedRead* read, char* readSeq, char* readQual, int readStart, int readOffset, int refOffset, int lenSeqToCheck, int minFlank):
        """
        Get all SNP candidates from a particular read segement.

        Args:
        char* readSeq -- The complete sequence of bases in the read.
        char* readQual -- The complete sequence of quality scores in the read.
        int readStart -- The starting position, in the reference sequence, of the read.
        int readOffset -- If we're not starting from the beginning of the read, then how far along the read sequence to start.
        int refOffset -- If the read contains indels, we need to offset our reference position accordingly, by this much.
        int lenSeqToCheck -- The number of bases to check.
        """
        cdef int index = 0
        cdef int baseQual = 0
        cdef int readIndex = 0
        cdef int refIndex = 0
        cdef char readChar
        cdef char refChar

        cdef bytes mSNPRef
        cdef bytes mSNPRead

        cdef int misMatchStartRef = -1
        cdef int misMatchEndRef = -1
        cdef int misMatchStartRead = -1
        cdef int misMatchEndRead = -1

        for index from 0 <= index < lenSeqToCheck:

            # Don't look at the first 'minFlank' bases for candidate generation
            if readOffset == 0 and index < minFlank:
                continue

            # Don't look at the last 'minFlank' bases for candidate generation
            if index + readOffset >= read.rlen - minFlank:
                continue

            #if index + readOffset + refOffset >= read.rlen:
            #    logger.info("read len = %s. index = %s. read offset = %s. ref offset = %s" %(read.rlen, index, readOffset, refOffset))

            readIndex = index + readOffset
            refIndex = (index + refOffset + readStart) - self.refSeqStart
            readChar = readSeq[readIndex]
            refChar = self.refSeq[refIndex]
            baseQual = readQual[readIndex]

            assert baseQual >= 0, "Something is very wrong. Base qual is %s" %(baseQual)
            assert baseQual <= 93, "Something is very wrong. Base qual is %s" %(baseQual)

            if readChar != refChar:

                if readChar != 'N' and refChar != 'N' and baseQual >= self.minBaseQual:

                    if misMatchStartRef == -1:
                        misMatchStartRef = refIndex
                        misMatchEndRef = refIndex
                        misMatchStartRead = readIndex
                        misMatchEndRead = readIndex

                    elif refIndex - misMatchEndRef <= minFlank:
                        misMatchEndRef = refIndex
                        misMatchEndRead = readIndex

                    else:

                        if self.verbosity >= 3:
                            logger.debug("Splitting long mis-match into two variants at %s." %(misMatchStartRef+self.refSeqStart))

                        mSNPRef = self.refSeq[misMatchStartRef: misMatchEndRef+1]
                        mSNPRead = readSeq[misMatchStartRead: misMatchEndRead+1]
                        #logger.info("Adding MNP. Added = %s. removed = %s. pos = %s" %(mSNPRead, mSNPRef, misMatchStartRef + self.refSeqStart))
                        self.addVariantToList(Variant(self.rname, misMatchStartRef + self.refSeqStart, mSNPRef, mSNPRead, 1, PLATYPUS_VAR))

                        misMatchStartRef = refIndex
                        misMatchEndRef = refIndex
                        misMatchStartRead = readIndex
                        misMatchEndRead = readIndex

            else:

                # We have a mis-match sequence, and now we're past the end of it, as there are
                # > minFlank matches at the end.
                if misMatchStartRef != -1 and refIndex - misMatchEndRef > minFlank:

                    mSNPRef = self.refSeq[misMatchStartRef: misMatchEndRef+1]
                    mSNPRead = readSeq[misMatchStartRead: misMatchEndRead+1]
                    #logger.info("Adding MNP. Added = %s. removed = %s. pos = %s" %(mSNPRead, mSNPRef, misMatchStartRef + self.refSeqStart))
                    self.addVariantToList(Variant(self.rname, misMatchStartRef + self.refSeqStart, mSNPRef, mSNPRead, 1, PLATYPUS_VAR))

                    misMatchStartRef = -1
                    misMatchEndRef = -1
                    misMatchStartRead = -1
                    misMatchEndRead = -1

        # Catch the last one... do I need to do this?
        if misMatchStartRef != -1:
            mSNPRef = self.refSeq[misMatchStartRef: misMatchEndRef+1]
            mSNPRead = readSeq[misMatchStartRead: misMatchEndRead+1]
            #logger.info("Adding MNP. Added = %s. removed = %s. pos = %s" %(mSNPRead, mSNPRef, misMatchStartRef + self.refSeqStart))
            self.addVariantToList(Variant(self.rname, misMatchStartRef + self.refSeqStart, mSNPRef, mSNPRead, 1, PLATYPUS_VAR))

    cdef void getVariantCandidatesFromSingleRead(self, cAlignedRead* read):
        """
        Check a single read for variant candidates. Variants are flagged by the CIGAR string. Pysam reports
        the CIGAR string information as a list of tuples, where each tuple is a pair, and the first element
        gives the type of feature (match, insertion or deletion), and the second element gives the number of
        nucleotides associated. For example, [(0, 1), (1, 2), (0, 1)] is a 1 base match, a 2 base insertion,
        and a 1 base match.
        """
        cdef int readStart = read.pos
        cdef int readLength = read.rlen
        cdef int mapQ = read.mapq
        cdef int flag = 0
        cdef int length = 0
        cdef int refOffset = 0
        cdef int readOffset = 0
        cdef int cigarIndex = 0
        cdef int cigarLength = read.cigarLen
        cdef char* readQual = read.qual
        cdef char* readSeq = read.seq
        cdef bytes insertedSequence = None
        cdef bytes deletedSequence = None

        for cigarIndex from 0 <= cigarIndex < cigarLength:

            flag = read.cigarOps[2*cigarIndex]
            length = read.cigarOps[ (2*cigarIndex) + 1]

            # An insertion take us further along the read, but not the reference
            if flag == self.CIGAR_I:

                # Skip this insertion if it isn't flanked by a matching sequence >= minFlank
                if cigarIndex > 0 and read.cigarOps[ (2*cigarIndex) -2] == self.CIGAR_M and read.cigarOps[ (2*cigarIndex)-1] >= self.minFlank:
                    pass
                elif cigarIndex < cigarLength-1 and read.cigarOps[ (2*cigarIndex) +2] == self.CIGAR_M and read.cigarOps[ (2*cigarIndex) +3] >= self.minFlank:
                    pass
                else:
                    readOffset += length
                    continue

                insertedSequence = readSeq[readOffset : readOffset+length]

                if insertedSequence.count("N") == 0 and self.genIndels:
                    #logger.info("Adding insertion. Length = %s. Added = %s. removed = %s. pos = %s" %(length, insertedSequence, "", readStart + refOffset -1))
                    self.addVariantToList(Variant(self.rname, readStart+refOffset-1, "", insertedSequence, 1, PLATYPUS_VAR))

                readOffset += length

            # A deletion take us further along the reference, but not the read
            elif flag == self.CIGAR_D:

                # Skip this deletion if it isn't flanked by a matching sequence >= minFlank
                if cigarIndex > 0 and read.cigarOps[2*cigarIndex-2] == self.CIGAR_M and read.cigarOps[2*cigarIndex-1] >= self.minFlank:
                    pass
                elif cigarIndex < cigarLength-1 and read.cigarOps[2*cigarIndex+2] == self.CIGAR_M and read.cigarOps[2*cigarIndex+3] >= self.minFlank:
                    pass
                else:
                    refOffset += length
                    continue

                deletedSequence = self.refFile.getSequence(self.rname, readStart+refOffset, (readStart+refOffset+length))

                # Don't look at deletions with Ns in them
                if deletedSequence.count("N") == 0 and self.genIndels:
                    #logger.info("Adding deletion. Length = %s. Added = %s. removed = %s. pos = %s" %("", len(deletedSequence), deletedSequence, readStart + refOffset -1))
                    self.addVariantToList(Variant(self.rname, readStart+refOffset-1, deletedSequence, "", 1, PLATYPUS_VAR))

                refOffset += length

            # A match take us further along the reference and the read
            elif flag == self.CIGAR_M:

                # Don't generate SNP candidates from matching sequences < minFlank
                if length < self.minFlank:
                    readOffset += length
                    refOffset += length
                    continue

                if self.genSNPs:
                    self.getSnpCandidatesFromReadSegment(read, readSeq, readQual, readStart, readOffset, refOffset, length, self.minFlank)

                readOffset += length
                refOffset += length

            # Skipped region from the reference.
            elif flag == self.CIGAR_N:
                #readOffset += length
                refOffset += length

            # Soft clipping. Sequence is present in read, but we should ignore it.
            elif flag == self.CIGAR_S:
                readOffset += length

                # We move back read positions when there is soft clipping at the beginning of
                # reads
                if cigarIndex == 0:
                    refOffset += length

            # Hard clipping. Sequence is not present in read.
            elif flag == self.CIGAR_H:
                continue

            # Padding. We do nothing here.
            elif flag == self.CIGAR_P:
                continue

            # Other kinds of flag.
            else:
                continue

    cdef void addCandidatesFromReads(self, cAlignedRead** readStart, cAlignedRead** readEnd):
        """
        Loop through all reads, and flag candidate variants.
        """
        cdef int nReads = 0
        cdef int compressed = 0

        while readStart != readEnd:

            compressed = 0

            if Read_IsQCFail(readStart[0]):
                readStart += 1
                continue

            if Read_IsCompressed(readStart[0]):
                compressed = 1
                uncompressRead(readStart[0], self.refSeq, self.refSeqStart, self.refSeqEnd, self.qualBinSize)

            self.getVariantCandidatesFromSingleRead(readStart[0])

            if compressed:
                compressRead(readStart[0], self.refSeq, self.refSeqStart, self.refSeqEnd, self.qualBinSize, 1)

            nReads += 1
            readStart += 1

    cdef list getCandidates(self, int minReads):
        """
        Return list of candidates.
        """
        return [x for x in sorted(self.variantHeap.values())]

###################################################################################################
