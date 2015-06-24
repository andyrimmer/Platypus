"""
FastaFile is a utility class used for reading the Fasta file format,
and facilitating access to reference sequences.
"""

import cython
import logging
import os
cimport cython

###################################################################################################

cdef extern from "stdlib.h":
    long long int atoll(char*)

###################################################################################################

logger = logging.getLogger("Log")

###################################################################################################

cdef class sequenceTuple:
    """
    Structure for storing data from line of fasta index file.
    """
    cdef public bytes SeqName
    cdef public long long int SeqLength
    cdef public long long int StartPos
    cdef public long long int LineLength
    cdef public long long int FullLineLength

    def __init__(self, bytes seqName, long long int seqLength, long long int startPos, long long int lineLength, long long int fullLength):
        """
        Constructor
        """
        self.SeqName = seqName
        self.SeqLength = seqLength
        self.StartPos = startPos
        self.LineLength = lineLength
        self.FullLineLength = fullLength

###################################################################################################

def expandedOpen(path, mode):
    try:
        return open(path, mode)
    except IOError:
        return open(os.path.expanduser(path), mode)

###################################################################################################

cdef class FastaIndex:
    """
    Index file for a FastaFile. Contains start and end positions of all
    sequences in the fasta file.
    """
    def __init__(self, fileName, mode='rb'):
        """
        Constructor. Takes the name of an Index file. For now, the whole file
        will be read into memory and the references stored in a set.
        """
        self.theFile = expandedOpen(fileName, mode)

    cdef dict getRefs(self, int parseNCBI):
        """
        Read next line
        """
        cdef dict theDict = {}

        for theLine in self.theFile:

            line = theLine.strip().split("\t")
            seqName = line[0].split()[0]

            if seqName.startswith('gi|') and parseNCBI:       # NCBI-formatted line
                ids = seqName.split('|')
                if len(ids)>=4 and ids[2] == "ref":
                    seqName = ids[3]

            theDict[seqName] = sequenceTuple(line[0], atoll(line[1]), atoll(line[2]), atoll(line[3]), atoll(line[4]))

        return theDict

###################################################################################################

cdef class FastaFile:
    """
    Utility for reading sequence from Fasta files.
    """
    def __init__(self, fileName, indexName, mode='rb', parseNCBI=True):
        """
        Constructor. Takes file-name and index file-name
        """
        self.theFile = expandedOpen(fileName, mode)
        self.theIndex = FastaIndex(indexName)
        self.refs = self.theIndex.getRefs( parseNCBI )
        self.cacheRefName = None
        self.cacheStartPos = -1
        self.cacheEndPos = -1
        self.cache = None

    def close(self):
        """
        Wrapper function to close self.theFile
        """
        self.theFile.close()

    def getTotalSequenceLength(self):
        """
        Return the accumulated lengths of all sequences in the
        reference index.
        """
        length = 0

        for name, seqTuple in self.refs.iteritems():
            length += seqTuple.SeqLength

        return length

    cdef bytes getCharacter(self, bytes seqName, long long int pos):
        """
        Returns the character at the specified (0-indexed) position
        of the specified sequence.
        """
        cdef sequenceTuple seqTuple = self.refs[seqName]
        cdef long long int seqLength = seqTuple.SeqLength
        cdef long long int seqStartPos = seqTuple.StartPos
        cdef long long int LineLength = seqTuple.LineLength
        cdef long long int FullLineLength = seqTuple.FullLineLength

        if pos >= seqLength or pos < 0:
            return <char*>"-"

        self.theFile.seek(seqStartPos + pos + (FullLineLength - LineLength)*(<long long int>(<double>(pos)/(LineLength))))

        try:
            return self.theFile.read(1).upper()
        except Exception:
            return <char*>"-"

    cdef void setCacheSequence(self, bytes seqName, long long int beginPos, long long int endPos):
        """
        """
        if seqName not in self.refs:
            raise StandardError, "Invalid contig name %s. Make sure your FASTA reference file and query regions have the same naming convention" %(seqName)

        cdef sequenceTuple seqTuple = self.refs[seqName]
        cdef long long int seqLength = seqTuple.SeqLength

        beginPos = max(0, beginPos)
        endPos = min(seqLength-1, endPos)

        cdef long long int seqStartPos = seqTuple.StartPos
        cdef long long int LineLength = seqTuple.LineLength
        cdef long long int FullLineLength = seqTuple.FullLineLength
        cdef long long int desiredSequenceStartPos = seqStartPos + beginPos + (FullLineLength - LineLength)*<long long int>(<double>beginPos/LineLength)
        cdef long long int desiredSequenceEndPos = seqStartPos + endPos + (FullLineLength - LineLength)*<long long int>(<double>endPos/LineLength)
        cdef long long int desiredSeqLength = (endPos - beginPos)
        cdef long long int desiredSequenceLengthInFile = (desiredSequenceEndPos - desiredSequenceStartPos)

        if endPos < beginPos:
            raise IndexError, "Cannot have beginPos = %s, endPos = %s" %(beginPos, endPos)

        if endPos > seqLength or beginPos < 0:
            raise IndexError, "Cannot return sequence from %s to %s. Ref seq length = %s" %(beginPos, endPos, seqLength)

        self.theFile.seek(desiredSequenceStartPos)
        self.cache = self.theFile.read(desiredSequenceLengthInFile).replace("\n", "").upper()
        self.cacheRefName = seqName
        self.cacheStartPos = beginPos
        self.cacheEndPos = endPos

    cdef bytes getSequence(self, bytes seqName, long long int beginPos, long long int endPos):
        """
        Returns the character sequence between the the specified (0-indexed) start
        and end positions. This returns a half-open sequence interval, i.e. the returned
        sequence includes the character at beginPos, but not the one at endPos. This is done
        in order to make down-stream sequence handling easier.
        """
        if self.cache is not None:
            if beginPos >= self.cacheStartPos and endPos < self.cacheEndPos:
                #logger.debug("Getting %s:%s-%s from cache. cache index = %s:%s" %(seqName, beginPos, endPos, beginPos - self.cacheStartPos, endPos -self.cacheStartPos))
                return self.cache[beginPos - self.cacheStartPos: endPos - self.cacheStartPos]

        cdef sequenceTuple seqTuple = self.refs[seqName]
        cdef long long int seqLength = seqTuple.SeqLength

        beginPos = max(0, beginPos)
        endPos = min(seqLength-1, endPos)

        cdef long long int seqStartPos = seqTuple.StartPos
        cdef long long int LineLength = seqTuple.LineLength
        cdef long long int FullLineLength = seqTuple.FullLineLength
        cdef long long int desiredSequenceStartPos = seqStartPos + beginPos + (FullLineLength - LineLength)*<long long int>(<double>beginPos/LineLength)
        cdef long long int desiredSequenceEndPos = seqStartPos + endPos + (FullLineLength - LineLength)*<long long int>(<double>endPos/LineLength)
        cdef long long int desiredSeqLength = (endPos - beginPos)
        cdef long long int desiredSequenceLengthInFile = (desiredSequenceEndPos - desiredSequenceStartPos)

        if endPos < beginPos:
            raise IndexError, "Cannot have beginPos = %s, endPos = %s" %(beginPos, endPos)

        if endPos > seqLength or beginPos < 0:
            raise IndexError, "Cannot return sequence from %s to %s. Ref seq length = %s" %(beginPos, endPos, seqLength)

        self.theFile.seek(desiredSequenceStartPos)
        cdef bytes seq = self.theFile.read(desiredSequenceLengthInFile)
        return seq.replace("\n", "").upper()

###################################################################################################
