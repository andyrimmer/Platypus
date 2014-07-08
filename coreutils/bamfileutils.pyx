"""
Utilities for handling multiple bam files.
"""

from __future__ import division

import logging
import os
import cython

cimport samtoolsWrapper
cimport fastafile
cimport cython

from heapq import heappush,heappop
from fastafile cimport FastaFile
from samtoolsWrapper cimport AlignedRead,Samfile

logger = logging.getLogger("Log")

###################################################################################################

cdef class BamFileIterator(object):
    """
    Wrapper class for Bam iterator. This class implements a comparison function,
    so that instances may be stored in a sorted container. As well as the iterator,
    a reference to the current value is stored.
    """
    def __init__(self, theIterator):
        """
        Constructor. Takes a reference to the iterator.
        """
        self.iterator = theIterator
        self.currentValue = None

    def __cmp__(self, BamFileIterator other):
        """
        Comparison function. Utilises the comparison function defined in
        the AlignedRead class.
        """
        if other.iterator == None or other.currentValue:
            return 1

        return cmp(self.currentValue, other.currentValue)

    cpdef next(self):
        cdef AlignedRead theRead = self.iterator.next()
        self.currentValue = theRead

###################################################################################################

cdef class MultiBamFileIterator(object):
    """
    Class to implement functionality for iterating over multiple
    bam files simultaneously.
    """
    def __init__(self, list files, char* region, int start, int end):
        """
        Constructor.
        """
        cdef BamFileIterator generator
        cdef Samfile file

        self.files = files
        self.queue = []
        self.nReadsInQueue = 0

        for file in self.files:
            generator = BamFileIterator(file.fetch(region, start, end))

            try:
                generator.next()
                heappush(self.queue, generator)
                self.nReadsInQueue += 1
            except StopIteration:
                continue

        self.__iter__()

    def __next__(self):
        """
        Increment the iterator and yield the new value. Also, store the
        current value for use in the comparison function.
        """
        if len(self.queue) == 0:
            raise StopIteration

        cdef AlignedRead read = None
        cdef AlignedRead nextRead = None

        while self.nReadsInQueue > 0 and not read:

            nextIterator = heappop(self.queue)
            self.nReadsInQueue -= 1
            nextRead = nextIterator.currentValue

            if not nextRead.is_unmapped():

                read = nextRead
                self.lastRead = nextRead

            try:
                nextIterator.next()
                heappush(self.queue, nextIterator)
                self.nReadsInQueue += 1
            except StopIteration:
                continue

        if read:
            return read
        else:
            raise StopIteration

    def __iter__(self):
        """
        """
        return self

###################################################################################################

cdef class MultiBamFileReader(object):
    """
    Utility class for reading multiple Bam files. A list of Bam file names is given, in the
    constructor, and reads are read, in order, through a generator, using a priority queue
    to sort the reads as they come in.

    Obviously, the Bam files must be individually sorted in the first place, or this will not
    work.
    """
    def __init__(self, fileNames):
        """
        Constructor. Takes a list of Bam file names
        """
        self.files = []

        for f in fileNames:
            try:
                self.files.append(samtoolsWrapper.Samfile(f, 'rb'))
            except Exception, e:
                logger.warning("Could not open file %s. Data from this file will not be used" %(f))
                continue

        if len(self.files) == 0:
            raise StandardError, "MultiBamFileReader Could not open any of the required BAM files"

    def __del__(self):
        """
        Make sure to clean up all the Bam files when this instance is
        deleted.
        """
        for f in self.files:
            f.close()

    cpdef getRName(self, int refId):
        """
        Returns the string reference name corresponding to the specified refId.
        """
        return self.files[0].getrname(refId)

    cpdef MultiBamFileIterator fetch(self, char* region, int start, int end):
        """
        A generator function. Returns reads in the specified region using pysam.
        """
        return MultiBamFileIterator(self.files, region, start, end)

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

        theFile = open(options.regions[0], 'r')

        for line in theFile:
            chr,region = line.split(":")
            start = int(region.split("-")[0])
            end = int(region.split("-")[1])

            regions.append( (chr,start,end) )

        return regions

    elif options.regions == None:

        try:
            header = file.header
            regions = [ (d['SN'], 1, d['LN']) for d in header['SQ'] ]
        except:

            for region,regionTuple in refFile.refs.iteritems():
                regions.append((region, 1, regionTuple.SeqLength))

    elif os.path.isfile(options.regions[0]): # Could be a file...
        regionFile = open( options.regions[0] )
        regions = regionFile.readlines()
        regionFile.close()
    else:
        for region in options.regions:

            split = region.split(":")
            chr = bytes(split[0])

            if len( split ) == 2 :
                [ start, end ] = split[1].split("-")
                regions.append((chr,int(start),int(end)))
            elif len(split) == 1:
                start = 1
                try:
                    header = file.header
                    pivot = dict(zip([d['SN'] for d in header['SQ']], [d['LN'] for d in header['SQ']]))
                    end = pivot[chr]
                    regions.append((chr,int(start),int(end)))
                except:
                    regions = []
                    for region,regionTuple in refFile.refs.iteritems():
                        if region == chr:
                            regions.append((region, 1, regionTuple.SeqLength))
            else:
                regions.append((chr,None,None))

    if len(regions) == 0:
        logger.error("Platypus found no regions to search. Check that you are using the correct reference FASTA file or have specified the 'regions' argument correctly")
    elif len(regions) < 1000:
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

def cigarJimKent(AlignedRead read, int g_start):
    """
    Convert cigar string into Jim Kent's format.
    """
    cdef char* cigarChars = "MID"
    cdef list read_blockstart = []
    cdef list genome_blockstart = []
    cdef list blocklens = []
    cdef list inslens = []
    cdef int r = 0
    cdef int g = g_start
    cdef int inblock = True
    cdef int cigarOp = 0
    cdef int cigarOpLen = 0
    cdef int index = 0
    cdef int cigarLen = read.getCigarLength()

    for index from 0 <= index < cigarLen:

        cigarOp = cigarChars[read.getCigarOpCode(index)]
        cigarOpLen = read.getCigarOpLength(index)

        if cigarOp == 'M':

            if not inblock:
                inslens.append(0)

            read_blockstart.append(r)
            genome_blockstart.append(g)
            blocklens.append(cigarOpLen)
            r += cigarOpLen
            g += cigarOpLen
            inblock = False

        else:

            if inblock:
                read_blockstart.append(r)
                genome_blockstart.append(g)
                blocklens.append(0)
            if cigarOp == 'I' or cigarOp == 'S':
                inslens.append(cigarOpLen)
                r += cigarOpLen
            elif cigarOp == 'D' or cigarOp == 'N':
                inslens.append(-cigarOpLen)
                g += cigarOpLen

            inblock = True

    if inblock:
        read_blockstart.append(r)
        genome_blockstart.append(g)
        blocklens.append(0)

    return read_blockstart, genome_blockstart, blocklens, inslens

###################################################################################################

def getReadInfo(Samfile bamFile, AlignedRead read):
    """
    Dumb utility function so that we can retrieve read data, without re-writing all
    the AlignedRead functions.
    """
    # obtain read's read group, or None if not present
    cdef dict tags = read.tags()

    cdef int flag = read.flag()
    cdef int tid = read.rname()
    cdef int pos = read.pos()
    cdef int isize = read.isize()

    cdef bytes seq = read.seq()
    cdef bytes qual = read.qual()
    cdef bytes label = read.fastQName()
    cdef bytes thisReadGroup = None

    try:
        thisReadGroup = tags['RG']
    except:
        pass

    cdef bytes chrom

    if tid != -1:
        chrom = bamFile.getrname(tid)
    else:
        chrom = bytes("*")

    return (thisReadGroup,chrom,flag,pos,isize,seq,qual,label)

###################################################################################################
