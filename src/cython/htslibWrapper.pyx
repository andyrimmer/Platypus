import os
import types
import itertools
import logging
cimport cython
import zlib

###################################################################################################

logger = logging.getLogger("Log")

###################################################################################################

# defines imported from samtools
DEF SEEK_SET = 0
DEF SEEK_CUR = 1
DEF SEEK_END = 2

## These are bits set in the flag.
## have to put these definitions here, in samtoolsWrapper.pxd they got ignored
## @abstract the read is paired in sequencing, no matter whether it is mapped in a pair */
DEF BAM_FPAIRED       =1
## @abstract the read is mapped in a proper pair */
DEF BAM_FPROPER_PAIR  =2
## @abstract the read itself is unmapped; conflictive with BAM_FPROPER_PAIR */
DEF BAM_FUNMAP        =4
## @abstract the mate is unmapped */
DEF BAM_FMUNMAP       =8
## @abstract the read is mapped to the reverse strand */
DEF BAM_FREVERSE      =16
## @abstract the mate is mapped to the reverse strand */
DEF BAM_FMREVERSE     =32
## @abstract this is read1 */
DEF BAM_FREAD1        =64
## @abstract this is read2 */
DEF BAM_FREAD2       =128
## @abstract not primary alignment */
DEF BAM_FSECONDARY   =256
## @abstract QC failure */
DEF BAM_FQCFAIL      =512
## @abstract optical or PCR duplicate */
DEF BAM_FDUP        =1024

DEF BAM_CIGAR_SHIFT=4
DEF BAM_CIGAR_MASK=((1 << BAM_CIGAR_SHIFT) - 1)

######################################################################
# valid types for sam headers
VALID_HEADER_TYPES = { "HD" : dict,
                       "SQ" : list,
                       "RG" : list,
                       "PG" : list,
                       "CO" : list }

# order of records within sam headers
VALID_HEADERS = ("HD", "SQ", "RG", "PG", "CO" )

# type conversions within sam header records
VALID_HEADER_FIELDS = { "HD" : { "VN" : str, "SO" : str, "GO" : str },
                        "SQ" : { "SN" : str, "LN" : int, "AS" : str, "M5" : str, "UR" : str, "SP" : str },
                        "RG" : { "ID" : str, "SM" : str, "LB" : str, "DS" : str, "PU" : str, "PI" : str, "CN" : str, "DT" : str, "PL" : str, "PG" : str},
                        "PG" : { "ID" : str, "VN" : str, "CL" : str , "PN" : str, "PP" : str, "DS": str},  }

# output order of fields within records
VALID_HEADER_ORDER = { "HD" : ( "VN", "SO", "GO" ),
                       "SQ" : ( "SN", "LN", "AS", "M5" , "UR" , "SP" ),
                       "RG" : ( "ID", "SM", "LB", "DS" , "PU" , "PI" , "CN" , "DT", "PL" , "PG"),
                       "PG" : ( "ID", "VN", "CL" ), }

######################################################################
## Public methods
######################################################################

cdef char* get_base(uint8_t *s, int i):
    cdef uint8_t idx = bam_seqi(s, i)
    return "=ACMGRSVTWYHKDBN"[idx]

def unpickleSamfile(samFilePointer, indexFilePointer, indexesByRegion, isbam, filename):
    """
    """
    #logger.debug("Calling Samfile unpickle")
    cdef Samfile theFile = Samfile(filename)

    if samFilePointer is None:
        theFile.samfile = NULL
    else:
        theFile.samfile = <samfile_t*>(samFilePointer)

    if indexFilePointer is None:
        theFile.index = NULL
    else:
        theFile.index = <bam_index_t*>(indexFilePointer)

    theFile.indexesByRegion = indexesByRegion
    theFile.filename = filename

    #logger.debug("Done Calling Samfile unpickle")

    return theFile

cdef class Samfile:
    """
    *(filename, mode='r', referencenames = None, referencelengths = None, text = NULL, header = None)*

    A *SAM* file. The file is automatically opened.

    *mode* should be ``r`` for reading or ``w`` for writing. The default is text mode so for binary
    (:term:`BAM`) I/O you should append ``b`` for compressed or ``u`` for uncompressed :term:`BAM` output.
    Use ``h`` to output header information  in text (:term:`TAM`)  mode.

    If ``b`` is present, it must immediately follow ``r`` or ``w``.
    Currently valid modes are ``r``, ``w``, ``wh``, ``rb``, ``wb`` and ``wbu``.

    so to open a :term:`BAM` file for reading::

        f=Samfile('ex1.bam','rb')


    For writing, the header of a :term:`TAM` file/:term:`BAM` file can be constituted from several
    sources:

        2. If *header* is given, the header is build from a multi-level dictionary. The first level are the four types ('HD', 'SQ', ...). The second level is then a list of lines, with each line being a list of tag-value pairs.

        3. If *text* is given, new header text is copied from raw text.

        4. The names (*referencenames*) and lengths (*referencelengths*) are supplied directly as lists.

    If an index for a BAM file exists (.bai), it will be opened automatically. Without an index random
    access to reads via :meth:`fetch` and :meth:`pileup` is disabled.
    """
    def __cinit__(self, fileName):
        """
        Constructor.
        """
        self.samfile = NULL
        self.header = NULL
        self.filename = fileName
        self.index = NULL
        self.indexesByRegion = {}
        self.lock = None
    
    def __dealloc__( self ):
        """
        Clean up. Here I've pasted the code from several clean-up
        methods, as calling methods is not guaranteed to work in
        this function.
        """
        # remember: dealloc cannot call other methods
        # Note that __del__ is not called.
        if self.samfile != NULL:
            sam_close(self.samfile)
            self.samfile = NULL
        
        if self.index != NULL:
            hts_idx_destroy(self.index);
            self.index = NULL
        
        if self.header != NULL:
            bam_hdr_destroy(self.header)
            self.header = NULL
    
    cdef void clearHeader(self):
        """
        Clear all the header data.
        """
        if self.header != NULL:
            bam_hdr_destroy(self.header)
            self.header = NULL
    
    cdef void clearIndex(self):
        """
        Clear all the index file data
        """
        if self.index != NULL:
            hts_idx_destroy(self.index);
            self.index = NULL
    
    cdef void clearCache(self):
        """
        Each BGZF file has 2 64kb buffers. If required, I clear these and set them to NULL,
        when not loading data. They are then allocated and freed each time data is loaded, to
        avoid memory problems when keeping large numbers of files open at the same time. The
        standard C i/o buffer is 4kb per file, which is small enough to not be a problem for now.
        """
        if self.samfile.fp.bgzf.uncompressed_block != NULL:
            free(self.samfile.fp.bgzf.uncompressed_block)
            self.samfile.fp.bgzf.uncompressed_block = NULL

        if self.samfile.fp.bgzf.compressed_block != NULL:
            free(self.samfile.fp.bgzf.compressed_block)
            self.samfile.fp.bgzf.compressed_block = NULL

    cdef void createCache(self):
        """
        If we have previously cleared the cache, then create it again.
        """
        if self.samfile != NULL:
            if self.samfile.fp.bgzf.uncompressed_block == NULL:
                self.samfile.fp.bgzf.uncompressed_block = malloc(64*1024)

                if self.samfile.fp.bgzf.uncompressed_block == NULL:
                    logger.error("Could not allocate BAM BGZF buffer in Samfile.createCache()")
                    raise StandardError, "Out of memory"

            if self.samfile.fp.bgzf.compressed_block == NULL:
                self.samfile.fp.bgzf.compressed_block = malloc(64*1024)

                if self.samfile.fp.bgzf.compressed_block == NULL:
                    logger.error("Could not allocate BAM BGZF buffer in Samfile.createCache()")
                    raise StandardError, "Out of memory"

    def __reduce__(self):
        """
        This function is called by Pickle/cPickle in order to grab the state of this
        object.
        """
        #logger.debug("Calling Samfile reduce")
        cdef int samFilePointer = <int>(self.samfile)
        cdef int indexFilePointer = <int>(self.index)
        #print "Pointer values are as follows: Sam file = %s. Index = %s" %(samFilePointer, indexFilePointer)

        data = (samFilePointer, indexFilePointer, self.indexesByRegion, self.is_bam(), self.filename)

        if self.samfile == NULL:
            data = (None,) + data[1:]

        if self.index == NULL:
            data = data[0:1] + (None,) + data[2:]

        #logger.debug("Done Calling Samfile reduce")
        return (unpickleSamfile, data)
    
    cdef bool is_bam(self):
        return self.samefile.is_bin
    
    cdef bool is_cram(self):
        return self.samefile.is_cram
    
    cdef int _isOpen(self):
        '''return true if samfile has been opened.'''
        return self.samfile != NULL
    
    cdef _hasIndex(self):
        '''return true if samfile has an existing (and opened) index.'''
        return self.index != NULL
    
    cdef void openBAMFile(self, mode):
        """
        """
        self.samfile = sam_open(self.filename, mode)
        self.header = bam_hdr_init()
        self.header = sam_hdr_read(self.samfile);
    
    cpdef _open(self, mode, loadIndex=True):
        """
        open a sam/bam/cram file.

        If _open is called on an existing bamfile, the current file will be
        closed and a new file will be opened.
        """
        if mode not in ("r", "rb", "rbh"):
            raise StandardError, "invalid file opening mode `%s`" % mode
        
        if self.samfile != NULL:
            if loadIndex and self.index == NULL:
                self.index = sam_index_load(self.samfile, self.filename)
            return

        if mode[0] == "r":
            self.openBAMFile(mode)
        else:
            raise StandardError, "BAM file is read-only"

        if self.samfile == NULL:
            raise IOError("Could not open file `%s`. Check that file/path exists." % self.filename)

        if mode[0] == "r" and self.is_bam():
            # returns NULL if there is no index or index could not be opened
            if loadIndex and self.index == NULL:
                self.index = sam_index_load(self.samfile, self.filename)

                if self.index == NULL:
                    raise IOError("Error while opening index for file `%s`. Check that index exists " % self.filename)
    
    cdef void loadOffsetsForRegions(self, list regions):
        """
        Load, and maintain an internal cache, of the file off-sets needed to access
        the beginning of the regions specified here. This will be used in preference to
        accessing the bam index, since loading the whole index is expensive. These offsets
        can be computed once at the beginning of the process.
        """
        cdef IteratorRow theIterator
        cdef samFile* fp = self.samfile.x.bam
        cdef int rtid = -1
        cdef int rstart = -1
        cdef int rend = -1
        cdef int n_off = -1
        cdef int index = 0

        if not self._isOpen():
            self._open('rbh', loadIndex=True) # Change to False

        if self.index == NULL:
            self.index = sam_index_load(self.samfile, self.filename)

        if self.index == NULL:
            raise IOError("Error while opening index for file `%s`. Check that index exists " % self.filename)

        cdef int nRegions = len(regions)

        for index,(chrom,start,end) in enumerate(regions):
            region, rtid, rstart, rend = self._parseRegion(chrom, start, end)
            theIterator = IteratorRow(self, rtid, rstart, rend, useIndex=1)
            #logger.info("Region index = %s (%s:%s-%s). Noff = %s" %(index, chrom, start, end, theIterator.bam_iter.n_off))
            self.indexesByRegion["%s:%s-%s" %(chrom,start,end)] = theIterator

        if self.index != NULL:
            hts_idx_destroy(self.index);
            self.index = NULL
    
    cdef char* getrname(self, int tid):
        '''
        Convert numerical :term:`tid` into :ref:`reference` name.
        '''
        if not 0 <= tid < self.header.n_targets:
            raise ValueError( "tid (%s) out of range 0<=tid<%i" % (tid, self.header.n_targets))
        
        return self.header.target_name[tid]
    
    cdef _parseRegion(self, reference=None, start=None, end=None, region=None):
        """
        parse region information.
        
        raise Value for for invalid regions.
        
        returns a tuple of region, tid, start and end. Region
        is a valid samtools :term:`region` or None if the region
        extends over the whole file.
        
        Note that regions are 1-based, while start,end are python coordinates.
        """
        cdef int rtid
        cdef int rstart
        cdef int rend
        cdef int max_pos
        
        max_pos = 2 << 29
        rtid = rstart = rend = 0
        
        # translate to a region
        if reference:
            if start != None and end != None:
                region = "%s:%i-%i" % (reference, max(1, start), end)
            else:
                region = reference
        
        if region:
            cdef const char *name_lim = hts_parse_reg(region, rstart, rend)
            cdef char *name = <char*>malloc(name_lim - region + 1)
            memcpy(name, reg, name_lim - reg)
            name[name_lim - reg] = '\0'
            cdef int tid = bam_name2id(header, name)
            free(name)
            
            if rtid < 0:
                raise ValueError("invalid region `%s`" % region)
            
            if rstart > rend:
                raise ValueError('invalid region: start (%i) > end (%i)' % (rstart, rend))
            
            if not 0 <= rstart < max_pos:
                raise ValueError('start out of range (%i)' % rstart)
            
            if not 0 <= rend < max_pos:
                raise ValueError('end out of range (%i)' % rend)
        
        return region, rtid, rstart, rend
    
    cpdef IteratorRow fetch(self, char* reference, int start, int end):
        '''
        Fetch reads from a specified region.
        '''
        cdef int rtid
        cdef int rstart
        cdef int rend
        cdef int cacheIndex = 0
        cdef bam_fetch_iterator_t* theIterator
        cdef IteratorRow theIt

        # Load from pre-existing iterator
        if "%s:%s-%s" %(reference, start, end) in self.indexesByRegion:
            theIt = self.indexesByRegion["%s:%s-%s" %(reference, start, end)]

            # Don't need to load the index, since we already have the pre-computed offsets
            if not self._isOpen():
                self._open('rb', loadIndex=False)

            # Re-create cache for BGZF file, if required.
            self.createCache()

            if self.is_bam():
                theIt.bam_iter.fp = self.samfile.x.bam
                return theIt
            else:
                raise StandardError, "Random access query only allowed for BAM files."

        # New query
        else:
            # Need to load the index for new queries.
            if not self._isOpen():
                self._open('rbh', loadIndex=True) # Change to False

            # Re-create cache for BGZF file, if required.
            self.createCache()

            if self.isbam:
                region, rtid, rstart, rend = self._parseRegion(reference, start, end)
                return IteratorRow(self, rtid, rstart, rend, useIndex=1)
            else:
                raise StandardError, "Random access query only allowed for BAM files."


cdef class SamIterator:
    """
    Iterates over mapped reads in a region.
    """
    def __cinit__(self, Samfile samfile, const char *region):
        """
        Constructor.
        """
        self.theIterator = NULL
        self.b = NULL
        
        if not samfile._isOpen():
            raise StandardError, "BAM file %s is not open. Cannot read from file." %(samfile.filename)

        self.theIterator = sam_itr_querys(samfile.index, samfile.header, region)
        self.b = bam_init1()
    