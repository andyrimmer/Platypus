import os
import types
import itertools
import logging
cimport cython
import zlib

###################################################################################################

logger = logging.getLogger("Log")

###################################################################################################

# Structs etc from samtools and C libtaries

cdef extern from "stdlib.h":
  void free(void *)
  void *alloca(size_t)
  void *calloc(size_t,size_t)
  int c_abs "abs" (int)

cdef extern from "math.h":
    double exp(double)
    double log(double)
    double log10(double)

cdef extern from "stdio.h":
    ctypedef struct FILE:
      pass
    FILE *fopen(char *,char *)
    FILE *freopen(char *path, char *mode, FILE *stream)
    int fileno(FILE *stream)
    int dup2(int oldfd, int newfd)
    int fflush(FILE *stream)
    FILE * stderr
    FILE * stdout
    int fclose(FILE *)
    int sscanf(char *str,char *fmt,...)
    int printf(char *str,char *fmt,...)
    int sprintf(char *str,char *fmt,...)
    int fprintf(FILE *ifile,char *fmt,...)
    char *fgets(char *str,int size,FILE *ifile)

cdef extern from "ctype.h":
  int toupper(int c)

cdef extern from "string.h":
  ctypedef int size_t
  void *memcpy(void *dst,void *src,size_t len)
  int strncmp(char *s1,char *s2,size_t len)
  char *strncpy(char *dest,char *src, size_t len)
  size_t strlen(char *s)
  int memcmp( void * s1, void *s2, size_t len )

cdef extern from "stdint.h":
    ctypedef int int64_t
    ctypedef int int32_t
    ctypedef int uint32_t
    ctypedef int uint8_t
    ctypedef int uint64_t

cdef extern from "bgzf.h":
    ctypedef struct BGZF:
        pass

cdef extern from "bam.h":

    ctypedef struct tamFile:
        pass

    ctypedef struct bamFile:
        void* uncompressed_block
        void* compressed_block

    bamFile* bam_open(char*, char*)

    ctypedef struct bam_header_t:
       int32_t n_targets
       char **target_name
       uint32_t *target_len
       void *hash
       void *rg2lib
       int l_text
       char *text

    ctypedef struct bam1_core_t:
        int32_t tid
        int32_t pos
        uint32_t bin
        uint32_t qual
        uint32_t l_qname
        uint32_t flag
        uint32_t n_cigar
        int32_t l_qseq
        int32_t mtid
        int32_t mpos
        int32_t isize

    ctypedef struct bam1_t:
      bam1_core_t core
      int l_aux
      int data_len
      int m_data
      uint8_t *data

    ctypedef struct bam_index_t:
        pass

    ctypedef int (*bam_fetch_f)(bam1_t *b, void *data)

    bamFile* razf_dopen(int data_fd, char *mode)
    bam1_t * bam_dup1( bam1_t *src )
    bam1_t * bam_copy1(bam1_t *bdst, bam1_t *bsrc)
    bam_index_t *bam_index_load(char *f )
    void bam_index_destroy(bam_index_t *idx)
    int bam_parse_region(bam_header_t *header, char *str, int *ref_id, int *begin, int *end)
    int bam_fetch(bamFile* fp, bam_index_t *idx, int tid, int beg, int end, void *data, bam_fetch_f func)
    int bam_read1(bamFile* fp, bam1_t *b)
    int bam_write1( bamFile* fp, bam1_t *b)
    bam_header_t *bam_header_init()
    int bam_header_write( bamFile* fp, bam_header_t *header)
    bam_header_t *bam_header_read( bamFile* fp )
    void bam_header_destroy(bam_header_t *header)
    uint8_t *bam_aux_get(bam1_t *b,  char tag[2])
    int bam_aux2i(uint8_t *s)
    float bam_aux2f(uint8_t *s)
    double bam_aux2d(uint8_t *s)
    char bam_aux2A( uint8_t *s)
    char *bam_aux2Z( uint8_t *s)
    int bam_reg2bin(uint32_t beg, uint32_t end)
    uint32_t bam_calend(bam1_core_t *c, uint32_t *cigar)

cdef extern from "sam.h":

  ctypedef struct samfile_t_un:
    tamFile tamr
    bamFile* bam
    FILE *tamw

  ctypedef struct samfile_t:
     int type
     samfile_t_un x
     bam_header_t *header

  samfile_t *samopen( char *fn, char * mode, void *aux)
  void samclose(samfile_t *fp)
  int samread(samfile_t *fp, bam1_t *b)
  int samwrite(samfile_t *fp, bam1_t *b)

cdef extern from "pysam_util.h":

    int pysam_dispatch(int argc, char *argv[] )

    # stand-in functions for samtools macros
    void pysam_bam_destroy1( bam1_t * b)

    # add *nbytes* into the variable length data of *src* at *pos*
    bam1_t * pysam_bam_update( bam1_t * b,
                               size_t nbytes_old,
                               size_t nbytes_new,
                               uint8_t * pos )

    # translate char to unsigned char
    unsigned char pysam_translate_sequence( char s )

    # stand-ins for samtools macros
    uint32_t* pysam_bam1_cigar(bam1_t * b)
    char * pysam_bam1_qname( bam1_t * b)
    uint8_t * pysam_bam1_seq( bam1_t * b)
    uint8_t * pysam_bam1_qual( bam1_t * b)
    uint8_t * pysam_bam1_aux( bam1_t * b)

    ctypedef struct pair64_t:
        uint64_t u
        uint64_t v

    # iterator implemenation
    ctypedef struct bam_fetch_iterator_t:
        bam1_t* b
        pair64_t* off
        int n_off
        uint64_t curr_off
        int curr_chunk
        bamFile* fp
        int tid
        int beg
        int end
        int n_seeks

    bam_fetch_iterator_t* bam_init_fetch_iterator(bamFile* fp, bam_index_t *idx, int tid, int beg, int end)
    bam1_t* bam_fetch_iterate(bam_fetch_iterator_t *iter)
    void bam_cleanup_fetch_iterator(bam_fetch_iterator_t *iter)

###################################################################################################

cdef int TYPE_BAM  = 1
cdef int TYPE_READ = 2

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
                        "RG" : { "ID" : str, "SM" : str, "LB" : str, "DS" : str, "PU" : str, "PI" : str, "CN" : str, "DT" : str, "PL" : str, },
                        "PG" : { "ID" : str, "VN" : str, "CL" : str , "PN" : str, "PP" : str}, }

# output order of fields within records
VALID_HEADER_ORDER = { "HD" : ( "VN", "SO", "GO" ),
                       "SQ" : ( "SN", "LN", "AS", "M5" , "UR" , "SP" ),
                       "RG" : ( "ID", "SM", "LB", "DS" , "PU" , "PI" , "CN" , "DT", "PL" ),
                       "PG" : ( "ID", "VN", "CL" ), }

######################################################################
## Public methods
######################################################################

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
    theFile.isbam = isbam
    theFile.filename = filename

    #logger.debug("Done Calling Samfile unpickle")

    return theFile

######################################################################

cdef bam_index_t* my_bam_index_load(char *f):
    """
    Wrapper for profiling.
    """
    return bam_index_load(f)

######################################################################

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
        self.isbam = False
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
            samclose(self.samfile)
            self.samfile = NULL

        if self.index != NULL:
            bam_index_destroy(self.index);
            self.index = NULL

    cdef void clearHeader(self):
        """
        Clear all the header data.
        """
        if self.samfile.header != NULL:
            bam_header_destroy(self.samfile.header)
            self.samfile.header = NULL

    cdef void clearIndex(self):
        """
        Clear all the index file data
        """
        if self.index != NULL:
            bam_index_destroy(self.index);
            self.index = NULL

    cdef void clearCache(self):
        """
        Each BGZF file has 2 64kb buffers. If required, I clear these and set them to NULL,
        when not loading data. They are then allocated and freed each time data is loaded, to
        avoid memory problems when keeping large numbers of files open at the same time. The
        standard C i/o buffer is 4kb per file, which is small enough to not be a problem for now.
        """
        if self.samfile.x.bam.uncompressed_block != NULL:
            free(self.samfile.x.bam.uncompressed_block)
            self.samfile.x.bam.uncompressed_block = NULL

        if self.samfile.x.bam.compressed_block != NULL:
            free(self.samfile.x.bam.compressed_block)
            self.samfile.x.bam.compressed_block = NULL

    cdef void createCache(self):
        """
        If we have previously cleared the cache, then create it again.
        """
        if self.samfile != NULL:
            if self.samfile.x.bam.uncompressed_block == NULL:
                self.samfile.x.bam.uncompressed_block = malloc(64*1024)

                if self.samfile.x.bam.uncompressed_block == NULL:
                    logger.error("Could not allocate BAM BGZF buffer in Samfile.createCache()")
                    raise StandardError, "Out of memory"

            if self.samfile.x.bam.compressed_block == NULL:
                self.samfile.x.bam.compressed_block = malloc(64*1024)

                if self.samfile.x.bam.compressed_block == NULL:
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

        data = (samFilePointer, indexFilePointer, self.indexesByRegion, self.isbam, self.filename)

        if self.samfile == NULL:
            data = (None,) + data[1:]

        if self.index == NULL:
            data = data[0:1] + (None,) + data[2:]

        #logger.debug("Done Calling Samfile reduce")
        return (unpickleSamfile, data)

    cdef int _isOpen( self ):
        '''return true if samfile has been opened.'''
        return self.samfile != NULL

    cdef _hasIndex( self ):
        '''return true if samfile has an existing (and opened) index.'''
        return self.index != NULL

    cdef void openBAMFile(self, mode):
        """
        """
        self.samfile = <samfile_t*>calloc(1, sizeof(samfile_t))

        if mode[0] == 'r':
            self.samfile.type |= TYPE_READ

            if len(mode) > 1 and mode[1] == 'b':
                self.samfile.type |= TYPE_BAM
                self.samfile.x.bam = bam_open(self.filename, "r")

                if <BGZF*>(self.samfile.x.bam) == NULL:
                    raise StandardError, "Could not open BAM file %s" %(self.filename)
        else:
            raise StandardError, "Writing to BAM files not supported"

        if len(mode) > 2 and mode[2] == 'h':
            self.samfile.header = bam_header_read(self.samfile.x.bam)
        else:
            self.samfile.header = NULL

    cpdef _open(self, mode, loadIndex=True):
        """
        open a sam/bam file.

        If _open is called on an existing bamfile, the current file will be
        closed and a new file will be opened.
        """
        if mode not in ("r","rb", "rbh"):
            raise StandardError, "invalid file opening mode `%s`" % mode

        # If file is already open, do nothing.
        if self.samfile != NULL:
            if loadIndex and self.index == NULL:
                self.index = my_bam_index_load(self.filename)

            return

        self.samfile = NULL
        self.isbam = len(mode) > 1 and mode[1] == 'b'

        if mode[0] == "r":
            self.openBAMFile(mode)
        else:
            raise StandardError, "BAM file is read-only"

        if self.samfile == NULL:
            raise IOError("Could not open file `%s`. Check that file/path exists." % self.filename )

        if mode[0] == "r" and self.isbam:
            # returns NULL if there is no index or index could not be opened
            if loadIndex and self.index == NULL:
                self.index = my_bam_index_load(self.filename)

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
        cdef bamFile* fp = self.samfile.x.bam
        cdef int rtid = -1
        cdef int rstart = -1
        cdef int rend = -1
        cdef int n_off = -1
        cdef int index = 0

        if not self._isOpen():
            self._open('rbh', loadIndex=True) # Change to False

        if self.index == NULL:
            self.index = my_bam_index_load(self.filename)

        if self.index == NULL:
            raise IOError("Error while opening index for file `%s`. Check that index exists " % self.filename)

        cdef int nRegions = len(regions)

        for index,(chrom,start,end) in enumerate(regions):
            region, rtid, rstart, rend = self._parseRegion(chrom, start, end)
            theIterator = IteratorRow(self, rtid, rstart, rend, useIndex=1)
            #logger.info("Region index = %s (%s:%s-%s). Noff = %s" %(index, chrom, start, end, theIterator.bam_iter.n_off))
            self.indexesByRegion["%s:%s-%s" %(chrom,start,end)] = theIterator

        if self.index != NULL:
            bam_index_destroy(self.index);
            self.index = NULL

    cdef char* getrname(self, int tid):
        '''
        Convert numerical :term:`tid` into :ref:`reference` name.
        '''
        if not 0 <= tid < self.samfile.header.n_targets:
            raise ValueError( "tid (%s) out of range 0<=tid<%i" % (tid, self.samfile.header.n_targets))

        return self.samfile.header.target_name[tid]

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
                region = "%s:%i-%i" % (reference, start+1, end)
            else:
                region = reference

        if region:
            bam_parse_region(self.samfile.header, region, &rtid, &rstart, &rend)

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

            if self.isbam:
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

    cpdef IteratorRowAll fetchAllReads(self):
        '''
        Return all reads from the file, including un-mapped reads.
        '''
        if not self._isOpen():
            raise ValueError("BAM file is not open. Cannot return reads.")

        if self.isbam:
            return IteratorRowAll(self)

    cpdef close(self):
        '''
        closes file.
        '''
        if self.samfile != NULL:
            samclose(self.samfile)
            self.samfile = NULL

        if self.index != NULL:
            bam_index_destroy(self.index);
            self.index = NULL

    property nreferences:
        '''number of :term:`reference` sequences in the file.'''
        def __get__(self):
            return self.samfile.header.n_targets

    property references:
        """tuple with the names of :term:`reference` sequences."""
        def __get__(self):
            t = []
            for x from 0 <= x < self.samfile.header.n_targets:
                t.append( self.samfile.header.target_name[x] )
            return tuple(t)

    property lengths:
        """tuple of the lengths of the :term:`reference` sequences. The lengths are in the same order as :attr:`pysam.Samfile.reference`
        """
        def __get__(self):
            t = []
            for x from 0 <= x < self.samfile.header.n_targets:
                t.append( self.samfile.header.target_len[x] )
            return tuple(t)

    property text:
        '''full contents of the :term:`sam file` header as a string.'''
        def __get__(self):
            # create a temporary 0-terminated copy
            cdef char * t
            t = <char*>calloc( self.samfile.header.l_text + 1, sizeof(char) )
            memcpy( t, self.samfile.header.text, self.samfile.header.l_text )
            cdef bytes result = t
            free(t)
            return result

    property header:
        '''header information within the :term:`sam file`. The records and fields are returned as
        a two-level dictionary.
        '''
        def __get__(self):
            result = {}

            if self.samfile.header.text != NULL:
                # convert to python string (note: call self.text to create 0-terminated string)
                t = self.text
                for line in t.split("\n"):
                    if not line.strip(): continue

                    if not line.startswith("@"):
                        raise StandardError, "Header line without '@': '%s. Total header text is %s'" % (line,t)

                    fields = line[1:].split("\t")
                    record = fields[0]
                    assert record in VALID_HEADER_TYPES, "header line with invalid type '%s': '%s'" % (record, line)

                    # treat comments
                    if record == "CO":
                        if record not in result: result[record] = []
                        result[record].append( "\t".join( fields[1:] ) )
                        continue

                    # the following is clumsy as generators do not work?
                    x = {}
                    for field in fields[1:]:
                        key, value = field.split(":",1)
                        if key not in VALID_HEADER_FIELDS[record]:
                            raise ValueError( "unknown field code '%s' in record '%s'" % (key, record) )
                        x[key] = VALID_HEADER_FIELDS[record][key](value)

                    if VALID_HEADER_TYPES[record] == dict:
                        if record in result:
                            raise ValueError( "multiple '%s' lines are not permitted" % record )
                        result[record] = x
                    elif VALID_HEADER_TYPES[record] == list:
                        if record not in result: result[record] = []
                        result[record].append( x )

            return result

###################################################################################################
## turning callbacks elegantly into iterators is an unsolved problem, see the following threads:
## http://groups.google.com/group/comp.lang.python/browse_frm/thread/0ce55373f128aa4e/1d27a78ca6408134?hl=en&pli=1
## http://www.velocityreviews.com/forums/t359277-turning-a-callback-function-into-a-generator.html
## Thus I chose to rewrite the functions requiring callbacks. The downside is that if the samtools C-API or code
## changes, the changes have to be manually entered.
###################################################################################################

def unpickleIteratorRow(samfile, bamIterPointer, useIndex):
    """
    """
    logger.debug("Calling IteratorRow unpickle")

    cdef IteratorRow theIt = IteratorRow(samfile, 0, 0, 0, 0)

    if bamIterPointer is None:
        theIt.bam_iter = NULL
    else:
        theIt.bam_iter = <bam_fetch_iterator_t*>(bamIterPointer)

    theIt.useIndex = useIndex
    logger.debug("Done Calling IteratorRow unpickle")

    return theIt

######################################################################

@cython.profile(False)
cdef class IteratorRow:
    """
    Iterates over mapped reads in a region.
    """
    def __cinit__(self, Samfile samfile, int tid, int beg, int end, int useIndex=1):
        """
        Constructor.
        """
        self.bam_iter = NULL
        self.b = NULL
        self.useIndex = useIndex

        if not samfile._isOpen():
            raise StandardError, "BAM file %s is not open. Cannot read from file." %(samfile.filename)

        cdef bamFile* fp = samfile.samfile.x.bam

        # Load from BAM by querying index for file off-set
        if useIndex:
            if not samfile._hasIndex():
                raise StandardError, "Cannot retrieve random region from BAM file %s, as it does not have a BAM index" %(samfile.filename)

            self.bam_iter = bam_init_fetch_iterator(fp, samfile.index, tid, beg, end)

        # Load from pre-existing BAM iterator. I can't pass that in here, so it will be shoved
        # in from outside after construction.
        else:
            pass

    def __reduce__(self):
        """
        This function is called by Pickle/cPickle in order to grab the state of this
        object.
        """
        logger.debug("Calling IteratorRow reduce")

        cdef int bamIter = <int>(self.bam_iter)
        cdef int useIndex = 0

        data = (self.samfile, bamIter, useIndex)

        if self.bam_iter == NULL:
            data = data[0:1] + (None,) + data[2:]

        logger.debug("Done Calling IteratorRow reduce")
        return (unpickleIteratorRow, data)

    cdef int cnext(self):
        '''
        cversion of iterator. Used by IteratorColumn
        '''
        self.b = bam_fetch_iterate(self.bam_iter)

        if self.b == NULL:
            return 0

        return 1

    def __dealloc__(self):
        '''
        remember: dealloc cannot call other methods!
        '''
        # Only want to free memory from iterators which were created
        # by this object, not those inserted from cache.
        if self.useIndex:
            if self.bam_iter:
                bam_cleanup_fetch_iterator(self.bam_iter)
                free(self.bam_iter)

###################################################################################################

cdef class IteratorRowAll:
    """
    iterates over all mapped reads
    """
    def __cinit__(self, Samfile samfile):

        if not samfile._isOpen():
            raise StandardError, "BAM file %s is not open. Cannot read from file." %(samfile.filename)

        self.fp = samfile.samfile
        self.b = <bam1_t*>calloc(1, sizeof(bam1_t))

    def __iter__(self):
        return self

    cdef bam1_t* getCurrent(self):
        return self.b

    cdef int cnext(self):
        '''
        cversion of iterator. Used by IteratorColumn
        '''
        return samread(self.fp, self.b)

    def __next__(self):
        """python version of next().

        pyrex uses this non-standard name instead of next()
        """
        cdef int ret = samread(self.fp, self.b)

        if ret > 0:
            return makeAlignedRead(self.b)
        else:
            raise StopIteration

    def __dealloc__(self):
        '''remember: dealloc cannot call other methods!'''
        pysam_bam_destroy1(self.b)

###################################################################################################

cdef void reOrientRead(cAlignedRead* theRead):
    """
    This function reverse-complements a read, as well as reversing its quality score sequence,
    and the cigar-string sequence.
    """
    cdef int i = 0
    cdef char tempSeq = 0
    cdef char tempQual = 0

    # Reverse sequence and quality scores
    for i in range((theRead.rlen-1)/2):
        tempSeq = theRead.seq[i]
        tempQual = theRead.qual[i]
        theRead.seq[i] = theRead.seq[(theRead.rlen-1) - i]
        theRead.qual[i] = theRead.qual[(theRead.rlen-1) - i]
        theRead.seq[(theRead.rlen-1)-i] = tempSeq
        theRead.qual[(theRead.rlen-1)-i] = tempQual

    cdef int tempOp = 0
    cdef int tempLen = 0

    # Reverse cigar string.
    for i in range((theRead.cigarLen-1)/2):
        tempOp = theRead.cigarOps[i]
        tempLen = theRead.cigarOps[i+1]
        theRead.cigarOps[i] = theRead.cigarOps[(theRead.cigarLen-2)-i]
        theRead.cigarOps[i+1] = theRead.cigarOps[(theRead.cigarLen-1)-i]
        theRead.cigarOps[(theRead.cigarLen-2)-i] = tempOp
        theRead.cigarOps[(theRead.cigarLen-1)-i] = tempLen

    # Toggle reverse bit
    if Read_IsReverse(theRead):
        Read_SetIsNotReverse(theRead)
    else:
        Read_SetIsReverse(theRead)

###################################################################################################

cdef cAlignedRead* createRead(bam1_t* src, int storeRgID, char** rgID):
    """
    Allocate memory for aligned read struct, and populate it.
    """
    cdef char* bamTable = "=ACMGRSVTWYHKDBN"
    cdef int characterIndex = 0
    cdef int lenSeq = src.core.l_qseq
    cdef int lenQual = src.core.l_qseq
    cdef double mLTOT = -0.23025850929940459    # Minus log ten over ten

    cdef uint8_t* p = pysam_bam1_seq(src)
    cdef uint8_t* q = pysam_bam1_qual(src)

    if lenSeq == 0:
        return NULL

    if lenQual == 0:
        return NULL

    if q[0] == 0xff:
        return NULL

    cdef cAlignedRead* theRead = <cAlignedRead*>(malloc(sizeof(cAlignedRead)))

    cdef char* seq = <char*>calloc(lenSeq + 1 , sizeof(char))
    cdef char* qual = <char*>calloc(lenQual + 1, sizeof(char))

    # Try to grab the read-group tag value
    cdef uint8_t* v = NULL
    cdef char* tempRgID = NULL
    cdef int lenRgID = 0

    if storeRgID:
        v = bam_aux_get(src, "RG")
        tempRgID = <char*>bam_aux2Z(v)
        lenRgID = strlen(tempRgID)
        rgID[0] = <char*>(calloc(lenRgID+1, sizeof(char)))
        strcpy(rgID[0], tempRgID)

    assert theRead != NULL
    assert seq != NULL
    assert qual != NULL

    cdef int index = 0
    #cdef int nQualsBelow20 = 0

    for index from 0 <= index < lenSeq:
        characterIndex = ((p)[(index) / 2] >> 4 * (1 - (index) % 2) & 0xf)
        seq[index] = bamTable[characterIndex]
        qual[index] = q[index]

        assert qual[index] <= 93
        assert qual[index] >= 0

        #if qual[index] < 20:
        #    nQualsBelow20 += 1

    theRead.seq = seq
    theRead.qual = qual
    theRead.chromID = src.core.tid
    theRead.pos = src.core.pos
    theRead.rlen = src.core.l_qseq
    theRead.mapq = src.core.qual
    theRead.cigarLen = src.core.n_cigar

    cdef int flag = src.core.flag
    cdef int end = theRead.pos
    cdef short* cigarOps = <short*>(malloc(2*theRead.cigarLen*sizeof(short)))
    cdef int cigarFlag = 0
    cdef int cigarFlagLen = 0

    assert cigarOps != NULL

    cdef uint32_t* cigar_p = pysam_bam1_cigar(src)

    for index from 0 <= index < theRead.cigarLen:

        cigarFlag = cigar_p[index] & BAM_CIGAR_MASK
        cigarFlagLen = cigar_p[index] >> BAM_CIGAR_SHIFT
        cigarOps[2*index] = cigarFlag
        cigarOps[ (2*index) + 1] = cigarFlagLen

        # 0 is match. 2 is deletion. 3 is reference skipping
        if cigarFlag == 0 or cigarFlag == 2 or cigarFlag == 3:
            end += cigarFlagLen

        # Soft-clipping of sequence at start of read changes the mapping
        # position. Recorded mapping pos is that of the first aligned (not soft-clipped)
        # base. I want to adjust this so that the read start refers to the first base in
        # the read.
        if index == 0 and cigarFlag == 4:
            theRead.pos -= cigarFlagLen
        # Don't want to pretend that the read extends to the end of the soft-clipped
        # sequence. Treat it in the same way as an insertion.
        elif cigarFlag == 4:
            pass

    theRead.end = end
    theRead.cigarOps = cigarOps
    theRead.bitFlag = flag
    theRead.hash = NULL
    theRead.mateChromID = src.core.mtid
    theRead.insertSize = src.core.isize
    theRead.matePos = src.core.mpos
    Read_SetUnCompressed(theRead)

    return theRead

###################################################################################################

cdef void destroyRead(cAlignedRead* theRead):
    """
    De-allocate memory for read.
    """

    if theRead.seq != NULL:
        free(theRead.seq)

    if theRead.qual != NULL:
        free(theRead.qual)

    if theRead.cigarOps != NULL:
        free(theRead.cigarOps)

    if theRead.hash != NULL:
        free(theRead.hash)

    free(theRead)

###################################################################################################

cdef void compressSeq(cAlignedRead* read, char* refSeq):
    """
    Does exactly what it says on the tin.
    """
    cdef char* seq = read.seq
    cdef char* newSeq = <char*>(alloca(sizeof(char)*read.rlen*2))
    #cdef char* newSeq = <char*>(malloc(sizeof(char)*read.rlen*2))
    cdef char* finalSeq = NULL
    cdef int i = 0
    cdef int nMatches = 0
    cdef int newSeqIndex = 0
    cdef int totalMatches = 0

    for i in range(read.rlen):

        # Ref match. Store count if > 40
        if seq[i] == refSeq[i]:
            totalMatches += 1

            if nMatches == 40:
                newSeq[newSeqIndex] = nMatches
                nMatches = 0
                newSeqIndex += 1

            nMatches += 1

        # Ref mis-match. Store base.
        else:
            if nMatches > 0:
                newSeq[newSeqIndex] = nMatches
                nMatches = 0
                newSeqIndex += 1

            newSeq[newSeqIndex] = seq[i]
            newSeqIndex += 1

    # If we finish off with a load of matches, catch those
    if nMatches > 0:
        newSeq[newSeqIndex] = nMatches
        nMatches = 0
        newSeqIndex += 1

    newSeq[newSeqIndex] = 0
    finalSeq = <char*>malloc(sizeof(char)*(newSeqIndex+1))
    strcpy(finalSeq, newSeq)
    finalSeq[newSeqIndex] = 0
    free(read.seq)
    #free(newSeq)
    read.seq = finalSeq

###################################################################################################

cdef void compressQual(cAlignedRead* read, int qualBinSize):
    """
    Does exactly what it says on the tin.
    """
    cdef int i = 0
    cdef char* qual = read.qual
    #cdef char* newQual = <char*>(malloc(sizeof(char)*2*read.rlen))
    cdef char* newQual = <char*>(alloca(sizeof(char)*2*read.rlen))
    cdef char* finalQual = NULL
    cdef char lastChar = 0
    cdef int lastCount = 0
    cdef int newQualIndex = 0

    if qualBinSize > 1:
        for i in range(read.rlen):
            qual[i] = (qual[i]/qualBinSize)*qualBinSize

    for i in range(read.rlen):

        assert qual[i] >= 0, "Shit 3. Qual = %s" %(qual[i])
        assert qual[i] <= 93, "Shit 4!. Qual = %s" %(qual[i])

        if i == 0:
            newQual[newQualIndex] = qual[i] + 33
            newQualIndex += 1
            lastChar = qual[i]
            lastCount = 1
        else:
            if qual[i] == lastChar:
                lastCount += 1
            else:
                newQual[newQualIndex] = lastCount
                newQualIndex += 1
                newQual[newQualIndex] = qual[i] + 33
                newQualIndex += 1
                lastChar = qual[i]
                lastCount = 1

    if lastCount > 0:
        newQual[newQualIndex] = lastCount
        newQualIndex += 1

    finalQual = <char*>malloc(sizeof(char)*(newQualIndex+1))
    newQual[newQualIndex] = 0
    strcpy(finalQual, newQual)
    finalQual[newQualIndex] = 0
    #free(newQual)
    free(read.qual)
    read.qual = finalQual

###################################################################################################

cdef void uncompressSeq(cAlignedRead* read, char* refSeq):
    """
    Does exactly what it says on the tin.
    """
    cdef int refIndex = 0
    cdef int newSeqIndex = 0
    cdef int i = 0
    cdef int j = 0
    cdef int seqLen = strlen(read.seq)
    cdef char* seq = read.seq
    cdef char* newSeq = <char*>(malloc(sizeof(char*)*(read.rlen+1)))

    for i in range(seqLen):
        if seq[i] <= 40:
            for j in range(seq[i]):
                newSeq[newSeqIndex] = refSeq[j + refIndex]
                newSeqIndex += 1
            refIndex += seq[i]
        else:
            newSeq[newSeqIndex] = seq[i]
            refIndex += 1
            newSeqIndex += 1

    free(read.seq)
    read.seq = newSeq
    read.seq[read.rlen] = 0

###################################################################################################

cdef void uncompressQual(cAlignedRead* read):
    """
    Does exactly what it says on the tin.
    """
    cdef char* newQual = <char*>(malloc(sizeof(char)*(read.rlen+1)))
    cdef char* qual = read.qual

    cdef int lenQual = strlen(read.qual)
    cdef int i = 0
    cdef int j = 0
    cdef int newQualIndex = 0

    #logger.info("Length = %s or %s" %(lenQual, len([x for x in qual])))
    #logger.info(",".join( [str(x) for x in qual] ) )

    for i in range(0, lenQual-1, 2):
        for j in range(qual[i+1]):
            newQual[newQualIndex] = qual[i] - 33
            assert qual[i] - 33 >= 0, "Shit 1. Qual = %s. Count = %s!" %(qual[i], qual[i+1])
            assert qual[i] - 33 <= 93, "Shit 2!. Qual = %s. Count = %s" %(qual[i], qual[i+1])
            newQualIndex += 1

    assert read.rlen == newQualIndex

    free(read.qual)
    read.qual = newQual
    read.qual[read.rlen] = 0

###################################################################################################

cdef void compressRead(cAlignedRead* read, char* refSeq, int refStart, int refEnd, int qualBinSize, int fullComp):
    """
    To save memory use, we compress reads. The sequence is compressed using reference-based
    compression. The qualities are compressed using zlib and (optionally) using course binning. We delete
    the read-group tag, and the cigar string information (this can only be done after the read has been used
    for candidate generation). A bit is set in the bit-field to label the read as compressed or uncompressed.
    """
    if Read_IsCompressed(read):
        return

    if read.seq != NULL:
        compressSeq(read, refSeq + (read.pos - refStart))

    if read.qual != NULL:
        compressQual(read, qualBinSize)

    #if fullComp == 1 and read.cigarOps != NULL:
    #    free(read.cigarOps)
    #    read.cigarOps = NULL

    if read.hash != NULL:
        free(read.hash)
        read.hash = NULL

    Read_SetCompressed(read)

###################################################################################################

cdef void uncompressRead(cAlignedRead* read, char* refSeq, int refStart, int refEnd, int qualBinSize):
    """
    Un-compress the read.
    """
    if not Read_IsCompressed(read):
        return

    if read.seq != NULL:
        uncompressSeq(read, refSeq + (read.pos - refStart))

    if read.qual != NULL:
        uncompressQual(read)

    Read_SetUnCompressed(read)

###################################################################################################


cdef class AlignedRead:
    '''
    Class representing an aligned read. see SAM format specification for meaning of
    fields (http://samtools.sourceforge.net/).

    This class stores a handle to the samtools C-structure representing
    an aligned read. Member read access is forwarded to the C-structure
    and converted into python objects. This implementation should be fast,
    as only the data needed is converted.

    For write access, the C-structure is updated in-place. This is
    not the most efficient way to build BAM entries, as the variable
    length data is concatenated and thus needs to resized if
    a field is updated. Furthermore, the BAM entry might be
    in an inconsistent state. The :meth:`~validate` method can
    be used to check if an entry is consistent.

    One issue to look out for is that the sequence should always
    be set *before* the quality scores. Setting the sequence will
    also erase any quality scores that were set previously.
    '''
    def __cinit__(self):
        self.hashValue = -1
        self.readEnd = -1

    def __dealloc__(self):
        '''
        clear up memory.
        '''
        if self._delegate:
            pysam_bam_destroy1(self._delegate)

        if self._seq:
            free(self._seq)

        if self._qual:
            free(self._qual)

    def __str__(self):
        """todo"""
        return "\t".join(map(str, (self.qname(),
                                   self.rname(),
                                   self.pos(),
                                   self.cigar(),
                                   self.seq(),
                                   self.qual(),
                                   self.mapq())))

    def __richcmp__(AlignedRead self, AlignedRead other, int opCode):
        """
        Cython's rich comparison method. All comparison operators (<, >, <=, >=, ==) are
        defined in this function.
        """
        cdef int thisChr = self.rname()
        cdef int otherChr = other.rname()
        cdef int thisPos = self.pos()
        cdef int otherPos = other.pos()
        cdef int thisLen = self.rlen()
        cdef int otherLen = other.rlen()
        cdef char* thisSeq = self.seq()
        cdef char* otherSeq = other.seq()
        cdef char* thisQual = self.qual()
        cdef char* otherQual = other.qual()

        # ==
        if opCode == 2:
            return thisChr == otherChr and thisPos == otherPos and thisLen == otherLen and thisSeq == otherSeq and thisQual == otherQual
        # <
        elif opCode == 0:
            if thisChr < otherChr:
                return True
            elif thisChr == otherChr and thisPos < otherPos:
                return True
            elif thisChr == otherChr and thisPos == otherPos and thisLen < otherLen:
                return True
            else:
                return False
        else:
            raise StandardError, "Argh. cmp code = %s" %(opCode)

    def __hash__(self):
        """
        Hash method for bam-file read object.
        """
        if self.hashValue == -1:
            self.hashValue = hash("%s_%s_%s_%s_%s_%s" %(self.pos(),self.mapq(),str(self.fastQName()),self.is_read1(),str(self.qual()),str(self.seq())))

        return self.hashValue

    cpdef object qname(AlignedRead self):
        """
        the query name (None if not present)
        """
        cdef bam1_t* src = self._delegate

        if src.core.l_qname == 0:
            return None

        return pysam_bam1_qname(src)

    cdef char* fastQName(AlignedRead self):
        """
        the query name (None if not present)
        """
        cdef bam1_t* src = self._delegate

        if src.core.l_qname == 0:
            return None

        return pysam_bam1_qname(src)

    cdef int getCigarOpCode(AlignedRead self, int index):
        """
        Return the operation code at the specified index in the
        CIGAR string.
        """
        cdef bam1_t* src = self._delegate
        cdef uint32_t* cigar_p = pysam_bam1_cigar(src)

        if index < src.core.n_cigar:
            return cigar_p[index] & BAM_CIGAR_MASK
        else:
            raise StandardError, "Index %s is >= cigar length (%s)" %(index, src.core.n_cigar)

    cdef int getCigarOpLength(AlignedRead self, int index):
        """
        Return the operation length at the specified index in the
        CIGAR string.
        """
        cdef bam1_t* src = self._delegate
        cdef uint32_t* cigar_p = pysam_bam1_cigar(src)

        if index < src.core.n_cigar:
            return cigar_p[index] >> BAM_CIGAR_SHIFT
        else:
            raise StandardError, "Index %s is >= cigar length (%s)" %(index, src.core.n_cigar)

    cdef int getCigarLength(AlignedRead self):
        """
        Return the length of the cigar string for this read.
        """
        return self._delegate.core.n_cigar

    cdef char* seq(AlignedRead self):
        """
        the query sequence (None if not present)
        """
        cdef uint8_t* p
        cdef int k = 0
        cdef char* bamTable = "=ACMGRSVTWYHKDBN"
        cdef bam1_t* src = self._delegate
        cdef int characterIndex = 0

        if not self._seq:

            # Parse qseq (bam1_seq)
            if src.core.l_qseq == 0:
                return None

            self._seq = <char*>calloc(src.core.l_qseq + 1 , sizeof(char))
            p = pysam_bam1_seq(src)

            # Equivalent to bam_nt16_rev_table[bam1_seqi(s, i)] (see bam.c)
            for k from 0 <= k < src.core.l_qseq:
                characterIndex = ((p)[(k) / 2] >> 4 * (1 - (k) % 2) & 0xf)
                self._seq[k] = bamTable[characterIndex]

        return self._seq

    cdef char* qual(AlignedRead self):
        """
        the base quality (None if not present)
        """
        cdef bam1_t* src = self._delegate
        cdef uint8_t* p

        if not self._qual:

            if src.core.l_qseq == 0:
                return None

            p = pysam_bam1_qual(src)

            if p[0] == 0xff:
                return None

            self._qual = <char*>calloc(src.core.l_qseq + 1 , sizeof(char))

            for k from 0 <= k < src.core.l_qseq:
                ## equivalent to t[i] + 33 (see bam.c)
                self._qual[k] = p[k] + 33

        return self._qual

    cdef dict tags(AlignedRead self):
        """
        the tags in the AUX field.
        """
        cdef char* ctag
        cdef bam1_t* src
        cdef uint8_t* s
        cdef char tpe
        cdef dict result = dict()

        src = self._delegate

        if src.l_aux == 0:
            return None

        s = pysam_bam1_aux( src )
        ctag = <char*>calloc( 3, sizeof(char) )
        cdef int x

        while s < (src.data + src.data_len):
            # get tag
            ctag[0] = s[0]
            ctag[1] = s[1]
            ctag[2] = 0
            pytag = ctag
            s += 2

            # get type and value
            # how do I do char literal comparison in cython?
            # the code below works (i.e, is C comparison)
            tpe = toupper(s[0])
            if tpe == <char>'S':
                value = <int>bam_aux2i(s)
                s += 2
            elif tpe == <char>'I':
                value = <int>bam_aux2i(s)
                s += 4
            elif tpe == <char>'F':
                value = <float>bam_aux2f(s)
                s += 4
            elif tpe == <char>'D':
                value = <double>bam_aux2d(s)
                s += 8
            elif tpe == <char>'C':
                value = <int>bam_aux2i(s)
                s += 1
            elif tpe == <char>'A':
                # there might a more efficient way
                # to convert a char into a string
                # yes there is:
                value = <char>bam_aux2A(s)
                s += 1
            elif tpe == <char>'Z' or tpe == <char>'H':
                value = <char*>bam_aux2Z(s)
                # +1 for NULL terminated string
                s += len(value) + 1
            else:
                print "Unrecognized type ",tpe
            s += 1

            # ignore type
            result[pytag] = value

        free( ctag )
        return result

    cdef int flag(AlignedRead self):
        """
        properties flag
        """
        return self._delegate.core.flag

    cdef int rname(AlignedRead self):
        """
        :term:`target` ID

        .. note::

            This field contains the index of the reference sequence
            in the sequence dictionary. To obtain the name
            of the reference sequence, use :meth:`pysam.Samfile.getrname()`

        """
        return self._delegate.core.tid

    cdef int pos(AlignedRead self):
        """
        0-based leftmost coordinate
        """
        return self._delegate.core.pos

    cdef int end(AlignedRead self):
        """
        0-based right-tmost coordinate, including
        cigar shifts (deletions, skips etc).
        """
        cdef int cigarLength = 0
        cdef int cigarIndex = 0
        cdef int flag = 0
        cdef int length = 0

        if self.readEnd == -1:
            self.readEnd = self.pos()
            cigarLength = self.getCigarLength()

            for cigarIndex from 0 <= cigarIndex < cigarLength:
                flag = self.getCigarOpCode(cigarIndex)
                length = self.getCigarOpLength(cigarIndex)

                if flag == 0 or flag == 2 or flag == 3:
                    self.readEnd += length

        return self.readEnd

    cdef bin(AlignedRead self):
        """
        properties bin
        """
        return self._delegate.core.bin

    cdef int rlen(AlignedRead self):
        '''
        length of the read (read only). Returns 0 if not given.
        '''
        return self._delegate.core.l_qseq

    cpdef int mapq(AlignedRead self):
        """
        mapping quality
        """
        return self._delegate.core.qual

    cdef int mrnm(AlignedRead self):
        """
        the :term:`reference` id of the mate
        """
        return self._delegate.core.mtid

    cdef int mpos(AlignedRead self):
        """
        the position of the mate
        """
        return self._delegate.core.mpos

    cdef int isize(AlignedRead self):
        """
        the insert size
        """
        return self._delegate.core.isize

    cdef int is_paired(AlignedRead self):
        """
        true if read is paired in sequencing
        """
        return (self._delegate.core.flag & BAM_FPAIRED) != 0

    cdef int is_proper_pair(AlignedRead self):
        """
        true if read is mapped in a proper pair
        """
        return (self.flag() & BAM_FPROPER_PAIR) != 0

    cdef int is_unmapped(AlignedRead self):
        """
        true if read itself is unmapped
        """
        return (self.flag() & BAM_FUNMAP) != 0

    cdef int mate_is_unmapped(AlignedRead self):
        """
        true if the mate is unmapped
        """
        return (self.flag() & BAM_FMUNMAP) != 0

    cdef int is_reverse(AlignedRead self):
        """
        true if read is mapped to reverse strand
        """
        return (self.flag() & BAM_FREVERSE) != 0

    cdef int mate_is_reverse(AlignedRead self):
        """
        true is read is mapped to reverse strand
        """
        return (self.flag() & BAM_FMREVERSE) != 0

    cdef int is_read1(AlignedRead self):
        """
        true if this is read1
        """
        return (self.flag() & BAM_FREAD1) != 0

    cdef int is_read2(AlignedRead self):
        """
        true if this is read2
        """
        return (self.flag() & BAM_FREAD2) != 0

    cdef int is_secondary(AlignedRead self):
        """
        true if not primary alignment
        """
        return (self.flag() & BAM_FSECONDARY) != 0

    cdef int is_qcfail(AlignedRead self):
        """
        true if QC failure
        """
        return (self.flag() & BAM_FQCFAIL) != 0

    cdef int is_duplicate(AlignedRead self):
        """
        true if optical or PCR duplicate
        """
        return (self.flag() & BAM_FDUP) != 0

    cdef opt(AlignedRead self, tag):
        """
        retrieves optional data given a two-letter *tag*
        """
        #see bam_aux.c: bam_aux_get() and bam_aux2i() etc
        cdef uint8_t * v
        v = bam_aux_get(self._delegate, tag)

        if v == NULL:
            raise KeyError( "tag '%s' not present" % tag )

        type = chr(v[0])

        if type == 'c' or type == 'C' or type == 's' or type == 'S' or type == 'i':
            return <int>bam_aux2i(v)
        elif type == 'f':
            return <float>bam_aux2f(v)
        elif type == 'd':
            return <double>bam_aux2d(v)
        elif type == 'A':
            # there might a more efficient way
            # to convert a char into a string
            return '%c' % <char>bam_aux2A(v)
        elif type == 'Z':
            return <char*>bam_aux2Z(v)

###################################################################################################

cdef AlignedRead makeAlignedRead(bam1_t * src):
    '''
    enter src into AlignedRead.
    '''
    cdef AlignedRead dest = AlignedRead()
    dest._delegate = bam_dup1(src)
    return dest

#####################################################################
