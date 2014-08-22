import cython

cdef extern from "string.h":
  ctypedef int size_t
  void *memcpy(void *dst,void *src,size_t len)
  void *memmove(void *dst,void *src,size_t len)
  void *memset(void *b,int c,size_t len)

cdef extern from "stdlib.h":
  void free(void *)
  void *malloc(size_t)
  void *calloc(size_t,size_t)
  void *realloc(void *,size_t)
  int c_abs "abs" (int)
  void qsort(void *base, size_t nmemb, size_t size,
             int (*compar)(void *,void *))

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
  int tolower(int c)
  
cdef extern from "unistd.h":
  char *ttyname(int fd)
  int isatty(int fd)  

cdef extern from "string.h":
  int strcmp(char *s1, char *s2)
  int strncmp(char *s1,char *s2,size_t len)
  char *strcpy(char *dest,char *src)
  char *strncpy(char *dest,char *src, size_t len)
  char *strdup(char *)
  char *strcat(char *,char *)
  size_t strlen(char *s)
  int memcmp( void * s1, void *s2, size_t len )

cdef extern from "razf.h":
  pass

cdef extern from "stdint.h":
  ctypedef int int64_t
  ctypedef int int32_t
  ctypedef int uint32_t
  ctypedef int uint8_t
  ctypedef int uint64_t


cdef extern from "bam.h":

  # IF _IOLIB=2, bamFile = BGZF, see bgzf.h
  # samtools uses KNETFILE, check how this works

    ctypedef struct tamFile:
        pass

    ctypedef struct bamFile:
        void* uncompressed_block
        void* compressed_block

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
    
    ctypedef int (*bam_fetch_f)(bam1_t *b, void *data)
    
    ctypedef struct bam_header_t:
       int32_t n_targets
       char **target_name
       uint32_t *target_len
       void *hash
       void *rg2lib
       int l_text
       char *text
    
    ctypedef struct bam_index_t:
        pass
    
    bamFile* razf_dopen(int data_fd, char *mode)
    
    # removed - macros not found
    
    # int64_t bam_seek( bamFile fp, uint64_t voffset, int where)
    # int64_t bam_tell( bamFile fp )
    # void bam_destroy1( bam1_t * b) 
    # void bam_init_header_hash(bam_header_t *header)
    
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
    
    bam1_t * bam_dup1( bam1_t *src ) 
    
    bam1_t * bam_copy1(bam1_t *bdst, bam1_t *bsrc)
    
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
    uint32_t * pysam_bam1_cigar( bam1_t * b)
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
        bamFile*	fp
        int	tid
        int	beg
        int	end

    bam_fetch_iterator_t* bam_init_fetch_iterator(bamFile* fp, bam_index_t *idx, int tid, int beg, int end)
    bam1_t * bam_fetch_iterate(bam_fetch_iterator_t *iter) nogil
    void bam_cleanup_fetch_iterator(bam_fetch_iterator_t *iter)

###################################################################################################
cdef class Samfile

@cython.final
cdef class IteratorRow:
    cdef bam_fetch_iterator_t*  bam_iter # iterator state object
    cdef bam1_t* b
    cdef int useIndex
    cdef int cnext(self) nogil

###################################################################################################

@cython.final
cdef class Samfile:
    cdef char* filename
    cdef samfile_t* samfile
    cdef bam_index_t *index
    cdef int isbam
    cdef void clearHeader(self)
    cdef void clearIndex(self)
    cdef void createCache(self)
    cdef void clearCache(self)
    cdef dict indexesByRegion
    cdef object lock
    cdef int _isOpen( self )
    cdef _hasIndex( self )
    cpdef _open(self, mode, loadIndex=*)
    cdef char* getrname(self, int tid)
    cdef _parseRegion(self, reference=*, start=*, end=*, region=*)
    cpdef IteratorRow fetch(self, char* reference, int start, int end)
    cpdef close(self)
    cdef void loadOffsetsForRegions(self, list regions)
    #cdef void cleanUpOffsetsForRegions(self)
    cdef void openBAMFile(self, mode)

###################################################################################################

ctypedef struct cAlignedRead:
    char* seq
    char* qual
    short* cigarOps
    short* hash
    short mateChromID
    short cigarLen
    short chromID
    short rlen
    int pos
    int end
    int insertSize
    int matePos
    int bitFlag
    char mapq

###################################################################################################

cdef cAlignedRead* createRead(bam1_t * src, int storeRgID, char** rgID)
cdef void destroyRead(cAlignedRead* theRead)

###################################################################################################

# These are bits set in the BAM bit-flag.
DEF BAM_FPAIRED = 1        # Read is paired in sequencing, no matter whether it is mapped in a pair
DEF BAM_FPROPER_PAIR = 2   # Read is mapped in a proper pair
DEF BAM_FUNMAP = 4         # Read itself is unmapped; conflictive with BAM_FPROPER_PAIR
DEF BAM_FMUNMAP = 8        # Mate is unmapped
DEF BAM_FREVERSE = 16      # Read is mapped to the reverse strand
DEF BAM_FMREVERSE = 32     # Mate is mapped to the reverse strand
DEF BAM_FREAD1 = 64        # This is read1
DEF BAM_FREAD2 = 128       # This is read2
DEF BAM_FSECONDARY = 256   # Not primary alignment
DEF BAM_FQCFAIL = 512      # QC failure
DEF BAM_FDUP = 1024        # Optical or PCR duplicate
DEF BAM_FCOMPRESSED = 2048 # Is the read compressed

###################################################################################################

# And here are accessor functions for the bit-fields
@cython.profile(False)
cdef inline int Read_IsReverse(cAlignedRead* theRead) nogil:
    return ( (theRead.bitFlag & BAM_FREVERSE) != 0)

@cython.profile(False)
cdef inline int Read_IsPaired(cAlignedRead* theRead) nogil:
    return ( (theRead.bitFlag & BAM_FPAIRED) != 0)

@cython.profile(False)
cdef inline int Read_IsProperPair(cAlignedRead* theRead) nogil:
    return ( (theRead.bitFlag & BAM_FPROPER_PAIR) != 0)

@cython.profile(False)
cdef inline int Read_IsDuplicate(cAlignedRead* theRead) nogil:
    return ( (theRead.bitFlag & BAM_FDUP) != 0)

@cython.profile(False)
cdef inline int Read_IsUnmapped(cAlignedRead* theRead) nogil:
    return ( (theRead.bitFlag & BAM_FUNMAP) != 0)

@cython.profile(False)
cdef inline int Read_MateIsUnmapped(cAlignedRead* theRead) nogil:
    return ( (theRead.bitFlag & BAM_FMUNMAP) != 0)

@cython.profile(False)
cdef inline int Read_MateIsReverse(cAlignedRead* theRead) nogil:
    return ( (theRead.bitFlag & BAM_FMREVERSE) != 0)

@cython.profile(False)
cdef inline int Read_IsQCFail(cAlignedRead* theRead) nogil:
    return ( (theRead.bitFlag & BAM_FQCFAIL) != 0)

@cython.profile(False)
cdef inline int Read_IsReadOne(cAlignedRead* theRead) nogil:
    return ( (theRead.bitFlag & BAM_FREAD1) != 0)

@cython.profile(False)
cdef inline int Read_IsSecondaryAlignment(cAlignedRead* theRead) nogil:
    return ( (theRead.bitFlag & BAM_FSECONDARY) != 0)

@cython.profile(False)
cdef inline int Read_IsCompressed(cAlignedRead* theRead) nogil:
    return ( (theRead.bitFlag & BAM_FCOMPRESSED) != 0)

@cython.profile(False)
cdef inline int Read_SetIsNotReverse(cAlignedRead* theRead) nogil:
    theRead.bitFlag &= (~BAM_FREVERSE)

@cython.profile(False)
cdef inline int Read_SetIsReverse(cAlignedRead* theRead) nogil:
    theRead.bitFlag |= BAM_FREVERSE

@cython.profile(False)
cdef inline void Read_SetQCFail(cAlignedRead* theRead) nogil:
    theRead.bitFlag |= BAM_FQCFAIL

@cython.profile(False)
cdef inline void Read_SetCompressed(cAlignedRead* theRead) nogil:
    theRead.bitFlag |= BAM_FCOMPRESSED

@cython.profile(False)
cdef inline void Read_SetUnCompressed(cAlignedRead* theRead) nogil:
    theRead.bitFlag &= (~BAM_FCOMPRESSED)

###################################################################################################

cdef void compressRead(cAlignedRead* read, char* refSeq, int refStart, int refEnd, int qualBinSize, int fullComp)
cdef void uncompressRead(cAlignedRead* read, char* refSeq, int refStart, int refEnd, int qualBinSize)
