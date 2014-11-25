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

from cpython cimport bool

cdef extern from "stdlib.h":
  void free(void *)
  void *alloca(size_t)
  void *calloc(size_t,size_t)
  void *malloc(size_t)
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
  int memcmp( void * s1, void *s2, size_t len)

cdef extern from "stdint.h":
    ctypedef int int8_t
    ctypedef int int64_t
    ctypedef int int32_t
    ctypedef int uint32_t
    ctypedef int uint8_t
    ctypedef int uint64_t

cdef extern from "htslib/bgzf.h":
    
    ctypedef struct z_stream:
        pass
    
    ctypedef struct hFILE:
        pass
    
    ctypedef struct bgzf_mtaux_t:
        pass
    
    ctypedef struct bgzidx_t:
        pass
    
    ctypedef struct BGZF:
        int errcode
        int is_write
        int is_be
        int compress_level
        int is_compressed
        int is_gzip
        int cache_size
        int block_length
        int block_offset
        int64_t block_address
        int64_t uncompressed_address
        void *uncompressed_block
        void *compressed_block
        void *cache
        hFILE *fp
        bgzf_mtaux_t *mt
        bgzidx_t *idx
        int idx_build_otf
        z_stream *gz_stream

cdef extern from "htslib/hts.h":
    
    ctypedef struct kstring_t:
        size_t l
        size_t m
        char *s
    
    ctypedef struct hts_idx_t:
        pass
    
    ctypedef struct hts_itr_t:
        pass
    
    void hts_idx_destroy(hts_idx_t *idx)
    const char *hts_parse_reg(const char *s, int *beg, int *end)

cdef extern from "htslib/sam.h":
    
    ctypedef struct bam_hdr_t:
        int32_t n_targets
        int32_t ignore_sam_err
        uint32_t l_text
        uint32_t *target_len
        int8_t *cigar_tab
        char **target_name
        char *text
        void *sdict

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
        int l_data
        int m_data
        uint8_t *data
        uint64_t id
    
    ctypedef struct cram_fd:
        pass
    
    ctypedef union samFileUnion:
        BGZF *bgzf
        cram_fd *cram
        hFILE *hfile
        void *voidp
    
    ctypedef struct samFile:
        uint32_t is_bin
        uint32_t is_write
        uint32_t is_be
        uint32_t is_cram
        uint32_t is_compressed
        uint32_t is_kstream
        uint32_t dummy
        int64_t lineno
        kstring_t line
        char *fn
        char *fn_aux
        samFileUnion fp
    
    samFile *sam_open(const char *fn, const char *mode)
    int bam_name2id(bam_hdr_t *h, const char *ref)
    hts_idx_t *sam_index_load(samFile *fp, const char *fn)
    bam_hdr_t *bam_hdr_init()
    bam_hdr_t *sam_hdr_read(samFile *fp)
    bam1_t *bam_init1()
    int sam_read1(samFile *fp, bam_hdr_t *h, bam1_t *b)
    bool bam_is_rev(const bam1_t *b)
    bool bam_is_mrev(const bam1_t *b)
    char *bam_get_qname(const bam1_t *b)
    uint8_t *bam_get_seq(const bam1_t *b)
    uint8_t *bam_get_qual(const bam1_t *b)
    uint32_t * bam_get_cigar(const bam1_t *b)
    hts_itr_t *sam_itr_queryi(const hts_idx_t *idx, int tid, int beg, int end)
    hts_itr_t *sam_itr_querys(const hts_idx_t *idx, bam_hdr_t *hdr, const char *region)
    int sam_itr_next(samFile *htsfp, hts_itr_t *itr, bam1_t *r)
    int sam_close(samFile *fp)
    uint8_t bam_seqi(uint8_t *s, int i)
    void bam_hdr_destroy(bam_hdr_t *h)
    void sam_itr_destroy(hts_idx_t *idx)
