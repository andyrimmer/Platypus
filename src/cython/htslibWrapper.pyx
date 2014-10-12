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
    
    char *bam_get_qname(bam1_t *b)
    uint32_t *bam_get_cigar(bam1_t *b)
    uint8_t *bam_get_seq(bam1_t *b)
    uint8_t *bam_get_qual(bam1_t *b)
    uint8_t *bam_get_aux(bam1_t *b)
    
    bam_hdr_t *bam_hdr_init(void)
    bam1_t *bam_init1(void)
    
    hts_idx_t *sam_index_load(htsFile *fp, const char *fn)
    
    void sam_itr_destroy(hts_idx_t *idx)
    hts_itr_t *sam_itr_queryi(const hts_idx_t *idx, int tid, int beg, int end)
    hts_itr_t *sam_itr_querys(const hts_idx_t *idx, bam_hdr_t *hdr, const char *region)
    int sam_itr_next(void *data, hts_itr_t *iter, void *r)
    
    htsFile sam_open(const char *fn, const char *mode)
    int sam_close(htsFile *fp)
    
    ctypedef htsFile samFile
    
    bam_hdr_t *sam_hdr_read(samFile *fp)
    
    int sam_read1(samFile *fp, bam_hdr_t *h, bam1_t *b)
    int sam_write1(samFile *fp, const bam_hdr_t *h, const bam1_t *b)
    
