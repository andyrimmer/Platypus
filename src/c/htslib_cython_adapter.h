#ifndef CRAM_htslib_cython_adapter_h
#define CRAM_htslib_cython_adapter_h

#include <stdint.h>

#include <htslib/sam.h>

htsFile *sam_open_(const char *fn, const char *mode);

hts_idx_t *sam_index_load_(htsFile *fp, const char *fn);

bam_hdr_t *bam_hdr_init_(void);

bam_hdr_t *sam_hdr_read_(samFile *fp);

bam1_t *bam_init1_(void);

int sam_read1_(samFile *fp, bam_hdr_t *h, bam1_t *b);

uint8_t *bam_get_seq_(const bam1_t *b);

hts_itr_t *sam_itr_querys_(const hts_idx_t *idx, bam_hdr_t *hdr, const char *region);

int sam_itr_next_(htsFile *htsfp, hts_itr_t *itr, bam1_t *r);

int sam_close_(htsFile *fp);

uint8_t bam_seqi_(uint8_t *s, int i);

#endif
