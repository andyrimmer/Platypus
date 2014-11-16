#include "htslib_cython_adapter.h"

htsFile *sam_open_(const char *fn, const char *mode)
{
    return sam_open(fn, mode);
}

hts_idx_t *sam_index_load_(htsFile *fp, const char *fn)
{
    return sam_index_load(fp, fn);
}

bam_hdr_t *bam_hdr_init_(void)
{
    return bam_hdr_init();
}

bam_hdr_t *sam_hdr_read_(samFile *fp)
{
    return sam_hdr_read(fp);
}

bam1_t *bam_init1_(void)
{
    return bam_init1();
}

int sam_read1_(samFile *fp, bam_hdr_t *h, bam1_t *b)
{
    return sam_read1(fp, h, b);
}

uint8_t *bam_get_seq_(const bam1_t *b)
{
    return bam_get_seq(b);
}

hts_itr_t *sam_itr_querys_(const hts_idx_t *idx, bam_hdr_t *hdr, const char *region)
{
    return sam_itr_querys(idx, hdr, region);
}

int sam_itr_next_(htsFile *htsfp, hts_itr_t *itr, bam1_t *r)
{
    return sam_itr_next(htsfp, itr, r);
}

int sam_close_(htsFile *fp)
{
    return sam_close(fp);
}

uint8_t bam_seqi_(uint8_t *s, int i)
{
    return bam_seqi(s, i);
}
