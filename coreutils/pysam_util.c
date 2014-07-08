#include "pysam_util.h"
#include "khash.h"
#include "ksort.h"
#include "bam_endian.h"
#include "knetfile.h"

// #######################################################
// utility routines to avoid using callbacks in bam_fetch
// taken from bam_index.c
// The order of the following declarations is important.
// #######################################################

//-------------------------------------------------------------------------------------------------

#define pair64_lt(a,b) ((a).u < (b).u)
KSORT_INIT(my_off, pair64_t, pair64_lt);

//-------------------------------------------------------------------------------------------------

typedef struct 
{
	uint32_t m, n;
	pair64_t *list;
} bam_binlist_t;

//-------------------------------------------------------------------------------------------------

typedef struct 
{
	int32_t n, m;
	uint64_t *offset;
} bam_lidx_t;

//-------------------------------------------------------------------------------------------------

KHASH_MAP_INIT_INT(my_i, bam_binlist_t);

//-------------------------------------------------------------------------------------------------

struct __bam_index_t
{
  int32_t n;
  khash_t(my_i) **index;
  bam_lidx_t *index2;
};

//-------------------------------------------------------------------------------------------------

typedef struct __linkbuf_t 
{
	bam1_t b;
	uint32_t beg, end;
	struct __linkbuf_t *next;
} lbnode_t;

//-------------------------------------------------------------------------------------------------

// standin for bam_destroy1 in bam.h
// deletes all variable length data
void pysam_bam_destroy1(bam1_t* b)
{
  if (b == NULL) return;
  if (b->data != NULL) free(b->data);
  free(b);
}

//-------------------------------------------------------------------------------------------------

unsigned char pysam_translate_sequence( const unsigned char s )
{
    // translate a nucleotide character to binary code
  return bam_nt16_table[s];
}

//-------------------------------------------------------------------------------------------------

char * pysam_bam1_qname( const bam1_t * b)
{
    // stand-ins for samtools macros in bam.h
  return (char*)b->data;
}

//-------------------------------------------------------------------------------------------------

uint32_t * pysam_bam1_cigar( const bam1_t * b) 
{
  return (uint32_t*)(b->data + b->core.l_qname);
}

//-------------------------------------------------------------------------------------------------

uint8_t * pysam_bam1_seq( const bam1_t * b) 
{
  return (uint8_t*)(b->data + b->core.n_cigar*4 + b->core.l_qname);
}

//-------------------------------------------------------------------------------------------------

uint8_t * pysam_bam1_qual( const bam1_t * b)
{
  return (uint8_t*)(b->data + b->core.n_cigar*4 + b->core.l_qname + (b->core.l_qseq + 1)/2);
}

//-------------------------------------------------------------------------------------------------

uint8_t * pysam_bam1_aux( const bam1_t * b)
{
  return (uint8_t*)(b->data + b->core.n_cigar*4 + b->core.l_qname + b->core.l_qseq + (b->core.l_qseq + 1)/2);
}

//-------------------------------------------------------------------------------------------------
//
//
//
// Iterator implementation follows
//
//
// functions defined in bam_index.c

extern pair64_t* get_chunk_coordinates(const bam_index_t *idx, int tid, int beg, int end, int* cnt_off);

//-------------------------------------------------------------------------------------------------

static inline int is_overlap(uint32_t beg, uint32_t end, const bam1_t *b)
{
	uint32_t rbeg = b->core.pos;
	uint32_t rend = b->core.n_cigar? bam_calend(&b->core, bam1_cigar(b)) : b->core.pos + 1;
	return (rend > beg && rbeg < end);
}

//-------------------------------------------------------------------------------------------------

bam_fetch_iterator_t* bam_init_fetch_iterator(bamFile* fp, const bam_index_t *idx, int tid, int beg, int end)
{
	// iterator contains current alignment position
	//      and will contain actual alignment during iterations
	bam_fetch_iterator_t* iter  = (bam_fetch_iterator_t*)calloc(1, sizeof(bam_fetch_iterator_t));
	iter->b                     = (bam1_t*)calloc(1, sizeof(bam1_t));
		
	// list of chunks containing our alignments
	iter->off = get_chunk_coordinates(idx, tid, beg, end, &iter->n_off);
	
	// initialise other state variables in iterator
	iter->fp                = fp;
	iter->curr_chunk        = -1;   
	iter->curr_off          =  0;
	iter->n_seeks           =  0;    
	iter->tid				= tid;
	iter->beg				= beg;
	iter->end				= end;
	return iter;
}

//-------------------------------------------------------------------------------------------------

bam1_t* bam_fetch_iterate(bam_fetch_iterator_t *iter)
{
	if (!iter->off) {
		return 0;
	}

	int ret;
	// iterate through all alignments in chunks
	for (;;) {
		if (iter->curr_off == 0 || iter->curr_off >= iter->off[iter->curr_chunk].v) { // then jump to the next chunk
			if (iter->curr_chunk == iter->n_off - 1) break; // no more chunks
			if (iter->curr_chunk >= 0) assert(iter->curr_off == iter->off[iter->curr_chunk].v); // otherwise bug
			if (iter->curr_chunk < 0 || iter->off[iter->curr_chunk].v != iter->off[iter->curr_chunk+1].u) { // not adjacent chunks; then seek
				bam_seek(iter->fp, iter->off[iter->curr_chunk+1].u, SEEK_SET);
				iter->curr_off = bam_tell(iter->fp);
				++iter->n_seeks;
			}
			++iter->curr_chunk;
		}
		if ((ret = bam_read1(iter->fp, iter->b)) > 0) {
			iter->curr_off = bam_tell(iter->fp);
			if (iter->b->core.tid != iter->tid || iter->b->core.pos >= iter->end) break; // no need to proceed
			else if (is_overlap(iter->beg, iter->end, iter->b)) 
				//
				//func(iter->b, data);
				//
				return iter->b;
		} else 
			return 0; // end of file
	}
	return 0;
}

//-------------------------------------------------------------------------------------------------

void bam_cleanup_fetch_iterator(bam_fetch_iterator_t *iter)
{
  //  fprintf(stderr, "[bam_fetch] # seek calls: %d\n", iter->n_seeks);
  bam_destroy1(iter->b);
  free(iter->off);
}

//-------------------------------------------------------------------------------------------------
