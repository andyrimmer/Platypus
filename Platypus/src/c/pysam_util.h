#ifndef PYSAM_UTIL_H
#define PYSAM_UTIL_H

#include <ctype.h>
#include <assert.h>
#include "bam.h"

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
// code for iterator

/*! @typedef
  @Structure for holding current state (current alignment etc.) for iterating through
  alignments overlapping a specified region.
  @field  b           pointer to the current alignment
  @field  off         pointer to an array of chunk loci (each with beg/end positions)
  @field  n_off       The number of chunks
  @field  curr_off    The current file positon
  @field  curr_chunk  The item in a list of chunk
  @discussion See also bam_fetch_iterate
*/

//-------------------------------------------------------------------------------------------------

typedef struct
{
  uint64_t u, v;
} pair64_t;

//-------------------------------------------------------------------------------------------------

typedef struct 
{
  bam1_t *        b;
  pair64_t *      off;
  int             n_off;
  uint64_t        curr_off;
  int             curr_chunk;
  bamFile* 		fp;
  int				tid;
  int				beg;
  int				end;
  int             n_seeks;
} bam_fetch_iterator_t;

//-------------------------------------------------------------------------------------------------



/*!
  @abstract Retrieve the alignments that are overlapped with the
  specified region.

  @discussion Returns iterator object to retrieve successive alignments ordered by
  start position.
  @param  fp    BAM file handler
  @param  idx   pointer to the alignment index
  @param  tid   chromosome ID as is defined in the header
  @param  beg   start coordinate, 0-based
  @param  end   end coordinate, 0-based
*/
bam_fetch_iterator_t* bam_init_fetch_iterator(bamFile* fp, const bam_index_t *idx, int tid, int beg, int end);
void bam_cleanup_fetch_iterator(bam_fetch_iterator_t *iter);


/*!
  @abstract Iterates through alignments overlapped the specified region.
  @discussion Returns pointer to successive alignments ordered by start position.
  Returns null pointer to signal the end of the iteration.
  The alignment data is nested within the iterator to avoid unnecessary allocations.
*/
bam1_t* bam_fetch_iterate(bam_fetch_iterator_t *iter);
bam_fetch_iterator_t* bam_init_fetchall_iterator(bamFile* fp, const bam_index_t *idx);
bam1_t* bam_fetchall_iterate(bam_fetch_iterator_t *iter);

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
// various helper functions

// stand-in for macro - not wrappable in pyrex
void pysam_bam_destroy1(bam1_t* b);

// stand-in for other samtools macros
uint32_t* pysam_bam1_cigar(const bam1_t* b);
char* pysam_bam1_qname(const bam1_t* b);
uint8_t* pysam_bam1_seq(const bam1_t* b);
uint8_t* pysam_bam1_qual(const bam1_t* b);
uint8_t* pysam_bam1_aux(const bam1_t* b);

// translate a nucleotide character to binary code
unsigned char pysam_translate_sequence(const unsigned char s);

#endif
