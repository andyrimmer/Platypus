#ifndef ALIGN_H
#define ALIGN_H

/*****************************************************************************************************************
 This code is copyright (c) Gerton Lunter, Jan 2009
 It may not be distributed, made public, or used in other software without the permission of the copyright holder
******************************************************************************************************************/

#include <emmintrin.h>

static inline short extract0(__m128i x) {return _mm_extract_epi16(x, 0);}
static inline short extract1(__m128i x) {return _mm_extract_epi16(x, 1);}
static inline short extract2(__m128i x) {return _mm_extract_epi16(x, 2);}
static inline short extract3(__m128i x) {return _mm_extract_epi16(x, 3);}
static inline short extract4(__m128i x) {return _mm_extract_epi16(x, 4);}
static inline short extract5(__m128i x) {return _mm_extract_epi16(x, 5);}
static inline short extract6(__m128i x) {return _mm_extract_epi16(x, 6);}
static inline short extract7(__m128i x) {return _mm_extract_epi16(x, 7);}

int fastAlignmentRoutine(char* seq1, char* seq2, char* qual2, int len1, int len2, int gapextend, int nucprior, short* localgapopen);

#endif
