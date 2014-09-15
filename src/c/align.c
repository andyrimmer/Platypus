#include "align.h"

/*****************************************************************************************************************
 This code is copyright (c) Gerton Lunter, Dec 2009
 It may not be distributed, made public, or used in other software without the permission of the copyright holder
******************************************************************************************************************/

//#include <emmintrin.h>
#include <stdio.h>
#include <assert.h>

//_________________________________________________________________________________________________

//inline short extract0(__m128i x) {return _mm_extract_epi16(x, 0);}
//inline short extract1(__m128i x) {return _mm_extract_epi16(x, 1);}
//inline short extract2(__m128i x) {return _mm_extract_epi16(x, 2);}
//inline short extract3(__m128i x) {return _mm_extract_epi16(x, 3);}
//inline short extract4(__m128i x) {return _mm_extract_epi16(x, 4);}
//inline short extract5(__m128i x) {return _mm_extract_epi16(x, 5);}
//inline short extract6(__m128i x) {return _mm_extract_epi16(x, 6);}
//inline short extract7(__m128i x) {return _mm_extract_epi16(x, 7);}
//
short (*extractors[8]) (__m128i x) = {extract0,extract1,extract2,extract3,extract4,extract5,extract6,extract7};

//_________________________________________________________________________________________________


int fastAlignmentRoutine(char* seq1, char* seq2, char* qual2, int len1, int len2, int gapextend, int nucprior,
			 short* localgapopen) {

  // seq2 is the read; the shorter of the sequences
  // no checks for overflow are done

  // the bottom-left and top-right corners of the DP table are just
  // included at the extreme ends of the diagonal, which measures
  // n=8 entries diagonally across.  This fixes the length of the
  // longer (horizontal) sequence to 14 (2*8-2) more than the shorter

  assert( len1 == len2 + 2*8 - 1);

  // make sure that special cases at beginning and end (initialization!)
  // do not mix.  (Todo: Check this assertion.)
  assert( len1 > 8);

  const short gap_extend = gapextend*4;
  const short nuc_prior = nucprior*4;
  const short n_score = 0*4;
  const short pos_inf = 0x7800;
  const int match_label = 0;
  const int insert_label = 1;
  const int delete_label = 3;

  register __m128i _m1;
  register __m128i _i1;
  register __m128i _d1;
  register __m128i _m2;
  register __m128i _i2;
  register __m128i _d2;
  __m128i _seq1win;
  __m128i _seq2win;
  __m128i _qual2win;
  __m128i _seq1nqual;   // 1 if N, pos_inf if not

  __m128i _gap_extend = _mm_set1_epi16( gap_extend );
  __m128i _nuc_prior = _mm_set1_epi16( nuc_prior );
  __m128i _three = _mm_set1_epi16( 3 );
  __m128i _initmask = _mm_set_epi16( 0,0,0,0,0,0,0,-1 );
  __m128i _initmask2 = _mm_set_epi16( 0,0,0,0,0,0,0,-0x8000 );
  __m128i _backpointers[ 2*(len1+8) ];

  // initialization
  _m1 = _mm_set1_epi16( pos_inf );
  _i1 = _m1;
  _d1 = _m1;
  _m2 = _m1;
  _i2 = _m1;
  _d2 = _m1;

  // sequence 1 is initialized with the n-long prefix, in forward direction
  // sequence 2 is initialized as empty; reverse direction
  _seq1win = _mm_set_epi16( seq1[7], seq1[6], seq1[5], seq1[4], seq1[3], seq1[2], seq1[1], seq1[0] );
  _seq2win = _m1;
  _qual2win = _mm_set1_epi16(64*4);

  // if N, make n_score; if != N, make pos_inf
  _seq1nqual = _mm_add_epi16( _mm_and_si128( _mm_cmpeq_epi16( _seq1win, 
							      _mm_set1_epi16( 'N' ) ), 
					     _mm_set1_epi16( n_score - pos_inf ) ), 
			      _mm_set1_epi16( pos_inf ) );


  __m128i _gap_open = _mm_set_epi16(localgapopen[7],localgapopen[6],localgapopen[5],localgapopen[4],
				    localgapopen[3],localgapopen[2],localgapopen[1],localgapopen[0]);

  short _score = 0;
  short minscore = pos_inf;
  short minscoreidx = -1;

  // main loop.  Do one extra iteration, with nucs from sequence 2 just moved out
  // of the seq2win/qual arrays, to simplify getting back pointers
  int s = 0;
  for (; s<2*(len2+8-1 + 1); s+=2) {

    // seq1 is current; seq2 needs updating
    _seq2win = _mm_slli_si128( _seq2win, 2 );
    _qual2win = _mm_slli_si128( _qual2win, 2 );
    if (s/2 < len2) {
      _seq2win = _mm_insert_epi16( _seq2win, seq2[ s/2 ], 0 );
      _qual2win = _mm_insert_epi16( _qual2win, 4*qual2[ s/2 ], 0 );
    } else {
      _seq2win = _mm_insert_epi16( _seq2win, '0', 0 );
      _qual2win = _mm_insert_epi16( _qual2win, 64*4, 0 );
    }

    //
    // S even
    //

    // initialize to -0x8000
    _m1 = _mm_or_si128( _initmask2, _mm_andnot_si128( _initmask, _m1 ) );
    // moved here, from below, to allow alignments to start with insertions or deletions
    _m2 = _mm_or_si128( _initmask2, _mm_andnot_si128( _initmask, _m2 ) );
    _m1 = _mm_min_epi16( _m1, _mm_min_epi16( _i1, _d1 ) );

    // at this point, extract minimum score.  Referred-to position must
    // be y==len2-1, so that current position has y==len2; i==0 so d=0 and y=s/2

    if (s/2 >= len2) {
        _score = (*extractors[s/2 - len2])(_m1);

      if (_score < minscore) {
          minscore = _score;
          minscoreidx = s;     // point back to the match state at this entry, so as not to
      }                      // have to store the state at s-2
    }

    _m1 = _mm_add_epi16( _m1, 
			 _mm_min_epi16( _mm_andnot_si128( _mm_cmpeq_epi16( _seq2win, 
									   _seq1win ), 
							  _qual2win ),
					_seq1nqual ) );
    _d1 = _mm_min_epi16( _mm_add_epi16( _d2,
					_gap_extend ),
			 _mm_add_epi16( _mm_min_epi16( _m2,
						       _i2 ),   // allow I->D
					_gap_open ) );

    _d1 = _mm_insert_epi16( _mm_slli_si128( _d1, 
					    2 ), 
			    pos_inf, 
			    0 );

    _i1 = _mm_add_epi16( _mm_min_epi16( _mm_add_epi16( _i2,
						       _gap_extend ),
					_mm_add_epi16( _m2,
						       _gap_open ) ),
			 _nuc_prior );

    //
    // S odd
    //

    // seq1 needs updating; seq2 is current
    char c = (8 + s/2 < len1) ? seq1[ 8+(s/2) ] : 'N';
    _seq1win   = _mm_insert_epi16( _mm_srli_si128( _seq1win,   2 ), c, 8-1 );
    _seq1nqual = _mm_insert_epi16( _mm_srli_si128( _seq1nqual, 2 ), (c=='N') ? n_score : pos_inf, 8-1 );
    _gap_open  = _mm_insert_epi16( _mm_slli_si128( _gap_open, 2 ),
				   localgapopen[ 8 + s/2 < len1 ? 8 + s/2 : len1-1 ],
				   0 );

    // Moved up:
    //_m2 = _mm_andnot_si128( _initmask, _m2 );
    _initmask = _mm_slli_si128( _initmask, 2 );
    _initmask2 = _mm_slli_si128( _initmask2, 2 );
    _m2 = _mm_min_epi16( _m2, _mm_min_epi16( _i2, _d2 ) );

    // at this point, extract minimum score.  Referred-to position must
    // be y==len2-1, so that current position has y==len2; i==0 so d=0 and y=s/2
    if (s/2 >= len2) {
        _score = (*extractors[s/2 - len2])(_m2);

      if (_score < minscore) {
          minscore = _score;
          minscoreidx = s+1;
      }
    }

    _m2 = _mm_add_epi16( _m2, 
			 _mm_min_epi16( _mm_andnot_si128( _mm_cmpeq_epi16( _seq2win, 
									   _seq1win ), 
							  _qual2win ),
					_seq1nqual ) );

    _d2 = _mm_min_epi16( _mm_add_epi16( _d1,
					_gap_extend ),
			 _mm_add_epi16( _mm_min_epi16( _m1,
						       _i1 ),  // allow I->D 
					_gap_open ) );


    _i2 = _mm_insert_epi16( _mm_srli_si128( _mm_add_epi16( _mm_min_epi16( _mm_add_epi16( _i1,
											 _gap_extend ),
									  _mm_add_epi16( _m1,
											 _gap_open ) ),
							   _nuc_prior ),
					    2 ),
			    pos_inf,
			    8-1 );

  }

  return (minscore + 0x8000) >> 2;
}

//_________________________________________________________________________________________________
