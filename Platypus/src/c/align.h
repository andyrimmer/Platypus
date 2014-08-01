/*****************************************************************************************************************
 This code is copyright (c) Gerton Lunter, Jan 2009
 It may not be distributed, made public, or used in other software without the permission of the copyright holder
******************************************************************************************************************/

int fastAlignmentRoutine(char* seq1, char* seq2, char* qual2, int len1, int len2, int gapextend, int nucprior, char* homopolgapq_or_localgapopen, char usehomopolgapq, char* aln1, char* aln2);
int fastAlignmentRoutine2(char* seq1, char* seq2, char* qual2, int len1, int len2, int gapextend, int nucprior, short* localGapOpen);
