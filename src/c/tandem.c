#include <stdlib.h>

#include <stdio.h>   // for testing
#include <string.h>  // for testing

static const int MAX_UNIT_LENGTH = 12;
static const int MIN_PARTIAL_MATCH = 5;

//=================================================================================================

unsigned char* twobit(char* sequence, int length, int offset)
{
    // Helper -- converts ACGT sequence into two-bit-per-nucleotide version,
    // with an offset, and pad to 64-nucleotide (128 bit) words

    // size in bytes, rounded up to whole number of 128-bit words, and make
    // sure that the result ends with an empty 128-bit word

    int buflen = (129 + (((length + MAX_UNIT_LENGTH) * 2 ) | 127)) / 8;
    unsigned char* buffer = malloc( (size_t)buflen );
    int current_nuc_index = offset;
    int i,j;

    // Loop over bytes
    for (i=0; i<buflen; i++)
    {
        // Loop over bit position within byte
        unsigned char byte = 0;
        for (j=0; j<8; j+=2)
        {
            switch (sequence[current_nuc_index] & 0xDF) {
                case 'A':
                    break;
                case 'C':
                    byte |= (1 << j);
                    break;
                case 'G':
                    byte |= (2 << j);
                    break;
                case 'T':
                    byte |= (3 << j);
                    break;
                case 0:
                    --current_nuc_index; // pad with As
                    break;
                default: // convert Ns into pseudo-random sequence
                    byte |= ((current_nuc_index % 257)*(1+current_nuc_index % 257)/2 + (current_nuc_index % 5)) % 4 << j;
                    break;
            }
            ++current_nuc_index;
        }
        buffer[i] = byte;
    }

    return buffer;
}

//=================================================================================================

inline int approximate_indel_rate_inline(int size, int displacement)
{
    // Helper -- returns guess of indel rate, in -10*phred units
    switch (displacement) {
        case 1: return -360 + 24*size;
        case 2: return -327 + 15*size;
        case 3: return -291 + 8*size;
        default: return -282 + 6*size;
    }
}

//=================================================================================================
// external version
int approximate_indel_rate(int size, int displacement)
{
    return approximate_indel_rate_inline(size, displacement);
}

//=================================================================================================

inline int min(int i, int j)
{
    // Helper -- minimum
    if (i<j) return i;
    return j;
}

//=================================================================================================

inline void foundmatch(char* sizes, char* displacements, int pos, int size, int displacement, int length)
{
    // Helper -- decides whether to accept (and store) a match
    char markfull = (length < 0);

    if (markfull)
    {
        length = -length;
    }

    if (pos + displacement + size > length)
    {
        size = length - displacement - pos;
    }

    // convert size to length of repetitive region
    size += displacement;
    // only accept true tandems and sufficiently long partial tandems

    if (size < displacement + min(MIN_PARTIAL_MATCH, displacement))
    {
        return;
    }

    if (approximate_indel_rate_inline(sizes[pos], displacements[pos]) < approximate_indel_rate_inline(size, displacement))
    {
        sizes[pos] = size;
        displacements[pos] = displacement;

        if (markfull)
        {
            int i = pos+1;
            for (; i<min(length, pos + size); i++)
            {
                sizes[i] = size;
                displacements[i] = displacement;
            }
        }
    }
}

//=================================================================================================

void annotate(char* sequence, char* sizes, char* displacements, int length)
{
    // Main function.  Annotates a sequence with the length (sizes) and unit size
    // (displacements) of the local repeat
    // Variable  length  is the length of sequence; sizes and displacements must be
    // at least the same length.
    // If  length  is negative, size/displacement will be marked along the full
    // repeat; if not, only the start (leftmost entry) will be marked.

    int original_length = length;

    if (length < 0)
    {
        length = -length;
    }

    unsigned char* seqs[4] = { twobit( sequence, length, 0 ),
        twobit( sequence, length, 1 ),
        twobit( sequence, length, 2 ),
        twobit( sequence, length, 3 ) };

    // DEBUG
    //int buflen = (129 + (((length + MAX_UNIT_LENGTH) * 2 ) | 127)) / 8;

    // initialize size and displacement arrays
    int pos, displacement;

    for (pos=0; pos < length; pos++)
    {
        sizes[pos] = 1;
        displacements[pos] = 1;
    }

    // loop over starting positions within the sequence
    for (pos=0; pos < length; pos+=4)
    {
        // get 128 bits, 64 nucleotides, at pos, in two 64-bit longs
        // printf("accessing pos/4=%i..%i buflen=%i length=%i size=%i address=%p\n",pos/4,pos/4+8,buflen,length,sizeof(long long),&((seqs[0])[pos/4]));
        unsigned long long original0 = *((long long*) & seqs[0][pos/4]);
        // printf("accessing(2nd) pos/4=%i..%i buflen=%i length=%i size=%i address=%p\n",(pos+32)/4,(pos+32)/4+8,buflen,length,sizeof(long long),&seqs[0][(pos+32)/4]);
        unsigned long long original1 = *((long long*) & seqs[0][(pos+32)/4]);

        for (displacement=1; displacement < MAX_UNIT_LENGTH; displacement+=1) {

            if (pos + displacement >= length)
            {
                break;
            }

            // get target
            unsigned long long target0 = *((long long*) & seqs[displacement % 4][(pos+displacement)/4]);

            // calculate match lengths
            target0 ^= original0;
            int size0 = 64;
            int size1 = 64;
            int size2 = 64;
            int size3 = 64;
            int nucscanleft = (__builtin_ffsll( target0 ) + 1)/2;  // 0 for no mismatches; 1 for nuc 0 mismatch, etc.

            if (nucscanleft == 1)
            {
                size0 = 0;                                         // found result for position 0
                target0 &= -4LL;                                   // remove blocking mismatch
                nucscanleft = (__builtin_ffsll( target0 ) + 1)/2;  // recompute
            }

            if (nucscanleft == 2)
            {
                size0 = min(size0,1);
                size1 = 0;
                target0 &= -16LL;
                nucscanleft = (__builtin_ffsll( target0 ) + 1)/2;
            }

            if (nucscanleft == 3)
            {
                size0 = min(size0, 2);
                size1 = min(size1, 1);
                size2 = 0;
                target0 &= -64LL;
                nucscanleft = (__builtin_ffsll( target0 ) + 1)/2;
            }

            if (nucscanleft == 0)
            {
                if (pos + displacement + 32 < length)
                {
                    unsigned long long target1 = *((long long*) & seqs[displacement % 4][(pos+displacement+32)/4]);
                    target1 ^= original1;
                    nucscanleft = (__builtin_ffsll( target1 ) +1)/2;

                    if (nucscanleft == 0)
                    {
                        nucscanleft = 32;
                    }
                    else
                    {
                        nucscanleft--;
                    }

                }
                else
                {
                    nucscanleft = 0;
                }

                size0 = min(size0, nucscanleft+32);
                size1 = min(size1, nucscanleft+31);
                size2 = min(size2, nucscanleft+30);
                size3 = nucscanleft + 29;
            }
            else
            {
                size0 = min(size0, nucscanleft-1);
                size1 = min(size1, nucscanleft-2);
                size2 = min(size2, nucscanleft-3);
                size3 = nucscanleft - 4;
            }

            foundmatch(sizes, displacements, pos, size0, displacement, original_length);
            foundmatch(sizes, displacements, pos+1, size1, displacement, original_length);
            foundmatch(sizes, displacements, pos+2, size2, displacement, original_length);
            foundmatch(sizes, displacements, pos+3, size3, displacement, original_length);
        }
    }

    free(seqs[0]);
    free(seqs[1]);
    free(seqs[2]);
    free(seqs[3]);
}

//=================================================================================================

int main()
{

    char* seq1 = "TATTTGCATGCGCTTTCGAGCTGTTGAAGAGACGTGTATTGGAATAAGTAATCACATAAGTGTTAGTAACTTATTTAAATACGTATAGAGTCGCCTATTTGCCTAGCCTTTTGGTTCTCAGATTTTTTAATTATTACATTGCTATAAGGGTGTAACTGTGTGATAGCCAAAATTTTAAGCTGCAAATGGTTTGTAAATATGATATATTACAAGCTTCATGAAAATCGGTTTATGACTGATCCGCGATTACGTTGAAAGGCGACTGGCAGAGATACTTTTGTTCAGATGTTTTTTCAGGTAGCGATTCCAATGAATAGGTAAAATACCTTGCAAGTTTTGTTGTTGTCGTTGGAGGAAATGTGGATGTGGTTGTTATTGTTGA";
    char* sizes = (char*)malloc( strlen(seq1)+1 );
    char* displacements = (char*)malloc( strlen(seq1)+1 );

    annotate( seq1, sizes, displacements, -strlen(seq1) );

    int i;

    for (i=0; i<strlen(seq1); i++)
    {
        printf ("%c\t%u\t%u\n",seq1[i], sizes[i], displacements[i]);
    }

    return 0;
}

//=================================================================================================
