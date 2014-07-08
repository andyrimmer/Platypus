import pysam
import sys
from palindrome import lcs as fastLCS

min_matching_stretch = 999    # require a dimer at minimum

###################################################################################################

def lcs(s, t, a=0, b=1e10):
    """
    longest common substring, which includes at least one character from s[a:b].  Returns length, starting positions in s and t
    """
    l0 = [0] * len(t)   # current row
    l1 = [0] * len(t)   # previous row
    z = 0               # lcs
    starts = -1
    startt = -1
    for i, sc in enumerate(s.upper()):
        for j, tc in enumerate(t.upper()):
            if sc == tc:
                if i==0 or j==0:
                    if i<b:
                        l0[j] = 1
                    else:
                        l0[j] = 0
                else:
                    if i<b or l1[j-1]>0:
                        l0[j] = l1[j-1] + 1
                    else:
                        l0[j] = 0
                if l0[j] >= z and i >= a:
                    if l0[j] > z or abs( startt + (z - len(t))//2 ) > abs( j-z+1 + (z - len(t)//2) ):
                        z = l0[j]
                        starts = i-z+1
                        startt = j-z+1
            else:
                l0[j] = 0
        l0, l1 = l1, l0

    return z, starts, startt

###################################################################################################

def get_max_palindrome(chrom, pos, fa, ref, alt, windowsize):
    """
    Returns length of longest palindromic match between the ref and alt alleles, which overlaps
    the longest allele by at least 1 nucleotide.  This is symmetric for insertions vs deletions.
    Choose ref == alt == '' for a palindromic match of the reference sequence
    """
    seq = fa.fetch(chrom, pos - windowsize, pos + windowsize + max(len(ref),len(alt))).upper()
    assert seq[ windowsize : windowsize + len(ref) ] == ref
    seq2 = seq[ :windowsize ] + alt + seq[ windowsize+len(ref) : ]
    if len(alt)>len(ref):
        # insertion
        #print "Comparing alt to -ref (ins); alt=",seq2[:windowsize+1].lower() + seq2[windowsize+1:windowsize+len(alt)].upper() + seq2[windowsize+len(alt):].lower()
        lng,strt1,strt2 = fastLCS( seq2, revcmp(seq), windowsize+1, windowsize + len(alt) )
        if strt2 > -1:
            return lng, pos - windowsize + (len(seq) - strt2 - lng)
        else:
            return lng, -1
    else:
        # deletion
        #print "Comparing ref to -alt (del); ref=",seq[:windowsize+1].lower() + seq[windowsize+1:windowsize+len(ref)].upper() + seq[windowsize+len(ref):].lower()
        lng,strt1,strt2 = fastLCS( seq, revcmp(seq2), windowsize+1, windowsize + len(ref) ) 
        if strt1 > -1:
            return lng, pos - windowsize + strt1
        else:
            return lng, -1

###################################################################################################

def revcmp( unit ):
    return ''.join(reversed([{'A':'T','T':'A','C':'G','G':'C'}.get(c,'N') for c in unit.upper()]))

###################################################################################################

fa = pysam.Fastafile(sys.argv[1])
palindrome = int(sys.argv[2])
outFile = open(sys.argv[3], 'w')

for line in sys.stdin:

    if line[0] == "#":
        outFile.write(line)
        continue

    cols = line.strip().split("\t")
    chrom =  cols[0]
    pos = int(cols[1]) - 1 # 0-indexed
    ref = cols[3]
    alt = cols[4]

    # Ig palindrome < 0, use only reference, not alt allele
    if palindrome > 0:
        pallen, palpos = get_max_palindrome(chrom, pos, fa, ref, alt, palindrome)
    else:
        pallen, palpos = get_max_palindrome(chrom, pos, fa, ref, ref, -palindrome)

    cols[7] = cols[7][0:-1] + ";PAL=%s" %(pallen)
    outFile.write("\t".join(cols) + "\n")
