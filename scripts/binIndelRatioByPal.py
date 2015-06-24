"""
This script produces summary plots and tables of the Platypus SNP and indel
calls, and comparison plots of Platypus calls with the validated 1000 genomes
calls.
"""
from __future__ import division

from math import sqrt,pow,log,exp,pi,log10
from collections import defaultdict

import sys
import gzip

###################################################################################################

def summariseVariantCalls(binSize):
    """
    Summarise the variant calls in a given vcf file.
    """
    nInsertions = defaultdict(int)
    nDeletions = defaultdict(int)

    for index,line in enumerate(sys.stdin):

        try:

            if line[0] == "#":
                continue

            cols = line.split("\t")
            chrom = cols[0]
            pos = int(cols[1])
            ref = cols[3]
            alt = cols[4]
            info = cols[7]
            pal = 0

            for infoVal in info.split(";"):
                name,value = infoVal.split("=")[0:2]

                if name == "PAL":
                    pal = int(value)

            if len(ref) > len(alt):
                nDeletions[pal//binSize] += 1
            else:
                nInsertions[pal//binSize] += 1
        except Exception:
            print "Error. Prob last line..."
            print line
            continue

    if binSize == 1:
        print "Max Palindrome Size\tnSnp\tnInd/nDels"

    for start in sorted(nInsertions.keys()):
        nIndels = nInsertions[start] + nDeletions[start]
        ratio  = None

        if nDeletions[start] > 0:
            ratio = nInsertions[start]/nDeletions[start]
        else:
            ratio = -1.0

        if binSize == 1:
            print "%s\t%s\t%1.2f" %(start, nIndels,ratio)
        else:
            print "For %s <= PAL < %s, nSNP = %s. nIns/nDels = %1.2f" %(start*binSize, (1+start)*binSize, nIndels,ratio)

###################################################################################################

if __name__ == "__main__":
    binSize = int(sys.argv[1])
    summariseVariantCalls(binSize)

###################################################################################################
