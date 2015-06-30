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
    nTs = defaultdict(int)
    nTv = defaultdict(int)

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

                if name == "HP":
                    pal = int(value)

            alleles = sorted([ref,alt])

            if alleles == ["A","G"] or alleles == ["C","T"]:
                nTs[pal//binSize] += 1
            else:
                nTv[pal//binSize] += 1
        except Exception:
            print "Error. Prob last line..."
            print line
            continue

    if binSize == 1:
        print "Max Palindrome Size\tnSnp\tTsTv"

    for start in sorted(nTs.keys()):
        nSnp = nTs[start] + nTv[start]
        TsTv  = None

        if nTv[start] > 0:
            TsTv = nTs[start]/nTv[start]
        else:
            TsTv = -1.0

        if binSize == 1:
            print "%s\t%s\t%1.2f" %(start, nSnp,TsTv)
        else:
            print "For %s <= PAL < %s, nSNP = %s. TsTv = %1.2f" %(start*binSize, (1+start)*binSize, nSnp,TsTv)

###################################################################################################

if __name__ == "__main__":
    binSize = int(sys.argv[1])
    summariseVariantCalls(binSize)

###################################################################################################
