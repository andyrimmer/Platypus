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
    totalIndels = 0
    totalIns = 0
    totalDels = 0

    for index,line in enumerate(sys.stdin):

        try:

            if line[0] == "#":
                continue

            cols = line.split("\t")
            chrom = cols[0]
            pos = int(cols[1])
            ref = cols[3]
            alts = cols[4]
            info = cols[7]
            hp = 0

            for infoVal in info.split(";"):
                name,value = infoVal.split("=")[0:2]

                if name == "HP":
                    hp = int(value)

            for alt in alts.split(","):

                if len(ref) > len(alt):
                    totalDels += 1
                    totalIndels += 1
                    nDeletions[hp//binSize] += 1
                elif len(alt) > len(ref):
                    totalIns += 1
                    totalIndels += 1
                    nInsertions[hp//binSize] += 1
                else:
                    pass
        except Exception:
            print "Error. Prob last line..."
            print line
            continue

    print "nInsertions = %s. nDeletions = %s. Total = %s. Ins/Dels = %s" %(totalIns,totalDels,totalIndels,totalIns/totalDels)
    if binSize == 1:
        print "HP Length\tnSnp\tnIns/nDel"

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
            print "For %s <= HP < %s, nIndels = %s. nIns/nDel = %1.2f" %(start*binSize, (1+start)*binSize, nIndels,ratio)

###################################################################################################

if __name__ == "__main__":
    binSize = int(sys.argv[1])
    summariseVariantCalls(binSize)

###################################################################################################
