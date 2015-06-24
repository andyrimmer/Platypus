"""
This script produces summary plots and tables of the Platypus SNP and indel
calls, and comparison plots of Platypus calls with the validated 1000 genomes
calls.
"""
from __future__ import division

from math import sqrt,pow,log,exp,pi,log10

import sys
import gzip

###################################################################################################

def summariseVariantCalls():
    """
    Summarise the variant calls in a given vcf file.
    """
    nSNPs = 0
    nTransitionSnps = 0
    nTransversionSnps = 0

    for index,line in enumerate(sys.stdin):

        try:

            if line[0] == "#":
                continue

            cols = line.split("\t")
            chrom = cols[0]
            pos = int(cols[1])
            ref = cols[3]
            alt = cols[4]
            qual = None
            alts = None

            nSNPs += 1
            alleles = sorted([ref,alt])

            if alleles == ["A","G"] or alleles == ["C","T"]:
                nTransitionSnps += 1
            else:
                nTransversionSnps += 1
        except Exception:
            #print line
            continue
            #raise

    print "nSNP = %s. \t TsTv = %s" %(nSNPs,nTransitionSnps/nTransversionSnps)

###################################################################################################

if __name__ == "__main__":
    summariseVariantCalls()

###################################################################################################
