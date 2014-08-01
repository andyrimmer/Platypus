"""
This script produces summary plots and tables of the Platypus SNP and indel
calls, and comparison plots of Platypus calls with the validated 1000 genomes
calls.
"""
from __future__ import division

from math import sqrt,pow,log,exp,pi,log10

import sys
import gzip
import pysam

###################################################################################################

def summariseVariantCalls():
    """
    Summarise the variant calls in a given vcf file.
    """
    kgSNPs = set()
    kgSNPFile = pysam.Tabixfile(sys.argv[1], 'r')

    for line in kgSNPFile.fetch():

        cols = line.strip().split("\t")
        chrom = cols[0]
        pos = int(cols[1])
        ref = cols[3]
        alt = cols[4]

        kgSNPs.add("%s:%s:%s:%s" %(chrom,pos,ref,alt))

    nSNPs = 0
    nSNPsIn1Kg = 0
    nSNPsNotIn1Kg = 0
    nPASSSNPs = 0
    nPASSSNPsIn1Kg = 0
    nPASSSNPsNotIn1Kg = 0
    nFAILSNPs = 0
    nFAILSNPsIn1Kg = 0
    nFAILSNPsNotIn1Kg = 0

    for index,line in enumerate(sys.stdin):

        if line[0] == "#":
            continue

        cols = line.split("\t")
        chrom = cols[0].replace("chr", "")
        pos = int(cols[1]) # 0 -indexed to match with pysam
        ref = cols[3]
        alts = cols[4]
        filters = cols[6]

        for alt in alts.split(","):

            if not (len(ref) == 1 and len(alt) == 1):
                continue

            nSNPs += 1

            if filters == "PASS":
                nPASSSNPs += 1
            else:
                nFAILSNPs += 1

            if "%s:%s:%s:%s" %(chrom,pos,ref,alt) in kgSNPs:
                nSNPsIn1Kg += 1

                if filters == "PASS":
                    nPASSSNPsIn1Kg += 1
                else:
                    nFAILSNPsIn1Kg += 1

            else:
                nSNPsNotIn1Kg += 1

                if filters == "PASS":
                    nPASSSNPsNotIn1Kg += 1
                else:
                    nFAILSNPsNotIn1Kg += 1

    print "%s SNPs in total" %(nSNPs)
    print "%s PASS SNPs in total (%s %% of all SNP calls)" %(nPASSSNPs, 100.0*nPASSSNPs/nSNPs)
    print "%s FAIL SNPs in total (%s %% of all SNP calls)" %(nFAILSNPs, 100.0*nFAILSNPs/nSNPs)
    print "%s SNPs (%s %%) in 1kg set" %(nSNPsIn1Kg, 100.0*nSNPsIn1Kg/nSNPs)
    print "%s SNPs (%s %%) not in 1kg set" %(nSNPsNotIn1Kg, 100.0*nSNPsNotIn1Kg/nSNPs)
    print "%s PASS SNPs (%s %%) in 1kg set" %(nPASSSNPsIn1Kg, 100.0*nPASSSNPsIn1Kg/nPASSSNPs)
    print "%s PASS SNPs (%s %%) not in 1kg set" %(nPASSSNPsNotIn1Kg, 100.0*nPASSSNPsNotIn1Kg/nPASSSNPs)
    print "%s FAIL SNPs (%s %%) in 1kg set" %(nFAILSNPsIn1Kg, 100.0*nFAILSNPsIn1Kg/nFAILSNPs)
    print "%s FAIL SNPs (%s %%) not in 1kg set" %(nFAILSNPsNotIn1Kg, 100.0*nFAILSNPsNotIn1Kg/nFAILSNPs)

###################################################################################################

if __name__ == "__main__":
    summariseVariantCalls()

###################################################################################################
