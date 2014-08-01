"""
Filter VCF and output only SNP calls in a sequence context where there
is a different homopolymer either side of the site, one of which is the same base as the ALT
allele, and one is the ref e.g. ref sequence is AAAAATTTTT and SNP is T-->A or A --> A.
"""

import sys
import pysam

a = pysam.Fastafile(sys.argv[1])
threshold = int(sys.argv[2])

for line in sys.stdin:

    if line[0] == "#":
        sys.stdout.write(line)
    else:

        cols = line.split("\t")
        chrom = cols[0]
        pos = int(cols[1]) - 1 # 0-indexed
        ref = cols[3]
        alt = cols[4]

        if len(ref) != 1 or len(alt) != 1:
            continue

        leftContext = a.fetch(chrom, pos-20, pos)
        rightContext = a.fetch(chrom, pos+1, pos+21)

        # The HPs must be different bases
        if leftContext[-1] == rightContext[0]:
            continue

        #if not ( (leftContext[-1] == alt and rightContext[0] == ref) or (leftContext[-1] == ref and rightContext[0] == alt) ):
        #    continue

        leftCount = 0
        rightCount = 0

        for i in range(1, 20):

            if leftContext[-i] == leftContext[-1]:
                leftCount += 1
            else:
                break

        for i in range(20):

            if rightContext[i] == rightContext[0]:
                rightCount += 1
            else:
                break

        if leftCount >= threshold and rightCount >= threshold:
            #sys.stdout.write(a.fetch(chrom, pos-20, pos+20) + "\n")
            sys.stdout.write(line)
