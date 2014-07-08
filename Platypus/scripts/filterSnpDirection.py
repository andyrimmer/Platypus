"""
Filter VCF and output only SNP calls in a sequence context where there
is a different homopolymer either side of the site, one of which is the same base as the ALT
allele, and one is the ref e.g. ref sequence is AAAAATTTTT and SNP is T-->A or A --> A.
"""

import sys
import pysam

flag = sys.argv[1]

for line in sys.stdin:

    if line[0] == "#":
        sys.stdout.write(line)
    else:

        cols = line.split("\t")
        chrom = cols[0]
        pos = int(cols[1]) - 1 # 0-indexed
        ref = cols[3]
        alt = cols[4]
        info = cols[7]
        context = None

        if len(ref) != 1 or len(alt) != 1:
            continue

        for infoVal in info.split(";"):
            name,value = infoVal.split("=")[0:2]

            if name == "SC":
                context = value
                break

        assert context[10] == ref

        if flag == "lr" and context[11] == alt:
            sys.stdout.write(line)
        elif flag == "rl" and context[9] == alt:
            sys.stdout.write(line)
        else:
            continue
