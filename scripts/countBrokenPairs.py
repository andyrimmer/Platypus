import sys
import pysam
from collections import defaultdict

bamFile = pysam.AlignmentFile(sys.argv[1], 'rb')


#for chrom in [str(x) for x in range(1,23)] + ["X"] + ["Y"]:
for chrom in [str(x) for x in range(20,21)]:
    nBroken = 0
    nReads = 0
    chromsWithBrokenPairs = {}
    for index,read in enumerate(bamFile.fetch(chrom)):
        nReads += 1
        if read.reference_id != read.next_reference_id:
            #print read.reference_id,read.next_reference_id
            chromsWithBrokenPairs[bamFile.get_reference_name(read.next_reference_id)] += 1
            nBroken += 1

        if index % 1000000 == 0:
            print("N broken pairs for chrom %s = %s out of %s reads" %(chrom, nBroken,nReads))

    print("List of chroms with broken pairs that map to chrom %s" %(chrom))
    for theChrom in sorted(chromsWithBrokenPairs.keys()):
        print(theChrom, chromsWithBrokenPairs[theChrom])
