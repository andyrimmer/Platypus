import pysam
import sys


f = pysam.Samfile(sys.argv[1], 'rb')

for index,read in enumerate(f.fetch("20", 0, 1000000)):

    if index < 5:
        print read.qual
