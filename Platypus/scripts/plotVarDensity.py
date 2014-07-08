import sys
from matplotlib import pylab

snpPos = []
indelPos = []

for line in sys.stdin:

    if line[0] == "#":
        continue

    cols = line.strip().split("\t")
    ref = cols[3]
    alt = cols[4]

    if "," in alt:
        continue
    elif len(alt) == len(ref):
        snpPos.append(int(cols[1]))
    else:
        indelPos.append(int(cols[1]))

pylab.subplot(121)
pylab.hist(snpPos, bins=500, label='snps')
pylab.subplot(122)
pylab.hist(indelPos, bins=1000)
pylab.hist(snpPos, bins=500, label='indels')
pylab.legend()
pylab.savefig("positions.png")
