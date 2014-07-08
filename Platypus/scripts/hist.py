import sys
from matplotlib import pylab

res = []
fileName = sys.argv[1]
nBins = int(sys.argv[2])

for line in sys.stdin:
    res.append(float(line.strip()))


if fileName == "-":
    pylab.hist(res, bins=nBins)
    pylab.show()
else:
    pylab.hist(res, bins=nBins)
    pylab.savefig(fileName)
