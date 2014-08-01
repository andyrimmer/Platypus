import sys
from matplotlib import pylab

gofs = []

for line in sys.stdin:

    if line[0] == "#":
        continue

    cols = line.strip().split("\t")
    gofs.append(int(cols[9].split(":")[-4]))

pylab.hist(gofs, bins=25)
pylab.savefig("gofs.png")
