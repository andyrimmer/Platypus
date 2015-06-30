from __future__ import division
import sys

for line in sys.stdin:

    if line[0] == "#":
        continue

    #try:
    cols = line.strip().split("\t")
    nTot = int(cols[9].split(":")[-1])

    if nTot >= 15 and nTot <= 50:
        print line.strip()
    #except Exception:
    #    print line.strip()
    #    continue
