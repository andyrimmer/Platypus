from __future__ import division
import sys

for line in sys.stdin:

    if line[0] == "#":
        continue

    try:
        cols = line.strip().split("\t")
        nVar = int(cols[9].split(":")[-1])
        nTot = int(cols[9].split(":")[-2])

        if nVar / nTot >= 0.30:
            print line.strip()
    except Exception:
        print line.strip()
