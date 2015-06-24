from __future__ import division
import sys

for line in sys.stdin:

    try:

        if line[0] == "#":
            print line.strip()
            continue

        cols = line.strip().split("\t")
        info = cols[7]
        TR = None
        TCR = None

        for infoVal in info.split(";"):
            name,value = infoVal.split("=")[0:2]

            if name == "TR":
               TR = int(value)

            if name == "TCR":
               TCR = int(value)

        if TR/TCR > 0.3:
            print line.strip()

    except Exception:
        continue
