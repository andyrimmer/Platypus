from __future__ import division
import sys

for line in sys.stdin:
    try:
        cols = line.split("\t")
        chars = cols[4].upper()
        lenChars = int(cols[3])
        nRef = chars.count(".") + chars.count(",")
        nNonRef = chars.count("A") + chars.count("C") + chars.count("T") + chars.count("G")
        print "N ref = %s (%s %%). N non-ref = %s (%s %%)" %(nRef, 100.0*(nRef/lenChars), nNonRef, 100.0*(nNonRef/lenChars))
    except Exception:
        continue
