import sys

for line in sys.stdin:

    cols = line.strip().split("\t")
    alt = cols[4]

    if "," in alt:
        continue

    nr = int(cols[-1].split(":")[-1])

    if nr == 0:
        print line.strip()
