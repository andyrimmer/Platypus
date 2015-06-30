import sys

threshold = int(sys.argv[1])

for line in sys.stdin:

    if line[0] == "#":
        continue

    try:
        cols = line.strip().split("\t")
        if int(cols[9].split(":")[-4]) < threshold:
            print line.strip()
    except Exception:
        print line.strip()
