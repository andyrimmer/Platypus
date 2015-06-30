import sys

threshold = int(sys.argv[1])

for line in sys.stdin:

    try:

        if line[0] == "#":
            print line.strip()
            continue

        cols = line.strip().split("\t")
        info = cols[7]

        for infoVal in info.split(";"):
            name,value = infoVal.split("=")[0:2]

            if name == "TU" and len(value) != threshold:
                print line.strip()

    except Exception:
        continue
