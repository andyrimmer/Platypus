import sys

for line in sys.stdin:

    try:

        if line[0] == "#":
            print line.strip()
            continue

        cols = line.strip().split("\t")
        info = cols[7]
        TR = None
        TU = None

        for infoVal in info.split(";"):
            name,value = infoVal.split("=")[0:2]

            if name == "TR":
                TR = value

            if name == "TU":
                TU = value

        if int(TR)*len(TU) > 10 and len(TU) > 1:
            continue

        if int(TR)*len(TU) > 5 and len(TU) == 1:
            continue

        print line.strip()

    except Exception:
        continue
