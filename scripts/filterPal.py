import sys

flag = sys.argv[1][0]
threshold = int(sys.argv[1][1:])

for line in sys.stdin:

    try:

        if line[0] == "#":
            print line.strip()
            continue

        cols = line.strip().split("\t")
        info = cols[7]

        for infoVal in info.split(";"):
            name,value = infoVal.split("=")[0:2]

            if flag == ">":
                if name == "PAL" and int(value) >= threshold:
                    print line.strip()
            elif flag == "<":
                if name == "PAL" and int(value) < threshold:
                    print line.strip()
            elif flag == "=":
                if name == "PAL" and int(value) == threshold:
                    print line.strip()
            else:
                raise StandardError, "Flag should be <,> or = and is %s" %(flag)

    except Exception:
        continue
