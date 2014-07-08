import sys

threshold = int(sys.argv[1])

for line in sys.stdin:

    if line[0] == "#":
        print line.strip()
        continue

    cols = line.strip().split("\t")
    info = cols[7]

    if "," in cols[4]:
        continue

    for infoVal in info.split(";"):
        name,value = infoVal.split("=")[0:2]

        if name == "TR" and int(value) >= threshold:
            print line.strip()
