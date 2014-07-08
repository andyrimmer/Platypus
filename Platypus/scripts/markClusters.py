import sys

lastPos = 0
lastLine =  None
threshold = 25
cluster = []

if len(sys.argv) > 1:
    threshold = int(sys.argv[1])

for line in sys.stdin:

    if line[0] == "#":
        print line.strip()
    else:

        cols = line.split("\t")
        chrom,pos = cols[0:2]
        pos = int(pos)

        if len(cluster) == 0:
            cluster.append( (chrom,pos,line) )
            continue
        else:
            if chrom == cluster[-1][0] and abs(pos - cluster[-1][1]) <= threshold:
                cluster.append( (chrom,pos,line) )

            elif len(cluster) > 1:

                for theTuple in cluster:
                    theLine = theTuple[-1]
                    theCols = theLine.split("\t")

                    if theCols[6] == "PASS":
                        theCols[6] = "clustered"
                    else:
                        theCols[6] += ";clustered"

                    print "\t".join(theCols).strip()

                cluster = []
                cluster.append( (chrom,pos,line) )

            else:
                print cluster[-1][-1].strip()
                cluster = []
                cluster.append( (chrom,pos,line) )

# Catch last one
if len(cluster) > 1:

    for theTuple in cluster:
        theLine = theTuple[-1]
        theCols = theLine.split("\t")

        if theCols[6] == "PASS":
            theCols[6] = "clustered"
        else:
            theCols[6] += ";clustered"
elif len(cluster) == 1:
    print cluster[-1][-1].strip()
else:
    pass
