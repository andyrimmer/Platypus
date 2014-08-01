import sys

lastPos = 0
lastLine =  None
threshold = int(sys.argv[1])
cluster = []

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
                    print theTuple[-1].strip()

                cluster = []
                cluster.append( (chrom,pos,line) )

            else:
                cluster = []
                cluster.append( (chrom,pos,line) )

# Catch last one
if len(cluster) > 1:
    for theTuple in cluster:
        print theTuple[-1].strip()
