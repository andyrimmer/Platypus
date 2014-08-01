import sys

lastPos = 0
lastLine =  None
threshold = int(sys.argv[1])
cluster = []
nCluster = 0

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

                nCluster += 1
                cluster = []
            else:
                cluster = []

print "There are %s clusters" %(nCluster)
