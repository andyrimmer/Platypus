import sys

###################################################################################################

def zopen(fileName, mode):
    if fileName.endswith(".gz"):
        return gzip.open(fileName, mode)
    else:
        return open(fileName, mode)

###################################################################################################

fileName = sys.argv[1]
theFile = zopen(fileName, 'r')

for line in theFile:
    if line.startswith("#"):
        continue
    else:
        chrom,pos,theId,ref,alt = line.split("\t")[0:5]
        alts = alt.split(",")

        for theAlt in alts:
            if theAlt[0] != ref[0]:
                break
            elif len(theAlt) != len(ref):
                break
        else:
            print "Error in VCF. Need to trim leading alt and ref allele sequence"
            print "Line is %s" %(line)
