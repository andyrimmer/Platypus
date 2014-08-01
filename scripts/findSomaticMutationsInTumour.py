"""
Simple program to find somatic mutations in tumour samples.
"""

import sys
from math import log,exp,log10

###################################################################################################

def computeSomaticPosterior(callQuality, normGLs, tumourGLs):
    """
    To compute the posterior for a specific variant being somatic, i.e. present in the
    tumour sample but not in the normal, we need to compute the bayesian normalisation which,
    in this model, is the sum of all (non-log) likelihood combinations across all possible
    genotypes:

    sum ( [0/0 * 0/0], [0/1 * 0/1], [0/0 * 0/1], [1/0 * 0/0], [1/1 * 1/1], [1/1 * 0/0], [0/0 * 1/1], [1/1 * 0/1], [0/1 * 1/1] )

    Here the 2 lists are ordered [homRef, het, homVar], and the likelihoods are in log10 space. Summing all the (10**GL1)*(10**GL2)
    values will result in underflows, so here we re-normalise by the largest value.

    The posterior can then be computed as follows:

    post = prior * sum(allDeNovoGLCombintations) / sum(allGLCombinations)

    then we cap this by the original variant-call quality:

    post = max(phredPost, callQuality)

    """
    somaticPrior = log10(1e-6) # Fudge. Depending on tumour type, could make this higher or lower
    normalPrior = log10(1.0- 1e-6)

    somaticLikelihoods = []
    allLikelihoods = []

    for normIndex in range(3):
        for tumIndex in range(3):

            # Somatic combination
            if normIndex == 0 and (tumIndex == 1 or tumIndex == 2):
                somaticLikelihoods.append( somaticPrior + normGLs[normIndex] + tumourGLs[tumIndex] )
                allLikelihoods.append( somaticPrior + normGLs[normIndex] + tumourGLs[tumIndex] )
            else:
                allLikelihoods.append( normalPrior + normGLs[normIndex] + tumourGLs[tumIndex] )

    # Now re-normalise likelihood combinations to avoid underflows
    maxLike = max(allLikelihoods + somaticLikelihoods)

    sumAllLikelihoods = sum( [10**(x - maxLike) for x in allLikelihoods] )
    sumSomaticLikelihoods = sum( [10**(x - maxLike) for x in somaticLikelihoods] )
    posterior = sumSomaticLikelihoods/sumAllLikelihoods

    phredPosterior = max(0, int( (-10*log10( max(1e-10, 1.0 - posterior))) + 0.5))
    #print normGLs, tumourGLs, sumSomaticLikelihoods,sumAllLikelihoods,posterior, min(callQuality, phredPosterior)
    return min(callQuality, phredPosterior)

###################################################################################################

vcfFile = sys.stdin
order = sys.argv[1]
threshold = int(sys.argv[2])
normalCol = None
tumourCol = None

if order == "TN":
    tumourCol = 9
    normalCol = 10
elif order == "NT":
    tumourCol = 10
    normalCol = 9
else:
    raise StandardError, "order argument must be 'TN' or 'NT' Input was %s" %(order)

for line in vcfFile:

    if line[0] == "#":
        print line.strip()
        continue

    cols = line.strip().split('\t')
    chrom,pos = cols[0],cols[1]
    alt = cols[4]
    tumour = cols[tumourCol].split(":")[0].split("/")
    normal = cols[normalCol].split(":")[0].split("/")
    callQuality = int(cols[5])

    # Exclude multi-allelic sites for now
    if "," in alt:
        continue

    # Skip sites with missing genotypes
    if "." in tumour or "." in normal:
        continue

    tumourGLs = [float(x) for x in cols[tumourCol].split(":")[1].split(",")]
    normalGLs = [float(x) for x in cols[normalCol].split(":")[1].split(",")]
    somaticPosterior = computeSomaticPosterior(callQuality, normalGLs, tumourGLs)

    if somaticPosterior >= threshold:
        cols[5] = str(somaticPosterior)
        print "\t".join(cols)

vcfFile.close()
