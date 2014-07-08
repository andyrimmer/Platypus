"""
Simple program to detect de-novo variant candidates in child-parent trios
"""

import sys
from math import log,exp,log10

###################################################################################################

def computeSomaticPosterior(callQuality, fatherGLs, motherGLs, childGLs):
    """
    To compute the posterior for a specific variant being somatic, i.e. present in the
    tumour sample but not in the normal, we need to compute the bayesian normalisation which,
    in this model, is the sum of all (non-log) likelihood combinations across all possible
    genotypes:

    sum_ijk ( GL_ijk )

    Here the 2 lists are ordered [homRef, het, homVar], and the likelihoods are in log10 space. Summing all the (10**GL1)*(10**GL2)
    values will result in underflows, so here we re-normalise by the largest value.

    The posterior can then be computed as follows:

    post = prior * sum(allDeNovoGLCombintations) / sum(allGLCombinations)

    then we cap this by the original variant-call quality:

    post = max(phredPost, callQuality)

    """
    deNovoPrior = log10(1e-6) # Fudge.
    normalPrior = log10(1.0- 1e-6)

    deNovoLikelihoods = []
    allLikelihoods = []

    for motherIndex in range(3):
        for fatherIndex in range(3):
            for childIndex in range(3):

                # DeNovo ref --> het combination
                if fatherIndex == 0 and motherIndex == 0 and (childIndex == 1 or childIndex == 2):
                    deNovoLikelihoods.append( deNovoPrior + childGLs[childIndex] + fatherGLs[fatherIndex] + motherGLs[motherIndex] )
                    allLikelihoods.append( deNovoPrior +  childGLs[childIndex] + fatherGLs[fatherIndex] + motherGLs[motherIndex] )

                # DeNovo het --> hom combination
                #elif (fatherIndex == 1 and motherIndex == 0) or (fatherIndex == 0 and motherIndex == 1) and childIndex == 2:
                #    deNovoLikelihoods.append( deNovoPrior + childGLs[childIndex] + fatherGLs[fatherIndex] + motherGLs[motherIndex] )
                #    allLikelihoods.append( deNovoPrior +  childGLs[childIndex] + fatherGLs[fatherIndex] + motherGLs[motherIndex] )

                # Normal combination
                else:
                    allLikelihoods.append( normalPrior + childGLs[childIndex] + fatherGLs[fatherIndex] + motherGLs[motherIndex] )

    # Now re-normalise likelihood combinations to avoid underflows
    maxLike = max(allLikelihoods + deNovoLikelihoods)

    sumAllLikelihoods = sum( [10**(x - maxLike) for x in allLikelihoods] )
    sumDeNovoLikelihoods = sum( [10**(x - maxLike) for x in deNovoLikelihoods] )
    posterior = sumDeNovoLikelihoods/sumAllLikelihoods

    phredPosterior = max(0,  int(0.5 + (-10*log10( max(1e-10, 1.0 - posterior)))) )
    return min(callQuality, phredPosterior)

###################################################################################################

childCol = int(sys.argv[1])
motherCol = int(sys.argv[2])
fatherCol = int(sys.argv[3])
threshold = int(sys.argv[4])

for index,line in enumerate(sys.stdin):

    if line[0] == "#":
        print line.strip()
        continue

    cols = line.strip().split('\t')
    chrom,pos,ref,alt = cols[0],cols[1],cols[2],cols[3]
    callQuality = int(cols[5])

    # Exclude multi-allelic sites for now
    if "," in alt:
        continue

    child = cols[childCol].split(":")[0].split("/")
    mother = cols[motherCol].split(":")[0].split("/")
    father = cols[fatherCol].split(":")[0].split("/")

    # Exclude sites with missing genotypes
    if "." in child or "." in mother or "." in father:
        continue

    childGLs = [float(x) for x in cols[childCol].split(":")[1].split(",")]
    fatherGLs = [float(x) for x in cols[fatherCol].split(":")[1].split(",")]
    motherGLs = [float(x) for x in cols[motherCol].split(":")[1].split(",")]

    deNovoPhredPosterior = computeSomaticPosterior(callQuality, fatherGLs, motherGLs, childGLs)

    if deNovoPhredPosterior >= threshold:
        cols[5] = str(deNovoPhredPosterior)
        print "\t".join(cols)
