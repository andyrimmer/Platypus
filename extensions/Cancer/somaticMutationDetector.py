"""
Simple program to find somatic mutations in tumour samples.
"""

import sys
import optparse
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

def parseOpts(args):
    """
    Run the Platypus variant-caller, with the specified arguments
    """
    parser = optparse.OptionParser()

    parser.add_option("--inputVCF", dest="inputVCF", help="Name of input VCF file", action='store', type='string', default=None)
    parser.add_option("--outputVCF", dest="outputVCF", help="Name of output VCF file", action='store', type='string', default=None)
    parser.add_option("--tumourSample", dest="tumourSample", help="Name of tumour sample", action='store', type='string', default=None)
    parser.add_option("--normalSample", dest="normalSample", help="Name of normal sample", action='store', type='string', default=None)
    parser.add_option("--minPosterior", dest="minPosterior", help="Minimum allowed value for somatic variant posterior", action='store', type='int', default=5)

    (options, args) = parser.parse_args(args)
    return options

###################################################################################################

options = parseOpts(sys.argv)

with open(options.inputVCF, 'r') as vcfFile:
    with open(options.outputVCF, 'w') as outputFile:
        threshold = options.minPosterior
        normalCol = None
        tumourCol = None

        for line in vcfFile:

            if line.startswith("##"):
                outputFile.write(line)
                continue

            elif line.startswith("#CHROM"):
                outputFile.write(line)
                cols = line.strip().split("\t")

                try:
                    normalCol = cols.index(options.normalSample)
                except ValueError:
                    print "\n\nCould not find normal sample %s in input VCF header\n\n" %(options.normalSample)
                    raise

                try:
                    tumourCol = cols.index(options.tumourSample)
                except ValueError:
                    print "\n\nCould not find tumour sample %s in input VCF header\n\n" %(options.tumourSample)
                    raise

            else:

                cols = line.strip().split("\t")
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
                    outputFile.write("\t".join(cols) + "\n")
