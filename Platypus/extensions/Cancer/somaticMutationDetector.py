"""
This program is used to detect and filter de novo SNP and indel mutations called by Platypus
(as well as other methods). This is the code used to provide the lists of de novos reported
in the Platypus paper.
"""
from __future__ import division

import sys
import gzip
import logging
import itertools
import optparse
from math import log10

###################################################################################################

# Set up logger
logger = logging.getLogger("TheLogger")
formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

# Log only INFO level for now
ch = logging.StreamHandler()
ch.setFormatter(formatter)
ch.setLevel(logging.INFO)

logger.addHandler(ch)
logger.setLevel(logging.INFO)

PRIOR_DENOVO = log10(2.0e-8) # Prior based on 2 de novos per 100MB
PRIOR_SNP = log10(1e-3) # 1 SNP per 1000 bases
PRIOR_NO_DENOVO = log10(1.0 - 10**PRIOR_DENOVO - 10**PRIOR_SNP) # The rest

MAX_PARENTAL_VAR_FRAC = 0.03
#MAX_PARENTAL_VAR_FRAC = 0.05
MIN_VAR_FRAC_IN_CHILD = 0.20 # Minimum fraction of variant-supporting reads in child
#MIN_VAR_FRAC_IN_CHILD = 0.25 # Minimum fraction of variant-supporting reads in child
#MIN_READS_IN_CHILD = 0 # Minimum number of supporting variant reads in child
#MIN_READS_IN_CHILD = 10 # Minimum number of supporting variant reads in child
MIN_READS_IN_CHILD = 8 # Minimum number of supporting variant reads in child
logBayesFactorThreshold = log10(1.0)
#logBayesFactorThreshold = log10(1.0/10.0)
#logBayesFactorThreshold = log10(1.0/50.0)
#logBayesFactorThreshold = log10(1.0/4.0)

###################################################################################################

def checkPloidy(chrom, pos, sex):
    """
    Return the appropriate ploidy for the specified sample, based on the chromosome
    and position, and sex of the sample. This attempts to deal with the pseudo-autosomal
    regions on X and Y, by treating these as always diploid.
    """
    if chrom == "X":

        # TODO: Adjust for the pseudo-autosomal regions.
        if sex == "M":
            return 1
        else:
            return 2

    elif chrom == "Y":
        # TODO: Adjust for the pseudo-autosomal regions.
        if sex == "M":
            return 1
        else:
            return 0

    else:
        return 2

###################################################################################################

class Variant(object):
    """
    Simple class to encapsulate a variant line in VCF
    """
    def __init__(self, line, samples):
        """
        Constructor.
        """
        self.line = line
        cols = line.strip().split("\t")

        self.chrom = cols[0].upper().strip("CHR")
        self.pos = int(cols[1])
        self.id = cols[2]
        self.ref = cols[3]
        self.alts = cols[4].split(",")
        self.qual = float(cols[5])
        self.filters = set(cols[6].split(";"))
        self.info = {}

        # Parse INFO field
        for item in cols[7].split(";"):
            try:
                key,value = item.split("=")
                self.info[key] = value.split(",")
            except:
                pass
                #print item
                #raise

        # Parse Sample data
        self.samples = {}
        formatData = cols[8].split(":")

        for sample,sampleData in zip(samples, cols[9:]):

            self.samples[sample] = {}
            sampleDataFields = sampleData.split(":")

            for fieldName,fieldValue in zip(formatData, sampleDataFields):
                self.samples[sample][fieldName] = fieldValue

        # Now extract genotype likelihoods. There are two standard cases here, when the GL or
        # PL tags are used. GL gives likelihoods in log10, PL in phred.
        for sample in samples:

            # Convert GQ into integer
            if "GQ" in self.samples[sample]:

                qualities = []

                for gq in self.samples[sample]['GQ'].split(","):
                    qualities.append(int(gq))

                self.samples[sample]['GQ'] = qualities

            else:
                self.samples[sample]['GQ'] = None

            if "GT" in self.samples[sample]:
                genotypeField = self.samples[sample]['GT']

                if "/" in genotypeField:
                    genotype = tuple(genotypeField.split("/"))
                elif "|" in genotypeField:
                    genotype = tuple(genotypeField.split("|"))
                else:
                    raise StandardError, "Invalid genotype %s in sample %s" %(genotypeField, sample)

                self.samples[sample]['GT'] = genotype
            else:
                raise StandardError, "No recognised genotype tag found in sample data for sample %s. Sample data is %s" %(sample, self.samples[sample])

            if "GL" in self.samples[sample]:
                likelihoods = []

                for gl in self.samples[sample]['GL'].split(","):
                    likelihoods.append(float(gl))

                self.samples[sample]['GL'] = likelihoods

            elif "PL" in self.samples[sample]:
                likelihoods = []

                for gl in self.samples[sample]['PL'].split(","):
                    likelihoods.append(-0.1*float(gl))

                self.samples[sample]['GL'] = likelihoods

            else:
                self.samples[sample]['GL'] = None
                #raise StandardError, "No recognised likelihood tag found in sample data for sample %s. Sample data is %s" %(sample, self.samples[sample])

            if "NV" in self.samples[sample]:
                nv = [int(x) for x in self.samples[sample]['NV'].split(",")]
                self.samples[sample]['NV'] = nv

            elif "AD" in self.samples[sample]:
                nv = [int(x) for x in self.samples[sample]['AD'].split(",")[1:]]
                self.samples[sample]['NV'] = nv

            else:
                self.samples[sample]['NV'] = None

            if "NR" in self.samples[sample]:
                nr = [int(x) for x in self.samples[sample]['NR'].split(",")]
                self.samples[sample]['NR'] = nr

            elif "AD" in self.samples[sample]:
                nr = [sum([int(x) for x in self.samples[sample]['AD'].split(",")])]
                self.samples[sample]['NR'] = nr
            else:
                self.samples[sample]['NR'] = None

    def adjustGenotypesAndLikelihoodsForPloidy(self, sampleSex):
        """
        Based on the location of the variant, adjust the genotypes and likelihoods
        to match the ploidy expected. For example, X is haploid in males.
        """
        for sample, in self.samples.keys():

            ploidyOfSample = checkPloidy(self.chrom, self.pos, sampleSex)

            if ploidyOfSample == 1:
                likelihoods = self.samples[sample]['GL']

                if likelihoods[2] > likelihoods[0]:
                    self.samples[sample]['GT'] = ("1",)
                else:
                    self.samples[sample]['GT'] = ("0",)

                self.samples[sample]['GL'] = (likelihoods[0], likelihoods[2])

            elif ploidyOfSample == 0:
                    self.samples[sample]['GT'] = ()
                    self.samples[sample]['GL'] = ()

            elif ploidyOfSample == 2:
                pass

            else:
                raise StandardError, "Unexpected ploidy %s" %(ploidyOfSample)

###################################################################################################

def parseSampleNamesVCF(vcfLine, normalSampleName, tumourSampleName):
    """
    Take a .ped file and parse it to extract sample names
    in the order "child father mother". The final column gives
    the sex of the child.
    """
    samples = [sample for sample in vcfLine.strip().split("\t")[9:] if sample == normalSampleName or sample == tumourSampleName]

    if len(samples) < 2:
        raise StandardError, "Not enough samples in VCF input file. Need at least 2 samples for tumour/normal comparison."

    if normalSampleName not in samples:
        raise StandardError, "Could not find normal sample (%s) in VCF. Samples found were %s" %(normalSampleName, samples)

    if tumourSampleName not in samples:
        raise StandardError, "Could not find tumour sample (%s) in VCF. Samples found were %s" %(tumourSampleName, samples)

    return samples

###################################################################################################

def isMendelError(variant, pedigree, sexOfChild):
    """
    Return True if the specified variant is inherited in a way not consistent with
    the mendelian pattern, and False otherwise.

    This works by looking at the genotypes of the 3 samples, extracted from a VCF
    file.

    NB. The genotype quality threshold was not being applied to non-Platypus calls,
    previously.
    """
    GQThreshold = 30 # Threshold on genotype quality

    childSample = pedigree['Child']
    motherSample = pedigree['Mother']
    fatherSample = pedigree['Father']

    childGenotype = variant.samples[childSample]['GT']
    fatherGenotype = variant.samples[fatherSample]['GT']
    motherGenotype = variant.samples[motherSample]['GT']

    if '.' in childGenotype:
        return False

    if '.' in motherGenotype:
        return False

    if '.' in fatherGenotype:
        return False

    GQChild = variant.samples[childSample]['GQ'][0]
    GQMother = variant.samples[motherSample]['GQ'][0]
    GQFather = variant.samples[fatherSample]['GQ'][0]

    # This is not optimal for X and Y chromosomes, as these
    # may be haploid. TODO: fix.
    if min(GQChild, GQFather, GQMother) < GQThreshold:
        return False


    # Genotypes that could be inherited by the child, according to the mendelian
    # pattern. If the observed genotype is not one of these then we have a mendelian
    # error.
    if variant.chrom == "X" and checkPloidy(variant.chrom, variant.pos, sexOfChild) == 1:
        for genotype in motherGenotype:
            if genotype == childGenotype[0]:
                return False
        else:
            return True

    elif variant.chrom == "Y" and checkPloidy(variant.chrom, variant.pos, sexOfChild) == 1:

        # TODO: This ignores the pseudo-autosomal Y-chrom regions in females. Fix this.
        if sexOfChild == "F":
            return False

        for genotype in fatherGenotype:
            if genotype == childGenotype[0]:
                return False
        else:
            return True

    else:
        for possibleChildGenotype in itertools.product(fatherGenotype, motherGenotype):
            if childGenotype == possibleChildGenotype or tuple(reversed(childGenotype)) == possibleChildGenotype:
                return False # Found matching mendelian inheritance pattern
        else:
            return True # Found no matching mendelian inheritance pattern

###################################################################################################

def isSomatic(variant, sampleSex, normalSampleName, tumourSampleName):
    """
    Return True if the specified variant represents a somatic mutation in the
    tumour sample, and False otherwise.
    """
    normalGenotype = variant.samples[normalSampleName]['GT']
    tumourGenotype = variant.samples[tumourSampleName]['GT']

    # Total coverage in each sample
    fatherNV = variant.samples[fatherSample]['NV'][0]
    motherNV = variant.samples[motherSample]['NV'][0]
    childNV = variant.samples[childSample]['NV'][0]

    # Variant coverage in each sample
    fatherNR = variant.samples[fatherSample]['NR'][0]
    motherNR = variant.samples[motherSample]['NR'][0]
    childNR = variant.samples[childSample]['NR'][0]

    if childNR == 0 or childNV / childNR < MIN_VAR_FRAC_IN_CHILD:
        return False

    if childNV < MIN_READS_IN_CHILD:
        return False

    if fatherNR == 0 or motherNR == 0:
        return False # No reads in one of the parents

    if fatherNR == 0 or fatherNV / fatherNR >= MAX_PARENTAL_VAR_FRAC:
        return False # Variant read frac in father is too high

    if motherNR == 0 or motherNV / motherNR >= MAX_PARENTAL_VAR_FRAC:
        return False # Variant read frac in mother is too high

    # Normal case, on autosomes
    if fatherGenotype ==  ("0", "0") and motherGenotype == ("0", "0"):

        if childGenotype != ("0", "0"):
            return True
        else:
            return False

    elif fatherGenotype == ("0",) and motherGenotype == ("0", "0"):

        if sexOfChild == "M":
            if childGenotype != ("0",):
                return True
            else:
                return False

        elif sexOfChild == "F":
            if childGenotype != ("0","0"):
                return True
            else:
                return False

        else:
            raise StandardError, "?"

    elif fatherGenotype == ("0",) and motherGenotype == ():

        if sexOfChild == "M":
            if childGenotype != ("0",):
                return True
            else:
                return False

    # We are currently only considering de-novos arising on a reference background
    elif "1" in fatherGenotype or "1" in motherGenotype:
        return False

    # Something I forgot?
    else:
        print "Should not be here.... variant = (%s:%s %s --> %s). Genotypes = child : %s, father: %s, mother %s" %(variant.chrom, variant.pos, variant.ref, variant.alts, childGenotype, fatherGenotype, motherGenotype)

    return False # Genotype pattern does not match de novo

###################################################################################################

def computeBayesFactor(childGLs, fatherGLs, motherGLs, variant, sexOfChild):
    """
    Compute and return the Bayes Factor for the de novo model.
    """
    noDeNovoPatterns = None
    deNovoPatterns = None

    if checkPloidy(variant.chrom, variant.pos, "M") == 2:
        # Genotype patterns are in order (Child, Mother, Father). Numbers are used to index into likelihod array. So a value of 1 for
        # the father is either het or hom alt depending on the location.
        noDeNovoPatterns = ((0,0,1), (0,1,0), (0,1,1), (1,0,1), (1,0,2), (1,1,0), (1,1,1), (1,1,2), (1,2,0), (1,2,1), (2,1,1), (2,1,2), (2,2,1), (2,2,2))
        deNovoPatterns = ((0,0,2), (0,1,2), (0,2,0) , (0,2,1), (0,2,2), (1,0,0), (1,2,2), (2,0,0), (2,0,1), (2,0,2), (2,1,0), (2,2,0))

    else:

        if sexOfChild == "F":
            # TODO: Mother is never homozygous here. Need to add patterns.
            noDeNovoPatterns = ((0,1,0), (1,1,0), (1,1,1), (1,0,1), (1,2,0), (2,1,1), (2,2,1))
            #noDeNovoPatterns = ((0,1,0),  (1,0,1), (1,1,0), (1,1,1), (2,1,1))
            deNovoPatterns = ((0,0,1), (0,1,1), (0,2,1), (0,2,0), (1,0,0), (1,2,1), (2,0,0), (2,1,0), (2,0,1))
            #deNovoPatterns = ((0,0,1), (0,1,1), (1,0,0), (2,0,0),  (2,0,1), (2,1,0))

        else:
            # TODO: Mother is never homozygous here. Need to add patterns.
            if variant.chrom == "X":
                # male child inherited from mother only
                noDeNovoPatterns = ((0,0,1), (0,1,0), (0,1,1), (1,0,0), (1,0,1), (1,1,0), (1,1,1), (1,2,0), (1,2,1))
                deNovoPatterns = ((0,2,0), (0,2,1), (1,0,0), (1,0,1))
            else:
                # male child inherited from father only
                noDeNovoPatterns = ()
                deNovoPatterns = ((0,0,1), (1,0,0))

            #noDeNovoPatterns = ((0,1,0),  (1,0,1), (1,1,0), (1,1,1))
            #deNovoPatterns = ((0,0,1), (0,1,1), (1,0,0))

    # Likelihoods
    likelihoodReference = 0
    likelihoodNoDeNovo = 0
    likelihoodDeNovo = 0

    # Likelihood of the homozygous reference genotype throughout
    likelihoodReference = 10**(childGLs[0] + fatherGLs[0] + motherGLs[0] + PRIOR_NO_DENOVO)

    # Sum the no de novo likelihood over all sets of genotypes consistent with
    # there being no de novo variant present
    for childGT,motherGT,fatherGT in noDeNovoPatterns:
        likelihoodNoDeNovo += 10**(childGLs[childGT] + fatherGLs[fatherGT] + motherGLs[motherGT] + PRIOR_SNP)
        #likelihoodNoDeNovo += 10**(childGLs[childGT] + fatherGLs[fatherGT] + motherGLs[motherGT] + PRIOR_NO_DENOVO)

    # Sum the de novo likelihood over all sets of genotypes consistent with
    # there being a de novo variant present
    for childGT,motherGT,fatherGT in deNovoPatterns:
        likelihoodDeNovo += 10**(childGLs[childGT] + fatherGLs[fatherGT] + motherGLs[motherGT] + PRIOR_DENOVO)

    # If likelihood is tiny, then this model is not supported
    logLikelihoodReference = log10(max(1e-300, likelihoodReference))
    logLikelihoodNoDeNovo = log10(max(1e-300, likelihoodNoDeNovo))
    logLikelihoodDeNovo =  log10(max(1e-300, likelihoodDeNovo))

    totalLogLikelihoodNoDeNovo = log10(10**logLikelihoodReference + 10**logLikelihoodNoDeNovo)

    #logger.info("No = %s (%s), Ref = %s (%s), TotNo = %s (%s)" %(logLikelihoodNoDeNovo, 10**logLikelihoodNoDeNovo, logLikelihoodReference, 10**logLikelihoodReference, totalLogLikelihoodNoDeNovo, 10**totalLogLikelihoodNoDeNovo))

    # This should be < 0 if the de novo is supported
    #logger.info("No = %s. Ref = %s. TotNo = %s. DeNovo = %s. BF = %s" %(logLikelihoodNoDeNovo, logLikelihoodReference, totalLogLikelihoodNoDeNovo, logLikelihoodDeNovo, totalLogLikelihoodNoDeNovo - logLikelihoodDeNovo))
    return totalLogLikelihoodNoDeNovo - logLikelihoodDeNovo

###################################################################################################

def passesBayesianFilter(variant, pedigree, sexOfChild):
    """
    Return True if the specified variant passes the de novo filter, and
    False otherwise.
    """
    childSample = pedigree['Child']
    motherSample = pedigree['Mother']
    fatherSample = pedigree['Father']

    childGenotype = variant.samples[childSample]['GT']
    fatherGenotype = variant.samples[fatherSample]['GT']
    motherGenotype = variant.samples[motherSample]['GT']

    childGLs = variant.samples[childSample]['GL']
    fatherGLs = variant.samples[fatherSample]['GL']
    motherGLs = variant.samples[motherSample]['GL']

    bayesFactor = computeBayesFactor(childGLs, fatherGLs, motherGLs, variant, sexOfChild)
    #print variant.pos, bayesFactor

    if bayesFactor < logBayesFactorThreshold:
        return True
    else:
        return False

###################################################################################################

def parseOptions(args):
    """
    Run the Platypus variant-caller, with the specified arguments
    """
    parser = optparse.OptionParser()

    # Input data and miscellaneous
    parser.add_option("--output", dest="output",  help="Output file of somatic variant calls", action='store', type='string', default="SomaticlVariants.vcf")
    parser.add_option("--input", dest="input",  help="Input file of variant calls", action='store', type='string', default=None, required=True)
    parser.add_option("--normalSampleName", dest="normalSampleName",  help="Name of sample to treat as normal", action='store', type='string', default=None, required=True)
    parser.add_option("--tumourSampleName", dest="tumourSampleName",  help="Name of sample to treat as tumour", action='store', type='string', default=None, required=True)
    parser.add_option("--sampleSex", dest="sampleSex",  help="Sex of sample", action='store', type='string', default=None, required=True)
    parser.add_option("--verbosity", dest="verbosity",  help="Verbosity", action='store', type='int', default=2, required=False)
    #parser.add_option("--refFile",dest="refFile", help="Fasta file of reference. Index must be in same directory", action='store', type='string', required=True)
    #parser.add_option("--regions", dest="regions", type="list", help = "region as comma-separated list of chr:start-end, or just list of chr, or nothing", default=None, action = 'store')

    (options, args) = parser.parse_args(args)
    return options

###################################################################################################

if __name__ == "__main__":

    options = parseOptions(sys.argv)

    logger.info("Detecting somatic mutations in file %s" %(options.input))
    logger.info("Normal sample is %s. Tumour sample is %s" %(options.normalSampleName, options.tumourSampleName))
    logger.info("Sample sex is %s." %(options.sampleSex))
    logger.info("Output will go to file %s" %(options.output))

    # Input VCF file of Platypus (or other) calls. Can be gzipped or
    # plain text.
    inVCFName = options.input
    inVCFFile = None

    if inVCFName.endswith("gz"):
        inVCFFile = gzip.open(inVCFName, 'r')
    else:
        inVCFFile = open(inVCFName, 'r')

    normalSample = options.normalSampleName
    tumourSample = options.tumourSampleName
    sampleSex = options.sampleSex
    outputFileName = options.output
    outputFile = open(outputFileName, 'w')

    nSomaticMutations = 0

    #badFilters = set(["strandBias", "Q20", "alleleBias"]) # Variants failing any of these filters are not considered
    badFilters = set(["strandBias", "Q20"]) # Variants failing any of these filters are not considered
    samples = None # All samples in the VCF

    for line in inVCFFile:

        if line.startswith("#"):

            if not line.startswith("##"):
                samples = parseSampleNamesVCF(line, normalSampleName, tumourSampleName)

            outputFile.write(line)

        else:
            assert samples is not None, "Invalid VCF header. No sample names available"

            line = line.strip()
            cols = line.split("\t")

            # Don't want to look at multi-allelic sites
            if "," in cols[4]:
                continue

            variant = Variant(line, samples)
            variant.adjustGenotypesAndLikelihoodsForPloidy(sampleSex)

            if len(variant.filters.intersection(badFilters)) > 0:
                continue

            if isSomaticCandidate(variant, sampleSex):
                outputFile.write(line + "\n")
                nSomaticMutations += 1

    # Print summary
    logger.info("Found %s somatic mutation candidates" %(nSomaticMutations))

    # Close all files
    outputFile.close()
    inVCFFile.close()
