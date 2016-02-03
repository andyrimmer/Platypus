"""
Utility code for Platypus. This is mostly concerned with processing needed to produce
VCF output from the various data structures used by Platypus.
"""

import logging
import vcf

cimport fastafile
cimport htslibWrapper

from cwindow cimport bamReadBuffer
from fastafile cimport FastaFile
from cgenotype cimport DiploidGenotype
from chaplotype cimport Haplotype
from variant cimport Variant
from variant cimport ASSEMBLER_VAR
from variant cimport PLATYPUS_VAR
from variant cimport FILE_VAR
from platypusutils cimport binomial
from platypusutils cimport betaBinomialCDF
from htslibWrapper cimport cAlignedRead
from htslibWrapper cimport Read_IsReverse
from htslibWrapper cimport Read_IsPaired
from htslibWrapper cimport Read_IsProperPair
from htslibWrapper cimport Read_IsDuplicate
from htslibWrapper cimport Read_IsUnmapped
from htslibWrapper cimport Read_MateIsUnmapped
from htslibWrapper cimport Read_MateIsReverse
from htslibWrapper cimport Read_IsQCFail
from htslibWrapper cimport Read_IsReadOne
from htslibWrapper cimport Read_IsSecondaryAlignment
from htslibWrapper cimport Read_SetQCFail

logger = logging.getLogger("Log")

canonicalBases = set(["A", "C", "T", "G"])

###################################################################################################

cdef extern from "math.h":
    double exp(double)
    double sqrt(double)
    double log(double)
    double log10(double)

###################################################################################################

cdef extern from "stdlib.h":
    void free(void *)
    void *malloc(size_t)
    void *calloc(size_t,size_t)
    void *realloc(void *,size_t)
    void *memset(void* ptr, int value, size_t num)
###################################################################################################

# Nasty global variables
cdef double LOG10_E = log10(exp(1.0))
cdef int CIGAR_M = 0 # Match
cdef int CIGAR_I = 1 # Insertion
cdef int CIGAR_D = 2 # Deletion
cdef int CIGAR_N = 3 # Skipped region from reference
cdef int CIGAR_S = 4 # Soft clipping. Sequence is present in read
cdef int CIGAR_H = 5 # Hard clipping. Sequence is not present in read
cdef int CIGAR_P = 6 # Padding. Used for padded alignment
cdef int CIGAR_EQ = 7 # Alignment match; sequence match
cdef int CIGAR_X = 8 # Alignment match; sequence mismatch

###################################################################################################

# Definititons for the individual and filter fields included in the header
vcfInfoSignature = {
    "FR":vcf.FORMAT('FR',1,".",'Float','Estimated population frequency of variant',-1),
    "PP":vcf.FORMAT('PP',1,".",'Float','Posterior probability (phred scaled) that this variant segregates',-1),
    "TC":vcf.FORMAT('TC',1,1,'Integer','Total coverage at this locus',-1),
    "WS":vcf.FORMAT('WS',1,1,'Integer','Starting position of calling window',-1),
    "WE":vcf.FORMAT('WE',1,1,'Integer','End position of calling window',-1),
    "TCR":vcf.FORMAT('TCR',1,1,'Integer','Total reverse strand coverage at this locus',-1),
    "TCF":vcf.FORMAT('TCF',1,1,'Integer','Total forward strand coverage at this locus',-1),
    "TR":vcf.FORMAT('TR',1,".",'Integer','Total number of reads containing this variant',-1),
    "NF":vcf.FORMAT('NF',1,".",'Integer','Total number of forward reads containing this variant',-1),
    "NR":vcf.FORMAT('NR',1,".",'Integer','Total number of reverse reads containing this variant',-1),
    "MGOF":vcf.FORMAT('MGOF',1,".",'Integer','Worst goodness-of-fit value reported across all samples',-1),
    "SC":vcf.FORMAT('SC',1,1,'String','Genomic sequence 10 bases either side of variant position',-1),
    "HP":vcf.FORMAT('HP',1,1,'Integer','Homopolymer run length around variant locus',-1),
    "BRF":vcf.FORMAT('BRF',1,1,'Float','Fraction of reads around this variant that failed filters',-1),
    "MMLQ":vcf.FORMAT('MMLQ',1,1,'Float','Median minimum base quality for bases around variant',-1),
    "QD":vcf.FORMAT('QD',1,1,'Float','Variant-quality/read-depth for this variant',-1),
    "Source":vcf.FORMAT('Source',1,".",'String','Was this variant suggested by Playtypus, Assembler, or from a VCF?',-1),
    "START":vcf.FORMAT('START',1,".",'Integer','Start position of reference call block', -1),
    "END":vcf.FORMAT('END',1,".",'Integer','End position of reference call block', -1),
    "Size":vcf.FORMAT('Size',1,".",'Integer','Size of reference call block', -1),
    "HapScore":vcf.FORMAT('HapScore',1,".",'Integer','Haplotype score measuring the number of haplotypes the variant is segregating into in a window', -1),
    "MQ":vcf.FORMAT('MQ',1,".",'Float','Root mean square of mapping qualities of reads at the variant position',-1),
    "FS":vcf.FORMAT('FS',1,".",'Float','Fisher\'s exact test for strand bias (Phred scale)',-1),
    "SbPval":vcf.FORMAT('SbPval',1,".",'Float','Binomial P-value for strand bias test',-1),
#   "MQRankSum":vcf.FORMAT('MQRankSum',1,".",'Float','Mann-Whitney Rank sum test for mapping quality difference between ref and alt (Phred scale)',-1),
    "ReadPosRankSum":vcf.FORMAT('ReadPosRankSum',1,".",'Float','Mann-Whitney Rank sum test for difference between in positions of variants in reads from ref and alt',-1),
}

vcfFilterSignature = {
    "alleleBias":vcf.FORMAT('alleleBias',1,0,'Flag','Variant frequency is lower than expected for het','.'),
    "strandBias":vcf.FORMAT('strandBias',1,0,'Flag','Variant fails strand-bias filter','.'),
    "badReads":vcf.FORMAT('badReads',1,0,'Flag','Variant supported only by reads with low quality bases close to variant position, and not present on both strands.','.'),
    "MQ":vcf.FORMAT('MQ',1,0,'Flag','Root-mean-square mapping quality across calling region is low.','.'),
    "Q20":vcf.FORMAT('Q20',1,0,'Flag','Variant quality is below 20.','.'),
    "QualDepth":vcf.FORMAT('HapScore',1,0,'Flag','Too many haplotypes are supported by the data in this region.','.'),
    "HapScore":vcf.FORMAT('QualDepth',1,0,'Flag','Variant quality/Read depth ratio is low.','.'),
    "GOF":vcf.FORMAT('GOF',1,0,'Flag','Variant fails goodness-of-fit test.','.'),
    "hp10":vcf.FORMAT('hp10',1,0,'Flag','Flanking sequence contains homopolymer of length 10 or greater','.'),
    "REFCALL":vcf.FORMAT('REFCALL',1,0,'Flag','This line represents a homozygous reference call','.'),
    "QD":vcf.FORMAT('QD',1,0,'Flag','Variants fail quality/depth filter.','.'),
    "SC":vcf.FORMAT('SC',1,0,'Flag','Variants fail sequence-context filter. Surrounding sequence is low-complexity','.'),
}

vcfFormatSignature = {
    "GT":vcf.FORMAT('GT',1,1,'String','Unphased genotypes','.'),
    "GL":vcf.FORMAT('GL',1,'.','Float','Genotype log10-likelihoods for AA,AB and BB genotypes, where A = ref and B = variant. Only applicable for bi-allelic sites','.'),
    "GQ":vcf.FORMAT('GQ',1,'.','Integer','Genotype quality as phred score','.'),
    "GOF":vcf.FORMAT('GOF',1,'.','Float','Goodness of fit value','.'),
    "NR":vcf.FORMAT('NR',1,'.','Integer','Number of reads covering variant location in this sample','.'),
    "NV":vcf.FORMAT('NV',1,'.','Integer','Number of reads containing variant in this sample','.'),
}

###################################################################################################

def outputSingleLineOfVCF(dict data, object outputFile, object vcfFile, options):
    """
    """
    #vcfDataLine = {'info':lineinfo,'format':format}
    # Write variant call info to VCF file.

    vcfFile.write_data(outputFile, data)

    #theChrom = data['chrom']
    #thePos = data['pos']
    #theId = data['id']
    #theRef = data['ref']
    #theAlts = data['alt']
    #theQual = data['qual']
    #theFilter = data['filter']

    #theLine = "\t".join( (theChrom,thePos,theId,theRef,",".join(theAlts), theQual,";".join(theFilter)) )

###################################################################################################

cdef int testGenotype(Variant var, int hapIsRefAtVarPos, set varsInHap):

    if var is None:

        if hapIsRefAtVarPos:
            return True
        else:
            return False
    else:
        if var in varsInHap:
            return True
        else:
            return False

###################################################################################################

cdef tuple computeGenotypeCallAndLikelihoods(int variantPosition, list haplotypes, list genotypes, int sampleIndex, double* haplotypeFrequencies, double** genotypeLikelihoods, double** gofValues, int** haplotypeIndexes, int** varInHap, list variantsThisPos, int* haplotypeIsRefAtThisPos, int nIndividuals, thisSample):
    """
    """
    # N.B. This may not deal correctly with multi-allelic sites, where there could
    # actually be an insertion and SNP in one haplotype, but reported at the same position
    # in VCF.

    cdef double marginalGenotypeLikelihood = 0.0
    cdef double currentLikelihood = 0.0
    cdef double sumLikelihoods = 0.0
    cdef double factor = 0.0
    cdef double bestGoodnessOfFitValue = 1e6
    cdef double haplotypeFrequency1 = 1.0 # Leave as 1.0 until we fix the EM calculations
    cdef double haplotypeFrequency2 = 1.0 # Leave as 1.0 until we fix the EM calculations
    cdef double bestLikelihood = -1.0
    cdef double nonRefPosterior = 0.0 # Sum of 1/0 and 1/1 posteriors for bi-allelic sites
    cdef double refPosterior = 0.0 # 0/0 posterior for bi-allelic sites
    cdef int bestVarIndex1 = -1
    cdef int bestVarIndex2 = -1
    cdef int index1,index2,genotypeIndex,hap1Index,hap2Index = -1
    cdef int hap1IsRefAtVarPos = -1
    cdef int hap2IsRefAtVarPos = -1

    cdef DiploidGenotype genotype = None
    cdef Variant var1,var2,variant = None
    cdef Haplotype haplotype,hap1,hap2 = None
    cdef list likelihoods = []
    cdef set varsHap1
    cdef set varsHap2

    cdef bint matchingGenotype = False
    cdef bint var1InHap1 = False
    cdef bint var1InHap2 = False
    cdef bint var2InHap1 = False
    cdef bint var2InHap2 = False

    cdef int nVariants = len(variantsThisPos)
    cdef int nGenotypes = len(genotypes)
    cdef int nHaplotypes = len(haplotypes)

    cdef int phasedIndex1 = -1
    cdef int phasedIndex2 = -1
    cdef double phasedMaxLike = -1e6

    for index1 in range(nVariants + 1):
        for index2 in range(nVariants + 1):

            if index2 > index1:
                break

            marginalGenotypeLikelihood = 0.0

            for genotypeIndex in range(nGenotypes):

                # Use pre-cached haplotype indices
                hap1Index = haplotypeIndexes[genotypeIndex][0]
                hap2Index = haplotypeIndexes[genotypeIndex][1]
                hap1IsRefAtVarPos = haplotypeIsRefAtThisPos[hap1Index]
                hap2IsRefAtVarPos = haplotypeIsRefAtThisPos[hap2Index]

                if hap1Index != hap2Index:
                    factor = 2.0
                else:
                    factor = 1.0

                # Use EM frequencies
                haplotypeFrequency1 = haplotypeFrequencies[hap1Index]
                haplotypeFrequency2 = haplotypeFrequencies[hap2Index]

                matchingGenotype = False

                # Ref/Ref genotype
                if index1 == 0 and index2 == 0:
                    if (hap1IsRefAtVarPos and hap2IsRefAtVarPos):
                        #logger.debug("hom ref genotype = %s. Pos = %s." %(genotypes[genotypeIndex], variantPosition))
                        matchingGenotype = True

                # Het genotype (always 1/0, 2/0 etc)
                elif index2 == 0:

                    var1InHap1 = varInHap[hap1Index][index1-1]
                    var1InHap2 = varInHap[hap2Index][index1-1]

                    if (hap2IsRefAtVarPos and var1InHap1) or (hap1IsRefAtVarPos and var1InHap2):
                        #logger.debug("het ref genotype = %s. Pos = %s." %(genotypes[genotypeIndex], variantPosition))
                        matchingGenotype = True

                # Both alleles are non-ref
                else:
                    var1InHap1 = varInHap[hap1Index][index1-1]
                    var1InHap2 = varInHap[hap2Index][index1-1]
                    var2InHap1 = varInHap[hap1Index][index2-1]
                    var2InHap2 = varInHap[hap2Index][index2-1]

                    if (var1InHap1 and var2InHap2) or (var2InHap1 and var1InHap2):
                        #logger.debug("hom var genotype = %s. Pos = %s." %(genotypes[genotypeIndex], variantPosition))
                        matchingGenotype = True

                # Accumulate likelihoods for matching genotypes
                if matchingGenotype:

                    # TODO: Fix this properly.
                    # Only use EM frequencies for large-ish populations.
                    if nIndividuals > 25:
                        currentLikelihood = (factor*haplotypeFrequency1*haplotypeFrequency2*genotypeLikelihoods[sampleIndex][genotypeIndex])
                    else:
                        currentLikelihood = (factor*genotypeLikelihoods[sampleIndex][genotypeIndex])

                    marginalGenotypeLikelihood += currentLikelihood

                    # Phase using maximum likelihood genotypes
                    if currentLikelihood > phasedMaxLike:
                        phasedMaxLike = currentLikelihood

                        # Hom Ref
                        if index1 == 0 and index2 == 0:
                            phasedIndex1 = index1
                            phasedIndex2 = index2

                        # Het. Make sure genotype call is phased correctly.
                        elif index2 == 0 and index1 != 0:

                            if var1InHap1:
                                phasedIndex1 = index1
                                phasedIndex2 = index2
                            elif var1InHap2:
                                phasedIndex1 = index2
                                phasedIndex2 = index1
                            # This will happen on the way
                            else:
                                pass

                        # Hom Var. Don't need to phase
                        elif index2 == index1 and index1 > 0:
                            phasedIndex1 = index1
                            phasedIndex2 = index2

                        # Multi-allelic het. Make sure call is phased correctly.
                        elif index2 > 0 and index1 > 0 and index2 != index1:

                            if var1InHap1 and var2InHap2:
                                phasedIndex1 = index1
                                phasedIndex2 = index2
                            elif var1InHap2 and var2InHap1:
                                phasedIndex1 = index2
                                phasedIndex2 = index1
                            # This will happen on the way
                            else:
                                pass

                    if gofValues[genotypeIndex][sampleIndex] < bestGoodnessOfFitValue:
                        bestGoodnessOfFitValue = gofValues[genotypeIndex][sampleIndex]

            # Pick the best variant pair for this site
            if marginalGenotypeLikelihood > bestLikelihood:

                bestLikelihood = marginalGenotypeLikelihood
                bestVarIndex1 = index1
                bestVarIndex2 = index2

            if (index1 == 1 and index2 == 0) or (index1 == 1 and index2 == 1):
                nonRefPosterior += marginalGenotypeLikelihood
            elif (index1 == 0 and index2 == 0):
                refPosterior += marginalGenotypeLikelihood

            sumLikelihoods += marginalGenotypeLikelihood
            likelihoods.append(marginalGenotypeLikelihood)

    #for genotypeIndex in range(nGenotypes):
    #    logger.debug("GL for sample %s, index %s = %s" %(sampleIndex, genotypeIndex, genotypeLikelihoods[sampleIndex][genotypeIndex]))

    return (phasedIndex1,phasedIndex2,likelihoods,bestLikelihood/sumLikelihoods,nonRefPosterior/sumLikelihoods,refPosterior/sumLikelihoods,bestGoodnessOfFitValue)

###################################################################################################

cdef void outputCallToVCF(dict varsByPos, dict vcfInfo, dict vcfFilter, list haplotypes, list genotypes, double* haplotypeFrequencies, double** genotypeLikelihoods, double** gofValues, int** haplotypeIndexes, list readBuffers, int nIndividuals, vcfFile, FastaFile refFile, outputFile, options, list allVariants, int windowStart, int windowEnd):
    """
    Output a call to a vcf file. Uses the signature of the callingModule.
    the file is opened and the header written when this object is initialised.
    The call object must have vcfINFO,vcfFILTER, genotypes, and variants methods.
    variants is a dictionary of POS-[array of variants]. Genotypes is a list of
    genotypes in the same order as the samples header of the vcf file.
    this code prints for each of the variants, the appropriate line.
    """
    cdef Haplotype refHap
    cdef Haplotype haplotype
    cdef bamReadBuffer theBuffer
    cdef cAlignedRead** startRead
    cdef cAlignedRead** endRead

    cdef dict info = vcfInfo
    cdef dict lineinfo

    cdef list readsThisIndividual
    cdef list genotypeLikelihoodsThisInd
    cdef list normalisedGLs
    cdef list FR
    cdef list PP
    cdef list NF
    cdef list NR
    cdef list TR
    cdef list GT
    cdef list linefilter

    cdef int gofThreshold = options.maxGOF
    cdef int i = 0
    cdef int variantIndex = 0
    cdef int nVariants = 0
    cdef int nVarsHap1 = 0
    cdef int nVarsHap2 = 0
    cdef int tnReadsThisIndividual = 0
    cdef int index1,index2,hapIndex = 0

    cdef double mLTOT = -0.23025850929940459 # Minus log ten over ten
    cdef double genotypePosterior = 0.0
    cdef double nonRefPosterior = 0.0
    cdef double refPosterior = 0.0
    cdef double maxLikelihood = 0.0
    cdef double gofValue = 0.0
    cdef double maxGofValue = 0.0
    cdef double minGofValue = 1e6

    cdef bytes thisSample
    cdef list sortedPositions = sorted(varsByPos.keys())
    cdef list variants
    cdef list likelihoods = None
    cdef list variantsThisPos

    cdef int** varThisPosInHap = <int**>(calloc(len(haplotypes), sizeof(int*)))
    cdef int* haplotypeIsRefAtThisPos = <int*>(calloc(len(haplotypes), sizeof(int)))
    cdef int hapIsRef = 0

    for hapIndex in range(len(haplotypes)):
        varThisPosInHap[hapIndex] = <int*>(calloc(len(allVariants), sizeof(int)))

    for POS in sortedPositions:

        #logger.info("Variants at position %s = %s" %(POS, varsByPos[POS]))

        maxGofValue = 0.0
        minGofValue = 1e6
        variants = varsByPos[POS]
        variantsThisPos = variants
        nVariants = len(variants)

        # Store boolean flag saying if each haplotype is reference at this position, for use in genotype/likelihood computation
        for hapIndex,haplotype in enumerate(haplotypes):

            haplotypeIsRefAtThisPos[hapIndex] = 1

            for variantIndex,variant in enumerate(variants):

                if variant in haplotype.variants:
                    varThisPosInHap[hapIndex][variantsThisPos.index(variant)] = 1

                    if variant.minRefPos <= POS <= variant.maxRefPos:
                        haplotypeIsRefAtThisPos[hapIndex] = 0
                else:
                    varThisPosInHap[hapIndex][variantsThisPos.index(variant)] = 0

            for variantIndex,variant in enumerate(allVariants):

                if variant in haplotype.variants:
                    if variant.minRefPos <= POS <= variant.maxRefPos:
                        haplotypeIsRefAtThisPos[hapIndex] = 0

        if options.verbosity >= 3:
            logger.debug("Haps are %s" %(haplotypes))
            logger.debug("Logging mutated sequences...")
            logger.debug("Lengths are %s" %([len(str(haplotype.getMutatedSequence())) for haplotype in haplotypes]))

            for haplotype in haplotypes:
                logger.debug(str(haplotype.getMutatedSequence())[50:-50])

            logger.debug("Unfiltered variants are %s" %(sorted(variants)))
            logger.debug("Var info = %s" %(info))
            logger.debug("Unfiltered variant info is %s" %([info[v] for v in sorted(variants)]))
            logger.debug("NVariants = %s" %(nVariants))

        chrom = variants[0].refName
        ref,alt = refAndAlt(chrom, POS, variants, refFile)
        theId = "."
        #id = "."#str(-1)
        qual = -1
        format = ['GT:GL:GOF:GQ:NR:NV']
        linefilter = []
        lineinfo = info[variants[0]]
        FR = []             # Can be generalised to all fields of variable length
        PP = []
        NF = []
        NR = []
        TR = []

        for var in variants:

            #logger.debug("Info for variant %s is %s" %(var, info[var]))
            linefilter.extend( [f for f in vcfFilter[var] if ( f in vcfFile.getfilter() )] )

            freq = info[var]['FR']
            posterior = info[var]['PP']
            varTR = info[var]['TR']
            varNF = info[var]['NF']
            varNR = info[var]['NR']

            FR.extend(freq)
            PP.extend(posterior)
            NR.extend(varNR)
            NF.extend(varNF)
            TR.extend(varTR)

        lineinfo['WS'] = [windowStart]
        lineinfo['WE'] = [windowEnd]
        lineinfo['FR'] = FR
        lineinfo['PP'] = PP
        lineinfo['NF'] = NF
        lineinfo['NR'] = NR
        lineinfo['TR'] = TR

        linefilter = list(set(linefilter))
        qual = max([int(pp) for pp in lineinfo['PP']])
        vcfDataLine = {'chrom':chrom,'pos':POS,'ref':ref,'alt':alt,'id':theId,'info':lineinfo,'filter':linefilter,'qual':qual,'format':format}

        # Now pull out the genotype informations for each sample, and compute genotype likelihoods and
        # Phred-based qualities scores.
        theVar = variants[0]
        gofsThisPos = []
        nNonRefCalls = 0

        for i from 0 <= i < nIndividuals:

            theBuffer = readBuffers[i]
            thisSample = theBuffer.sample

            # Missing genotype call. Probably due to zero coverage for this individual.
            if theBuffer.reads.windowEnd - theBuffer.reads.windowStart == 0:
                vcfDataLine[thisSample] = dict(GT=[[".", "/", "."]], GL=[0,0,0], GQ=[0],GOF=[0],NR=[0], NV=[0])
                continue

            # Have genotype. Need to construct call string ("0/1" etc), and compute likelihoods
            # and quality.
            else:
                #logger.info("nReads for sample %s = %s" %(thisSample, theBuffer.reads.windowEnd - theBuffer.reads.windowStart))
                (index1,index2,likelihoods,genotypePosterior,nonRefPosterior,refPosterior,gofValue) = computeGenotypeCallAndLikelihoods(POS, haplotypes, genotypes, i, haplotypeFrequencies, genotypeLikelihoods, gofValues, haplotypeIndexes, varThisPosInHap, variantsThisPos, haplotypeIsRefAtThisPos, nIndividuals, thisSample)

                if not (index1 == 0 and index2 == 0):
                    gofsThisPos.append(gofValue)
                    nNonRefCalls += 1

                GT = [str(index1), "/", str(index2)]
                phredPosterior = int(min(99, round(-10.0*log10(max(1e-10, 1.0-genotypePosterior)))))
                phredNonRefPosterior = int(min(99, round(-10.0*log10(max(1e-10, 1.0-nonRefPosterior)))))
                phredRefPosterior = int(min(99, round(-10.0*log10(max(1e-10, 1.0-refPosterior)))))

                if options.verbosity >= 3:
                    logger.debug("Sample %s has genotype posterior = %s. non-ref posterior = %s. Ref posterior = %s." %(thisSample,phredPosterior, phredNonRefPosterior, phredRefPosterior))

                # Don't make a call if the non-ref posterior is too low
                if nVariants == 1 and phredNonRefPosterior < options.minPosterior and phredRefPosterior < options.minPosterior:
                    GT = [".", "/", "."]
                    maxLikelihood = max(likelihoods)
                    normalisedGLs = [round(log10(max(x / maxLikelihood, 1e-300)),2) for x in likelihoods]

                elif nVariants == 1 and phredNonRefPosterior < options.minPosterior:
                    GT = ["0", "/", "0"]
                    maxLikelihood = max(likelihoods)
                    normalisedGLs = [round(log10(max(x / maxLikelihood, 1e-300)),2) for x in likelihoods]

                # Only normalise for sites which are not mult-allelic
                elif nVariants == 1:
                    maxLikelihood = max(likelihoods)
                    normalisedGLs = [round(log10(max(x / maxLikelihood, 1e-300)),2) for x in likelihoods]

                else:
                    # Multi-allelic site. Don't compute likelihoods here.
                    normalisedGLs = [-1,-1,-1]

                readsPerSample = [vcfInfo[v]['nReadsPerSample'][i] for v in variants]
                varReadsPerSample = [vcfInfo[v]['nVarReadsPerSample'][i] for v in variants]

                #vcfDataLine[thisSample] = dict(GT=[GT], GL=normalisedGLs, GQ=[phredPosterior], GOF=[min(99,int(0.5 + -10.0*log10(max(1e-10, 1.0 - round(gofValue,2)))))],NR=readsPerSample, NV=varReadsPerSample)

                # Actual nReads covering variant locus
                if nVariants == 1 and readsPerSample[0] < options.minReads:
                    vcfDataLine[thisSample] = dict(GT=[[".", "/", "."]], GL=normalisedGLs, GQ=[phredPosterior],GOF=[int(gofValue)],NR=readsPerSample, NV=varReadsPerSample)
                else:
                    vcfDataLine[thisSample] = dict(GT=[GT], GL=normalisedGLs, GQ=[phredPosterior], GOF=[int(gofValue)],NR=readsPerSample, NV=varReadsPerSample)

                if gofValue > maxGofValue:
                    maxGofValue = gofValue

                if gofValue < minGofValue:
                    minGofValue = gofValue

        # These are just temporary values. Don't output them
        vcfDataLine['info'].pop('nReadsPerSample')
        vcfDataLine['info'].pop('nVarReadsPerSample')
        vcfDataLine['info'].pop('ABPV')
        vcfDataLine['info']['MGOF'] = [int(round(maxGofValue, 2))]

        gofsThisPos.sort()

        #if minGofValue > gofThreshold:
        #if len(gofsThisPos) > 0 and gofsThisPos[ len(gofsThisPos) // 2 ] > gofThreshold:
        #    vcfDataLine['filter'].append("GOF")

        #vcfDataLine['info'].pop('SBPV')

        # Write variant call info to VCF file.
        if nNonRefCalls > 0 or options.minPosterior == 0 or options.outputRefCalls == 1:

            # Slight hack. We pad variant alleles to the left with reference sequence, when
            # they overlap with other variants, so these all get reported at the same position.
            # We need to check, at this point, if the alleles to be output actually need padding. That is,
            # we may, by now, have filtered all the overlapping variants, and so the padded variants which are
            # left should be un-padded.
            #logger.info("Variants = %s" %(variants))
            trimLeftPadding(vcfDataLine)

            # Don't output if the reference bases are non-canonical
            for c in vcfDataLine['ref']:

                if c not in canonicalBases:
                    theChrom = vcfDataLine['chrom']
                    thePos = vcfDataLine['pos']
                    theRef = vcfDataLine['ref']
                    theAlt = vcfDataLine['alt']
                    logger.warning("Skipping output for call %s:%s %s --> %s, as the reference sequence contains non-canonical bases (i.e. not A,C,T,G)" %(theChrom,thePos,theRef,theAlt))
                    break
            else:
                outputSingleLineOfVCF(vcfDataLine, outputFile, vcfFile, options)

    # Free memory. (TODO: move this somewhere else).
    for hapIndex in range(len(haplotypes)):
        free(varThisPosInHap[hapIndex])

    free(varThisPosInHap)
    free(haplotypeIsRefAtThisPos)

###################################################################################################
cdef void outputHLACallToVCF(list haplotypes, list readBuffers, int nIndividuals, FastaFile refFile, outputFile, options, int windowStart, int windowEnd):
    """
    Output a call to a vcf file. Uses the signature of the callingModule.
    the file is opened and the header written when this object is initialised.
    The call object must have genotypes, and variants methods. Genotypes is a list of
    genotypes in the same order as the samples header of the vcf file.
    this code prints for each of the variants, the appropriate line.
    """
    cdef Haplotype refHap
    cdef Haplotype haplotype
    cdef bamReadBuffer theBuffer
    cdef cAlignedRead** startRead
    cdef cAlignedRead** endRead
    cdef bytes thisSample
    cdef Variant longVar
    haplotype = haplotypes[0]
    longVar = haplotype.longVar  
 
    chrom = longVar.refName

    alts = []
    
    cdef: 
        int index = 0
        list indexOfBestGenotype=[]
        double likelihoodThisGenotype =0.0
        double maxEMLikelihood = 0.0
        DiploidGenotype gt
        int nHaplotypes = len(haplotypes)
        list gtIdx
        int hapIndex1, hapIndex2
        
        list genotypeCalls = [] # Genotype data type genotype calls of all individuals, stored by index in the "genotypes"  list
        int nReads =0           
        list allGTs =[]   # all genotypes of individuals, indexing  by the order of alts, 
        list thisGTs = [] # genotype this individual, indexing by the order in the alts list 
        list GLs = []     # genotype likelihood of the called genotype
        list NRs = []     # number of reads supporting first haplotype
        list NV1s = []     # number of reads supporting 
        list NV2s = []     # number of reads supporting 
        list CFs = []     # confidence score = (maximumLikelihood - secondMaximumLikelihood) / haplotypeLen
        list secondMaxCandidates = []
        double secondMaxLikelihood = 0.0
        int nv1, nv2, readIndex, totalReads, nBadReads, nBrokenPairReads
        double* alignScoreArr1 = NULL
        double* alignScoreArr2 = NULL
        cAlignedRead* thisRead
        cdef Haplotype hap1, hap2
        list thisNV1, thisNV2
        double like1, like2, likelihood
        double* gof = <double*>(calloc(nIndividuals+10, sizeof(double)))
        double* arr = NULL

    varSource = set()    
    logger.debug("Start output to vcf")
    ref = '' 
    memset(gof, 0, sizeof(double)* (nIndividuals+3))
    fo = None
    if options.alignScoreFile != "":
        fo = open(options.alignScoreFile, 'a')

    for index from 0<= index < nIndividuals:
        theBuffer = readBuffers[index]
        nReads = theBuffer.reads.windowEnd - theBuffer.reads.windowStart
        nBadReads = theBuffer.badReads.windowEnd - theBuffer.badReads.windowStart
        nBrokenPairReads = theBuffer.brokenMates.windowEnd - theBuffer.brokenMates.windowStart
        totalReads = nReads + nBadReads + nBrokenPairReads
        if options.alignScoreFile != "":
            fo.write("Individual\t%d\t%d\t%d:%d-%d\n" %(index, len(haplotypes), nReads,windowStart, windowEnd))
        secondMaxCandidates = []

        thisGTs = []
        NRs.append(nReads) 
        thisNV1, thisNV2 = [], []
        if nReads ==0:
            genotypeCalls.append([])
            allGTs.append(thisGTs)
            GLs.append(0.0)
            NV1s.append([])       
            NV2s.append([])
            CFs.append(0.0)
            continue

        indexOfBestGenotype = []
        maxEMLikelihood = 0.0
        likelihoodThisGenotype = -1e7
        if options.alignScoreFile != "":
            for hap1 in haplotypes:
                thisStr = hap1.shortHaplotypeSequence
                fo.write("%d %d %s\n" %(hap1.startPos+1, hap1.endPos, thisStr))

        for hapIndex1 from 0<= hapIndex1 < nHaplotypes:
            LKs = []
            hap1 = haplotypes[hapIndex1]
            arr = hap1.alignReads(index, theBuffer.reads.windowStart, theBuffer.reads.windowEnd, theBuffer.badReads.windowStart, theBuffer.badReads.windowEnd, theBuffer.brokenMates.windowStart, theBuffer.brokenMates.windowEnd,1)

        for hapIndex1 from 0<= hapIndex1 < nHaplotypes:
            LKs = []
            hap1 = haplotypes[hapIndex1]
            for hapIndex2 from 0<=hapIndex2 < nHaplotypes:
                hap2 = haplotypes[hapIndex2]
                gt = DiploidGenotype(hap1, hap2)
                likelihoodThisGenotype= gt.calculateDataLikelihood(theBuffer.reads.windowStart, theBuffer.reads.windowEnd, theBuffer.badReads.windowStart, theBuffer.badReads.windowEnd, theBuffer.brokenMates.windowStart, theBuffer.brokenMates.windowEnd, index, nIndividuals, gof, options.HLATyping)
                LKs.append(likelihoodThisGenotype)
                tempIndex = sorted([hapIndex1, hapIndex2])
                if indexOfBestGenotype == [] or likelihoodThisGenotype > maxEMLikelihood:
                    maxEMLikelihood = likelihoodThisGenotype
                    indexOfBestGenotype = [tempIndex]
                    secondMaxCandidates.append(likelihoodThisGenotype)
                elif likelihoodThisGenotype == maxEMLikelihood and tempIndex not in indexOfBestGenotype:
                    indexOfBestGenotype.append(tempIndex)
            if options.alignScoreFile != "":
                fo.write("%s\n" %("\t".join(map(str,LKs))))
        #find the secondMaxLikelihood value to calculate confidence score
        if len(secondMaxCandidates) > 1:
            secondMaxLikelihood = sorted(secondMaxCandidates, reverse = True)[1]
        else:
            secondMaxLikelihood = maxEMLikelihood -100.0
        
        for tempIdx, gtIndex in enumerate(indexOfBestGenotype):
            hap1, hap2 = haplotypes[gtIndex[0]], haplotypes[gtIndex[1]]
            alt1, alt2 = hap1.getShortHaplotypeSequence(), hap2.getShortHaplotypeSequence()
            varSource.add(hap1.longVar.varSource)
            varSource.add(hap2.longVar.varSource)
            if tempIdx ==0:
                ref = hap1.shortReferenceSequence
                       
            if alt1 != ref and alt1 not in alts:
               alts.append(alt1)
            if alt2 != ref and alt2 not in alts:
               alts.append(alt2)
            altIndex1, altIndex2 = 0 , 0
            if alt1 in alts:
                altIndex1 = alts.index(alt1) + 1
            if alt2 in alts:
                altIndex2 = alts.index(alt2) + 1
            thisGTs.append(str(altIndex1) +  "/" +str( altIndex2))
            nv1, nv2 = 0, 0
            alignScoreArr1 = hap1.alignReads(index, theBuffer.reads.windowStart, theBuffer.reads.windowEnd, theBuffer.badReads.windowStart, theBuffer.badReads.windowEnd, theBuffer.brokenMates.windowStart, theBuffer.brokenMates.windowEnd, True)            
            alignScoreArr2 = hap2.alignReads(index, theBuffer.reads.windowStart, theBuffer.reads.windowEnd, theBuffer.badReads.windowStart, theBuffer.badReads.windowEnd, theBuffer.brokenMates.windowStart, theBuffer.brokenMates.windowEnd, True)            
            for readIndex from 0<= readIndex < nReads:
                if alignScoreArr1[readIndex]==999 and alignScoreArr2[readIndex] ==999:
                    break
                score1 = -10 *  alignScoreArr1[readIndex]
                score2 = -10 *  alignScoreArr2[readIndex]
                if score1 < 5:  nv1 +=1
                if score2 < 5:  nv2 +=1
            thisNV1.append(nv1)
            thisNV2.append(nv2)
        if maxEMLikelihood!=0.0:
            confidenceScore = -(maxEMLikelihood - secondMaxLikelihood) *  (windowEnd - windowStart)/maxEMLikelihood
        elif len(haplotypes) > 1:
            confidenceScore = maxEMLikelihood - secondMaxLikelihood
        else:
            confidenceScore = 100

        genotypeCalls.append(indexOfBestGenotype)
        allGTs.append(thisGTs)
        GLs.append(maxEMLikelihood)
        NV1s.append(thisNV1)
        NV2s.append(thisNV2)
        CFs.append(confidenceScore)
    if options.alignScoreFile != "":
        fo.close()

    theId = "."
    qual = -1
    format = 'GT:GL:NR:NV1:NV2'
    linefilter = ["PASS"]
    linefilter = ""
    qual = max([int(cf) for cf in CFs])
    if alts== []:
        alts = ["."]
        linefilter = "REFCALL"
    else:
        linefilter = "PASS"

    infoLine = "WS=" + str(windowStart+1) + ";WE=" + str(windowEnd) + ";Size=" + str(windowEnd - windowStart + 1) + ";varSource=" + ",".join(map(str, list(varSource)))
    theLine = "\t".join([chrom, str(windowStart+1), theId, ref, ",".join(alts), str(qual), linefilter,  infoLine, format]) 

    # Now pull out the genotype informations for each sample, and compute genotype likelihoods and

    for i from 0 <= i < nIndividuals:
        theBuffer = readBuffers[i]
        thisSample = theBuffer.sample
        sampleLine = ":".join([",".join(allGTs[i]), str(GLs[i]), str(NRs[i]), ",".join(map(str, NV1s[i])), ",".join(map(str, NV2s[i]))])
        theLine += "\t" + sampleLine

    outputFile.write("%s\n" %(theLine))    
    free(gof)
    logger.debug('Out vcf writing')
###################################################################################################


def trimLeftPadding(vcfDataLine):
    """
    This is a small hack to make sure that no un-needed left padding is left
    in the VCF output.
    """
    ref = vcfDataLine['ref']
    alt = vcfDataLine['alt']

    #logger.info("ref = %s. alt = %s. pos = %s" %(ref, alt, vcfDataLine['pos']))

    if alt == []:
        pass
    else:
        # Only loop up to length of shortest of all ref/alt sequences
        minLenRefAlt = min(len(ref), min(map(len,alt)))
        maxLenDiff = max([abs(len(ref) - len(allele)) for allele in alt])

        for i in range(1, minLenRefAlt):

            firstAltChars = list(set(allele[0].upper() for allele in alt))
            secondAltChars = list(set(allele[1].upper() for allele in alt if len(allele) > 1))
            firstRefChar = ref[0].upper()

            #logger.info("Len first alt chars = %s. First Ref = %s. First alt = %s" %(len(firstAltChars), firstRefChar, firstAltChars[0]))
            #logger.info("Max len diff = %s. ref[1] = %s. alt[0][1] = %s" %(maxLenDiff, ref[1], alt[0][1]))

            # Stop if alt alleles don't have same first character, or don't match ref
            if len(firstAltChars) > 1 or firstRefChar != firstAltChars[0]:
                break

            # Stop if there are indels present, and the 2nd base disagrees. This is because VCF
            # represents indels with a leading reference base. N.B. Need to check the second
            # base on all alt alleles, not just the first
            elif maxLenDiff > 0 and (len(secondAltChars) > 1 or ref[1] != secondAltChars[0]):
                break

            # Trim one base from ref and all alts, and adjust position
            ref, alt = ref[1:], [allele[1:] for allele in alt]
            vcfDataLine['pos'] += 1

    #logger.info("now ref = %s. alt = %s. pos = %s" %(ref, alt, vcfDataLine['pos']))

    vcfDataLine['ref'] = ref
    vcfDataLine['alt'] = alt

###################################################################################################

cdef tuple refAndAlt(char* chrom, int POS, list variants, FastaFile refFile):
    """
    Calculates the right REF and ALT fields for a set of variants.
    In general, this is pretty difficult - it's ok if there's just
    one, but when there's more then one, it's a total nightmare.
    http://1000genomes.org/wiki/doku.php?do=show&id=1000_genomes%3Aanalysis%3Avcf4.0
    """
    cdef int nonSnpPresent = 0
    cdef int indelPresent = 0
    cdef Variant v
    cdef list ALT
    cdef list seq
    cdef bytes REF

    for v in variants:
        nAdded = v.nAdded
        nRemoved = v.nRemoved

        if nRemoved != 1 or nAdded != 1:
            nonSnpPresent = 1

            if nRemoved != nAdded:
                indelPresent = 1

    if not nonSnpPresent:
        # Just SNPs
        REF = refFile.getCharacter(chrom, POS)
        ALT = [ v.added for v in variants ]
        return REF,ALT
    else:
        # If there are indels around, then the reference base is the base before the insertion
        # or the base before the first deleted base. This means that the alternative allele for
        # snps will have two bases

        rlen = max([v.nRemoved for v in variants]) # number of extra bases in REF
        REF = None

        if indelPresent:
            REF = refFile.getSequence(chrom, POS, (POS+rlen) + 1)
        else:
            REF = refFile.getSequence(chrom, POS, (POS+rlen))

        ALT = []

        for v in variants:
            seq = list(REF)

            if v.nRemoved == v.nAdded:
                seq[0: len(v.added)] = v.added
            else:
                seq[1:1+v.nRemoved] = v.added

            ALT.append("".join(seq))

        return REF,ALT

###################################################################################################

cdef int readOverlapsVariant(cAlignedRead* pRead, int varBAMMinPos, int varBAMMaxPos):
    """
    Check if a read overlaps the variant position
    """
    cdef int readStart = pRead.pos
    cdef int readEnd = pRead.end

    if readStart <= varBAMMaxPos and readEnd > varBAMMinPos:
        #logger.info("Read %s - %s overlaps var %s - %s" %(readStart, readEnd, varBAMMinPos, varBAMMaxPos))
        return True
    else:
        #logger.info("Read %s -%s does not overlap var %s - %s" %(readStart, readEnd, varBAMMinPos, varBAMMaxPos))
        return False

###################################################################################################

cdef int readQualIsGoodVariantPosition(cAlignedRead* pRead, int varBAMMinPos, int varBAMMaxPos, int minBaseQual):
    """
    If a read overlaps the variant position, then check the qualities at that position. If quals
    are all < minBaseQual, then return False, otherwise True
    """
    if pRead.qual == NULL:
        return True

    cdef int readStart = pRead.pos
    cdef int readEnd = pRead.end

    cdef char* qualStart = pRead.qual + max(0, min(pRead.rlen, varBAMMinPos - pRead.pos))
    cdef char* qualEnd = pRead.qual + max(0, min(pRead.rlen, varBAMMaxPos - pRead.pos))

    # Conservative. I want good coverage all the way across
    while qualStart != qualEnd:

        #logger.info("Base at pos %s has qual %s. Var min = %s. Var max = %s" %(readStart + (qualStart - pRead.qual), qualStart[0], varBAMMinPos, varBAMMaxPos))

        # Ok, this is a slight fudge, as I still want to catch those cases when there's a mis-match supporting
        # a variant with a low-ish quality. I'm just not counting coverage which is very low quality.
        if qualStart[0] < 5:
            return False

        qualStart += 1

    return True

###################################################################################################

cdef int overlap(int varStart, int varEnd, int seqStart, int seqEnd):
    """
    Compute and return the number of bases by which a read overlaps the haplotype of interest.
    """
    cdef int overlapStart = max(varStart, seqStart)
    cdef int overlapEnd = min(varEnd, seqEnd)

    if overlapEnd > overlapStart:
        return overlapEnd - overlapStart
    else:
        return -1

###################################################################################################

cdef int variantSupportedByRead(cAlignedRead* pRead, int varBAMMinPos, int varBAMMaxPos, Variant variant, int countOnlyExactIndelMatches):
    """
    Return True if the specified variant is supported by the specified read
    """
    cdef int refOffset = 0
    cdef int readOffset = 0
    cdef int readStart = pRead.pos
    cdef int match = 0
    cdef int cigarIndex = 0
    cdef int cigarLength = pRead.cigarLen
    cdef int flag = 0
    cdef int length = 0
    cdef int startPosOfVarInRead = 0
    cdef int varPos = variant.refPos
    cdef int lenVarAdded = variant.nAdded
    cdef int lenVarRemoved = variant.nRemoved
    cdef char* varAdded = variant.added
    cdef char* varRemoved = variant.removed

    for cigarIndex from 0 <= cigarIndex < cigarLength:

        flag = pRead.cigarOps[2*cigarIndex]
        length = pRead.cigarOps[2*cigarIndex+1]

        # An insertion take us further along the read, but not the reference
        if flag == CIGAR_I:
            startPosOfVarInRead = readOffset

            # Read overlaps variant site and has an insertion.
            if lenVarAdded != lenVarRemoved:

                # Are we looking for exact matches?
                if countOnlyExactIndelMatches:
                    # Want exact length match
                    if variant.nAdded - variant.nRemoved == length:
                        # And exact sequence match
                        if pRead.seq[startPosOfVarInRead:startPosOfVarInRead+lenVarAdded] == varAdded:
                            return True
                    # Otherwise match is not valid
                    return False

                # If not looking for exact matches then any insertion counts
                else:
                    return True

            readOffset += length

        elif flag == CIGAR_D:
            startPosOfVarInRead = readOffset

            # Counting support if the read overlaps variant site and has an indel.
            if variant.nAdded != variant.nRemoved:

                # Are we looking for exact matches?
                if countOnlyExactIndelMatches:

                    # Want exact length match
                    if variant.nRemoved - variant.nAdded == length:
                        return True
                    # Otherwise match is not valid
                    else:
                        return False

                # If not looking for exact matches then any deletion counts
                else:
                    return True

            refOffset += length

        # A match take us further along the reference and the read
        elif flag == CIGAR_M or flag == CIGAR_EQ or flag == CIGAR_X:

            # Check if read base matches variant

            # If there are no insertions/deletions, then the variant position in the read is just the reference
            # position - the start position of the read.
            startPosOfVarInRead = varPos - readStart

            # An insertion before the variant means that it will occur later in the read by insertionSize
            startPosOfVarInRead += readOffset

            # A deletion before the variant means that it will occur earlier in the read by deletion Size
            startPosOfVarInRead -= refOffset
            #logger.info("Read start = %s. var pos = %s. read seq = %s. var added = %s" %(readStart, varPos, pRead.seq[startPosOfVarInRead:startPosOfVarInRead+lenVarAdded], varAdded))

            if refOffset + readStart <= varPos and refOffset + readStart + length > varPos and lenVarAdded == lenVarRemoved:
                if startPosOfVarInRead + lenVarAdded <= pRead.rlen:
                    if pRead.seq[startPosOfVarInRead:startPosOfVarInRead+lenVarAdded] == varAdded:
                        return True

            readOffset += length
            refOffset += length

        # Skipped region from the reference.
        elif flag == CIGAR_N:
            readOffset += length
            refOffset += length

        # Soft clipping. Sequence is present in read, but we should ignore it.
        elif flag == CIGAR_S:
            readOffset += length

            # We moved back read starts when there is soft clipping at the beginning of the read
            if cigarIndex == 0:
                refOffset += length

        # Other kinds of flag. including H (hard clipping) and P (padding)
        else:
            continue

    # No match found
    return False

###################################################################################################

cdef int computeHaplotypeScore(list genotypes):
    """
    Calculate Haplotype Score: The number of haplotypes that a variant segregate into within the window.
    Compute the likelihood of all haplotypes, perform a simple clustering of haplotype likelihoods by a hard
    boudary of 20 from the last element, then return the size of either the first cluster if dist(cluster1, cluster2) > 50
    or return the size of (cluster1 U cluster2 ) if dist <= 50
    """
    cdef DiploidGenotype genotype
    cdef Haplotype hap1,hap2
    hapScores = {}

    for genotype in genotypes:
        hap1 = genotype.hap1
        hap2 = genotype.hap2
        hapScores[hap1] = - genotype.hap1Like
        hapScores[hap2] = - genotype.hap2Like

    HapScores = sorted(hapScores.values())
    hapScoreClusters = [[HapScores[0]]]
    dist = 0

    for i in range(1, len(HapScores)):

        if HapScores[i] -  HapScores[i-1] >20:
           if len(hapScoreClusters) ==1:
              dist = HapScores[i] - HapScores[i-1]
           if len(hapScoreClusters)==2 :
              break
           hapScoreClusters.append([HapScores[i]])

        else:
           hapScoreClusters[-1].append(HapScores[i])

    HapScore = len(hapScoreClusters[0])

    if dist < 50 and dist >0:
       HapScore += len(hapScoreClusters[1])

    return HapScore

###################################################################################################

cdef dict getHaplotypeInfo(list haplotypes, dict variantPosteriors, double* haplotypeFrequencies, int nHaplotypes):
    """
    Some info is computed at the haplotype level. Sort this
    out here.
    """
    cdef Haplotype haplotype
    cdef Variant var
    cdef int hapIndex = 0
    cdef dict hapInfo
    cdef dict INFO = {}
    cdef str PP
    cdef double FR = 0.0
    cdef double posteriorProbThisVar

    # Start INFO calculations
    for hapIndex from 0 <= hapIndex < nHaplotypes:

        haplotype = haplotypes[hapIndex]
        hapINFO = haplotype.vcfINFO()

        for var,value in hapINFO.iteritems():

            # Skip low-posterior variants
            if var not in variantPosteriors:
                continue

            if var not in INFO:
                posteriorProbThisVar = variantPosteriors[var]
                PP = "%.0f" %(posteriorProbThisVar)
                FR = haplotypeFrequencies[hapIndex]
                INFO[var] = dict(HP=value['HP'], PP=[PP], FR=[FR], SC=value['SC'])
            else:
                INFO[var]['FR'][0] += haplotypeFrequencies[hapIndex]

    return INFO

###################################################################################################

cdef double computeAlleleBiasPValue(int totalReads, int variantReads):
    """
    Compute and return allele-bias p-value
    """
    cdef double pValue = 0.0

    # Frequency is higher than het.
    if totalReads > 0 and <double>variantReads / <double>(totalReads) >= 0.5:
        return 1.0

    # No reads touching variant
    elif totalReads == 0:
        return 1.0

    # Compare total variant reads to total coverage, and computed frequency
    else:
        pValue = betaBinomialCDF(variantReads, totalReads, 20, 20)
        return min(pValue, 1.0-pValue)

###################################################################################################

cdef double computeStrandBiasPValue(int nFwdReads, int nRevReads, int nFwdVarReads, int nRevVarReads):
    """
    Compute and return a binomial p-value for the number of reads supporting a
    variant on the forward/reverse strand.
    """
    cdef int alpha = 0
    cdef int beta = 0
    cdef double pVal = 0.0
    cdef double freq = 0.0
    cdef int useForward = True

    if nFwdReads == 0 or nRevReads == 0:
        return 1.0

    if nFwdReads < nRevReads:
        useForward = False

    if nFwdReads + nRevReads > 0 and nFwdVarReads + nRevVarReads > 0:

        if useForward:
            freq = <double>(nFwdReads) / <double>(nFwdReads + nRevReads)
        else:
            freq = <double>(nRevReads) / <double>(nFwdReads + nRevReads)

        #print "Freq = %s. nFwd = %s. nRev = %s. nVarFwd = %s. nVarRev = %s" %(freq, nFwdReads, nRevReads, nFwdVarReads, nRevVarReads)

        if freq < 0.5:
            alpha = 20
            beta = int((<double>alpha)/freq - alpha)

        elif freq > 0.5:
            beta = 20
            alpha = int(beta*freq/(1.0-freq))

        else:
            alpha = 20
            beta = 20

        if useForward:
            pVal = betaBinomialCDF(nFwdVarReads, nFwdVarReads+nRevVarReads, alpha, beta)
        else:
            pVal = betaBinomialCDF(nRevVarReads, nFwdVarReads+nRevVarReads, alpha, beta)

        return pVal
    else:
        return 1.0

###################################################################################################

cdef dict vcfINFO(double* haplotypeFrequencies, dict variantPosteriors, list genotypeCalls, list genotypes, list haplotypes, list readBuffers, int nHaplotypes, options, FastaFile refFile):
    """
    INFO field for vcf output, two level dictionary of variant - info field - value
    """
    cdef Variant variant
    cdef list nReadsPerSample
    cdef list nVarReadsPerSample
    cdef list listOfMinBaseQuals = []
    cdef double BRF = 0.0
    cdef int badReadsWindow = options.badReadsWindow
    cdef int verbosity = options.verbosity
    cdef int minMapQual = options.minMapQual
    cdef int minBaseQual = options.minBaseQual
    cdef int nReadsThisSample = 0
    cdef int nVarReadsThisSample = 0
    cdef int nBadReads = 0
    cdef int nGoodReads = 0
    cdef int index = 0
    cdef int varBAMMinPos = 0
    cdef int varBAMMaxPos = 0
    cdef int TC = 0
    cdef int TC_bad = 0
    cdef int TR = 0
    cdef int TC_ab = 0
    cdef int TR_ab = 0
    cdef int NR_sb = 0
    cdef int NF_sb = 0
    cdef int TCR = 0
    cdef int TCF = 0
    cdef int TCR_sb = 0
    cdef int TCF_sb = 0
    cdef int minBaseQualInWindow = 0
    cdef int windowStart = 0
    cdef int windowSize = badReadsWindow
    cdef int windowEnd = 0
    cdef int windowIndex = 0
    cdef int startPosOfVarInRead = 0
    cdef int countOnlyExactIndelMatches = options.countOnlyExactIndelMatches
    cdef float RMSMQ = 0
    cdef bint varInGenotype = False
    cdef cAlignedRead** pStartRead
    cdef cAlignedRead** pEndRead
    cdef cAlignedRead** pStartBadRead
    cdef cAlignedRead** pEndBadRead
    cdef cAlignedRead* pRead
    cdef bamReadBuffer readBuffer
    
    # Make sure to do this in the right order
    cdef int HapScore = computeHaplotypeScore(genotypes)
    cdef dict INFO = getHaplotypeInfo(haplotypes, variantPosteriors, haplotypeFrequencies, nHaplotypes)
    cdef list sortedVars = sorted(INFO.keys())
    
    for variant in INFO:

        listOfMinBaseQuals = []
        nReadsPerSample = []
        nVarReadsPerSample = []
        nGoodReads = 0
        nBadReads = 0
        RMSMQ = 0
        TC = 0
        TC_bad = 0
        TR = 0
        TC_ab = 0
        TR_ab = 0
        TCR = 0
        TCF = 0
        NR = 0
        NF = 0
        NR_sb = 0
        NF_sb = 0
        TCR_sb = 0
        TCF_sb = 0
        varBAMMinPos = variant.bamMinPos
        varBAMMaxPos = variant.bamMaxPos

        for index,genotype in enumerate(genotypeCalls):

            varInGenotype = genotype is not None and variant in genotype
            readBuffer = readBuffers[index]
            pStartRead = readBuffer.reads.windowStart
            pEndRead = readBuffer.reads.windowEnd
            pStartBadRead = readBuffer.badReads.windowStart
            pEndBadRead = readBuffer.badReads.windowEnd
            nGoodReads += (readBuffer.reads.windowEnd - readBuffer.reads.windowStart)
            nBadReads += (readBuffer.badReads.windowEnd - readBuffer.badReads.windowStart)

            nReadsThisSample = 0
            nVarReadsThisSample = 0

            while pStartBadRead != pEndBadRead:

                pRead = pStartBadRead[0]

                # Does read overlap variant site
                if not readOverlapsVariant(pRead, varBAMMinPos, varBAMMaxPos):
                    pStartBadRead += 1
                    continue

                # Does read have good quality across variant site
                if not readQualIsGoodVariantPosition(pRead, varBAMMinPos, varBAMMaxPos, minBaseQual):
                    pStartBadRead += 1
                    continue

                TC_bad += 1
                RMSMQ += (pRead.mapq*pRead.mapq)
                pStartBadRead += 1

            while pStartRead != pEndRead:

                pRead = pStartRead[0]

                # Does read overlap variant site
                if not readOverlapsVariant(pRead, varBAMMinPos, varBAMMaxPos):
                    pStartRead += 1
                    continue

                # Does read have good quality across variant site
                if not readQualIsGoodVariantPosition(pRead, varBAMMinPos, varBAMMaxPos, minBaseQual):
                    pStartRead += 1
                    continue

                nReadsThisSample += 1
                TC += 1
                RMSMQ += (pRead.mapq*pRead.mapq)

                if varInGenotype:
                    TC_ab += 1

                    if Read_IsReverse(pRead):
                        TCR_sb += 1
                    else:
                        TCF_sb += 1

                if Read_IsReverse(pRead):
                    TCR += 1
                else:
                    TCF += 1

                # Check if read supports variant.
                if variantSupportedByRead(pRead, varBAMMinPos, varBAMMaxPos, variant, countOnlyExactIndelMatches):

                    TR += 1
                    nVarReadsThisSample += 1

                    if varInGenotype:
                        TR_ab += 1

                        if Read_IsReverse(pRead):
                            NR_sb += 1
                        else:
                            NF_sb += 1

                    if Read_IsReverse(pRead):
                        NR += 1
                    else:
                        NF += 1

                    # Loop over window around variant position, and store smallest quality
                    # windowEnd is always 1 past the end

                    # Only compute MMLQ for samples which support the variant
                    if varInGenotype:
                        windowStart = max(0, varBAMMinPos - pRead.pos - (windowSize-1)//2)
                        windowEnd = min(pRead.rlen, varBAMMaxPos - pRead.pos + (windowSize-1)//2)
                        minBaseQualInWindow = 0

                        for windowIndex in range(windowStart, windowEnd):
                            if windowIndex == windowStart:
                                minBaseQualInWindow = pRead.qual[windowIndex]
                            else:
                                minBaseQualInWindow = min(minBaseQualInWindow, pRead.qual[windowIndex])
                        
                        listOfMinBaseQuals.append(minBaseQualInWindow)
                    
                    #logger.debug('minBaseQuals list %s' %(';'.join(map(str, listOfMinBaseQuals))))

                pStartRead += 1

            nReadsPerSample.append(nReadsThisSample)
            nVarReadsPerSample.append(nVarReadsThisSample)
        
        # Compute per-sample strand-bias
        INFO[variant]['ABPV'] = [round(computeAlleleBiasPValue(TC_ab, TR_ab), 2)]
        INFO[variant]['SbPval'] = [round(computeStrandBiasPValue(TCF_sb, TCR_sb, NF_sb, NR_sb), 2)]
        INFO[variant]['TR'] = [TR]
        INFO[variant]['NF'] = [NF]
        INFO[variant]['NR'] = [NR]
        
        if TR > 0:
            qual = float(INFO[variant]['PP'][0])
            if qual > 2500:
                INFO[variant]['QD'] = [options.qdThreshold + 10]
            else:
                INFO[variant]['QD'] = [(qual + (-10*log10(variant.calculatePrior(refFile))) )/TR]
                #logger.debug("Old QD = %s. new QD = %s. TR = %s. Qual = %s. Prior = %s. Var = %s" %(qual/TR, (qual + (-10*log10(variant.calculatePrior(refFile))))/TR, TR, qual, variant.calculatePrior(refFile), variant))
        else:
            #logger.warning("Something is wrong. TR == 0 for variant %s" %(variant))
            INFO[variant]['QD'] = [0]

        INFO[variant]['BRF'] = [round(nBadReads/<double>(nGoodReads+nBadReads), 2)]
        INFO[variant]['TC'] = [TC]
        INFO[variant]['TCR'] = [TCR]
        INFO[variant]['TCF'] = [TCF]

        if TC + TC_bad > 0 and RMSMQ > 0:
            INFO[variant]['MQ'] = [round(sqrt(RMSMQ/(TC+TC_bad)), 2)]
        else:
            INFO[variant]['MQ'] = [0]

        INFO[variant]['nReadsPerSample'] = nReadsPerSample
        INFO[variant]['nVarReadsPerSample'] = nVarReadsPerSample
        INFO[variant]['FR'][0] = "%1.4f" %(INFO[variant]['FR'][0])
        INFO[variant]['HapScore'] = [HapScore]

        listOfMinBaseQuals.sort()
        
        if len(listOfMinBaseQuals) > 0:
            INFO[variant]['MMLQ'] = [listOfMinBaseQuals[ len(listOfMinBaseQuals) // 2 ]]
        else:
            INFO[variant]['MMLQ'] = [100]
        
        INFO[variant]["Source"] = []
        
        if variant.varSource & PLATYPUS_VAR:
            INFO[variant]["Source"].append("Platypus")
        
        if variant.varSource & ASSEMBLER_VAR:
            INFO[variant]["Source"].append("Assembler")
        
        if variant.varSource & FILE_VAR:
            INFO[variant]["Source"].append("File")
    
    return INFO

###################################################################################################

cdef double andersonDarlingTest(list readPositions):
    """
    Compute the Anderson-Darling test statistic, to see if the positions
    within reads of variant bases are uniformly distributed.
    """
    # Sort in ascending order.
    readPositions.sort()
    nReads = len(readPositions)
    Asquared = -nReads

    for i in xrange(nReads):
        Asquared -= (1.0/nReads)*( (2.0*i - 1.0)*log(readPositions[i]) + (2.0*nReads + 1.0  - 2.0*i)*log(1.0 - readPositions[i]) )

    return Asquared

###################################################################################################

cdef double computeSCValue(sequence):
    """
    Compute the value for Platypus' SC filter. This is simply the
    fraction of the surrounding sequence which can be made of of
    any 2 bases so, e.g. a di-nucleotide repeat would have an SC
    value of 1.0, as would a homopolymer.
    """
    counter = {}

    for char in sequence:
        if char not in counter:
            counter[char] = 1
        else:
            counter[char] += 1

    cdef int lenSeq = len(sequence)
    cdef int nCharsFromTop2Bases = sum([x[0] for x in sorted( [(nChars,base) for base,nChars in counter.iteritems()], reverse=True)[0:2] ])
    cdef double SC = <double>(nCharsFromTop2Bases)/<double>(lenSeq)
    return SC

###################################################################################################

cdef dict vcfFILTER(list genotypeCalls, list haplotypes, dict vcfInfo, dict varsByPos, options):
    """
    FILTER field for vcf output, two level dictionary of variant - info field - value
    """
    cdef Variant v
    cdef dict FILTER = {}
    cdef dict hapFILTER
    cdef dict infoThisVar
    cdef int nVars = 0
    cdef int nVarsQD = 0
    cdef int nVarsHapScore = 0
    cdef int nVarsSB = 0
    cdef int nVarsAB = 0
    cdef int nVarsRMSMQ = 0
    cdef int totalVarReads = 0
    cdef int hapScore = 0
    cdef int totalReads = 0
    cdef int totalFwdReads = 0
    cdef int totalRevReads = 0
    cdef int nVarsFailingMMLQFilter = 0
    cdef int medMinQualBases = 0
    cdef double abPVal = 0.0
    cdef double sbPVal = 0.0
    cdef double BRF = 0.0
    cdef double RMSMQ = 0.0
    cdef double QD = 0.0
    cdef int bestQual = 0
    cdef int thisQual = 0
    cdef int badReadsThreshold = options.badReadsThreshold
    cdef double abThreshold = options.abThreshold
    cdef double sbThreshold = options.sbThreshold
    cdef double rmsmqThreshold = options.rmsmqThreshold
    cdef double filteredReadsFrac = options.filteredReadsFrac
    cdef double qdThreshold = options.qdThreshold
    cdef double hapScoreThreshold = options.hapScoreThreshold
    cdef double scThreshold = options.scThreshold
    cdef int failsSCFilter = 0

    for pos,varsAtPos in varsByPos.iteritems():

        nVars = len(varsAtPos)
        nVarsSB = 0
        nVarsAB = 0
        nVarsQD = 0
        nVarsHapScore = 0
        nVarsRMSMQ = 0
        nVarsFailingMMLQFilter = 0
        hapScore = 0
        QD = 0.0
        bestQual = 0

        sequence = vcfInfo[varsAtPos[0]]['SC'][0]

        if computeSCValue(sequence) > scThreshold:
            failsSCFilter = 1
        else:
            failsSCFilter = 0

        for v in varsAtPos:

            FILTER[v] = []
            infoThisVar = vcfInfo[v]
            hapScore = int(infoThisVar['HapScore'][0])
            QD = float(infoThisVar['QD'][0])
            BRF = float(infoThisVar['BRF'][0])
            RMSMQ = float(infoThisVar['MQ'][0])
            totalForwardVarReads = int(infoThisVar['NF'][0])
            totalReverseVarReads = int(infoThisVar['NR'][0])
            totalVarReads = int(infoThisVar['TR'][0])
            totalReads = int(infoThisVar['TC'][0])
            totalFwdReads = int(infoThisVar['TCF'][0])
            totalRevReads = int(infoThisVar['TCR'][0])
            medMinQualBases = int(infoThisVar.get('MMLQ', [100])[0])
            thisQual = int(infoThisVar.get('PP', [0])[0])

            if thisQual > bestQual:
                bestQual = thisQual

            if medMinQualBases < badReadsThreshold:
                nVarsFailingMMLQFilter += 1

            abPVal = float(infoThisVar['ABPV'][0])
            sbPVal = float(infoThisVar['SbPval'][0])

            if QD < qdThreshold:
                nVarsQD += 1

            if hapScore > hapScoreThreshold:
                nVarsHapScore += 1

            if totalReads > 0 and abPVal < abThreshold:
                nVarsAB += 1

            if sbPVal < sbThreshold:
                nVarsSB += 1

            if RMSMQ < rmsmqThreshold:
                nVarsRMSMQ += 1

            if failsSCFilter:
                FILTER[v].append('SC')

        for v in varsAtPos:

            if nVarsQD == nVars:
                FILTER[v].append('QD')

            if nVarsHapScore == nVars:
                FILTER[v].append('HapScore')

            if nVarsRMSMQ == nVars:
                FILTER[v].append('MQ')

            if nVarsSB == nVars:
                FILTER[v].append('strandBias')

            if nVarsAB == nVars:
                FILTER[v].append('alleleBias')

            if nVarsFailingMMLQFilter == nVars or BRF >= filteredReadsFrac:
                FILTER[v].append('badReads')

            if bestQual < 20:
                FILTER[v].append('Q20')

    return FILTER

###################################################################################################
