"""
Code for identifying variants in illumina reads, based on Gerton's haplotype realignment
algorithm and initial implementation.
"""

from __future__ import division

import multiprocessing
import variantcaller
import extendedoptparse
import sys
import os
import random
import heapq
import math
import ast
import logging
import filez
import logging.handlers
import platypusutils
import time

from variantcaller import PlatypusSingleProcess
from variantcaller import PlatypusMultiProcess
from platypusutils import Open

###################################################################################################

class FileForQueueing(object):
    """
    """
    def __init__(self, theFile, line):
        """
        Store the file, and initialise the current value
        """
        self.theFile = theFile
        self.finishedReadingFile = False
        self.heap = []

        line = line
        cols = line.strip().split("\t")
        chrom = cols[0]

        # Where possible, convert chromosome names into
        # integers for sorting. If not possible, use
        # original names.
        try:
            chrom = int(chrom.upper().strip("CHR"))
        except Exception:
            pass

        pos = int(cols[1])
        heapq.heappush(self.heap, (chrom, pos, line))

        while not self.finishedReadingFile and len(self.heap) < 100:

            try:
                line = self.theFile.next()
                cols = line.strip().split("\t")
                chrom = cols[0]

                try:
                    chrom = int(chrom.upper().strip("CHR"))
                except Exception:
                    pass

                pos = int(cols[1])
            except StopIteration:
                self.finishedReadingFile = True
                break

            heapq.heappush(self.heap, (chrom, pos, line))

        # Now take the top line
        self.chrom, self.pos, self.line = heapq.heappop(self.heap)

    def __cmp__(self, other):
        """
        Comparison function. Utilises the comparison function defined in
        the AlignedRead class.
        """
        return cmp(self.chrom, other.chrom) or cmp(self.pos, other.pos)

    def __del__(self):
        """
        Destructor
        """
        self.theFile.close()
        os.remove(self.theFile.name)

    def next(self):
        """
        Increment the iterator and yield the new value. Also, store the
        current value for use in the comparison function.
        """
        if not self.finishedReadingFile:

            try:
                line = self.theFile.next()
                cols = line.strip().split("\t")
                chrom = cols[0]

                # Where possible, convert chromosome names into
                # integers for sorting. If not possible, use
                # original names.
                try:
                    chrom = int(chrom.upper().strip("CHR"))
                except Exception:
                    pass

                pos = int(cols[1])
                heapq.heappush(self.heap, (chrom, pos, line))

            except StopIteration:
                self.finishedReadingFile = True

        if len(self.heap) != 0:
            # Now take the top line
            self.chrom, self.pos, self.line = heapq.heappop(self.heap)
        else:
            raise StopIteration

###################################################################################################

def regionSort(x, y):
    """
    Sort chromosomal regions
    """
    chrom1 = x[0]
    chrom2 = y[0]
    pos1 = int(x[1])
    pos2 = int(y[1])

    try:
        chrom1 = int(chrom1.replace("chr", ""))
        chrom2 = int(chrom2.replace("chr", ""))
    except ValueError:
        pass

    return cmp(chrom1, chrom2) or cmp(pos1, pos2)

###################################################################################################

def chromAndPosSort(x, y):
    """
    Comparison function for use in sort routines. Compares strings of the form
    chr10:0-100. Sorting is done first by chromosome, in alphabetical order, and then
    by start position in numerical order.
    """
    xChrom = x.split("_")[-1].split(":")[0]
    yChrom = y.split("_")[-1].split(":")[0]
    xStart = int(x.split(":")[1].split("-")[0])
    yStart = int(y.split(":")[1].split("-")[0])

    try:
        xChrom = int(xChrom.replace("chr", ""))
        yChrom = int(yChrom.replace("chr", ""))
    except ValueError:
        pass

    return cmp(xChrom, yChrom) or cmp(xStart, yStart)

###################################################################################################

def parseTypeFromString(value):
    """
    Parse a string representation of a variable into a true, typed, python variable
    """
    return ast.literal_eval(value)

###################################################################################################

def parsePlatypusOptionsFromVCFHeader(line):
    """
    """
    class fakeClass:
        pass

    optionsStr = line.split("=")[1].replace("{", "").replace("}","")
    theOptions = fakeClass()

    for option in optionsStr.split(","):
        name,value = option.split(":", 1)

        # Get rid of extra quotes and white-space
        name = name.strip().strip("'")
        value = value.strip()

        # Get correct type, and set attribute
        value = parseTypeFromString(value)
        setattr(theOptions, name, value)

    return theOptions

###################################################################################################

def continueCalling(args):
    """
    This function allows the user to re-start Platypus from the partially completed output of
    a previous job. This takes a single argument: the VCF file of a previous incomplete job. Platypus
    then picks up all the options for the previous job from the VCF header, and restarts calling from the latest
    sensible position (the last integer multipls of --bufferSize on the last chromosome in the VCF).
    """
    # Create a logger
    logger = logging.getLogger("ATemporaryLog")
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    ch = logging.StreamHandler()
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    ch.setLevel(logging.DEBUG)
    logger.setLevel(logging.DEBUG)


    # Seed the Python random number generator
    random.seed("Yet acquiescingly I did turn as he pointed: neither pride nor hope rekindling at the end descried, so much as gladness that some end might be.")
    parser = extendedoptparse.OptionParser()
    parser.add_option("--vcfFile", dest="vcfFile", help="Platypus will start again from the nearest possible co-ordinate to the end of this VCF. This must be a VCF produced by Platypus", action='store', type='string')
    (options, args) = parser.parse_args(args)

    newOutputFileName = options.vcfFile. replace(".vcf", "_ContinuedFromFailedProcess.vcf")

    logger.info("Platypus will now attempt to finish running a failed process, from the VCF output in file %s" %(options.vcfFile))
    logger.info("Complete output (old + new) will go to file %s" %(newOutputFileName))

    theVCF = Open(options.vcfFile, 'r')
    lastLine = None
    platypusOptions = None

    for line in theVCF:

        if "platypusOptions=" in line:
            platypusOptions = parsePlatypusOptionsFromVCFHeader(line)

        lastLine = line

    if platypusOptions is None:
        logger.error("Could not parse old platypus options from VCF %s" %(options.vcfFile))
        logger.error("Check that VCF file is a valid platypus output file")
        logger.error("Quitting now.")
        return

    cols = lastLine.strip().split("\t")

    lastChrom = cols[0]
    realLastPos = int(cols[1]) - 1
    lastPos = (realLastPos//platypusOptions.bufferSize)*platypusOptions.bufferSize

    if platypusOptions.nCPU != 1:
        logger.error("Platypus can only currently continue from single process jobs")
        logger.error("The VCF you specified was produced from a multi-process Platypus job (--nCPU != 1).")
        logger.error("Quitting now.")

    logger.info("Previous job failed at %s:%s. Job will be re-run from %s:%s" %(lastChrom,realLastPos,lastChrom,lastPos))
    allRegions = sorted(platypusutils.getRegions(platypusOptions), cmp=regionSort)
    theIndex = -1

    for index,region in enumerate(allRegions):
        if region[0] == lastChrom and region[2] == lastPos:
            theIndex = index + 1

    if theIndex == -1:
        raise StandardError, "Could not find region which was unfinished in input VCF"

    logger.info("Platypus will continue calling. Output will go to file %s." %(options.vcfFile))

    doneRegions = allRegions[:theIndex]
    doneChroms = set([x[0] for x in doneRegions if x[0] != lastChrom])

    # Reset input VCF file
    theVCF.seek(0,0)

    # Make new file to store complete output
    outputVCF = Open(newOutputFileName, "w")

    # Copy old, unfinished VCF into new VCF
    for line in theVCF:

        if line[0] == "#":
            outputVCF.write(line)
        else:
            cols = line.split("\t")
            chrom = cols[0]
            pos = int(cols[1]) - 1

            if chrom in doneChroms:
                outputVCF.write(line)

            elif chrom == lastChrom and pos < lastPos:
                outputVCF.write(line)

            else:
                break

    outputVCF.close()
    setattr(platypusOptions, "unfinishedRegions", allRegions[theIndex:])
    platypusOptions.output = newOutputFileName
    runVariantCaller(platypusOptions, continuing=True)

###################################################################################################

def mergeVCFFiles(tempFileNames, finalFileName, log):
    """
    """
    log.info("Merging output VCF file(s) into final file %s" %(finalFileName))

    # Final output file
    if finalFileName == "-":
        outputVCF = sys.stdout
    else:
        outputVCF = Open(finalFileName, 'wb')
    theHeap = []

    # Initialise queue
    for index, fileName in enumerate(tempFileNames):
        theFile = Open(fileName, 'rb')

        for line in theFile:

            # End of this file
            if line[0] == "#":
                if index == 0:
                    outputVCF.write(line)
                else:
                    continue
            else:
                theFileForQueueing = FileForQueueing(theFile, line)
                heapq.heappush(theHeap, theFileForQueueing)
                break
        # If there are no calls in the temp file, we still want to
        # remove it.
        else:
            theFile.close()
            os.remove(fileName)

    # Merge-sort the output using a priority queue
    while len(theHeap) != 0:

        # Get file from heap in right order
        nextFile = heapq.heappop(theHeap)
        outputVCF.write(nextFile.line)

        # Put file back on heap
        try:
            nextFile.next()
            heapq.heappush(theHeap, nextFile)
        except StopIteration:
            continue

    # Close final output file
    if finalFileName != "-":
        outputVCF.close()
    log.info("Finished merging VCF file(s)")

###################################################################################################

def expandPaths(options):
    if not os.path.exists(options.refFile):
        options.refFile = os.path.expanduser(options.refFile)
    
    expandedBamPaths = []
    for bamPath in options.bamFiles:
        if os.path.exists(bamPath):
            expandedBamPaths.append(bamPath)
        else:
            expandedBamPaths.append(os.path.expanduser(bamPath))
    options.bamFiles = expandedBamPaths
    
    if not os.path.exists(os.path.dirname(options.output)):
        options.output = os.path.expanduser(options.output)
    
    if options.sourceFile:
        expandedSourcePaths = []
        for sourcePath in options.sourceFile:
            if os.path.exists(sourcePath):
                expandedSourcePaths.append(sourcePath)
            else:
                expandedSourcePaths.append(os.path.expanduser(sourcePath))
        options.sourceFile = expandedSourcePaths
    
    if options.logFileName and not os.path.exists(os.path.dirname(options.logFileName)):
        options.logFileName = os.path.expanduser(options.logFileName)
    
    if options.regions and not os.path.exists(options.regions[0]) and os.path.exists(os.path.expanduser(options.regions[0])):
        options.regions[0] = os.path.expanduser(options.regions[0])
    
    if options.skipRegionsFile and not os.path.exists(options.skipRegionsFile) and os.path.exists(os.path.expanduser(options.skipRegionsFile)):
        options.skipRegionsFile = os.path.expanduser(options.skipRegionsFile)
    
    return options

def runVariantCaller(options, continuing=False):
    """
    Run the variant caller. If continuing == True, then we are picking up a failed job from
    where it left off.
    """
    
    options = expandPaths(options)
    
    # Seed the Python random number generator
    random.seed("Full many a flower is born to blush unseen and waste its sweetness on the desert air")

    # Set up basic logging

    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

    log = logging.getLogger('Log')

    fh = None
    ch = logging.StreamHandler()

    if continuing:
        fh = logging.FileHandler(options.logFileName, 'a')
    else:
        fh = logging.FileHandler(options.logFileName, 'w')

    ch.setFormatter(formatter)
    fh.setFormatter(formatter)

    if options.verbosity == 0:
        log.setLevel(logging.DEBUG)
        ch.setLevel(logging.ERROR)
        fh.setLevel(logging.DEBUG)
    elif options.verbosity == 1:
        log.setLevel(logging.DEBUG)
        ch.setLevel(logging.WARNING)
        fh.setLevel(logging.DEBUG)
    elif options.verbosity == 2:
        log.setLevel(logging.DEBUG)
        ch.setLevel(logging.INFO)
        fh.setLevel(logging.DEBUG)
    elif options.verbosity >= 3:
        # Debug goes to file only.
        log.setLevel(logging.DEBUG)
        ch.setLevel(logging.INFO)
        fh.setLevel(logging.DEBUG)
    else:
        raise StandardError, "Value of 'verbosity' input parameter must be between 0 and 3 inclusive"

    log.addHandler(ch)
    log.addHandler(fh)

    if continuing:
        log.info("Continuing variant calling from where we left off.")
    else:
        log.info("Beginning variant calling")
    
    log.info("Output will go to %s" %(options.output))

    regions = None

    if continuing:
        regions = options.unfinishedRegions
    else:
        regions = sorted(platypusutils.getRegions(options), cmp=regionSort)

    # Always create process manager even if nCPU=1, so that we can listen for signals from the main thread
    fileNames = set()
    processes = []
    regionsForEachProcess = []

    # In this case, create all the BAM files here, before splitting into separate processes. The files will be left open until
    # the end of the parent process, and all child processes will share the same open files via pointers.
    bamFileNames = None
    samples = None
    samplesByID = None
    samplesByBAM = None
    bamFiles = None
    theLocks = None

    for i in range(options.nCPU):
        regionsForEachProcess.append([])

    for index,region in enumerate(regions):
        regionsForEachProcess[index % options.nCPU].append(region)
    
    if options.nCPU == 1 and options.output == "-":
        processes.append(PlatypusMultiProcess("-", options, regionsForEachProcess[0]))
    else:
        for index in range(options.nCPU):
            fileName = options.output + "_temp_%s" %(index)
            fileNames.add(fileName)
            processes.append(PlatypusMultiProcess(fileName, options, regionsForEachProcess[index]))

    for process in processes:
        process.start()

    # listen for signals while any process is alive
    while True in [process.is_alive() for process in processes]:
        try:
            time.sleep(1)
        except KeyboardInterrupt:
            print "KeyboardInterrupt detected, terminating all processes..."
            for process in processes:
                process.terminate()
            log.error("Variant calling aborted due to keyboard interrupt")
            sys.exit(1)

    # make sure all processes are finished
    for process in processes:
        process.join()

    # Final output file
    if options.output != "-":
        mergeVCFFiles(fileNames, options.output, log)

    # All done. Write a message to the log, so that it's clear when the
    # program has actually finished, and not crashed.
    log.info("Finished variant calling")

###################################################################################################

def callVariants(args):
    """
    Run the Platypus variant-caller, with the specified arguments
    """
    parser = extendedoptparse.OptionParser()
    
    # Input data and miscellaneous
    parser.add_option("-o", "--output", dest="output",  help="Output SNP data file", action='store', type='string', default="AllVariants.vcf")
    parser.add_option("--refFile",dest="refFile", help="Fasta file of reference. Index must be in same directory", action='store', type='string', required=True)
    parser.add_option("--regions", dest="regions", type="list", help = "region as comma-separated list of chr:start-end, or just list of chr, or nothing", default=None, action = 'store')
    parser.add_option("--skipRegionsFile", dest="skipRegionsFile", type="string", help = "region as comma-separated list of chr:start-end, or just list of chr, or nothing", default=None, action = 'store')
    parser.add_option("--bamFiles", dest="bamFiles", type="list", help = "Comma-delimited list of bam or cram file names", default=None, required=True)
    parser.add_option("--bufferSize", dest="bufferSize", type="int", help = "Data will be buffered in regions of this size", default=100000, required=False)
    parser.add_option("--minReads", dest="minReads", help="Minimum number of supporting reads required before a variant candidate will be considered.", action='store', type='int', default=2)
    parser.add_option("--maxReads", dest="maxReads", help="Maximium coverage in window", action='store', type='float', default=5000000)
    parser.add_option("--verbosity", dest="verbosity", help="Level of logging", action='store', type='int', default=2)
    parser.add_option("--maxReadLength", dest="rlen", help="Maximum read length", action='store', type = 'int', default=150)
    parser.add_option("--logFileName", dest="logFileName", help="Name of log file", action='store', type='string', default="log.txt")
    parser.add_option("--source", dest="sourceFile", help="vcf file(s) to get candidates from", action='store', type='list', default=None)
    parser.add_option("--nCPU", dest="nCPU", help="Number of processors to use", action='store', type='int', default=1)
    parser.add_option("--parseNCBI", dest="parseNCBI", help="", type=int, action='store', default=0)
    parser.add_option("--longHaps", dest="longHaps", help="If this is set to 1, then don't trim replacement variants from input VCFs.", type='int', action='store', default=0)
    parser.add_option("--alignScoreFile", dest="alignScoreFile", help="If this is set to a string, then alignment scores of reads to haplotypes will be writen to this file. This only work when --HLATyping flag is on", type='string', action='store', default="")
    parser.add_option("--HLATyping", dest="HLATyping", help="If this is set to 1, then run HLA genotyping mode which require a source file containing HLA haplotypes", type='int', action='store', default=0)
    parser.add_option("--compressReads", dest="compressReads", help="If this is set to 1, then all reads will be compressed, and decompressd on demand. This will slow things down, but reduce memory usage.", type='int', action='store', default=0)
    parser.add_option("--qualBinSize", dest="qualBinSize", help="This sets the granularity used when compressing quality scores. If > 1 then quality compression is lossy", type='int', action='store', default=1)
    parser.add_option("--fileCaching", dest="fileCaching", help="Sets file caching level. 0: BAM/CRAM files cached. 1: CRAM files cached. 2: No file caching.", type=int, action='store', default=0)
    
    # Calling Parameters
    parser.add_option("--maxSize", dest="maxSize", help="Largest variant to consider", action='store', type='int', default=1500)
    parser.add_option("--largeWindows", dest="largeWindows", help="If set to 1, window size can be up to 'maxSize'", action='store', type='int', default=0)
    parser.add_option("--maxVariants", dest="maxVariants", help="Maximium variants to consider in a given window", action='store', type='int', default=8)
    parser.add_option("--coverageSamplingLevel", dest="coverageSamplingLevel", help="Downsample to this level of coverage when filtering haplotypes in divergent regions.", action='store', type='int', default=30)
    parser.add_option("--maxHaplotypes", dest="maxHaplotypes", help="Maximium haplotypes to consider in a given window", action='store', type='int', default=50)
    parser.add_option("--skipDifficultWindows", dest="skipDifficultWindows", help="If set to 1, skip windows with > maxVariants candidates", action='store', type='int', default=0)
    parser.add_option("--getVariantsFromBAMs", dest="getVariantsFromBAMs", help="If set to TRUE (default), variant candidates will be generated from BAMs as well as any other inputs", action='store', type='int', default=1)
    parser.add_option("--genSNPs", dest="genSNPs", help="If set to TRUE (default), SNP candidates will be considered", action='store', type='int', default=1)
    parser.add_option("--genIndels", dest="genIndels", help="If set to TRUE (default), Indel candidates will be considered", action='store', type='int', default=1)
    parser.add_option("--mergeClusteredVariants", dest="mergeClusteredVariants", help="If set to 1, variant-containing windows which are close together will be merged, resulting in slower, more accurate variant calls in diverse regions", action='store', type='int', default=1)
    parser.add_option("--minFlank", dest="minFlank", help="Ignore base-changes closer than minFlank bases to the end of reads. Also, merge SNPs within this distance into MNPs or complex replacements", action='store', type = 'int', default=10)
    parser.add_option("--trimReadFlank", dest="trimReadFlank", help="Set base-qualities to 0 within 'trimReadFlank' bases of the end of reads", action='store', type = 'int', default=0)
    parser.add_option("--filterVarsByCoverage", dest="filterVarsByCoverage", help="If 1, Platypus filters variants in difficult regions by the number of times each variant is seen.", action='store', type='int', default=1)
    parser.add_option("--filteredReadsFrac", dest="filteredReadsFrac", help="If > this fraction of reads are filtered in a given window, the 'badReads filter is triggered.", action='store', type='float', default=0.7)
    parser.add_option("--maxVarDist", dest="maxVarDist", help="Max distance between variants to be considered in the same window", action='store', type='int', default=15) # 9 is 1 base longer than the max possible alignment shift
    parser.add_option("--minVarDist", dest="minVarDist", help="Min distance allowed between windows", action='store', type='int', default=9) # 9 is 1 base longer than the max possible alignment shift
    parser.add_option("--useEMLikelihoods", dest="useEMLikelihoods", help="If 1, likelihoods computed from EM algorithm will be used to call genotypes for each sample, otherwise likelihoods from individual sample will be used.", action='store', type='int', default=0)
    parser.add_option("--countOnlyExactIndelMatches", dest="countOnlyExactIndelMatches", help="If 1, only exactly matching indels will be counted in the NV field", action='store', type='int', default=0)
    parser.add_option("--calculateFlankScore", dest="calculateFlankScore", help="If 1, an additional alignment routine is used to calculate scores from flanks outside windows (EXPERIMENTAL).", action='store', type='int', default=0)
    
    # Assembly parameters
    parser.add_option("--assemble", dest="assemble", help="If 1, Cortex will be used to assemble variant candidates for Platypus to call.", action='store', type='int', default=0)
    parser.add_option("--assembleAll", dest="assembleAll", help="If 1 then Platypus will assemble all regions.'.", action='store', type='int', default=1)
    parser.add_option("--assemblyRegionSize", dest="assemblyRegionSize", help="Size of region to assemble with Cortex", action='store', type='int', default=1500)
    parser.add_option("--assembleBadReads", dest="assembleBadReads", help="If 1, then use filtered 'bad' reads for local assembly", action='store', type='int', default=1)
    parser.add_option("--assemblerKmerSize", dest="assemblerKmerSize", help="Kmer size to use for cortex assembly'.", action='store', type='int', default=15)
    parser.add_option("--assembleBrokenPairs", dest="assembleBrokenPairs", help="If 1, then use broken read pairs for local assembly", action='store', type='int', default=0)
    parser.add_option("--noCycles", dest="noCycles", help="If 1, then don't allow cycles in the graph", action='store', type='int', default=0)
    
    # QC Parameters
    parser.add_option("--minMapQual", dest="minMapQual", help="Minimum mapping quality of read. Any reads with map qual below this are ignored", action='store', type = 'int', default=20, required=False)
    parser.add_option("--minBaseQual", dest="minBaseQual", help="Minimum allowed base-calling quality. Any bases with qual below this are ignored in SNP-calling", action='store', type = 'int', default=20, required=False)
    parser.add_option("--minGoodQualBases", dest="minGoodQualBases", help="Min bases per read that must have base-quality >= 20.", action='store', type = 'int', default=20, required=False)
    parser.add_option("--filterDuplicates", dest="filterDuplicates", help="If set to 1, duplicate reads will be removed based on the read-pair start and end", action='store', type = 'int', default=1, required=False)
    parser.add_option("--filterReadsWithUnmappedMates", dest="filterReadsWithUnmappedMates", help="If set to 1, reads with un-mapped mates will be removed", action='store', type = 'int', default=1, required=False)
    parser.add_option("--filterReadsWithDistantMates", dest="filterReadsWithDistantMates", help="If set to 1, reads with mates mapped far away will be removed", action='store', type = 'int', default=1, required=False)
    parser.add_option("--filterReadPairsWithSmallInserts", dest="filterReadPairsWithSmallInserts", help="If set to 1, read pairs with insert sizes < one read length will be removed", action='store', type = 'int', default=1, required=False)
    parser.add_option("--trimOverlapping", dest="trimOverlapping", help="If set to 1, overlapping paired reads have overlap set to qual 0", action='store', type = 'int', default=1, required=False)
    parser.add_option("--trimAdapter", dest="trimAdapter", help="If set to 1, then sets to qual 0 any part of read which exceeds the mapped fragment length. This is mainly useful for trimming adapter sequences", action='store', type = 'int', default=1, required=False)
    parser.add_option("--trimSoftClipped", dest="trimSoftClipped", help="If set to 1, then sets to qual 0 any soft clipped parts of the read.", action='store', type = 'int', default=1, required=False)
    
    # Variant-calling Filter Parameters
    parser.add_option("--maxGOF", dest="maxGOF", help="Max allowed value for goodness-of-fit test. Higher than this triggers GOF filter (Phred-scaled).", action='store', type='int', default=30)
    parser.add_option("--minPosterior", dest="minPosterior", help="Only variants with posterior >= this will be outpu to the VCF. Value is a Phred-score.", action='store', type='int', default=5)
    parser.add_option("--sbThreshold", dest="sbThreshold", help="P-value for strand-bias filtering..", action='store', type='float', default=1e-3)
    parser.add_option("--scThreshold", dest="scThreshold", help="Cut-off for SC filter.", action='store', type='float', default=0.95)
    parser.add_option("--abThreshold", dest="abThreshold", help="P-value for allele-bias filtering..", action='store', type='float', default=1e-3)
    parser.add_option("--minVarFreq", dest="minVarFreq", help="Variants below this frequency will be flagged as allele-biased", action='store', type='float', default=0.05)
    parser.add_option("--badReadsWindow", dest="badReadsWindow", help="Size of window around variant to look for low-quality bases.", action='store', type='int', default=11)
    parser.add_option("--badReadsThreshold", dest="badReadsThreshold", help="Variants where the median minimum quality in a window of badReadsWindow around the variant position falls below this value will be filtered with the flag 'badReads'.", action='store', type='int', default=15)
    parser.add_option("--rmsmqThreshold", dest="rmsmqThreshold", help="RMSMQ filter triggers when root-mean-square mapping quality across region containing variant is below this.", action='store', type='int', default=40)
    parser.add_option("--qdThreshold", dest="qdThreshold", help="QD filter triggers quality/depth for variant is below this.", action='store', type='int', default=10)
    parser.add_option("--hapScoreThreshold", dest="hapScoreThreshold", help="HapScore filter triggers HapScore for variant is above this.", action='store', type='int', default=4)
    
    # Genome VCF parameters
    parser.add_option("--outputRefCalls", dest="outputRefCalls", help="If 1, output block reference calls.", action='store', type='int', default=0)
    parser.add_option("--refCallBlockSize", dest="refCallBlockSize", help="Max size of reference call block.", action='store', type='int', default=1000)

    (options, args) = parser.parse_args(args)
    runVariantCaller(options)

###################################################################################################
