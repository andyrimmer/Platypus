"""
This module contains various functions and classes for handling
variant information, including utilities for generating combinations
of variants, haplotypes, and genotypes.
"""

from __future__ import division

import logging
import os
import ctabix

from variant import Variant

FILE_VAR = 2 # Duplication of definition in variant.pxd

logger = logging.getLogger("Log")

###################################################################################################

class VariantCandidateReader(object):
    """
    A class to read variant information from a vcf file, and return
    a stream of variant instances.
    """
    def __init__(self, fileNames, options):
        """
        Constructor. Takes the name of a vcf file.
        """
        self.options = options
        self.vcfFiles = []

        for fileName in fileNames:

            if ".gz" not in fileName:

                logger.error("")
                logger.error("Source File %s does not look like a bgzipped VCF (i.e. name does not end in .gz)" %(fileName))
                logger.error("You must supply a VCF that has been compressed with bgzip, and indexed with tabix")
                logger.error("Compress the VCF using bgzip: bgzip %s --> %s.gz" %(fileName, fileName))
                logger.error("Remember to use the 'tabix -p vcf %s.gz' command to index the compressed file" %(fileName))
                logger.error("This should create an index file: %s.gz.tbi" %(fileName))
                logger.error("")
                raise StandardError, "Input VCF source file %s was not compressed and indexed" %(fileName)

            # Compressed VCF. Tabix will complain if there's no index.
            else:
                self.vcfFiles.append(ctabix.Tabixfile(fileName))

    def __del__(self):
        """
        Destructor. Make sure to close files.
        """
        pass

    def Variants(self, chromosome, start, end):
        """
        Generator funtion. Yields variants in order of
        genomic co-ordinate.
        """
        vcfLines = None
        varList = []
        maxSize = self.options.maxSize

        for vcfFile in self.vcfFiles:
            try:
                vcfLines = vcfFile.fetch(chromosome, start, end, parser=ctabix.asVCF())
            except Exception, e:
                logger.warning("Could not retrieve variants from source file in region %s:%s-%s. Error was %s" %(chromosome,start,end,e))
                continue

            for line in vcfLines:

                if not isValidVcfLine(line):
                    continue

                # Get the components of the VCF line
                chrom = line.contig
                pos = line.pos
                ref = line.ref
                altCol = line.alt
                alts = altCol.split(",")

                lenRef = len(ref)

                for alt in alts:
                    lenAlt = len(alt)
                    varSize = abs(lenAlt - lenRef)

                    if varSize > maxSize:
                        logger.debug("Skipping large variant of size %s in source file. Maximum allowed variant size is %s" %(varSize, maxSize))
                        continue

                    # SNP
                    if lenRef == 1 and lenAlt == 1:
                        var = Variant(chromosome, pos, ref, alt, 0, FILE_VAR)
                        varList.append(var)

                    # MNP
                    elif lenRef == lenAlt:
                        # MNPs may leading and/or trailing bases trimming
                        #var = Variant(chromosome, pos, ref, alt, 0, FILE_VAR)
                        #varList.append(var)

                        tempRef = ref
                        tempAlt = alt
                        tempPos = pos
                        removed = tempRef
                        added = tempAlt

                        # Trim leading bases
                        while len(tempRef) > 0 and len(tempAlt) > 0 and tempRef[0] == tempAlt[0]:
                            tempRef = tempRef[1:]
                            tempAlt = tempAlt[1:]
                            removed = tempRef
                            added = tempAlt
                            tempPos +=1

                        # Trim trailing bases
                        while len(tempRef) > 0 and len(tempAlt) > 0 and tempRef[-1] == tempAlt[-1]:
                            tempRef = tempRef[:-1]
                            tempAlt = tempAlt[:-1]
                            removed = tempRef
                            added = tempAlt

                        var = Variant(chromosome, tempPos, removed, added, 0, FILE_VAR)
                        varList.append(var)

                    # Anything else
                    else:
                        if self.options.longHaps ==1:
                            var = Variant(chromosome, pos, ref, alt, 0, FILE_VAR)
                            varList.append(var)
                            continue

                        # VCF4 is -1 indexed for indels, so trim off first base
                        tempRef = ref[1:]
                        tempAlt = alt[1:]
                        tempPos = pos
                        removed = tempRef
                        added = tempAlt

                        # Trim the matching bits off and shift position. This will decompose
                        # multi-variant sites into individual alleles at different positions.
                        while len(tempRef) > 0 and len(tempAlt) > 0 and tempRef[0] == tempAlt[0]:
                            tempRef = tempRef[1:]
                            tempAlt = tempAlt[1:]
                            removed = tempRef
                            added = tempAlt
                            tempPos +=1

                        # Skip weird cases for now
                        #if len(removed) != 0 and len(added) != 0:
                        #    continue
                            #logger.error("Dodgy variant found at %s:%s, with ref=%s, alt = %s" %(chrom,pos,ref,alt))
                            #logger.error("This will probably break something later on...")

                        var = Variant(chromosome, tempPos, removed, added, 0, FILE_VAR)
                        varList.append(var)

        varList = sorted(list(set(varList)))
        logger.debug("Found %s variants in region %s in source file" %(len(varList), "%s:%s-%s" %(chromosome,start,end)))
        return varList

###################################################################################################

# Performs basic validation checks on a single VCF file line.
def isValidVcfLine(line):

    chromosome = line.contig # TODO any possible checks on chromosome?
    position   = line.pos
    reference  = line.ref
    variants   = line.alt.split(",")

    try:
        if int(position) < 0:
            return False
    except ValueError:
        logger.warning("Non inetgral position at chromosome " + chromosome)
        return False

    validBases = set(['A', 'C', 'G', 'T', 'N'])

    invalidBasesInReference = set(reference) - validBases
    if len(invalidBasesInReference) > 0:
        logger.warning("Invalid reference sequence at chromosome " + chromosome)
        return False

    for variant in variants:
        invalidBasesInVariant = set(variant) - validBases
        if len(invalidBasesInVariant) > 0:
            logger.warning("Invalid alternative at chromosome " + chromosome)
            return False

    return True

###################################################################################################
