"""
This module contains utilities for generating a list windows, and the reads
and variant candidates they contain, for processing by the Platypus variant
caller.
"""
from __future__ import division

import logging
import variant
import variantutils

from variant import Variant,VariantCandidateGenerator

logger = logging.getLogger("Log")

###################################################################################################

class WindowGenerator(object):
    """
    This class is used to generate a list of all possible variants for a given region.
    """
    def __init__(self):
        """
        Constructor. Currently does nothing.
        """
        pass

    def getVariantsByPos(self, chromosome, start, end, sortedVariants):
        """
        Returns a list of lists of variant candidates. Each list contains all candidates
        at that position. Only variants which start within the region are allowed.
        """
        varsByPos = {}

        for v in sortedVariants:
            if v.refName == chromosome and v.refPos >= start and v.refPos < end:
                try:
                    varsByPos[v.refPos].append(v)
                except KeyError:
                    varsByPos[v.refPos] = [v]

        listOfVarsByPos = []

        for pos in sorted(varsByPos.keys()):
            listOfVarsByPos.append(varsByPos[pos])

        return listOfVarsByPos

    def getBunchesOfInteractingVariants(self, varsByPos, options):
        """
        Go through the list of lists or variants sorted by position, and
        concatenate neighbouring lists, if any of the variants in those
        lists interact.
        """
        bunchedVars = []
        maxWindowSize = options.rlen

        for varList in varsByPos:

            if len(bunchedVars) == 0:
                bunchedVars.append(varList)
            else:
                minLastPos = min([x.minRefPos for x in bunchedVars[-1]])
                maxLastPos = max([x.maxRefPos for x in bunchedVars[-1]])
                minThisPos = min([x.minRefPos for x in varList])
                maxThisPos = max([x.maxRefPos for x in varList])

                # Always put interacting variants in the same window
                if maxLastPos >= minThisPos:
                    bunchedVars[-1].extend(varList)

                # Try to merge windows, if edge variants are close together, and resulting window
                # size is not too large, and the mergeClusteredVariants option is set
                elif options.mergeClusteredVariants:

                    thisWindowSize = maxThisPos - minLastPos

                    # Limit the window size. Normally this is limited to 1 read length, but
                    # can be optionally made much larger (max variant size is ~1kb).
                    maxWindowSize = None

                    if options.largeWindows == 1:
                        maxWindowSize = options.maxSize
                    else:
                        maxWindowSize = options.rlen

                    edgeVarDist = minThisPos - maxLastPos
                    nVarsThisWindow = len(varList)
                    nVarsLastWindow = len(bunchedVars[-1])

                    # If adjacent windows are close together, try to merge them
                    if edgeVarDist < options.maxVarDist:

                        # I want to keep windows <= one read length
                        if thisWindowSize <= maxWindowSize:

                            # If there aren't too many variants, then just add to the existing window
                            if nVarsLastWindow + nVarsThisWindow <= options.maxVariants:
                                bunchedVars[-1].extend(varList)

                            # Too many variants.
                            else:

                                # If distance to next variant is big enough (probably >= 9), then split
                                if edgeVarDist >= options.minVarDist:
                                    bunchedVars.append(varList)

                                # Otherwise, merge the two windows
                                else:
                                    bunchedVars[-1].extend(varList)

                        # Large window. Don't split, as Cortex will produce large variants
                        else:
                            #logger.info("Splitting due to large window. edge var dist = %s" %(edgeVarDist))
                            bunchedVars.append(varList)
                            #bunchedVars[-1].extend(varList)

                    # Variants are sufficiently far apart. Split.
                    else:
                        bunchedVars.append(varList)

                # If "--mergeClusteredVariants==0" then always start a new window if variants do not
                # directly overlap
                else:
                    bunchedVars.append(varList)

        return bunchedVars

    def getWindowVariants(self, chromosome, start, end, sortedVariants, options):
        """
        Iterate over the list of bunches of interacting vars, and concatenate
        neighbouring bunches, if the total number of variants in the two neighbouring
        bunches is <= maxVarsInWindow. Return a list of lists, where each sub-list
        is a set of variants to be considered in the same window.
        """
        varsByPos = self.getVariantsByPos(chromosome, start, end, sortedVariants)
        bunchedVars = self.getBunchesOfInteractingVariants(varsByPos, options)
        return bunchedVars

    def WindowsAndVariants(self, chromosome, start, end, maxContigPos, sortedVariants, options):
        """
        Generator: return a series of dictionaries, each containing the start and end
        positions of a haplotype window, plus a list of all the variants in that window.
        Basically takes a stream of variant candidates, and figures out where to break
        them up into windows and then gives that to the haplotype generator.

        Yields a dictionary.

        The basic procedure is as follows: loop through all variants. If the next variant is outside the
        window, then either extend the window, or, if we already have enough variants in the current
        window, then yield this window, and start a new one which will contain the next variant.

        The maximum allowed number of variants in a particular window is currently set to 3.
        Once we hit the end of the variants, we catch the StopIteration, and return whatever is in the the variant buffer.
        """
        windowVars = self.getWindowVariants(chromosome, start, end, sortedVariants, options)

        # Debug output
        if options.verbosity >= 3:
            varsThisRegion = [v for v in sortedVariants if v.refPos >= start and v.refPos <= end]
            nVar = len(varsThisRegion)
            nSnp = len([v for v in varsThisRegion if abs(v.nAdded - v.nRemoved) == 0])
            logger.debug("There are %s vars in total in this region" %(nVar))
            logger.debug("There are %s SNPs and %s indels" %(nSnp, nVar-nSnp))
            logger.debug("There will be %s windows used in this region" %(len(windowVars)))

        nWindows = len(windowVars)

        for index,varsThisWindow in enumerate(windowVars):

            winStart = max(min([v.minRefPos for v in varsThisWindow]) - options.minVarDist, start)

            # Don't go past the end of the chromosome/contig. Don't mind here if we overlap slightly into
            # the next region, as this allows better calling at region boundaries. In any case we should
            # never have variants being reported in more than one region, even if called independently, as
            # only variants which start inside the region are considered for calling.
            winEnd = min(max([v.maxRefPos for v in varsThisWindow]) + options.minVarDist, maxContigPos)


            # Check for gaps between windows and, if required, yield windows spanning the gaps
            # and containing no variants
            if options.outputRefCalls:

                if index == 0:
                    firstVarPos = max(min([v.minRefPos for v in varsThisWindow]) + 1, start)

                    if firstVarPos - start >= 1:

                        for refBlockStart in xrange(start, firstVarPos, options.refCallBlockSize):
                            refBlockEnd = min(refBlockStart + options.refCallBlockSize, firstVarPos-1)

                            if refBlockStart == refBlockEnd:
                                continue

                            thisWindow = dict(chromosome=chromosome,startPos=refBlockStart, endPos=refBlockEnd,  variants=[], nVar=0)

                            if options.verbosity >= 3:
                                logger.debug("Window = %s" %(thisWindow))

                            yield thisWindow

                else:
                    lastVarPos = max([v.maxRefPos for v in windowVars[index-1]])
                    nextVarPos = min([v.minRefPos for v in varsThisWindow]) + 1 # To account for VCF 1-indexing vs Variant 0-indexing

                    if nextVarPos - lastVarPos > 1:

                        for refBlockStart in xrange(lastVarPos+1, nextVarPos, options.refCallBlockSize):
                            refBlockEnd = min(refBlockStart + options.refCallBlockSize, nextVarPos-1)

                            if refBlockStart == refBlockEnd:
                                continue

                            thisWindow = dict(chromosome=chromosome,startPos=refBlockStart, endPos=refBlockEnd,  variants=[], nVar=0)

                            if options.verbosity >= 3:
                                logger.debug("Window = %s" %(thisWindow))

                            yield thisWindow

            # Now yield the window which actually contains variants
            thisWindow = dict(chromosome=chromosome,startPos=winStart, endPos=winEnd, variants=varsThisWindow, nVar=len(varsThisWindow))

            if options.verbosity >= 3:
                logger.debug("")
                logger.debug("#####################################################################")
                logger.debug("Current window spans %s:%s-%s and has %s variants" %(chromosome, winStart, winEnd, len(varsThisWindow)))
                logger.debug("Printing all variants in window...")
                for v in varsThisWindow:
                    logger.debug(v)
                logger.debug("")

            #if thisWindow['endPos'] <= thisWindow['variants'][-1].maxRefPos:
            #    raise StandardError, "Shit. Problem with window %s" %(thisWindow)

            #logger.info(thisWindow)
            #logger.debug("Window size = %s" %(thisWindow['endPos'] - thisWindow['startPos']))
            yield thisWindow

###################################################################################################
