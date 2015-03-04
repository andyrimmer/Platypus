cimport fastafile
cimport chaplotype

from fastafile cimport FastaFile
from chaplotype cimport Haplotype
from variant cimport Variant
from cwindow cimport bamReadBuffer

cdef list filterVariants(list varList, FastaFile refFile, int maxReadLength, int minSupport, int maxDiff, int verbosity, options)
cdef void filterVariantsInWindow(dict thisWindow, bytes chrom, int windowStart, int windowEnd, FastaFile refFile, options, list variants, Haplotype refHaplotype, list readBuffers)
cdef list getFilteredHaplotypes(dict thisWindow, bytes chrom, int windowStart, int windowEnd, FastaFile refFile, options, list variants, Haplotype refHaplotype, list readBuffers)
cdef void filterVariantsByCoverage(dict thisWindow, bytes chrom, int windowStart, int windowEnd, FastaFile refFile, options, list variants, Haplotype refHaplotype, list readBuffers)
cdef list getHaplotypesInWindow(dict window, int nReads, FastaFile refFile, int maxCoverage, int minMapQual, int minReadQual, int maxHaplotypes, int maxVariants, int maxReadLength, int verbosity, list readBuffers, options)
cdef list padVariants(list sortedVariants, FastaFile refFile, bytes chrom)
cdef double computeVariantReadSupportFrac(Variant variant, bamReadBuffer readBuffer)

cdef list getAllHLAHaplotypesInRegion(bytes chrom, int windowStart, int windowEnd, FastaFile refFile, options, list variants, Haplotype refHaplotype, list readBuffers)
cdef Variant normaliseVar(Variant v)
cdef Variant trimLongVar(Variant v, int windowStart, int windowEnd)
cdef list getAllAssemblerHaplotypesInRegion(bytes chrom, int windowStart, int windowEnd, FastaFile refFile, options, list variants, Haplotype refHaplotype, list readBuffers)