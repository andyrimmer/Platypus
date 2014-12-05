from htslibWrapper cimport cAlignedRead

cdef list assembleReadsAndDetectVariants(char* chrom, int assemStart, int assemEnd, int refStart, int refEnd, list readBuffers, char* refSeq, options)
