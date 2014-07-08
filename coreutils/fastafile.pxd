cdef class FastaIndex:
    cdef object theFile
    cdef dict getRefs(self, int parseNCBI)

cdef class FastaFile(object):
    cdef object theFile
    cdef FastaIndex theIndex
    cdef dict refs
    cdef bytes cache
    cdef bytes cacheRefName
    cdef long long int cacheStartPos
    cdef long long int cacheEndPos
    cdef bytes getCharacter(self, bytes seqName, long long int pos)
    cdef bytes getSequence(self, bytes seqName, long long int beginPos, long long int endPos)
    cdef void setCacheSequence(self, bytes seqName, long long int beginPos, long long int endPos)
