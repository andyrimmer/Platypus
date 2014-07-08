import cython

###################################################################################################

cdef extern from "stdlib.h":
    void free(void *)
    void *malloc(size_t)
    void *calloc(size_t,size_t)
    void *realloc(void *,size_t)

###################################################################################################

cdef extern from "math.h":
    double exp(double)
    double log(double)
    double log10(double)
    double fabs(double)
    double round(double)

###################################################################################################

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.nonecheck(False)
def lcs(char* s, char* t, int a=0, double b=1e10):
    """
    longest common substring, which includes at least one character from s[a:b].  Returns length
    and starting positions in s and t.
    """
    cdef int lenS = len(s)
    cdef int lenT = len(t)
    cdef int* l0 = <int*>(calloc(sizeof(int), lenT)) # Current row
    cdef int* l1 = <int*>(calloc(sizeof(int), lenT)) # Previous row
    cdef int* temp
    cdef int z = 0               # lcs
    cdef int starts = -1
    cdef int startt = -1
    cdef int i,j
    cdef char sc,tc

    for i in range(lenS):
        for j in range(lenT):
            sc = s[i]
            tc = t[j]

            if sc == tc:
                if i==0 or j==0:
                    if i<b:
                        l0[j] = 1
                    else:
                        l0[j] = 0
                else:
                    if i<b or l1[j-1]>0:
                        l0[j] = l1[j-1] + 1
                    else:
                        l0[j] = 0
                if l0[j] >= z and i >= a:
                    if l0[j] > z or fabs( startt + (z - lenT)/2 ) > fabs( j-z+1 + (z - lenT/2) ):
                        z = l0[j]
                        starts = i-z+1
                        startt = j-z+1
            else:
                l0[j] = 0

        temp = l1
        l1 = l0
        l0 = temp
        #l0, l1 = l1, l0

    free(l0)
    free(l1)
    return z, starts, startt

###################################################################################################
