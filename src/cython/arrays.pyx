"""
Module containing simple array classes, for storing grow-able arrays
of ints and doubles. These can be used from Python or Cython code.
"""

cimport cython

###################################################################################################

cdef extern from "stdlib.h":
    void free(void *)
    void *malloc(size_t)
    void *calloc(size_t,size_t)
    void *realloc(void *,size_t)

###################################################################################################

@cython.profile(False)
cdef class IntArray:
    """
    Simple structure to wrap a raw C array, with some bounds
    checking.
    """
    def __init__(self, int size, int init):
        """
        Allocate an array of size 'size', with initial values
        'init'.
        """
        self.array = <int*> malloc(size*sizeof(int))
        self.__size = size
        self.__capacity = size

        for i in range(size):
            self.array[i] = init

    cdef int getSize(self):
        """
        Return the size of the array
        """
        return self.__size

    cdef int getLast(self):
        """
        Return the last element of the array
        """
        return self.array[self.__size - 1]

    cdef int at(self, int index):
        """
        Unchecked buffer access.
        """
        return self.array[index]

    cdef void addValue(self, int index, int value):
        """
        Unchecked buffer access.
        """
        self.array[index] += value

    cdef void append(self, int value):
        """
        Append a new value to the array, re-allocating if necessary.
        """
        cdef int* temp = NULL

        if self.__size == self.__capacity:
            temp = <int*>(realloc(self.array, 2*sizeof(int)*self.__capacity))

            if temp == NULL:
                raise StandardError, "Could not re-allocate IntArray"
            else:
                self.array = temp
                self.__capacity *= 2

        self.array[self.__size] = value
        self.__size += 1

    def __len__(self):
        """
        Return length of array
        """
        return self.__size

    def __getitem__(self, int idx):
        """
        For Python's [] operator
        """
        if idx < self.__size:
            return self.array[idx]
        else:
            raise StandardError, "DoubleArray: Out of bounds. %s > %s" %(idx, self.__size)

    def __setitem__(self, int idx, int value):
        """
        For Python's [] operator
        """
        if idx < self.__size:
            self.array[idx] = value
        else:
            raise StandardError, "DoubleArray: Out of bounds. %s > %s" %(idx, self.__size)

    def __getslice__(self, int i, int j):
        """
        This allows the user to use python's slice syntax on instances of
        the IntArray class.
        """
        return [self.__getitem__(k) for k in range(i,j)]

###################################################################################################

@cython.profile(False)
cdef class DoubleArray:
    """
    Simple structure to wrap a raw C array, with some bounds
    checking.
    """
    def __init__(self, int size, double init):
        """
        Allocate an array of size 'size', with initial values
        'init'.
        """
        self.array = <double*> malloc(size*sizeof(double))
        self.__size = size
        self.__capacity = size

        for i in range(size):
            self.array[i] = init

    cdef int getSize(self):
        """
        Return the size of the array
        """
        return self.__size

    cdef double getLast(self):
        """
        Return the last element of the array
        """
        return self.array[self.__size - 1]

    cdef double at(self, int index):
        """
        Unchecked buffer access.
        """
        return self.array[index]

    cdef void append(self, double value):
        """
        Append a new value to the array, re-allocating if necessary.
        """
        cdef double* temp = NULL

        if self.__size == self.__capacity:
            temp = <double*>(realloc(self.array, 2*sizeof(double)*self.__capacity))

            if temp == NULL:
                raise StandardError, "Could not re-allocate DoubleArray"
            else:
                self.array = temp
                self.__capacity *= 2

        self.array[self.__size] = value
        self.__size += 1

    cdef void addValue(self, int index, double value):
        """
        Unchecked buffer access.
        """
        self.array[index] += value

    def __len__(self):
        """
        Return length of array
        """
        return self.__size

    def __getitem__(self, int idx):
        """
        For Python's [] operator
        """
        if idx < self.__size:
            return self.array[idx]
        else:
            raise StandardError, "DoubleArray: Out of bounds. %s > %s" %(idx, self.__size)

    def __setitem__(self, int idx, double value):
        """
        For Python's [] operator
        """
        if idx < self.__size:
            self.array[idx] = value
        else:
            raise StandardError, "DoubleArray: Out of bounds. %s > %s" %(idx, self.__size)

    def __getslice__(self, int i, int j):
        """
        This allows the user to use python's slice syntax on instances of
        the DoubleArray class.
        """
        return [self.__getitem__(k) for k in range(i,j)]

###################################################################################################
