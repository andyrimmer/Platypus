"""
The histogram module contains utilty functions and classes for histogramming
data.
"""

from bisect import bisect

cimport cython
cimport arrays

from arrays cimport IntArray
from arrays cimport DoubleArray

###################################################################################################

def make_label(label, lane, mate):
    """
    Create and return a string label to go on a line of the QC
    output file. This will consist of the lane number, the number
    of the mate pair (1 or 2), a flag to say if this value should be
    filtered (always 'N' at the moment), and the name of the value stored
    on this line.
    """
    filtered = 'N'
    return "%s\t%s\t%s\t%s" % (lane, mate, filtered, label)

###################################################################################################

def output(stream, label, value, lane, mate):
    """
    Output a value to the specified stream.
    """
    stream.write("%s\t%s\n" % (make_label(label, lane, mate), value))

###################################################################################################

def row_output(stream, label, counts, format, lane, mate, selection=None):
    """
    Output a row to the specified stream.
    """
    r = make_label(label, lane, mate)

    if selection is None:
        for i in range(len(counts)):
            r += "\t" + (format %  counts[i])
    else:
        for i in range(len(counts)):
            if i in selection:
                r += "\t" + (format %  counts[i])

    stream.write(r + "\n")

###################################################################################################

def histogram_output(stream, label, bins, counts, sub1, sub2, lane, mate):
    """
    Output a histogram to the specified stream.
    """
    row_output(stream, label+sub1, bins, "%s", lane, mate)
    row_output(stream, label+sub2, counts, "%i", lane, mate)

###################################################################################################

cdef class Histogrammer:
    """
    The Histogrammer class keeps an unbounded histogram of a continuous variable,
    and keep some information on where the extreme outliers have occurred
    """
    def __init__(self, granularity = 20, max = 1, bins = None, counts = None):
        """
        Constructor for histogram.

        granularity: histogram bin width; <0 for a standard histogram (width 1)
        max: initial maximum for histogram
        """
        if bins:
            self.bins = IntArray(len(bins), 0)
            self.counts = DoubleArray(len(counts), 0)

            for i in range(len(bins)):
                self.bins.array[i] = bins[i]

            for i in range(len(counts)):
                self.counts.array[i] = counts[i]
        else:
            self.bins = IntArray(1, 1)
            self.counts = DoubleArray(1, 0)

        self.granularity = granularity   # 20 = about a deciban
        self._histbin(max)

    cdef int _histbin(self, int value):
        """
        Convert value into histogram bin
        """
        cdef int lastBinValue = self.bins.getLast()
        cdef int b = 0

        # Simple histogram. Bins are integer valued
        if self.granularity <= 0:

            while lastBinValue <= value:
                lastBinValue += 1
                self.bins.append(lastBinValue)
                self.counts.append(0)

            return value

        # Bins are still integer-valued here.
        else:
            b = bisect(self.bins, value)

            if b < self.bins.getSize():
                return b

            while lastBinValue <= value:
                lastBinValue = 1 + lastBinValue + (lastBinValue -1) // self.granularity
                self.bins.append(lastBinValue)
                self.counts.append(0)

            return bisect(self.bins, value)

    cdef void enter(self, int value, double number=1):
        """
        Enter a number of values into histogram, in the
        relevant bin.
        """
        cdef int curbin = self._histbin(value)
        self.counts.addValue(curbin, number)

    @cython.profile(False)
    cdef void fast_enter(self, int value, int number):
        """
        A fast version of the '_histbin' function from histogram.py.
        Use only when granularity <= 0
        """
        cdef int i = self.bins.getLast()

        while i <= value:
            i += 1
            self.bins.append(i)
            self.counts.append(0)

        self.counts.addValue(value, number)

    def __getitem__(self, int idx):
        """
        Emulate a smooth histogram, by linearly interpolating
        between bins.
        """
        cdef int b = bisect(self.bins, idx)
        cdef int left = 0
        cdef int right = 0

        if b >= len(self.bins):
            return 0

        if b == 0:
            left = 0
        else:
            left = self.bins[b-1]

        right = self.bins[b]
        return (self.counts[b] + right - 1 - idx) / (right-left)

    def __getslice__(self, int i, int j):
        """
        This allows the user to use python's slice syntax on instances of
        the histogram class.
        """
        return [ self.__getitem__(k) for k in range(i,j) ]

    def __len__(self):
        """
        This allows the user to call len(x) where x is an instance of
        the histogram class. Length is always equal to the value of the
        last bin.
        """
        return self.bins.getLast()

    def total(self):
        """
        Return the sum of all values in this histogram
        """
        cdef double theSum = 0.0

        for i from 0 <= i < self.counts.getSize():
            theSum += self.counts.array[i]

        return theSum

    def mean(self):
        """
        Compute and return the mean value of this
        histogram.
        """
        cdef double theSum = 0.0
        cdef int i = 0
        cdef int mid = 0

        for i from 0 <= i < self.bins.getSize():

            if i == 0:
                mid = (0 + self.bins[i]) // 2
            else:
                mid = (self.bins[i-1] + self.bins[i]) // 2

            theSum += mid * self.counts[i]

        cdef double n = self.total()

        if n > 0:
            return theSum/n
        else:
            return 0

    def percentile(self, double pct):
        """
        Compute and return the specified percentile value of this
        histogram.
        """
        cdef int i = 0
        cdef double n = self.total()
        cdef double m = 0

        if n == 0:
            return 0

        for i from 0 <= i < self.bins.getSize():
            m += self.counts[i]

            if m >= n*pct:
                break

        if i == 0:
            return (0 + self.bins[i])//2

        return (self.bins[i-1] + self.bins[i]) // 2

    def median(self):
        """
        Compute and return the median value of this
        histogram.
        """
        return self.percentile(0.5)

    def output(self, stream, label, max, lane, mate, selection=None):
        """
        Output function for histogram. Writes to the specified
        stream.
        """
        bins = IntArray(1,0)

        for i in range(1, len(self.bins)):
            bins.append((self.bins[i-1] + self.bins[i]) // 2)

        maxidx = bisect(self.bins, max)
        bns, cnts = bins[:maxidx], [ c for c in self.counts[:maxidx] ]

        if selection is not None:
            i = 0
            while i < len(bns):
                if bns[i] in selection:
                    i += 1
                else:
                    del bns[i], cnts[i]
        else:
            # remove runs of 0s
            zstart, runs = -1, []

            for i in range(len(bns)):
                if cnts[i] == 0:
                    if zstart == -1: zstart = i
                else:
                    if zstart > -1: runs.append( (zstart, i) )
                    zstart = -1

            if zstart > -1:
                runs.append( (zstart, len(bns)) )

            runs.reverse()

            for zstart, zend in runs:
                if zend-zstart > 2:
                    del bns[zstart+1:zend-1]
                    del cnts[zstart+1:zend-1]

        histogram_output(stream, label, bns, cnts, "_bins", "_counts", lane, mate)
        return bns

###################################################################################################
