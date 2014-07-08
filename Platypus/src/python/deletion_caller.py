"""
"""
from __future__ import division

import pysam
import sys
import numpy as np
import matplotlib

matplotlib.interactive(False)
matplotlib.use('PDF')

from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages

from itertools import izip
from bisect import bisect_left

infile = pysam.Samfile(sys.argv[1], 'rb')
region = sys.argv[2]
start = int(sys.argv[3])
stop = int(sys.argv[4])

print infile.header.keys()
print "Region length is ",infile.lengths[infile.references.index(region)]
fetched = infile.fetch(region, start, stop)
print "Extracting ",region,start,stop,"..."

readPos = []
readMQ = []
insertSizes = []

for read in fetched:

    if read.mapq < 20:
        continue

    if read.is_paired:
        if read.is_proper_pair:

            if read.pos < read.mpos:
                insertSizes.append(read.isize)
                readPos.append(read.pos)
                readMQ.append(10**(-0.1*(read.mapq)))
                #pread[pi].append([read.pos,read.pos+read.rlen,read.mpos,read.mpos+read.rlen])
                #rd2.append([read.pos,read.pos+read.rlen,read.mpos,read.mpos+read.rlen,pi,abs(read.isize),255,read.mapq])
            else:
                pass
                #pread[pi].append([read.mpos,read.mpos+read.rlen,read.pos,read.pos+read.rlen])
                #rd2.append([read.mpos,read.mpos+read.rlen,read.pos,read.pos+read.rlen,pi,abs(read.isize),read.mapq,255])

iSizeMedian = np.median(insertSizes)
iSizeSD = np.std(insertSizes)

print "Insert size median = %s. SD = %s." %(iSizeMedian, iSizeSD)

upperBound = iSizeMedian + 2*iSizeSD
lowerBound = iSizeMedian - 2*iSizeSD

#for i in range(len(insertSizes)):
#    insertSizes[i] = max(0,insertSizes[i] - iSizeMedian)


windowInsertSizes = []
windowMidPoints = []

windowSlide = 1
windowSize = 2*iSizeMedian
windowStart = readPos[0]
windowEnd = windowStart + windowSize

lastReadPos = readPos[-1]
index = 0

while windowStart < lastReadPos:

    index += 1
    windowStart += windowSlide
    windowEnd += windowSlide
    windowStartIndex = bisect_left(readPos, windowStart)
    windowEndIndex = bisect_left(readPos, windowEnd)
    windowInsertSizeTot = 0

    if windowStartIndex == windowEndIndex:
        continue

    if windowStartIndex >= len(readPos):
        break

    if index % 10000 == 0:
        print "Index = %s. Window start = %s. Window End = %s" %(index, windowStart, windowEnd)

    nReads = 0

    for i in range(windowStartIndex, windowEndIndex):
        nReads += 1
        pos = readPos[i]
        iSize = insertSizes[i]
        mapRightProb = 1.0 - readMQ[i]
        #windowInsertSizeTot += (mapRightProb*iSize)
        windowInsertSizeTot += (iSize)

    windowInsertSizeMean = windowInsertSizeTot / nReads
    windowInsertSizes.append(max(0,windowInsertSizeMean - iSizeMedian))
    windowMidPoints.append((windowStart + windowEnd) / 2.0)


pdfFile = PdfPages("InsertSizes.pdf")
fig = pyplot.figure(figsize=(12,12))
fig.suptitle("Insert Sizes")
ax = fig.add_subplot(111)
#ax.plot(readPos, insertSizes, label='Insert Sizes')
ax.plot(windowMidPoints, windowInsertSizes, label='Window Insert Sizes')
ax.legend()
pdfFile.savefig(fig)
pyplot.clf()
pyplot.cla()
pdfFile.close()



#for pos,iSize in izip(readPos,insertSizes):
#    if iSize > upperBound:
#        print pos,iSize

infile.close()
