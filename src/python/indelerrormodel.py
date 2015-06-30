import math
import heapq
import random
from collections import defaultdict

import pysam
import cerrormodel
import chaplotype
import platypusutils

#
# TODO:
#
# Calculate error models per read group, rather than entire BAM
#
# Do left alignment of indels (since Stampy doesn't do this!)
#
# Filter for read pos existing, is_qc_fail
#
# Repackage so that "utils" module is not stand alone
#
# Write BAM not SAM
#
# Refactor PreprocessBam so that it's not multiple pages long, and
#  to separate the looping, the error model code, and the broken-read code
#

# SHELVED:
#
# Compute likelihood of read under some useful default model (or re-align after model has been computed); put model in BAM
#  (Needs bit of rewiring of Platypus code - do when we get to assembly)
#
# Accumulate coverage along read, for per-cycle errors
#  (Not for now; per-cycle errors not a priority)
#
# Still a coverage accounting error (assertion fails sometimes; too many not quite full hom alt counts) 
#  (Does not affect results)
#


#########################################################################################
#
# Some math functions
#
#########################################################################################

def stirling(n):
    # http://en.wikipedia.org/wiki/Stirling%27s_approximation
    if n<=20: return math.factorial(n)
    return math.sqrt(2*math.pi*n)*(n/math.e)**n

def logstirling(n):
    if n<=20: return math.log(math.factorial(n))
    return 0.5 * math.log(2*math.pi*n) + n * (math.log(n) - 1)

def npr(n,r):
    return stirling(n)/stirling(n-r)

def ncr(n,r):    
    return stirling(n)/(stirling(r)*stirling(n-r))

def sample_binomial(p, n):
    """ trivial implementation """
    assert 0 <= p <= 1.0
    k = 0
    for i in range(n):
        if random.random() < p:
            k += 1
    return k


#########################################################################################
#
# Definitions
#
#########################################################################################

CIGAR_M = 0
CIGAR_I = 1
CIGAR_D = 2
CIGAR_N = 3
CIGAR_S = 4
CIGAR_MA = 7
CIGAR_X = 8

MAX_POS = 3000000000

MAX_UNIT_LENGTH = 12    # from tandem.c.  Maximum repeat unit length that is considered
MIN_TANDEM_LENGTH = 4   # minimum length of repeat tract 

MAX_INDELS_PER_READ = 1 # reads with more indels are considered mapping errors
MIN_INDEL_DISTANCE = 10 # nearby indels are lumped together

MAX_READ_LENGTH = 200

USEBFGS = False         # my own simple stepwise adaptive optimizer works better

BASEQ = 33              # base of the Q score encoding

DEFAULT_INDEL_ERROR_MODEL = chaplotype.default_indel_error_model

EXTRAPOLATION_EXPONENTS = [2.8, 2.0, 1.5, 1.4, 1.4]    # per-base increase for unit length 1,2,...; used to extrapolate error models
EXTRAPOLATION_RATE_THRESHOLD = 0.05                    # use extrapolation when error rates exceed this
EXTRAPOLATION_HET_THRESHOLD = 0.15                     # use extrapolation when heterozygosity exceed this

if USEBFGS:
    import scipy.optimize

#########################################################################################
#
# Main loop
#
#########################################################################################

def PreprocessBam( bamFileName, fastaFileName, orphanBamFileName=None,
                   minCoverage=6, maxCoverage=150, minErrors=5, maxNumReads=1e10,
                   minMapQ=20, minReadQ=10, minAnchor=4, maxReadStretch=200, min_isize=500, min_tandem_length=MIN_TANDEM_LENGTH ):
    """ Computes indel error model from BAM file, and extracts orphan reads.
        The min/maxCoverage variables should be set to reasonable limits to exclude iffy regions.
        maxReadStretch is the maximum aligned extent of a read. """

    MAX_REPEAT_LENGTH = 64    # fixed maximum due to ffsll instruction (coreutils/tandem.c)

    bamfile = pysam.Samfile( bamFileName, mode='rb' )

    # extract list of read group identifiers
    readgroups = bamfile.header.get('RG',None)
    if type(readgroups) == dict: readgroups = [readgroups]
    if readgroups: 
        readgroupdict = {}
        for idx, ident in enumerate( rg.get('ID',"") for rg in readgroups ):
            readgroupdict[ident] = idx

    # open orphan bam - write SAM file for now
    if orphanBamFileName:
        #orphanBamFile = htslibWrapper.Samfile( orphanBamFileName, mode='wb', template=bamfile )
        orphanSamFile = open( orphanBamFileName, 'w' )
        orphanSamFile.write( bamfile.text )
    else:
        #orphanBamFile = None
        orphanSamFile = None

    # open fasta file
    fastafile = pysam.Fastafile( fastaFileName )

    # make a coverage ring buffer
    covbuflen = max( maxReadStretch*2, 1000 )
    coveragebuf = CoverageBuffer( covbuflen, minReadQ )

    # make read buffer, to implement limited lookback
    readbuffer = Readbuffer( bamfile, size = maxReadStretch + MAX_UNIT_LENGTH + MAX_REPEAT_LENGTH )

    # make indel queue
    indelqueue = Indelqueue( minAnchor, minReadQ )

    # make repetitive region queue
    repeatqueue = Repeatqueue( fastafile, covbuflen, maxReadStretch + MAX_UNIT_LENGTH + MAX_REPEAT_LENGTH,
                               min_tandem_length = min_tandem_length )

    # make object to collect results
    indelhistogram = IndelHistogram( minCoverage, maxCoverage, minErrors )
    
    # setup main loop
    curchrom = -1      # int!
    curchromS = ""     # string
    curpos = -1
    repeatqueue.reset( curchromS, curpos )

    num_reads = 0
    for read in readbuffer:

        num_reads += 1

        if num_reads % 100000 == 0:
            print "Processed %s reads" % num_reads
            if num_reads > maxNumReads:
                break

        # extract index of read group
        if readgroups:
            readgroup = read.opt('RG')
            readgroupidx = readgroupdict[readgroup]
        else:
            readgroupidx = 0

        # decide whether to write as orphan read
        if ( orphanSamFile and
             read.is_paired and 
             #(not read.is_qc_fail) and
             (not read.is_duplicate) and
             (not read.is_secondary) and
             read.mapq >= minMapQ and
             (not read.is_proper_pair) and 
             (not read.mate_is_unmapped) and 
             (not read.is_unmapped) and 
             (read.rname != read.mrnm or abs(read.isize) > min_isize) ):
            
            isize = read.isize
            if isize == 0: isizeS = "*"
            else:          isizeS = str(-isize)

            mq = read.opt('MQ')
            uq = read.opt('UQ')              # TODO: this should be re-computed with a sensible indel error model

            optS = "\tMQ:i:%s" % read.mapq   # store original mapping quality
            optS += "\tUQ:i:%s" % uq         # store likelihood of aligning read to original location

            if readgroups: optS += "\tRG:Z:%s" % readgroup

            flag = read.flag

            # swap reverse strand bits
            if flag & 0x30 in [0x10, 0x20]:
                flag ^= 0x30

            # swap 1st/last bits
            if flag & 0xc0 in [0x40, 0x80]:
                flag ^= 0xc0

            orphanSamFile.write( "%s\t%i\t%s\t%i\t%i\t%iM\t%s\t%i\t%i\t%s\t%s%s\n" % (read.qname,
                                                                                      flag,
                                                                                      bamfile.getrname( read.mrnm ),   # mate
                                                                                      read.mpos,                       # mate
                                                                                      read.mapq,
                                                                                      read.rlen,
                                                                                      bamfile.getrname( read.tid ),    # this
                                                                                      read.pos,                        # this
                                                                                      isize,
                                                                                      read.seq,
                                                                                      read.qual,
                                                                                      optS) )

        # filter
        if (read.is_unmapped or 
            (read.is_paired and not read.is_proper_pair) or
            #read.is_qc_fail() or
            read.is_duplicate or
            read.is_secondary or
            read.mapq < minMapQ):
            continue

        # enter any indels in queue
        indelqueue.enter( read )

        # prepare to enter read into coverage buffer -- flush buffer
        readpos = read.pos
        readchrom = read.tid

        # calculate position to flush to
        if readchrom == curchrom:
            assert readpos >= curpos                      # BAM should be sorted
            flushpos = min(readpos, curpos +  covbuflen)
        else:
            if curchrom > -1:
                flushpos = curpos + covbuflen
            else:
                flushpos = curpos
                
        # get positions of next repetitive region, and next indel;
        # these are guaranteed to be on curchrom
        indelpos, indelallele, indelreadpos, indelread1 = indelqueue.nextindel()
        repeatqueue.setpos( flushpos )

        repeatend = curpos
        while repeatend <= curpos:
            repeatstart, repeatend, repeatunit, repeatlength = repeatqueue.nextsegment()

        # process region up to flushpos
        curcov = coveragebuf.getcurcov()
        while curpos < flushpos:

            # 1.  update coverage
            curcov += coveragebuf.getupdate( curpos )

            # 2.  process any indel
            alleles = []
            while indelpos <= curpos:
                
                # only process indels outside repetitive regions
                if curpos < repeatstart:
                    alleles.append( (curpos, indelallele, indelreadpos, indelread1) )

                indelpos, indelallele, indelreadpos, indelread1 = indelqueue.nextindel()

            if minCoverage <= curcov <= maxCoverage:
                indelhistogram.add( curcov, 1, 1, alleles )

            # 3.  process repetitive regions
            if curpos == repeatend:

                # skip bits of repeats that are overlapped by more error-prone repeats
                if repeatlength != -1:

                    # retrieve reads that fully overlap repetitive region, plus anchoring
                    reads = readbuffer.retrieve( repeatstart - minAnchor, repeatend + minAnchor, minReadQ )

                    repeatcov = len(reads)
                    if minCoverage <= repeatcov <= maxCoverage:

                        alleles = haplotypes( reads, repeatstart, repeatend+1, minAnchor, minReadQ )
                        indelhistogram.add( repeatcov, repeatunit, repeatlength, alleles )
                        indelhistogram.prune()

                repeatstart, repeatend, repeatunit, repeatlength = repeatqueue.nextsegment()

            curpos += 1

        # push back unprocessed indels and repeats
        repeatqueue.pushback( repeatstart, repeatend, repeatunit, repeatlength )
        indelqueue.pushback( indelpos, indelallele, indelreadpos, indelread1 )

        # process jumps and change of chromosome
        if readchrom != curchrom:
            curchrom = readchrom
            curchromS = bamfile.getrname( readchrom )
            curpos = readpos
            coveragebuf.reset( curchrom, curpos )
            repeatqueue.reset( curchromS, curpos )
        elif flushpos < readpos:
            curpos = readpos
            repeatqueue.setpos( curpos )
        
        # set current coverage
        coveragebuf.setcurcov( curcov )

        # enter read into coverage buffer, accounting for anchorage required
        coveragebuf.enter( read, minAnchor )

    # Done        
    print "Processed",num_reads,"reads"
    indelhistogram.computeModels()
    indelhistogram.output()


#########################################################################################
#
# Readbuffer, allows iteration as well as limited look-back
#
#########################################################################################


class Readbuffer():
    """ Allows sequential block access to reads """

    def __init__(self, bamfile, size, region=None):
        self.size = size
        self.generator = bamfile.fetch( until_eof=True, region=region )
        self.buffer = []
        self._bufend = 0
        self._current = None
        self._lastrname = -2  # -1 is used for unmapped reads
        self._resetbuf = False
        try:
            self.next()
        except StopIteration:
            pass

    def __iter__(self):
        return self

    def next(self):
        current = self._current
        if self._resetbuf:
            self._buffer = [self._firstread]
            self._bufend = 0            # first free location
            self._lastrname = self._firstread.tid
            self._resetbuf = False
        try:
            self._current = self.generator.next()
        except StopIteration:
            self._current = None
        if current:
            if current.tid != self._lastrname:
                self._resetbuf = True
                self._firstread = current  # first read in new buffer
            elif not current.is_unmapped:
                # store mapped reads
                if self._buffer[self._bufend].pos < current.pos - self.size:
                    # old read is stale; we can store read
                    self._buffer[self._bufend] = current
                    self._bufend = (self._bufend + 1) % len(self._buffer)
                else:
                    # old read is not stale; we need to extend buffer
                    # put the last-entered read at the end of the buffer, to allow
                    # append (repeatedly, if necessary)
                    if self._bufend != 0:
                        # TODO: this should be done by copying, rather than by building a new list
                        self._buffer = self._buffer[self._bufend:] + self._buffer[:self._bufend]
                        self._bufend = 0
                    self._buffer.append( current )
        else:
            raise StopIteration
        return current


    def retrieve(self, start, end, minreadq):
        """ Returns reads fully overlapping [start,end) """
        
        # start at end going backwards, for efficiency
        idx = self._bufend
        reads = []
        
        while True:
            idx = (idx - 1) % len(self._buffer)

            # for efficiency, do a quick test first
            if self._buffer[idx].pos <= start:

                ##DEBUG
                #rstart0, rend0 = pruned_ref_start_end(self._buffer[idx], minreadq)
                rstart, rend = platypusutils.pruned_ref_start_end(self._buffer[idx], minreadq, 0)
                #assert rstart0 == rstart
                #assert rend0 == rend

                if rend <= start:
                    # it is theoretically possible that subsequent reads, that start leftward of
                    # the current read, do still fully overlap [start,end), even though the current
                    # read is disjoint with [start,end).  In practice we will miss very few indel
                    # errors this way, as long error indels tend to occur in long tandems, and
                    # errors outside tandems are short.
                    break

                if rstart <= start and rend >= end:
                    reads.append( self._buffer[idx] )

            if idx == self._bufend:
                break

        return reads


#########################################################################################
#
# Indelqueue
#
#########################################################################################

        
class Indelqueue:

    def __init__(self, minAnchor, minReadQ ):
        self._heap = []
        self.minAnchor = minAnchor
        self.minReadQ = minReadQ

    def enter(self, read):
        indels = get_indels( read, self.minAnchor, self.minReadQ )
        for indel in indels:
            heapq.heappush( self._heap, indel )

    def nextindel(self):
        """ Returns tuple (genome_pos, length_and_type, read_pos, is_read_one) """
        if len(self._heap) == 0:
            return MAX_POS, 0, 0, True
        return heapq.heappop(self._heap)

    def pushback(self, indelpos, indelallele, readpos, isread1):
        heapq.heappush( self._heap, (indelpos, indelallele, readpos, isread1) )


#########################################################################################
#
# Repeat region queue
#
#########################################################################################


class Repeatqueue:

    def __init__(self, fastafile, lookahead, maxtandemsize, min_tandem_length):
        self.fastafile = fastafile
        self.lookahead = lookahead
        self.maxtandemsize = maxtandemsize
        self.queue = []
        self.min_tandem_length = min_tandem_length
        self.reset( None, -MAX_POS )
        self._warn_chromosomes = {}

    def reset(self, chromS, pos):
        self.chrom = chromS
        self.curpos = -MAX_POS
        self.setpos( pos )

    def setpos(self, pos):
        if pos < self.curpos: return
        if not self.chrom: return
        # read from self.curpos to pos+lookahead
        # but do not process more than 2*lookahead
        topos = pos + self.lookahead
        frompos = max( self.curpos, pos, 0 )
        # obtain sequence
        sequence = self.read( self.chrom, frompos - self.maxtandemsize, topos + self.maxtandemsize )
        # get repeats as (pos, size, unit) elements
        repeats = cerrormodel.get_repeats( sequence, self.min_tandem_length, frompos - self.maxtandemsize )
        # enter into queue
        newreps = []
        for reppos, size, unit in repeats:
            unit = unit.upper()   # Ignore strandedness
            if self.curpos <= reppos < topos:
                # ensure repeats do not overlap; choose most error prone if they do,
                # and ignore rest of segment (as it will have a high indel rate, and
                #  may contaminate the background counts)
                newreps.append( (reppos, reppos+size, unit, size) )
                if len(newreps) >= 2:
                    if newreps[-2][1] > newreps[-1][0]:
                        if type(newreps[-2][2]) == type(""):
                            unit2 = len(newreps[-2][2])
                        else:
                            unit2 = newreps[-2][2]
                        if type(newreps[-1][2]) == type(""):
                            unit1 = len(newreps[-1][2])
                        else:
                            unit1 = newreps[-1][2]
                        if (cerrormodel.approx_indel_error_rate( newreps[-2][3], unit2) >
                            cerrormodel.approx_indel_error_rate( newreps[-1][3], unit1) ):
                            newreps[-1] = ( (newreps[-2][1], newreps[-1][1], -1, -1) )
                        else:
                            newreps[-2] = ( (newreps[-2][0], newreps[-1][0], -1, -1) )
        for repeat in newreps:
            self.pushback( *repeat )
        self.curpos = topos

    def read(self, chrom, start, end):
        if start<0:
            if end<=0: return "N" * (end-start)
            return ("N" * (-start)) + self.read( chrom, 0, end )
        seq = self.fastafile.fetch( chrom, start, end )
        if len(seq) == 0 and chrom not in self._warn_chromosomes:
            if len(self.fastafile.fetch( chrom, 0, 1 )) == 0:
                print "*** Warning: chromosome '%s' appears to be absent from reference file '%s'" % (chrom, self.fastafile.filename)
                self._warn_chromosomes[chrom] = 1
        if len(seq) < end-start:
            seq += "N"*(end-start-len(seq))
        return seq

    def nextsegment(self):
        next = None
        if len(self.queue):
            next = heapq.heappop( self.queue )
        if next:
            return next
        return MAX_POS, MAX_POS, 0, 0

    def pushback(self, repeatstart, repeatend, repeatunit, repeatlength):
        if repeatstart == MAX_POS:
            return
        heapq.heappush( self.queue, (repeatstart, repeatend, repeatunit, repeatlength) )


#########################################################################################
#
# Coverage buffer
#
#########################################################################################


class CoverageBuffer:

    def __init__(self, covbuflen, minReadQ ):
        self.covbuflen = covbuflen
        self.minReadQ = minReadQ
        self.reset( None, 0 )

    def reset( self, chromId, pos):
        self.covbuf = [0] * self.covbuflen
        self.curcov = 0

    def getcurcov(self):
        return self.curcov

    def setcurcov(self, curcov):
        self.curcov = curcov

    def getupdate(self, pos):
        update = self.covbuf[ pos % self.covbuflen ]
        self.covbuf[ pos % self.covbuflen ] = 0
        return update

    def enter(self, read, minAnchor):
        ##DEBUG
        #pos0, end0 = pruned_ref_start_end( read, self.minReadQ, minAnchor )
        pos, end = platypusutils.pruned_ref_start_end( read, self.minReadQ, minAnchor )
        #assert pos0 == pos
        #assert end0 == end
        if end - pos > self.covbuflen:
            # read too long
            return
        self.covbuf[pos % self.covbuflen] += 1
        self.covbuf[end % self.covbuflen] -= 1


#########################################################################################
#
# Indel histogram
#
#########################################################################################


class IndelHistogram:

    def __init__(self, minCoverage, maxCoverage, minErrCount, maxUnits=25000):
        self.minCoverage = minCoverage
        self.maxCoverage = maxCoverage
        self.minErrCount = minErrCount
        self.maxUnits = maxUnits
        self.numNoData = 0
        # Main data store: repeatunit -> length -> (length histogram, read -> unit count read pos histogram)
        self.histograms = {}

    
    def prune(self):
        """ Permanently remove repeat units with very little data """

        if len(self.histograms) - self.numNoData < self.maxUnits:
            return

        # remove units with fewer than minErrCount ticks
        for unit in self.histograms.keys():
            if not self.histograms[unit]: 
                continue
            ticks = sum( sum( sum(hist[0]) for hist in self.histograms[unit][replen].values() ) for replen in self.histograms[unit] )
            if ticks < self.minErrCount:
                self.histograms[unit] = None
                self.numNoData += 1


    def report(self):
        try:
            self.i += 1
        except Exception:
            self.i = 0
        if self.i % 10000 == 0:
            print "Histogram report:"
            print "Num units: ",len(self.histograms) - self.numNoData
            print "Num deleted:",self.numNoData
            print "Num unit/replens: ", sum( len(self.histograms[repeatunit]) 
                                             for repeatunit in self.histograms 
                                             if self.histograms[repeatunit] != None )
            print "Num unit/replen/coverage:", sum( sum( len( hist_unit[replen] ) for replen in hist_unit ) 
                                                    for hist_unit in self.histograms.values()
                                                    if hist_unit != None )
                

    def add( self, coverage, repeatunit, repeatlength, alleles ):
        """ Adds counts for the particular repeat unit, and for units this size """
        if type(repeatunit) == type(""):
            self._add( coverage, len(repeatunit), repeatlength, alleles )
        self._add( coverage, repeatunit, repeatlength, alleles )


    def _add( self, coverage, repeatunit, repeatlength, alleles ):
        # First select the histograms to add the observation to; and add new
        # histograms if necessary

        if not repeatunit in self.histograms:
            # Add a dictionary for the repeat length
            self.histograms[repeatunit] = {}
        hist_unit = self.histograms[repeatunit]
        if hist_unit == None:
            # Previously marked as having very little data
            return
        if not repeatlength in hist_unit:
            # Add a dictionary for the coverage
            hist_unit[repeatlength] = {}
        hist_unit_replen = hist_unit[repeatlength]
        if not coverage in hist_unit_replen:
            # Add a tuple (length histogram, read dictionary to unit counts hsitogram by cycle number)
            hist_unit_replen[coverage] = self._create_histogram( coverage )
        hurc = hist_unit_replen[coverage]

        # Now add the observation
        non_ref_allele_count, read_pos, is_read1 = self._nonrefcount( alleles, coverage )
        hurc[ 0 ][ non_ref_allele_count ] += 1

        #if non_ref_allele_count == 1:
        #    hurc[ 1 ][ is_read1 ][ read_pos ] += 1


    def _create_histogram(self, coverage):
        #return ( [0]*(coverage+1), {True: [0]*MAX_READ_LENGTH, False: [0]*MAX_READ_LENGTH} )
        return ( [0]*(coverage+1), None )


    def _rescale_data( self, newcounts, old_coverage, old_histogram, new_coverage ):

        """ Moves data from a higher coverage to a lowr coverage bin, for efficiency, and
            to improve estimates from very high coverage data """

        # ensure a histogram exists to receive the rescaled data
        if new_coverage not in newcounts:
            newcounts[new_coverage] = self._create_histogram(new_coverage)
        if new_coverage == old_coverage:
            # just copy
            for ac in range(old_coverage+1):
                newcounts[new_coverage][0][ac] += old_histogram[0][ac]
                #for read in [True,False]:
                #    for cycle in range(MAX_READ_LENGTH):
                #        newcounts[new_coverage][1][read][cycle] += old_histogram[1][read][cycle]
            return
        assert new_coverage < old_coverage
        # copy hom ref and hom alt (high counts)
        newcounts[new_coverage][0][0] += old_histogram[0][0]
        newcounts[new_coverage][0][new_coverage] += old_histogram[0][old_coverage]
        # move counts over
        for ac in range(1, old_coverage):
            if old_coverage > new_coverage*4:
                # the relative error (in ac) in the new histogram is more than 2* that of the old histogram,
                # so stochastically add other allele counts, and accept the additional error caused by 
                # adding the variances together
                for site in range(old_histogram[0][ac]):
                    newcounts[new_coverage][0][ sample_binomial( float(ac)/old_coverage, new_coverage ) ] += 1
            else:
                # the histograms are quite similar in size, so simply copy the counts to the
                # best fitting new position
                new_ac = ((ac-1)*(new_coverage-2))/(old_coverage-2) + 1
                newcounts[new_coverage][0][ new_ac ] += old_histogram[0][ac]
        ## move cycle-specific counts for the 1-bin over
        #for read in [True,False]:
        #    for cycle in range(MAX_READ_LENGTH):
        #        old_count = old_histogram[1][read][cycle]
        #        new_count = sample_binomial( float(new_coverage)/old_coverage, old_count )
        #        newcounts[new_coverage][1][read][cycle] += new_count


    def _nonrefcount(self, data, coverage ):
        # Data is a list of (pos, length_type, read_pos, is_read1) tuples.  
        # First turn into a dictionary of counts

        haps = defaultdict(int)
        annotations = {}
        coverage_remaining = coverage
        for pos, length_type, read_pos, is_read1 in data:
            haps[ length_type ] += 1
            coverage_remaining -= 1
            annotations[ length_type ] = (read_pos, is_read1)
            # see below
            if coverage_remaining == 0:
                break

        # This fails very infrequently, and is not very reproducible (it depends on the
        #  history rather than just local reads).  A very small number of errors will not
        #  influence final results, so there.
        #assert coverage_remaining >= 0

        if coverage_remaining > 0:
            haps[ 0 ] += coverage_remaining

        if len(haps)>2:
            # aggregate the minor alleles
            alleles = sorted( [(count,hap) for hap,count in haps.iteritems()] )
            minors = sum( count for (count,hap) in alleles[:-1] )
            haps = { alleles[-2][1]: minors, 
                     alleles[-1][1]: alleles[-1][0] }

        # if two alleles are present, but not the reference allele, map the major allele to the reference
        if len(haps)==2 and 0 not in haps:
            alleles = sorted( [(count,hap) for hap,count in haps.iteritems()] )
            minorallele = alleles[1][1]
            haps = { alleles[0][1]: alleles[0][0], 
                     0: alleles[1][0] }
            annotations[ 0 ] = annotations[ alleles[1][1] ]
                
        # calculate non-ref allele count
        count = sum( count for hap,count in haps.iteritems() if hap != 0 )

        if count != 1:
            return count, None, None

        nonrefallele = [ hap for hap,count in haps.iteritems() if hap != 0 ][0]
        return count, annotations[nonrefallele][0], annotations[nonrefallele][1]


    def computeModels(self):

        results = {}
        het_estimates = {}
        for repeat_unit in sorted( self.histograms.keys() ):
            hist_unit = self.histograms[repeat_unit]
            if not hist_unit:
                continue

            observations = sum( sum( sum(hist_unit[k][c][0]) for c in hist_unit[k]) for k in hist_unit )
            nonrefs = sum( sum( sum(hist_unit[k][c][0][1:]) for c in hist_unit[k]) for k in hist_unit )
            if observations < 10*self.minErrCount or nonrefs < 2*self.minErrCount:
                continue

            error_rate_guess = 1e-3
            for repeat_length in sorted( hist_unit.keys() ):
                hist_unit_replen = hist_unit[repeat_length]
                add_neighbours = False

                while True:
                    # rescale data for high coverage, where the expected numbers of errors is >0.5
                    # also sparsen the data for higher coverage, to reduce computation time
                    newcounts = {}
                    for cov in sorted(hist_unit_replen.keys()):
                        if cov * error_rate_guess < 0.5:
                            density = 1 + (cov // 15)
                            rescaled_cov = cov - (cov % density)  # may be ==cov, in which case _rescale_data copies
                        else:
                            rescaled_cov = max(6,int(0.5 / error_rate_guess))
                        self._rescale_data( newcounts, cov, hist_unit_replen[cov], rescaled_cov )
                        if add_neighbours:
                            neighbour_hurc = hist_unit.get( repeat_length-1, {} ).get( cov, None )
                            if neighbour_hurc != None: 
                                self._rescale_data( newcounts, cov, neighbour_hurc, rescaled_cov )
                            neighbour_hurc = hist_unit.get( repeat_length+1, {} ).get( cov, None )
                            if neighbour_hurc != None: 
                                self._rescale_data( newcounts, cov, neighbour_hurc, rescaled_cov )

                    err_counts = sum( newcounts[cov][0][1] for cov in newcounts.keys() )

                    if err_counts < self.minErrCount:
                        # low counts -- add counts of neighbouring repeat lengths and try again;
                        # if that doesn't solve the issue, don't fit a model
                        if not add_neighbours:
                            add_neighbours = True
                        else:
                            break

                    # fit model to counts
                    N00, N01, N11, epsilon, beta = fitmodel( newcounts, errorguess = error_rate_guess, usebfgs=USEBFGS )

                    # if the estimated error rate is substantially larger than the guessed rate,
                    # redo the rescaling as this may further increase the error rate
                    if epsilon < error_rate_guess * 1.5:
                        break
                    
                    # set guess for next iteration for this data
                    error_rate_guess = epsilon
                
                    
                if err_counts < self.minErrCount:
                    # too little data.
                    continue

                # set guess for next iteration for NEXT data
                error_rate_guess = epsilon

                # store results
                if repeat_unit not in results:
                    results[repeat_unit] = {}

                # store model as phred score -- convert to characters later
                results[repeat_unit][repeat_length] = -10.0*math.log(epsilon) / math.log(10.0)

                # store het estimate
                het_estimates[ (repeat_unit, repeat_length) ] = N01/(N00+N01+N11+1e-10)

                ##DEBUG
                counts = sum( sum( hist_unit_replen[cov][0] ) for cov in hist_unit_replen.keys() )
                output = (repeat_unit, repeat_length, "%s\t%s\t%s\t%s\t%s\t%1.6f\t%1.6f\t%1.6f" % (repeat_unit,
                                                                                                   repeat_length,
                                                                                                   observations,
                                                                                                   counts,
                                                                                                   err_counts,
                                                                                                   het_estimates[ (repeat_unit, repeat_length) ],
                                                                                                   beta,
                                                                                                   epsilon ) )
                print output[2]

        # build models
        models = {}

        try:
            base = results[1][1]  # if this fails, there was no data
        except Exception:
            raise ValueError("Too few observations to build error model!")

        for repeat_unit in results:
            if type(repeat_unit) == type(""):
                unitlen = len(repeat_unit)
            else:
                unitlen = repeat_unit
            extrapolate_exp = 10 * math.log( EXTRAPOLATION_EXPONENTS[ min(unitlen, len(EXTRAPOLATION_EXPONENTS)-1) ] ) / math.log(0.1)
            model = [base] # dummy zero-th entry is removed afterwards
            for repeat_length in sorted(results[repeat_unit]):
                while repeat_length > len(model):
                    # extrapolate in the absence of any data
                    if len(model) < 5 + unitlen:
                        model.append( model[-1] )
                    else:
                        ##DEBUG
                        print " -- (repunit %s replen %s: extrapolating without data)" % (repeat_unit, len(model))
                        model.append( model[-1] + extrapolate_exp )
                # we have data; but extrapolate when the error rate or heterozygosity become too high
                if (het_estimates[ (repeat_unit, repeat_length) ] > EXTRAPOLATION_HET_THRESHOLD or
                    results[repeat_unit][repeat_length] < -10.0*math.log(EXTRAPOLATION_RATE_THRESHOLD) / math.log(10.0)):
                    print " -- repunit %s replen %s: prepare to extrapolating with data" % (repeat_unit, repeat_length)
                    extrapolated_rate = model[-1] + extrapolate_exp
                else:
                    extrapolated_rate = model[-1]
                model.append( min( extrapolated_rate, results[repeat_unit][repeat_length] ) )
                ##DEBUG
                if model[-1] == extrapolated_rate:
                    print " -- repunit %s replen %s: extrapolating with data" % (repeat_unit, repeat_length)
            while model[-1] > 0:
                model.append( model[-1] + extrapolate_exp )
            # convert to a string
            models[ repeat_unit ] = ''.join( chr(BASEQ + max(0,int(phred + 0.5))) for phred in model[1:] )

        # prune superfluous models
        for repeat_unit in models.keys():
            if type(repeat_unit) == type(""):
                unit_len = len(repeat_unit)
                generic_model = models[ unit_len ]
                for idx, c in enumerate( generic_model ):
                    if idx < len(models[repeat_unit]) and ord(models[repeat_unit][idx]) < ord(c)-1:
                        break
                else:
                    # specific model isn't predicting substantially larger indel rates -- superfluous
                    print " Removing model for ",repeat_unit
                    del models[repeat_unit]

        self.models = models


    def output(self):
        print self.models


#########################################################################################
#
# Indel parse utilities
#
#########################################################################################


def get_indels( read, minAnchor, minReadQ ):
    """ Returns (position, length_type, read position, is_read_1) tuples """
    
    indels = []
    cigar = read.cigar
    if not cigar or len(cigar) == 1: return indels

    # the algorithm to determine if an indel is properly anchored is a
    # mix of two flawed algorithms: one determines only if a read contains
    # sufficiently high-quality bases beyond and before the indel; and one
    # determines only if a read is aligned sufficiently long beyond and before
    # the indel.  The correct algorithm would check both simultaneously.
    # however in practice there will be very few corner cases, 
    # and I can't currently be bothered
    
    #readstart, readend = pruned_read_start_end(read, minReadQ, minAnchor)
    readstart, readend = platypusutils.pruned_read_start_end(read, minReadQ, minAnchor)
    putative_indels = []
    is_anchored = False
    nucpos = read.pos
    if read.is_reverse:
        # we start at rlen, not rlen-1; this ensures deletions
        # get placed correctly, and for insertions we need to 
        # adjust by the insertion length
        readpos = read.rlen
        readinc = -1
        # convert to sequenced read coordinates
        readstart, readend = readpos - readend, readpos - readstart
    else:
        readpos = 0
        readinc = 1
    for op,le in cigar:
        if op == CIGAR_M:
            if le >= minAnchor:
                if is_anchored:
                    # anchored left and right -- accept indels
                    for indel in putative_indels:
                        indels.append( indel )
                    putative_indels = []
                is_anchored = True
            nucpos += le
            readpos += readinc*le
        elif op == CIGAR_D:
            if len(putative_indels)==0 or (putative_indels[-1][0]-nucpos > MIN_INDEL_DISTANCE and 
                                           putative_indels[-1][2]-readpos > MIN_INDEL_DISTANCE):
                # ensure indel is anchored by high-quality bases
                if readstart <= readpos < readend:
                    putative_indels.append( (nucpos, -le, readpos, read.is_read1) )
            nucpos += le
        elif op == CIGAR_I:
            if readinc == 1: 
                inspos = readpos
            else:           
                inspos = readpos - le
            if len(putative_indels)==0 or (putative_indels[-1][0]-nucpos > MIN_INDEL_DISTANCE and 
                                           putative_indels[-1][2]-readpos > MIN_INDEL_DISTANCE):
                # ensure indel is anchored by high-quality bases
                if readstart <= inspos and inspos + le < readend:
                    # further ensure that the insertion is not at the very right end of the read in
                    # the forward direction.  For insertions there are otherwise 1 more allowed positions
                    # than there are bases in the read, and this gives problems with the coverage.
                    # for read already in the forward direction, this condition is implied in the 
                    # condition inspos + le < readend; but for reads mapped in the reverse direction
                    # it has to be explicitly enforced:
                    if readinc == 1 or readstart < inspos:
                        putative_indels.append( (nucpos, le, inspos, read.is_read1) )
            readpos += readinc*le
        elif op == CIGAR_S:
            readpos += readinc*le
        elif op == CIGAR_N:
            nucpos += le
    if len(indels) > MAX_INDELS_PER_READ:
        return []
    return indels


def pruned_read_start_end(read, minq, minAnchor=0):
    """ Calculates the start and end of the read, pruning low-quality bases at either end """
    qual = read.qual
    readlen = len(qual)
    # find pruned start and end on read
    readstart = 0
    readend = readlen
    while readend > 0 and ord(qual[readend-1]) - BASEQ < minq:
        readend -= 1
    readend -= minAnchor
    while readstart < readend and ord(qual[readstart]) - BASEQ < minq:
        readstart += 1
    readstart += minAnchor
    return readstart, readend


def pruned_ref_start_end(read, minq, minAnchor=0):
    """ Calculates the reference start and end of the read, pruning low-quality bases at either end,
        and removing an anchor sequence at either end """
    cigar = read.cigar
    refpos = read.pos
    readstart, readend = pruned_read_start_end(read, minq, minAnchor)
    # deal with very low quality reads
    if readstart >= readend or not cigar: 
        return refpos, refpos
    # now convert start/end to refstart/refend
    readpos = 0
    refstart, refend = -1, -1
    for op, le in cigar:
        if op == CIGAR_M or op == CIGAR_I or op == CIGAR_S:
            # does start lie in this segment?
            if readpos <= readstart < readpos+le:
                # compute ref position and store
                if op == CIGAR_M:
                    refstart = refpos + (readstart - readpos)
                else:
                    refstart = refpos
            # does end lie in this segment?
            if readpos < readend <= readpos+le:
                # compute ref position and store
                if op == CIGAR_M:
                    refend = refpos + (readend - readpos)
                else:
                    refend = refpos
            readpos += le
        if op == CIGAR_M or op == CIGAR_D or op == CIGAR_N:
            refpos += le
            # deal with the case of a deletion at the right boundary
            if readpos == readend:
                refend = refpos
    # problem
    assert refend != -1 and refstart != -1
    return refstart, refend

                
def haplotypes( reads, repeatstart, repeatend, minAnchor, minReadQ ):
    haps = []
    for read in reads:
        indels = get_indels( read, minAnchor, minReadQ )
        for indel in indels:
            pos = indel[0]
            if repeatstart <= pos < repeatend:
                haps.append( indel )
    return haps





#########################################################################################
#
# Actual indel error model
#
#########################################################################################


def model( pars, counts ):
    """ computes rate, lambda of the possible outcomes --- to score a Poisson approximation. """
    
    cov = len(counts)-1    # for coverage C, counts array includes counts for 0,..,C
    lambdas_noerr = [0] * (cov+1)
    lambdas = [0] * (cov+1)
    
    N00, N01K, N11K, epsilonK, beta = pars
    N00 = max(0, N00)
    N01 = max(0, N01K / 1000.0)
    N11 = max(0, N11K / 1000.0)
    epsilon = max(0, epsilonK / 1000.0)

    ## rescale so N00+N01+N11 = sum(counts)
    # N = N00 + N01 + N11
    ## No, instead weigh with counts, but do not scale: reduces degrees of freedom in maximization
    N = 1
    cov_scaling = sum(counts)/N
    N00 *= cov_scaling
    N01 *= cov_scaling
    N11 *= cov_scaling

    # hom ref model
    lambdas_noerr[0] = N00
    
    # het model
    for k in range(0,cov+1):
        lambdas_noerr[k] += N01 * ncr(cov,k) * math.pow(beta,k) * math.pow(1-beta,cov-k)

    # hom alt model
    lambdas_noerr[cov] += N11

    # error model
    for k in range(0, cov+1):
        # 0 errors
        error_1_prob = cov*epsilon*math.pow(1-epsilon,cov-1)
        error_2_prob = cov*(cov-1)*0.5*epsilon*epsilon*math.pow(1-epsilon,cov-2)
        error_0_prob = 1.0 - error_1_prob - error_2_prob

        lambdas[k] += error_0_prob * lambdas_noerr[k]

        # 1 error
        if k>0:
            lambdas[k-1] += k*error_1_prob*lambdas_noerr[k]/cov
        if k<cov:
            lambdas[k+1] += (cov-k)*error_1_prob*lambdas_noerr[k]/cov

        # 2 errors
        if k>1:
            lambdas[k-2] += k*(k-1)*error_2_prob*lambdas_noerr[k]/(cov*(cov-1))
        if k<cov-1:
            lambdas[k+2] += (cov-k)*(cov-1-k)*error_2_prob*lambdas_noerr[k]/(cov*(cov-1))
        # one error on a ref background, one on an alt background
        lambdas[k] += 2*k*(cov-k)*error_2_prob*lambdas_noerr[k]/(cov*(cov-1))   

    # likelihood
    ll = 0
    if cov < 10:
        rnge = range(0,cov+1)
    else:
        rnge = range(0,cov-2) + [cov]
        # do not include cov-1 and cov-2; these are often populated by mismapped reads in a high-coverage
        # region, rather than true errors mutating alt to ref
    for k in rnge:
        lmda = lambdas[k]
        ll += counts[k]*math.log(lmda + 1e-10) - lmda - logstirling(counts[k])

    return ll


def multimodel( pars, countdict ):
    """ Returns minus the log likelihood, for minimization """
    return -sum( model( pars, countdict[cov][0] ) for cov in countdict )


def fitmodel( countdict, errorguess = 0.001, usebfgs=False ):
    """ fits a (homref + error) + het model to count data """
    
    N00, N01, N11 = 1e-10,1e-10,1e-10
    for cov in countdict:
        N00 += float( countdict[cov][0][0] )
        N01 += float( sum(countdict[cov][0][1:-1]) )
        N11 += float( countdict[cov][0][-1] )
    N = N00+N01+N11
    N00 /= N
    N01 /= N
    N11 /= N
    eps = errorguess
    beta = 0.5

    # scale N01, N11 and eps by 1000
    pars = [N00,N01*1000,N11*1000,eps*1000,beta]
    minpars = [1e-9,1e-9*1000,1e-9*1000,1e-6*1000,0.35]
    maxpars = [1e10,1e10*1000,1e10*1000,0.5*1000,0.65]

    if usebfgs:
        x, f, d = scipy.optimize.fmin_l_bfgs_b(multimodel, 
                                               pars, 
                                               fprime=None, 
                                               args=(countdict,), 
                                               approx_grad=True, 
                                               bounds=zip(minpars, maxpars),
                                               m=25, 
                                               factr=1e12,
                                               pgtol=1e-04, 
                                               epsilon=1e-05, 
                                               iprint=-1, 
                                               maxfun=15000)
        """
        x = scipy.optimize.fmin_cobyla(multimodel,
                                       pars,
                                       cons = [ lambda x:x[idx]-minpar for idx,minpar in enumerate(minpars) ] + \
                                              [ lambda x:maxpar-x[idx] for idx,maxpar in enumerate(maxpars) ],
                                       args=(countdict,),
                                       consargs = (),
                                       rhobeg = 0.05,
                                       rhoend = 1e-3,
                                       iprint=1,
                                       maxfun=1000)

        """
        pars = x

    else:

        dpars = [0.25,0.25,0.25,0.25,0.25]
        ddpars = 0.5
        ddparsplus = 1.2
        k = 0
        change = True
        ll = multimodel( pars, countdict )
        while sum(dpars)>0.01 and (k>0 or change):
            if k == 0: change = False
            parsplus = pars[:]
            parsminus = pars[:]
            parsplus[k] *= 1.0 + dpars[k]
            parsminus[k] /= 1.0 + dpars[k]
            if parsplus[k] < maxpars[k]:
                llplus = multimodel( parsplus, countdict )
            else:
                llplus = ll + abs(ll/2)
            if parsminus[k] > minpars[k]:
                llminus = multimodel( parsminus, countdict )
            else:
                llminus = ll + abs(ll/2)
            if ll <= min(llplus, llminus):
                dpars[k] *= ddpars
                change = True
            elif llplus < min(ll, llminus):
                pars[k] = parsplus[k]
                ll = llplus
                dpars[k] *= ddparsplus
                change = True
            else:
                pars[k] = parsminus[k]
                ll = llminus
                dpars[k] *= ddparsplus
                change = True
            k = (k+1) % len(pars)
    
    # undo scaling
    pars = [ pars[0], pars[1]/1000.0, pars[2]/1000.0, pars[3]/1000.0, pars[4] ]
    return pars
 




#########################################################################################
#
# Main
#
#########################################################################################



if __name__ == "__main__":

    import sys

    if len(sys.argv) != 4:
        print "Usage: %s bamfile fastafile orphanfile" % sys.argv[0]
        sys.exit(1)

    bamname = sys.argv[1]
    fastaname = sys.argv[2]
    orphanname = sys.argv[3]

    PreprocessBam( bamname, fastaname, orphanname, maxNumReads=1e11, minErrors=5 )
