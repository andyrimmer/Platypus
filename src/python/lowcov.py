import math
import heapq

import htslibWrapper
import pysam

#
# TODO:
#
# Current algorithm lumps reads from all read groups together
#


#########################################################################################
#
# Some math functions
#
#########################################################################################

def stirling(n):
    # http://en.wikipedia.org/wiki/Stirling%27s_approximation
    return math.sqrt(2*math.pi*n)*(n/math.e)**n

def logstirling(n):
    if n == 0: return 0
    return 0.5 * math.log(2*math.pi*n) + n * (math.log(n) - 1)

def npr(n,r):
    return (stirling(n)/stirling(n-r) if n>20 else
            math.factorial(n)/math.factorial(n-r))

def ncr(n,r):    
    return (stirling(n)/stirling(r)/stirling(n-r) if n>20 else
            math.factorial(n)/math.factorial(r)/math.factorial(n-r))


#########################################################################################
#
# Definitions and utilities
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

MAX_UNIT_LENGTH = 12    # from tandem.c
MIN_UNIT_LENGTH = 5 

#########################################################################################
#
# Main loop
#
#########################################################################################

def PreprocessBam( bamFileName, fastaFileName, orphanBamFileName=None,
                   minCoverage=6, maxCoverage=50, 
                   minMapQ=20, minAnchor=5, maxReadStretch=200, min_isize=500, min_tandemlength=5 ):
    """ Computes indel error model from BAM file, and extracts orphan reads.
        The min/maxCoverage variables should be set to reasonable limits to exclude iffy regions.
        maxReadStretch is the maximum aligned extent of a read. """

    MAX_REPEAT_LENGTH = 64    # fixed maximum due to ffsll instruction (coreutils/tandem.c)

    bamfile = htslibWrapper.Samfile( bamFileName, mode='rb' )

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
    coveragebuf = CoverageBuffer( covbuflen )

    # make read buffer, to implement limited lookback
    readbuffer = Readbuffer( bamfile, size = maxReadStretch + MAX_UNIT_LENGTH + MAX_REPEAT_LENGTH )

    # make indel queue
    indelqueue = Indelqueue()

    # make repetitive region queue
    repeatqueue = Repeatqueue( fastafile, covbuflen, maxReadStretch + MAX_UNIT_LENGTH + MAX_REPEAT_LENGTH,
                               min_tandemlength = min_tandemlength )

    # make object to collect results
    indelhistogram = IndelHistogram( minCoverage, maxCoverage )
    
    # setup main loop
    curchrom = -1      # int!
    chrchromS = ""     # string
    curpos = -1
    repeatqueue.reset( curchromS, curpos )

    for read in readbuffer:

        # extract index of read group
        if readgroups:
            readgroup = read.opt('RG')
            readgroupidx = readgroupdict[readgroup]
        else:
            readgroupidx = 0

        # decide whether to write as orphan read
        if ( orphanSamFile and
             read.is_paired() and 
             (not read.is_qc_fail()) and
             (not read.is_duplicate()) and
             read.mapq() >= minMapQ and
             (not read.is_proper_pair()) and 
             (not read.mate_is_unmapped()) and 
             (not read.is_unmapped()) and 
             (read.rname() != read.mrnm() or abs(read.isize()) > min_isize) ):
            
            isize = read.isize()
            if isize == 0: isizeS = "*"
            else:          isizeS = str(-isize)

            if readgroups: rgS = "\tRG:Z:" + readgroup
            else:          rgS = ""

            flag = read.flag()   # TODO: swap 1st/2nd, strand/matestrand, etc.

            orphanSamFile.write( "%s\t%i\t%s\t%i\t%i\t*\t%s\t%i\t%i\t%s\t%s\t%s%s\n" % (read.fastQName(),
                                                                                        flag,
                                                                                        bamfile.getrname( read.mrnm() ),
                                                                                        read.npos(),
                                                                                        read.mapq(),
                                                                                        bamfile.getrname( read.rname() ),
                                                                                        read.pos(),
                                                                                        isize,
                                                                                        read.seq(),
                                                                                        read.qual(),
                                                                                        rgS) )

        # filter
        if (read.is_unmapped() or 
            (read.is_paired() and not read.is_proper_pair()) or
            read.is_qc_fail() or
            read.is_duplicate() or
            read.mapq() < minMapQ):
            continue

        # enter any indels in queue
        indelqueue.enter( read )

        # prepare to enter read into coverage buffer -- flush buffer
        readpos = read.pos()
        readchrom = read.rname()  # int!

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
        indelpos, indelallele = indelqueue.nextindel()
        repeatqueue.setpos( flushpos )
        repeatstart, repeatend, repeatunit, repeatlength = repeatqueue.nextsegment()

        # process region up to flushpos
        curcov = coveragebuf.getcurcov()
        while curpos < flushpos:
            
            # 1.  update coverage
            curcov += coveragebuf.getupdate( curpos )

            # 2.  process any indel
            alleles = defaultdict(int)
            while indelpos <= curpos:
                
                # only process indels outside repetitive regions
                if curpos < repeatstart:
                    alleles[indelallele] += 1

                indelpos, indelallele = indelqueue.nextindel()

            if minCoverage <= curcov <= maxCoverage:
                indelhistogram.add( curcov, 1, 1, alleles )

            # 3.  process repetitive regions
            if curpos == repeatend:

                # retrieve reads that fully overlap repetitive region,
                # plus anchoring
                reads = readbuffer.retrieve( repeatstart - minAnchor, repeatend + minAnchor )

                repeatcov = len(reads)
                if minCoverage <= repeatcov <= maxCoverage:

                    alleles = haplotypes( reads, repeatstart, repeatend )
                    indelhistogram.add( repeatcov, repeatunit, repeatlength, alleles )

                repeatstart, repeatend, repeatunit, repeatlength = repeatqueue.nextsegment()

            curpos += 1

        # push back unprocessed indels and repeats
        repeatqueue.pushback( repeatstart, repeatend, repeatunit, repeatlength )
        indelqueue.pushback( indelpos, indelallele )

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

        # enter read into coverage buffer
        coveragebuf.enter( read )

    # Done        
    indelhistogram.computeModels()
    indelhistogram.output()


#########################################################################################
#
# Readbuffer, allows iteration as well as limited look-back
#
#########################################################################################


class ReadBuffer():
    """ Allows sequential block access to reads """

    def __init__(self, bamfile, size):
        self.size = size
        self.generator = bamfile.fetchAll()
        self.buffer = []
        self._bufend = 0
        self._current = None
        self._laststart = -1
        self._lastrname = -1
        self._resetbuf = False
        self._next()

    def _next(self):
        current = self._current
        if self._resetbuf:
            self.buffer = [self._firstread]
            self._bufend = 0            # first free location
            self._laststart = -1
            self._resetbuf = False
        try:
            self._current = self.generator.next()
        except StopIteration:
            self._current = None
        if current:
            if current.rname() != self._lastrname:
                self._resetbuf = True
                self._firstread = current  # first read in new buffer
            else:
                if self._buffer[self._bufend].pos() < current.pos() - self.size:
                    # old read is stale; we can store read
                    self._buffer[self._bufend] = current
                    self._bufend = (self._bufend + 1) % self.size
                else:
                    # old read is not stale; we need to extend buffer
                    # put the last-entered read at the end of the buffer, to allow
                    # append (repeatedly, if necessary)
                    if self._bufend != 0:
                        self._buffer = self._buffer[self._bufend:] + self._buffer[:self._bufend]
                        self._bufend = 0
                    self._buffer.append( current )
        return current


    def retrieve(self, start, end):
        """ Returns reads fully overlapping [start,end) """
        
        # start at end going backwards, for efficiency
        idx = self._bufend
        reads = []
        
        while True:
            idx = (idx - 1) % self.size

            if self._buffer[idx].pos() < start:
                break

            if self._buffer[idx].end() > end:
                idx = (idx - 1) % self.size
                continue

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

    def __init__(self, minAnchor):
        self._heap = []
        self.minAnchor = minAnchor

    def enter(self, read):
        n = read.getCigarLength()
        if n == 0: return

        putative_indels = []
        is_anchored = False
        nucpos = read.pos()
        for i in range(n):
            op = read.getCigarOpCode( i )
            le = read.getCigarOpLength( i )
            if op == CIGAR_M:
                if le >= minAnchor:
                    if is_anchored:
                        # anchored left and right -- accept indels
                        for indel in putative_indels:
                            heapq.heappush( self._heap, indel )
                        putative_indels = []
                    is_anchored = True
                nucpos += le
            elif op == CIGAR_D:
                putative_indels.append( (pos, -le) )
                nucpos += le
            elif op == CIGAR_I:
                putative_indels.append( (pos, le) )
            elif op == CIGAR_N:
                nucpos += le

    def nextindel(self):
        if len(self._heap) == 0:
            return MAX_POS, 0
        return heapq.heappop(self._heap)

    def pushback(self, indelpos, indelallele):
        heapq.heappush( self._heap, (indelpos, indelallele) )


#########################################################################################
#
# Repeat region queue
#
#########################################################################################


class Repeatqueue:

    def __init__(self, fastafile, lookahead, maxtandemsize, min_tandem_length = MIN_UNIT_LENGTH):
        self.fastafile = fastafile
        self.lookahead = lookahead
        self.maxtandemsize = maxtandemsize
        self.queue = []
        self.min_tandem_length = min_tandem_length
        self.reset( None, -MAX_POS )

    def reset(self, chromS, pos):
        self.chrom = chromS
        self.curpos = -MAX_POS
        self.setpos( pos )

    def setpos(self, pos):
        if pos < self.curpos: return
        if not self.chrom: return
        # read from self.curpos to pos+lookahead
        # but do not process more than 2*lookahead
        topos = pos + lookahead
        frompos = max( self.curpos, pos - lookahead, 0 )
        # obtain sequence
        sequence = self.read( self.chrom, frompos - maxtandemsize, topos + maxtandemsize )
        # get repeats as (pos, size, unit) elements
        repeats = cerrormodel.get_repeats( sequence, self.min_tandem_length, frompos - maxtandemsize )
        # enter into queue
        for reppos, size, unit in repeats:
            if self.curpos <= reppos < topos:
                self.pushback( reppos, reppos + size, unit, size )

    def read(self, chrom, start, end):
        if start<0:
            if end<=0: return "N" * (end-start)
            return ("N" * (-start)) + self.read( chrom, 0, end )
        seq = self.fastafile.fetch( chrom, start, end )
        # TODO: check length of reference, and add Ns myself if required
        assert len(seq) > 0
        if len(seq) < end-start:
            seq += "N"*(end-start-len(seq))
        return seq

    def nextsegment(self):
        return repeatstart, repeatend, repeatunit, repeatlength

    def pushback(self, repeatstart, repeatend, repeatunit, repeatlength):
        heapq.heappush( self.queue, (repeatstart, repeatend, repeatunit, repeatlength) )


#########################################################################################
#
# Repeat region queue
#
#########################################################################################


class CoverageBuffer:

    def __init__(self, covbuflen):
        self.covbuflen = covbuflen
        self.reset( None, 0 )

    def reset( self, chromId, pos):
        self.covbuf = [0] * self.covbuflen
        self.curpos = pos
        self.curcov = 0

    def getcurcov(self):
        return self.curcov

    def setcurcov(self, curcov):
        self.curcov = curcov

    def getupdate(self, pos):
        update = self.covbuf[ self.curpos % self.covbuflen ]
        self.covbuf[ self.curpos % self.covbuflen ] = 0
        return update

    def enter(self, read):
        pos = read.pos() % self.covbuflen
        end = read.end() % self.covbuflen
        self.covbuf[pos] += 1
        self.covbuf[end] -= 1


#########################################################################################
#
# Indel histogram
#
#########################################################################################


class IndelHistogram:

    def __init__(self, minCoverage, maxCoverage):
        self.minCoverage = minCoverage
        self.maxCoverage = maxCoverage
        self.minTandem = MIN_UNIT_LENGTH
        self.maxTandem = MAX_UNIT_LENGTH
        pass

    def add( self, coverage, repeatunit, repeatlength, alleles ):
        pass

    def computeModels(self):
        pass

    def output(self):
        pass









def haplotypes( reads, repeatstart, repeatend ):
    return alleles




#
# Assume that indels with low support, in otherwise well covered regions, are errors?
#

maxprocessedmotifs = 1e10
thinner = 4/15.0                    # fraction of reads considered
min_tot_count = 10                  # minimum count to report a motif

Read = collections.namedtuple('Read','flag, rname, pos, mapq, cigar, mname, mpos, isize, seq, qual')

CIGAR_M = 0
CIGAR_I = 1
CIGAR_D = 2
CIGAR_N = 3
CIGAR_S = 4
CIGAR_MA = 7
CIGAR_X = 8



def readLengthOnGenome( read ):
    """ returns length of read on genome """
    if read.cigar: return sum( [leng for op,leng in read.cigar if op in (CIGAR_M,CIGAR_D,CIGAR_N,CIGAR_MA,CIGAR_X)] )
    if read.seq: return len(read.seq)
    return 0



class ControlledPopen:

    def __init__(self, arguments):
        self.args = arguments

    def __enter__(self):
        self.pipe = subprocess.Popen( self.args, stdout = subprocess.PIPE)
        return self.pipe

    def __exit__(self, type, value, traceback):
        self.pipe.kill()
        del self.pipe


def readGenerator( bamfilename, chromosome ):
    """ Generates reads """

    arguments = ["samtools","view",bamfilename,chromosome]
    with ControlledPopen( arguments) as pipe:
        for line in pipe.stdout:
            if line.startswith('@'):
                continue
            if random.random() > thinner:
                continue
            label, flagS, rname, posS, mapqS, cigarS, mname, mposS, isizeS, seq, qual = line[:-1].split('\t')[:11]
            flag = int(flagS)
            mapq = int(mapqS)
            if posS == "*":
                pos = -1
            else:
                pos = int(posS)
            if mposS == "*":
                mpos = -1
            else:
                mpos = int(mposS)
            if cigarS == "*":
                cigar = None
            else:
                cigar = []
                p = -1
                for arg in re.split('[MIDNSHP=X]',cigarS)[:-1]:
                    p += len(arg) + 1
                    op = 'MIDNSHP=X'.find(cigarS[p])
                    cigar.append( (op, int(arg)) )
            yield Read(flag=flag,
                       rname=rname,
                       pos=int(posS),
                       mapq=mapq,
                       cigar=cigar,
                       mname=mname,
                       mpos=mpos,
                       isize=int(isizeS),
                       seq=seq,
                       qual=qual)


def filterReads( generator, minmapq=10, minanchor=5 ):
    """ Filters for properly mapped reads, and imposes minimum mapQ.
        Remove any read with >1 gap
        Remove reads not anchored with minimum anchor length """
    for read in generator:

        if (read.flag & 4 == 1) or (read.flag & 2 != 2) or read.pos == -1 or read.mapq < minmapq:
            continue

        gaps = 0
        matches = []
        for op,arg in read.cigar:
            if op==CIGAR_I or op==CIGAR_D:
                gaps += 1
            elif op==CIGAR_M:
                matches.append(arg)
        if gaps > 1:
            continue
        if len(matches)==0 or matches[0] < minanchor or matches[-1] < minanchor:
            continue

        yield read





def motifGenerator( infile, chromosome ):
    """ reads motif file """
    for line in infile:
        chrom, pos, hlen, hom, tlen, tandem = line[:-1].split('\t')
        assert chrom == chromosome
        yield (chrom, int(pos), int(tlen), tandem)


def haplotypes( reads, start, end ):
    """ Generates list of indel haplotypes (+:insertion, -:deletion, 0:reference) supported by the reads, and their multiplicity """
    haps = {}
    for read in reads:
        indel = 0
        pos = read.pos
        for op,arg in read.cigar:
            if op == CIGAR_M or op == CIGAR_N:
                pos += arg
            elif op == CIGAR_I:
                if start <= pos <= end+1:
                    indel = arg
            elif op == CIGAR_D:
                if start <= pos <= end+1:
                    indel = -arg
                pos += arg
        haps[indel] = haps.get(indel,0) + 1
    return haps




def estimateErrorRate( chromosome, motiffile, lowcovbam,
                       minmapq=10, minanchor=5, coverage=5 ):
    
    motifgenerator = motifGenerator( motiffile, chromosome )

    lowcovBuf = ReadBuffer( filterReads( readGenerator( lowcovbam, chromosome ), 
                                          minmapq=minmapq, minanchor=minanchor ) )

    counts = {}  # the distribution of non-ref allele coverage, for single-allele sites

    num = 0
    numtried = 0

    try:
        for motif in motifgenerator:

            chrom, pos, tlen, tunit = motif

            start = pos - minanchor
            end = pos + tlen + minanchor + 1

            reads = lowcovBuf.retrieve(start,end)

            numtried += 1

            cov = len(reads)
            if cov < 4 or cov > coverage: continue

            haps = haplotypes( reads, pos, pos+tlen )
            key = "%s:%s" % (tunit, tlen)
            
            if len(haps)>2:
                # aggregate the minor alleles
                alleles = sorted( [(count,hap) for hap,count in haps.iteritems()] )
                minors = sum( count for (count,hap) in alleles[:-1] )
                haps = { alleles[0][1]: minors, 
                         alleles[-1][1]: alleles[-1][0] }

            # if two alleles are present, but not the reference allele, map the major allele to the reference
            if len(haps)==2 and 0 not in haps:
                alleles = sorted( [(count,hap) for hap,count in haps.iteritems()] )
                haps = { alleles[0][1]: alleles[0][0], 
                         0: alleles[1][0] }
                
            # calculate non-ref allele count
            count = sum( count for hap,count in haps.iteritems() if hap != 0 )

            histogram = counts.get(key,[0]*((coverage-3)*(coverage+1)))    # include hom alt counts
            histogram[(cov-4)*(coverage+1)+count] += 1
            counts[key] = histogram

            num += 1
            if num > maxprocessedmotifs: break

    except Exception:
        raise

    #  create aggregate counts
    for key in set( counts.keys() ):
        tunit, tlen = key.split(':')
        newkey = "%s:%s" % (len(tunit),tlen)
        hist = counts.get(newkey, [0]*len(counts[key]))
        for idx,c in enumerate(counts[key]):
            hist[idx] += c
        counts[newkey] = hist

    # remove keys with small counts
    for key in set( counts.keys() ):
        if sum(counts[key]) < min_tot_count:
            del counts[key]

    report( counts, coverage )



def report( counts, coverage ):

    # build output
    output = []
    for key in counts.keys():
        tunit, tlen = key.split(':')
        tlen = int(tlen)
        try:      tunit = int(tunit)
        except Exception:   pass

        # fit model to counts
        N00, N01, N11, epsilon, beta = fitmodel( counts[key], coverage )

        output.append( (tunit, tlen, "%s\t%s\t%s\t%1.6f\t%1.6f\t%1.6f" % (tunit,
                                                                          tlen,
                                                                          counts[key],
                                                                          N01/(N00+N01+N11+1e-10),
                                                                          beta,
                                                                          epsilon ) ) )
    output.sort()
    for u,l,line in output:
        print line


def parse_counts( infile ):

    counts = {}
    for line in infile:
        if line.startswith('#'):
            continue
        elts = line[:-1].split('\t')
        if len(elts) < 3:
            continue
        unit = elts[0]
        length = int(elts[1])
        count = eval(elts[2])
        try:
            # convert numbers into ints; leave tandem units alone
            unit = int(unit)
        except Exception:
            pass

        key = "%s:%s" % (unit,length)
        count0 = counts.get(key,None)
        if count0:
            for idx,c in enumerate(count0):
                count[idx] += c
        counts[key] = count

    return counts


def multimodel(pars, counts, maxcoverage):

    ll = 0
    N = float(sum(counts))
    for i in range(0, len(counts), maxcoverage+1):
        cov = (i // (maxcoverage+1)) + 4
        ll += model( pars, counts[i:i+cov+1], N )
    return ll
                        

def model( pars, counts, N ):
    """ computes rate, lambda of the possible outcomes --- to score a Poisson approximation. """
    
    cov = len(counts)-1    # for coverage C, counts array includes counts for 0,..,C
    lambdas_noerr = [0] * (cov+1)
    lambdas = [0] * (cov+1)
    
    cov_scaling = sum(counts)/N
    N00, N01, N11, epsilon, beta = pars
    N00 *= cov_scaling
    N01 *= cov_scaling
    N11 *= cov_scaling

    # hom ref model
    lambdas_noerr[0] = N00
    
    # het model
    for k in range(0,cov+1):
        lambdas_noerr[k] += N01 * ncr(cov,k) * math.pow(beta,k) * math.pow(1-beta,cov-k)

    # hom alt model
    lambdas_noerr[cov] = N11

    # error model
    for k in range(0, cov+1):
        # 0 errors
        lambdas[k] += math.pow(1-epsilon,cov) * lambdas_noerr[k]

        # 1 error
        factor = cov*epsilon*math.pow(1-epsilon,cov-1)
        if k>0:
            lambdas[k-1] += k*factor*lambdas_noerr[k]/cov
        if k<cov:
            lambdas[k+1] += (cov-k)*factor*lambdas_noerr[k]/cov

        # 2 errors
        factor = cov*(cov-1)*0.5*epsilon*epsilon*math.pow(1-epsilon,cov-2)
        if k>1:
            lambdas[k-2] += k*(k-1)*factor*lambdas_noerr[k]/(cov*(cov-1))
        if k<cov-1:
            lambdas[k+2] += (cov-k)*(cov-1-k)*factor*lambdas_noerr[k]/(cov*(cov-1))
        # one error on a ref background, one on an alt background
        lambdas[k] += 2*k*(cov-k)*factor*lambdas_noerr[k]/(cov*(cov-1))   
    
    # likelihood
    ll = 0
    for k in range(0,cov+1):
        lmda = lambdas[k]
        ll += counts[k]*math.log(lmda + 1e-10) - lmda - logstirling(counts[k])

    return ll
        


def fitmodel( counts, coverage ):
    """ fits a (homref + error) + het model to count data """
        
    N00, N01, N11 = 0,0,0
    for i in range(0, len(counts), coverage+1):
        N00 += float(counts[i])
        N01 += float(sum(counts[i+1:i+coverage]))
        N11 += float(counts[i+coverage])
    eps = 0.001
    beta = 0.5

    pars = [N00,N01,N11,eps,beta]
    dpars = [0.05,0.05,0.05,0.05,0.05]
    minpars = [0.01,0.01,0.01,1e-8,0.35]
    maxpars = [1e10,1e10,1e10,0.2,0.65]
    ddpars = 0.9

    k = 0
    ll = multimodel( pars, counts, coverage )
    change = 1
    while sum(dpars)>0.001 and (change + k)>0:
        if k == 0: change = 0
        parsplus = pars[:]
        parsminus = pars[:]
        parsplus[k] *= 1.0 + dpars[k]
        parsminus[k] /= 1.0 + dpars[k]
        if parsplus[k] < maxpars[k]:
            llplus = multimodel( parsplus, counts, coverage )
        else:
            llplus = ll
        if parsminus[k] > minpars[k]:
            llminus = multimodel( parsminus, counts, coverage )
        else:
            llminus = ll
        if ll >= max(llplus, llminus):
            dpars[k] *= ddpars
            change += 1
        elif llplus > max(ll, llminus):
            pars[k] = parsplus[k]
            ll = llplus
            change += 1
        else:
            pars[k] = parsminus[k]
            ll = llminus
            change += 1
        k = (k+1) % len(pars)

    return pars
        

#
# main
#

if len(sys.argv) not in [1,2,4,5]:
    print "Usage: %s chromosome motiffile lowcovbam [maxcoverage]" % sys.argv[0]
    print "Usage: %s [maxcoverage] < output"
    sys.exit(1)

coverage = 5
minmapq = 30
minanchor = 5

if len(sys.argv) in [1,2]:
    if len(sys.argv) == 2:
        coverage= int(sys.argv[1])
    counts = parse_counts( sys.stdin )
    report( counts, coverage )
    sys.exit(0)
    

chromosome, motiffilename, lowcovbam = sys.argv[1:4]
if len(sys.argv) == 5:
    coverage = int(sys.argv[4])
    if coverage > 10:
        thinner = 1.0

motiffile = filez.open(motiffilename)

print "# chromosome      \t",chromosome
print "# bamfile         \t",lowcovbam
print "# motifs          \t",motiffilename
print "# maxcoverage     \t",coverage
print "# processed motifs\t",maxprocessedmotifs
print "# thinner         \t",thinner
print "# min_tot_count   \t",min_tot_count
print "# minanchor       \t",minanchor
print "# minmapq         \t",minmapq

estimateErrorRate( chromosome, motiffile, lowcovbam,
                   minmapq=minmapq, minanchor=minanchor, coverage=coverage)

