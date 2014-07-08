"""
annotate calls that are present in a secondary VCF file
todo: filter for genotype quality
"""

import vcf
import filez
import sys
import getopt

###################################################################################################

def wrap(line, length=80):
    """
    Free function to wrap a long line. If the line is longer than 'length', then
    the remainder is wrapped into another line. I think this can be replaced by the
    textwrap.wrap function.
    """
    if len(line) <= length:
        return line

    if len(line) == 0:
        return ""

    while length>=0 and line[length] != " ":
        length -= 1

    start = length

    while length>=0 and line[length] == " ":
        length -= 1

    if length>=0:
        return line[:length+1] + "\n" + wrap(line[start+1:])

    while length<len(line) and line[length] != " ":
        length += 1

    return line[:length] + "\n" + wrap(line[length+1:])

###################################################################################################

class VCForder():
    """
    Class to compare the position of records from a
    sorted VCF file.
    """
    def __init__(self, chromosomes):
        """
        Initialize with a sorted list of chromosomes
        """
        self._maxchromsize = 1000000000   # 1 Gb should do
        self._chromdict = {}

        for i,c in enumerate(chromosomes):
            self._chromdict[c] = i

    def cmp(self, x, y):
        """
        cmp function for sorting
        """
        c1 = x["chrom"]
        c2 = y["chrom"]

        if c1 != c2:
            return self._chromdict[c1] - self._chromdict[c2]
        else:
            return x["pos"] - y["pos"]

    def pos(self, x):
        """
        Generalized position
        """
        return self._maxchromsize*self._chromdict[x["chrom"]] + x["pos"]

    def distance(self, x, y):
        """
        Distance
        """
        return abs(self.pos(x) - self.pos(y))

###################################################################################################

def filter( vcfdata, filter ):
    for v in vcfdata['filter']:
        if v not in filter:
            return True
        else:
            return False

###################################################################################################

# merges records from two VCF files
# yields tuples (n,rec) where n=0, 1 for stream1, stream2
def VCFmerge( stream1, stream2, vcforder, filter1=[], filter2=[] ):
    try:
        r1, r2 = None, None  # to identify stream causing the exception
        r1 = stream1.next()
        r2 = stream2.next()
        while True:
            if vcforder.cmp( r1, r2 ) < 0:
                if not filter(r1, filter1): yield (0,r1)
                r1 = None
                r1 = stream1.next()
            else:
                if not filter(r2, filter2): yield (1,r2)
                r2 = None
                r2 = stream2.next()
    except StopIteration:
        if r1 == None:
            if r2 == None: raise StopIteration
            if not filter(r2, filter2): yield (1,r2)
            for r2 in stream2:
                if not filter(r2, filter2): yield (1,r2)
        else:
            if not filter(r1, filter1): yield (0,r1)
            for r1 in stream1:
                if not filter(r1, filter1): yield (0,r1)

###################################################################################################

# returns float encoding how well two calls match
# also returns the relevant pair of alternate alleles (1=allele1, etc.)
#
# 0: sequence, length, type, pos matches
# 1: length, type, pos matches
# 2: type, pos matches
# 3: only pos matches
# above + k/(maxdiff+1): position mismatch by k nucleotides
#
# 4: mismatch
#
# Based on an idea by Iain Mathieson
#
def cmpcall( call1, call2, maxdistance, vcforder ):
    r1, r2 = call1["ref"], call2["ref"]
    posdiff = vcforder.pos( call2 ) - vcforder.pos( call1 )
    distance = abs(posdiff)
    if distance > maxdistance: return 4, (r1, r2)
    # try all pairs of alternate alleles
    mincmp, alleles = 4, None
    for i1, a1 in enumerate(call1["alt"]):
        for i2, a2 in enumerate(call2["alt"]):
            code = distance / (maxdistance+1.0) + _cmpcall(r1,a1,r2,a2,posdiff)
            if code < mincmp:
                mincmp = code
                alleles = i1+1, i2+1
    return mincmp, alleles

###################################################################################################

# bring original and replacement sequence into normal form, apart from rotation
# due to position mismatch.  For insertions or deletions one of r, a will be empty,
# but for replacement calls both will be nonempty
def _truncate( r, a, pos ):
    while len(r)>0 and len(a)>0 and r[-1] == a[-1]:
        r = r[:-1]
        a = a[:-1]
    while len(r)>0 and len(a)>0 and r[0] == a[0]:
        r = r[1:]
        a = a[1:]
        pos = pos+1
    return r, a, pos

###################################################################################################

# rotate variant call assuming a repetitive reference.  Diff = relative position of current call
def _rotate( r, a, diff ):
    if len(r)>0 and len(a)>0:
        if diff == 0: return r, a
        return "X"*len(r), "Y"*len(a)  # do not rotate replacement calls; ensure mismatch
    if len(r)>0:
        a, r = _rotate(a, r, diff)
        return r, a
    # now len(r) == 0
    return r, a[-(diff % len(a)):] + a[:-(diff % len(a))]

###################################################################################################

def _cmpcall( r1, a1, r2, a2, diff ):
    r1, a1, p1 = _truncate(r1,a1,0)
    r2, a2, p2 = _truncate(r2,a2,diff)
    r2, a2 = _rotate( r2,a2,p2-p1)
    if r1 == r2 and a1 == a2: return 0  # perfect match
    if len(r1) == len(r2) and len(a1) == len(a2): return 1 # sequence mismatch
    if cmp(len(r1),len(a1)) == cmp(len(r2),len(a2)): return 2 # type match (ins/del/snp)
    return 3

###################################################################################################

# forms pairs of best-matching records out of a merged stream
def VCFpair( stream, vcforder, maxdistance=0, ignoretype=False, ignorelength=False, ignoresequence=False ):
    window = []
    # encode criteria into match code
    if ignoretype: matchcode = 3.9999
    elif ignorelength: matchcode = 2.9999
    elif ignoresequence: matchcode = 1.9999
    else: matchcode = 0.9999
    streamgood = True
    while streamgood:
        # invariant: window contains calls over no more than 2*maxdistance nucleotides
        # (and if the stream has run out, it is empty)
        try:
            record = stream.next()
            window.append(record)
        except StopIteration:
            streamgood = False
        while (len(window)>0 and
               ((len(window)>1 and vcforder.distance( window[0][1], window[-1][1] ) > 2*maxdistance) or
                (not streamgood))):
            # Make invariant true again by removing calls at head of list.
            # Any call that matches by current criteria with another in the window is paired
            # However, find reciprocal best match to avoid a spurious matching taking out a good pair
            curpair = [0]
            oldpair = []
            while set(curpair) != set(oldpair):
                oldpair = curpair
                curpair = [curpair[-1]]
                bestmatch = None
                bestcode = 4
                for i in range(len(window)):
                    j = curpair[0]
                    if i == j: continue
                    # check that we're comparing between VCF files
                    if window[j][0] != window[i][0]:
                        code, alleles = cmpcall( window[j][1], window[i][1], maxdistance, vcforder )
                        if code < min(matchcode, bestcode):
                            bestmatch = i
                            bestcode = code
                            bestalleles = alleles
                if bestmatch != None: curpair.append( bestmatch )
            result = [None, None]
            if len(curpair)==2:
                # found a best pair.  Store VCF line, and matching allele
                result[ window[oldpair[0]][0] ] = (window[oldpair[0]][1], bestalleles[0])
                result[ window[oldpair[1]][0] ] = (window[oldpair[1]][1], bestalleles[1])
                yield result
                del window[ max(oldpair) ]
                del window[ min(oldpair) ]
            else:
                # call at head does not match another
                result[ window[oldpair[0]][0] ] = (window[oldpair[0]][1], -1)
                yield result
                del window[oldpair[0]]
        # invariant is true again

###################################################################################################

# 0=r/r, 1=r/a, 2=a/a, 3=other combination or no call
def allelecode( g, allele ):
    if 'GT' not in g: raise ValueError("Allele call does not have GT field")
    g = g['GT']
    if len(g[0]) == 1:
        if g[0][0] == 0: return 0,1
        if g[0][0] == allele: return 1,1
        return 3,1
    if not (g[0][0] in [0,allele] and g[0][2] in [0,allele]): return 3,2
    return (g[0][0] == allele) + (g[0][2] == allele), 2

###################################################################################################

# number of occurrences of an allele among a set of shared samples
def countalleles( call, sharedsamples, allele, bins, homopolbin, microsatbin ):
    count, popsize = 0, 0
    for s in sharedsamples:
        if s not in call: raise ValueError("Expected sample %s in call" % s)
        code, ploidy = allelecode( call[s], allele )
        if code in [0,1,2]:
            count += code
            popsize += ploidy
    # calculate histogram bin
    if popsize>0: bin = int( bins*count / (popsize+1.0) )
    else:         bin = 0
    hr = call["info"].get("HR",[0])[0]
    tr = call["info"].get("TR",[0])[0]
    if hr>=tr:
        idx = len([c for c in homopolbin if c < hr])
        if idx>0: return count, (bin, "r"+str(idx-1))
        return count, (bin, "0")
    idx = len([c for c in microsatbin if c < tr])
    if idx>0: return count, (bin, "u"+str(idx-1))
    return count, (bin, "0")

###################################################################################################

def genotypeconcordance( call1, call2, allele1, allele2, sharedsamples ):
    gtconcordance = [ [0,0,0,0] for i in range(4) ]
    for s in sharedsamples:
        if s not in call1 or s not in call2: raise ValueError("Expected sample %s in call" % s)
        code1, ploidy = allelecode( call1[s], allele1 )
        code2, ploidy = allelecode( call2[s], allele2 )
        gtconcordance[ code1 ][ code2 ] += 1
    return gtconcordance

###################################################################################################

def updatepairconcordance( call1, call2, allele1, allele2, samples1, samples2, pairconcordance ):
    gt1 = [ allelecode( call1[s], allele1 )[0] for s in samples1 ]
    gt2 = [ allelecode( call2[s], allele2 )[0] for s in samples2 ]
    for i1 in range(len(gt1)):
        for i2 in range(len(gt2)):
            if gt1[i1]<3 and gt2[i2]<3:
                pairconcordance[i1][i2][gt1[i1]][gt2[i2]] += 1

###################################################################################################

def concordances( s1, s2, pairconcordance ):
    pc = pairconcordance[s1][s2]
    homcc = pc[0][0] / (0.01+pc[0][0]+pc[0][1]+pc[0][2])
    hetcc = pc[1][1] / (0.01+pc[1][0]+pc[1][1]+pc[1][2])
    altcc = pc[2][2] / (0.01+pc[2][0]+pc[2][1]+pc[2][2])
    conc = (pc[0][0] + pc[1][1] + pc[2][2]) / (0.01+sum(pc[0][:3])+sum(pc[1][:3])+sum(pc[2][:3]))
    hh1 = (pc[1][0] + pc[1][1] + pc[1][2]) / (0.01+pc[0][0]+pc[0][1]+pc[0][2]+pc[2][0]+pc[2][1]+pc[2][2])
    hh2 = (pc[0][1] + pc[1][1] + pc[2][1]) / (0.01+pc[0][0]+pc[1][0]+pc[2][0]+pc[0][2]+pc[1][2]+pc[2][2])
    return conc, homcc, hetcc, altcc, hh1, hh2

###################################################################################################

def pair_conc_results_1( samples1, samples2, pairconcordance ):
    s2dict = dict( ((s,n) for n,s in enumerate(samples2)) )
    for n1, s1 in enumerate(samples1):
        result = []
        n2 = s2dict.get(s1,-1)
        if n2>-1:
            conc, homcc, hetcc, altcc, hh, dummy = concordances(n1,n2,pairconcordance)
            concs2 = s1
            result += map(lambda x: "%0.3f"%x, [conc,homcc,hetcc,altcc])
        else:
            conc, homcc, hetcc, altcc, hh = [0.0]*5
            concs2 = "*"
            result += ["*"]*4
        for n2, s2 in enumerate(samples2):
            conc0, homcc0, hetcc0, altcc0, hh0, dummy = concordances(n1,n2,pairconcordance)
            if conc0 > conc or (conc0 == conc and s2 == s1):
                concs2, conc, homcc, hetcc, altcc, hh = s2, conc0, homcc0, hetcc0, altcc0, hh0
        if concs2 == s1:
            result += ["*"]*5
        else:
            result += [concs2] + map(lambda x: "%0.3f"%x, [conc, homcc, hetcc, altcc])
        yield [s1,"%0.3f"%hh] + result

###################################################################################################

def pair_conc_results_2( samples1, samples2, pairconcordance ):
    s1dict = dict( ((s,n) for n,s in enumerate(samples1)) )
    for n2, s2 in enumerate(samples2):
        result = []
        n1 = s1dict.get(s2,-1)
        if n1>-1:
            conc, homcc, hetcc, altcc, dummy, hh = concordances(n1,n2,pairconcordance)
            concs1 = s2
            result += map(lambda x: "%0.3f"%x, [conc,homcc,hetcc,altcc])
        else:
            conc, homcc, hetcc, altcc, hh = [0.0]*5
            concs1 = "*"
            result += ["*"]*4
        for n1, s1 in enumerate(samples1):
            conc0, homcc0, hetcc0, altcc0, dummy, hh0 = concordances(n1,n2,pairconcordance)
            if conc0 > conc or (conc0 == conc and s1 == s2):
                concs1, conc, homcc, hetcc, altcc, hh = s1, conc0, homcc0, hetcc0, altcc0, hh0
        if concs1 == s2:
            result += ["*"]*5
        else:
            result += [concs1] + map(lambda x: "%0.3f"%x, [conc, homcc, hetcc, altcc])
        yield [s2,"%0.3f"%hh] + result

###################################################################################################

def cross_concordances( samples1, samples2, pairconcordance ):

    print "VCF1, cross-concordance with VCF2:"
    print "\t".join( ["Sample","H/H","var","hom","het","alt","Match","var","hom","het","alt"] )
    for line in pair_conc_results_1( samples1, samples2, pairconcordance ):
        print "\t".join(line)
    print
    print "VCF2, cross-concordance with VCF1:"
    print "\t".join( ["Sample","H/H","var","hom","het","alt","Match","var","hom","het","alt"] )
    for line in pair_conc_results_2( samples1, samples2, pairconcordance ):
        print "\t".join(line)
    print
    print "Columns: H/H=het/hom ratio (primary VCF)"
    print "         var=variant concordance; hom/hom + het/het + var/var / total"
    print "             * if no matching sample in secondary VCF"
    print "         hom=homozygous concordance; hom/hom / (hom/ref + hom/het + hom/alt)"
    print "         het=heterozygous concordance"
    print "         alt=homozygous alternative concordance"
    print "         Match=best matching sample; * if matches same sample"
    print "         var,hom,het,alt: concordances with best-matching sample"
    print

###################################################################################################

def pretty( b, homopolbin, microsatbin ):
    j, contextbin = b
    if contextbin == "0": return "f"+str(j)
    if contextbin.startswith("r"): return "f%s;HR>=%s" % (j, homopolbin[int(contextbin[1:])])
    return "f%s;TR>=%s" % (j,microsatbin[int(contextbin[1:])])

def help():
    print "Compare two VCF files"
    print "Usage: vcf-compare.py [options]"
    print
    print " -h, --help           This help"
    print " -oF, --output F      Send output (calls from set 1 annotated by match status) to file F"
    print " -iF, --input F       Read input 1 from file F (default stdin)"
    print " -cF, --input2 F      Read input 2 from file F (required)"
    print " -n, --nosamples      Ignore sample information, just compare calls"
    print " -sR, --select R      Select region or regions to compare (chr:start-end,chr:start-end,...)"
    print " -mK, --match K       Match stringency (0 [default] exact; 1: ignore sequence; 2: ignore length; 3: ignore type)"
    print " -pN, --position N    Allow position to mismatch N positions (default 0)"
    print " -fF1,F2,F3           Allow FILTER values F1, F2, F3 for VCF1 (default None)"
    print " -gF1,F2,F3           Allow FILTER values F1, F2, F3 for VCF2 (default None)"
    print " -bN, --bins N        Bin by N frequency bins, using input 1 allele counts (default 1)"
    print " -rN[,N...]           Bin by homopolymer run length using thresholds N,N,..."
    print " -uN[,N...]           Bin by microsatellite length using thresholds N,N,..."
    print " -t, --samplestats    Compute per-sample statistics"
    print " -a, --callstats      Compute per-call statistics"
    print " -3                   Set default input version to 3.3"
    print " -F                   Fast and sloppy GT parsing"
    print " -xE, --ignore E      Ignore error/warning E"
    print
    print " Note: When positions do not match exactly, indel sequences are rotated assuming that the reference is repetitive."
    # todo: select single chromosome
    # todo: implement filtering
    print
    print "Possible errors are:"
    v = vcf.VCF()
    print wrap( ", ".join([v.split(':')[0] for v in v._errors.values()]) )
    sys.exit(1)

###################################################################################################

def main():
    infile, fname = sys.stdin, "<stdin>"
    secondaryin = None
    outfile = None
    inversion = 40
    nosamples = False
    stringency = 0
    maxdistance = 0
    regions = None
    fastGT = False
    samplestats = False
    callstats = False
    homopolbin = []
    microsatbin = []
    allowfilter1, allowfilter2 = [], []
    bins = 1
    v = vcf.sortedVCF()
    try:
        if len(sys.argv)==1: sys.argv.append("-h")
        opts, args = getopt.getopt(sys.argv[1:], 
                                   "ho:i:x:c:l:d:3nm:p:f:g:s:tr:u:b:aF", 
                                   ["help","output","input","input2","ignore","nosamples","match","position","select","samplestats","bins=","callstats"])
    except getopt.GetoptError, err:
        print str(err)
        help()
    for o, a in opts:
        if o in ["-h","--help"]:
            help()
            sys.exit()
        elif o in ["-n","--nosamples"]:         nosamples = True
        elif o in ["-m","--match"]:             stringency = int(a)
        elif o in ["-p","--position"]:          maxdistance = int(a)
        elif o in ["-x","--ignore"]:            v.ignoreerror(a)
        elif o in ["-o","--output"]:            outfile = open(a,'w')
        elif o in ["-b","--bins"]:              bins = int(a)
        elif o in ["-s","--select"]:            regions = vcf.parse_regions(a)
        elif o in ["-t","--samplestats"]:       samplestats = True
        elif o in ["-r"]:                       homopolbin = map(int,a.split(','))
        elif o in ["-u"]:                       microsatbin = map(int,a.split(','))
        elif o in ["-a","--callstats"]:         callstats = True
        elif o in ["-F"]:                       fastGT = True
        elif o in ["-i","--input"]:             
            infile = filez.open(a,'r')
            fname = a
        elif o in ["-c","--input2"]:
            secondaryin = filez.open(a,'r')
            secondaryfname = a
        elif o == "-3":                         inversion = 33
        elif o == "-f":                         allowfilter1 = a.split(',')
        elif o == "-g":                         allowfilter2 = a.split(',')

    if not secondaryin:
        raise ValueError("Need two input VCF files")

    v.setversion(inversion)
    v._fastGT = fastGT

    # pre-parse first 1000 lines of each file, to catch formatting problems early
    v0 = vcf.VCF( v, lines=1000, _fastGT=fastGT )
    for dummy in v0.parse( infile ): pass
    infile = filez.open(fname,'r')

    v0 = vcf.VCF( v, lines=1000, _fastGT=fastGT )
    for dummy in v0.parse( secondaryin ): pass
    secondaryin = filez.open(secondaryfname, 'r')

    # process data
    v.setregions(regions)
    secondary = vcf.sortedVCF(v)
    secondary._fastGT = fastGT
    vcfstream1 = v.parse( infile )
    vcfstream2 = secondary.parse( secondaryin )

    samples1 = v.getsamples()
    samples2 = secondary.getsamples()
    sharedsamples = list( set.intersection( set(samples1), set(samples2) ) )

    if len(sharedsamples) == 0 and not nosamples:
        raise ValueError("Files do not share samples (consider option -n)")

    for f in allowfilter1:
        if f not in v.getfilter(): print "#WARNING: filter value %s not defined in header" % f

    for f in allowfilter2:
        if f not in secondary.getfilter(): print "#WARNING: filter value %s not defined in header" % f

    if samplestats:
        pairsamples1 = samples1
        pairsamples2 = samples2
        pairconcordance = [ [ [ [0,0,0] for i in range(3) ] for s2 in samples2 ] for s1 in samples1 ]

    if not nosamples:
        samples1 = sharedsamples
        samples2 = sharedsamples

    if outfile:
        vcfout = vcf.VCF(v)
        vcfout.setversion(40)
        for i in ["1","2"]:
            vcfout.getinfo()['CI'+i] = \
                vcf.FORMAT('CI'+i,vcfout.NT_NUMBER,1,"String",
                           "Concordance info VCF%s: NC/NA/AP = Not Called/No Allele in shared samples/Allele Present" % i, "")
        vcfout.getinfo()['GCC'] = \
            vcf.FORMAT('GCC',vcfout.NT_NUMBER,16,"Integer",
                       "Genotype concordance counts: HomRef(1)/HomRef(2), Het(1)/HomRef(2), HomAlt(1)/HomRef(2), Other(1)/HomRef(2), ...",
                       [0]*16)
        vcfout.getheader().append( ("comment","Secondary input file (for CI1/CI2/GCC): %s" % secondaryfname) )
        vcfout.getheader().append( ("source"," ".join(sys.argv)) )
        vcfout.writeheader( outfile )
    else:
        vcfout = None

    chromosomes = set(v.chr_order())
    chromosomes.update(secondary.chr_order())
    chromosomes = v.chr_order( list(chromosomes) )
    vcforder = VCForder(chromosomes)

    mergedstream = VCFmerge(vcfstream1, vcfstream2, vcforder, filter1=allowfilter1, filter2=allowfilter2 )

    pairstream = VCFpair( mergedstream, vcforder, maxdistance,
                          ignoretype=(stringency>=3), ignorelength=(stringency>=2), ignoresequence=(stringency>=1) )
    # setup results counters
    #  callconcordance[i][j]  i,j are 0,1 or 2 encoding:
    #  0 = site not called;
    #  1 = site called but allele not called in intersection
    #  2 = site and allele called
    #
    #  gtconcordance[i][j]   callconcordance[2][2], and i,j encode
    #  0 = ref/ref
    #  1 = alt/ref
    #  2 = alt/alt
    #  3 = other  (
    #  Here "alt" = the matching non-reference allele, "ref" is the
    #  reference, "other" a non-call or non-matching non-ref allele

    contextbins = ["0"] + ["r"+str(i) for i in range(len(homopolbin))] + ["u"+str(i) for i in range(len(microsatbin))]
    allbins = []
    for cb in contextbins:
        for j in range(bins):
            allbins.append( (j,cb) )

    callconcordance = dict([(j,[ [0,0,0] for i in range(3) ]) for j in allbins ])
    gtconcordance = dict([(j,[ [0,0,0,0] for i in range(4) ]) for j in allbins ])

    for pair in pairstream:

        if pair[1] == None:
            # VCF1 call not matched in VCF2
            call = pair[0][0]
            call["info"]["CI2"] = ["NC"]
            numalleles, bin = countalleles( call, samples1, 1, bins, homopolbin, microsatbin )
            if numalleles == 0:
                # VCF1 call with no alleles in intersection
                callconcordance[bin][1][0] += 1
                call["info"]["CI1"] = ["NA"]
            else:
                callconcordance[bin][2][0] += 1
                call["info"]["CI1"] = ["AP"]
            if outfile: vcfout.write_data( outfile, call )
        elif pair[0] == None:
            # VCF2 call not matched in VCF1
            call = pair[1][0]
            numalleles, dummy = countalleles( call, samples2, 1, bins, homopolbin, microsatbin )
            if numalleles == 0:
                # VCF2 call with no alleles in intersection
                callconcordance[(0,"0")][0][1] += 1
            else:
                callconcordance[(0,"0")][0][2] += 1
        else:
            # Paired call
            call1, call2 = pair[0][0], pair[1][0]
            allele1, allele2 = pair[0][1], pair[1][1]
            numalleles1, bin = countalleles( call1, samples1, allele1, bins, homopolbin, microsatbin )
            numalleles2, dummy = countalleles( call2, samples2, allele2, bins, homopolbin, microsatbin )
            callconcordance[bin][ 1+(numalleles1>0) ][ 1+(numalleles2>0) ] += 1
            call1["info"]["CI1"] = [["NA","AP"][numalleles1>0]]
            call1["info"]["CI2"] = [["NA","AP"][numalleles2>0]]
            gc = genotypeconcordance( call1, call2, allele1, allele2, sharedsamples )
            call1["info"]["GCC"] = gc[0]+gc[1]+gc[2]+gc[3]
            for i1, gc1 in enumerate(gc):
                for i2, n in enumerate(gc1):
                    gtconcordance[bin][i1][i2] += n
            if samplestats: updatepairconcordance( call1, call2, allele1, allele2, pairsamples1, pairsamples2, pairconcordance )
            if outfile: vcfout.write_data( outfile, call1 )

    print "\nConcordance summary"
    print "-------------------\n"
    print "\tVCF1: ",fname
    print "\tVCF2: ",secondaryfname
    print "\n\nCall concordance:"
    print "-----------------\n"
    for b in allbins:
        print "\t%s\t%s\t%s\t%s\t%s" % ("(bin %s)"%pretty(b,homopolbin,microsatbin),"VCF2:","NC","NA","AP")
        for i in range(3):
            print "\t%s\t%s\t%s\t%s\t%s" % (("VCF1:",["NC","NA","AP"][i])+tuple(callconcordance[b][i]))
        print
    print "\tNC=no call, NA=no alleles in sample intersection, AP=allele present in intersection"
    print "\tWith -n, NA and AP report allele presence or absence within all samples rather than the intersection"
    print "\n\nGenotype concordance (AP/AP calls):"
    print "-----------------------------------\n"
    for b in allbins:
        print "\t%s\t%s\t%s\t%s\t%s\t%s" % ("(bin %s)"%pretty(b,homopolbin,microsatbin),"VCF2:","HomRef","Het","HomAlt","Other")
        for i in range(4):
            print "\t%s\t%s\t%s\t%s\t%s\t%s" % (("VCF1:",["HomRef","Het","HomAlt","Other"][i])+tuple(gtconcordance[b][i]))
        print
    print "\t(HomRef: homozygous reference; Het: ref/alt where alt=matching non-ref allele,"
    print "\t HomAlt: alt/alt, Other: no call, or unmatched non-ref allele)"
    print
    print "Genotype quality (AP/AP calls)"
    print "------------------------------\n"
    print "Freqbin\tHomRef\t\t\tHet\t\t\tHomAlt"
    print "\tTrue\tTotal\tFrac\tTrue\tTotal\tFrac\tTrue\tTotal\tFrac"
    for b in allbins:
        print ("\t".join(["%s"]*10)) % (pretty(b,homopolbin,microsatbin),
                                        gtconcordance[b][0][0],sum(gtconcordance[b][0][0:3]),
                                        "%1.3f" % (gtconcordance[b][0][0]/(0.01+sum(gtconcordance[b][0][0:3]))),
                                        gtconcordance[b][1][1],sum(gtconcordance[b][1][0:3]),
                                        "%1.3f" % (gtconcordance[b][1][1]/(0.01+sum(gtconcordance[b][1][0:3]))),
                                        gtconcordance[b][2][2],sum(gtconcordance[b][2][0:3]),
                                        "%1.3f" % (gtconcordance[b][2][2]/(0.01+sum(gtconcordance[b][2][0:3]))))
    print
    print "(Assumes VCF1 is true.  Non-calls and unmatched non-ref alleles are not considered.)\n"
    if samplestats:
        print
        cross_concordances( pairsamples1, pairsamples2, pairconcordance )

###################################################################################################

if __name__ == "__main__":
    main()
