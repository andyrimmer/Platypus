#
# Polarizes indel in a VCF file, using indels called in an outgroup species alignment
#
# Usage: script.py indels.vcf outgroupIndels1.vcf [outgroupIndels2.vcf [outgroupIndels3.vcf]]
#
# A polymorphic indel wrt the reference, which is also found by aligning the reference to
# one or more outgroups, is assumed to represent the ancestral state.
#
# More precisely, the criterion for calling ancestral state for an indel is:
#  The site aligns to at least 'mincallable' of the outgroups (default 1,2,2 for 1,2,3 outgroups)
#  Either of these events is true:
#   1. None of the alignable outgroups show indels in a window around the reference indel site
#   2. All of the alignable outgroups show an indel matching the reference indel within the window,
#      and no other indels in that window
#  The window is defined as 'match_shoulder' (default 5) bp left of the reference site, and
#  'match_shoulder' plus the maximum of the homopolymer run length and the indel length right of
#  the reference site.
#  Indels are considered to match if their size and type matches; their sequence is not compared.
#
# The VCF files must be sorted by chromosome and position
#
# The input vcf file can be either vcf3.3 or vcf4.0.
# WARNING: only works with vcf4.0 files for ancestral calls
#
# Fields added: AA (0 if REF=ANC, 1 if ALT=ANC, . if not called); SZ (size; negative is deletion)
#

import sys
import bisect
from collections import namedtuple
import filez

# to add prefix to chromosome identifiers in input
chromprefix = ""

# to account for offset -1 in vcf 3.3 files
vcf_offset = -1

# match positions within window +- match_shoulder
match_shoulder = 5

# minimum number of callable outgorups
mincallable = 2


Indel = namedtuple('Indel','chrom pos ref alt start end hr')



def readVCF( filename, chrom ):

    vcf = []
    callablelist = []
    print >>sys.stderr, "#### Loading vcf for choromosome",chrom
    for line in filez.open(filename,'r'):
        if line.startswith('#'): continue
        if not line.startswith( chrom + "\t"): continue
        elts = line[:-1].split('\t')
        chrom, pos, id, ref, alt, qual, filter, info = elts[:8]
        start, end, callable, hr = -1, -1, False, 1
        for e in info.split(';'):
            if e.startswith('START'): start = int(e.split('=')[1])
            if e.startswith('END'): end = int(e.split('=')[1])
            if e.startswith('CALLABLE'): callable = True
            if e.startswith('HR') or e.startswith('TR'): hr = max(hr,int(e.split('=')[1]))
        if callable:
            if len(callablelist) >0 and callablelist[-1]+1 == start:
                callablelist[-1] = end
            else:
                callablelist.append( start )
                callablelist.append( end )
        else:
            vcf.append( Indel(chrom, int(pos),ref,alt,start,end,hr) )
    print>>sys.stderr,  "#### Loaded ",len(vcf)," calls"
    vcf.sort()
    return vcf, callablelist


def getVCFcalls( vcf, first, last ):

    if len(vcf) == 0: return []
    chrom = vcf[0].chrom
    left = bisect.bisect_left( vcf, Indel(chrom,first,'','', -1, -1, -1) )
    right = bisect.bisect_right( vcf, Indel(chrom,last,'','', -1, -1, -1) )
    return vcf[left:right]


def isCallable( pos, callablelist ):

    idx = bisect.bisect_left( callablelist, pos )
    return (idx % 2) == 1


# leaves vcf4.0 calls alone
# deals with vcf3.3 calls
# also handles D<seq> calls, and with ref=='.'
def convert( ref, alt ):
    if ref == '.': ref = 'N'
    if alt.startswith('I'):
        alt = ref + alt[1:]
        return ref, alt
    if alt.startswith('D'):
        try:
            ref, alt = ref + 'N' * int(alt[1:]), ref
        except:
            ref, alt = ref + alt[1:], ref
        return ref, alt
    return ref, alt


def findMatch( anccalls, pos, typesize, hr):
    found = False
    othercall = None
    for call in anccalls:
        anctypesize = len(call.alt) - len(call.ref)
        if typesize == anctypesize:
            if typesize != 0:
                found = (found or
                         call.pos==pos or
                         call.start - match_shoulder <= pos <= max(call.end, call.start + hr) + match_shoulder)
            else:
                found = found or (call.alt == alt and call.pos == pos)
        else:
            if (typesize != 0) and (anctypesize != 0):
                if call.pos == pos or call.start - match_shoulder <= pos <= max(call.end, call.start + hr) + match_shoulder:
                    othercall = anctypesize
    return found, othercall


def processVCF( invcf, ancvcf, ancvcf2, ancvcf3, ancvcf4, window=50 ):

    invcf = filez.open(invcf,'r')
    curchrom = None
    for line in invcf:
        # copy header
        if line.startswith('#'):
            if line.startswith('#CHROM'):
                ## add info fields
                print "##INFO=<ID=AA, Number=1, Type=String, Description=""Ancestral Allele"">"
                print "##INFO=<ID=SZ, Number=1, Type=Integer, Description=""Indel size and type"">"
                print "##source=vcf-indel-polarize.py %s %s %s %s %s" % (invcf, ancvcf, ancvcf2, ancvcf3, ancvcf4)
            print line[:-1]
            continue
        # get data
        elts = line[:-1].split('\t')
        chrom, pos, id, ref, alt, qual, filter, info = elts[:8]
        pos = int(pos)
        chrom = chromprefix + chrom
        ref, alt = convert(ref,alt)
        typesize = len(alt)-len(ref)
        if typesize != 0:
            # not a snp -- adjust coordinate
            pos -= vcf_offset
        hr = 1
        for e in info.split(';'):
            if e.startswith('HR') or e.startswith('TR'): hr = max(hr,int(e.split('=')[1]))
        # load ancestral chromosome calls if required
        if curchrom != chrom:
            vcf, callablelist = readVCF( ancvcf, chrom )
            if ancvcf2: vcf2, callablelist2 = readVCF( ancvcf2, chrom )
            if ancvcf3: vcf3, callablelist3 = readVCF( ancvcf3, chrom )
            if ancvcf4: vcf4, callablelist4 = readVCF( ancvcf4, chrom )
            curchrom = chrom
        # get calls in vicinity, find matching calls
        anccalls = getVCFcalls( vcf, pos-window, pos+window )
        found, other = findMatch( anccalls, pos, typesize, hr )
        callable = isCallable( pos, callablelist )

        # do same for ancvcf2
        if ancvcf2:
            anccalls2 = getVCFcalls( vcf2, pos-window, pos+window )
            found2, other2 = findMatch( anccalls2, pos, typesize, hr )
            callable2 = isCallable( pos, callablelist2 )
        else:
            found2 = found
            callable2 = callable
            other2 = other

        # and ancvcf3
        if ancvcf3:
            anccalls3 = getVCFcalls( vcf3, pos-window, pos+window )
            found3, other3 = findMatch( anccalls3, pos, typesize, hr )
            callable3 = isCallable( pos, callablelist3 )
        else:
            found3 = found
            callable3 = callable
            other3 = other

        # and ancvcf4
        if ancvcf4:
            anccalls4 = getVCFcalls( vcf4, pos-window, pos+window )
            found4, other4 = findMatch( anccalls4, pos, typesize, hr )
            callable4 = isCallable( pos, callablelist4 )
        else:
            found4 = found
            callable4 = callable
            other4 = other

        numcallable = callable + callable2 + callable3 + callable4

        # get status from one of the called alleles
        if callable:
            found0, other0 = found, other
        elif callable2:
            found0, other0 = found2, other2
        elif callable3:
            found0, other0 = found3, other3
        elif callable4:
            found0, other0 = found4, other4
        else:
            found0, other0 = False, None

        # fill in the uncallable states with one of the others
        if not callable:
            found, other = found0, other0
        if not callable2:
            found2, other2 = found0, other0
        if not callable3:
            found3, other3 = found0, other0
        if not callable4:
            found4, other4 = found0, other0

        # if sufficiently many outgroups are callable, and found
        # status is consistent, and no other indels were found, polarize.

        if numcallable >= mincallable and (found == found2) and (found == found3) and (found == found4) and \
               (other == None) and (other2 == None) and (other3 == None) and (other4 == None):
            if found:
               info = "AA=1;SZ=%s;%s" % (-typesize, info)
            else:
                info = "AA=0;SZ=%s;%s" % (typesize, info)
        else:
            info = "AA=.;SZ=%s;%s" % (typesize, info)

        # make output
        print "\t".join([chrom, str(pos), id, ref, alt, qual, filter, info] + elts[8:])


vcfin = sys.argv[1]
vcfrefvsout = sys.argv[2]
try:
    vcfrefvsout2 = sys.argv[3]
    try:
        vcfrefvsout3 = sys.argv[4]
        try:
            vcfrefvsout4 = sys.argv[5]
        except:
            vcfrefvsout4 = None
    except:
        vcfrefvsout3 = None
except:
    vcfrefvsout2 = None
    mincallable = 1
processVCF( vcfin, vcfrefvsout, vcfrefvsout2,vcfrefvsout3, vcfrefvsout4, window=100 )

        

    
