# obtains indels from axt files

import sys
from collections import namedtuple

# set to 33 or 40
vcfversion = 40

# set to True if point mutations also included
snp = True

# maximum number of differences accepted when shifting the gap
maxdiffs = 2

# indel and alignment types
Indel = namedtuple('Indel','chrom pos ref alt start end')
Alignment = namedtuple('Alignment','chr1 pos1 end1 chr2 pos2 end2 strand score align1 align2')


# get v3.3 representation
def encode(indel):
    # point to first non-matching column from left
    i = 0
    while i<min(len(indel.alt), len(indel.ref)) and indel.alt[i] == indel.ref[i]: 
        i += 1
    # point to first non-matching coumn from right
    j = 0
    while len(indel.alt)-j-1 >= i and len(indel.ref)-j-1 >= i and indel.alt[-j-1] == indel.ref[-j-1]: 
        j += 1
    # pure indel?
    if i+j == min(len(indel.ref),len(indel.alt)):
        # success
        length = max(len(indel.ref),len(indel.alt)) - i - j
        if len(indel.alt) < len(indel.ref):
            return indel.pos + i, "D" + str( length )
        else:
            return indel.pos + i, "I" + indel.alt[i:len(indel.alt)-j]
    else:
        # mixed type
        return indel.pos + i, "R" + indel.ref[i:len(indel.ref)-j] + ":" + indel.alt[i:len(indel.alt)-j]

def testencode():
    
    indels = [ Indel("chr1",1000,"AA",""),
               Indel("chr1",1000,"AT","A"),
               Indel("chr1",1000,"TA","A"),
               Indel("chr1",1000,"TAT","A"),
               Indel("chr1",1000,"","AT"),
               Indel("chr1",1000,"A","TA"),
               Indel("chr1",1000,"A","AT"),
               Indel("chr1",1000,"A","TAT") ]
    for i in indels:
        print i, encode(i)

def axt2alignment(stream):

    while True:
        line = stream.readline()
        if line.startswith('#'): continue
        if line == "": raise StopIteration
        idx, chrom, pos, posend, chrom2, pos2, posend2, strand, score = line[:-1].split(' ')
        line1 = stream.readline()[:-1].upper()
        line2 = stream.readline()[:-1].upper()
        line3 = stream.readline()
        yield Alignment(chrom,int(pos),int(posend),chrom2,int(pos2),int(posend2),strand,int(score),line1,line2)

def alignment2indel(stream):

    for alignment in stream:

        # first generate a "genotypeable" record
        yield Indel(alignment.chr1, alignment.pos1, 'N', 'N', alignment.pos1, alignment.end1)

        i = 0
        pos1 = alignment.pos1
        while i < len(alignment.align1):
            if alignment.align1[i] != '-' and alignment.align2[i] != '-':
                if snp:
                    if alignment.align1[i] != alignment.align2[i]:
                        var = Indel(alignment.chr1, pos1, alignment.align1[i], alignment.align2[i], -1, -1)
                        yield var
                i += 1
                pos1 += 1
                continue

            j = i
            pos = pos1
            while i < len(alignment.align1) and (alignment.align1[i] == '-' or alignment.align2[i] == '-'):
                if alignment.align1[i] != '-': pos1 += 1
                i += 1

            # j:i is the gap
            jmax, imax, a1max, a2max = movegap( alignment.align1, alignment.align2, j, i, 1, maxdiffs)
            jmin, imin, a1min, a2min = movegap( alignment.align1, alignment.align2, j, i, -1, maxdiffs)

            # add anchor (VCF V4.0)
            if j>0:
                j -= 1
                pos -= 1

            indel = Indel(alignment.chr1, pos,
                          alignment.align1[j:i].replace('-',''),
                          alignment.align2[j:i].replace('-',''),
                          jmin-j+pos, jmax-j+pos)
                
            yield indel

def indel2vcf(stream, inputfilename = "<stdin>"):

    if vcfversion == 33: 
        filter = "0"
        print "##fileformat=VCFv3.3"
    else:         
        filter = "PASS"
        print "##fileformat=VCFv4.0"
    print "##source=axt2indel.py " + inputfilename
    print "##INFO=<ID=START,Number=1,Type=Integer,Description=\"Leftmost placement allowing %s mismatches\">" % maxdiffs
    print "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Rightmost placement allowing %s mismatches\">" % maxdiffs
    print "##INFO=<ID=CALLABLE,Number=0,Type=Flag,Description=\"Genotypable region\">"
    print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + inputfilename
    
    for indel in stream:
        
        # edit the 2nd allele to remove potential SNP at first position
        # however, do this only for indels, not for actual SNPs
        if len(indel.ref) > 1 or len(indel.alt) > 1:
            alt = indel.ref[0] + indel.alt[1:]
            indel = indel._replace( alt=alt )

        if vcfversion == 33:
            pos, alt = encode(indel)
            ref = indel.ref[0]
        else:
            pos, ref, alt = indel.pos, indel.ref, indel.alt

        info = ""
        if indel.start>=0: info += "START=%s;" % indel.start
        if indel.end>=0: info += "END=%s;" % indel.end
        if ref == "N" and alt == "N": info += "CALLABLE;"
        if len(info)==0: info = "."
        else:            info = info[:-1]

        print "%s\t%s\t.\t%s\t%s\t99\t%s\t%s\tGT\t1" % (indel.chrom,pos,ref,alt,filter,info)

def movegap( a1, a2, i, j, direction, maxdiff=0):

    a1 = [c for c in a1]
    a2 = [c for c in a2]
    if direction>0:
        frm = j
        to = i
    else:
        frm = i-1
        to = j-1
    cont = True
    diff = 0
    while cont and 0 <= frm < len(a1) and 0 <= to < len(a2):
        if a1[to] == '-':
            if a2[to] != a1[frm]: diff += 1
        else:
            if a1[to] != a2[frm]: diff += 1
        if a1[to] == '-' and diff <= maxdiff:
            a1[to] = a1[frm]
            a1[frm] = '-'
        elif a2[to] == '-' and diff <= maxdiff:
            a2[to] = a2[frm]
            a2[frm] = '-'
        else:
            # done
            cont = False
        if cont:
            frm += direction
            to += direction
    if direction>0:
        return to, frm, ''.join( a1 ), ''.join( a2 )
    return frm+1, to+1, ''.join( a1 ), ''.join( a2 )


def testmovegap():

    print movegap( 'aa--aa','aaaaaa', 2, 4, 1 )
    print movegap( 'aa--aa','aaaaaa', 2, 4, -1 )
    print movegap( 'aaaaaa','aa--aa', 2, 4, 1 )
    print movegap( 'aaaaaa','aa--aa', 2, 4, -1 )

    print movegap( 'at--at','atatat', 2, 4, 1 )
    print movegap( 'at--at','atatat', 2, 4, -1 )
    print movegap( 'atatat','at--at', 2, 4, 1 )
    print movegap( 'atatat','at--at', 2, 4, -1 )

    print movegap( 'at--at','atagat', 2, 4, 1 )
    print movegap( 'at--at','atagat', 2, 4, -1 )
    print movegap( 'atagat','at--at', 2, 4, 1 )
    print movegap( 'atagat','at--at', 2, 4, -1 )

    print movegap( 'at--at','atagat', 2, 4, 1, 1 )
    print movegap( 'at--at','atagat', 2, 4, -1, 1 )
    print movegap( 'atagat','at--at', 2, 4, 1, 1 )
    print movegap( 'atagat','at--at', 2, 4, -1, 1 )


# main

filename = sys.argv[1]
if len(sys.argv)>2 and sys.argv[2] == "--nosnps": snp = False
file = open(filename,'r')
align = axt2alignment(file)
indel = alignment2indel(align)
indel2vcf(indel, filename)
            



        
        
            

