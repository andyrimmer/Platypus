#!/usr/bin/env python

#
# Converts coordinates between several coordinate systems
#
# Todo: check for unique names
#

import sys
import os


#
# Coordinate-wise access to fasta-formatted chromosome files using the slice operator
#
class Chromosome:

    def __init__(self, filename):
        self.file = open(filename,'rb')
        self.chrom = self.file.readline()[1:-1]
        self.line0 = self.file.tell()
        self.numnucs = len(self.file.readline()[:-1])
        self.linelen = self.file.tell() - self.line0
        self.file.seek(0,2)
        totdata = self.file.tell() - self.line0
        lines = totdata // self.linelen
        self.totlen = lines * self.numnucs
        if totdata % self.linelen != 0:
            self.totlen += totdata % self.linelen - self.linelen + self.numnucs
        print "%s\t%s" % (self.chrom,len(self))

    def __getslice__(self,i,j):
        if 0 > i or j > self.totlen: raise IndexError, "Trying to access slice ["+str(i)+":"+str(j)+"], have [0:"+str(self.totlen)+"]"
        if i>=j: return ''
        startline = i//self.numnucs
        endline = (j-1)//self.numnucs
        leftover = (endline+1)*self.numnucs - j
        if leftover == 0:
            datalen = (endline-startline+1)*self.linelen
        else:
            datalen = (endline-startline+1)*self.linelen - self.linelen + self.numnucs - leftover
        self.file.seek(self.line0 + startline*self.linelen)
        return ''.join(self.file.read(datalen).splitlines())[i % self.numnucs:]

    def __getitem__(self,i):
        return self.__getslice__(i,i+1)

    def __len__(self):
        return self.totlen


#
# Open all chromosomes in a directory
#
def getchromosomes(directory):
    f = os.popen("ls -1 --color=never "+directory+"*").readlines()
    chroms = []
    for extension in [".fa",".fa.gz",".fa.dz"]:
        chroms += [ c[:-len(extension)-1] for c in f if c[-len(extension)-1:-1] == extension]   # ignore '\n' at end
    chromdict = {}
    regions = []
    for chrfile in chroms:
        chrom = chrfile.split('/')[-1]
        chromdict[chrom] = Chromosome( directory + '/' + chrom + ".fa" )
        # add region; the fifth entry makes this a 1-based (own) coordinates system
        regions.append( (chrom, 0, len(chromdict[chrom]), '+', chrom, 1) )
    return regions, chromdict


#
# Sequence utilities
#
def revcomp( sequence ):
    d = {'A':'T','C':'G','G':'C','T':'A','N':'N','-':'-'}
    s = [ d[c] for c in sequence ]
    s.reverse()
    return ''.join(s)

def codon2aminoacid(codon):
    return aa2aminoacid( codon2aa( codon ))

def aa2aminoacid(aa):
    return {'A':'Ala','C':'Cys','D':'Asp','E':'Glu','F':'Phe','G':'Gly','H':'His','I':'Ile','K':'Lys','L':'Leu',
            'M':'Met','N':'Asn','P':'Pro','Q':'Gln','R':'Arg','S':'Ser','T':'Thr','V':'Val','W':'Trp','Y':'Tyr','X':'Stop'}[aa]


def codon2aa(codon):
    return {
        'GCA':'A','GCC':'A','GCG':'A','GCT':'A','TGC':'C','TGT':'C','GAC':'D','GAT':'D',
        'GAA':'E','GAG':'E','TTC':'F','TTT':'F','GGA':'G','GGC':'G','GGG':'G','GGT':'G',
        'CAC':'H','CAT':'H','ATA':'I','ATC':'I','ATT':'I','AAA':'K','AAG':'K','TTA':'L',
        'TTG':'L','CTA':'L','CTC':'L','CTG':'L','CTT':'L','ATG':'M','AAC':'N','AAT':'N',
        'CCA':'P','CCC':'P','CCG':'P','CCT':'P','CAA':'Q','CAG':'Q','AGA':'R','AGG':'R',
        'CGA':'R','CGC':'R','CGG':'R','CGT':'R','AGC':'S','AGT':'S','TCA':'S','TCC':'S',
        'TCG':'S','TCT':'S','ACA':'T','ACC':'T','ACG':'T','ACT':'T','GTA':'V','GTC':'V',
        'GTG':'V','GTT':'V','TGG':'W','TAC':'Y','TAT':'Y','TAA':'X','TGA':'X','TAG':'X'}[codon.upper()];

#
# Read transcript file
#
# The transcript is converted into regions located on the genome.  The location of the regions are
# represented in two ways: first by the native coordinate system, (chrom, start, end), 0-based, half-open,
# and in the forward direction.  Second, in the "regional" coordinate system identified by a name (e.g. the
# transcript), an offset, and a strand relative to the native system.  The region over which this coordinate
# system is valid is implicitly defined by the segment on the chromosome.  Multiple regions sharing an
# identifier may exist, but these are supposed to be non-overlapping.  The data for these two coordinate
# systems is represented in a six-tuple:
#
#  (chrom, start, end, strand, name, offset)
#
# For transcripts, the first (coding) exon has offset 1; this also defined the frame.
#
def readtxfile( filename ):
#
# Example entry: (from ucsc tables)
#
# 2107	ENST00000360372	chr1	-	199594764	199609005	199594960	199609005	14
# 199594764,199595373,199597029,199597663,199598136,199599046,199600048,199600941,199601360,199602588,199603521,199603912,199604168,199608964,
# 199595006,199595414,199597120,199597773,199598145,199599157,199600126,199601058,199601421,199602622,199603557,199603978,199604179,199609005,
# 0	ENSG00000118194	cmpl	cmpl	2,0,2,0,0,0,0,0,2,1,1,1,2,0,
#
# Coordinates are 0-based half-open
#

    regions = []

    f = open(filename,'r')
    for line in f:
        elts = line[:-1].split('\t')

        bin, name, chrom, strand, start, end, cstart, cend, exons, exonstart, exonend, id, name2 = elts[:13]
        start, end, cstart, cend, exons = map(int, (start,end,cstart,cend,exons))
        exonstart = map(int, exonstart.split(',')[:-1])
        exonend = map(int, exonend.split(',')[:-1])

        if strand == '+': exonidx = range(len(exonstart))
        else:             exonidx = range(len(exonstart)-1,-1,-1)

        # set initial position within transcript
        txpos = 1

        # loop over all exons, in the forward transcript direction
        for i in exonidx:

            # skip over non-coding exons
            if exonend[i] <= cstart: continue
            if exonstart[i] >= cend: continue

            # skip UTR for partly coding exons
            if exonstart[i] < cstart: exonstart[i] = cstart
            if exonend[i] > cend: exonend[i] = cend

            assert exonstart[i] < exonend[i]

            regions.append( (chrom, exonstart[i], exonend[i], strand, name, txpos) )
            
            txpos += exonend[i] - exonstart[i]

    return regions

                
def readampfile( filename ):
#
# Example entry:
# chr1 1000000 1000500 + amp_1
#
# Coordinates are 1-based, closed
#

    regions = []

    for line in open(filename,'r'):
        elts = line[:-1].split('\t')
        chrom, start, end, strand, name = elts
        start = int(start) - 1
        end = int(end)

        regions.append( (chrom, start, end, strand, name, 1) )

    return regions


# Takes a position in a regional coordinate system, and returns the native coordinates
# After mapping to nearest region using pos, adjust by adjustment (in regional coordinates)
def to_native( pos, name, regions, adjustment = 0 ):

    distances = []
    mindist = sys.maxint
    minidx = -1

    for i in range(len(regions)):
        if regions[i][4] != name:
            distance = sys.maxint
        else:
            # get start and end of region in regional coordinates (half-open)
            rstart = regions[i][5]
            rend = regions[i][5] + (regions[i][2] - regions[i][1])
            if rstart <= pos < rend:
                # within region: distance <= 0
                distance = -min( pos - rstart, rend - pos - 1 )
            elif pos < rstart:
                # not within region: distance > 0
                distance = rstart - pos
            else:
                # not within region: distance > 0
                distance = pos - rend + 1
        distances.append( distance )
        if distance < mindist:
            minidx = i
            mindist = distance

    if minidx == -1:
        sys.stderr.write("# Warning: locus %s:%s could not be converted (identifier not found)\n" % (name,pos))
        return None, 0

    # transform to 0-based coordinates, relative to start of region, and in regional direction
    newcoord = pos - regions[minidx][5] + adjustment

    # transform to native coordinates
    if regions[minidx][3] == "-":
        regionlen = regions[minidx][2] - regions[minidx][1]
        newcoord = regionlen - 1 - newcoord
    newcoord += regions[minidx][1]

    # return region (which includes chromosome), and native coordinate
    return regions[minidx], newcoord



# Takes a position in the native coordinate system, and calculates region-based coordinates
def to_region( pos, name, regions ):

    distances = []
    mindist = sys.maxint
    minidx = -1

    for i in range(len(regions)):
        if regions[i][0] != name:
            distance = sys.maxint
        elif regions[i][1] <= pos < regions[i][2]:
            distance = -min( pos - regions[i][1], regions[i][2] - 1 - pos )
        else:
            if pos < regions[i][1]: distance = regions[i][1] - pos
            else:                   distance = pos - regions[i][2] + 1
        distances.append( distance )
        if distance < mindist:
            minidx = i
            mindist = distance

    if minidx == -1:
        sys.stderr.write("# Warning: locus %s:%s could not be converted (identifier not found)\n" % (name,pos+1))
        return None, 0

    newcoord = pos - regions[minidx][1]

    if regions[minidx][3] == "-":
        newcoord = (regions[minidx][2] - regions[minidx][1]) - 1 - newcoord
    newcoord += regions[minidx][5]

    return regions[minidx], newcoord


# reverse-complementing variant data
def have_snp( data ):
    return len(data)>=4 and len(data[2])==1 and len(data[3])==1


def have_indel( data ):
    try:
        indel = int( data[2] )
    except:
        return data
    return len(data)>=4 and len(data[3]) == abs(indel)


def revcomp_data( data ):
    if have_snp( data ):
        data[2] = revcomp( data[2].upper() )
        data[3] = revcomp( data[3].upper() )
    elif have_indel( data ):
        data[3] = revomp( data[3].upper() )
    return data

    

def readsnps( handle, regions ):
#
# Parses a SNP file, and returns the native coordinates
#
# Example entry:
#
# chr/amplicon/gene  position[+-position]  oldnuc/indel_id  newnuc/more_indel_id  [even_more_indel_id]
#
    variants = []
    
    for line in handle:
        elts = line[:-1].split('\t')
        
        name = elts[0]
        
        posstring = elts[1]
        if posstring.count('+') == 1:
            pos, adjustment = map(int, posstring.split('+'))
        elif posstring.count('-') == 1:
            pos, adjustment = map(int, posstring.split('-'))
            adjustment = -adjustment
        else:
            pos = int(elts[1])
            adjustment = 0
            
        region,  nativepos = to_native( pos, name, regions, adjustment )

        # continue if no mapping found
        if region == None: continue

        chrom = region[0]

        newdata = [chrom,nativepos] + elts[2:]
        if region[3] == '-': revcomp_data( newdata )

        variants.append( newdata )

    return variants




#
# Writes out variants in region-based coordinates
#
def writesnps( snps, handle, regions, coding, chromdict ):

    for data in snps:

        chrom, pos = data[:2]
        region, newpos = to_region( pos, chrom, regions )
        newdata = data[:]

        within = region[5] <= newpos < region[5] + region[2]-region[1]

        if not within:
            if newpos < region[5]: adjustment = newpos - region[5]
            else:                  adjustment = newpos - (region[5] + region[2]-region[1] - 1)
            # obtain new basal coordinates just within the region
            if region[3] == "+": adjpos = pos - adjustment
            else:                adjpos = pos + adjustment
            region, newpos = to_region( adjpos, chrom, regions )
        else:
            adjustment = 0
        
        havesnps = have_snp( data )

        if region[3] == '-': revcomp_data( newdata )

        if coding and havesnps and within:

            frame = (newpos - 1) % 3          # this assumes that tx are 1-based
            region0, codonnativepos0 = to_native( newpos - frame, region[4], regions )
            region1, codonnativepos1 = to_native( newpos - frame+1, region[4], regions )
            region2, codonnativepos2 = to_native( newpos - frame+2, region[4], regions )
            assert region0[0] == chrom
            assert region1[0] == chrom
            assert region2[0] == chrom
            # obtain codon on forward strand
            codonseq = ( chromdict[chrom][ min(codonnativepos0,codonnativepos2) ] +
                         chromdict[chrom][ codonnativepos1 ] +
                         chromdict[chrom][ max(codonnativepos0,codonnativepos2) ] )
            # change to exon strand; build alternative codon
            if region[3] == '-':
                codonseq = revcomp( codonseq )
                if codonseq[frame] != revcomp( data[2] ):
                    sys.stderr.write("# Warning: reference for SNP %s:%s has %s (in reverse direction), but SNP file says %s" % (chrom,pos+1,codonseq[frame],revcomp(data[2])))
                altcodonseq = list(codonseq)
                altcodonseq[frame] = revcomp( data[3] )
            else:
                if codonseq[frame] != data[2]:
                    sys.stderr.write("# Warning: reference for SNP %s:%s has %s, but SNP file says %s" % (chrom,pos+1,codonseq[frame],data[2]))
                
                altcodonseq = list(codonseq)
                altcodonseq[frame] = data[3]
            altcodonseq = ''.join(altcodonseq)
            a2a = "%s%s>%s" % ((newpos+2) // 3, codon2aminoacid(codonseq), codon2aminoacid(altcodonseq))
            newdata.append( a2a )

        elif havesnps:

            snp = chromdict[ chrom ][ pos : pos+1 ]
            if snp.upper() != data[2].upper():
                sys.stderr.write("# Warning: reference for SNP %s:%s has %s, but SNP file says %s\n" % (chrom,pos+1,snp,data[2]))

        newdata[0] = region[4]
        newdata[1] = str(newpos)
        if adjustment > 0:
            newdata[1] += "+" + str(adjustment)
        elif adjustment < 0:
            newdata[1] += str(adjustment)

        handle.write( "\t".join(newdata) + "\n" )

        



def parseopts( argv ):

    help = """
Converts SNP/indel calls on stdin to another coordinate system

Usage:  convertcoords.py [options]

Options:

  --from genomic|amplicon|cdna   Set coordinate system to convert from
  --to   genomic|amplicon|cdna   Set coordinate system to convert to

  --ref directory                Directory with .fa[.gz] chromosome files; names correspond to chromosomes
  --amp file                     File with amplicon definitions (chromosome/start/end/strand/name)
  --tx file                      File with transcript definitions (ucsc ensGene / refGene format)


File formats:

  This utility uses 2 types of files: one to define a coordinate system, and one to define variants (SNPs or
  indels).  For each type, there are specific file formats for each of the 3 coordinate systems supported:
  genomic, amplicon-based, or cdna (transcript) based.


  1. Coordinate systems:

  A genomic coordinate system is defined by providing a directory with chromosome .fa files for the reference.
  A reference genomic coordinate system is required in all cases.  The chromosome names are taken to be the
  first word on the fasta comment (">") line.  Each .fa file must contain exactly one chromosome.

  An amplicon coordinate system is defined by providing the chromosome, start and end coordinates, read
  directions and names of a set of amplicons, in a tab-delimited format.  Coordinates here are relative to
  the reference, and are 1-based closed-segment.  For example,

      chr1 1000000 1000500 + amplicon_1
      chr1 1001000 1001500 - amplicon_2

  A transcript coordinate system is defined by providing gene models, in the ucsc ensGene.txt / refGene.txt
  format.  Note that the coordinates in these files use the 0-based open-segment convention.


  2. Variants

  The file formats for SNPs in the various coordinate systems are similar, and follow this pattern:

      <name> <position> <reference> <alternative>

  where fields are tab-delimited.  Here, 'name' refers to the chromosome, amplicon, or transcript.  Position
  is 1-based in all cases.  When a locus falls outside all of the regions, the displacement relative to the
  start or end of the nearest region is given in the usual notation, e.g.

      amplicon4 1-10 A C
      amplicon4 123+1 A G

  means that the A-to-C variant (in the direction of the amplicon) is located 10 bases before the start of
  the amplicon, and the A-to-G variant is located 1 bp beyond the end of the amplicon (which is 123 bp long).

  On output, when translation is to a transcript, the reference and variant codon, and their position within
  the amino acid sequence, is provided on column 5.

  Indel variants are represented by a number signifying the (type and) length of the indel in column 3, and
  the inserted or deleted sequence in column 4.  When the absolute value of column 3 equals the length of the
  sequence in column 4, the data is assumed to represent an indel.  In this case, the indel sequence is
  reverse-complemented when necessary.

  Both for SNPs and indels, additional data may be present in columns 5 and onwards, and will be copied without
  change.

"""

    coorddict = {'genomic':[None],'amplicon':[None],'cdna':[None]}

    fromcoord = None
    tocoord = None
    chromdict = None

    if len(argv) == 1:
        print help
        sys.exit(1)

    i = 1
    while i < len(argv):

        if argv[i].startswith('--help') or argv[i].startswith('-?'):
            print help
            sys.exit(0)

        if not argv[i].startswith('--') or i + 1 >= len(argv):
            print "Error parsing command line: '%s'" % argv[i]
            sys.exit(1)

        if argv[i].startswith('--from'):
            fromcoord = coorddict.get( argv[i+1], None )
        elif argv[i].startswith('--to'):
            tocoord = coorddict.get( argv[i+1], None )
        elif argv[i].startswith('--ref'):
            coorddict['genomic'][0], chromdict = getchromosomes( argv[i+1] )
        elif argv[i].startswith('--amp'):
            coorddict['amplicon'][0] = readampfile( argv[i+1] )
        elif argv[i].startswith('--tx'):
            coorddict['cdna'][0] = readtxfile( argv[i+1] )
        else:
            raise ValueError("Don't understand option %s" % argv[i])

        i += 2

    if fromcoord == None: raise ValueError("Please specify system to convert from")
    if fromcoord[0] == None: raise ValueError("Please provide data for original coordinate system")
    if tocoord == None: raise ValueError("Please specify system to convert to")
    if tocoord[0] == None: raise ValueError("Please provide data for target coordinate system")

    return fromcoord[0], tocoord[0], coorddict, chromdict





if __name__ == "__main__":

    fromcoord, tocoord, coorddict, chromdict = parseopts( sys.argv )
    snps = readsnps( sys.stdin, fromcoord )
    writesnps( snps, sys.stdout, tocoord, tocoord == coorddict['cdna'][0], chromdict )
    
        
