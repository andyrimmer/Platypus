#
# naively call variants from BAM
#

import vcf
import filez
import pysam
import sys
import getopt
import re

def wrap(line, length=80):
    if len(line)<=length: return line
    if len(line)==0: return ""
    while length>=0 and line[length] != " ": length -= 1
    start = length
    while length>=0 and line[length] == " ": length -= 1
    if length>=0: return line[:length+1] + "\n" + wrap(line[start+1:])
    while length<len(line) and line[length] != " ": length += 1
    return line[:length] + "\n" + wrap(line[length+1:])


def help():
    print "Naively call variants from BAM"
    print
    print " -h, --help         This help"
    print " -oF, --output F    Send output to file F (default stdout)"
    print " -iF, --input F     Read input from file F (default stdin)"
    print " -rF, --reference F Use reference fasta file; must be samtools faidx'ed (required)"
    print " --sample=S         Use S as sample label"
    print " -sN                Report variants if <N SNPs per 1000 bp"
    print " -nN                Report variants if <N indels per 1000 bp"




def main():
    infile = sys.stdin
    outfile = sys.stdout
    reference = None
    sample = "sample"
    maxsnps = 10
    maxindels = 5
    maxwindow = 1000
    v = vcf.VCF( _fastGT = False )
    try:
        opts, args = getopt.getopt(sys.argv[1:], "ho:i:r:s:n:", ["help","output=","input=","reference=","sample="])
    except:
        help()
        raise
    for o, a in opts:
        if o in ["-h","--help"]:
            help()
            sys.exit()
        elif o in ["-o","--output"]:
            outfile = open(a,'w')
        elif o in ["-i","--input"]:
            infile = filez.open(a,'r')
        elif o in ["-r","--reference"]:
            reference = a
        elif o in ["-s"]:
            maxsnps = int(a)
        elif o in ["-n"]:
            maxindels = int(a)
        elif o in ["--sample"]:
            sample = a
            
    if not reference:
        print "Reference required"
        help()
        sys.exit()

    # open reference
    fa = pysam.Fastafile( reference )
    v.setreference(fa)
    
    # instantiate vcfout from v, to include header and definitions.
    vcfout = vcf.VCF(v)
    vcfout.setversion(41)

    vcfout.getheader().append( ("source",' '.join(sys.argv)))

    vcfout.getinfo()['RN'] = vcf.FORMAT('RN', vcfout.NT_NUMBER, 1, "String", "Name of read supporting the variant", -1)
    vcfout.getfilter()['HighCallDensity'] = vcf.FORMAT('mask',vcfout.NT_NUMBER,0,"Flag","Number of indels/snps exceeding %s/%s in %s bp window" % (maxindels,maxsnps,maxwindow),".")

    vcfout.setsamples( [sample] )

    vcfout.writeheader( outfile )

    for data in sys.stdin:
        
        if data.startswith("#"):
            continue

        readname, flag, chrom, pos, mapq, cigar, mchrom, mpos, isize, seq, qual = data[:-1].split('\t')[:11]

        vcflist = []
        vcfdata = {'chrom':chrom, 'id':'.', 'qual':100, 'filter': ["HighCallDensity"], 'info': {'RN': [readname]}, 'format': ['GT'], sample:{'GT':["1"]}}

        rpos = int(pos)-1
        spos = 0

        for i in re.finditer( "[0-9]*[MIDSNHP]", cigar ):
            ash = i.group(0)
            num = int( ash[:-1] )
            if ash[-1] == 'P' or ash[-1] == 'H':
                pass
            elif ash[-1] == 'M':
                refseq = fa.fetch(chrom,rpos,rpos+num)
                for idx,c in enumerate(refseq):
                    c = c.upper()
                    if c != 'N' and seq[spos+idx] != 'N' and c != seq[spos+idx]:
                        # found a SNP
                        vcfdata['pos'] = rpos + idx
                        vcfdata['ref'] = c
                        vcfdata['alt'] = [seq[spos+idx]]
                        vcflist.append( vcfdata.copy() )
                spos += num
                rpos += num
            elif ash[-1] == 'S':
                spos += num
            elif ash[-1] == 'I':
                vcfdata['pos'] = rpos - 1
                vcfdata['ref'] = fa.fetch(chrom, rpos-1, rpos)
                vcfdata['alt'] = [fa.fetch(chrom, rpos-1, rpos) + seq[spos:spos+num]]
                vcflist.append( vcfdata.copy() )
                spos += num
            elif ash[-1] == 'D':
                vcfdata['pos'] = rpos - 1
                vcfdata['ref'] = fa.fetch(chrom, rpos-1, rpos + num)
                vcfdata['alt'] = [fa.fetch(chrom, rpos-1, rpos)]
                vcflist.append( vcfdata.copy() )
                rpos += num

        # see if snp or indel density is too high
        windowidx = 0
        for idx,call in enumerate(vcflist):
            while vcflist[windowidx]['pos'] < call['pos'] - maxwindow:
                windowidx += 1
            indels = len( [ c for c in vcflist[windowidx:idx+1] if len(c['ref']+c['alt'][0])>2 ] )
            snps = len( [ c for c in vcflist[windowidx:idx+1] if len(c['ref']+c['alt'][0])==2 ] )
            if indels > maxindels or snps > maxsnps: 
                break
        else:
            for call in vcflist:
                call['filter'] = None

        for call in vcflist:
            vcfout.write_data( outfile, call )

if __name__ == "__main__": main()


