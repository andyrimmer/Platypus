#
# add allele count field
#

import vcf
import filez
import sys
import getopt


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
    print "Add allele counts (AC field) to INFO column"
    print "Usage: vcf-add-ac.py [options]"
    print
    print " -h, --help        This help"
    print " -oF, --output F   Send output to file F (default stdout)"
    print " -iF, --input F    Read input from file F (default stdin)"
    print " -3                Set default input version to 3.3"
    print " -F                Use fast and sloppy genotype parsing"
    print " -xE, --ignore E   Ignore error/warning E"
    print
    print "Possible errors are:"
    v = vcf.VCF()
    print wrap( ", ".join([v.split(':')[0] for v in v._errors.values()]) )

def main():
    infile = sys.stdin
    outfile = sys.stdout
    inversion = 40
    fastGT = False
    v = vcf.VCF( leftalign=True, _fastGT=True )
    try:
        opts, args = getopt.getopt(sys.argv[1:], "ho:i:x:3F", ["help","output","input","ignore"])
    except:
        help()
    for o, a in opts:
        if o in ["-h","--help"]:
            help()
            sys.exit()
        elif o in ["-o","--output"]:
            outfile = open(a,'w')
        elif o in ["-i","--input"]:
            infile = filez.open(a,'r')
        elif o == "-3":
            inversion = 33
        elif o == "-F":
            fastGT = True
        elif o in ["-x","--ignore"]:
            v.ignoreerror(a)
    
    # process data
    v.setversion(inversion)
    v._fastGT = fastGT
    vcfstream = v.parse( infile )

    # instantiate vcfout from v, to include header and definitions.
    vcfout = vcf.VCF(v)  
    vcfout.setversion(40)
    vcfout.getinfo()['AC'] = vcf.FORMAT('AC',vcfout.NT_ALLELES,-1,"Integer","Allele counts",-1)

    vcfout.writeheader( outfile )

    for data in vcfstream:
        samples = v.getsamples()
        ac = {}
        for s in samples:
            sampledata = data[s]
            genotype = sampledata.get('GT',['.'])[0]
            if len(genotype) == 1: continue    # haploid, or missing data
            for idx in 0,2:
                if genotype[idx] == ".": continue
                ac[genotype[idx]] = ac.get(genotype[idx],0) + 1
            if len(ac.keys()) > 0:
                data['info']['AC'] = [ ac.get(i,0) for i in range(len(data['alt'])+1) ]
        vcfout.write_data( outfile, data )


if __name__ == "__main__": main()
