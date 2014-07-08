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
    print "Marks Mendelian errors"
    print "Usage: vcf-mendelian.py [options]"
    print
    print " -h, --help        This help"
    print " -oF, --output F   Send output to file F (default stdout)"
    print " -iF, --input F    Read input from file F (default stdin)"
    print " -3                Set default input version to 3.3"
    print " -F                Use fast and sloppy genotype parsing"
    print " -cF, --child F    Set child label to F"
    print " -qN               Set minimum genotype quality to N (default 0)"
    print " -xE, --ignore E   Ignore error/warning E"
    print
    print "Possible errors are:"
    v = vcf.VCF()
    print wrap( ", ".join([v.split(':')[0] for v in v._errors.values()]) )

def main():
    infile = sys.stdin
    outfile = sys.stdout
    inversion = 40
    child = None
    mingq = 0
    fastGT = False
    v = vcf.VCF( leftalign=True, _fastGT=True )
    try:
        opts, args = getopt.getopt(sys.argv[1:], "ho:i:x:3Fc:q:", ["help","output","input","ignore","child"])
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
        elif o in ["-c","--child"]:
            child = a
        elif o in ["-q"]:
            mingq = int(a)
        elif o in ["-x","--ignore"]:
            v.ignoreerror(a)

    # check that we have a child
    if child == None: raise ValueError("Need to set child label")
    
    # process data
    v.setversion(inversion)
    v._fastGT = fastGT
    vcfstream = v.parse( infile )

    # instantiate vcfout from v, to include header and definitions.
    vcfout = vcf.VCF(v)  
    vcfout.setversion(40)
    vcfout.getinfo()['MENDELERROR'] = vcf.FORMAT('MENDELERROR',vcfout.NT_NUMBER,0,'Flag','Non-Mendelian segregation','')

    vcfout.writeheader( outfile )

    for data in vcfstream:
        samples = v.getsamples()
        if child not in samples:
            raise ValueError("Label for child (%s) not found among labels in file (%s)" % (child, ",".join(samples)))
        if len(samples) != 3:
            raise ValueError("Expect exactly 3 samples")
        parent_gts, child_gt = [], None
        gqs = []
        for s in samples:
            sampledata = data[s]
            genotype = sampledata.get('GT',['.'])[0]
            gq = sampledata.get('GQ',[0])[0]
            if gq < mingq: continue  # low genotype quality
            if len(genotype) == 1: continue    # haploid, or missing data
            if genotype[0] == "." or genotype[2] == ".": continue  # missing data
            if s == child:
                child_gt = genotype
            else:
                parent_gts.append(genotype)

        # check for mendel error
        if child_gt != None and len(parent_gts) == 2:
            if ( (child_gt[0] in [parent_gts[0][0], parent_gts[0][2]] and child_gt[2] in [parent_gts[1][0], parent_gts[1][2]]) or
                 (child_gt[2] in [parent_gts[0][0], parent_gts[0][2]] and child_gt[0] in [parent_gts[1][0], parent_gts[1][2]]) ):
                # no error
                pass
            else:
                # check genotype qualities
                
                data['info']['MENDELERROR'] = []

        vcfout.write_data( outfile, data )


if __name__ == "__main__": main()
