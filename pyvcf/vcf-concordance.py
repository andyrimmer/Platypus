#
# annotate calls that are present in a secondary VCF file
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
    print "Annotate calls by their presence in a secondary VCF"
    print "Usage: vcf-concordance.py [options]"
    print
    print " -h, --help           This help"
    print " -oF, --output F      Send output to file F (default stdout)"
    print " -iF, --input F       Read input from file F (default stdin)"
    print " -cF, --concordance F Read secondary calls from file F"
    print " -lL, --label L       Set INFO label to L"
    print " -dD, --description D Set description field to D"
    print " -3                   Set default input version to 3.3"
    print " -xE, --ignore E      Ignore error/warning E"
    print
    print "Possible errors are:"
    v = vcf.VCF()
    print wrap( ", ".join([v.split(':')[0] for v in v._errors.values()]) )

def main():
    infile = sys.stdin
    secondaryin = None
    outfile = sys.stdout
    label = None
    description = None
    inversion = 40
    adjust = 0      # hack, to cope with Quang's off-by-one VCFs
    v = vcf.VCF()
    try:
        opts, args = getopt.getopt(sys.argv[1:], "ho:i:x:c:l:d:3", ["help","output","input","ignore","concordance","label","description","quang"])
    except getopt.GetOptError, err:
        print str(err)
        help()
    for o, a in opts:
        if o in ["-h","--help"]:
            help()
            sys.exit()
        elif o in ["-o","--output"]:            outfile = open(a,'w')
        elif o in ["-i","--input"]:             infile = filez.open(a,'r')
        elif o in ["-c","--concordance"]:       secondaryin = filez.open(a,'r')
        elif o in ["-l","--label"]:             label = a
        elif o in ["-d","--description"]:       description = a
        elif o == "-3":                         inversion = 33
        elif o in ["-x","--ignore"]:            v.ignoreerror(a)
        elif o == "--quang":                    adjust = -1

    if not description or not label or not secondaryin:
        raise ValueError("Need concordance file; label; and description")
    
    # process data
    v.setversion(inversion)
    secondary = vcf.VCF(v)
    vcfstream = v.parse( infile )

    # read secondary file, and store calls in a dictionary
    secondarystream = secondary.parse( secondaryin )
    secondarydata = {}
    for data in secondarystream:
        if data['chrom'] not in secondarydata: secondarydata[data['chrom']] = {}
        d = secondarydata[data['chrom']]
        for a in data['alt']:
            if data['pos'] not in d: d[data['pos']] = data['ref']+"|"+a
            else:                    d[data['pos']] = d[data['pos']] + "|" + data['ref']+"|"+a

    # instantiate vcfout from v, to include header and definitions.
    vcfout = vcf.VCF(v)  
    vcfout.setversion(40)
    vcfout.getinfo()[label] = vcfout.FORMAT(label,vcfout.NT_NUMBER,0,"Flag",description,".")
    vcfout.writeheader( outfile )

    for data in vcfstream:
        # mark position if ANY of the alleles exist in dbSNP
        p, r = data['pos'] + adjust, data['ref']
        data['pos'] = p
        call = False
        secondary = secondarydata.get(data['chrom'],{}).get(p,None)
        if secondary:
            calls = secondary.split('|')
            for a in data['alt']:
                for idx in range(0,len(calls),2):
                    if vcfout.compare_calls(p,r,a,p,calls[idx],calls[idx+1]):
                        call = True
                        break
        if call:
            data['info'][label] = []
        vcfout.write_data( outfile, data )


if __name__ == "__main__": main()
