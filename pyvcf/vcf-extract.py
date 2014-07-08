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
    print "Extract information from a VCF file in a column-oriented format"
    print "Usage: vcf-add-hr.py [options]"
    print
    print " -h, --help         This help"
    print " -oF, --output=F    Send output to file F (default stdout)"
    print " -iF, --input=F     Read input from file F (default stdin)"
    print " -eS, --extract=S   Specify information to extract (see below)"
    print " --keepheader       Copy the VCF header to output"
    print " -3                 Set default input version to 3.3"
    print " -parsegt           Disable fast genotype parsing"
    print " -xE, --ignore=E    Ignore error/warning E"
    print
    print "The specification string S is a comma-separated string of field specifiers, or *"
    print "for all fields.  Field specifiers can be any column heading, possibly followed by"
    print "a subfield and an index, separated by colons.  For example, the string"
    print 
    print "    CHROM,POS,INFO:DP,NA12878:GT,NA12878:HQ:1,NA12878:HQ:2"
    print
    print "would specify the chromosome and position columns, the depth of coverage subfield"
    print "from the INFO column, and the two haplotype quality fields for the sample NA12878"
    print "in separate columns.  The three elements of the genotype field (first haplotype,"
    print "haplotype/genotype indicator, and second haplotype) can also be extracted"
    print
    print "Possible errors are:"
    v = vcf.VCF()
    print wrap( ", ".join([v.split(':')[0] for v in v._errors.values()]) )


def parse_extract(e):
    v = vcf.VCF()
    cols = []
    for col in e.split(','):
        elts = col.split(':')
        # change standard columns to lower case
        if elts[0] in v._required: elts[0] = elts[0].lower()
        # parse
        if len(elts) == 1:
            cols.append( (col, elts[0], None, None) )
        elif len(elts) == 2:
            cols.append( (col, elts[0], elts[1], None) )
        elif len(elts) == 3:
            try:
                idx = int(elts[2])
            except ValueError:
                raise ValueError("Problem parsing '%s': expected integer in 3rd position" % col)
            cols.append( (col, elts[0], elts[1], int(elts[2])-1) )
        else:
            raise ValueError("Problem parsing '%s': expected 1, 2 or 3 elements" % col)
    return cols


def allfields(v):
    cols = [ (i, i.lower(), None, None) for i in v._required[:6] + v._required[8:9] ]     # required fields, except FILTER, INFO
    # add FILTER fields
    for i in v.getfilter().itervalues(): cols.append( ("FILTER:"+i.id, "filter", i.id, None) )
    # add INFO fields
    for i in v.getinfo().itervalues():
        if   i.numbertype == v.NT_NUMBER: n = i.number
        elif i.numbertype == v.NT_ALLELES: n = 2         # this should really be data-dependent
        elif i.numbertype == v.NT_GENOTYPES: n = 3       # this should really be data-dependent
        else:
            print "found numbertype",i.numbertype,i
        if n <= 1: cols.append( ("INFO:"+i.id, "info", i.id, None) )
        else:
            cols += [ ("INFO:"+i.id+":"+str(j+1), "info", i.id, j) for j in range(n) ]
    # add sample fields
    cols += [ (s, s, None, None) for s in v.getsamples() ]
    return cols


def main():
    infile = sys.stdin
    outfile = sys.stdout
    inversion = 40
    extract = None
    keepheader = False
    v = vcf.VCF( _fastGT = True )
    try:
        opts, args = getopt.getopt(sys.argv[1:], "ho:i:x:e:3", ["help","output=","input=","ignore=","extract=","keepheader","parsegt"])
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
        elif o == "-3":
            inversion = 33
        elif o in ["-x","--ignore"]:
            print a
            v.ignoreerror(a)
        elif o in ["-e","--extract"]:
            extract = a
        elif o in ["--keepheader"]:
            keepheader = True
        elif o in ["--parsegt"]:
            v._fastGT = False
    if not extract:
        print "Specification string required"
        help()
        sys.exit()

    if extract != "*": columns = parse_extract(extract)

    # process data
    v.setversion(inversion)
    vcfstream = v.parse( infile )

    # copy vcf header
    if keepheader: v.writeheader(outfile)

    # deal with 'all fields'
    if extract == "*": columns = allfields(v)

    # write header line
    outfile.write("\t".join( c[0] for c in columns ) + "\n")

    for data in vcfstream:
        cols = []
        for e in columns:
            if e[1] not in data:
                raise ValueError("No column found for '%s'" % e[0])
            col = data[e[1]]
            if e[1] == "pos": col += 1  # use 1-based coordinates
            if e[2] != None:
                if e[1] == "filter":
                    # see if key exists in list
                    if e[2] in col: col = [e[2]]
                    else: col = ["."]
                else:
                    # extract from dictionary
                    if type(col) != type({}): 
                        print col
                        raise ValueError("Cannot extract '%s' from column '%s'" % (e[2],e[1]))
                    col = col.get(e[2],["."])
                # Allow extraction from the single GT element
                if e[2] == "GT": col = col[0]
                if e[3] != None:
                    if len(col) > e[3]:
                        col = [col[e[3]]]
                    else:
                        col = ["."]
            # format the various types
            if type(col) in [type(""),type(0),type(0.0)]:
                col = str(col)
            elif type(col) == type({}):
                if e[0] == "INFO":
                    col = v.format_formatdata( col, v._info, separator=";" )
                else:
                    col = v.format_formatdata( col, v._format, key=False )
            elif type(col) == type([]):
                # change unextracted GT back into text format; other unextracted fields are comma-separated
                if e[2] == "GT" and len(col)==3:     col = ''.join(map(str,col))  
                elif e[1] == "ref" or e[1] == "alt": col = ','.join(map(str,col))
                else:                                col = ":".join(map(str,col))
            else:
                print "Unexpected type found: ",type(col)
                print "Value:",col
                raise ValueError("")

            cols.append(col)

        print "\t".join(cols)

if __name__ == "__main__": main()
