"""

Converts .bam file to .fastq format

Gerton Lunter, Jan 2010

"""

import optparse, subprocess, os, gzip

FLAG_PAIRED=1
FLAG_REVERSE = 16
FLAG_FIRST=64
FLAG_SECOND=128

parser = optparse.OptionParser(usage="usage: %prog [options] INPUT.BAM")
parser.add_option("-m", "--maxmem", dest="maxmem", type="int", help="(approximate) maximum amount of memory to use")
parser.add_option("-o", "--output", dest="outfile", help="write output to FILE[.|_1.|_2.]fastq.gz (default INPUT)", metavar="FILE")
parser.add_option("-r", "--readgroup", dest="readgroup", help="only process reads from READGROUP", metavar="READGROUP")
parser.add_option("-t", "--tempdir", dest="tempdir", help="use DIR for temporary file (default %default)", metavar="DIR", default="./")
parser.add_option("-F", "--force", dest="force", action="store_true", help="ignore errors during conversion - use at own risk")
parser.add_option("-v", "--verbose", dest="verbose", action="store_true" )

options, args = parser.parse_args()
if len(args) != 1: parser.parse_args(["-h"])

if options.readgroup:
    tempfileprefix = options.tempdir + os.path.split(args[0])[1] + "." + options.readgroup + ".sorted"
else:
    tempfileprefix = options.tempdir + os.path.split(args[0])[1] + ".sorted"

try:
    if options.verbose: print "Sorting",args[0]
    sortcommand = "samtools sort -n "
    if options.maxmem: sortcommand += "-m " + str(options.maxmem) + " "
    if not options.readgroup:
        sortcommand += args[0] + " " + tempfileprefix
        retcode = subprocess.call(args=sortcommand.split())
    else:
        p1 = subprocess.Popen(["samtools","view","-u","-r",options.readgroup,args[0]],stdout=subprocess.PIPE)
        sortcommand += "- " + tempfileprefix
        p2 = subprocess.Popen(sortcommand.split(),stdin=p1.stdout)
        retcode = os.waitpid(p2.pid,0)[1]
except OSError:
    raise ValueError("Cannot execute samtools")
# unfortunately samtools does not report errors by return code
if retcode>0: raise ValueError("Error executing '%s'" % sortcommand)


def revcomp( seq ):
    d = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
    try:
        s = [ d[c] for c in seq.upper() ]
    except:
        raise ValueError("Unknown character encountered in sequence '%s'" % seq)
    s.reverse()
    return ''.join(s)


def reverse( quals ):
    q = [ c for c in quals ]
    q.reverse()
    return ''.join(q)


def reads( inputfile ):
    global errors
    try:
        viewcommand = "samtools view " + tempfileprefix + ".bam"
        viewer = subprocess.Popen(args=viewcommand.split(), bufsize=-1, stdout=subprocess.PIPE)
        reads = 0
        for line in viewer.stdout:
            if line.startswith('@'): continue
            cols = line[:-1].split('\t')
            if len(cols) < 11: 
                errors += 1
                if options.force: continue
                raise ValueError("Cannot parse output from 'samtools view': %s" % line[:-1])
            label, flag, seq, qual = cols[0], cols[1], cols[9], cols[10]
            for opts in cols[11:]:
                if opts.startswith('OQ:Z:'): qual = opts[5:]
            try: intflag = int(flag)
            except: 
                errors += 1
                if options.force: continue
                raise ValueError("Cannot parse output from 'samtools view': %s" % line[:-1])
            if len(qual) != len(seq): raise ValueError("Sequence and quality not the same length: %s" % line[:-1])
            if intflag & FLAG_PAIRED:
                if intflag & FLAG_FIRST: pairflag=1
                elif intflag & FLAG_SECOND: pairflag=2
                else: 
                    errors += 1
                    if options.force: continue
                    else: raise ValueError("Inconsistent flag encountered (paired, but not 1st or 2nd): %s" % line[:-1])
            else:
                pairflag = 0
            if intflag & FLAG_REVERSE:
                seq = revcomp(seq)
                qual = reverse(qual)
            reads += 1
            yield (label, pairflag, seq, qual)
    except ValueError:
        raise
    except:
        raise ValueError("Error executing '%s'" % viewcommand)
    if viewer.returncode > 0:
        raise ValueError("Command '%s' returned with return code %s" % (viewcommand,returncode))


def write(fnum, label, seq, qual):
    if not outputfiles[fnum]:
        try:
            outputfiles[fnum] = gzip.GzipFile(outputfilenames[fnum],'w')
            if options.verbose: print "Writing to %s" % outputfilenames[fnum]
        except:
            raise ValueError("Could not open '%s' for writing" % outputfilenames[fnum])
    counts[fnum] += 1
    outputfiles[fnum].write("@%s\n%s\n+%s\n%s\n" % (label,seq,'',qual))


def sameread(first,second):
    if first==second: return True
    if first.endswith('/1') and second.endswith('/2') and first[:-2] == second[:-2]: return True
    return False


def addsuffix(label,suffix):
    if label.endswith("/1") or label.endswith("/2"): label = label[:-2]
    return label + suffix


outprefix = options.outfile
if not outprefix: outprefix, ext = os.path.splitext(args[0])
outputfilenames = [outprefix+".fastq.gz", outprefix+"_1.fastq.gz", outprefix+"_2.fastq.gz"]
outputfiles, counts = [None,None,None], [0,0,0]
first, second, unpaired = None, None, 0
errors = 0

for label, pairflag, seq, qual in reads( tempfileprefix + ".bam" ):
    if pairflag == 0: write( 0, label, seq, qual )
    elif pairflag == 1:
        if first:
            write( 0, first[0], first[1], first[2])
            unpaired += 1
        first = label, seq, qual
    elif pairflag == 2:
        if second:
            write( 0, second[0], second[1], second[2])
            unpaired += 1
        second = label, seq, qual
    if first and second:
        if sameread(first[0],second[0]):
            write( 1, addsuffix(first[0],"/1"), first[1], first[2] )
            write( 2, addsuffix(second[0],"/2"), second[1], second[2] )
        else:
            write( 0, first[0],first[1],first[2] )
            write( 0, second[0],second[1],second[2] )
        first, second = None, None

if first:
    write( 1, first[0], first[1], first[2] )
    unpaired += 1
if second:
    write( 2, second[0], second[1], second[2] )
    unpaired += 1

for idx,f in enumerate(outputfiles):
    try:
        if f: f.close()
    except:
        raise ValueError("Problem closing '%s'" % outputfilenames[idx])

os.unlink( tempfileprefix + ".bam" )

if unpaired > 0: print "Warning: found %s mateless pairs" % unpaired
if errors > 0: print "Ignored %s reads with errors" % errors
if options.verbose: 
    print "Written %s single reads and %s mate pairs" % (counts[0],counts[1])
