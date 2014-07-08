# 
# Diff two vcf files with various options for diff options and output
# 
# Iain Mathieson, 28 July 2010
# 
# TODO
# 1) Make the outputbuffer collect stats and possibly point to different files etc...
# 2) Option to relax the alt comparison
# 3) Collect better stats, print them more nicely

import vcf
import filez
import sys
import getopt

from collections import deque, defaultdict
from heapq import heappush, heappop

# Match codes - 0 means no match, and higher codes are other sorts of match - in increasing order of preference
# Should be able to change these codes arbitrarily. 


MATCH_MATCH  = 2 # Lines match apart from position
MATCH_GTYPE  = 1 # Lines match apart from genotype and position
MATCH_NO     = 0 # Lines do not match

def help():
    print "Diff two vcf files"
    print "By default compares position and variant exactly and ignores genotype"
    print ""
    print "Usage: python2.6 vcf-diff.py file1.vcf file2.vcf [ -options ]"
    print ""
    print "Output:"
    print "+[line] Line added in file2 relative to file1"
    print "-[line] Line removed in file2 relative to file1"
    print "*[line] Line matched but genotypes are different"
    print ""
    print "Options ( --command as [command] ):"
    print "-h              Print this [help] text"
    print "-s              Just print a [sumamry] of the output"
    print "-o outfile.txt  Write [output] to outfile.txt"
    print "-r n            [relax] position - match if within n bases"
#TODO:     print "-a n            Relax [alt]ernative - match if edit distance < n"

def parse_options():
    """
    First two arguments are file1 and file2, which must be readable vcf files.
    Other options are described by the help() function
    """
    options ={
        "outfile":   sys.stdout,
        "summary":   False,
        "relax_pos": 0,
        "relax_alt": 0,
        "ignore_genotypes": False
        }

   
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hso:r:ga:", ["help","summary","output","relax","ignore_genotypes","alt"])
        options["file1"] = open(args[0],'r')
        options["file2"] = open(args[1],'r')
    except Exception as err:
        print str(err)
        help()


    for o, a in opts:
        if o in ["-h","--help"]:
            help()
            sys.exit()
        elif o in ["-o","--output"]:            options["outfile"] = open(a,'w')
        elif o in ["-s","--summary"]:           options["summary"] = True
        elif o in ["-r","--relax"]:             options["relax_pos"] = int(a)        
        elif o in ["-g","--ignore_genotypes"]:  options["ignore_genotypes"] = True                
#         elif o in ["-a","--alt"]:               options["relax_alt"] = int(a)

    return options

class match_engine:
    """
    An object that can work out how well two lines from two vcf files match
    """

    def __init__(self, vcf1, vcf2, options):
        self.ignore_genotypes = options["ignore_genotypes"]
        self.samples = vcf1.getsamples()
        if self.samples != vcf2.getsamples(): 
            print "Samples are different, ignoring genotypes"
            self.ignore_genotypes = True

    def match(self, line1, line2):
        alt_match = (line1["alt"] == line2["alt"])
        gt_match = self.ignore_genotypes or all([ line1[s] == line2[s] for s in self.samples ] )

        if alt_match and gt_match: return MATCH_MATCH
        elif alt_match: return MATCH_GTYPE
        else: return MATCH_NO

class outputbuffer:
    """
    Buffer the output - keeps output in a heap so that we always print in the right order
    """

    def __init__(self,vcf,options):
        self.max_size = 100
        self.buffer = []
        self.stats = defaultdict(int)
        self.vcf = vcf
        if options["summary"]:
            self.out = open( "/dev/null", "w" )
        else:
            self.out = options["outfile"]

    def add(self, element ):
        """
        element is a tuple of ( position, prefix, line )
        """
        heappush(self.buffer, element)
        if len(self.buffer) > self.max_size: self.lineout()

    def lineout(self):
        if len(self.buffer) > 0:
            element = heappop(self.buffer)
            self.out.write(element[1] + " ")
            self.vcf.write_data(self.out,element[2])
    
    def flush(self):
        while len(self.buffer) > 0:
            self.lineout()
    

def main():
    """
    We iterate over each line of file1. For each line, we build a list of all the lines in file2 which are within a 
    given distance of that line. Then we check if any of them match, and work out which lines to output. The lines 
    are output to the output buffer which attempts to output them in the correct order. 
    """
    options = parse_options()
    vcf1 = vcf.VCF()
    vcf2 = vcf.VCF()
    line_gen1 = vcf1.parse(options["file1"])
    line_gen2 = vcf2.parse(options["file2"])
    
    # This will get moved to the outputbuffer class
    out = outputbuffer(vcf1,options)
    offset = options["relax_pos"]
    matcher = match_engine(vcf1, vcf2, options)
    stats = defaultdict(int)

    f2_buffer = deque()
    line2 = line_gen2.next()

    for line1 in line_gen1:        
        # First throw out all the lines in buffer 2 which are too now too far away. 
        # Move this to some outputbuffer logic. 
        for line in f2_buffer:
            if line1["pos"] - line["pos"] > offset:
                stats["added"] += 1
                out.add((line["pos"],"+",line))
        f2_buffer = [ line for line in f2_buffer if line1["pos"] - line["pos"] <= offset ]

        # Now add any new lines to the buffer. 
        while line2["pos"] - line1["pos"] <= offset:
            f2_buffer.append(line2)
            try: line2 = line_gen2.next()
            except StopIteration: break

        # now compare line1 to all the lines in the buffer and find the best match
        # See match codes defined above.
        matches = [ matcher.match( line1, line ) for line in f2_buffer ]    
        best_match = max(matches + [0])
        if best_match == MATCH_NO:
            stats["removed"] += 1
            out.add((line1["pos"],"-",line1))
        elif best_match == MATCH_MATCH:
            stats["matched"] += 1
            f2_buffer.pop(matches.index(best_match))
        elif best_match == MATCH_GTYPE:
            stats["matched/genotype"] += 1
            f2_buffer.pop(matches.index(best_match))
            out.add((line["pos"],"*",line1))
        else:
            raise Exception( "Don't know how to handle match code %d"%(best_match) )

    # Finally deal with anything that's left in the buffer - must be added in file2
    for line in f2_buffer:
        stats["added"] += 1
        out.add((line["pos"],"+",line))

    out.flush()
    print stats.items()

if __name__ == "__main__" : main()
