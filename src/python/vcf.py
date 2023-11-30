#
# Code to read, write and edit VCF files
#
# VCF lines are encoded as a dictionary with these keys (note: all lowercase):
# 'chrom':  string
# 'pos':    integer
# 'id':     string
# 'ref':    string
# 'alt':    list of strings
# 'qual':   integer
# 'filter': None (missing value), or list of keys (strings); empty list parsed as ["PASS"]
# 'info':   dictionary of values (see below)
# 'format': list of keys (strings)
# sample keys: dictionary of values (see below)
#
# The sample keys are accessible through vcf.getsamples()
#
# A dictionary of values contains value keys (defined in ##INFO or ##FORMAT lines) which map
# to a list, containign integers, floats, strings, or characters.  Missing values are replaced
# by a particular value, often -1 or .
#
# Genotypes are not stored as a string, but as a list of 1 or 3 elements (for haploid and diploid samples),
# the first (and last) the integer representing an allele, and the second the separation character.
# Note that there is just one genotype per sample, but for consistency the single element is stored in a list.
#
# Header lines other than ##INFO, ##FORMAT and ##FILTER are stored as (key, value) pairs and are accessible
# through getheader()
#
# The VCF class can be instantiated with a 'regions' variable consisting of tuples (chrom,start,end) encoding
# 0-based half-open segments.  Only variants with a position inside the segment will be parsed.  A regions
# parser is available under parse_regions.
#
# When instantiated, a reference can be passed to the VCF class.  This may be any class that supports a
# fetch(chrom, start, end) method.
#
#
#
# NOTE: the position that is returned to Python is 0-based, NOT 1-based as in the VCF file.
#
#
#
# TODO:
#  only v4.0 writing is complete; alleles are not converted to v3.3 format
#

from collections import namedtuple, defaultdict
from operator import itemgetter
import sys, re, copy, bisect

gtsRegEx = re.compile("[|/\\\\]")
alleleRegEx = re.compile('^[ACGTN]+$')
cleanupEx = re.compile('[^ACGTN]')

# Utility function.  Uses 0-based coordinates
def get_sequence(chrom, start, end, fa):
    # obtain sequence from .fa file, without truncation
    if end<=start: return ""
    if not fa: return "N"*(end-start)
    if start<0: return "N"*(-start) + get_sequence(chrom, 0, end, fa).upper()
    sequence = fa.fetch(chrom, start, end).upper()
    sequence = cleanupEx.sub('N',sequence)
    if len(sequence) < end-start: sequence += "N"*(end-start-len(sequence))
    return sequence

# Utility function.  Parses a region string
def parse_regions( string ):
    result = []
    for r in string.split(','):
        elts = r.split(':')
        chrom, start, end = elts[0], 0, 3000000000
        if len(elts)==1: pass
        elif len(elts)==2:
            if len(elts[1])>0:
                ielts = elts[1].split('-')
                if len(ielts) != 2: ValueError("Don't understand region string '%s'" % r)
                try:    start, end = int(ielts[0])-1, int(ielts[1])
                except Exception: raise ValueError("Don't understand region string '%s'" % r)
        else:
            raise ValueError("Don't understand region string '%s'" % r)
        result.append( (chrom,start,end) )
    return result


FORMAT = namedtuple('FORMAT','id numbertype number type description missingvalue')

###########################################################################################################
#
# New class
#
###########################################################################################################

class VCF:

    # types
    NT_UNKNOWN = 0
    NT_NUMBER = 1
    NT_ALLELES = 2
    NT_NR_ALLELES = 3
    NT_GENOTYPES = 4
    NT_PHASED_GENOTYPES = 5

    _errors = { 0:"UNKNOWN_FORMAT_STRING:Unknown file format identifier",
                1:"BADLY_FORMATTED_FORMAT_STRING:Formatting error in the format string",
                2:"BADLY_FORMATTED_HEADING:Did not find 9 required headings (CHROM, POS, ..., FORMAT) %s",
                3:"BAD_NUMBER_OF_COLUMNS:Wrong number of columns found (%s)",
                4:"POS_NOT_NUMERICAL:Position column is not numerical",
                5:"UNKNOWN_CHAR_IN_REF:Unknown character in reference field",
                6:"V33_BAD_REF:Reference should be single-character in v3.3 VCF",
                7:"V33_BAD_ALLELE:Cannot interpret allele for v3.3 VCF",
                8:"POS_NOT_POSITIVE:Position field must be >0",
                9:"QUAL_NOT_NUMERICAL:Quality field must be numerical, or '.'",
               10:"ERROR_INFO_STRING:Error while parsing info field",
               11:"ERROR_UNKNOWN_KEY:Unknown key (%s) found in formatted field (info; format; or filter)",
               12:"ERROR_FORMAT_NOT_NUMERICAL:Expected integer or float in formatted field; got %s",
               13:"ERROR_FORMAT_NOT_CHAR:Eexpected character in formatted field; got string",
               14:"FILTER_NOT_DEFINED:Identifier (%s) in filter found which was not defined in header",
               15:"FORMAT_NOT_DEFINED:Identifier (%s) in format found which was not defined in header",
               16:"BAD_NUMBER_OF_VALUES:Found too many of values in sample column (%s)",
               17:"BAD_NUMBER_OF_PARAMETERS:Found unexpected number of parameters (%s)",
               18:"BAD_GENOTYPE:Cannot parse genotype (%s)",
               19:"V40_BAD_ALLELE:Bad allele found for v4.0 VCF (%s)",
               20:"MISSING_REF:Reference allele missing",
               21:"V33_UNMATCHED_DELETION:Deleted sequence does not match reference (%s)",
               22:"V40_MISSING_ANGLE_BRACKETS:Format definition is not deliminted by angular brackets",
               23:"FORMAT_MISSING_QUOTES:Description field in format definition is not surrounded by quotes",
               24:"V40_FORMAT_MUST_HAVE_NAMED_FIELDS:Fields in v4.0 VCF format definition must have named fields",
               25:"HEADING_NOT_SEPARATED_BY_TABS:Heading line appears separated by spaces, not tabs",
               26:"WRONG_REF:Wrong reference %s",
               27:"ERROR_TRAILING_DATA:Numerical field ('%s') has semicolon-separated trailing data",
               28:"BAD_CHR_TAG:Error calculating chr tag for %s",
               29:"ZERO_LENGTH_ALLELE:Found zero-length allele",
               30:"MISSING_INDEL_ALLELE_REF_BASE:Indel alleles must begin with single reference base",
               31:"ERROR_NON_FLAG_WITHOUT_VALUE: Formatted field is not a flag, but does not have a value specified"
                }

    # required columns
    _required = ["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]

    def __init__(self, _copy=None, reference=None, regions=None, lines=None, leftalign=False, _fastGT=False):
        # make error identifiers accessible by name
        for id in list(self._errors.keys()): self.__dict__[self._errors[id].split(':')[0]] = id
        self._header = []         # tag-value pairs; tags are not unique; does not include fileformat, INFO, FILTER or FORMAT fields
        self._version = 40        # version number; 33=v3.3; 40=v4.0
        self._info = {}           # info, filter and format data
        self._filter = {}
        self._format = {}
        self._samples = []
        self._ignored_errors = set([11])   # ERROR_UNKNOWN_KEY
        self._warn_errors = set([])
        self._leftalign = False
        self._fastGT = False
        self._reference = None    # reference sequence
        self._regions = None      # regions to include; None includes everything
        self._lineno = -1
        self._line = None
        self._lines = None
        if _copy != None:
            self._header = _copy._header[:]
            self._version = _copy._version
            self._info = copy.deepcopy(_copy._info)
            self._filter = copy.deepcopy(_copy._filter)
            self._format = copy.deepcopy(_copy._format)
            self._samples = _copy._samples[:]
            self._ignored_errors = copy.deepcopy(_copy._ignored_errors)
            self._warn_errors = copy.deepcopy(_copy._warn_errors)
            self._leftalign = _copy._leftalign
            self._reference = _copy._reference
            self._regions = _copy._regions
        if reference: self._reference = reference
        if regions: self._regions = regions
        if leftalign: self._leftalign = leftalign
        if _fastGT: self._fastGT = _fastGT
        self._lines = lines

    def error(self,line,error,opt=None):
        if error in self._ignored_errors: return
        errorlabel, errorstring = self._errors[error].split(':')
        if opt: errorstring = errorstring % opt
        errwarn = ["Error","Warning"][error in self._warn_errors]
        sys.stderr.write("Line %s: '%s'\n%s %s: %s\n" % (self._lineno,line,errwarn,errorlabel,errorstring))
        if error in self._warn_errors: return
        raise ValueError(errorstring)

    def parse_format(self,line,format,filter=False):
        if self._version >= 40:
            if not format.startswith('<'):
                self.error(line,self.V40_MISSING_ANGLE_BRACKETS)
                format = "<"+format
            if not format.endswith('>'):
                self.error(line,self.V40_MISSING_ANGLE_BRACKETS)
                format += ">"
            format = format[1:-1]
        data = {'id':None,'number':None,'type':None,'descr':None}
        idx = 0
        while len(format.strip())>0:
            elts = format.strip().split(',')
            first, rest = elts[0], ','.join(elts[1:])
            if first.find('=') == -1 or (first.find('"')>=0 and first.find('=') > first.find('"')):
                if self._version >= 40: self.error(line,self.V40_FORMAT_MUST_HAVE_NAMED_FIELDS)
                if idx == 4: self.error(line,self.BADLY_FORMATTED_FORMAT_STRING)
                first = ["ID=","Number=","Type=","Description="][idx] + first
            if first.startswith('ID='):            data['id'] = first.split('=')[1]
            elif first.startswith('Number='):      data['number'] = first.split('=')[1]
            elif first.startswith('Type='):        data['type'] = first.split('=')[1]
            elif first.startswith('Description='):
                elts = format.split('"')
                if len(elts)<3:
                    self.error(line,self.FORMAT_MISSING_QUOTES)
                    elts = first.split('=') + [rest]
                data['descr'] = elts[1]
                rest = '"'.join(elts[2:])
                if rest.startswith(','): rest = rest[1:]
            else:
                self.error(line,self.BADLY_FORMATTED_FORMAT_STRING)
            format = rest
            idx += 1
            if filter and idx==1: idx=3  # skip number and type fields for FILTER format strings
        if not data['id']: self.error(line,self.BADLY_FORMATTED_FORMAT_STRING)
        if data['descr'] == None:
            self.error(line,self.BADLY_FORMATTED_FORMAT_STRING)
            data['descr'] = '<none>'
        if not data['type'] and not data['number']:
            # fine, ##filter format
            return FORMAT(data['id'],self.NT_NUMBER,0,"Flag",data['descr'],'.')
        if not data['type'] in ["Integer","Float","Character","String","Flag"]:
            self.error(line,self.BADLY_FORMATTED_FORMAT_STRING)
        # I would like a missing-value field, but it isn't there
        if data['type'] in ['Integer','Float']: data['missing'] = None    # Do NOT use arbitrary int/float as missing value
        else:                                   data['missing'] = '.'
        if not data['number']: self.error(line,self.BADLY_FORMATTED_FORMAT_STRING)
        try:
            n = int(data['number'])
            t = self.NT_NUMBER
        except ValueError:
            n = -1
            if data['number'] == '.':                   t = self.NT_UNKNOWN
            elif data['number'] == '#alleles':          t = self.NT_ALLELES
            elif data['number'] == 'A':                 t = self.NT_ALLELES
            elif data['number'] == '#nonref_alleles':   t = self.NT_NR_ALLELES
            elif data['number'] == '#genotypes':        t = self.NT_GENOTYPES
            elif data['number'] == 'G':                 t = self.NT_GENOTYPES
            elif data['number'] == '#phased_genotypes': t = self.NT_PHASED_GENOTYPES
            else:
                self.error(line,self.BADLY_FORMATTED_FORMAT_STRING)
        return FORMAT(data['id'],t,n,data['type'],data['descr'],data['missing'])


    def format_format( self, fmt, filter=False ):
        values = [('ID',fmt.id)]
        if fmt.number != None and not filter:
            if fmt.numbertype == self.NT_UNKNOWN: nmb = "."
            elif fmt.numbertype == self.NT_NUMBER: nmb = str(fmt.number)
            elif fmt.numbertype == self.NT_ALLELES: nmb = "#alleles"
            elif fmt.numbertype == self.NT_NR_ALLELES: nmb = "#nonref_alleles"
            elif fmt.numbertype == self.NT_GENOTYPES: nmb = "#genotypes"
            elif fmt.numbertype == self.NT_PHASED_GENOTYPES: nmb = "#phased_genotypes"
            else:
                raise ValueError("Unknown number type encountered: %s" % fmt.numbertype)
            values.append( ('Number',nmb) )
            values.append( ('Type', fmt.type) )
        values.append( ('Description', '"' + fmt.description + '"') )
        if self._version == 33:
            format = ",".join(v for k,v in values)
        else:
            format = "<" + (",".join( "%s=%s" % (k,v) for (k,v) in values )) + ">"
        return format

    def get_expected(self, format, formatdict, alt):
        fmt = formatdict.get(format,None)
        if fmt == None: return -1
        if fmt.numbertype == self.NT_UNKNOWN: return -1
        if fmt.numbertype == self.NT_NUMBER: return fmt.number
        if fmt.numbertype == self.NT_ALLELES: return len(alt)+1
        if fmt.numbertype == self.NT_NR_ALLELES: return len(alt)
        if fmt.numbertype == self.NT_GENOTYPES: return ((len(alt)+1)*(len(alt)+2)) // 2
        if fmt.numbertype == self.NT_PHASED_GENOTYPES: return (len(alt)+1)*(len(alt)+1)
        return 0


    def _add_definition(self, formatdict, key, data, line ):
        if key in formatdict: return
        self.error(line,self.ERROR_UNKNOWN_KEY,key)
        if data == None:
            formatdict[key] = FORMAT(key,self.NT_NUMBER,0,"Flag","(Undefined tag)",".")
            return
        if data == []: data = [""]             # unsure what type -- say string
        if type(data[0]) == type(0.0):
            formatdict[key] = FORMAT(key,self.NT_UNKNOWN,-1,"Float","(Undefined tag)",None)
            return
        if type(data[0]) == type(0):
            formatdict[key] = FORMAT(key,self.NT_UNKNOWN,-1,"Integer","(Undefined tag)",None)
            return
        formatdict[key] = FORMAT(key,self.NT_UNKNOWN,-1,"String","(Undefined tag)",".")


    # todo: trim trailing missing values
    def format_formatdata( self, data, format, key=True, value=True, separator=":" ):
        output, sdata = [], []
        if type(data) == type([]): # for FORMAT field, make data with dummy values
            d = {}
            for k in data: d[k] = []
            data = d
        # convert missing values; and silently add definitions if required
        for k in data:
            self._add_definition( format, k, data[k], "(output)" )
            for idx,v in enumerate(data[k]):
                if v == format[k].missingvalue: data[k][idx] = "."
        # make sure GT comes first; and ensure fixed ordering; also convert GT data back to string
        for k in data:
            if k != 'GT': sdata.append( (k,data[k]) )
        sdata.sort()
        if 'GT' in data:
            sdata = [('GT',list(map(self.convertGTback,data['GT'])))] + sdata
        for k,v in sdata:
            if v == []: v = None
            if key and value:
                if v != None: output.append( k+"="+','.join(map(str,v)) )
                else: output.append( k )
            elif key: output.append(k)
            elif value:
                if v != None: output.append( ','.join(map(str,v)) )
                else: output.append( "." )                    # should not happen
        # snip off trailing missing data
        while len(output) > 1:
            last = output[-1].replace(',','').replace('.','')
            if len(last)>0: break
            output = output[:-1]
        return separator.join(output)


    def enter_default_format(self):
        for f in [FORMAT('GT',self.NT_NUMBER,1,'String','Genotype','.'),
                  FORMAT('GQ',self.NT_NUMBER,1,'Integer','Genotype Quality',-1),
                  FORMAT('DP',self.NT_NUMBER,1,'Integer','Read depth at this position for this sample',-1),
                  FORMAT('HQ',self.NT_UNKNOWN,-1,'Integer','Haplotype Quality',-1),    # unknown number, since may be haploid
                  FORMAT('FT',self.NT_NUMBER,1,'String','Sample Genotype Filter','.')]:
            if f.id not in self._format:
                self._format[f.id] = f


    def parse_header( self, line ):
        assert line.startswith('##')
        elts = line[2:].split('=')
        key = elts[0].strip()
        value = '='.join(elts[1:]).strip()
        if key == "fileformat":
            if value == "VCFv3.3":
                self._version = 33
            elif value == "VCFv4.0":
                self._version = 40
            # Temp fudge. TODO: fix this properly
            elif value == "VCFv4.1":
                self._version = 41
            else:
                print(value)
                self.error(line,self.UNKNOWN_FORMAT_STRING)
        elif key == "INFO":
            f = self.parse_format(line, value)
            self._info[ f.id ] = f
        elif key == "FILTER":
            f = self.parse_format(line, value, filter=True)
            self._filter[ f.id ] = f
        elif key == "FORMAT":
            f = self.parse_format(line, value)
            self._format[ f.id ] = f
        else:
            # keep other keys in the header field
            self._header.append( (key,value) )


    def write_header( self, stream ):
        stream.write("##fileformat=VCFv%s.%s\n" % (self._version // 10, self._version % 10))
        for key,value in self._header: stream.write("##%s=%s\n" % (key,value))
        for var,label in [(self._info,"INFO"),(self._filter,"FILTER"),(self._format,"FORMAT")]:
            for f in var.values(): stream.write("##%s=%s\n" % (label,self.format_format(f,filter=(label=="FILTER"))))


    def parse_heading( self, line ):
        assert line.startswith('#')
        assert not line.startswith('##')
        headings = line[1:].split('\t')
        if len(headings)==1 and len(line[1:].split()) >= 9:
            self.error(line,self.HEADING_NOT_SEPARATED_BY_TABS)
            headings = line[1:].split()

        for i,s in enumerate(self._required):

            if len(headings)<=i or headings[i] != s:

                if len(headings) <= i:
                    err = "(%sth entry not found)" % (i+1)
                else:
                    err = "(found %s, expected %s)" % (headings[i],s)

                #self.error(line,self.BADLY_FORMATTED_HEADING,err)

                # allow FORMAT column to be absent
                if len(headings) == 8:
                    headings.append("FORMAT")
                else:
                    self.error(line,self.BADLY_FORMATTED_HEADING,err)

        self._samples = headings[9:]


    def write_heading( self, stream ):
        stream.write("#" + "\t".join(self._required + self._samples) + "\n")


    def convertGT(self, GTstring):
        if GTstring == "." or GTstring == "1" or GTstring == "" or GTstring.startswith(".:"): return ["."]
        if self._fastGT:
            try:
                if GTstring[0] == '.' or GTstring[2] == '.': return ["."]
                return [int(GTstring[0]),GTstring[1],int(GTstring[2])]
            except Exception:
                # hack, to accept "0" genotypes
                return [int(GTstring[0])]
        try:
            gts = gtsRegEx.split(GTstring)
            if len(gts) == 1: return [int(gts[0])]
            if len(gts) != 2: raise ValueError()
            if gts[0] == "." and gts[1] == ".": return [gts[0],GTstring[len(gts[0]):-len(gts[1])],gts[1]]
            return [int(gts[0]),GTstring[len(gts[0]):-len(gts[1])],int(gts[1])]
        except ValueError:
            self.error(self._line,self.BAD_GENOTYPE,GTstring)
            return [".","|","."]


    def convertGTback(self, GTdata):
        return ''.join(map(str,GTdata))


    def parse_formatdata( self, key, value, formatdict, line ):
        # To do: check that the right number of values is present
        f = formatdict.get(key,None)
        if f == None:
            self._add_definition(formatdict, key, value, line )
            f = formatdict[key]
        if f.type == "Flag":
            if value is not None: self.error(line,self.ERROR_FLAG_HAS_VALUE)
            return []
        elif value is None: 
            self.error(line, self.ERROR_NON_FLAG_WITHOUT_VALUE)
            return []
        values = value.split(',')
        # deal with trailing data in some early VCF files
        if f.type in ["Float","Integer"] and len(values)>0 and values[-1].find(';') > -1:
            self.error(line,self.ERROR_TRAILING_DATA,values[-1])
            values[-1] = values[-1].split(';')[0]
        if f.type == "Integer":
            try:
                for idx,v in enumerate(values):
                    if v == ".": values[idx] = f.missingvalue
                    else:        values[idx] = int(v)
            except Exception:
                self.error(line,self.ERROR_FORMAT_NOT_NUMERICAL,values)
                return [0] * len(values)
            return values
        elif f.type == "String":
            self._line = line
            if f.id == "GT": values = list(map( self.convertGT, values ))
            return values
        elif f.type == "Character":
            for v in values:
                if len(v) != 1: self.error(line,self.ERROR_FORMAT_NOT_CHAR)
            return values
        elif f.type == "Float":
            try:
                for idx,v in enumerate(values):
                    if v == ".":  values[idx] = f.missingvalue
                    else:         values[idx] = float(values[idx])
            except Exception:
                self.error(line,self.ERROR_FORMAT_NOT_NUMERICAL,values)
                return [0.0] * len(values)
            return values
        else:
            # can't happen
            self.error(line,self.ERROR_INFO_STRING)


    def inregion(self, chrom, pos):
        if not self._regions: return True
        for r in self._regions:
            if r[0] == chrom and r[1] <= pos < r[2]: return True
        return False


    def parse_data( self, line, lineparse=False, parseGenotypes=True):
        cols = line.split('\t')
        if len(cols) != len(self._samples)+9:
            # gracefully deal with absent FORMAT column
            if len(cols) == 8 and len(self._samples)==0:
                cols.append("")
            else:
                self.error(line,
                           self.BAD_NUMBER_OF_COLUMNS,
                           "expected %s for %s samples (%s), got %s" % (len(self._samples)+9, len(self._samples), self._samples, len(cols)))

        chrom = cols[0]

        # get 0-based position
        try:    pos = int(cols[1])-1
        except Exception: self.error(line,self.POS_NOT_NUMERICAL)
        if pos < 0: self.error(line,self.POS_NOT_POSITIVE)

        # implement filtering
        if not self.inregion(chrom,pos): return None

        # end of first-pass parse for sortedVCF
        if lineparse: return chrom, pos, line

        id = cols[2]

        ref = cols[3].upper()
        if ref == ".":
            self.error(line,self.MISSING_REF)
            if self._version == 33: ref = get_sequence(chrom,pos,pos+1,self._reference)
            else:                   ref = ""
        else:
            for c in ref:
                if c not in "ACGTN": self.error(line,self.UNKNOWN_CHAR_IN_REF)
            if "N" in ref: ref = get_sequence(chrom,pos,pos+len(ref),self._reference)

        # make sure reference is sane
        if self._reference:
            left = max(0,pos-100)
            faref_leftflank = get_sequence(chrom,left,pos+len(ref),self._reference)
            faref = faref_leftflank[pos-left:]
            if faref != ref: self.error(line,self.WRONG_REF,"(reference is %s, VCF says %s)" % (faref,ref))
            ref = faref

        # convert v3.3 to v4.0 alleles below
        if cols[4] == ".": alt = []
        else: alt = cols[4].upper().split(',')

        if cols[5] == ".": qual = -1
        else:
            try:    qual = float(cols[5])
            except Exception: self.error(line,self.QUAL_NOT_NUMERICAL)

        # postpone checking that filters exist.  Encode missing filter or no filtering as empty list
        if cols[6] == "." or cols[6] == "PASS" or cols[6] == "0": filter = []
        else: filter = cols[6].split(';')

        # dictionary of keys, and list of values
        info = {}
        if cols[7] != "." and cols[7] != "":
            for blurp in cols[7].split(';'):
                elts = blurp.split('=')
                if len(elts) == 1: v = None
                elif len(elts) == 2: v = elts[1]
                else: self.error(line,self.ERROR_INFO_STRING)
                info[elts[0]] = self.parse_formatdata(elts[0], v, self._info, line)

        # Gracefully deal with absent FORMAT column
        if cols[8] == "": format = []
        else: format = cols[8].split(':')

        # hack: only keep genotype if doing fast parsing
        if self._fastGT:
            format = format[:1]

        # check: all filters are defined
        if filter is not None:
            for f in filter:
                if f not in self._filter: self.error(line,self.FILTER_NOT_DEFINED, f)

        # check: format fields are defined
        for f in format:
            if f not in self._format: self.error(line,self.FORMAT_NOT_DEFINED, f)

        # convert v3.3 alleles
        if self._version == 33:
            if len(ref) != 1: self.error(line,self.V33_BAD_REF)
            newalts = []
            have_deletions = False
            for a in alt:
                if len(a) == 1: a = a + ref[1:]                       # SNP; add trailing reference
                elif a.startswith('I'): a = ref[0] + a[1:] + ref[1:]  # insertion just beyond pos; add first and trailing reference
                elif a.startswith('D'): # allow D<seq> and D<num>
                    have_deletions = True
                    try:
                        l = int(a[1:])          # throws ValueError if sequence
                        if len(ref) < l:        # add to reference if necessary
                            addns = get_sequence(chrom,pos+len(ref),pos+l,self._reference)
                            ref += addns
                            for i,na in enumerate(newalts): newalts[i] = na+addns
                        a = ref[l:]             # new deletion, deleting pos...pos+l
                    except ValueError:
                        s = a[1:]
                        if len(ref) < len(s):   # add Ns to reference if necessary
                            addns = get_sequence(chrom,pos+len(ref),pos+len(s),self._reference)
                            if not s.endswith(addns) and addns != 'N'*len(addns):
                                self.error(line,self.V33_UNMATCHED_DELETION,
                                           "(deletion is %s, reference is %s)" % (a,get_sequence(chrom,pos,pos+len(s),self._reference)))
                            ref += addns
                            for i,na in enumerate(newalts): newalts[i] = na+addns
                        a = ref[len(s):]        # new deletion, deleting from pos
                else:
                    self.error(line,self.V33_BAD_ALLELE)
                newalts.append(a)
            alt = newalts
            # deletion alleles exist, add dummy 1st reference allele, and account for leading base
            if have_deletions:
                if pos == 0:
                    # Petr Danacek's: we can't have a leading nucleotide at (1-based) position 1
                    addn = get_sequence(chrom,pos+len(ref),pos+len(ref)+1,self._reference)
                    ref += addn
                    alt = [allele+addn for allele in alt]
                else:
                    addn = get_sequence(chrom,pos-1,pos,self._reference)
                    ref = addn + ref
                    alt = [addn + allele for allele in alt]
                    pos -= 1
        else:
            # format v4.0 -- just check for nucleotides
            for allele in alt:
                # hack -- allowed ALT alleles are defined in the header
                if self._version == 41 and allele == "<DEL>":
                    continue
                if not alleleRegEx.match(allele):
                    self.error(line,self.V40_BAD_ALLELE,allele)

        # check for leading nucleotide in indel calls
        for allele in alt:
            if allele == "<DEL>": continue
            if len(allele) != len(ref):
                if len(allele) == 0: self.error(line,self.ZERO_LENGTH_ALLELE)
                if ref[0].upper() != allele[0].upper() and "N" not in (ref[0]+allele[0]).upper():
                    self.error(line,self.MISSING_INDEL_ALLELE_REF_BASE)

        # trim trailing bases in alleles
        if alt == []:
            pass
        else:
            for i in range(1,min(len(ref),min(list(map(len,alt))))):
                if len(set(allele[-1].upper() for allele in alt)) > 1 or ref[-1].upper() != alt[0][-1].upper():
                    break
                ref, alt = ref[:-1], [allele[:-1] for allele in alt]

        # left-align alleles, if a reference is available
        if self._leftalign and self._reference:
            while left < pos:
                movable = True
                for allele in alt:
                    if len(allele) > len(ref):
                        longest, shortest = allele, ref
                    else:
                        longest, shortest = ref, allele
                    if len(longest) == len(shortest) or longest[:len(shortest)].upper() != shortest.upper():
                        movable = False
                    if longest[-1].upper() != longest[len(shortest)-1].upper():
                        movable = False
                if not movable:
                    break
                ref = ref[:-1]
                alt = [allele[:-1] for allele in alt]
                if min(len(allele) for allele in alt) == 0 or len(ref) == 0:
                    ref = faref_leftflank[pos-left-1] + ref
                    alt = [faref_leftflank[pos-left-1] + allele for allele in alt]
                    pos -= 1

        # parse sample columns
        if parseGenotypes:
            samples = []
            for sample in cols[9:]:

                if self._fastGT:
                    # hack: ignore all information except genotype; assume diploid; do not check anything
                    samples.append( {format[0]: [self.convertGT(sample)]} )
                    continue

                dict = {}
                values = sample.split(':')
                if len(values) > len(format):
                    self.error(line,self.BAD_NUMBER_OF_VALUES,"(found %s values in element %s; expected %s)" % (len(values),sample,len(format)))
                for idx in range(len(format)):
                    expected = self.get_expected(format[idx], self._format, alt)
                    if idx < len(values): value = values[idx]
                    else:
                        if expected == -1: value = "."
                        else: value = ",".join(["."]*expected)
                    dict[format[idx]] = self.parse_formatdata(format[idx], value, self._format, line)
                    if expected != -1 and len(dict[format[idx]]) != expected:
                        self.error(line,self.BAD_NUMBER_OF_PARAMETERS,
                                   "id=%s, expected %s parameters, got %s" % (format[idx],expected,dict[format[idx]]))
                        if len(dict[format[idx]] ) < expected: dict[format[idx]] += [dict[format[idx]][-1]]*(expected-len(dict[format[idx]]))
                        dict[format[idx]] = dict[format[idx]][:expected]
                samples.append( dict )

        # done
        d = {'chrom':chrom,
             'pos':pos,      # return 0-based position
             'id':id,
             'ref':ref,
             'alt':alt,
             'qual':qual,
             'filter':filter,
             'info':info,
             'format':format}

        if parseGenotypes:
            for key,value in zip(self._samples,samples):
                d[key] = value

        return d


    def write_data(self, stream, data):
        required = ['chrom','pos','id','ref','alt','qual','filter','info','format'] + self._samples
        for k in required:
            if k not in data: raise ValueError("Required key %s not found in data" % str(k))
        if data['alt'] == []: alt = "."
        else: alt = ",".join(data['alt'])
        if data['filter'] == None: filter = "."
        elif data['filter'] == []:
            if self._version == 33: filter = "0"
            else: filter = "PASS"
        else: filter = ';'.join(data['filter'])
        if data['qual'] == -1: qual = "."
        else:
            qual = str(data['qual'])
            if qual.endswith(".0"): qual = qual[:-2]    # report as integer if possible

        output = [data['chrom'],
                  str(data['pos']+1),   # change to 1-based position
                  data['id'],
                  data['ref'],
                  alt,
                  qual,
                  filter,
                  self.format_formatdata( data['info'], self._info, separator=";" ),
                  self.format_formatdata( data['format'], self._format, value=False ) ]

        for s in self._samples:
            output.append( self.format_formatdata( data[s], self._format, key=False ) )

        stream.write( "\t".join(output) + "\n" )

    def _parse_header(self, stream):
        self._lineno = 0
        line = None

        for line in stream:
            self._lineno += 1
            if line.startswith('##'):
                self.parse_header( line.strip() )
            elif line.startswith('#'):
                self.parse_heading( line.strip() )
                self.enter_default_format()
            else:
                break
        return line

    def _parse(self, line, stream, parseGenotypes):
        if len(line.strip()) > 0:
            d = self.parse_data( line.strip(), parseGenotypes=parseGenotypes)
            if d: yield d
        for line in stream:
            self._lineno += 1
            if self._lines and self._lineno > self._lines: raise StopIteration
            d = self.parse_data( line.strip(), parseGenotypes=parseGenotypes)
            if d: yield d

    ######################################################################################################
    #
    # API follows
    #
    ######################################################################################################

    def getsamples(self):
        """ List of samples in VCF file """
        return self._samples

    def setsamples(self,samples):
        """ List of samples in VCF file """
        self._samples = samples

    def getheader(self):
        """ List of header key-value pairs (strings) """
        return self._header

    def setheader(self,header):
        """ List of header key-value pairs (strings) """
        self._header = header

    def getinfo(self):
        """ Dictionary of ##INFO tags, as VCF.FORMAT values """
        return self._info

    def setinfo(self,info):
        """ Dictionary of ##INFO tags, as VCF.FORMAT values """
        self._info = info

    def getformat(self):
        """ Dictionary of ##FORMAT tags, as VCF.FORMAT values """
        return self._format

    def setformat(self,format):
        """ Dictionary of ##FORMAT tags, as VCF.FORMAT values """
        self._format = format

    def getfilter(self):
        """ Dictionary of ##FILTER tags, as VCF.FORMAT values """
        return self._filter

    def setfilter(self,filter):
        """ Dictionary of ##FILTER tags, as VCF.FORMAT values """
        self._filter = filter

    def setversion(self, version):
        if version != 33 and version != 40 and version != 41: raise ValueError("Can only handle v3.3 and v4.0/v4.1 VCF files")
        self._version = version

    def setregions(self, regions):
        self._regions = regions

    def setreference(self, ref):
        """ Provide a reference sequence; a Python class supporting a fetch(chromosome, start, end) method, e.g. PySam.FastaFile """
        self._reference = ref

    def ignoreerror(self, errorstring):
        try:             self._ignored_errors.add(self.__dict__[errorstring])
        except KeyError: raise ValueError("Invalid error string: %s" % errorstring)

    def warnerror(self, errorstring):
        try:             self._warn_errors.add(self.__dict__[errorstring])
        except KeyError: raise ValueError("Invalid error string: %s" % errorstring)

    def parse(self, stream, parseGenotypes=True):
        """ Parse a stream of VCF-formatted lines.  Initializes class instance and return generator """
        last_line = self._parse_header(stream)
        # now return a generator that does the actual work.  In this way the pre-processing is done
        # before the first piece of data is yielded
        return self._parse(last_line, stream, parseGenotypes)

    def write(self, stream, datagenerator):
        """ Writes a VCF file to a stream, using a data generator (or list) """
        self.write_header(stream)
        self.write_heading(stream)
        for data in datagenerator: self.write_data(stream,data)

    def writeheader(self, stream):
        """ Writes a VCF header """
        self.write_header(stream)
        self.write_heading(stream)

    def compare_calls(self, pos1, ref1, alt1, pos2, ref2, alt2):
        """ Utility function: compares two calls for equality """
        # a variant should always be assigned to a unique position, one base before
        # the leftmost position of the alignment gap.  If this rule is implemented
        # correctly, the two positions must be equal for the calls to be identical.
        if pos1 != pos2: return False
        # from both calls, trim rightmost bases when identical.  Do this safely, i.e.
        # only when the reference bases are not Ns
        while len(ref1)>0 and len(alt1)>0 and ref1[-1] == alt1[-1]:
            ref1 = ref1[:-1]
            alt1 = alt1[:-1]
        while len(ref2)>0 and len(alt2)>0 and ref2[-1] == alt2[-1]:
            ref2 = ref2[:-1]
            alt2 = alt2[:-1]
        # now, the alternative alleles must be identical
        return alt1 == alt2

###########################################################################################################
#
# New class
#
###########################################################################################################

class sortedVCF(VCF):
    """
    A vcf reader/write which reads the entire VCF file into memory and sorts it,
    storing it in a dictionary of the form {chrom: [ lines ]}. It can be accesed
    using parse as in the VCF class, in which case it will yield the lines sorted
    within chromosome and with some determinisitc ordering on the chromosomes.
    Alternatively, there is another interface which just gives you the data
    """
    def __init__(self):
        """
        """
        self._sorted_lines = None
        VCF.__init__(self)

    def load_and_sort_data(self, line, stream):
        """
        Read in the data from stream - not including the header, and store it
        in the _sorted_lines cache. This is a dictionary of chromosome to array
        of lines, where the lines are sorted. For tedious reasons, needs to use
        the first line separately.
        """
        if len(line.strip()) > 0:
            line_data = self.parse_data( line.strip(), lineparse=True, parseGenotypes=True)
            if line_data: self._sorted_lines[ line_data[0] ].append( (line_data[1],line_data[2]) )

        for line in stream:
            line_data = self.parse_data( line.strip(), lineparse=True, parseGenotypes=True)
            if line_data: self._sorted_lines[ line_data[0] ].append( (line_data[1],line_data[2]) )

        for key in self._sorted_lines.keys():
            self._sorted_lines[key].sort()

    def chr_tag(self, chrom):
        """
        remove initial chr, then either get an int or a string for the rest
        fortunately int < string always so this gives us the ordering we want
        returns a tuple so we can sort but get back the original keys.
        """
        if chrom[0:3].upper() == "CHR": val = chrom[3:]
        else: val = chrom

        try:
            return (chrom,int(val))
        except ValueError:
            return (chrom,val)

        self.error(None, self.BAD_CHR_TAG, chrom)

    def chr_order(self, chroms=None):
        """
        The order in which to return the chromosomes. Should sort by chromosome
        number then by letter i.e 1-21, M, X, Y... etc.
        If a list of chromosomes is provided, this list is sorted instead.
        """
        if chroms == None: chroms = list(self._sorted_lines.keys())
        chrs = list(map(self.chr_tag, chroms))
        chrs.sort(key=itemgetter(1))
        return list(map(itemgetter(0), chrs))

    def _parse(self, region=None):
        if region != None:
            (chr, start, end)=region
            poss=[dat[0] for dat in self._sorted_lines[chr]]
            left = bisect.bisect_left(poss, start)
            right = bisect.bisect_right(poss, end)
            for pos, line in self._sorted_lines[chr][left:right]:
                yield self.parse_data( line, parseGenotypes=True)
        else:
            for chrom in self.chr_order():
                for pos, line in self._sorted_lines[chrom]:
                    yield self.parse_data( line, parseGenotypes=True)

    def getdata(self):
        return self._sorted_lines

    def parse(self, stream, region=None):
        """ Parse a stream of VCF-formatted lines.  Initializes class instance and return generator """
        if self._sorted_lines is None:
            self._sorted_lines = defaultdict(list)
            last_line = self._parse_header(stream)
            self.load_and_sort_data(last_line, stream)

        return self._parse(region)

def Test():

    vcf33 = """##fileformat=VCFv3.3
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##INFO=NS,1,Integer,"Number of Samples With Data"
##INFO=DP,1,Integer,"Total Depth"
##INFO=AF,-1,Float,"Allele Frequency"
##INFO=AA,1,String,"Ancestral Allele"
##INFO=DB,0,Flag,"dbSNP membership, build 129"
##INFO=H2,0,Flag,"HapMap2 membership"
##FILTER=q10,"Quality below 10"
##FILTER=s50,"Less than 50% of samples have data"
##FORMAT=GT,1,String,"Genotype"
##FORMAT=GQ,1,Integer,"Genotype Quality"
##FORMAT=DP,1,Integer,"Read Depth"
##FORMAT=HQ,2,Integer,"Haplotype Quality"
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003
20\t14370\trs6054257\tG\tA\t29\t0\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:-1,-1
17\t17330\t.\tT\tA\t3\tq10\tNS=3;DP=11;AF=0.017\tGT:GQ:DP:HQ\t0|0:49:3:58,50\t0|1:3:5:65,3\t0/0:41:3:-1,-1
20\t1110696\trs6040355\tA\tG,T\t67\t0\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t1|2:21:6:23,27\t2|1:2:0:18,2\t2/2:35:4:-1,-1
17\t1230237\t.\tT\t.\t47\t0\tNS=3;DP=13;AA=T\tGT:GQ:DP:HQ\t0|0:54:7:56,60\t0|0:48:4:51,51\t0/0:61:2:-1,-1
20\t1234567\tmicrosat1\tG\tD4,IGA\t50\t0\tNS=3;DP=9;AA=G\tGT:GQ:DP\t0/1:35:4\t0/2:17:2\t1/1:40:3"""

    vcf40 = """##fileformat=VCFv4.0
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003
M\t1230237\t.\tT\t.\t47\tPASS\tNS=3;DP=13;AA=T\tGT:GQ:DP:HQ\t0|0:54:7:56,60\t0|0:48:4:51,51\t0/0:61:2
20\t1234567\tmicrosat1\tGTCT\tG,GTACT\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t0/1:35:4\t0/2:17:2\t1/1:40:3
17\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,.
20\t17330\t.\tT\tA\t3\tq10\tNS=3;DP=11;AF=0.017\tGT:GQ:DP:HQ\t0|0:49:3:58,50\t0|1:3:5:65,3\t0/0:41:3
20\t1110696\trs6040355\tA\tG,T\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t1|2:21:6:23,27\t2|1:2:0:18,2\t2/2:35:4"""

    if False:
        print("Parsing v3.3 file:")
        print(vcf33)
        vcf = VCF()
        lines = [data for data in vcf.parse( (line+"\n" for line in vcf33.split('\n') ) )]
        print("Writing v3.3 file:")
        vcf.write( sys.stdout, lines )

    if False:
        print("Parsing v4.0 file:")
        print(vcf40)
        vcf = VCF()
        lines = [data for data in vcf.parse( (line+"\n" for line in vcf40.split('\n') ) )]
        print("Writing v4.0 file:")
        vcf.write( sys.stdout, lines )

    if True:
        print("Parsing v3.3 file:")
        print(vcf33)
        vcf = sortedVCF()
        lines = [data for data in vcf.parse( (line+"\n" for line in vcf33.split('\n') ) )]
        print("Writing v3.3 file:")
        vcf.write( sys.stdout, lines )

    if True:
        print("Parsing v4.0 file:")
        print(vcf40)
        vcf = sortedVCF()
        lines = [data for data in vcf.parse( (line+"\n" for line in vcf40.split('\n') ) )]
        print("Writing v4.0 file:")
        vcf.write( sys.stdout, lines )

if __name__ == "__main__":
    Test()
