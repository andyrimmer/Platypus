###################################################################################################
#
# Platypus readme file
#
###################################################################################################


1. Requirements
===============

Platypus requires Python 2.6 or greater. It has not tested with the Python 3.X
series.

htslib 1.2.1 or greater: See the htslib website (http://www.htslib.org/download/) for details on how to download and install htslib.

2. Installation
===============

To build Platypus, run the 'buildPlatypus.sh' script included in the distribution:

cd Platypus_X.X.X
./buildPlatypus.sh

and this will compile all the required libraries. Platypus can then be run from the Platypus_X.X.X
directory as follows:


3. Running Platypus
===================

The easiest way to run Platypus is as follows:

python Platypus.py callVariants --bamFiles=LIST_OF_BAMS --refFile=REF.fa --output=Calls.vcf


You can see a list of all the possible input options by running the following comand:

python Platypus.py callVariants --help

However, in most cases the default parameter values should be fine, and you will only need to specify the --bamFiles
and --refFile and --output arguments. By default, if you do not specify a region or regions or interest, Platypus will
run through all the data in your BAM files. The --regions argument can be used to specify regions of interest.


3. Running in Variant-Calling Mode
==================================

The standard way of running Platypus is to use it to detect variants in one or more BAM files. Variants are detected by
comparing the BAM reads with a reference sequence. This can be done using the following command:

python Platypus.py callVariants --bamFiles=DATA.bam --regions=chr20 --output=test.vcf --refFile=GENOME.fa

where the input BAM files, and the genome reference must be indexed using samtools, or a program that produces compatible
index files.


4. Variant calling with additional Reference calling
====================================================

Platypus can also output reference calls. When a region is well covered by reads, and there is no evidence
of variation from the reference, a 'REFCALL' block will be output. This can be useful if you want to exclude
the possibility of any variation in a specific region. Each REFCALL block comes with an associated quality score,
in the 'QUAL' column. If this is high, then there is good support for the reference sequence; if this score is low,
then there is some evidence for variation, but not enough for Platypus to make an explicit variant call.

To enable reference-calling, use the '--outputRefCalls=1' flag:

python Platypus.py callVariants --bamFiles=DATA.bam --regions=chr20 --output=test.vcf --refFile=GENOME.fa --outputRefCalls=1


5. Running in Genotyping Mode
=============================

Platypus can take as input a compressed, indexed VCF file, and genotype the BAMs for all alleles in the compressed, indexed VCF.
To run Platypus in genotyping mode, use the following command:

python Platypus.py callVariants --bamFiles=DATA.bam --regions=chr20 --output=test.vcf --refFile=GENOME.fa --source=INPUT.vcf.gz --minPosterior=0 --getVariantsFromBAMs=0

You must use as input a VCF that has been compressed, with bgzip, and indexed with tabix. To create this, do the following:

bgzip file.vcf # Converts file.vcf into file.vcf.gz
tabix -p vcf file.vcf.gz # Creates the index file.vcf.gz.tbi


5.1 Known issues with genotyping
================================

Sometimes variants in the input VCF don't get genotypes in the output VCF. This generally only happens in complex regions, with many
variant candidates, when the particular variant is not well supported by the data, typically when there are lots of indels close together.


6. Running in Combined Genotyping/Calling Mode
==============================================

If required, Platypus can detect variants from the input BAM files, as well as using an input allele list. This can be done as
follows:

python Platypus.py callVariants --bamFiles=DATA.bam --regions=chr20 --output=test.vcf --refFile=GENOME.fa --source=INPUT.vcf.gz

For cases 3 and 4, if you want '0/0', i.e. reference genotypes to be reported, then set the argument --minPosterior=0, and the output VCF
will then also contain records for variants candidates which were genotyped '0/0'.


7. Main command-line arguments to Platypus
==========================================

# General variant-calling arguments
--output                   # Name of output VCF file
--logFileName              # Name of output log file (default is log.txt)
--refFile                  # Name of indexed FASTA reference file
--regions                  # List of regions to search. The format is e.g. chrX:1000-100000. This argument can be a comma-separated list of regions, or a text file of regions in the same format.
--bamFiles                 # List of BAM files. Currently must be either a comma-separated list of indexed BAMs, or the name of a text file with a list of input BAMs, one per line.
--bufferSize               # Specifies how much (as a genomic region) of the BAMs to buffer into memory at any time (default is 100kb)
--minReads                 # The minimum number of reads required to support a variants (default is 2)
--maxReads                 # The maximium allowed coverage in a region of "bufferSize". If there are more reads than this, the region will be skipped, and a warning issued (default is 5,000,000)
--maxVariants              # The maximium number of variants to consider in a given window (default is 8)
--verbosity                # Level of logging (default is 2. Set to 3 for debug output)
--source                   # Name of input VCF file to use for genotyping (default is None)
--nCPU                     # Number of processors to use (default is 1. If > 1, then multiple processes are run in parallel, and the output is combined at the end)
--getVariantsFromBAMs      # If set to 1 (default), variant candidates will be generated from BAMs as well as any other inputs
--genSNPs                  # If set to 1 (default), SNP candidates will be considered
--genIndels                # If set to 1 (default), Indel candidates will be considered
--minPosterior             # Only variants with posterior >= minPosterior will be output to the VCF. (default is 5; phred-scaled)
--maxSize                  # Only variant candidates smaller than this will be considered. Anything larger is filtered out (default is 1500 bases)
--minFlank                 # Variant candidates must be > minFlank bases from the end of a read to be considered (default is 10 bases).


# Arguments for BAM data filtering.
--minMapQual                          # Minimum mapping quality of reads to consider. Any reads with map qual below this are ignored (default is 20)
--minBaseQual                         # Minimum allowed base-calling quality. Any bases with qual below this are ignored in SNP-calling (default is 20)
--minGoodQualBases                    # Minimum number of bases with quality above --minBaseQual that must be present for a read to be used (default is 20).
--filterDuplicates                    # If set to 1 then Platypus will skip duplicate read-pairs based on the start and end position of both reads (default is 1).
--filterReadsWithUnmappedMates        # If set to 1, Platypus filters reads whose mates are un-mapped (default is 1).
--filterReadsWithDistantMates         # If set to 1, Platypus filters reads whose mates are mapped far away, using the IS_PROPER_PAIR flag in the BAM record. (Default is 1).
--filterReadPairsWithSmallInserts     # If set to 1, Platypus filters read-pairs with insert sizes < one read length (default is 1).


# Arguments for local assembly
--assemble                 # Set to 1 to turn on local assembly. Default is 0

# Arguments for referene calling
--outputRefCalls          # If set to 1, Platypus will output reference call blocks (default is 0).


These are not all the possible arguments. To get a complete list, run "python Platypus.py callVariants --help". But anything not listed here
should generally be left alone.


8. VCF output
=============

The VCF files output by Platypus contain a number of annotations, some of which are deprecated and will be removed in future releases.


8.1 INFO fields
===============

FR        # Estimated haplotype population frequency
TC        # Total coverage at this locus
TCR       # Total reverse strand coverage at this locus
TCF       # Total forward strand coverage at this locus
NR        # Total number of reverse reads containing this variant
NF        # Total number of forward reads containing this variant
TR        # Total number of reads containing this variant
HP        # Homopolmer run length in 20 bases either side of variant position
PP        # Posterior probability (phred scaled) that this variant segregates in the data.
SC        # Genomic sequence 10 bases either side of variant position
MMLQ      # Median minimum base quality for bases around variant. If this is low (<=10) then the variant is only supported by low-quality reads.
WS        # Start position of window in which variant was called
WE        # End position of window in which variant was called
SOURCE    # Flag to say if the variant was found by Platypus, the Assembler, or an input file
END       # End position of reference call-block. Only used for reference calling
Sb        # Pval Binomial P-value for strand bias test
MQ        # Root mean square of mapping qualities of reads at the variant position
QD        # Variant-quality/read-depth for this variant
SC        # Genomic sequence 10 bases either side of variant position
BRF       # Fraction of reads around this variant that failed filters
HapScore  # The number of haplotypes supported in the calling window
Size      # Size of reference call block. Only used for reference calling


8.2 FILTER fields
=================

strandBias     # Variant fails strand-bias filter
alleleBias     # Variant fails allele-bias filter
badReads       # Variant is supported only by low-quality reads
Q20            # Variant-call has low posterior score ( < phred-score of 20)
GOF            # Variant-call fails goodness-of-fit test
PASS           # Variant passes all filters
QD             # Ratio of variant quality to number of supporting reads is low
SC             # Sequence context surrounding variant has low complexity

For multi-allelic sites, the filter field is based on information from the best variant
call at that site, so there may be variants which fail filters at these sites. VCF does not
allow allele-specific filter values.


8.3 FORMAT fields
=================

GT   # Un-phased genotype calls
GL   # Genotype log-likelihoods (natural log) for AA,AB and BB genotypes, where A = ref and B = variant. Only applicable for bi-allelic sites
GOF  # Phred-scaled goodness-of-fit score for the genotype call
GQ   # Phred-scaled quality score for the genotype call
NR   # Number of reads covering variant position in this sample
NV   # Number of reads at variant position which support the called variant in this sample


9 Known Issues
==============

When genotyping, if there are many variants close together, then not all will be reported in the output VCF. This is only true for variants
which are not supported by the data. Variants that are well supported will have genotype calls.


10 Release History
=================

0.1.5
===== 

Released in November 2011. First stable release of Platypus

0.1.6
===== 

Released on 23 Feb 2012. Various bug-fixes and tweaks, the most significant being
the following:

- Improved output when genotyping. FILTER and INFO fields are now correctly computed
- Fixed multi-sample genotyping, which previously had an underflow error in high-coverage regions
- Fixed small bug in posterior calculation, which was sometimes slightly wrong for high-confidence calls
- Improved data-filtering for duplicate reads and broken read pairs
- Added NV and NR fields to the per-sample FORMAT columns. NV is the number of variant-supporting reads, and NR is the total number of reads covering that position
- Fixed off-by-one error in VCF parsing, by adding pysam version check
- Added experimental flag --mergeClusteredVariants for calling in divergent samples. This is not fully tested, and does slow down calling substantially.
- Numerous other small fixes and performance tweaks

0.1.7
===== 

Released on 14 March 2012.

- Fixed error caused by zero coverage in window
- Improved read filtering: pairs with small inserts no longer lead to spurious insertion calls
- Improved haplotype generation: adjacent SNP/deletion haplotypes not generated correctly
- Switched to using single sample likelihoods by default, not EM likelihoods: this fixes various
  inconsistencies between the genotype calls and likelihoods
- Fixed bug when genotyping complex regions, which lead to incorrect genotypes being assigned
- Improved window algorithm when the --mergeClusteredVariants flag is set to 1
- Removed cap on likelihoods in alignment function, to improve calling in divergent regions
- Switched off EM filtering. Now uses simple coverage filter to remove noise.
- Fixed small bug in VCF INFO computations
- Fixed bug in left-normalisation. If Platypus cannot left-normalise an indel, it reports it in the position it was originally
  seen in the BAM.

0.1.8
===== 

Released on 19 March 2012.

- Various small bug-fixes
- Extended badReads filter to cover regions where many reads are filtered
- Improved Haplotype representation. Padded Haplotypes with Ns around window
- Improved strand bias and allele bias filtering

0.1.9
===== 

Released on 12 April 2012.

- Numerous small bug-fixes

- Platypus can now automatically handle samples spread across multiple BAMs. The BAMs are merged on the fly based on the
  value of the SM tag


0.2.0
===== 

Released on 25 July 2012

- First release to include (experimental) use of Cortex for local haplotype assembly


0.2.1
===== 

Released on 25 Jan 2013.

- Numerous small bug-fixes and minor improvements to filtering


0.2.3
===== 

- Introduced MNP and complex-replacement calling. Fixed various small bugs. Improved efficiency of EM algorithm.

0.2.4
=====

- Platypus now correctly splits/merges BAM files based on the read-group tag information. Multi-sample BAMs are now
  supported as well as single samples split over multiple BAMs.


0.4.0
=====

- Generally improved variant-calling and filtering
- Assembly with Cortex replaced with new implementation
- Improved memory usage

0.5.1
=====

- Added reference calling
- Fixed several small bugs in calling and filtering algorithms

0.5.2
=====

- Fixed bug when computing NR sample info field.

0.6.0
=====

- Released on 4 March 2014
- Includes per-sample variant-frequency filtering to reduce candidate list and run-time
- Fixed minor bug in MMLQ calculation


0.7.4
=====

- Released on 23 July 2014
- Major memory leak fixed when skipping regions due to high coverage.

0.8
=====

- Released on 3 March 2015
- CRAM support 
- Various bug fixes.
