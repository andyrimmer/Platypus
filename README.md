
The Platypus variant caller. To build Platypus, do the following:

make

Then to run do

python bin/Platypus.py --bamFiles=BAM.bam --refFile=REF.fa --output=variants.vcf

Platypus has been tested with Python 2.6 and 2.7, and requires Cython 0.20.1 or later
to build.
