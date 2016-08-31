Platypus
========

The Platypus variant caller.

Installation and execution
--------------------------

To build Platypus, do the following:

    make

Then to run do

    python bin/Platypus.py --bamFiles=BAM.bam --refFile=REF.fa --output=variants.vcf

Platypus has been tested with Python 2.6 and 2.7, and requires Cython 0.20.2 or later
to build.

### Prerequisites ###

Platypus requires htslib 1.2.1 or greater. htslib can be downloaded from the htslib website http://www.htslib.org/download/.

To build and install htslib, cd into htslib source and type `make install`. This will install htslib under `/usr/local/` (see note below). To install htslib in any other directory use `make install prefix=/path/to/dir`.

NOTE: htslib should be installed in a standard location (e.g. `/usr/local/`). If not installed in a standard location, you will need to set your library paths:

For linux:

    export C_INCLUDE_PATH=/path/to/dir/include
    export LIBRARY_PATH=/path/to/dir/lib (only for making)
    export LD_LIBRARY_PATH=/path/to/dir/lib

Note the `/include` and `/lib` sub-directories. e.g. if you installed htslib under `/Users/me/htslib` then set

    export C_INCLUDE_PATH=/Users/me/htslib/include
    export LIBRARY_PATH=/Users/me/htslib/lib
    export LD_LIBRARY_PATH=/Users/me/htslib/lib

htslib will automatically make the `include` and `lib` directories on install.

For OSX:

    export C_INCLUDE_PATH=/path/to/dir/include
    export LIBRARY_PATH=/path/to/dir/lib (only for making)
    export DYLD_FALLBACK_LIBRARY_PATH=/path/to/dir/lib
