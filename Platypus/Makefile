
PYTHON := python
HEADERS := src/cython/cgenotype.pxd
SOURCES := src/cython/cpopulation.pyx src/cython/cgenotype.pyx src/cython/cwindow.pyx
COREUTILS_C := calign.c cerrormodel.c chaplotype.c histogram.c samtoolsWrapper.c variant.c fastafile.c
PLATYPUS_C := cpopulation.c cgenotype.c cwindow.c vcfutils.c platypusutils.c variantcaller.c variantFilter.c assembler.c
COREUTILS_SO := $(COREUTILS_C:.c=.so)
PLATYPUS_SO := $(PLATYPUS_C:.c=.so)

# Python source files needed for building release
PYTHONSRC := src/python/Platypus.py src/python/platypusexceptions.py src/python/runner.py src/python/variantutils.py\
src/python/window.py src/python/vcf.py ../coreutils/filez.py src/python/extendedoptparse.py

VERSION := 0.7.5
TARGET := Platypus_${VERSION}

# C source files needed for building release (these are not compiled, only distributed)
CYTHONSRC:= src/cython/platypusutils.pyx src/cython/cgenotype.pyx src/cython/cpopulation.pyx\
src/cython/cwindow.pyx src/cython/variantcaller.pyx src/cython/vcfutils.pyx src/cython/variantFilter.pyx\
src/cython/assembler.pyx ../coreutils/calign.pyx ../coreutils/chaplotype.pyx\
../coreutils/fastafile.pyx ../coreutils/samtoolsWrapper.pyx ../coreutils/variant.pyx ../coreutils/cerrormodel.pyx\
src/pysam/ctabix.pyx src/pysam/ctabix.pxd src/pysam/TabProxies.pxd src/pysam/TabProxies.pyx

# C source files needed for building release (these will be compiled when building the released code)
CSRC := src/cython/platypusutils.c src/cython/cgenotype.c src/cython/cpopulation.c\
src/cython/cwindow.c src/cython/variantcaller.c src/cython/vcfutils.c src/cython/variantFilter.c\
src/cython/assembler.c ../coreutils/align.c ../coreutils/calign.c ../coreutils/chaplotype.c ../coreutils/fastafile.c\
../coreutils/pysam_util.c ../coreutils/samtoolsWrapper.c ../coreutils/variant.c ../coreutils/pysam_util.h ../coreutils/align.h\
../coreutils/cerrormodel.c ../coreutils/tandem.h ../coreutils/tandem.c src/pysam/ctabix.c src/pysam/TabProxies.c src/pysam/tabix_util.c

OTHER := misc/README.txt misc/LICENSE release/setup.py release/buildPlatypus.sh

platypus: ${HEADERS} ${SOURCES}
	echo 'Building Platypus'
	cd ../coreutils; ${PYTHON} setup.py build 
	cd src; ${PYTHON} setup.py build
	mkdir -p bin
	cp ../coreutils/build/*/*.so bin/
	cp ../coreutils/build/*/*.py bin/
	cp src/build/*/*.so bin/
	cp src/build/*/python/*.py bin/

clean:
	rm -rf bin
	cd src; rm -rf build; cd cython; rm -f ${PLATYPUS_C}
	cd ../coreutils; rm -rf build; rm ${COREUTILS_C}

.PHONY:
releasedir: platypus
	echo ''
	echo 'Building new Platypus release'
	rm -rf ${TARGET}
	mkdir ${TARGET}
	mkdir ${TARGET}/src
	mkdir ${TARGET}/src/pysam
	echo 'Copying Platypus C source files to ' ${TARGET}
	cp -f ${CSRC} ${TARGET}
	echo 'Copying Platypus Python source files to ' ${TARGET}
	cp -f ${PYTHONSRC} ${TARGET}
	echo 'Copying Platypus Cython source files to ' ${TARGET}/src
	cp -f ${CYTHONSRC} ${TARGET}/src
	echo 'Copying Platypus README, LICENCE etc files to ' ${TARGET}
	cp -f ${OTHER} ${TARGET}
	echo 'Copying Samtools files'
	rsync -a --exclude='.svn' ../coreutils/samtools ${TARGET}/
	echo 'Copying pysam *.c files'
	cp -f src/pysam/*.c ${TARGET}/src/pysam/
	cp -f src/pysam/COPYING ${TARGET}/src/pysam/
	echo 'Copying tabix files'
	rsync -a --exclude='.svn' src/tabix ${TARGET}/
	echo 'Cleaning up other unwanted files'

.PHONY:
tarball: releasedir
	tar cvzf ${TARGET}.tgz ${TARGET}
	rm -rf ${TARGET}

.PHONY:
release: tarball
	echo 'Finished building release for Platypus_${VERSION}'


.PHONY: test
test: install
	python sanityChecks.py test.vcf
