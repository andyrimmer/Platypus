PYTHON := python
HEADERS := src/cython/cgenotype.pxd
SOURCES := src/cython/cpopulation.pyx src/cython/cgenotype.pyx src/cython/cwindow.pyx

PLATYPUS_C := cpopulation.c cgenotype.c cwindow.c vcfutils.c platypusutils.c variantcaller.c variantFilter.c assembler.c\
calign.c cerrormodel.c chaplotype.c samtoolsWrapper.c variant.c fastafile.c

PLATYPUS_SO := $(PLATYPUS_C:.c=.so)

# Python source files needed for building release
PYTHONSRC := src/python/Platypus.py src/python/platypusexceptions.py src/python/runner.py src/python/variantutils.py\
src/python/window.py src/python/vcf.py src/python/filez.py src/python/extendedoptparse.py

VERSION := 0.7.9.3
TARGET := Platypus_${VERSION}

# Cython source files needed for building release (these are not compiled, only distributed)
CYTHONSRC:= src/cython/platypusutils.pyx src/cython/cgenotype.pyx src/cython/cpopulation.pyx\
src/cython/cwindow.pyx src/cython/variantcaller.pyx src/cython/vcfutils.pyx src/cython/variantFilter.pyx\
src/cython/assembler.pyx src/cython/calign.pyx src/cython/chaplotype.pyx\
src/cython/fastafile.pyx src/cython/samtoolsWrapper.pyx src/cython/variant.pyx src/cython/cerrormodel.pyx\
src/pysam/ctabix.pyx src/pysam/ctabix.pxd src/pysam/TabProxies.pxd src/pysam/TabProxies.pyx

# C source files needed for building release (these will be compiled when building the released code)
CSRC := src/cython/platypusutils.c src/cython/cgenotype.c src/cython/cpopulation.c\
src/cython/cwindow.c src/cython/variantcaller.c src/cython/vcfutils.c src/cython/variantFilter.c\
src/cython/assembler.c src/c/align.c src/cython/calign.c src/cython/chaplotype.c src/cython/fastafile.c\
src/c/pysam_util.c src/cython/samtoolsWrapper.c src/cython/variant.c src/c/pysam_util.h src/c/align.h\
src/cython/cerrormodel.c src/c/tandem.h src/c/tandem.c src/pysam/ctabix.c src/pysam/TabProxies.c src/pysam/tabix_util.c

OTHER := misc/README.txt LICENSE release/setup.py release/buildPlatypus.sh

platypus: ${HEADERS} ${SOURCES}
	echo 'Building Platypus'
	cd src; ${PYTHON} setup.py build
	mkdir -p bin
	cp src/build/*/*.so bin/
	cp src/build/*/python/*.py bin/

clean:
	rm -rf bin
	cd src; rm -rf build; cd cython; rm -f ${PLATYPUS_C}

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
	rsync -a --exclude='.svn' src/samtools ${TARGET}/
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
