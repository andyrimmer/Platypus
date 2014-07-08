
#PYTHON := python
#HEADERS := src/cython/cgenotype.pxd src/cython/cfilter.pxd
#SOURCES := src/cython/cpopulation.pyx src/cython/cgenotype.pyx src/cython/cwindow.pyx src/cython/cfilter.pyx

#platypus: ${HEADERS} ${SOURCES}
#	echo 'C Libraries will be installed to ' ${PLATYPUS}
#	cd src/cython; ${PYTHON} setup.py install --prefix=${PLATYPUS}
#	cd src/python; ${PYTHON} setup.py install --prefix=${PLATYPUS}


SAMTOOLSSRC := $(wildcard samtools/*.c)
SAMTOOLSOBJ := $(subst .c,.o,${SAMTOOLSSRC})

coreutils: ${SAMTOOLSOBJ}

${SAMTOOLSOBJ}: %.o: %.c
	${CC} -c $< -o $@

clean:
	${RM} ${SAMTOOLSOBJ}

#clean:
#	cd src/cython; rm -rf build; cpopulation.c cgenotype.c cwindow.c cfilter.c
#	echo ${PLATYPUS}
#	rm ${PLATYPUS}/lib/${PYTHON}/site-packages/cpopulation.so
#	rm ${PLATYPUS}/lib/${PYTHON}/site-packages/cgenotype.so
#	rm ${PLATYPUS}/lib/${PYTHON}/site-packages/cwindow.so
#	rm ${PLATYPUS}/lib/${PYTHON}/site-packages/cfilter.so
#	rm ${PLATYPUS}/lib/${PYTHON}/site-packages/platypusexceptions.py
#	rm ${PLATYPUS}/lib/${PYTHON}/site-packages/indel.py
#	rm ${PLATYPUS}/lib/${PYTHON}/site-packages/window.py
#	rm ${PLATYPUS}/lib/${PYTHON}/site-packages/variantcaller.py
#	rm ${PLATYPUS}/lib/${PYTHON}/site-packages/variantutils.py
