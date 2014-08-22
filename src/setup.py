import os, sys, glob

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext as build_pyx

#
# Build Python modules first
#

setup(name="window", py_modules=['python/window'])
setup(name="filez", py_modules=['python/filez'])
setup(name="variantutils", py_modules=['python/variantutils'])
setup(name="platypusexceptions", py_modules=['python/platypusexceptions'])
setup(name="extendedoptparse", py_modules=['python/extendedoptparse'])
setup(name="vcf", py_modules=['python/vcf'])
setup(name="runner", py_modules=['python/runner'])
setup(name="Platypus", py_modules=['python/Platypus'])


#
# Then deal with Cython extensions
#

extModules = []
corMods = ['cython/chaplotype.pxd', 'cython/variant.pxd', 'cython/fastafile.pxd', 'cython/calign.pxd', 'cython/samtoolsWrapper.pxd']
incDirs = ["samtools", "./", "c"]

# Debug for Valgrind
cFlags = ["-msse2", "-msse3", "-funroll-loops", "-D_LARGEFILE64_SOURCE", "-D_FILE_OFFSET_BITS=64" ,"-g"]

# No debug. I don't think this affects performance?
#cFlags = ["-msse2", "-msse3", "-funroll-loops", "-D_LARGEFILE64_SOURCE", "-D_FILE_OFFSET_BITS=64"]

extModules.append(Extension(name='samtoolsWrapper', sources=['cython/samtoolsWrapper.pyx', 'c/pysam_util.c'] + glob.glob(os.path.join("samtools", "*.c")), include_dirs=incDirs,libraries=['z'], language='c', extra_compile_args=cFlags))
extModules.append(Extension(name='arrays', sources=['cython/arrays.pyx'], include_dirs=incDirs, language='c', extra_compile_args=cFlags))
extModules.append(Extension(name='fastafile', sources=['cython/fastafile.pyx'], include_dirs=incDirs,  extra_compile_args=cFlags))
extModules.append(Extension(name='variant', sources=['cython/variant.pyx', 'cython/samtoolsWrapper.pxd'], include_dirs=incDirs, extra_compile_args=cFlags))
extModules.append(Extension(name='cerrormodel', sources=['cython/cerrormodel.pyx', 'c/tandem.c'], include_dirs=incDirs, language='c', extra_compile_args=cFlags))
extModules.append(Extension(name='calign', sources=['cython/calign.pyx', 'c/align.c', 'cython/samtoolsWrapper.pxd', 'cython/cerrormodel.pxd'], include_dirs=incDirs, extra_compile_args=cFlags))
extModules.append(Extension(name='chaplotype', sources=['cython/chaplotype.pyx', 'cython/variant.pxd','cython/samtoolsWrapper.pxd', 'cython/calign.pxd', 'c/align.c'], include_dirs=incDirs, language='c', extra_compile_args=cFlags))

extModules.append(Extension(name="ctabix", sources=["pysam/ctabix.pyx"] + ["pysam/tabix_util.c"] + glob.glob("tabix/*.pysam.c"), include_dirs=["tabix", "pysam"], libraries=["z"], language="c"))
extModules.append(Extension(name="TabProxies", sources=[ "pysam/TabProxies.pyx", ], libraries=[ "z", ], language="c"))
extModules.append(Extension(name='cwindow', sources=corMods+['cython/cwindow.pyx'], include_dirs=incDirs,  extra_compile_args=cFlags, language='c'))
extModules.append(Extension(name='assembler', sources=corMods+['cython/assembler.pyx','cython/cwindow.pxd'], include_dirs=incDirs, extra_compile_args=cFlags, language='c'))
extModules.append(Extension(name='cgenotype', sources=corMods+['cython/cgenotype.pyx','cython/cwindow.pxd'], include_dirs=incDirs,  extra_compile_args=cFlags, language='c'))
extModules.append(Extension(name='platypusutils', sources=corMods+['cython/platypusutils.pyx','cython/cwindow.pxd'], include_dirs=incDirs,  extra_compile_args=cFlags, language='c'))
extModules.append(Extension(name='vcfutils', sources=corMods+['cython/vcfutils.pyx', 'cython/platypusutils.pxd', 'cython/cwindow.pxd'], include_dirs=incDirs,  extra_compile_args=cFlags, language='c'))
extModules.append(Extension(name='cpopulation', sources=corMods+['cython/cpopulation.pyx', 'cython/platypusutils.pxd', 'cython/vcfutils.pxd', 'cython/cwindow.pxd'],  include_dirs=incDirs, extra_compile_args=cFlags, language='c'))
extModules.append(Extension(name='variantFilter', sources=corMods+['cython/variantFilter.pyx', 'cython/cpopulation.pxd', 'cython/cwindow.pxd', 'cython/platypusutils.pxd'], include_dirs=incDirs,  extra_compile_args=cFlags, language='c'))
extModules.append(Extension(name='variantcaller', sources=corMods+['cython/variantcaller.pyx', 'cython/cpopulation.pxd', 'cython/platypusutils.pxd', 'cython/vcfutils.pxd', 'cython/cwindow.pxd'], include_dirs=incDirs,  extra_compile_args=cFlags, language='c'))

# Set Cython directives for each module
for module in extModules:
    module.cython_directives = {"boundscheck": False, "nonecheck" : False, "cdivision" : True, "profile" : False, "initializedcheck" : False, "wraparound" : True}

setup(name = "PlatypusCythonModules", ext_modules=extModules, cmdclass = {'build_ext': build_pyx})
