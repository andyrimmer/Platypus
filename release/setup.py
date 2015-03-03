import os, sys, glob

from distutils.core import setup
from distutils.extension import Extension

cFlags = ["-msse2", "-msse3", "-funroll-loops", "-D_LARGEFILE64_SOURCE", "-D_FILE_OFFSET_BITS=64", "-fPIC"]
tabix_exclude = ( "main.c", )
tabix_dest = os.path.abspath( "tabix" )


extModules = []
extModules.append(Extension(name="ctabix", sources=["ctabix.c"] + ["tabix_util.c"] + glob.glob(os.path.join( "tabix", "*.pysam.c" )), include_dirs=["tabix", "pysam"], libraries=["z"], language="c"))
extModules.append(Extension(name="TabProxies", sources=[ "TabProxies.c", ], libraries=[ "z", ], language="c"))
extModules.append(Extension(name='htslibWrapper', sources=['htslibWrapper.c'], libraries=['hts'], extra_compile_args=cFlags))
extModules.append(Extension(name='fastafile', sources=['fastafile.c'], include_dirs=["samtools", "./"],  extra_compile_args=cFlags))
extModules.append(Extension(name='variant', sources=['variant.c'], include_dirs=["samtools", "./"], extra_compile_args=cFlags))
extModules.append(Extension(name='cerrormodel', sources=['cerrormodel.c', 'tandem.c'], include_dirs=["samtools", "./"], extra_compile_args=cFlags))
extModules.append(Extension(name='calign', sources=['calign.c', 'align.c'], include_dirs=["samtools", "./"], extra_compile_args=cFlags))
extModules.append(Extension(name='chaplotype', sources=['chaplotype.c', 'align.c'], include_dirs=["samtools", "./"], language='c', extra_compile_args=cFlags))
extModules.append(Extension(name='assembler', sources=['assembler.c'], include_dirs=["samtools", "./"], extra_compile_args=cFlags, language='c'))
extModules.append(Extension(name='platypusutils', sources=['platypusutils.c'],  extra_compile_args=cFlags))
extModules.append(Extension(name='cgenotype', sources=['cgenotype.c'],  extra_compile_args=cFlags))
extModules.append(Extension(name='vcfutils', sources=['vcfutils.c'],  extra_compile_args=cFlags))
extModules.append(Extension(name='cpopulation', sources=['cpopulation.c'],  extra_compile_args=cFlags))
extModules.append(Extension(name='cwindow', sources=['cwindow.c'],  extra_compile_args=cFlags))
extModules.append(Extension(name='vcfutils', sources=['vcfutils.c'],  extra_compile_args=cFlags))
extModules.append(Extension(name='variantFilter', sources=['variantFilter.c'],  extra_compile_args=cFlags))
extModules.append(Extension(name='variantcaller', sources=['variantcaller.c'], extra_compile_args=cFlags))

setup(name="vcf", py_modules=['vcf'])
setup(name="extendedoptparse", py_modules=['extendedoptparse'])
setup(name="window", py_modules=['window'])
setup(name="variantutils", py_modules=['variantutils'])
setup(name="platypusexceptions", py_modules=['platypusexceptions'])
setup(name="extendedoptparse", py_modules=['extendedoptparse'])
setup(name="filez", py_modules=['filez'])
setup(name="runner", py_modules=['runner'])
setup(name="PlatypusCythonModules", ext_modules=extModules)
