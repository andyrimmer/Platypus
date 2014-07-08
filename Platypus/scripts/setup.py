import os, sys, glob

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext as build_pyx

extModules = []
extModules.append(Extension(name='palindrome', sources=['palindrome.pyx',]))
setup(name = "Palindrome", ext_modules=extModules, cmdclass={'build_ext': build_pyx})
