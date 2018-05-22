#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
test with:
    python setup.py build_ext --inplace

build wheels with:
    python setup.py bdist_wheel --plat-name win_amd64 --python-tag cp36
    python setup.py bdist_wheel --plat-name macosx_10_10_x86_64 --python-tag cp36
    python setup.py sdist --formats=gztar
    twine upload dist/*

https://pypi.python.org/pypi/ssw-py
'''
DESCRIPTION = ("Complete Striped Smith-Waterman Library for Python")
LONG_DESCRIPTION = '''
**ssw-py** is a Python package
Cythonized wrapped version of:

https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library

Original C library authors and paper should be cited if used.

Support for Python 3 at the moment

License is MIT
'''

try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup, Extension

import numpy.distutils.misc_util
import os
import sys
import shutil
import re
import ast

pjoin = os.path.join
rpath = os.path.relpath

PACKAGE_PATH =      os.path.abspath(os.path.dirname(__file__))
MODULE_PATH =       pjoin(PACKAGE_PATH, 'ssw')

common_include = ['lib/CSSWL/src', 'lib']

if sys.platform == 'win32':
    extra_compile_args = ['']
else:
    extra_compile_args = ['-Wno-unused-function']

# fasta dataset files to include in installation
# ssw_files = [rpath(pjoin(root, f), MODULE_PATH) for root, _, files in
#                  os.walk(DATASETS_PATH) for f in files if '.fa' in f]
ssw_files = []

ssw_ext = Extension(
    'ssw.sswpy',
    sources=['ssw/sswpy.pyx',
             'lib/CSSWL/src/ssw.c',
             'lib/str_util.c'],
    include_dirs=common_include + [numpy.get_include()],
    extra_compile_args=extra_compile_args
)

# Initialize subtree with
# git subtree add --prefix=lib/CSSWL git@github.com:mengyao/Complete-Striped-Smith-Waterman-Library.git master

# Begin modified code from Flask's version getter
# BSD license
# Copyright (c) 2015 by Armin Ronacher and contributors.
# https://github.com/pallets/flask
_version_re = re.compile(r'__version__\s+=\s+(.*)')

with open('ssw/__init__.py', 'rb') as initfile:
    VERSION = str(ast.literal_eval(_version_re.search(
                                   initfile.read().decode('utf-8')).group(1)))

DISTNAME = 'ssw-py'
LICENSE = 'MIT'
AUTHORS = "Nick Conway"
EMAIL = "nick.conway@wyss.harvard.edu"
CLASSIFIERS = [
    'Development Status :: 4 - Beta',
    'Environment :: Console',
    'Intended Audience :: Science/Research',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Cython',
    'Topic :: Scientific/Engineering',
]

setup(
    name=DISTNAME,
    version=VERSION,
    author=AUTHORS,
    author_email=EMAIL,
    url='https://github.com/Wyss/ssw-py',
    packages=['ssw'],
    ext_modules=[ssw_ext],
    include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs(),
    package_data={'ssw': ssw_files},
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    license=LICENSE,
    classifiers=CLASSIFIERS,
    zip_safe=False
)
