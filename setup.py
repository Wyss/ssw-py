#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
test with:
    python setup.py build_ext --inplace
"""
DESCRIPTION = ("Complete Complete Striped Smith-Waterman Library")
LONG_DESCRIPTION = """
**sswpy** is a Python package
Cythonized wrapped version of 

https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library

Original authors should be cited if used.

License is MIT
"""

DISTNAME = 'sswpy'
LICENSE = 'MIT'
AUTHORS = "Nick Conway"
EMAIL = "nick.conway@wyss.harvard.edu"
URL = ""
DOWNLOAD_URL = ''
CLASSIFIERS = [
    'Development Status :: 1 - Beta',
    'Environment :: Console',
    'Intended Audience :: Science/Research',
    'Programming Language :: Python',
    'Programming Language :: Python :: 2',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3.3',
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Cython',
    'Topic :: Scientific/Engineering',
]

try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup, Extension

from Cython.Build import cythonize
import numpy.distutils.misc_util
import os
import sys
import shutil

pjoin = os.path.join
rpath = os.path.relpath

PACKAGE_PATH =      os.path.abspath(os.path.dirname(__file__))
MODULE_PATH =       pjoin(PACKAGE_PATH, 'ssw')
DATASETS_PATH =     pjoin(MODULE_PATH, 'datasets')

common_include = ['src']

if sys.platform == 'win32':
    extra_compile_args = ['']
else:
    extra_compile_args = ['-Wno-unused-function']

# fasta dataset files to include in installation
ssw_files = [rpath(pjoin(root, f), MODULE_PATH) for root, _, files in
                 os.walk(DATASETS_PATH) for f in files if '.fa' in f]

ssw_ext = Extension(
    'ssw.sswpy',
    depends=[],
    sources=['ssw/sswpy.pyx',
             'src/ssw.c',
             'src/str_util.c'],
    include_dirs=common_include + [numpy.get_include()],
    extra_compile_args=extra_compile_args
)
# ssw_files.append('sswpy.pxd')

is_py_3 = int(sys.version_info[0] > 2)

cython_extensions = [
    ssw_ext
]
cython_ext_list = cythonize(cython_extensions,
                            compile_time_env={'IS_PY_THREE': is_py_3})
setup(
    name=DISTNAME,
    maintainer=AUTHORS,
    packages=['ssw', 'ssw.sswpy'],
    ext_modules=cython_ext_list,
    include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs(),
    package_data={'ssw': ssw_files},
    maintainer_email=EMAIL,
    description=DESCRIPTION,
    license=LICENSE,
    url=URL,
    download_url=DOWNLOAD_URL,
    long_description=LONG_DESCRIPTION,
    classifiers=CLASSIFIERS,
    zip_safe=False
)
