#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
=========================================================================
ssw-py: SIMD Smith-Waterman Python Module for Use in Genomic Applications
=========================================================================

Installation
------------

**ssw-py** has no external runtime library dependencies and should compile
on easily most Linux and MacOS systems

https://pypi.python.org/pypi/ssw-py

To build **ssw-py** within the package directory run::

  $ python setup.py build_ext --inplace

If you would like to install ssw-py in your local Python environment
you may do so using either pip or the setup.py script::

  $ pip install ssw-py
            or
  $ python setup.py install

Create wheels with wheel installed and tar.gz::

  $ python setup.py bdist_wheel
  $ python setup.py sdist --formats=gztar

upload wheels with:
    twine upload dist/*

'''

try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup, Extension

import os
import os.path as op
import sys

PACKAGE_PATH =      op.abspath(op.dirname(__file__))
MODULE_PATH =       op.join(PACKAGE_PATH, 'ssw')
LIB_PATH =          op.join(MODULE_PATH, 'lib')

COMMON_INCLUDE = ['ssw/lib/CSSWL/src', 'ssw/lib']

if sys.platform == 'win32':
    EXTRA_COMPILE_ARGS = ['']
else:
    EXTRA_COMPILE_ARGS = [
        '-Wno-error=declaration-after-statement',
        '-Wno-unused-function',
    ]

# source files to include in installation for tar.gz
SSW_FILES = [
    op.relpath(op.join(_root, _f), MODULE_PATH) 
    for _root, _, _files in os.walk(LIB_PATH) 
    for _f in _files 
    if (
        ('.h' in _f) or ('.c' in _f) or ('.cpp' in _f)
    )
]

ssw_ext = Extension(
    'ssw.sswpy',
    sources=[
        'ssw/sswpy.pyx',
        'ssw/lib/CSSWL/src/ssw.c',
        'ssw/lib/str_util.c'
    ],
    include_dirs=COMMON_INCLUDE,
    extra_compile_args=EXTRA_COMPILE_ARGS
)

# Initialize subtree with
# git subtree add --prefix=lib/CSSWL git@github.com:mengyao/Complete-Striped-Smith-Waterman-Library.git master


CLASSIFIERS = [
    'Development Status :: 5 - Production/Stable',
    'Environment :: Console',
    'Intended Audience :: Science/Research',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3',
    'Programming Language :: Cython',
    'Topic :: Scientific/Engineering',
]

import ssw

with open('README.md') as fd:
    LONG_DESCRIPTION = fd.read()

def try_cythonize(extension_list, *args, **kwargs):
    '''
    Light cythonize wrapper
    '''
    try:
        from Cython.Build import cythonize
    except (ImportError, ModuleNotFoundError):
        def cythonize(x, **kwargs):
            return x

    cython_compiler_directives = dict(
        language_level='3',
        c_string_encoding='utf-8',
    )

    return cythonize(
        extension_list,
        compiler_directives=cython_compiler_directives,
    )


setup(
    name='ssw-py',
    version=ssw.__version__,
    author=ssw.__author__,
    author_email='a.grinner@gmail.com',
    url='https://github.com/Wyss/ssw-py',
    packages=['ssw'],
    ext_modules=try_cythonize([ssw_ext]),
    package_data={'ssw': SSW_FILES},
    description=ssw.DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    license=ssw.__license__,
    classifiers=CLASSIFIERS,
    setup_requires=['Cython', 'setuptools>=65.6.3'],
    zip_safe=False
)
