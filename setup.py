#!/usr/bin/env python
# Copyright (C) 2023 Nick Conway
# Copyright (C) 2014-2018, Nick Conway; Wyss Institute Harvard University
#
# The MIT License

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

# http://www.opensource.org/licenses/mit-license.php
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
    from setuptools import (
        Extension,
        setup,
    )
except ImportError:
    from distutils.core import (
        setup,
        Extension,
    )

import os
import os.path as op
import sys

PACKAGE_PATH = op.abspath(op.dirname(__file__))
MODULE_PATH = op.join(PACKAGE_PATH, 'ssw')
LIB_PATH = op.join(MODULE_PATH, 'lib')
TESTS_PATH = op.join(PACKAGE_PATH, 'tests')


COMMON_INCLUDE = ['ssw/lib/CSSWL/src', 'ssw/lib']

if sys.platform == 'win32':
    EXTRA_COMPILE_ARGS = ['']
else:
    EXTRA_COMPILE_ARGS = [
        '-Wno-error=declaration-after-statement',
        '-Wno-unused-function',
    ]

# source files to include in installation for tar.gz
SSWPY_FILES = [
    op.relpath(op.join(_root, _f), MODULE_PATH)
    for _root, _, _files in os.walk(LIB_PATH)
    for _f in _files
    if (
        ('.h' in _f) or
        ('.c' in _f) or
        ('.cpp' in _f) or
        ('.md' in _f) or
        ('.pyx' in _f) or
        ('.pxd' in _f) or
        ('Makefile' in _f)
    )
]

SSWPY_TEST_FPS = [
    op.relpath(op.join(root, fp), MODULE_PATH) for root, _, fps in
    os.walk(TESTS_PATH) for fp in fps
]

ALIGNMENTMGR_EXT = Extension(
    'ssw.alignmentmgr',
    sources=[
        op.join('ssw', 'alignmentmgr.pyx'),
        op.join('ssw', 'lib', 'CSSWL', 'src', 'ssw.c'),
        op.join('ssw', 'lib', 'str_util.c'),
    ],
    include_dirs=COMMON_INCLUDE,
    extra_compile_args=EXTRA_COMPILE_ARGS,
)


CLASSIFIERS = [
    'Development Status :: 5 - Production/Stable',
    'Environment :: Console',
    'Intended Audience :: Science/Research',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3',
    'Programming Language :: Cython',
    'Topic :: Scientific/Engineering',
]


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
        # embedsignature=True,
        binding=True,
    )

    return cythonize(
        extension_list,
        compiler_directives=cython_compiler_directives,
    )


import ssw

setup(
    name='ssw-py',
    version=ssw.__version__,
    author=ssw.__author__,
    author_email='a.grinner@gmail.com',
    url='https://github.com/libnano/ssw-py',
    packages=['ssw'],
    ext_modules=try_cythonize([ALIGNMENTMGR_EXT]),
    package_data={'ssw': SSWPY_FILES + SSWPY_TEST_FPS},
    description=ssw.DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    license=ssw.__license__,
    classifiers=CLASSIFIERS,
    setup_requires=['Cython', 'setuptools>=65.6.3'],
    zip_safe=False,
)
