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
ssw-py
~~~~~~~~~~

Python bindings for Complete-Striped-Smith-Waterman-Library (SSW) project

'''
__author__ = 'Nick Conway'
__copyright__ = (
    'Copyright 2023 Nick Conway; Copyright 2018, Nick Conway; Wyss Institute'
    'Harvard University'
)
__license__ = 'MIT'
__version__ = '1.0.0'
DESCRIPTION = (
    'Python bindings for Complete-Striped-Smith-Waterman-Library '
    '(SSW) project'
)


# `try` block here allows __version__, etc to be available prior to
# Cython extension building
from .alignmenttuple import Alignment

try:
    from .alignmentmgr import (  # type: ignore
        AlignmentMgr,
        BitwiseAlignmentFlag,
        force_align,
        format_force_align,
    )
    SSW = AlignmentMgr  # deprecated alias
except BaseException:
    pass
