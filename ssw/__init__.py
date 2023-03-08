# Copyright (C) 2014-2023. Ben Pruitt & Nick Conway; 2014-2018 Wyss Institute
# See LICENSE.md for full MIT license.
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
__version__ = '1.0.0a2'
DESCRIPTION = (
    'Python bindings for Complete-Striped-Smith-Waterman-Library '
    '(SSW) project'
)


# `try` block here allows __version__, etc to be available prior to
# Cython extension building
from .alignment import Alignment

try:
    from .alignmentmgr import (  # type: ignore
        AlignmentMgr,
        force_align,
        format_force_align,
    )
    SSW = AlignmentMgr  # deprecated alias
except BaseException:
    pass
