# ssw-py: SIMD Smith-Waterman Python Module for Use in Genomic Applications

This is derived from:
[original source repository](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library)

Please cite this [PLOS ONE paper](http://dx.plos.org/10.1371/journal.pone.0082138) by Zhao et al. 2013

## Overview

SSW is a fast implementation of the Smith-Waterman algorithm, which uses the
Single-Instruction Multiple-Data (SIMD) instructions to parallelize the
algorithm at the instruction level. It can return the Smith-Waterman score,
alignment location and traceback path (cigar) of the optimal alignment
accurately; and return the sub-optimal alignment score and location
heuristically.

Note: When SSW opens a gap, the gap open penalty alone is applied.

## ssw-py specific updates
`sse2neon.h` from the amazing https://github.com/DLTcollab/sse2neon project
was used to enable builds on arm 64 targets as a conditional macro target if
`__ARM_ARCH_ISA_A64` is defined

## Installation

from [PyPi](https://pypi.org/project/ssw-py/)

    pip install ssw-py


or from source

    python setup.py install
