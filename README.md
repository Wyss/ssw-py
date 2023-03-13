# ssw-py: Striped Smith-Waterman SIMD accelerated Python Package for Use in Genomic Applications

<a href="https://github.com/libnano/ssw-py/actions/" rel="actions">![Actions](https://github.com/libnano/ssw-py/actions/workflows/ssw-py-ci-github-action.yml/badge.svg)</a>
<a href="http://www.gnu.org/licenses/gpl-2.0.html" rel="license">![License](https://img.shields.io/pypi/l/ssw-py.png)</a>
<a href="https://pypi.python.org/pypi/ssw-py" rel="pypi">![PyPi](https://img.shields.io/pypi/v/ssw-py.png)</a>


This library uses the excellent source code from this is
[original source repository](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library)

Please cite this [PLOS ONE paper](http://dx.plos.org/10.1371/journal.pone.0082138) by Zhao et al. 2013 when using this software.

## Overview

**ssw-py** provides a fast implementation of the Smith-Waterman algorithm,
which uses the Single-Instruction Multiple-Data (SIMD) instructions to parallelize
the algorithm at the instruction level.

Using `ssw.AlignmentMgr`, you can compute the Smith-Waterman score, alignment location and traceback path
([CIGAR](https://genome.sph.umich.edu/wiki/SAM#What_is_a_CIGAR.3F)) of the
optimal alignment accurately; and return the sub-optimal alignment score and
location heuristically.

Note: When Striped Smith-Waterman opens a gap, the gap open penalty alone is applied.

## Installation

from [PyPi](https://pypi.org/project/ssw-py/)

    $ pip install ssw-py


or from source

    $ python setup.py install

## Documentation
See [documentation](https://libnano.github.io/ssw-py/) for help on using these
bindings.
