# Copyright (C) 2023, Nick Conway
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
ssw.alignmenttuple
~~~~~~~~~~~~~~~~~~

Alignment data structure for ssw-py

NOTE: See this link for a info on the CIGAR format:
`What is a CIGAR? <http://genome.sph.umich.edu/wiki/SAM#What_is_a_CIGAR.3F>`_

'''

from typing import (
    List,
    NamedTuple,
    Tuple,
)


class Alignment(NamedTuple):
    '''Alignment namedtuple structure.  Used for storing results of
    :meth:`ssw.alignmentmgr.AlignmentMgr.align`

    Attributes:
        CIGAR: The CIGAR string
        optimal_score: The optimal alignment score (SSW primary score)
        sub_optimal_score: The sub-optimal alignment score (SSW secondary score)
        reference_start: The reference start index
        reference_end: The reference end index
        read_start: The read start index
        read_end: The read end index
        cigar_pair_list: Convenience list of cigar length code pairs

    '''
    CIGAR: str
    optimal_score: int
    sub_optimal_score: int
    reference_start: int
    reference_end: int
    read_start: int
    read_end: int
    cigar_pair_list: List[Tuple[int, str]]

    def len_cigar(self) -> int:
        '''
        Returns:
            Number of CIGAR pairs
        '''
        return len(self.cigar_pair_list)

    def get_cigar_values(self, idx: int) -> Tuple[int, str]:
        '''Fetch CIGAR value at an index in the CIGAR string

        Args:
            idx: Index to fetch

        Returns:
            Tuple of the form:

                <CIGAR length>, <CIGAR code>

        '''
        return self.cigar_pair_list[idx]

    def pprint(self):
        print(
            'Alignment(\n'
            f'  CIGAR={self.CIGAR}\n'
            f'  optimal_score={self.optimal_score}\n'
            f'  sub_optimal_score={self.sub_optimal_score}\n'
            f'  reference_start={self.reference_start}\n'
            f'  reference_end={self.reference_end}\n'
            f'  read_start={self.read_start}\n'
            f'  read_end={self.read_end}\n'
            f'  cigar_pair_list={self.cigar_pair_list}\n'
            ')',
        )
