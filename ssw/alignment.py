
# Copyright 2023 Nick Conway;
#
# See LICENSE.md for full MIT license.
'''
ssw.alignment
~~~~~~~~~~~~~~

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
    '''Alignment namedtuple structure

    Attributes:
        CIGAR: The CIGAR string
        cigar_pair_list: list of cigar pairs
        optimal_score: The optimal alignment score
        sub_optimal_score: The sub-optimal alignment score
        reference_start: The reference start index
        reference_end: The reference end index
        read_start: The read start index
        read_end: The read end index

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
