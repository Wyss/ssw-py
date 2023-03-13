# Copyright 2023 Nick Conway;
#
# See LICENSE.md for full MIT license.
'''
example.blast

Example derivation of upstream C libraries ``example.c``

Calls Python ``sswpy_write`` to display an :meth:`ssw.alignmenttuple.Alignment`
result

Also call :meth:`ssw.alignmentmgr.AlignmentMgr.print_result` code to print out
same result as as a confirmation
'''
import os.path as op
import sys
from typing import List

PACKAGE_DIR = op.dirname(op.dirname(op.abspath(__file__)))
sys.path.append(PACKAGE_DIR)

from ssw import (
    Alignment,
    AlignmentMgr,
)


def sswpy_write(
        alignment: Alignment,
        ref_seq: str,
        read_seq: str,
):
    '''Print a BLAST-like output

    Args:
        alignment: :class:`ssw.alignmenttuple.Alignment` named tuple
        ref_seq: Reference sequence python string or bytes-string
        read_seq: Read sequence python string or bytes-string

    '''
    print(f'CIGAR start index {alignment.reference_start}: {alignment.CIGAR}')
    print(
        f'optimal_alignment_score: {alignment.optimal_score}\n'
        f'sub-optimal_alignment_score: {alignment.sub_optimal_score}\n',
        end='',
    )
    if alignment.reference_start + 1:
        print(f'target_begin: {alignment.reference_start}')

    print(f'target_end: {alignment.reference_end}')

    if alignment.read_start + 1:
        print(f'query_begin: {alignment.read_start}')

    print(f'query_end: {alignment.read_end}\n')

    if alignment.CIGAR:
        c = 0
        left = 0
        e = 0
        qb = alignment.reference_start
        pb = alignment.read_start
        cigar_len = len(alignment.cigar_pair_list)
        break_to_next_step = False
        while (e < cigar_len) or (left > 0):
            count = 0
            q = qb
            p = pb
            print('Target: %8d    ' % (q), end='')
            for c in range(e, cigar_len):
                clength, letter = alignment.cigar_pair_list[c]
                lrange = left if ((count == 0) and (left > 0)) else clength
                for i in range(lrange):
                    if letter == 'I':
                        print('-', end='')
                    else:
                        print(ref_seq[q], end='')
                        q += 1
                    count += 1
                    if count == 60:
                        # goto step2
                        # Need to break out of the for loop to get to step2
                        break_to_next_step = True
                        break
                if break_to_next_step:
                    # Reset break flag and goto to step2
                    break_to_next_step = False
                    break
# step2:
            print(f'    {q - 1}\n                    ', end='')
            q = qb
            count = 0
            for c in range(e, cigar_len):
                clength, letter = alignment.cigar_pair_list[c]
                lrange = left if ((count == 0) and (left > 0)) else clength
                for i in range(lrange):
                    if letter == 'M':
                        if ref_seq[q] == read_seq[p]:
                            print('|', end='')
                        else:
                            print('*', end='')
                        q += 1
                        p += 1
                    else:
                        print('*', end='')
                        if letter == 'I':
                            p += 1
                        else:
                            q += 1
                    if count == 60:
                        qb = q
                        # goto step3
                        # Need to break out of the for loop to get to step3
                        break_to_next_step = True
                        break

                if break_to_next_step:
                    # Reset break flag and goto to step3
                    break_to_next_step = False
                    break
# step3:
            p = pb
            print('\nQuery:  %8d    ' % (p + 1), end='')
            count = 0
            for c in range(e, cigar_len):
                clength, letter = alignment.cigar_pair_list[c]
                lrange = left if ((count == 0) and (left > 0)) else clength
                for i in range(lrange):
                    if letter == 'D':
                        print('-', end='')
                    else:
                        print(read_seq[p], end='')
                        p += 1
                    count += 1
                    if count == 60:
                        pb = p
                        left = lrange - i - 1
                        e = (c + 1) if (left == 0) else c
                        # goto end_step
                        # Need to break out of the while loop to get to end_step
                        break_to_next_step = True
                        break
                if break_to_next_step:
                    # Need to break out of the while loop to get to end_step
                    break
            c += 1
            if not break_to_next_step:
                e = c
                left = 0
# end_step:
            # Reset break flag and goto to end_step
            break_to_next_step = False
            print(f'    {p}\n\n', end='')
        # end while


# Align a pair of genome sequences.
def main():
    gap_open = 3
    gap_extension = 1  # default parameters for genome sequence alignment

    # Reference sequence
    # Should align to indices 8 to 21 inclusive
    ref_seq = ''.join([
        'C', 'A', 'G', 'C', 'C', 'T', 'T', 'T',
        'C', 'T', 'G', 'A', 'C', 'C', 'C', 'G',
        'G', 'A', 'A', 'A', 'T', 'C', 'A', 'A',
        'A', 'A', 'T', 'A', 'G', 'G', 'C', 'A',
        'C', 'A', 'A', 'C', 'A', 'A', 'A',
    ])

    # Read sequence
    read_seq = ''.join([
        'C', 'T', 'G', 'A', 'G', 'C', 'C', 'G',
        'G', 'T', 'A', 'A', 'A', 'T', 'C',
    ])

    align_mgr = AlignmentMgr(
        match_score=2,
        mismatch_penalty=2,
    )
    align_mgr.set_read(read_seq)
    align_mgr.set_reference(ref_seq)

    alignment = align_mgr.align(
        gap_open=gap_open,
        gap_extension=gap_extension,
    )
    print('> 1. sswpy_write() call\n')
    sswpy_write(
        alignment,
        ref_seq,
        read_seq,
    )

    # Confirm with print_result call
    print('\n> 2. AlignmentMgr.print_result()')
    align_mgr.print_result(
        alignment,
    )


if __name__ == '__main__':
    main()
