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
tests.test_ssw
~~~~~~~~~~~~~~

Tests for the ssw-py
'''
import unittest

from ssw import (  # type: ignore
    AlignmentMgr,
    BitwiseAlignmentFlag,
    force_align,
    format_force_align,
)

# Reference sequence
REF_SEQ = ''.join([
    'C', 'A', 'G', 'C', 'C', 'T', 'T', 'T',
    'C', 'T', 'G', 'A', 'C', 'C', 'C', 'G',
    'G', 'A', 'A', 'A', 'T', 'C', 'A', 'A',
    'A', 'A', 'T', 'A', 'G', 'G', 'C', 'A',
    'C', 'A', 'A', 'C', 'A', 'A', 'A',
])

# Read sequence
READ_SEQ = ''.join([
    'C', 'T', 'G', 'A', 'G', 'C', 'C', 'G',
    'G', 'T', 'A', 'A', 'A', 'T', 'C',
])


class TestAlignmentMgr(unittest.TestCase):
    '''Tests for the AlignmentMgr class
    '''

    def setUp(self):
        self.align_mgr = AlignmentMgr(
            match_score=2,
            mismatch_penalty=2,
        )

    def test_deprecation(self):
        '''Test camelCase Deprecation'''
        with self.assertWarns(UserWarning):
            self.align_mgr.setRead(b'ACGT')

        with self.assertWarns(UserWarning):
            self.align_mgr.setReference(b'ACGT')

    def test_exact(self):
        a = self.align_mgr
        a.set_read(b'ACGT')
        ref = b'TTTTACGTCCCCC'
        a.set_reference(b'TTTTACGTCCCCC')
        res = a.align()

        self.assertEqual(res.optimal_score, 8)
        self.assertEqual(res.sub_optimal_score, 0)

        self.assertEqual(res.reference_start, 4)
        self.assertEqual(res.reference_end, 7)
        self.assertEqual(res.read_start, 0)
        self.assertEqual(res.read_end, 3)

        # a.print_result(res)

        with self.assertWarns(UserWarning):
            a.setReference(ref)

    def test_deletion(self):
        a = self.align_mgr
        a.set_read(b'ACGT')
        a.set_reference(b'TTTTACAGTCCCCC')
        res = a.align()

        self.assertEqual(res.optimal_score, 5)
        self.assertEqual(res.sub_optimal_score, 0)

        self.assertEqual(res.reference_start, 4)
        self.assertEqual(res.reference_end, 8)
        self.assertEqual(res.read_start, 0)
        self.assertEqual(res.read_end, 3)

        # a.print_result(res)

    def test_insertion_with_gap_open(self):
        a = self.align_mgr
        a.set_read(b'ACGT')
        a.set_reference(b'TTTTACTCCCCC')
        res = a.align(gap_open=3)

        self.assertEqual(res.optimal_score, 4)
        self.assertEqual(res.sub_optimal_score, 0)

        self.assertEqual(res.reference_start, 4)
        self.assertEqual(res.reference_end, 5)
        self.assertEqual(res.read_start, 0)
        self.assertEqual(res.read_end, 1)

        # a.print_result(res)

    def test_insertion_with_no_gap_open_penalty(self):
        a = self.align_mgr
        a.set_read(b'ACGT')
        a.set_reference(b'TTTTACTCCCCC')
        res = a.align(gap_open=0)

        self.assertEqual(res.optimal_score, 6)
        self.assertEqual(res.sub_optimal_score, 0)

        self.assertEqual(res.reference_start, 4)
        self.assertEqual(res.reference_end, 6)
        self.assertEqual(res.read_start, 0)
        self.assertEqual(res.read_end, 3)

        # a.print_result(res)

    def test_start_idx_test(self):
        a = self.align_mgr
        a.set_read(b'ACTG')
        a.set_reference(b'ACTCACTG')
        res = a.align(start_idx=4)

        self.assertEqual(res.optimal_score, 8)
        self.assertEqual(res.sub_optimal_score, 0)

        self.assertEqual(res.reference_start, 0)
        self.assertEqual(res.reference_end, 3)
        self.assertEqual(res.read_start, 0)
        self.assertEqual(res.read_end, 3)

        # a.print_result(res, start_idx=4)

    def test_force_align(self):
        read = b'ACTG'
        ref = b'TTTTCTGCCCCCACG'
        res = force_align(read, ref)

        self.assertEqual(res.optimal_score, 6)
        self.assertEqual(res.sub_optimal_score, 0)

        self.assertEqual(res.reference_start, 4)
        self.assertEqual(res.reference_end, 6)
        self.assertEqual(res.read_start, 1)
        self.assertEqual(res.read_end, 3)

        # self.a.set_read(read)
        # self.a.set_reference(ref)
        # self.a.print_result(res)

        out = format_force_align(read, ref, res)
        self.assertEqual(out, ('TTTTCTGCCCCCACG', '   ACTG'))

    def test_bitwise_flag(self):
        a = self.align_mgr
        a.set_read(READ_SEQ)
        a.set_reference(REF_SEQ)

        # 1. Test BitwiseAlignmentFlag.best_idxs
        res = a.align(
            gap_open=3,
            gap_extension=1,
            bitwise_flag=BitwiseAlignmentFlag.best_idxs,
        )

        self.assertEqual(res.CIGAR, '9M1I5M')

        self.assertEqual(res.optimal_score, 21)
        self.assertEqual(res.sub_optimal_score, 8)

        self.assertEqual(res.reference_start, 8)
        self.assertEqual(res.reference_end, 21)

        self.assertEqual(res.read_start, 0)
        self.assertEqual(res.read_end, 14)

        # 2. Test BitwiseAlignmentFlag.best_idxs_no_cigar
        res = a.align(
            gap_open=3,
            gap_extension=1,
            bitwise_flag=BitwiseAlignmentFlag.best_idxs_no_cigar,
        )

        self.assertEqual(res.CIGAR, '')

        self.assertEqual(res.optimal_score, 21)
        self.assertEqual(res.sub_optimal_score, 8)

        self.assertEqual(res.reference_start, 8)
        self.assertEqual(res.reference_end, 21)

        self.assertEqual(res.read_start, 0)
        self.assertEqual(res.read_end, 14)

        # 3. Test BitwiseAlignmentFlag.end_idxs_only_no_cigar
        res = a.align(
            gap_open=3,
            gap_extension=1,
            bitwise_flag=BitwiseAlignmentFlag.end_idxs_only_no_cigar,
        )

        self.assertEqual(res.CIGAR, '')

        self.assertEqual(res.optimal_score, 21)
        self.assertEqual(res.sub_optimal_score, 8)

        self.assertEqual(res.reference_start, -1)
        self.assertEqual(res.reference_end, 21)

        self.assertEqual(res.read_start, -1)
        self.assertEqual(res.read_end, 14)

        # 4. Test BitwiseAlignmentFlag.score_filter
        # 4.a change ``score_cutoff`` to 0
        res = a.align(
            gap_open=3,
            gap_extension=1,
            bitwise_flag=BitwiseAlignmentFlag.score_filter,
            score_cutoff=0,
        )

        self.assertEqual(res.CIGAR, '9M1I5M')

        self.assertEqual(res.optimal_score, 21)
        self.assertEqual(res.sub_optimal_score, 8)

        self.assertEqual(res.reference_start, 8)
        self.assertEqual(res.reference_end, 21)

        self.assertEqual(res.read_start, 0)
        self.assertEqual(res.read_end, 14)

        # 4.b change ``score_cutoff`` to 22
        res = a.align(
            gap_open=3,
            gap_extension=1,
            bitwise_flag=BitwiseAlignmentFlag.score_filter,
            score_cutoff=22,
        )

        self.assertEqual(res.CIGAR, '')

        self.assertEqual(res.optimal_score, 21)
        self.assertEqual(res.sub_optimal_score, 8)

        self.assertEqual(res.reference_start, -1)
        self.assertEqual(res.reference_end, 21)

        self.assertEqual(res.read_start, -1)
        self.assertEqual(res.read_end, 14)

        # 5. Test BitwiseAlignmentFlag.distance_filter
        # 5.a change ``distance_cutoff`` to 0
        res = a.align(
            gap_open=3,
            gap_extension=1,
            bitwise_flag=BitwiseAlignmentFlag.distance_filter,
            distance_cutoff=0,
        )

        self.assertEqual(res.CIGAR, '')

        self.assertEqual(res.optimal_score, 21)
        self.assertEqual(res.sub_optimal_score, 8)

        self.assertEqual(res.reference_start, 8)
        self.assertEqual(res.reference_end, 21)

        self.assertEqual(res.read_start, 0)
        self.assertEqual(res.read_end, 14)

        # 5.b change ``distance_cutoff`` to 13
        # (equal to the READ_SEQ length - 2)
        res = a.align(
            gap_open=3,
            gap_extension=1,
            bitwise_flag=BitwiseAlignmentFlag.distance_filter,
            distance_cutoff=13,
        )

        self.assertEqual(res.CIGAR, '')

        self.assertEqual(res.optimal_score, 21)
        self.assertEqual(res.sub_optimal_score, 8)

        self.assertEqual(res.reference_start, 8)
        self.assertEqual(res.reference_end, 21)

        self.assertEqual(res.read_start, 0)
        self.assertEqual(res.read_end, 14)

        # 5.c change ``distance_cutoff`` to 14
        # (equal to the READ_SEQ length - 1)
        res = a.align(
            gap_open=3,
            gap_extension=1,
            bitwise_flag=BitwiseAlignmentFlag.distance_filter,
            distance_cutoff=len(READ_SEQ) - 1,
        )

        self.assertEqual(res.CIGAR, '9M1I5M')
        self.assertEqual(res.optimal_score, 21)
        self.assertEqual(res.sub_optimal_score, 8)

        self.assertEqual(res.reference_start, 8)
        self.assertEqual(res.reference_end, 21)

        self.assertEqual(res.read_start, 0)
        self.assertEqual(res.read_end, 14)
