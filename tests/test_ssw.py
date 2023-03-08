# -*- coding: utf-8 -*-
# Copyright 2023 Nick Conway; Copyright 2018, Nick Conway; Wyss Institute
# Harvard University
#
# See LICENSE.md for full MIT license.
'''
tests.test_ssw
~~~~~~~~~

Tests for the ssw-py
'''
import unittest

from ssw import (  # type: ignore
    AlignmentMgr,
    force_align,
    format_force_align,
)


class TestAlignmentMgr(unittest.TestCase):
    '''Tests for the AlignmentMgr class
    '''

    def setUp(self):
        self.align_mgr = AlignmentMgr()
        with self.assertWarns(UserWarning):
            self.align_mgr.setRead(b'ACGT')

        self.align_mgr.set_read(b'ACGT')

    def test_exact(self):
        a = self.align_mgr
        ref = b'TTTTACGTCCCCC'
        a.set_reference(ref)
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
        ref = b'TTTTACAGTCCCCC'
        a.set_reference(ref)
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
        ref = b'TTTTACTCCCCC'
        a.set_reference(ref)
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
        ref = b'TTTTACTCCCCC'
        a.set_reference(ref)
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
        a.set_read('ACTG')
        a.set_reference('ACTCACTG')
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
