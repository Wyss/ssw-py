# cython: language_level=3, boundscheck=False, wraparound=False
# Copyright 2023 Nick Conway; Copyright 2018, Nick Conway; Wyss Institute
# Harvard University
#
# See LICENSE.md for full MIT license.
'''
ssw.alignmentmgr
~~~~~~~~~~~~~~~~

Alignment Manager Cython wrapper for ssw library code.

NOTE: See this link for a info on the CIGAR format:
`What is a CIGAR? <http://genome.sph.umich.edu/wiki/SAM#What_is_a_CIGAR.3F>`_

'''
from cpython.bytes cimport (
    PyBytes_AsStringAndSize,
    PyBytes_Check,
)
from cpython.mem cimport (
    PyMem_Free,
    PyMem_Malloc,
)
from cpython.unicode cimport PyUnicode_AsUTF8AndSize
from cython.operator cimport postincrement as inc
from libc.stdint cimport (
    int8_t,
    int32_t,
    uint8_t,
    uint16_t,
    uint32_t,
)


cdef extern from "str_util.h":
    void dna_to_int8(const char*, int8_t*, int32_t)
    void ssw_write_cigar(const s_align*)
    void ssw_writer(const s_align*, const char*, const char*)

cdef extern from "ssw.h":
    # leave out a few members
    ctypedef struct s_profile:
        const int8_t* read
        const int8_t* mat
        int32_t readLen
        int32_t n
        uint8_t bias

    ctypedef struct s_align:
        uint16_t score1
        uint16_t score2
        int32_t ref_begin1
        int32_t ref_end1
        int32_t read_begin1
        int32_t read_end1
        int32_t ref_end2
        uint32_t* cigar
        int32_t cigarLen

    s_profile* ssw_init(const int8_t*, const int32_t, const int8_t*, const int32_t, const int8_t)
    void init_destroy (s_profile*)
    s_align* ssw_align (const s_profile*, const int8_t*, int32_t, const uint8_t, const uint8_t, const uint8_t, const uint16_t, const int32_t, const int32_t)
    void align_destroy (s_align*)
    char cigar_int_to_op (uint32_t)
    uint32_t cigar_int_to_len(uint32_t)
    uint32_t to_cigar_int(uint32_t, char)


import warnings as pywarnings
from typing import (
    Optional,
    Tuple,
    Union,
)

from .alignment import Alignment

SNAKE_CASE_DEPRECATED_MSG = 'Function deprecated please use "%s" instead'


cdef inline str convert_bytes_to_str(s):
    if isinstance(s, bytes):
        # encode to the specific encoding used inside of the module
        return (<bytes>s).decode('utf8')
    else:
        return s


cdef inline char* obj_to_cstr_len(
        object o1,
        Py_ssize_t *length,
) except NULL:
    '''Convert a Python string or bytes-string object to a C string

    Args:
        o1 - python string or bytes-string object

    Returns:
        C char* pointer to the internal string of the o1
    '''
    cdef char* c_str1
    if PyBytes_Check(o1):
        if PyBytes_AsStringAndSize(o1, &(c_str1), length) == -1:
            raise TypeError("obj_to_cstr: PyBytes_AsStringAndSize error")
        return <char*> c_str1
    else:
        c_str1 = <char*> PyUnicode_AsUTF8AndSize(o1, length)
        if c_str1 == NULL:
            raise OSError("obj_to_cstr: PyUnicode_AsUTF8AndSize error")
    return c_str1


STR_T = Union[str, bytes]


cdef class AlignmentMgr:
    '''Class to manage SSW-based alignment
    '''
    cdef:
        int8_t* score_matrix
        s_profile* profile
        int8_t* read_arr
        Py_ssize_t read_length

        int8_t* ref_arr
        Py_ssize_t ref_length

    cdef:
        object read       # type: STR_T
        object reference  # type: STR_T

    def __cinit__(
            self,
            int match_score=2,
            int mismatch_penalty=2,
    ):
        self.score_matrix = NULL
        self.profile = NULL
        self.read_arr = NULL
        self.ref_arr = NULL
    # end def

    def __init__(
            self,
            int match_score=2,
            int mismatch_penalty=2,
    ):
        '''Class initialization method

        Args:
            match_score: for scoring matches
            mismatch_penalty: for scoring mismatches

        '''
        self.score_matrix = <int8_t*> PyMem_Malloc(25 * sizeof(int8_t))
        self.build_dna_score_matrix(
            <uint8_t>match_score,
            <uint8_t> mismatch_penalty,
            self.score_matrix,
        )
        self.read = None
        self.reference = None
    # end def

    def __dealloc__(self):
        PyMem_Free(self.score_matrix)

        if self.profile != NULL:
            init_destroy(self.profile)
            self.profile = NULL

        if self.read_arr != NULL:
            PyMem_Free(self.read_arr)

        if self.ref_arr != NULL:
            PyMem_Free(self.ref_arr)
    # end def

    cdef int print_result_c(
            self,
            s_align* result,
            Py_ssize_t start_idx,
    ) except -1:
        '''Print the ``result`` argument and ``self.read``

        Args:
            result: ``s_align`` pointer to result
            start_idx: start index to print from

        Returns:
            0

        Raises:
            ValueError: result argument is NULL

        '''
        cdef:
            const char* read_cstr = NULL
            Py_ssize_t read_length, ref_length
            const char* ref_cstr = NULL

        print(self.read)
        if result != NULL:
            read_cstr = obj_to_cstr_len(self.read, &read_length)
            ref_cstr = obj_to_cstr_len(self.reference, &ref_length)

            ssw_write_cigar(result)
            ssw_writer(
                result,
                <const char*> &ref_cstr[start_idx],
                <const char*> read_cstr,
            )
        else:
            raise ValueError('result argument is NULL')
        return 0


    def print_result(self, result: Alignment, start_idx: int = 0):
        '''Rebuild an ``s_align`` struct from a result dictionary
        so as to be able to call ssw terminal print functions

        Args:
            result: :class:`ssw.alignment.Alignment` containing the result from a
                call to :meth:`align`
            start_idx: index to start printing from. defaults to 0

        Raises:
            OSError: Memory allocation issue

        '''
        cdef:
            uint32_t* cigar_array = NULL
            s_align *res_align = NULL
            char letter
            uint32_t length
            Py_ssize_t i, j

        try:
            res_align = <s_align*> PyMem_Malloc(
                sizeof(s_align),
            )
            if res_align == NULL:
                raise OSError('Memory allocation issue')
            cigar = result.CIGAR
            if cigar is not None:
                cigar_array = <uint32_t*> PyMem_Malloc(
                    len(cigar) // 2*sizeof(uint32_t)
                )
                if cigar_array == NULL:
                    raise OSError('Memory allocation issue')
                j = 0
                for i in range(0, len(cigar), 2):
                    length = int(cigar[i])
                    letter = ord(cigar[i+1])
                    cigar_array[j] = to_cigar_int(length, letter)
                    j += 1
                res_align.cigar = cigar_array
                res_align.cigarLen = len(cigar)//2
            else:
                res_align.cigar = NULL

            res_align.score1 = result.optimal_score
            res_align.score2 = result.sub_optimal_score
            res_align.ref_begin1 = result.reference_start
            res_align.ref_end1 = result.reference_end
            res_align.read_begin1 = result.read_start
            res_align.read_end1 = result.read_end

            self.print_result_c(res_align, start_idx)
        finally:
            if cigar is not None:
                PyMem_Free(cigar_array)
            PyMem_Free(res_align)

    def printResult(self, result: Alignment, start_idx: int = 0):
        '''.. deprecated:: 1.0.0. Choose :meth:`print_result` method instead

        '''
        pywarnings.warn(SNAKE_CASE_DEPRECATED_MSG % 'print_result')
        self.print_result(result, start_idx)


    def set_read(self, read: STR_T):
        '''Set the query read string

        Args:
            read: String-like (str or bytestring) that represents the read.
                Must be set for a prpoer run

        '''
        cdef:
            Py_ssize_t read_length
            const char* read_cstr = obj_to_cstr_len(
                read,
                &read_length,
            )
            int8_t* read_arr = <int8_t*> PyMem_Malloc(
                read_length * sizeof(char)
            )

        if self.profile != NULL:
            init_destroy(self.profile)
            self.profile = NULL
        if self.read_arr != NULL:
            PyMem_Free(self.read_arr)
            self.read_arr = NULL

        self.read = read
        self.read_arr = read_arr
        self.read_length = read_length

        dna_to_int8(
            read_cstr,
            read_arr,
            read_length,
        )

        self.profile = ssw_init(
            read_arr,
            <int32_t> read_length,
            self.score_matrix,
            5,
            2, # don't know best score size
        )
    # end def

    def setRead(self, read: STR_T):
        '''.. deprecated:: 1.0.0. Choose :meth:`set_read` method instead

        '''
        pywarnings.warn(SNAKE_CASE_DEPRECATED_MSG % 'set_read')
        self.set_read(read)

    def set_reference(self, reference: STR_T):
        '''Set the query reference string

        Args:
            reference: String-like (str or bytestring) that represents the
                reference sequence must be set

        '''
        cdef:
            Py_ssize_t ref_length
            const char* ref_cstr = obj_to_cstr_len(
                reference,
                &ref_length
            )
            int8_t* ref_arr = <int8_t*> PyMem_Malloc(
                ref_length * sizeof(char)
            )

        dna_to_int8(ref_cstr, ref_arr, ref_length)
        self.reference = reference
        if self.ref_arr != NULL:
            PyMem_Free(self.ref_arr)
            self.ref_arr = NULL
        self.ref_arr = ref_arr
        self.ref_length = ref_length
    # end def

    def setReference(self, reference: STR_T):
        '''.. deprecated:: 1.0.0. Choose :meth:`set_reference` method instead

        '''
        pywarnings.warn(SNAKE_CASE_DEPRECATED_MSG % 'set_reference')
        self.set_reference(reference)

    cdef s_align* align_c(
            self,
            int gap_open,
            int gap_extension,
            Py_ssize_t start_idx,
            int32_t mod_ref_length,
    ) except NULL:
        '''C version of the alignment code

        Args:
            gap_open:
            gap_extension:
            start_idx:
            mod_ref_length:

        Returns:
            ``s_align`` pointer results

        Raises:
            ValueError: Must set profile first
            ValueError: Problem Running alignment, see ``stdout``

        '''
        cdef:
            Py_ssize_t read_length
            const char* _ = obj_to_cstr_len(
                self.read,
                &read_length,
            )
            s_align* result = NULL
            int32_t mask_len = read_length // 2

        mask_len = 15 if mask_len < 15 else mask_len

        if self.profile != NULL:
            result = ssw_align(
                self.profile,
                &self.ref_arr[start_idx],
                mod_ref_length,
                gap_open,
                gap_extension,
                1,
                0,
                0,
                mask_len,
            )
        else:
            raise ValueError('Must set profile first via `set_read`')
        if result == NULL:
            raise ValueError('Problem Running alignment, see stdout')
        return result
    # end def

    def align(
            self,
            gap_open: int = 3,
            gap_extension: int = 1,
            start_idx: int = 0,
            end_idx: int = 0,
    ) -> Alignment:
        '''Align a read to the reference with optional index offseting

        returns a dictionary no matter what as align_c can't return
        NULL

        Args:
            gap_open: Penalty for gap_open. default 3
            gap_extension: Penalty for gap_extension. default 1
            start_idx: Index to start search. default 0
            end_idx: Index to end search (trying to avoid a
                target region). default 0 means use whole reference length

        Returns:
            Alignment with keys `CIGAR`,        <for depicting alignment>
                                `optimal_score`,
                                `sub-optimal_score`,
                                `reference_start`, <index into reference>
                                `reference_end`,   <index into reference>
                                `read_start`,  <index into read>
                                `read_end`     <index into read>

        Raises:
            ValueError: Negative indexing not supported
            ValueError: ``start_idx`` or ``end_idx`` error
            ValueError: Call :meth:`set_reference` first

        '''
        cdef:
            Py_ssize_t c
            char letter
            int letter_int
            uint32_t length
            int32_t search_length
            Py_ssize_t end_idx_final

        if start_idx < 0 or end_idx < 0:
            raise ValueError('Negative indexing not supported')
        if end_idx > self.ref_length or start_idx > self.ref_length:
            err = (
                f'start_idx: {start_idx} or end_idx: {end_idx} cannot be '
                f'greater than ref_length: {self.ref_length}'
            )
            raise ValueError(err)
        if end_idx == 0:
            end_idx_final = self.ref_length
        else:
            end_idx_final = end_idx
        search_length = end_idx_final - start_idx

        cigar = None

        if self.reference is None:
            raise ValueError('Call set_reference first')

        cdef s_align* result =  self.align_c(
            gap_open,
            gap_extension,
            start_idx,
            search_length,
        )
        if result.cigar != NULL:
            cigar = ""
            for c in range(result.cigarLen):
                letter = cigar_int_to_op(result.cigar[c])
                letter_int = letter
                length = cigar_int_to_len(result.cigar[c])
                cigar += "%d%s" % (<int>length, chr(letter_int))
        out = Alignment(
            cigar,
            result.score1,
            result.score2,
            result.ref_begin1,
            result.ref_end1,
            result.read_begin1,
            result.read_end1,
        )
        align_destroy(result)
        return out
    # end def

    cdef int build_dna_score_matrix(
            self,
            const uint8_t match_score,
            const uint8_t mismatch_penalty,
            int8_t* matrix,
    ) except -1:
        '''Mismatch_penalty should be positive

        The score matrix looks like
                            A,  C,  G,  T,  N
        score_matrix  = {   2, -2, -2, -2,  0, // A
                           -2,  2, -2, -2,  0, // C
                           -2, -2,  2, -2,  0, // G
                           -2, -2, -2,  2,  0, // T
                            0,  0,  0,  0,  0  // N
                        }

        Args:
            match_score: Match score
            mismatch_penalty: Mismatch penalty
            matrix: Pointer to matrix to populate

        Returns:
            0

        '''
        cdef:
            Py_ssize_t i, j
            Py_ssize_t idx = 0

        for i in range(4):
            for j in range(4):
                if i == j:
                    matrix[idx] =  <int8_t> match_score
                else:
                    matrix[idx] = <int8_t> (-mismatch_penalty)
                inc(idx)
            matrix[idx] = 0;
            inc(idx)
        for i in range(5):
            matrix[inc(idx)] = 0
        return 0


def force_align(
        read: STR_T,
        reference: STR_T,
        force_overhang: bool = False,
        aligner: Optional[AlignmentMgr] = None,
) -> Alignment:
    '''Enforces no gaps by raising the ``gap_open`` penalty

    Args:
        read: Read sequence python string or bytes-string
        reference: Reference sequence python string or bytes-string
        force_overhang: Make sure only one end overhangs
        aligner: pass an existing :class:`AlignmentMgr` object

    Returns:
        :class:`ssw.alignment.Alignment` result

    Raises:
        ValueError: No solution found
        ValueError: Read does not align to one overhang

    '''
    a: AlignmentMgr = AlignmentMgr() if aligner is None else aligner
    a.set_read(read)
    a.set_reference(reference)
    len_x: int = len(read)
    # Set the gap_open penalty high to drop all splits in alignment
    # pick the first max hit
    res: Alignment = a.align(gap_open=len_x)
    if res.optimal_score < 4:
        raise ValueError('No solution found')
    if force_overhang:
        # Read must align to either the beginning or end of the reference string
        if (
            res.reference_start != 0 or
            res.reference_end != len(reference) - 1
        ):
            raise ValueError('Read does not align to one overhang')
    return res


def format_force_align(
        read: STR_T,
        reference: STR_T,
        alignment: Alignment,
        do_print: bool = False,
) -> Tuple[str, str]:
    '''Does not truncate strings. Optionally prints these formatted strings

    Args:
        read: Read sequence python string or bytes-string
        reference: Reference sequence python string or bytes-string
        alignment: :class:`ssw.alignment.Alignment` named tuple
        do_print: Default is ``False``. If ``True``, print output

    Returns:
        (<reference output>, <read output>)

    '''
    start_ref: int = alignment.reference_start
    start_read: int = alignment.read_start
    buffer_ref: str = ''
    buffer_read: str = ''
    if start_ref < start_read:
        buffer_ref = ' '*(start_read - start_ref)
    else:
        buffer_read = ' '*(start_ref - start_read)

    ref_out = f'{buffer_ref}{convert_bytes_to_str(reference)}'
    read_out = f'{buffer_read}{convert_bytes_to_str(read)}'
    if do_print:
        print(ref_out)
        print(read_out)
    return ref_out, read_out
