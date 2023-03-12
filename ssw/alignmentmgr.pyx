# cython: language_level=3, boundscheck=False, wraparound=False
# Copyright (C) 2023, Nick Conway
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
ssw.alignmentmgr
~~~~~~~~~~~~~~~~

Alignment Manager Cython wrapper for ssw library code.

NOTE: See this link for a info on the CIGAR format:
`What is a CIGAR? <http://genome.sph.umich.edu/wiki/SAM#What_is_a_CIGAR.3F>`_

Refer to
`Wiki Smith–Waterman_algorithm <https://en.wikipedia.org/wiki/Smith–Waterman_algorithm>`_
to understand Smith-Waterman scoring

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

    s_profile* ssw_init(
        const int8_t*, const int32_t,
        const int8_t*, const int32_t,
        const int8_t
    )

    void init_destroy(s_profile*)

    s_align* ssw_align(
        const s_profile*, const int8_t*, int32_t, const uint8_t,
        const uint8_t, const uint8_t,
        const uint16_t, const int32_t, const int32_t
    ) nogil

    void align_destroy(s_align*)
    char cigar_int_to_op(uint32_t)
    uint32_t cigar_int_to_len(uint32_t)
    uint32_t to_cigar_int(uint32_t, char)


import warnings as pywarnings
from typing import (
    Optional,
    Tuple,
    Union,
)

from .alignmenttuple import Alignment

SNAKE_CASE_DEPRECATED_MSG = 'Function deprecated please use "%s" instead'


class BitwiseAlignmentFlag:
    '''``bitwise_flag`` (from high- to low-bit) for :meth:`AlignmentMgr.align`

    Attributes:
        best_idxs_no_cigar: Bit 5. When set as 1, ``ssw_align`` will return
            the best alignment beginning position;
            NOTE: this is setting ``bitwise_flag == 8``

        distance_filter: Bit 6.  When set as 1::

                if (
                    ((reference_end - reference_start) < distance_filter) and
                    ((read_end - read_start) < distance_filter)
                )

            (whatever bit 5 is set) the function will return the best
            alignment beginning position and cigar;
            NOTE: this is setting ``bitwise_flag == 4``

        score_filter: Bit 7. When set as 1, if the best
            ``alignment_score >= score_filter``, (whatever bit 5 is set as) the
            function  will return the best alignment beginning position and
            cigar; NOTE: this is setting ``bitwise_flag == 2``

        best_idxs: Bit 8. When set as 1, (whatever bit 5, 6 or 7 are set as) the
            function will always return the best alignment beginning position
            and cigar.
            NOTE: this is setting ``bitwise_flag = 1``

        end_idxs_only_no_cigar: When ``bitwise_flag == 0``, only the optimal
            and sub-optimal scores and the optimal alignment ending position
            will be returned.
    '''
    end_idxs_only_no_cigar =    0
    best_idxs =                 1  # includes cigar
    score_filter =              2  # cutoff at ``score_filter`` arg
    distance_filter =           4  # cutoff at ``distance_filter`` arg distance
    best_idxs_no_cigar =        8


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

    Attributes:
        read: Read sequence python string or bytes-string
        reference: Reference sequence python string or bytes-string
        match_score: 0 to 255 value for scoring matches
        mismatch_penalty: 0 to 255 value penalty for mismatches
    '''
    cdef:
        int8_t* score_matrix
        s_profile* profile
        int8_t* read_arr
        Py_ssize_t read_length
        int32_t num_elements_sq_root

        int8_t* ref_arr
        Py_ssize_t ref_length


    cdef:
        object read       # type: STR_T
        object reference  # type: STR_T
        readonly int _match_score
        readonly int _mismatch_penalty

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
            match_score: 0 to 255 value for scoring matches
            mismatch_penalty: 0 to 255 value penalty for mismatches

        '''
        # Hard coded to a 5 x 5 matrix
        self.score_matrix = <int8_t*> PyMem_Malloc(25 * sizeof(int8_t))

        # the square root of the number of elements in `self.score_matrix`
        # (``self.score_matrix`` has
        # ``num_elements_sq_root*num_elements_sq_root`` elements)
        self.num_elements_sq_root = 5

        self.match_score = match_score
        self.mismatch_penalty = mismatch_penalty

        self.build_dna_score_matrix()

        self.read = ''
        self.reference = ''
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

    @property
    def match_score(self) -> int:
        '''0 to 255 value for scoring matches

        '''
        return self._match_score

    @match_score.setter
    def match_score(self, value: int):
        '''Set the :attr:`match_score` attribute
        Must be in range 0 to 255

        Raises:
            ValueError: Must be in range(0, 255)
        '''
        if value > 255 or value < 0:
            raise ValueError(f'{value} must be in range(0, 255)')
        self._match_score = value

    @property
    def mismatch_penalty(self) -> int:
        '''0 to 255 penalty value for scoring mismatches

        '''
        return self._mismatch_penalty

    @mismatch_penalty.setter
    def mismatch_penalty(self, value: int):
        '''Set the :attr:`mismatch_penalty` attribute
        Must be in range 0 to 255

        Raises:
            ValueError: Must be in range(0, 255)
        '''
        if value > 255 or value < 0:
            raise ValueError(f'{value} must be in range(0, 255)')
        self._mismatch_penalty = value

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
            result: :class:`ssw.alignmenttuple.Alignment` containing the result from
                a call to :meth:`align`
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
                Must be set for a call to :meth:`align`

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

            # ``score_size`` is the estimated Smith-Waterman score; if your
            # estimated best alignment score is surely < 255 please set 0;
            # if your estimated best alignment score >= 255, please set 1; if
            # you don't know, please
            # set 2
            int8_t score_size = 2  # Don't know best score size

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
            self.num_elements_sq_root,
            score_size,
        )
    # end def

    def setRead(self, read: STR_T):
        '''.. deprecated:: 1.0.0. Choose :meth:`set_read` method instead

        '''
        pywarnings.warn(SNAKE_CASE_DEPRECATED_MSG % 'set_read')
        self.set_read(read)

    def set_reference(self, reference: STR_T):
        '''Set the query reference string.

        Args:
            reference: String-like (str or bytestring) that represents the
                reference sequence. Must be set a call to :meth:`align`

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
            uint8_t bitwise_flag,
            uint16_t distance_filter,
            uint16_t score_filter,
    ) except NULL:
        '''C version of the alignment code

        Args:
            gap_open: The absolute value of gap open penalty
            gap_extension: The absolute value of gap extension penalty
            start_idx: The start index of the reference array to begin the s
                rearch
            mod_ref_length: Length of the target reference sequence
            bitwise_flag: Flag using a member of :class:`BitwiseAlignmentFlag`
            distance_cutoff: Filter by index start to end distances less than
                this value.  see :class:`BitwiseAlignmentFlag` docs for details
                ``bitwise_flag=BitwiseAlignmentFlag.distance_filter``
            score_cutoff: Filter by fits with a score greater score value cutoff
                see :class:`BitwiseAlignmentFlag` docs for details
                ``bitwise_flag=BitwiseAlignmentFlag.score_filter``

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

            # Should be read_length divided by 2
            int32_t mask_len = read_length // 2

            s_profile* profile = self.profile
            int8_t* ref_arr = self.ref_arr

        # Guard the ``mask_len`` parameter
        mask_len = 15 if mask_len < 15 else mask_len

        if profile != NULL:
            with nogil:
                result = ssw_align(
                    self.profile,
                    &ref_arr[start_idx],
                    mod_ref_length,
                    gap_open,
                    gap_extension,
                    bitwise_flag,
                    score_filter,
                    distance_filter,
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
            bitwise_flag: int = BitwiseAlignmentFlag.best_idxs,
            distance_cutoff: int = 0,
            score_cutoff: int = 0,
    ) -> Alignment:
        '''Align a read to the reference with optional index offseting

        Returns a dictionary no matter what as `align_c` can't return ``NULL``

        Args:
            gap_open: Penalty for ``gap_open``. default 3
            gap_extension: Penalty for ``gap_extension``. default 1
            start_idx: Index to start search. Default 0
            end_idx: Index to end search (trying to avoid a
                target region). Default 0 means use whole reference length
            bitwise_flag: Flag using a member of :class:`BitwiseAlignmentFlag`
            distance_cutoff: Filter by index start to end distances less than
                this value.  see :class:`BitwiseAlignmentFlag` docs for details
                ``bitwise_flag=BitwiseAlignmentFlag.distance_filter``
            score_cutoff: Filter by fits with a score greater score value cutoff
                see :class:`BitwiseAlignmentFlag` docs for details
                ``bitwise_flag=BitwiseAlignmentFlag.score_filter``

        Returns:
            :class:`ssw.alignmenttuple.Alignment` instance

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

        # NOTE: .. deprecated ``None`` type for Alignment.CIGAR
        # cigar_str = None
        cigar_str = ''
        cigar_str_list = []
        if self.reference is None:
            raise ValueError('Call set_reference first')


        #`distance_filter`: when bit 6 of flag is set as 1 and bit 8 is
        # set as 0, `distance_filter` will be used (Please check the
        # decription of the flag parameter for detailed usage.)
        cdef uint16_t _distance_cutoff = distance_cutoff

        #`score_cutoff`: when bit 7 of flag is set as 1 and bit 8 is
        # set as 0, `score_cutoff` will be used (Please check the
        # decription of the flag parameter for detailed usage.)
        cdef uint16_t _score_cutoff = score_cutoff

        cdef s_align* result =  self.align_c(
            gap_open,
            gap_extension,
            start_idx,
            search_length,
            <uint8_t> bitwise_flag,
            _distance_cutoff,
            _score_cutoff,
        )
        cigar_tuple_list = []
        if result.cigar != NULL:
            for c in range(result.cigarLen):
                letter = cigar_int_to_op(result.cigar[c])
                letter_int = letter
                length = cigar_int_to_len(result.cigar[c])
                idx_tuple = (<int>length, chr(letter_int))
                cigar_str_list.append('%d%s' % idx_tuple)
                cigar_tuple_list.append(idx_tuple)
            cigar_str = ''.join(cigar_str_list)
        out = Alignment(
            cigar_str,
            result.score1,
            result.score2,
            result.ref_begin1,
            result.ref_end1,
            result.read_begin1,
            result.read_end1,
            cigar_tuple_list,
        )
        align_destroy(result)
        return out
    # end def

    def build_dna_score_matrix(
            self,
    ):
        '''Mismatch_penalty should be positive

        The score matrix looks like:

                            A,  C,  G,  T,  N
        score_matrix  = {   2, -2, -2, -2,  0, // A
                           -2,  2, -2, -2,  0, // C
                           -2, -2,  2, -2,  0, // G
                           -2, -2, -2,  2,  0, // T
                            0,  0,  0,  0,  0  // N
                        }

        Hard coded to a 5 x 5 matrix so this should not be called externally
        as of now

        Uses:
            :attr:`match_score`: Match score
            :attr:`mismatch_penalty`: Mismatch penalty
            :attr:`score_matrix`: Pointer to matrix to populate

        '''
        cdef:
            Py_ssize_t i, j
            Py_ssize_t idx = 0
            int8_t* matrix = self.score_matrix

            int8_t match_score = <uint8_t>   self._match_score
            int8_t neg_mismatch_penalty = <int8_t> (-self._mismatch_penalty)

        for i in range(4):
            for j in range(4):
                if i == j:
                    matrix[idx] = match_score
                else:
                    matrix[idx] = neg_mismatch_penalty
                inc(idx)
            matrix[idx] = 0
            inc(idx)
        for i in range(5):
            matrix[inc(idx)] = 0


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
        :class:`ssw.alignmenttuple.Alignment` result

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
        alignment: :class:`ssw.alignmenttuple.Alignment` named tuple
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
