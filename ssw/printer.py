# -*- coding: utf-8 -*-
from typing import Union

from ssw import (
    SSW,
    Alignment
)

STR_T = Union[str, bytes]

def printer(alignment: Alignment, reference: STR_T, read: STR_T):
    print("optimal_score: %d\tsub-optimal_score: %d\t" %
        (alignment.optimal_score, alignment.sub_optimal_score))

    reference_start: int = alignment.reference_start

    if reference_start >= 0:
        print("reference_start: %d\t" % (reference_start), end='')
    print("reference_end: %d\t" % (alignment.reference_end) )

    read_start: int = alignment.read_start

    if read_start >= 0:
        print("read_start: %d\t" % (read_start), end='')
    print("read_end: %d\n" % (alignment.read_end))

    cigar: str = alignment.CIGAR

    if cigar is not None:
        e: int      = 0
        left: int   = 0
        pb: int = read_start
        qb: int = reference_start
        STEP2: bool = False
        STEP3: bool = False
        STEPEND: bool = False
        lim: int = len(cigar)

        while e < lim or left > 0:
            q: int = qb
            p: int = pb
            count: int = 0

            print("Refer:  %8d    " % (q), end='')
            for c in range(e, lim, 2):
                length = int(cigar[c])
                letter = cigar[c + 1]
                L1 = left if (count == 0 and left > 0) else length
                for j in range(L1):
                    if letter == 'I':
                        print("-", end='')
                    else:
                        print(reference[q], end='')
                        q += 1
                    count += 1
                    if count == 60:
                        STEP2 = True
                        break # goto STEP 2
                if STEP2:
                    STEP2 = False
                    break
            # end for
            # STEP 2
            print("    %d\n                    " % (q - 1), end='')
            q = qb
            count = 0
            for c in range(e, lim, 2):
                length = int(cigar[c])
                letter = cigar[c + 1]
                L1 = left if (count == 0 and left > 0) else length
                for j in range(L1):
                    if letter == 'M':
                        if reference[q] == read[p]:
                            print('|', end='')
                        else:
                            print("*", end='')
                        q += 1
                        p += 1
                    else:
                        print('*', end='')
                        if letter == 'I':
                            p += 1
                        else:
                            q += 1
                    count += 1
                    if count == 60:
                        qb = q
                        STEP3 = True
                        break # goto STEP 3
                if STEP3:
                    STEP3 = False
                    break
            # end for
            # STEP 3
            p = pb
            print("\nRead:   %8d    " % (p), end='')
            count = 0
            for c in range(e, lim, 2):
                length = int(cigar[c])
                letter = cigar[c + 1]
                L1 = left if (count == 0 and left > 0) else length
                for j in range(L1):
                    if letter == 'D':
                        print('-', end='')
                    else:
                        print(read[p], end='')
                        p += 1
                    count += 1
                    if count == 60:
                        pb = p
                        left = L1 - j - 1
                        e = c + 2 if left == 0 else c
                        STEPEND = True
                        break # goto STEPEND
                if STEPEND:
                    STEPEND = False
                    break
            e = c + 2
            left = 0
            # STEPEND
            print("    %d\n" % (p - 1))
        # end while
    #end if
# end def