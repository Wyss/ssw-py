
#include "str_util.h"
#include <Python.h>
#include <stdio.h>

static const int8_t DNA_BASE_LUT[128] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  //   A     C            G
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
  //             T
    4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  //   a     c            g
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
  //             t
    4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

void dnaToInt8(const char* c_str, int8_t* arr, uint32_t len) {
    /* Assume ascii */
    const char* str_lim  = c_str + len;
    while(c_str < str_lim) {
        *arr++ = DNA_BASE_LUT[(unsigned char) *c_str++];
    }
}

void ssw_write_cigar(const s_align* a) {
    int32_t c;
    if (a->cigar) {
        printf("CIGAR start index %d: ", a->ref_begin1);
        for (c = 0; c < a->cigarLen; ++c) {
            char letter = cigar_int_to_op(a->cigar[c]);
            uint32_t length = cigar_int_to_len(a->cigar[c]);
            printf("%d%c", length, letter);
        }
        printf("\n");
    } else {
        printf("no CIGAR\n");
    }
}

//  Print the BLAST like output.
void ssw_writer(const s_align* a,
      const char* ref_seq,
      const char* read_seq) {

    const int8_t* table = DNA_BASE_LUT;
    fprintf(stdout, "optimal_score: %d\tsub-optimal_score: %d\t\n", a->score1, a->score2);
    if (a->ref_begin1 >= 0) {
        fprintf(stdout, "target_begin: %d\t", a->ref_begin1);
    }
    fprintf(stdout, "target_end: %d\t\n", a->ref_end1);
    if (a->read_begin1 >= 0) {
        fprintf(stdout, "query_begin: %d\t", a->read_begin1);
    }
    fprintf(stdout, "query_end: %d\n\n", a->read_end1);
    if (a->cigar) {
        int32_t c = 0, left = 0, e = 0, qb = a->ref_begin1, pb = a->read_begin1;
        uint32_t i;
        while (e < a->cigarLen || left > 0) {
            int32_t count = 0;
            int32_t q = qb;
            int32_t p = pb;
            fprintf(stdout, "Target: %8d    ", q);
            for (c = e; c < a->cigarLen; ++c) {
                char letter = cigar_int_to_op(a->cigar[c]);
                uint32_t length = cigar_int_to_len(a->cigar[c]);
                uint32_t l = (count == 0 && left > 0) ? left: length;
                for (i = 0; i < l; ++i) {
                  if (letter == 'I') {
                      fprintf(stdout, "-");
                  } else {
                      fprintf(stdout, "%c", *(ref_seq + q));
                      ++q;
                  }
                  ++count;
                  if (count == 60) {
                      goto step2;
                  }
                }
            }
step2:
            fprintf(stdout, "    %d\n                    ", q-1);
            q = qb;
            count = 0;
            for (c = e; c < a->cigarLen; ++c) {
                char letter = cigar_int_to_op(a->cigar[c]);
                uint32_t length = cigar_int_to_len(a->cigar[c]);
                uint32_t l = (count == 0 && left > 0) ? left: length;
                for (i = 0; i < l; ++i) {
                    if (letter == 'M') {
                        if (table[(int)*(ref_seq + q)] == table[(int)*(read_seq + p)]) {
                            fprintf(stdout, "|");
                        } else {
                            fprintf(stdout, "*");
                        }
                        ++q;
                        ++p;
                    } else {
                      fprintf(stdout, "*");
                      if (letter == 'I') { 
                          ++p;
                      } else {
                          ++q;
                      }
                    }
                    ++count;
                    if (count == 60) {
                        qb = q;
                        goto step3;
                    }
                }
            }
step3:
            p = pb;
            fprintf(stdout, "\nQuery:  %8d    ", p);
            count = 0;
            for (c = e; c < a->cigarLen; ++c) {
                char letter = cigar_int_to_op(a->cigar[c]);
                uint32_t length = cigar_int_to_len(a->cigar[c]);
                uint32_t l = (count == 0 && left > 0) ? left: length;
                for (i = 0; i < l; ++i) {
                    if (letter == 'D') {
                        fprintf(stdout, "-");
                    } else {
                       fprintf(stdout, "%c", *(read_seq + p));
                       ++p;
                    }
                    ++ count;
                    if (count == 60) {
                      pb = p;
                      left = l - i - 1;
                      e = (left == 0) ? (c + 1) : c;
                      goto end;
                    }
                }
            }
            e = c;
            left = 0;
end:
            fprintf(stdout, "    %d\n\n", p-1);
        }
    }
}