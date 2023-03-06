/*
* Copyright 2023 Nick Conway; Copyright 2018, Nick Conway; Wyss Institute
* Harvard University
*
* lib/str_util.h
*
* Utilities for assisting DNA C-string to integer conversion
* and CIGAR format IO
*/
#ifndef SSWPY_STR_UTIL_H
#define SSWPY_STR_UTIL_H

#include <inttypes.h>
#include "ssw.h"

void dna_to_int8(const char* c_str, int8_t* arr, uint32_t len);
void ssw_write_cigar(const s_align* a);
void ssw_writer(const s_align* a, const char* ref_seq, const char* read_seq);

#endif
