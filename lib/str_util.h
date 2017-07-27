#include <inttypes.h>
#include "ssw.h"

void dnaToInt8(const char* c_str, int8_t* arr, uint32_t len);
void ssw_write_cigar(const s_align* a);
void ssw_writer(const s_align* a, const char* ref_seq, const char* read_seq);