#include <stdint.h>

#define stid(tid, strand) ((((char)strand)=='-')?(INT32_MIN + (int32_t)tid):((int32_t)tid))
#define stidpos(tid, strand, pos) ((((uint64_t)stid(tid, strand))<<32u)+(uint32_t) (pos))
