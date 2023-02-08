/* The MIT License (MIT)

   Copyright (c) 2023 Anrui Liu <liuar6@gmail.com>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   “Software”), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
 */

#ifndef __BIOIDX_H
#define __BIOIDX_H

#include <stdint.h>

#define BIOIDX_VERSION "1.0.0"
static inline int32_t bioidx_key(int32_t tid, char strand){
    return ((tid >= 0 && (strand)=='-')?(INT32_MIN+tid):(tid));
}
typedef struct bioidx_t bioidx_t;
typedef struct binidx_itr_t{
    void *bidx;
    int32_t start;
    int32_t end;
    int32_t bin_start;
    int32_t bin_end;
    void *item;
    void *prev_item;
    int l;
} binidx_itr_t;
typedef binidx_itr_t bioidx_itr_t;
bioidx_t *bioidx_init();
void bioidx_destroy(bioidx_t *bioidx);
int bioidx_chrom_insert(bioidx_t *bioidx, int32_t bioidx_key, uint32_t min_shift, uint32_t step);
int bioidx_insert(bioidx_t *bioidx, int32_t bioidx_key, int32_t start, int32_t end, void *data);
int bioidx_search(bioidx_t *bioidx, bioidx_itr_t *itr, int32_t bioidx_key, int32_t start, int32_t end);
bioidx_itr_t *bioidx_itr_init();
void bioidx_itr_destroy(bioidx_itr_t *itr);
void *binidx_itr_next(void *itr);
int binidx_itr_remove(void *itr);
static inline void *bioidx_itr_next(bioidx_itr_t *itr){
    if (!itr->bidx) return NULL;
    return binidx_itr_next(itr);
}
static inline int bioidx_itr_remove(bioidx_itr_t *itr) {
    if (!itr->bidx) return -1;
    return binidx_itr_remove(itr);
}
#endif /* __BIOIDX_H */