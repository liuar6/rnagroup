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

#include "khash.h"

typedef struct bin_item_t{
    struct bin_item_t *next;
    int32_t start;
    int32_t end;
    void *data;
} bin_item_t;

typedef struct bin_t{
    bin_item_t *item;
} bin_t;

KHASH_MAP_INIT_INT(bin, bin_t*)
typedef struct binidx_t{
    uint32_t min_shift;
    uint32_t step;
    uint32_t n_level;
    khash_t(bin) **bh;
    struct fspool_s *bp;
    struct fspool_s *bip;
} binidx_t;

typedef struct binidx_itr_t{
    binidx_t *bidx;
    int32_t start;
    int32_t end;
    int32_t bin_start;
    int32_t bin_end;
    bin_item_t *item;
    bin_item_t *prev_item;
    int l;
} binidx_itr_t;

void *binidx_init(uint32_t min_shift, uint32_t step);
void binidx_destroy(void *_bidx);
int binidx_insert(void *_bidx, int32_t start, int32_t end, void *data);
int binidx_search(void *_bidx, void *_itr, int32_t start, int32_t end);
void *binidx_itr_next(void *_itr);
int binidx_itr_remove(void *_itr);