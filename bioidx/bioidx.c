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

#include <stdint.h>
#include "khash.h"
#include "binidx.h"

KHASH_MAP_INIT_INT(idx, void *)
typedef struct bioidx_t{
    khash_t(idx) *idx;
    uint32_t min_shift;
    uint32_t step;
} bioidx_t;

typedef binidx_itr_t bioidx_itr_t;

bioidx_t *bioidx_init(){
    bioidx_t *bioidx;
    bioidx = malloc(sizeof(*bioidx));
    if (!bioidx) return NULL;
    bioidx->idx = kh_init(idx);
    if (!bioidx->idx) {free(bioidx); return NULL;}
    bioidx->min_shift = 12;
    bioidx->step = 3;
    return bioidx;
}

void bioidx_destroy(bioidx_t *bioidx){
    khiter_t k;
    khash_t (idx) *h = bioidx->idx;
    for (k = kh_begin(h); k != kh_end(h); ++k)
        if (kh_exist(h, k))
            binidx_destroy(kh_val(h, k));
    kh_destroy(idx, h);
    free(bioidx);
}

int bioidx_chrom_insert(bioidx_t *bioidx, int32_t bioidx_key, uint32_t min_shift, uint32_t step){
    khiter_t k;
    int ret;
    binidx_t *binidx = binidx_init(min_shift, step);
    if (!binidx) return -1;
    k = kh_put(idx, bioidx->idx, bioidx_key, &ret);
    if (ret < 1) {binidx_destroy(binidx); return -1;} /*operation failed or key already present*/
    kh_val(bioidx->idx, k) = binidx;
    return 0;
}

int bioidx_insert(bioidx_t *bioidx, int32_t bioidx_key, int32_t start, int32_t end, void *data){
    khiter_t k;
    k = kh_get(idx, bioidx->idx, bioidx_key);
    if (k == kh_end(bioidx->idx)) {
        if (bioidx_key == -1 || bioidx_chrom_insert(bioidx, bioidx_key, bioidx->min_shift, bioidx->step) < 0) return -1;
        k = kh_get(idx, bioidx->idx, bioidx_key);
    }
    return binidx_insert(kh_val(bioidx->idx, k), start, end, data);
}

int bioidx_search(bioidx_t *bioidx, bioidx_itr_t *itr, int32_t bioidx_key, int32_t start, int32_t end){
    khiter_t k;
    k = kh_get(idx, bioidx->idx, bioidx_key);
    if (k == kh_end(bioidx->idx)){
        itr->bidx = NULL;
        return 0;
    }
    return binidx_search(kh_val(bioidx->idx, k), itr, start, end);
}

bioidx_itr_t *bioidx_itr_init(){
    return calloc(1, sizeof(bioidx_itr_t));
}

void bioidx_itr_destroy(bioidx_itr_t *itr){
    free(itr);
}
