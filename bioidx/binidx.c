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
#include "mempool.h"
#include "binidx.h"

#define min(a, b) (((a)>(b))?(b):(a))
#define max(a, b) (((a)>(b))?(a):(b))
#define RANGE_INTERSECT(start1, end1, start2, end2) (min(end1, end2)-max(start1, start2) > 0)

static inline void reg2bin(uint32_t beg, uint32_t end, uint32_t min_shift, uint32_t step, uint32_t *level, int *bin){
    uint32_t s = min_shift, l = 0;
    end--;
    while (beg>>s < end>>s){s += step; l++;}
    *level = l;
    *bin = beg>>s;
}

void *binidx_init(uint32_t min_shift, uint32_t step){
    binidx_t *bidx;
    bidx = malloc(sizeof(*bidx));
    if (!bidx) return NULL;
    bidx->min_shift = min_shift;
    bidx->step = step;
    bidx-> n_level = 0;
    bidx->bh = NULL;
    bidx->bp = fspool_init(sizeof(bin_t));
    if (!bidx->bp) {free(bidx); return NULL;}
    bidx->bip = fspool_init(sizeof(bin_item_t));
    if (!bidx->bip) {free(bidx); fspool_destroy(bidx->bp); return NULL;}
    return bidx;
}

void binidx_destroy(void *_bidx){
    binidx_t * bidx = _bidx;
    khash_t(bin) **bh = bidx->bh;
    khash_t(bin) *h;
    bin_t *b;
    bin_item_t *item;
    khiter_t k;
    int i;
    for (i = 0; i < bidx->n_level; ++i){
        if (!(h = bh[i])) continue;
        for (k = kh_begin(h); k != kh_end(h); ++k)
            if (kh_exist(h, k)){
                b = kh_val(h, k);
                while (b->item){
                    item = b->item;
                    b->item = item->next;
                    fsfree(bidx->bp, item);
                }
                fsfree(bidx->bip, b);
            }
        kh_destroy(bin, h);
    }
    free(bh);
    fspool_destroy(bidx->bp);
    fspool_destroy(bidx->bip);
    free(bidx);
}

static int bin_insert(bin_t *b, bin_item_t *new_item, int32_t start, int32_t end, void *data){
    if (!new_item) return -1;
    new_item->start = start;
    new_item->end = end;
    new_item->data = data;
    new_item->next = b->item;
    b->item = new_item;
    return 0;
}

int binidx_insert(void *_bidx, int32_t start, int32_t end, void *data){
    binidx_t *bidx = _bidx;
    khash_t(bin) *h;
    bin_t *b;
    uint32_t level;
    int bin, ret;
    khint_t k;
    if (start < 0 || end <= start) return -1;
    reg2bin(start, end, bidx->min_shift, bidx->step, &level, &bin);
    if (level >= bidx->n_level) {
        khash_t(bin) **new_bh;
        new_bh = realloc(bidx->bh, (level + 1) * sizeof(*bidx->bh));
        if (!new_bh) return -1;
        for (int i = bidx->n_level; i <= level; ++i) new_bh[i] = NULL;
        bidx->bh = new_bh;
        bidx->n_level = level + 1;
    }
    h = bidx->bh[level];
    if (!h){
        if (!(h = kh_init(bin))) return -1;
        bidx->bh[level] = h;
    }
    if ((k = kh_get(bin, h, bin)) == kh_end(h)){
        b = fsalloc(bidx->bp);
        if (!b) return -1;
        b->item = NULL;
        k = kh_put(bin, h, bin, &ret);
        if (ret < 0) {fsfree(bidx->bp, b); return -1;}
        kh_val(h, k) = b;
    } else b = kh_val(h, k);
    return bin_insert(b, fsalloc(bidx->bip), start, end, data);
}

int binidx_search(void *_bidx, void *_itr, int32_t start, int32_t end){
    if (start < 0 || end <= start) return -1;
    binidx_itr_t *itr = _itr;
    itr->bidx = _bidx;
    itr->start = start;
    itr->end = end;
    itr->bin_start = 0;
    itr->bin_end = 0;
    itr->l = -1;
    itr->item = NULL;
    return 0;
}

void *binidx_itr_next(void *_itr) {
    binidx_itr_t *itr = _itr;
    do {
        if (!(itr->prev_item = itr->item) || !(itr->item = itr->prev_item->next)) {
            khash_t (bin) *h = itr->bidx->bh[itr->l];
            khiter_t k;
            do {
                itr->bin_start++;
                if (itr->bin_start > itr->bin_end) {
                    binidx_t *bidx = itr->bidx;
                    khash_t (bin) **bh = bidx->bh;
                    uint32_t min_shift = bidx->min_shift;
                    uint32_t step = bidx->step;
                    uint32_t n_level = bidx->n_level;
                    int new_l;
                    for (new_l = itr->l + 1; new_l < n_level && !(h = bh[new_l]); new_l++);
                    if (new_l == n_level) { itr->bin_start--; return NULL;}
                    itr->l = new_l;
                    itr->bin_start = (uint32_t)(itr->start) >> (min_shift + step * new_l);
                    itr->bin_end = (uint32_t)(itr->end - 1) >> (min_shift + step * new_l);
                }
            } while ((k = kh_get(bin, h, itr->bin_start)) == kh_end(h));
            itr->item = kh_val(h, k)->item;;
        }
    } while (!RANGE_INTERSECT(itr->start, itr->end, itr->item->start, itr->item->end));
    return itr->item->data;
}

int binidx_itr_remove(void *_itr){
    binidx_itr_t *itr = _itr;
    binidx_t *bidx = itr->bidx;
    bin_item_t *item = itr->item;
    /* item == NULL only happen when the binidx_itr is just initialized or when the first item of a bin is removed */
    /* item == itr->prev_item only happen when the non-first item of a bin is removed */
    if (item == NULL || item == itr->prev_item) return -1;
    /* prev_item == NULL only happen after call of binidx_search_next after initialization or frist item of bin removed */
    if (itr->prev_item == NULL || itr->prev_item->next == NULL){ /* first item of the bin */
        khash_t (bin) *h;
        khiter_t k;
        bin_t *b;
        h = bidx->bh[itr->l];
        k =kh_get(bin, h, itr->bin_start);
        b = kh_val(h, k);
        if (item->next == NULL){ /* bin */
            fsfree(bidx->bp, b);
            kh_del(bin, h, k);
        } else b->item = item->next;
        fsfree(bidx->bip, item);
        /*upon next call of binidx_next, this bin will be rechecked */
        itr->item = NULL;
        itr->bin_start--;
    } else {    /* the prev_item is the previous item in the bin */
        /*upon next call of binidx next, the item will goes to the next item of the current item*/
        itr->prev_item->next = item->next;
        itr->item = itr->prev_item;
        fsfree(bidx->bip, item);
    }
    return 0;
}

























