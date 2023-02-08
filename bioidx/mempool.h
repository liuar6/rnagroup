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

#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#ifndef __MEMPOOL_H
#define __MEMPOOL_H

#define MEMPOOL_VERSION "0.1.0"

#ifndef BPOOL_DEFAULT_MAX_FREE
#define BPOOL_DEFAULT_MAX_FREE 256
#endif

#ifndef VSPOOL_DEFAULT_EXTENSION_FACTOR
#define VSPOOL_DEFAULT_EXTENSION_FACTOR 2
#endif

#ifndef VSPOOL_DEFAULT_MAX_DSIZE
#define VSPOOL_DEFAULT_MAX_DSIZE 1u<<24u
#endif

#ifndef FSPOOL_DEFAULT_EXTENSION_FACTOR
#define FSPOOL_DEFAULT_EXTENSION_FACTOR 2
#endif

typedef struct memb_s{
    void *next;
    void *prev;
    size_t size;
    size_t alloc_size;
} memb_t;

typedef struct bpool_s{
    memb_t *free;
    memb_t *free_end;
    size_t n_free;
    size_t max_free;
    int ref_count;
} bpool_t;

typedef struct vspool_s{
    bpool_t *bp;
    bpool_t *extra_bp;
    size_t max_dsize;
    void *block;
    void *extra;
    size_t block_size;
    double block_extension_factor;
    void *block_avail;
    void *block_end;
    int n_bin;
    size_t *bin_capacity;
    void **free;
    size_t n_alloc;
    size_t n_extra_alloc;
} vspool_t;

typedef struct fspool_s{
    bpool_t *bp;
    void *block;
    size_t dsize;
    size_t block_size;
    double block_extension_factor;
    void *block_avail;
    void *block_end;
    void *free;
    size_t n_alloc;
} fspool_t;

#define memb_size(s) (((memb_t *)(s)-1)->size)
#define memb_alloc_size(s) (((memb_t *)(s)-1)->alloc_size)
#define memb_detach(b, begin, end)  \
do { \
    if (b->next) ((memb_t *)b->next)->prev = b->prev; \
    else end = b->prev; \
    if (b->prev) ((memb_t *)b->prev)->next = b->next; \
    else begin = b->next; \
} while (0)
#define memb_attach(b, begin, end) \
do {    \
    b->next = begin;  \
    b->prev = NULL; \
    if (begin) ((memb_t *)begin)->prev = b; \
    else end = b;   \
    begin = b;    \
} while (0)

#define BPOOL_SET_MAX_FREE 1
static int bpool_set(bpool_t *p, int option, ...){
    va_list ap;
    switch(option){
        case BPOOL_SET_MAX_FREE:
            va_start(ap, option);
            size_t max_free = va_arg(ap, size_t);
            va_end(ap);
            p->max_free = max_free;
            return 0;
        default:
            return -1;
    }
}

bpool_t *bpool_init(){
    bpool_t *p;
    p = malloc(sizeof(*p));
    if (!p) return NULL;
    p->free = NULL;
    p->free_end = NULL;
    p->n_free = 0;
    p->max_free = BPOOL_DEFAULT_MAX_FREE;
    p->ref_count = 1;
    return p;
}

void bpool_deref(bpool_t *p){
    memb_t *b, *b1;
    p->ref_count--;
    if (p->ref_count == 0) {
        for (b = p->free; b; b1 = b, b = b->next, free(b1));
        free(p);
    }
}

void bpool_ref(bpool_t *p){
    p->ref_count++;
}

void bpool_destroy(bpool_t *p){
    bpool_deref(p);
}

void bpool_join(bpool_t *p1, bpool_t *p2){
    p1->free_end->next = p2->free;
    p2->free->prev = p1->free_end;
    p1->free_end = p2->free;
    p2->free = NULL;
    p2->free_end = NULL;
    p1->n_free += p2->n_free;
}

void bpool_shrink(bpool_t *p, size_t size){
    memb_t *b, *b1;
    size_t tsize = 0;
    for (b = p->free; b && (tsize += b->size) <= size; b = b->next);
    if (!b) return;
    if (b->prev) {
        ((memb_t*)b->prev)->next = NULL;
        p->free_end = b->prev;
    }
    while (b) {b1 = b; b = b->next; free(b1), p->n_free--;}
}

void *balloc_range(bpool_t *p, size_t size, size_t max_size){
    memb_t *b = NULL;
    for (b = p->free; b; b = b->next) if (b->size >= size && b->size <= max_size) break;
    if (!b){
        b = malloc(sizeof(*b) + size);
        if (!b) return NULL;
        b->size = size;
    } else {
        memb_detach(b, p->free, p->free_end);
        p->n_free--;
    }
    b->alloc_size = size;
    return b + 1;
}

void *balloc(bpool_t *p, size_t size){
    return balloc_range(p, size, size + (size<<1u));
}

void bdrop(void *s){
    memb_t *b = (memb_t *)s - 1;
    free(b);
}

void bfree(bpool_t *p, void *s){
    if (!p) bdrop(s);
    memb_t *b = (memb_t *)s - 1;
    memb_attach(b, p->free, p->free_end);
    p->n_free++;
    while (p->n_free > p->max_free){
        b = p->free_end;
        memb_detach(b, p->free, p->free_end);
        p->n_free--;
        free(b);
    }
}

void *brealloc_range(bpool_t *p, void *s, size_t size, size_t max_size){
    if (size <= memb_size(s)) {
        memb_alloc_size(s) = size;
        return s;
    }
    void *new_s = balloc_range(p, size, max_size);
    if (!new_s) return NULL;
    memcpy(new_s, s, memb_alloc_size(s));
    bfree(p, s);
    memb_alloc_size(new_s) = size;
    return new_s;
}

static void *brealloc(bpool_t *p, void *s, size_t size){
    return brealloc_range(p, s, size, size + (size<<1u));
}

int vsbin_power16[16] = {0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3};
static inline int vsbin_power(size_t x){
    int power = 0;
    if (x >= (1lu<<32u)) {power += 32; x >>= 32u;}
    if (x >= (1lu<<16u)) {power += 16; x >>= 16u;}
    if (x >= (1lu<<8u)) {power += 8; x >>= 8u;}
    if (x >= (1lu<<4u)) {power += 4; x >>= 4u;}
    power += vsbin_power16[x];
    return power;
}
int vsbin_void_pointer_power;
static inline int vsbin(size_t *dsize){
    size_t udsize;
    int power, bin;
    if (*dsize <= sizeof(void *)){
        *dsize = sizeof(void *);
        return 0;
    }
    power = vsbin_power(*dsize - 1) + 1;
    if (power < vsbin_void_pointer_power + 2){
        bin = power - vsbin_void_pointer_power;
        udsize = 1u<<power;
    } else {
        bin = (power - vsbin_void_pointer_power) * 2 - 1;
        udsize = 1u<<power;
        if (udsize - (udsize>>2u) >= *dsize){
            bin--;
            udsize -= udsize>>2u;
        }
    }
    *dsize = udsize;
    return bin;
}

static inline size_t vsbin_capacity(int bin){
    if (bin < 2) return 1u<<(bin + vsbin_void_pointer_power);
    size_t ret = 1u<<(bin / 2 + vsbin_void_pointer_power + 1);
    if (!(bin % 2)) ret -= ret>>2u;
    return ret;
}

#define VSPOOL_SET_MAX_DSIZE 1
#define VSPOOL_SET_BPOOL 2
#define VSPOOL_SET_EXTRA_BPOOL 3
static int vspool_set(vspool_t *p, int option, ...){
    va_list ap;
    switch(option){
        case VSPOOL_SET_MAX_DSIZE:
            va_start(ap, option);
            size_t max_dsize = va_arg(ap, size_t);
            va_end(ap);
            if (p->block || p->extra_bp) return -1;
            int n_bin = vsbin(&max_dsize) + 1;
            void *free_ = calloc(n_bin, sizeof(void *));
            if (!free_) return -1;
            void *bin_capacity = malloc(n_bin * sizeof(size_t));
            if (!bin_capacity) {free(free_); return -1;}
            for (int i = 0; i < n_bin; ++i) ((size_t *)bin_capacity)[i] = vsbin_capacity(i);
            if (p->free) free(p->free);
            if (p->bin_capacity) free(p->bin_capacity);
            p->max_dsize = max_dsize;
            p->n_bin = n_bin;
            p->free = free_;
            p->bin_capacity = bin_capacity;
            return 0;
        case VSPOOL_SET_BPOOL:
            va_start(ap, option);
            bpool_t * bp = va_arg(ap, bpool_t *);
            va_end(ap);
            if (p->bp) bpool_deref(p->bp);
            bpool_ref(bp);
            p->bp = bp;
            return 0;
        case VSPOOL_SET_EXTRA_BPOOL:
            va_start(ap, option);
            bpool_t * extra = va_arg(ap, bpool_t *);
            va_end(ap);
            if (p->extra_bp) bpool_deref(p->extra_bp);
            bpool_ref(extra);
            p->extra_bp = extra;
            return 0;
        default:
            return -1;
    }
}

static vspool_t *vspool_init(){
    vspool_t *p;
    p = malloc(sizeof(*p));
    if (!p) return NULL;
    p->bp = NULL;
    p->extra_bp = NULL;
    p->block = NULL;
    p->extra = NULL;
    p->block_size = 0;
    p->block_extension_factor = VSPOOL_DEFAULT_EXTENSION_FACTOR;
    p->block_avail = NULL;
    p->block_end = NULL;
    p->free = NULL;
    p->bin_capacity = NULL;
    if (vspool_set(p, VSPOOL_SET_MAX_DSIZE, VSPOOL_DEFAULT_MAX_DSIZE) != 0) {free(p); return NULL;}
    p->n_alloc = 0;
    p->n_extra_alloc = 0;
    vsbin_void_pointer_power = vsbin_power(sizeof(void *));
    return p;
}

static void vspool_destroy(vspool_t *p){
    if (!p) return;
    void *b, *b1;
    for (b = p->block; b; b1 = b, b = *(void **)b, (memb_size(b1)<(1u<<20u))? bdrop(b1): bfree(p->bp, b1));
    if (p->bp) bpool_deref(p->bp);
    for (b = p->extra; b; b1 = b, b = *(void **)b, (memb_size(b1)<(1u<<20u))? bdrop(b1): bfree(p->extra_bp, b1));
    if (p->extra_bp) bpool_deref(p->extra_bp);
    free(p->free);
    free(p->bin_capacity);
    free(p);
}

static void vspool_rewind(vspool_t *p){
    void *b, *b1;
    if (!p->block) return;
    for (b = *(void **)(p->block); b; b1 = b, b = *(void **)b, (memb_size(b1)<(1u<<20u))? bdrop(b1): bfree(p->bp, b1));
    *(void **)(p->block) = NULL;
    p->block_avail = p->block + sizeof(void *);
    for (int i = 0; i < p->n_bin; ++i) p->free[i] = NULL;
}

static void *vsbrk(vspool_t *p, size_t dsize){
    void *ret;
    dsize = (dsize + sizeof(void *) - 1) & ~(sizeof(void *) - 1);
    if (dsize > p->block_end - p->block_avail){
        if (!p->bp) p->bp = bpool_init();
        size_t extend;
        extend = (size_t) (p->block_size * p->block_extension_factor);
        if (extend < dsize) extend = dsize;
        extend = (extend + sizeof(void *) - 1) & ~(sizeof(void *) - 1);
        void *new_block;
        new_block = balloc(p->bp, (extend + sizeof(void *)));
        if (!new_block) return NULL;
        *(void **)new_block = p->block;
        p->block = new_block;
        p->block_size = memb_size(new_block);
        p->block_avail = new_block + sizeof(void *);
        p->block_end = p->block_avail + extend;
    }
    ret = p->block_avail;
    p->block_avail += dsize;
    return ret;
}

#define vs_extra_next(e) (*((void **)e))
#define vs_extra_prev(e) (*((void **)e + 1))
static void *vsalloc(vspool_t *p, size_t dsize){
    if (dsize > p->max_dsize){
        if (!p->extra_bp) p->extra_bp = bpool_init();
        void *new_extra = balloc(p->extra_bp, sizeof(void *[3]) + dsize);
        if (!new_extra) return NULL;
        vs_extra_next(new_extra) = p->extra;
        vs_extra_prev(new_extra) = NULL;
        if (p->extra) vs_extra_prev(p->extra) = new_extra;
        p->extra = new_extra;
        *(int *)(new_extra + sizeof(void *[2]))  = -1;
        p->n_extra_alloc++;
        return  new_extra + sizeof(void *[3]);
    }
    void *pt;
    size_t udsize = dsize;
    int bin = vsbin(&udsize);
    if (p->free[bin]){
        pt = p->free[bin];
        p->free[bin] = *(void **)pt;
        p->n_alloc++;
        return pt;
    } else {
        pt = vsbrk(p, sizeof(void *) + udsize);
        if (!pt) return NULL;
        *(int *)pt = bin;
        p->n_alloc++;
        return (void **)pt + 1;
    }
}

static void vsfree(vspool_t *p, void *s){
    int bin = *(int *)((void **)s - 1);
    if (bin < 0) {
        p->n_extra_alloc--;
        void *extra = s - sizeof(void *[3]);
        if (vs_extra_next(extra)) vs_extra_prev(vs_extra_next(extra)) = vs_extra_prev(extra);
        if (vs_extra_prev(extra)) vs_extra_next(vs_extra_prev(extra)) = vs_extra_next(extra);
        else p->extra = vs_extra_next(extra);
        bfree(p->extra_bp, extra);
    } else {
        p->n_alloc--;
        *(void **)s = p->free[bin];
        p->free[bin] = s;
    }
}

static void *vsrealloc(vspool_t *p, void *s, size_t dsize){
    void *ret;
    int bin = *(int *)((void **)s - 1);
    if (bin < 0) {
        if (dsize > p->max_dsize) {
            void *new_extra = brealloc(p->extra_bp, s - sizeof(void *[3]), sizeof(void *[3]) + dsize);
            if (vs_extra_prev(new_extra)) vs_extra_next(vs_extra_prev(new_extra)) = new_extra;
            else p->extra = new_extra;
            if (vs_extra_next(new_extra)) vs_extra_prev(vs_extra_next(new_extra)) = new_extra;
            ret = new_extra + sizeof(void *[3]);
        } else {
            ret = vsalloc(p, dsize);
            if (!ret) return NULL;
            memcpy(ret, s, dsize);
            vsfree(p, s);
        }
    } else {
        size_t capacity = p->bin_capacity[bin];
        if (dsize <= capacity) ret = s;
        else {
            ret = vsalloc(p, dsize);
            memcpy(ret, s, capacity);
            vsfree(p, s);
        }
    }
    return ret;
}

#define FSPOOL_SET_BPOOL 2
static int fspool_set(fspool_t *p, int option, ...){
    va_list ap;
    switch(option){
        case FSPOOL_SET_BPOOL:
            va_start(ap, option);
            bpool_t * bp = va_arg(ap, bpool_t *);
            va_end(ap);
            if (p->bp) bpool_deref(p->bp);
            bpool_ref(bp);
            p->bp = bp;
            return 0;
        default:
            return -1;
    }
}

static fspool_t *fspool_init(size_t dsize){
    fspool_t *p;
    p = malloc(sizeof(*p));
    if (!p) return NULL;
    p->bp = NULL;
    if (dsize == 0) dsize = sizeof(void *);
    p->dsize = (dsize + sizeof(void *) - 1) & ~(sizeof(void *) - 1);
    p->block = NULL;
    p->block_size = 0;
    p->block_extension_factor = FSPOOL_DEFAULT_EXTENSION_FACTOR;
    p->block_avail = NULL;
    p->block_end = NULL;
    p->free = NULL;
    p->n_alloc = 0;
    return p;
}

static void fspool_destroy(fspool_t *p){
    if (!p) return;
    void *b, *b1;
    for (b = p->block; b; b1 = b, b = *(void **)b, (memb_size(b1)<(1u<<20u)) ? bdrop(b1): bfree(p->bp, b1));
    if (p->bp) bpool_deref(p->bp);
    free(p);
}

static void *fsalloc(fspool_t *p){
    void *ret;
    if (p->free){
        ret = p->free;
        p->free = *(void **)ret;
        p->n_alloc++;
        return ret;
    }
    if (p->block_avail == p->block_end){
        if (!p->bp) p->bp = bpool_init();
        size_t extend;
        extend = (size_t) (p->block_size * p->block_extension_factor);
        if (extend < p->dsize) extend = p->dsize;
        extend = (extend + sizeof(void *) - 1) & ~(sizeof(void *) - 1);
        void *new_block;
        new_block = balloc(p->bp, (extend + sizeof(void *)));
        if (!new_block) return NULL;
        *(void **)new_block = p->block;
        p->block = new_block;
        p->block_size = memb_size(new_block);
        p->block_avail = new_block + sizeof(void *);
        p->block_end = p->block_avail + extend - (extend % p->dsize);
    }
    ret = p->block_avail;
    p->block_avail += p->dsize;
    p->n_alloc++;
    return ret;
}

static void fsfree(fspool_t *p, void *s){
    *(void **)s = p->free;
    p->free = s;
}
#endif /*__MEMPOOL_H */
