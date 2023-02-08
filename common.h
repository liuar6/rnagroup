#include <stdio.h>
#include <string.h>
#include "bioidx/bioidx.h"

#include "htslib/khash.h"
#include "vector.h"
#include "biomacro.h"

typedef struct position_t{
    int32_t pos;
    double count;
    double norm_count;
} position_t;
static int position_compare(const void *a, const void *b){
    position_t *r1=(position_t *)a;
    position_t *r2=(position_t *)b;
    return r1->pos - r2->pos;
}
static int position_compare1(const void *a, const void *b){
    position_t *r1=*(position_t **)a;
    position_t *r2=*(position_t **)b;
    return r2->pos - r1->pos;
}


VEC_INIT(position, position_t)
typedef struct record_t{
    int32_t tid;
    int8_t strand;
    position_t *start;
    vec_t (position) *end;
} record_t;

static int record_compare(const void *a, const void *b){
    int ret = 0;
    record_t *r1=*(record_t **)a;
    record_t *r2=*(record_t **)b;
    if (r1->tid != r2->tid) ret = ((int) r1->tid) - ((int) r2->tid);
    else if (r1->strand != r2->strand) ret = r1->strand - r2->strand;
    else if (r1->strand == '+' && r1->start->pos != r2->start->pos) ret = ((int) r1->start->pos) - ((int) r2->start->pos);
    else if (r1->strand == '-' && r1->start->pos != r2->start->pos) ret = ((int) r2->start->pos) - ((int) r1->start->pos);
    return ret;
}

KHASH_MAP_INIT_INT64(record, record_t*)
VEC_INIT(record, record_t*)


typedef struct transcript_cluster{
    int32_t tid;
    char strand;
    int32_t fstart; /* zero-based */
    int32_t fend; /* zero-based */
    int32_t tstart;
    int32_t tend;
    vec_t (position) *start;
    vec_t (position) *end;
    double count;
    double norm_count;
    double bg_count;
    double bg_norm_count;
} transcript_cluster;

typedef transcript_cluster rg_t;

VEC_INIT(rg, rg_t*)
KHASH_MAP_INIT_STR(str2int32, int32_t)
void rg_read(const char* input, vec_t(rg) *rg, khash_t(str2int32) *chrom2id);
void rg_count(vec_t(rg) *rg, vec_t(record) *rv, int n_targets, uint32_t *target_len);
void rg_free(vec_t(rg) *rg);
void rg_output(const char*output, vec_t(rg) *rg, char **id2chrom, const char *comment);
