#include "common.h"
#include "file.h"
#include "htslib/sam.h"

static char ** strsplit(char * line, char ** results, int length, char c){
    char *start=line;
    char *end=NULL;
    int i=0;
    while ((end=strchr(start, c))!=NULL && i<length){
        end[0]='\0';
        results[i]=start;
        start=end+1;
        i=i+1;
    }
    if (i<length && start[0]!='\0') {
        results[i]=start;
        i=i+1;
    }
    for (;i<length;++i) results[i]=NULL;
    return results;
}

void bam_parser_process(bam_parser_t *bp, khash_t (record) *rh, vec_t (record) *rv){
    bam1_t *b1, *b2;
    int32_t tid;
    int8_t strand;
    int32_t pos;
    for (int i = 0; i < bp->read_index; ++i) {
        if (bp->read_valid[i] == 0) continue;
        b1 = bp->read1[i];
        b2 = bp->read2[i];

        tid = b1->core.tid;
        strand = (b1->core.flag & BAM_FREVERSE) ? '-' : '+';
        pos = (b1->core.flag & BAM_FREVERSE) ? bam_endpos(b1) - 1 : b1->core.pos;
        khint_t key_i = kh_get(record, rh, stidpos(tid, strand, pos));
        if (key_i == kh_end(rh)) {
            int ret;
            record_t *new_record = malloc(sizeof(record_t));
            new_record->tid = tid;
            new_record->strand = strand;
            new_record->start = malloc(sizeof(*new_record->start));
            new_record->start->pos = pos;
            new_record->start->count = 0;
            new_record->start->norm_count = 0;
            new_record->end = vec_init(position);
            key_i = kh_put(record, rh, stidpos(tid, strand, pos), &ret);
            kh_val(rh, key_i) = new_record;
            vec_add(record, rv, new_record);
        }
        record_t *r = kh_val(rh, key_i);
        r->start->count++;
        r->start->norm_count += 1.0 / bp->valid_count;

        int32_t end = (r->strand == '-') ? b2->core.pos : bam_endpos(b2) - 1;
        position_t *p = NULL;
        for (int j = 0; j < r->end->size; ++j)
            if (end == r->end->data[j].pos) p = r->end->data + j;
        if (!p){
            position_t new_pos;
            new_pos.pos = end;
            new_pos.count = 0;
            new_pos.norm_count = 0;
            vec_add(position, r->end, new_pos);
            p = r->end->data + r->end->size - 1;
        }
        p->count++;
        p->norm_count += 1.0 / bp->valid_count;
    }
}

void rg_read(const char* input, vec_t(rg) *rg, khash_t(str2int32) *chrom2id){
    FILE *fp = fopen(input, "r");
    char *buffer = malloc(1u<<20u);
    char *field[7];
    rg_t *new_rg;
    while (fgets(buffer, 1u<<20u, fp)){
        if (buffer[0] == '#') continue;
        strsplit(buffer, field, 7, '\t');
        new_rg = calloc(1, sizeof(*new_rg));
        new_rg->tid = kh_val(chrom2id, kh_get(str2int32, chrom2id, field[0]));
        new_rg->strand = field[1][0];
        new_rg->fstart = strtol(field[2], NULL, 10);
        new_rg->fend = strtol(field[3], NULL, 10);
        new_rg->tstart = strtol(field[4], NULL, 10);
        new_rg->tend = strtol(field[5], NULL, 10);
        vec_add(rg, rg, new_rg);
    }
    free(buffer);
}

void rg_count(vec_t(rg) *rg, vec_t(record) *rv, int n_targets, uint32_t *target_len){
    int i, j;
    khiter_t k;
    int32_t tid;
    int8_t strand;
    int32_t region_start, region_end, start, end;
    double count, norm_count;
    rg_t *t;
    position_t *p;
    bioidx_t *bidx1 = bioidx_init();
    bioidx_t *bidx2 = bioidx_init();
    bioidx_itr_t *bitr = calloc(1, sizeof(*bitr));
    for (i = 0; i < rg->size; ++i){
        t = rg->data[i];
        t->start = vec_init(position);
        t->end = vec_init(position);
        bioidx_insert(bidx1, bioidx_key(t->tid, t->strand), t->fstart, t->fend+1, t);
        bioidx_insert(bidx2, bioidx_key(t->tid, t->strand), (t->strand=='+')?t->fstart:t->tstart,((t->strand=='+')?t->tend:t->fend)+1, t);
    }

    for (i = 0; i < rv->size; ++i){
        record_t *r = rv->data[i];
        tid = r->tid;
        strand = r->strand;
        start = r->start->pos;
        for (j = 0; j < r->end->size; ++j){
            p = r->end->data + j;
            end = p->pos;
            count = p->count;
            norm_count = p->norm_count;
            region_start = (r->strand == '+')?start:end;
            region_end = (r->strand == '+')?(end+1):(start+1);
            bioidx_search(bidx1, bitr, bioidx_key(tid, strand), start, start+1);
            while ((t = bioidx_itr_next(bitr)) != NULL){
                if (end >= t->tstart && end <= t->tend){
                    t->count+=count;
                    t->norm_count += norm_count;
                    for (k = 0; k < t->start->size; ++k)
                        if ((p = t->start->data + k)->pos == start){
                            p->count+=count;
                            p->norm_count += norm_count;
                            break;
                        }
                    if (k == t->start->size){
                        position_t new_pos;
                        new_pos.pos = start;
                        new_pos.count = count;
                        new_pos.norm_count = norm_count;
                        vec_add(position, t->start, new_pos);
                    }
                    for (k = 0; k < t->end->size; ++k)
                        if ((p = t->end->data + k)->pos == end){
                            p->count+=count;
                            p->norm_count += norm_count;
                            break;
                        }
                    if (k == t->end->size){
                        position_t new_pos;
                        new_pos.pos = end;
                        new_pos.count = count;
                        new_pos.norm_count = norm_count;
                        vec_add(position, t->end, new_pos);
                    }
                }
            }
            bioidx_search(bidx2, bitr, bioidx_key(tid, strand), region_start, region_end);
            while ((t = bioidx_itr_next(bitr)) != NULL){
                t->bg_count+=count;
                t->bg_norm_count+=norm_count;
            }

        }
    }
    bioidx_destroy(bidx1);
    bioidx_destroy(bidx2);
    free(bitr);
}

void rg_output(const char*output, vec_t(rg) *rg, char **id2chrom, const char *comment){
    int i, j;
    rg_t *t;
    FILE *out = fopen(output, "w");
    if (!out) fprintf(stderr, "[rg_output] Error: unable to open the output file\n\n");
    if (comment) fprintf(out, "#%s\n", comment);
    for (i = 0; i < rg->size; ++i){
        t = rg->data[i];
        fprintf(out, "%s\t%c\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t", id2chrom[t->tid], t->strand, t->fstart, t->fend, t->tstart, t->tend, t->count, t->norm_count, t->bg_count, t->bg_norm_count);

        if (t->start == NULL) {
            fprintf(out, "\n");
            continue;
        } else if (t->count == 0){
            fprintf(out, "\t\t\t\t\t\n");
            continue;
        }

        qsort(t->start->data, t->start->size, sizeof(position_t), position_compare);
        qsort(t->end->data, t->end->size, sizeof(position_t), position_compare);
        fprintf(out, "%d", (t->start->data + 0)->pos);
        for (j = 1; j < t->start->size; ++j) fprintf(out, ",%d", (t->start->data + j)->pos);
        fprintf(out, "\t");
        fprintf(out, "%f", (t->start->data + 0)->count);
        for (j = 1; j < t->start->size; ++j) fprintf(out, ",%f", (t->start->data + j)->count);
        fprintf(out, "\t");
        fprintf(out, "%f", (t->start->data + 0)->norm_count);
        for (j = 1; j < t->start->size; ++j) fprintf(out, ",%f", (t->start->data + j)->norm_count);
        fprintf(out, "\t");
        fprintf(out, "%d", (t->end->data + 0)->pos);
        for (j = 1; j < t->end->size; ++j) fprintf(out, ",%d", (t->end->data + j)->pos);
        fprintf(out, "\t");
        fprintf(out, "%f", (t->end->data + 0)->count);
        for (j = 1; j < t->end->size; ++j) fprintf(out, ",%f", (t->end->data + j)->count);
        fprintf(out, "\t");
        fprintf(out, "%f", (t->end->data + 0)->norm_count);
        for (j = 1; j < t->end->size; ++j) fprintf(out, ",%f", (t->end->data + j)->norm_count);
        fprintf(out, "\n");
    }
    fclose(out);
}

void rg_free(vec_t(rg) *rg){
    for (int i = 0; i < rg->size; ++i){
        rg_t *t = rg->data[i];
        if (t->start) vec_destroy(position, t->start);
        if (t->end) vec_destroy(position, t->end);
        free(t);
    }
    vec_destroy(rg, rg);
}