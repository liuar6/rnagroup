#include <getopt.h>
#include <math.h>

#include "htslib/khash.h"
#include "htslib/bgzf.h"
#include "htslib/sam.h"

#include "vector.h"

#include "file.h"
#include "common.h"

#define RGFIND_VERSION "1.0.0"

VEC_INIT(tc, transcript_cluster*);

void infer_rna_group(record_t **input, int size, double shrink_prop, vec_t (tc) *tc);

struct rgfind_option{
    char **input;
    int input_count;
    char *output;
    double shrink_prop;
};


void rgfind_usage(const char *msg){
    const char *usage_info = "\
rgfind: identify the transcript clusters from alignments.\n\
Usage:  rgfind [options] --fi <alignment file> ... --fo <output file>\n\
[options]\n\
-i/--fi             : input bam file (should be unsorted). [required]\n\
-o/--fo             : output rg file. [required]\n\
-h/--help           : show help informations.\n\
--shrink            : the extent ([0, 1)) of trimming the begin range and end range. default: 0.05.\n\n";
    if (msg==NULL || msg[0] == '\0') fprintf(stderr, "%s", usage_info);
    else fprintf(stderr, "%s\n\n%s", msg, usage_info);
    exit(1);
}

void rgfind_version(){
    fprintf(stderr, "rgfind-%s\n\n", RGFIND_VERSION);
}
void rgfind_option(struct rgfind_option *options, int argc, char *argv[]){
    char c;
    options->input = NULL;
    options->input_count = 0;
    options->output = NULL;
    options->shrink_prop = 0.05;
    if (argc == 1) rgfind_usage("");
    const char *short_options = "hvs:i:o:";
    const struct option long_options[] =
            {
                    { "help" , no_argument , NULL, 'h' },
                    { "version" , no_argument , NULL, 'v' },
                    { "fi" , required_argument , NULL, 'i' },
                    { "fo" , required_argument, NULL, 'o' },
                    { "shrink" , required_argument, NULL, 's' },
                    { NULL, 0, NULL, 0} ,
            };

    while ((c = getopt_long(argc, argv, short_options, long_options, NULL)) >= 0)
    {
        switch (c)
        {
            case 'h':
                rgfind_usage("");
                exit(0);
            case 'v':
                rgfind_version();
                exit(0);
            case 'o':
                options->output = optarg;
                break;
            case 'i':
                options->input = &argv[optind-1];
                options->input_count = 1;
                while (optind < argc && argv[optind][0] != '-')  {
                    options->input_count++;
                    optind++;
                }
                break;
            case 's':
                options->shrink_prop = strtof(optarg, NULL);;
                break;
            default:
                rgfind_usage("[rgfind] Error:unrecognized parameter");
        }
    }
    if (argc != optind) rgfind_usage("[rgfind] Error:unrecognized parameter");
};

int main(int argc, char **argv) {
    struct rgfind_option options;
    rgfind_option(&options, argc, argv);
    bam_parser_t *bp = NULL;

    khash_t (record) *rh = kh_init(record);
    vec_t (record) *rv = vec_init(record);

    double n_valid_read_pair = 0;

    for (int fidx = 0; fidx < options.input_count; ++fidx) {
        if (bp) bam_parser_close(bp);
        bp = bam_parser_open(options.input[fidx]);
        if (!bp) fprintf(stderr, "[rgfind] Error: unable to open the bam file\n\n");
        while (!bam_parser_read(bp)) {
            bam_parser_process(bp, rh, rv);
        }
        n_valid_read_pair+=bp->n_valid_read_pair;
    }

    /* statistics for debugging */
    double count1 = 0, count2 = 0;
    double norm1 = 0, norm2 = 0;
    for (int i = 0; i < rv->size; ++i) {
        record_t *r =rv->data[i];
        count1+=r->start->count;
        norm1+=r->start->norm_count;
        for (int j = 0 ; j < r->end->size; ++j){
            count2+=(r->end->data+j)->count;
            norm2+=(r->end->data+j)->norm_count;
        }
    }
    //fprintf(stderr, "%f\t%f\t%f\t%f\n", count1, count2, norm1, norm2);
    //fprintf(stderr, "%f\n", n_valid_read_pair);

    int i = 0, j = 0;
    int index_start, index_end;
    int loci_count = 0;
    qsort(rv->data, rv->size, sizeof(record_t *), record_compare);
    vec_t (tc) *tc = vec_init(tc);
    while (i < rv->size){
        index_start = i;
        record_t *record_start = rv->data[index_start];
        int32_t range_start = record_start->start->pos;
        int32_t range_end = record_start->end->data[0].pos;
        double range_count = 0;
        while (1){
            if (i >= rv->size) break;
            record_t *r1 =  rv->data[i];
            if (r1->tid != record_start->tid) break;
            if (r1->strand != record_start->strand) break;
            if (r1->strand == '+' && r1->start->pos > range_end) break;
            if (r1->strand == '-' && r1->start->pos < range_end) break;
            for (j = 0; j < r1->end->size; ++j){
                position_t *p = r1->end->data + j;
                if (r1->strand == '+' && p->pos > range_end) range_end = p->pos;
                if (r1->strand == '-' && p->pos < range_end) range_end = p->pos;
                range_count+=p->count;
            }
            ++i;
        }
        if (range_count >= 10) {
            //fprintf(stderr, "%s\t%d\t%d\n", bp->hdr->target_name[record_start->tid], range_start, range_end);
            infer_rna_group(rv->data + index_start, i - index_start, options.shrink_prop, tc);
        }
    }

    vec_clear(record, rv);
    kh_clear(record, rh);
    n_valid_read_pair = 0;
    for (int fidx = 0; fidx < options.input_count; ++fidx) {
        if (bp) bam_parser_close(bp);
        bp = bam_parser_open(options.input[fidx]);
        while (!bam_parser_read(bp)) {
            bam_parser_process(bp, rh, rv);
        }
        n_valid_read_pair+=bp->n_valid_read_pair;
    }

    rg_count((vec_t (rg) *)tc, rv, bp->hdr->n_targets, bp->hdr->target_len);

    char comment[100];
    sprintf(comment, "%f", n_valid_read_pair);
    rg_output(options.output, (vec_t(rg) *)tc, bp->hdr->target_name, comment);
    rg_free((vec_t (rg) *)tc);

    bam_parser_close(bp);

    vec_destroy(record, rv);
    kh_destroy(record, rh);
    return 0;
}

void infer_position_cluster(int32_t *pos, double *count, int size, int *cluster, int *cluster_count){
    double sd_fold = 3;
    double sd_extend = 1;
    double max_count;
    int index_start = 0 , last_index_start = 0;
    int index_end = 0, last_index_end = 0;
    int i = 0, j = 0;
    double current_count;
    double mean;
    double sd;
    double total_count;
    *cluster_count = 0;

    while (1){
        max_count = 0;
        i = 0;
        while(1) {
            current_count = 0;
            for (j = i; j < size && abs(pos[j] - pos[i]) < 7 && count[j] != -1; ++j)
                current_count += count[j];
            if (current_count > max_count) {
                max_count = current_count;
                index_start = i;
                index_end = j;
            }
            if (abs(pos[i++] - pos[size - 1]) < 7) break;
        }

        if (max_count < 10) break;

        while(1) {

            last_index_start = index_start;
            last_index_end = index_end;

            total_count = 0;
            mean = 0;
            sd = 0;
            for (i = index_start; i < index_end; ++i) {
                mean += pos[i] * count[i];
                total_count += count[i];
            }
            mean /= total_count;
            for (i = index_start; i < index_end; ++i) {
                sd += pow(pos[i] - mean, 2) * count[i];
            }
            if (total_count > 1) {
                sd /= (total_count - 1);
                sd = pow(sd, 0.5);
            } else sd = 0;

            while (1) {
                if (fabs(pos[index_start] - mean) <= sd_fold * sd + sd_extend) {
                    if (index_start == 0) break;
                    else if (count[index_start - 1] == -1) break;
                    else if (fabs(pos[index_start - 1] - mean) <= sd_fold * sd + sd_extend) index_start--;
                    else if (fabs(pos[index_start - 1] - mean) > sd_fold * sd + sd_extend) break;
                } else index_start++;
            }
            while(1){ /* check here please */
                if (fabs(pos[index_end - 1] - mean) <= sd_fold * sd + sd_extend) {
                    if (index_end == size) break;
                    else if (count[index_end] == -1) break;
                    else if (fabs(pos[index_end] - mean) > sd_fold * sd + sd_extend) break;
                    else if (fabs(pos[index_end] - mean) <= sd_fold * sd + sd_extend) index_end++;
                } else index_end--;
            }
            if (index_start == last_index_start && index_end == last_index_end) break;
        }
        (*cluster_count)++;
        for (i = index_start; i < index_end; ++i) {
            count[i] = -1;
            cluster[i] = *cluster_count;
        }
    }
};




KHASH_INIT(db, int32_t, double, 1, kh_int_hash_func, kh_int_hash_equal);
KHASH_INIT(int32, int32_t, int, 1, kh_int_hash_func, kh_int_hash_equal);

int int32_compare(const void *a, const void *b){
    return *(int32_t *)b - *(int32_t *)a;
}

void *shrink_rna_group(int32_t *pos, double *count, int32_t *shrink_index, int32_t *shrink_size, int32_t size, double drop){
    double curr_count, shrink_count = 0;
    int32_t curr_range, shrink_range = INT32_MAX;
    double total_count = 0;
    *shrink_index = 0, *shrink_size = 0;
    int i, j, s;
    for (i = 0; i < size; ++i) total_count+=count[i];
    for (s = size; s >= 1; --s) {
        for (i = 0; i + s <= size; ++i) {
            curr_count = 0;
            for (j = i; j < i + s; ++j) curr_count += count[j];
            curr_range = abs(pos[i+s-1]-pos[i])+1;
            if (total_count - curr_count <= drop &&  curr_range < shrink_range){
                *shrink_index = i;
                *shrink_size = s;
                shrink_range = curr_range;
            }
        }
    }
}

void infer_rna_group(record_t **input, int size, double shrink_prop, vec_t (tc) *tc){
    int32_t *pos, *pos2;
    double *count, *count2;
    int *cluster, *cluster2;
    int cluster_count, member_count;
    int cluster_count2, member_count2;
    record_t **r = calloc(sizeof(*r), size);

    int32_t key_i;
    int i, j, k, c1, c2, ret;

    kh_db_t *kh = kh_init(db);
    khash_t (int32) *kh1 = kh_init(int32);
    khash_t (int32) *kh2 = kh_init(int32);
    khash_t (int32) *kh1s = kh_init(int32);
    khash_t (int32) *kh2s = kh_init(int32);

    pos = calloc(sizeof(*pos), size);
    count = calloc(sizeof(count), size);
    cluster = calloc(sizeof(*cluster), size);
    cluster_count = 0;

    for (i = 0; i < size; ++i){
        pos[i] = input[i]->start->pos;
        count[i] = input[i]->start->count;
    }

    infer_position_cluster(pos, count, size, cluster, &cluster_count);
    for (c1 = 1; c1 <= cluster_count; ++c1){
        kh_clear(int32, kh1);
        member_count = 0;
        for (i = 0; i < size; ++i)
            if (cluster[i] == c1) {
                r[member_count] = input[i];
                key_i = kh_put(int32, kh1, pos[i], &ret);
                kh_val(kh1, key_i) = member_count;
                member_count++;
            }
        kh_clear(db, kh);
        int size2 = 0;
        for (i = 0; i < member_count; ++i){
            for (j = 0; j < r[i]->end->size; ++j){
                position_t *p = r[i]->end->data + j;
                key_i=kh_get(db, kh, p->pos);
                if (key_i == kh_end(kh)){
                    key_i = kh_put(db, kh, p->pos, &ret);
                    kh_val(kh, key_i) = 0;
                    size2++;
                }
                kh_val(kh, key_i)+=p->count;
            }
        }

        pos2 = calloc(size2, sizeof(*pos2));
        count2 = calloc(size2, sizeof(*count2));
        cluster2 = calloc(size2, sizeof(int));
        cluster_count2 = 0;

        for (i = kh_begin(kh), j = 0; i != kh_end(kh); ++i) if (kh_exist(kh, i)) pos2[j++] = kh_key(kh, i);
        qsort(pos2, size2, sizeof(*pos2), int32_compare);
        for (i = 0; i < size2; ++i){
            key_i = kh_get(db, kh, pos2[i]);
            count2[i] = kh_val(kh, key_i);
        }

        infer_position_cluster(pos2, count2, size2, cluster2, &cluster_count2);

        for (c2 = 1; c2 <= cluster_count2; ++c2){
            kh_clear(int32, kh2);
            member_count2 = 0;
            for (i = 0; i < size2; ++i)
                if (cluster2[i] == c2){
                    key_i = kh_put(int32, kh2, pos2[i], &ret);
                    kh_val(kh2, key_i) = member_count2;
                    member_count2++;
                }

            kh_clear(int32, kh1s);
            kh_clear(int32, kh2s);

            int32_t *c1_pos = calloc(member_count, sizeof(int32_t));
            double *c1_count = calloc(member_count, sizeof(double));
            khiter_t key1, key2;
            int32_t shrink_index, shrink_size, shrink_index2, shrink_size2;
            position_t *start, *end;
            double total_count = 0, total_count2 = 0;
            int32_t member_count_s = 0, member_count2_s = 0;


            for (i = 0; i < size; ++i){
                start = input[i]->start;
                if ((key1 = kh_get(int32, kh1, start->pos))  != kh_end(kh1)){
                    c1_pos[kh_val(kh1, key1)] = start->pos;
                    for (j = 0; j < input[i]->end->size; ++j){
                        end = input[i]->end->data + j;
                        if ((key2 = kh_get(int32, kh2, end->pos)) != kh_end(kh2))
                            c1_count[kh_val(kh1, key1)] += end->count;
                    }
                }
            }
            for (i = 0; i < member_count; ++i) total_count+=c1_count[i];
            shrink_rna_group(c1_pos, c1_count, &shrink_index, &shrink_size, member_count, total_count * shrink_prop);
            for (i = shrink_index; i < shrink_index + shrink_size; ++i){
                key_i = kh_put(int32, kh1s, c1_pos[i], &ret);
                kh_val(kh1s, key_i) = member_count_s++;
            }

            int32_t *c2_pos = calloc(member_count2, sizeof(int32_t));
            double *c2_count = calloc(member_count2, sizeof(double));
            for (i = kh_begin(kh2); i != kh_end(kh2); ++i){
                if (kh_exist(kh2, i)) c2_pos[kh_val(kh2, i)] = kh_key(kh2, i);
            }
            for (i = 0; i < size; ++i){
                start = input[i]->start;
                if ((key1 = kh_get(int32, kh1s, start->pos))  != kh_end(kh1s)){
                    for (j = 0; j < input[i]->end->size; ++j){
                        end = input[i]->end->data + j;
                        if ((key2 = kh_get(int32, kh2, end->pos)) != kh_end(kh2)){
                            c2_count[kh_val(kh2, key2)] += end->count;
                        }
                    }
                }
            }
            for (i = 0; i < member_count2; ++i) total_count2+=c2_count[i];
            shrink_rna_group(c2_pos, c2_count, &shrink_index2, &shrink_size2, member_count2, total_count * shrink_prop);
            for (i = shrink_index2; i < shrink_index2 + shrink_size2; ++i){
                key_i = kh_put(int32, kh2s, c2_pos[i], &ret);
                kh_val(kh2s, key_i) = member_count2_s++;
            }

            transcript_cluster* new_tc;
            new_tc = calloc(1, sizeof(*new_tc));
            new_tc->tid = input[0]->tid;
            new_tc->strand = input[0]->strand;
            new_tc->fstart = c1_pos[shrink_index];
            new_tc->fend = c1_pos[shrink_index+shrink_size-1];
            if (new_tc->fend < new_tc->fstart) {
                int32_t v;
                v = new_tc->fstart;
                new_tc->fstart = new_tc->fend;
                new_tc->fend = v;
            }
            new_tc->tstart = c2_pos[shrink_index2];
            new_tc->tend = c2_pos[shrink_index2+shrink_size2-1];
            if (new_tc->tend < new_tc->tstart) {
                int32_t v;
                v = new_tc->tstart;
                new_tc->tstart = new_tc->tend;
                new_tc->tend = v;
            }
            vec_add(tc, tc, new_tc);

            free(c2_pos);
            free(c2_count);
            free(c1_pos);
            free(c1_count);
        }
        free(pos2);
        free(count2);
        free(cluster2);
    }

    kh_destroy(db, kh);
    kh_destroy(int32, kh1);
    kh_destroy(int32, kh2);
    kh_destroy(int32, kh1s);
    kh_destroy(int32, kh2s);
    free(cluster);
    free(pos);
    free(count);
}

