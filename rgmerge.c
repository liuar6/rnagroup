#include <getopt.h>
#include "file.h"
#include "common.h"

#define min(a, b) (((a)>(b))?(b):(a))
#define max(a, b) (((a)>(b))?(a):(b))
#define RANGE_INTERSECT(start1, end1, start2, end2) (min(end1, end2)-max(start1, start2) > 0)

struct rgmerge_option{
    char **input;
    int input_count;
    char *ref;
    char *output;
};

#define RGMERGE_VERSION "1.0.0"
void rgmerge_usage(const char *msg){
    const char *usage_info = "\
rgmerge: merge the transcript clusters from runs of rgfind.\n\
Usage:  rgmerge [options] --rg <rg file> ...\n\
[options]\n\
-r/--rg             : rg files to be merged. [required]\n\
-b/--bam            : bam file (only provide metadata). [required]\n\
-o/--fo             : output rg file. [required]\n\
-h/--help           : show help informations.\n\n";
    if (msg==NULL || msg[0] == '\0') fprintf(stderr, "%s", usage_info);
    else fprintf(stderr, "%s\n\n%s", msg, usage_info);
    exit(1);
}

void rgmerge_version(){
    fprintf(stderr, "rgmerge-%s\n\n", RGMERGE_VERSION);
}
void rgmerge_option(struct rgmerge_option *options, int argc, char *argv[]){
    char c;
    options->input = NULL;
    options->input_count = 0;
    if (argc == 1) rgmerge_usage("");
    const char *short_options = "hvb:r:o:";
    const struct option long_options[] =
            {
                    { "help" , no_argument , NULL, 'h' },
                    { "version" , no_argument , NULL, 'v' },
                    { "rg" , required_argument , NULL, 'r' },
                    { "bam" , required_argument , NULL, 'b' },
                    { "fo" , required_argument , NULL, 'o' },
                    { NULL, 0, NULL, 0} ,
            };

    while ((c = getopt_long(argc, argv, short_options, long_options, NULL)) >= 0)
    {
        switch (c)
        {
            case 'h':
                rgmerge_usage("");
                exit(0);
            case 'v':
                rgmerge_version();
                exit(0);
            case 'o':
                options->output = optarg;
                break;
            case 'b':
                options->ref = optarg;
                break;
            case 'r':
                options->input = &argv[optind-1];
                options->input_count = 1;
                while (optind < argc && argv[optind][0] != '-')  {
                    options->input_count++;
                    optind++;
                }
                break;
            default:
                rgmerge_usage("[rgmerge] Error:unrecognized parameter");
        }
    }
    if (argc != optind) rgmerge_usage("[rgmerge] Error:unrecognized parameter");
};
int main(int argc, char **argv) {
    struct rgmerge_option options;
    rgmerge_option(&options, argc, argv);
    const char *output = options.output;

    bam_parser_t *bp;
    int i;
    bp = bam_parser_open(options.ref);
    khash_t (str2int32) *chrom2id = kh_init(str2int32);
    for (i = 0; i < bp->hdr->n_targets; ++i){
        int ret;
        khiter_t k = kh_put(str2int32, chrom2id, bp->hdr->target_name[i], &ret);
        kh_val(chrom2id, k) = i;
    }
    vec_t (rg) *rg = vec_init(rg);
    vec_t (rg) *rg1 = vec_init(rg);
    for (int fidx = 0; fidx < options.input_count; ++fidx) {
        rg_read(options.input[fidx], rg, chrom2id);
    }
    bioidx_t *bidx1 = bioidx_init();
    bioidx_itr_t *bitr = bioidx_itr_init();
    int last_x_count = INT_MAX;
    while (1) {
        int x_count = 0;
        for (i = 0; i < rg->size; ++i) {
            rg_t *t = rg->data[i];
            if (t->strand == 'x') {
                x_count++;
                continue;
            }
            rg_t *t1;
            bioidx_search(bidx1, bitr, bioidx_key(t->tid, t->strand), t->fstart, t->fend + 1);
            while ((t1 = bioidx_itr_next(bitr))) {
                if (RANGE_INTERSECT(t->tstart, t->tend + 1, t1->tstart, t1->tend + 1)) {
                    bioidx_itr_remove(bitr);
                    if (t == t1) continue;
                    t->fstart = min(t->fstart, t1->fstart);
                    t->fend = max(t->fend, t1->fend);
                    t->tstart = min(t->tstart, t1->tstart);
                    t->tend = max(t->tend, t1->tend);
                    t1->strand ='x';
                }
            }
            bioidx_insert(bidx1, bioidx_key(t->tid, t->strand), t->fstart, t->fend + 1, t);
        }
        fprintf(stderr, "%d\n", x_count);
        if (x_count == last_x_count) break;
        last_x_count = x_count;
    }
    for (i = 0; i < rg->size; ++i) {
        rg_t *t = rg->data[i];
        if (t->strand != 'x')
            vec_add(rg, rg1, t);
    }
    rg_output(output, rg1, bp->hdr->target_name, NULL);
    bam_parser_close(bp);
    bioidx_destroy(bidx1);
    free(bitr);
    rg_free(rg);
}