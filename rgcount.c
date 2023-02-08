#include <stdint.h>
#include <getopt.h>
#include "htslib/khash.h"

#include "vector.h"


#include "file.h"
#include "common.h"

struct rgcount_option{
    char **input;
    int input_count;
    char *rgfile;
    char *output;
};

#define RGCOUNT_VERSION "1.0.0"
void rgcount_usage(const char *msg){
    const char *usage_info = "\
rgcount: evaluate the abundance of transcript clusters given alignments.\n\
Usage:  rgcount [options] --fi <alignment file> ... --rg <rg file> --fo <output file>\n\
[options]\n\
-i/--fi             : input bam file (should be unsorted). [required]\n\
-r/--rg             : transcript clusters to be evaluated. [required]\n\
-o/--fo             : output rg file. [required]\n\
-h/--help           : show help informations.\n\n";
    if (msg==NULL || msg[0] == '\0') fprintf(stderr, "%s", usage_info);
    else fprintf(stderr, "%s\n\n%s", msg, usage_info);
    exit(1);
}

void rgcount_version(){
    fprintf(stderr, "rgcount-%s\n\n", RGCOUNT_VERSION);
}
void rgcount_option(struct rgcount_option *options, int argc, char *argv[]){
    char c;
    options->input = NULL;
    options->input_count = 0;
    options->rgfile = NULL;
    options->output = NULL;
    if (argc == 1) rgcount_usage("");
    const char *short_options = "hvi:o:r:";
    const struct option long_options[] =
            {
                    { "help" , no_argument , NULL, 'h' },
                    { "version" , no_argument , NULL, 'v' },
                    { "fi" , required_argument , NULL, 'i' },
                    { "rg" , required_argument , NULL, 'r' },
                    { "fo" , required_argument, NULL, 'o' },
                    { NULL, 0, NULL, 0} ,
            };

    while ((c = getopt_long(argc, argv, short_options, long_options, NULL)) >= 0)
    {
        switch (c)
        {
            case 'h':
                rgcount_usage("");
                exit(0);
            case 'v':
                rgcount_version();
                exit(0);
            case 'o':
                options->output = optarg;
                break;
            case 'r':
                options->rgfile = optarg;
                break;
            case 'i':
                options->input = &argv[optind-1];
                options->input_count = 1;
                while (optind < argc && argv[optind][0] != '-')  {
                    options->input_count++;
                    optind++;
                }
                break;
            default:
                rgcount_usage("[rgcount] Error:unrecognized parameter");
        }
    }
    if (argc != optind) rgcount_usage("[rgcount] Error:unrecognized parameter");
};

int main(int argc, char **argv) {
    struct rgcount_option options;
    rgcount_option(&options, argc, argv);

    bam_parser_t *bp = NULL;
    int i = 0;

    bp = bam_parser_open(options.input[0]);
    khash_t (str2int32) *chrom2id = kh_init(str2int32);
    for (i = 0; i < bp->hdr->n_targets; ++i){
        int ret;
        khiter_t k = kh_put(str2int32, chrom2id, bp->hdr->target_name[i], &ret);
        kh_val(chrom2id, k) = i;
    }
    vec_t (rg) *rg = vec_init(rg);
    rg_read(options.rgfile, rg, chrom2id);

    khash_t (record) *rh = kh_init(record);
    vec_t (record) *rv = vec_init(record);
    double n_valid_read_pair = 0;
    for (int fidx = 0; fidx < options.input_count; ++fidx) {
        bp = bam_parser_open(options.input[fidx]);
        while (!bam_parser_read(bp)) {
            bam_parser_process(bp, rh, rv);
        }
        n_valid_read_pair+=bp->n_valid_read_pair;
    }
    rg_count(rg, rv, bp->hdr->n_targets, bp->hdr->target_len);

    char comment[100];
    sprintf(comment, "%f", n_valid_read_pair);
    rg_output(options.output, rg, bp->hdr->target_name, comment);

    bam_parser_close(bp);
    kh_destroy(record, rh);
    vec_destroy(record, rv);
    rg_free(rg);

    return 0;
}
