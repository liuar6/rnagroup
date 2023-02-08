#include "htslib/bgzf.h"
#include "htslib/sam.h"

typedef struct bam_parser_t{
    BGZF *fp;
    sam_hdr_t *hdr;
    bam1_t *alignment[201];
    bam1_t* read1[100];
    bam1_t* read2[100];
    int read_valid[100];
    int read_index;
    int valid_count;
    int ret;
    int n_valid_read_pair;

} bam_parser_t;

bam_parser_t *bam_parser_open(const char* input);
int bam_parser_read(bam_parser_t *p);
void bam_parser_close(bam_parser_t *p);