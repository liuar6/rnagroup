#include "htslib/sam.h"
#include "htslib/bgzf.h"

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

bam_parser_t *bam_parser_open(const char* input){
    bam_parser_t *p;
    p = malloc(sizeof(*p));
    p->fp = bgzf_open(input, "r");
    p->hdr = bam_hdr_read(p->fp);
    for (int i = 0; i < 201; ++i) p->alignment[i]=bam_init1();
    p->n_valid_read_pair = 0;
    p->ret = bam_read1(p->fp, p->alignment[0]);
    return p;
}
void bam_parser_close(bam_parser_t *p){
    int i;
    for (i = 0; i < 201; ++i) bam_destroy1(p->alignment[i]);
    bam_hdr_destroy(p->hdr);
    bgzf_close(p->fp);
    free(p);
}
int bam_parser_read(bam_parser_t *p){
    if (p->ret <= 0) return -1;
    BGZF *fp = p->fp;
    bam1_t **alignment = p->alignment;
    bam1_t **read1 = p->read1;
    bam1_t **read2 = p->read2;
    bam1_t *b;
    int *read_valid = p->read_valid;
    int ret;

    int record_count = 0;
    int record_index = 0;

    p->read_index = 0;
    while ((ret = bam_read1(fp, alignment[++record_count])) > 0 &&
           strcmp(bam_get_qname(alignment[0]), bam_get_qname(alignment[record_count])) == 0);
    while (record_index < record_count) {
        if ((alignment[record_index]->core.flag & BAM_FPAIRED) &&
            (alignment[record_index]->core.flag & BAM_FPROPER_PAIR)) {
            read1[p->read_index] = alignment[record_index];
            read2[p->read_index] = alignment[record_index + 1];
            record_index += 2;
        } else {
            read1[p->read_index] = alignment[record_index];
            read2[p->read_index] = NULL;
            record_index++;
        }
        if (read1[p->read_index]->core.flag & BAM_FREAD2) {
            b = read1[p->read_index];
            read1[p->read_index] = read2[p->read_index];
            read2[p->read_index] = b;
        }
        p->read_index++;
    }
    p->valid_count = 0;
    for (int i = 0; i < p->read_index; ++i) {
        read_valid[i] = 1;
        if (read1[i] == NULL || read2[i] == NULL) {
            read_valid[i] = 0;
            continue;
        }
        if (!(read1[i]->core.flag & BAM_FPAIRED) && !(read1[i]->core.flag & BAM_FPROPER_PAIR)) {
            read_valid[i] = 0;
            continue;
        }
        if (read1[i]->core.isize > 700 || read1[i]->core.isize < -700) {
            read_valid[i] = 0;
            continue;
        }
        if (read_valid[i]) p->valid_count++;
    }
    if (p->valid_count > 0) p->n_valid_read_pair++;

    b = alignment[0];
    alignment[0] = alignment[record_count];
    alignment[record_count] = b;

    p->ret = ret;
    return 0;
}


