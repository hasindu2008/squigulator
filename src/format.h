#ifndef SQ_FORMAT_H
#define SQ_FORMAT_H

#include <stdint.h>
#include <stdio.h>
#include "ref.h"

typedef struct {
    char *read_id;
    int64_t len_raw_signal;
    int64_t sig_start;
    int64_t sig_end;
    char strand;
    char *tid;
    int64_t tlen;
    int64_t t_st;
    int64_t t_end;

    int64_t si_st_ref; //for sam
    int64_t si_end_ref;
    int32_t *ss;
    int64_t ss_n;
    int64_t ss_c;

    int64_t prefix_end;
} aln_t;

aln_t *init_aln();
void free_aln(aln_t *aln);
char *paf_str(aln_t *aln);
void sam_hdr_wr(FILE *fp, ref_t *ref);
char *sam_str(aln_t *aln, char *seq, char *rname, int32_t ref_pos_st);

#endif