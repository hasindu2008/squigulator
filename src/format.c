/* @file  format.c
** squigulator, a nanopore signal simulator
**
** @@
******************************************************************************/
/*
MIT License

Copyright (c) 2023 Hasindu Gamaarachchi

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "format.h"
#include "error.h"
#include "str.h"
#include "seq.h"
#include "version.h"

aln_t *init_aln() {
    aln_t *aln = (aln_t*)malloc(sizeof(aln_t));
    MALLOC_CHK(aln);
    aln->ss_n = 0;
    aln->ss_c = 1000; //todo can be the readlen
    aln->ss = (int32_t*)malloc(aln->ss_c * sizeof(int32_t));
    MALLOC_CHK(aln->ss);

    aln->prefix_end = 0;
    return aln;
}

void free_aln(aln_t *aln) {
    free(aln->ss);
    free(aln);
}

char *paf_str(aln_t *aln) {
    kstring_t str;
    kstring_t *sp = &str;
    str_init(sp, aln->tlen*3+1000);

    assert(aln->sig_end > aln->sig_start);
    sprintf_append(sp, "%s\t%ld\t%ld\t%ld\t%c\t", aln->read_id, (long)aln->len_raw_signal,(long)aln->sig_start, (long)aln->sig_end, aln->strand);
    //target: name, start, end
    sprintf_append(sp, "%s\t%ld\t%ld\t%ld\t", aln->tid, aln->tlen, (long)aln->t_st, (long)aln->t_end);
    //residue matches, block len, mapq
    int64_t blocklen = aln->t_end - aln->t_st;
    blocklen = blocklen > 0 ? blocklen : -blocklen;
    sprintf_append(sp, "%ld\t%ld\t%d\t",blocklen,blocklen,255);
    sprintf_append(sp, "sc:f:%f\t",1.0);
    sprintf_append(sp, "sh:f:%f\t",0.0);
    sprintf_append(sp, "ss:Z:");

    int8_t rna = aln->t_st > aln->t_end ? 1 : 0;
    for(int i=0; i<aln->ss_n; i++) {
        int idx = rna ? aln->ss_n - i - 1 : i;
        sprintf_append(sp, "%d,", aln->ss[idx]);
    }

    sprintf_append(sp, "\n");
    str.s[str.l] = '\0';
    return sp->s;
}


void sam_hdr_wr(FILE *fp, ref_t *ref) {
    for(int i=0; i<ref->num_ref; i++) {
        fprintf(fp, "@SQ\tSN:%s\tLN:%ld\n", ref->ref_names[i], (long)ref->ref_lengths[i]);
    }
    fprintf(fp, "@PG\tID:squigulator\tPN:squigulator\tVN:%s\n", SQ_VERSION);

}

char *sam_str(aln_t *aln, char *seq, char *rname, int32_t ref_pos_st) {
    kstring_t str;
    kstring_t *sp = &str;
    str_init(sp, aln->tlen*3+1000);

    assert(aln->sig_end > aln->sig_start);

    int flag = aln->strand == '+' ? 0 : 16;
    sprintf_append(sp, "%s\t%d\t", aln->read_id, flag); //qname, flag
    sprintf_append(sp, "%s\t%ld\t%d\t", rname, (long)ref_pos_st+1, 255); //rname, pos, mapq
    sprintf_append(sp, "%ldM\t%c\t%d\t%d\t",strlen(seq), '*', 0, 0); //cigar, rnext, pnext, tlen
    if (aln->strand == '+'){
        str_cat(sp, seq, strlen(seq));//seq
    } else {
        char *rc = reverse_complement(seq);
        str_cat(sp, rc, strlen(rc));
        free(rc);
    }
    sprintf_append(sp, "\t%c\t",'*'); //seq, qual

    sprintf_append(sp, "si:Z:%ld,%ld,%ld,%ld\t",(long)aln->sig_start, (long)aln->sig_end, (long)aln->si_st_ref, (long)aln->si_end_ref);
    sprintf_append(sp, "ss:Z:");

    int8_t rna = aln->t_st > aln->t_end ? 1 : 0;
    for(int i=0; i<aln->ss_n; i++) {
        int idx = rna ? aln->ss_n - i - 1 : i;
        sprintf_append(sp, "%d,", aln->ss[idx]);
    }

    sprintf_append(sp, "\n");
    str.s[str.l] = '\0';
    return sp->s;
}

