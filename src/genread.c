/* @file  genread.c
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

#include <assert.h>
#include "seq.h"
#include "sq.h"
#include "format.h"


const char* polya = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
const char *adaptor_dna = "GGCGTCTGCTTGGGTGTTTAACCTTTTTTTTTTAATGTACTTCGTTCAGTTACGTATTGCT";
const char *adaptor_rna = "TGATGATGAGGGATAGACGATGGTTGTTTCTGTTGGTGCTGATATTGCTTTTTTTTTTTTTATGATGCAAGATACGCAC";

int8_t short_warn = 0;

int16_t * gen_sig_core_seq(core_t *core, int16_t *raw_signal, int64_t* n, int64_t *c, double offset, const char *read, int32_t len, int tid, aln_t *aln);

void gen_prefix_dna(int16_t *raw_signal, int64_t* n, int64_t *c, profile_t *profile, double offset){

    //pre
    float s = 80;
    for(int i=0;i<500;i++){
        raw_signal[*n] = s*(profile->digitisation)/(profile->range)-(offset);
        *n = *n+1;
    }

    //jump
    s = 140;
    for(int i=0;i<500;i++){
        raw_signal[*n] = s*(profile->digitisation)/(profile->range)-(offset);
        *n = *n+1;
    }

    //adaptor
    s=90;
    for(int i=0;i<1000;i++){
        raw_signal[*n] = s*(profile->digitisation)/(profile->range)-(offset);
        *n = *n+1;
    }

}


int16_t *gen_prefix_rna(core_t *core, int16_t *raw_signal, int64_t* n, int64_t *c, double offset, int tid, aln_t *aln){

    profile_t *profile = &core->profile;
    //float s;
    //polya
    // raw_signal=gen_sig_core_seq(core, raw_signal, n, c, offset, polya, strlen(polya));

    //adaptor
    int st = *n-(strlen(adaptor_rna)*(int)profile->dwell_mean);
    // raw_signal=gen_sig_core_seq(core, raw_signal, n, c, offset, adaptor_rna, strlen(adaptor_rna));
    int end = *n;
    //fprintf(stderr, "adaptor_rna: %d %d\n", st, end);
    int16_t off = 30*(profile->digitisation)/(profile->range);
    for(int i=st; i<end; i++){
        raw_signal[i] = raw_signal[i]-off;
    }

    const char *stall = "AAAAAGAAAAAACCCCCCCCCCCCCCCCCC";
    raw_signal=gen_sig_core_seq(core, raw_signal, n, c, offset, stall, strlen(stall), tid, aln);


    return raw_signal;
}

char *attach_prefix(core_t *core, const char *read, int32_t *len, aln_t *aln){

    char *s = (char *)read;
    if(core->opt.flag & SQ_RNA){
        int polya_len = strlen(polya);
        int adapt_len = strlen(adaptor_rna);
        char *seq = malloc((*len+polya_len+adapt_len+1)*sizeof(char));
        MALLOC_CHK(seq);
        strncpy(seq, read, *len);
        strncpy(seq+*len, polya, polya_len);
        strncpy(seq+*len+polya_len, adaptor_rna, adapt_len);
        *len = *len+polya_len+adapt_len;
        seq[*len] = '\0';
        s = seq;
    } else{
        const char *stall = "TTTTTTTTTTTTTTTTTTAATCAA";
        int stall_len = strlen(stall);
        int adapt_len = strlen(adaptor_dna);
        char *seq = malloc((*len+stall_len+adapt_len+1)*sizeof(char));
        MALLOC_CHK(seq);
        strncpy(seq, stall, stall_len);
        strncpy(seq+stall_len, adaptor_dna, adapt_len);
        strncpy(seq+stall_len+adapt_len, read, *len);
        *len = *len+stall_len+adapt_len;
        seq[*len] = '\0';
        s = seq;
    }
    return s;
}

static inline int32_t is_bad_read(char *seq, int32_t len){
    if (len < 200){
        return -1;
    }
    if (len >= UINT32_MAX){
        return -2;
    }
    int64_t r = 100;
    int nc = 0;
    for (int i=0; i<len; i++){
        if (seq[i] == 'N'){
            nc++;
            int n = round(rng(&r)*3);
            seq[i] = n==0 ? 'A' : n==1 ? 'C' : n==2 ? 'G' : 'T';
            if(nc > 0.1*len){
                return nc;
            }
        }
    }

    return 0;
}

static char *gen_read_dna(core_t *core, ref_t *ref, char **ref_id, int32_t *ref_len, int32_t *ref_pos, int32_t *rlen, char *c, int tid){

    char *seq = NULL;
    while(1){

        int len = grng(core->rand_rlen[tid]);
        int64_t ref_sum_pos = round(rng(&core->ref_pos[tid])*ref->sum); //serialised pos
        assert(ref_sum_pos <= ref->sum);
        int64_t s = 0;
        int seq_i = 0;
        for(seq_i=0; seq_i<ref->num_ref; seq_i++){ //check manually if logic is right
            s+=ref->ref_lengths[seq_i];
            if(s>=ref_sum_pos){
                *ref_id = ref->ref_names[seq_i];
                *ref_pos = ref_sum_pos-s+ref->ref_lengths[seq_i];
                *ref_len = ref->ref_lengths[seq_i];
                break;
            }
        }
        assert(s<=ref->sum);

        int64_t strand = round(rng(&core->rand_strand[tid]));
        *c = strand ? '+' : '-';

        seq= (char *)malloc((len+1)*sizeof(char));
        MALLOC_CHK(seq);
        strncpy(seq,(ref->ref_seq[seq_i])+(*ref_pos),len);
        seq[len] = '\0';
        *rlen = strlen(seq);
        LOG_DEBUG("%d\t%d\t%d",*rlen, len, seq_i);
        assert(*rlen <= len);
        int nc = 0;
        if((nc = is_bad_read(seq, *rlen))==0){
           break;
        } else {
            if(nc == -1) {
                LOG_TRACE("Too short read: <%d. %s:%d-%d. Trying again!",200,*ref_id,*ref_pos,*ref_pos+*rlen);
                if(short_warn==0 && ref->ref_lengths[seq_i]<200){
                    WARNING("Reference sequence is too short: %d. Expected to be >=200. Open a pull request if you need support for such tiny references.",ref->ref_lengths[seq_i]);
                    short_warn = 1;
                }
            }else if (nc == -2){
                WARNING("Too long read: >=%d. %s:%d-%d. Trying again!",UINT32_MAX,*ref_id,*ref_pos,*ref_pos+*rlen);
            }
            else{
                LOG_TRACE("Too many Ns in read: %d. %s:%d-%d. Trying again!",nc,*ref_id,*ref_pos,*ref_pos+*rlen);
            }
            free(seq);
        }
    }

    if(*c == '-'){
        char *r = reverse_complement(seq);
        r[*rlen] = '\0';
        free(seq);
        seq = r;
    }

    return seq;
}

static inline int transcript_idx(core_t *core,int tid){
    int seq_i = 0;
    trans_t *trans = core->ref->trans_counts;
    if(trans == NULL){
        seq_i = round(rng(&core->ref_pos[tid])*(core->ref->num_ref-1));
    } else {
        float r = rng(&core->ref_pos[tid]);
        //fprintf(stderr, "r: %f, ", r);
        for(int i=0; i<trans->n; i++){
            if(r <= trans->trans_csum[i]){
                seq_i = trans->trans_idx[i];
                //fprintf(stderr, "seq_i: %d\n", seq_i);
                break;
            }
        }
    }
    return seq_i;
}

static char *gen_read_rna(core_t *core, ref_t *ref, char **ref_id, int32_t *ref_len, int32_t *ref_pos, int32_t *rlen, char *c, int tid){

    char *seq = NULL;
    while(1){

        //int len = grng(core->rand_rlen); //for now the whole transcript is is simulated

        int seq_i = transcript_idx(core,tid); //random transcript
        int len = ref->ref_lengths[seq_i];
        *ref_pos=0;
        *ref_id = ref->ref_names[seq_i];
        *ref_len = len;

        //int64_t strand = round(rng(&core->rand_strand));
        *c = '+'; //strand is alwats plus

        seq= (char *)malloc((len+1)*sizeof(char));
        MALLOC_CHK(seq);
        strncpy(seq,(ref->ref_seq[seq_i])+(*ref_pos),len);
        seq[len] = '\0';
        *rlen = strlen(seq);
        LOG_DEBUG("%d\t%d\t%d",*rlen, len, seq_i);
        assert(*rlen == len);
        int nc = 0;
        if((nc = is_bad_read(seq, *rlen))==0){
           break;
        } else {
            if(nc == -1) {
                LOG_TRACE("Too short read: <%d. %s:%d-%d. Trying again!",200,*ref_id,*ref_pos,*ref_pos+*rlen);
                if(short_warn ==0 && ref->ref_lengths[seq_i]<200){
                    WARNING("Reference sequence is too short: %d. Expected to be >=200. Open a pull request if you need support for such tiny references.",ref->ref_lengths[seq_i]);
                    short_warn = 1;
                }
            } else if (nc == -2){
                WARNING("Too long read: >=%d. %s:%d-%d. Trying again!",UINT32_MAX,*ref_id,*ref_pos,*ref_pos+*rlen);
            }
            else{
                LOG_TRACE("Too many Ns in read: %d. %s:%d-%d. Trying again!",nc,*ref_id,*ref_pos,*ref_pos+*rlen);
            }
            free(seq);
        }

    }

    return seq;

}

char *gen_read(core_t *core, ref_t *ref, char **ref_id, int32_t *ref_len, int32_t *ref_pos, int32_t *rlen, char *c, int8_t rna, int tid){
    return rna ? gen_read_rna(core, ref, ref_id, ref_len, ref_pos, rlen, c, tid): gen_read_dna(core, ref, ref_id, ref_len, ref_pos, rlen, c, tid);
}