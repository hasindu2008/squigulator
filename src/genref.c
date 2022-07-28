
/* @file genref.c
**

** @@
******************************************************************************/

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <zlib.h>
#include <math.h>

#include "error.h"
#include "sigsim.h"
#include "ref.h"

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

//#define NORMALISE_ALL 1

static inline void normalise(float *rawptr, uint64_t n){

    uint64_t start_idx =  0;
    uint64_t end_idx = n;

    float event_mean = 0;
    float event_var = 0;
    float event_stdv = 0;
    float num_samples = end_idx-start_idx;

    for(uint64_t j=start_idx; j<end_idx; j++){
        event_mean += rawptr[j];
    }
    event_mean /= num_samples;
    for(uint64_t j=start_idx; j<end_idx; j++){
        event_var += (rawptr[j]-event_mean)*(rawptr[j]-event_mean);
    }
    event_var /= num_samples;
    event_stdv = sqrt(event_var);

    for(uint64_t j=start_idx; j<end_idx; j++){
        rawptr[j] = (rawptr[j]-event_mean)/event_stdv;
    }

}

static inline void normalise_ref(refsynth_t *ref, int8_t rna){

    uint64_t n = 0;
    float mean = 0;
    float var = 0;
    float stdv = 0;
    for(int i=0; i<ref->num_ref; i++){
        n += ref->ref_lengths[i];
        for(int32_t j= 0; j<ref->ref_lengths[i]; j++){
            mean += ref->forward[i][j];
        }
    }
    mean /= n;

    for(int i=0; i<ref->num_ref; i++){
        for(int32_t j= 0; j<ref->ref_lengths[i]; j++){
            float a = ref->forward[i][j]-mean;
            var += a*a;
        }
    }
    var /= n;
    stdv =  sqrt(var);

    for(int i=0; i<ref->num_ref; i++){
        for(int32_t j= 0; j<ref->ref_lengths[i]; j++){
            ref->forward[i][j] = (ref->forward[i][j]-mean)/stdv;
            if (!rna){
                ref->reverse[i][j] = (ref->reverse[i][j]-mean)/stdv;
            }
        }

    }


}

refsynth_t *gen_ref(const char *genome, model_t *pore_model, uint32_t kmer_size, uint32_t flag, int32_t query_size){

    gzFile fp;
    kseq_t *seq;
    int l;
    fp = gzopen(genome, "r");
    F_CHK(fp,genome);
    seq = kseq_init(fp);
    MALLOC_CHK(seq);
    int8_t rna = flag & SIGSIM_RNA;

    refsynth_t *ref = (refsynth_t *) malloc(sizeof(refsynth_t));

    int c = 1;
    ref->ref_seq_lengths = (int32_t *) malloc(sizeof(int32_t));
    ref->ref_lengths = (int32_t *) malloc(sizeof(int32_t));
    ref->ref_names = (char **) malloc(sizeof(char *));
    ref->forward = (float **) malloc(sizeof(float *));
    if(!rna){
        ref->reverse = (float **) malloc(sizeof(float *));
    }
    else{
        ref->reverse = NULL;
    }

    int i = 0;
    while ((l = kseq_read(seq)) >= 0) {
        assert(l==(int)strlen(seq->seq.s));

        if(i+1 > c){
            c *= 2;
            ref->ref_lengths = (int32_t *) realloc(ref->ref_lengths, c*sizeof(int32_t));
            ref->ref_seq_lengths = (int32_t *) realloc(ref->ref_seq_lengths, c*sizeof(int32_t));
            ref->ref_names = (char **) realloc(ref->ref_names, c*sizeof(char *));
            ref->forward = (float **) realloc(ref->forward, c*sizeof(float *));
            if(!rna){
                ref->reverse = (float **) realloc(ref->reverse, c*sizeof(float *));
            }
        }

        int32_t ref_len;
        if(!rna || flag & SIGSIM_REF){ //dna or use full reference
            ref_len = l+1-kmer_size;
        }
        else{ //rna
            uint32_t rlen_heu=query_size*1.5;
            ref_len = (rlen_heu > l+1-kmer_size ? l+1-kmer_size : rlen_heu);
            LOG_TRACE("Only %d bases of %d bases in reference sequence will be used\n", ref_len, l);
        }
        //int32_t ref_len = rna ? ((uint32_t)query_size > l+1-kmer_size ? l+1-kmer_size : query_size) : l+1-kmer_size;
        //fprintf(stderr,"%s\t%d\n",seq->name.s,ref_len);
        //int32_t ref_len =  l+1-kmer_size;
        ref->ref_lengths[i] = ref_len;
        ref->ref_seq_lengths[i] = l;
        ref->ref_names[i] = (char *) malloc(strlen(seq->name.s)+1);
        strcpy(ref->ref_names[i], seq->name.s);
        ref->forward[i] = (float *) malloc(ref_len*sizeof(float));


        char *rc = NULL;
        if(!rna){
            ref->reverse[i] = (float *) malloc(ref_len*sizeof(float));
            rc = reverse_complement(seq->seq.s);
        }

        // fprintf(stderr,"%s\n",seq->seq.s);
        // fprintf(stderr,"%s\n",rc);

        if(!rna){ //dna (if DNA we just use the full reference)
            for (int j=0; j< ref_len; j++){
                uint32_t kmer_rank = get_kmer_rank(seq->seq.s+j, kmer_size);
                ref->forward[i][j] = pore_model[kmer_rank].level_mean;

                kmer_rank = get_kmer_rank(rc+j, kmer_size);
                ref->reverse[i][j] = pore_model[kmer_rank].level_mean;
            }
        }else{ //rna
            if(flag & SIGSIM_INV){ //would be incorrect - have not tested recently
                VERBOSE("%s","Reversing the reference to be 5' -> 3'\n");
                // char *f = seq->seq.s;
                // char *reverse = (char *) malloc(strlen(f)*sizeof(char));
                // for(int j=0; j<(int)strlen(f); j++){
                //     reverse[j] = f[strlen(f)-j-1];
                // }
                char *seq_end = seq->seq.s+l-ref_len-(kmer_size-1);
                for(int j=0; j<ref_len; j++){
                    uint32_t kmer_rank = get_kmer_rank(seq_end+j, kmer_size);
                    ref->forward[i][ref_len-j-1] = pore_model[kmer_rank].level_mean;
                }
                // for(int j=0; j<ref_len; j++){
                //     uint32_t kmer_rank = get_kmer_rank(seq->seq.s+j, kmer_size);
                //     ref->forward[i][ref_len-j-1] = pore_model[kmer_rank].level_mean;
                // }
            }
                // free(reverse);
            else{
                char *st;
                if(flag &  SIGSIM_END){ //if the end of query then it is the beginning of the reference in RNA
                    st = seq->seq.s;

                }
                else{ //if the beginning of query then it is the end of the reference in RNA
                    st = seq->seq.s+l-ref_len-(kmer_size-1);
                }
                LOG_TRACE("%s:%ld-%ld\n",seq->name.s,st-seq->seq.s,st-seq->seq.s+ref_len);
                for (int j=0; j< ref_len; j++){
                    uint32_t kmer_rank = get_kmer_rank(st+j, kmer_size);
                    ref->forward[i][j] = pore_model[kmer_rank].level_mean;
                }
            }
        }

        // for (int j=0; j< ref->ref_lengths[i]; j++){
        //     fprintf(stderr,"%f,",ref->forward[i][j]);
        // }
        // fprintf(stderr,"\n");
        // for (int j=0; j< ref->ref_lengths[i]; j++){
        //     fprintf(stderr,"%f,",ref->reverse[i][j]);
        // }
        // fprintf(stderr,"\n");

        #ifndef NORMALISE_ALL
        normalise(ref->forward[i], ref->ref_lengths[i]);
        #endif
            if (!rna){
            free(rc);
            #ifndef NORMALISE_ALL
            normalise(ref->reverse[i], ref->ref_lengths[i]);
            #endif
        }


        // for (int j=0; j< ref->ref_lengths[i]; j++){
        //     fprintf(stderr,"%f,",ref->forward[i][j]);
        // }
        // fprintf(stderr,"\n");

        i++;

    }

    ref->num_ref = i;
    #ifdef NORMALISE_ALL
    normalise_ref(ref,rna);
    #endif


    kseq_destroy(seq);
    gzclose(fp);

    return ref;

}

void free_ref(refsynth_t *ref){

    for(int i=0;i<ref->num_ref;i++){
        free(ref->ref_names[i]);
        free(ref->forward[i]);
        if(ref->reverse){
            free(ref->reverse[i]);
        }
    }

    free(ref->ref_lengths);
    free(ref->ref_seq_lengths);
    free(ref->ref_names);
    free(ref->forward);
    free(ref->reverse);

    free(ref);
}
