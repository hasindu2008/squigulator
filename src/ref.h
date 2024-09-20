#ifndef SQ_REF_H
#define SQ_REF_H

#include <stdint.h>

typedef struct {
    float *trans_csum;
    int32_t *trans_idx;
    int32_t n;
} trans_t;

typedef struct{
    char **ref_names;
    int32_t *ref_lengths;
    char **ref_seq;
    int num_ref;
    int64_t sum;
    trans_t *trans_counts;
    uint8_t **ref_meth;
} ref_t;

ref_t *load_ref(const char *genome);
void free_ref_sim(ref_t *ref);
void load_trans_count(const char *trans_count, ref_t *ref);
void load_meth_freq(const char *meth_freq, ref_t *ref);

#endif