#ifndef SQ_REF_H
#define SQ_REF_H

#include <stdint.h>

typedef struct{
    char **ref_names;
    int32_t *ref_lengths;
    char **ref_seq;
    int num_ref;
    int64_t sum;
} ref_t;

ref_t *load_ref(const char *genome);
void free_ref_sim(ref_t *ref);

#endif