/* @file sigsim.h
**
******************************************************************************/

#ifndef SIGSIM_H
#define SIGSIM_H

#include <stdint.h>
#include "slow5/slow5.h"

#define SIGSIM_VERSION "0.1.0"

//model types
#define MODEL_TYPE_NUCLEOTIDE 1
#define MODEL_TYPE_METH 2

#define MAX_KMER_SIZE 6 //maximum k-mer size
#define MAX_NUM_KMER 4096   //maximum number of k-mers in nucleotide model
#define MAX_NUM_KMER_METH 15625 //maximum number of k-mers in methylated model

//default model IDs
#define MODEL_ID_DNA_NUCLEOTIDE 1
#define MODEL_ID_RNA_NUCLEOTIDE 2

/*******************************************************
 * flags related to the user specified options (opt_t) *
 *******************************************************/

#define SIGSIM_RNA 0x001 //if RNA or not
#define SIGSIM_DTW 0x002 //if dtw-std
#define SIGSIM_INV 0x004 //if set, reverse reference events instead of query events
#define SIGSIM_SEC 0x008 //if secondaries are printed
#define SIGSIM_REF 0x010 //map to the whole reference
#define SIGSIM_END 0x020 //map the end of the query

typedef struct {
    int32_t num_ref;
    char **ref_names;
    int32_t *ref_lengths;
    int32_t *ref_seq_lengths;

    float **forward;
    float **reverse;
} refsynth_t;

/* k-mer model */
typedef struct {
    float level_mean;
    float level_stdv;

#ifdef CACHED_LOG
    float level_log_stdv;     //pre-calculated for efficiency
#endif

#ifdef LOAD_SD_MEANSSTDV
    //float sd_mean;
    //float sd_stdv;
    //float weight;
#endif
} model_t;



#endif