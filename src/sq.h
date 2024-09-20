/* @file sq.h
**
******************************************************************************/

#ifndef SQ_H
#define SQ_H

#include <stdint.h>
#include "ref.h"
#include "rand.h"
#include "slow5/slow5.h"

//model types
#define MODEL_TYPE_NUCLEOTIDE 1
#define MODEL_TYPE_METH 2

#define MAX_KMER_SIZE 9 //maximum k-mer size
#define MAX_NUM_KMER 262144   //maximum number of k-mers in nucleotide model
#define MAX_NUM_KMER_METH 1953125 //maximum number of k-mers in methylated model

//default model IDs
#define MODEL_ID_DNA_R9_NUCLEOTIDE 1
#define MODEL_ID_RNA_R9_NUCLEOTIDE 2
#define MODEL_ID_DNA_R10_NUCLEOTIDE 3
#define MODEL_ID_RNA_RNA004_NUCLEOTIDE 4
#define MODEL_ID_DNA_R9_CPG 5
#define MODEL_ID_DNA_R10_CPG 6

/*******************************************************
 * flags related to the user specified options (opt_t) *
 *******************************************************/
#define SQ_RNA 0x001 //if RNA or not
#define SQ_FULL_CONTIG 0x002 //if fullcontigs
#define SQ_IDEAL 0x004 //if set, ideal signals with no noise
#define SQ_IDEAL_TIME 0x008 //signal with no time domain noise
#define SQ_IDEAL_AMP 0x010 //signal with no time amplitude domain noise
#define SQ_PREFIX 0x020 //generate prefix or not
#define SQ_R10 0x040 //R10 or R9
#define SQ_PAF_REF 0x080 //in paf output, use ref as target
#define SQ_TRANS_TRUNC 0x100 //trans-trunc
#define SQ_CDNA 0x200 //CDNA
#define SQ_ONT 0x400 //ont friendly

#define WORK_STEAL 1 //simple work stealing enabled or not (no work stealing mean no load balancing)
#define STEAL_THRESH 1 //stealing threshold

typedef struct {
    double digitisation;
    double sample_rate;
    double bps;
    double range;
    double offset_mean;
    double offset_std;
    double median_before_mean;
    double median_before_std;
    double dwell_mean;
    double dwell_std;
} profile_t;

/* k-mer model */
typedef struct {
    float level_mean;
    float level_stdv;

#ifdef CACHED_LOG
    float level_log_stdv;     //pre-calculated for efficiency
#endif
} model_t;

typedef struct{
    int32_t rlen;
    int64_t seed;

    //int8_t rna;
    const char *model_file;
    //int8_t prefix;

    uint32_t flag;

    int32_t num_thread; //t
    int32_t batch_size; //K

    float amp_noise;

    const char *meth_freq;
    const char *meth_model_file;
} opt_t;

typedef struct {
    int64_t *ref_pos;
    int64_t *rand_strand;
    nrng_t **rand_time;
    grng_t **rand_rlen;
    nrng_t **rand_offset;
    nrng_t **rand_median_before;
    int64_t *rand_meth;

    profile_t profile;
    model_t *model;
    nrng_t ***kmer_gen;
    uint32_t kmer_size;
    uint32_t num_kmer;

    model_t *cpgmodel;
    //opt
    opt_t opt;

    //files
    FILE *fp_fasta;
    FILE *fp_paf;
    FILE *fp_sam;
    slow5_file_t *sp;

    //reference
    ref_t *ref;

    double process_db_time;
    double output_time;

    //stats //set by output_db
    int64_t n_samples;
    int64_t total_reads;
} core_t;


/* a batch of read data (dynamic data based on the reads) */
typedef struct {
    int32_t n_rec;
    int32_t capacity_rec;

    char **mem_records;
    size_t *mem_bytes;

    char **fasta;
    char **paf;
    char **sam;

    int64_t n_samples;
} db_t;


/* argument wrapper for the multithreaded framework used for data processing */
typedef struct {
    core_t* core;
    db_t* db;
    int32_t starti;
    int32_t endi;
    void (*func)(core_t*,db_t*,int,int);
    int32_t thread_index;
#ifdef WORK_STEAL
    void *all_pthread_args;
#endif
} pthread_arg_t;

#endif
