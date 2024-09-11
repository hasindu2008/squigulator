/* @file  ref.c
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

#define _XOPEN_SOURCE 700
#include <zlib.h>
#include <math.h>
#include "ref.h"
#include "error.h"
#include "str.h"
#include "kseq.h"
#include "ksort.h"
#include "khash.h"

KSEQ_INIT(gzFile, gzread)

typedef struct {
    char *name;
    int32_t count;
} trans_counts_t;

#define sq_pair_lt(a, b) ((a).count > (b).count)
KSORT_INIT(sq_trans_ct, trans_counts_t, sq_pair_lt)

KHASH_MAP_INIT_STR(sq_s2ui32, uint32_t)

//todo can be optimised for memory by better encoding if necessary
ref_t *load_ref(const char *genome){

    gzFile fp;
    kseq_t *seq;
    int l;
    fp = gzopen(genome, "r");
    F_CHK(fp,genome);
    seq = kseq_init(fp);
    MALLOC_CHK(seq);

    ref_t *ref = (ref_t *) malloc(sizeof(ref_t));
    MALLOC_CHK(ref);

    int c = 1;
    ref->ref_lengths = (int32_t *) malloc(sizeof(int32_t));
    MALLOC_CHK(ref->ref_lengths);
    ref->ref_names = (char **) malloc(sizeof(char *));
    MALLOC_CHK(ref->ref_names);
    ref->ref_seq = (char **) malloc(sizeof(char *));
    MALLOC_CHK(ref->ref_seq);
    ref->sum = 0;

    int8_t large_warn = 1;
    int i = 0;
    while ((l = kseq_read(seq)) >= 0) {
        assert(l==(int)strlen(seq->seq.s));

        if(i+1 > c){
            c *= 2;
            ref->ref_lengths = (int32_t *) realloc(ref->ref_lengths, c*sizeof(int32_t));
            ref->ref_names = (char **) realloc(ref->ref_names, c*sizeof(char *));
            ref->ref_seq = (char **) realloc(ref->ref_seq, c*sizeof(char *));
        }

        ref->ref_lengths[i] = l;
        ref->ref_names[i] = (char *) malloc(strlen(seq->name.s)+1);
        MALLOC_CHK(ref->ref_names[i]);
        strcpy(ref->ref_names[i], seq->name.s);
        ref->ref_seq[i] = (char *) malloc((l+1)*sizeof(char));
        MALLOC_CHK(ref->ref_seq[i]);
        strcpy(ref->ref_seq[i], seq->seq.s);

        ref->sum += l;
        i++;

        if(large_warn && ref->sum > 20000000000){
            WARNING("%s","The input FASTA/FASTQ file seems >20 Gbases. You are seeing this warning because part by part loading is not implemented. If you input file is larger than available RAM, terminate the programme and open an issue on github. The feature will then be implemented as soon as possible.");
            large_warn = 0;
        }

    }

    ref->num_ref = i;
    ref->trans_counts = NULL;
    ref->ref_meth = NULL;

    kseq_destroy(seq);
    gzclose(fp);

    VERBOSE("Loaded %d reference sequences with total length %f Mbases",ref->num_ref,ref->sum/1000000.0);

    return ref;

}

void free_trans_count(trans_t *t){
    free(t->trans_csum);
    free(t->trans_idx);
    free(t);
}

void free_ref_sim(ref_t *ref){

    for(int i=0;i<ref->num_ref;i++){
        free(ref->ref_names[i]);
        free(ref->ref_seq[i]);
        if (ref->ref_meth) free(ref->ref_meth[i]);
    }
    if(ref->trans_counts != NULL){
        free_trans_count(ref->trans_counts);
    }

    if(ref->ref_meth != NULL){
        free(ref->ref_meth);
    }

    free(ref->ref_lengths);
    free(ref->ref_names);
    free(ref->ref_seq);
    free(ref);
}

trans_counts_t *load_sorted_trans_table(const char *trans_count, int32_t *table_n_p, int32_t *table_c_p){

    FILE *fp = fopen(trans_count, "r");
    F_CHK(fp, trans_count);

    //buffers for getline
    char* buffer = (char*)malloc(sizeof(char) * (1000));
    MALLOC_CHK(buffer);
    size_t bufferSize = 1000;
    ssize_t readlinebytes = 0;

    int32_t table_n = 0;
    int32_t table_c = 1000;
    trans_counts_t *table = malloc(table_c * sizeof(trans_counts_t));
    MALLOC_CHK(table);

    while ((readlinebytes = getline(&buffer, &bufferSize, fp)) != -1) {
        if(buffer[0] == '#'){
            continue;
        }
        char *name = strtok(buffer, "\t");
        char *count = strtok(NULL, "\t");
        if(count == NULL){
            ERROR("Invalid format in file %s. Each line should have two tab separated columns. The first column should be the transcript name and the second column should be the transcipt count", trans_count);
            exit(EXIT_FAILURE);
        }
        int32_t c = atoi(count);
        if(count<=0){
            ERROR("Invalid count in file %s. The count should be a positive integer", trans_count);
            exit(EXIT_FAILURE);
        }
        if(table_n == table_c){
            table_c *= 2;
            table = realloc(table, table_c * sizeof(trans_counts_t));
            MALLOC_CHK(table);
        }
        table[table_n].name = malloc(strlen(name)+1);
        MALLOC_CHK(table[table_n].name);
        strcpy(table[table_n].name, name);
        table[table_n].count = c;
        table_n++;

        assert(table_n <= table_c);
        assert(table_n > 0);
        assert(table_c > 0);

    }

    ks_mergesort(sq_trans_ct, table_n, table, 0);

    free(buffer);
    fclose(fp);

    *table_n_p = table_n;
    *table_c_p = table_c;
    return table;

}


void load_trans_count(const char *trans_count, ref_t *ref){

    int32_t table_n = 0;
    int32_t table_c = 0;
    trans_counts_t *table = load_sorted_trans_table(trans_count, &table_n, &table_c);

    khash_t(sq_s2ui32) *h = kh_init(sq_s2ui32);

    for(int i=0; i<ref->num_ref; i++){
        int absent;
        khint_t k = kh_put(sq_s2ui32, h, ref->ref_names[i], &absent);
        if (absent == -1 || absent == 0) {
            ERROR("Transcript '%s' is duplicated in the provided reference", ref->ref_names[i]);
            exit(EXIT_FAILURE);
        }
        kh_value(h, k) = i;
    }

    trans_t *trans = (trans_t *) malloc(sizeof(trans_t));
    MALLOC_CHK(trans);

    trans->trans_csum = (float *) malloc(table_n*sizeof(float));
    MALLOC_CHK(trans->trans_csum);

    trans->trans_idx = (int32_t *) malloc(table_n*sizeof(int32_t));
    MALLOC_CHK(trans->trans_idx);

    trans->n = table_n;

    double count_sum = 0;
    for(int i=0; i<table_n; i++){
        count_sum += table[i].count;
    }
    double cumsum = 0;
    for(int i=0; i<table_n; i++){
        cumsum += table[i].count/count_sum;
        trans->trans_csum[i] = cumsum;

        khint_t k = kh_get(sq_s2ui32, h, table[i].name);
        if (k == kh_end(h)) {
            ERROR("Transcript '%s' in the transcript count file is not found in the reference", table[i].name);
            exit(EXIT_FAILURE);
        }
        int j = trans->trans_idx[i] = kh_value(h, k);
        assert(j < ref->num_ref);
        if(strcmp(ref->ref_names[j], table[i].name) != 0){
            ERROR("ASSERTION failed for %s", table[i].name);
            exit(EXIT_FAILURE);
        }
    }

    // for(int i=0; i<table_n; i++){
    //     fprintf(stderr, "%f %d\n", trans_csum[i], trans_idx[i]);
    // }

    for(int i=0; i<table_n; i++){
        //fprintf(stderr, "%s\t%d\n", table[i].name, table[i].count);
        free(table[i].name);
    }
    free(table);

    kh_destroy(sq_s2ui32, h);

    ref->trans_counts = trans;

    return;

}

static inline int get_ref_idx(const char *name, ref_t *ref){
    for(int i=0; i< ref->num_ref; i++){
        if(strcmp(ref->ref_names[i],name)==0){
            return i;
        }
    }
    return -1;
}

static inline void check_meth_line(void *ptr, const char *file){
    if(ptr == NULL){
        ERROR("Invalid format in file %s. Each line should have three tab separated values: col1: chromosome, col2: 0-based pos; col3: methylation frequency", file);
        exit(EXIT_FAILURE);
    }
}

void load_meth_freq(const char *meth_freq, ref_t *ref){

    ref->ref_meth = calloc(ref->num_ref, sizeof(uint8_t *));
    MALLOC_CHK(ref->ref_meth);

    FILE *fp = fopen(meth_freq, "r");
    F_CHK(fp, meth_freq);

    //buffers for getline
    char* buffer = (char*)malloc(sizeof(char) * (1000));
    MALLOC_CHK(buffer);
    size_t bufferSize = 1000;
    ssize_t readlinebytes = 0;

    int line = 0;
    while ((readlinebytes = getline(&buffer, &bufferSize, fp)) != -1) {
        if(buffer[0] == '#'){
            line++;
            continue;
        }
        //fprintf(stderr,"%s\n",buffer);


        char *name = strtok(buffer, "\t");
        check_meth_line(name, meth_freq);
        int ref_idx = get_ref_idx(name, ref);
        //fprintf(stderr,"Refidx %d\n",ref_idx);
        if(ref_idx < 0){
            ERROR("There was no such chromosome in the reference. Check line %d value %s of the input methy-freq file.", line, name);
            exit(EXIT_FAILURE);
        }
        char *position = strtok(NULL, "\t");
        check_meth_line(position, meth_freq);
        int32_t pos = atoi(position);
        //fprintf(stderr,"%s,%d\n",position,pos);
        if(pos < 0){
            ERROR("Chromosome position cannot be negative. Check line %d value %d of the input methy-freq file.",line,pos);
            exit(EXIT_FAILURE);
        } else if (pos >= ref->ref_lengths[ref_idx]){
            ERROR("Chromosome %s position must be less than the length %d. Check line %d value %d of the input methy-freq file.", ref->ref_names[ref_idx], ref->ref_lengths[ref_idx], line, pos);
            exit(EXIT_FAILURE);
        }
        char c = ref->ref_seq[ref_idx][pos];
        if(!(c=='C' || c=='c')){
            ERROR("The chromosome %s position %d in the reference was a %c. How can it be methylated C? Check line %d of the input methy-freq file.",name,pos,c,line);
            exit(EXIT_FAILURE);
        }

        char *frequency = strtok(NULL, "\t");
        check_meth_line(frequency, meth_freq);
        float freq = atof(frequency);
        if(freq<0 || freq >1){
            ERROR("Methylation frequency must be between 0 to 1. Check line %d value %f of the input methy-freq file.",line, freq);
            exit(EXIT_FAILURE);
        }
        uint8_t f = roundf(freq * 255);
        if(ref->ref_meth[ref_idx]==NULL){
            ref->ref_meth[ref_idx] = calloc(ref->ref_lengths[ref_idx], sizeof(uint8_t));
            MALLOC_CHK(ref->ref_meth[ref_idx]);
        }
        ref->ref_meth[ref_idx][pos] = f;

        line ++;
        //fprintf(stderr,"%s,%d,%f,%d\n",name,pos,freq,f);

    }

    free(buffer);
    fclose(fp);
    return;
}