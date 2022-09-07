/* @file misc.h
**
** miscellaneous definitions and function prototypes
** @@
******************************************************************************/

#ifndef MISC_H
#define MISC_H

#include <sys/resource.h>
#include <sys/time.h>
#include <stdint.h>
#include <math.h>
#include "error.h"

/* models */

uint32_t read_model(model_t* model, const char* file, uint32_t type);
uint32_t set_model(model_t* model, uint32_t model_id);

//refsynth_t *gen_ref(const char *genome, model_t *pore_model, uint32_t kmer_size, uint32_t rna, int32_t query_size);
void free_ref(refsynth_t *ref);

// taken from minimap2/misc
static inline double realtime(void) {
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return tp.tv_sec + tp.tv_usec * 1e-6;
}

// taken from minimap2/misc
static inline double cputime(void) {
    struct rusage r;
    getrusage(RUSAGE_SELF, &r);
    return r.ru_utime.tv_sec + r.ru_stime.tv_sec +
           1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

//taken from minimap2
static inline long peakrss(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
#ifdef __linux__
	return r.ru_maxrss * 1024;
#else
	return r.ru_maxrss;
#endif
}

// Prints to the provided buffer a nice number of bytes (KB, MB, GB, etc)
//from https://www.mbeckler.org/blog/?p=114
static inline void print_size(const char* name, uint64_t bytes)
{
    const char* suffixes[7];
    suffixes[0] = "B";
    suffixes[1] = "KB";
    suffixes[2] = "MB";
    suffixes[3] = "GB";
    suffixes[4] = "TB";
    suffixes[5] = "PB";
    suffixes[6] = "EB";
    uint64_t s = 0; // which suffix to use
    double count = bytes;
    while (count >= 1024 && s < 7)
    {
        s++;
        count /= 1024;
    }
    if (count - floor(count) == 0.0)
        fprintf(stderr, "[%s] %s : %d %s\n", __func__ , name, (int)count, suffixes[s]);
    else
        fprintf(stderr, "[%s] %s : %.1f %s\n", __func__, name, count, suffixes[s]);
}

//replace u with t in a string
static inline void replace_char(char *str, char u, char t){
    while(*str){
        if(*str == u){
            *str = t;
        }
        str++;
    }
}

static inline float *signal_in_picoamps(slow5_rec_t *rec){
    int16_t* rawptr = rec->raw_signal;
    float range = rec->range;
    float digitisation = rec->digitisation;
    float offset = rec->offset;
    int32_t nsample = rec->len_raw_signal;

    // convert to pA
    float *current_signal = (float*)malloc(sizeof(float) * nsample);
    MALLOC_CHK(current_signal);

    float raw_unit = range / digitisation;
    for (int32_t j = 0; j < nsample; j++) {
        current_signal[j] = ((float)rawptr[j] + offset) * raw_unit;
    }

    return current_signal;
}

static inline int8_t drna_detect(slow5_file_t *sp){

    const slow5_hdr_t* hdr = sp->header;
    int8_t rna = 0;
    char *exp =slow5_hdr_get("experiment_type", 0, hdr);
    if(exp==NULL){
        WARNING("%s","experiment_type not found in SLOW5 header. Assuming genomic_dna");
        return 0;
    }
    if (strcmp(exp,"genomic_dna")==0){
        rna = 0;
    }else if (strcmp(exp,"rna")==0){
        rna = 1;
    } else {
        WARNING("Unknown experiment type: %s. Assuming genomic_dna", exp);
    }

    for(uint32_t  i=1; i < hdr->num_read_groups; i++){
        char *curr =slow5_hdr_get("experiment_type", i, hdr);
        if (strcmp(curr, exp)){
            WARNING("Experiment type mismatch: %s != %s in read group %d. Defaulted to %s", curr, exp, i, exp);
        }
    }
    return rna;
}

static inline void drna_mismatch(slow5_file_t *sp, int8_t rna){
    const char *expected = rna ? "rna" : "genomic_dna";
    const slow5_hdr_t* hdr = sp->header;
    for(uint32_t  i=0; i < hdr->num_read_groups; i++){
        char *curr =slow5_hdr_get("experiment_type", i, hdr);
        if (strcmp(curr, expected)){
            WARNING("Experiment type mismatch: %s != %s in read group %d. Double check for --rna.", curr, expected, i);
        }
    }
}


#define TO_PICOAMPS(RAW_VAL,DIGITISATION,OFFSET,RANGE) (((RAW_VAL)+(OFFSET))*((RANGE)/(DIGITISATION)))

#define SIGSIM_MEAN_VAL 104.6
#define SIGSIM_STDV_VAL 20.39
#define SIGSIM_WINDOW_SIZE 2000
#define SIGSIM_SIZE 1000

typedef struct {
    int64_t x;
    int64_t y;
} pair_t;

pair_t find_polya(float *raw, int64_t nsample, float top, float bot);
pair_t find_adaptor(slow5_rec_t *rec);
float meanf(float *x, int n);
float stdvf(float *x, int n);
float medianf(float *x, int n);
#endif
