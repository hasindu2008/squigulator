/* @file  sim.c
** Stupidly simple signal simulator
**
** @@
******************************************************************************/
#include <assert.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include "sigsim.h"
#include "error.h"
#include "ref.h"
#include "misc.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

typedef struct{
    char **ref_names;
    int32_t *ref_lengths;
    char **ref_seq;
    int num_ref;
    int64_t sum;
} ref_t;

typedef struct{
    double m;
    double s;
    int64_t x;
} nrng_t;

typedef struct{
    double a;
    double b;
    int64_t x;
} grng_t;

typedef struct {
    double digitisation;
    double sample_rate;
    //double bases_per_second;
    double range;
    double offset_mean;
    double offset_std;
    double median_before_mean;
    double median_before_std;
    double drift_mean;
    double drift_std;
} profile_t;

profile_t minion_r9_dna_prof = {
    .digitisation = 8192,
    .sample_rate = 4000,
    //.bases_per_second = 450,
    .range = 1443.030273,
    .offset_mean=13.7222605,
    .offset_std=10.25279688,
    .median_before_mean=200.815801,
    .median_before_std=20.48933762,
    .drift_mean=9.0, //this must be sample_rate/bases_per_second for now
    .drift_std=5.0
};
profile_t prom_r9_dna_prof = {
    .digitisation = 2048,
    .sample_rate = 4000,
    //.bases_per_second = 450,
    .range = 748.5801,
    .offset_mean=-237.4102,
    .offset_std=14.1575,
    .median_before_mean=214.2890337,
    .median_before_std=18.0127916,
    .drift_mean=9.0, //this must be sample_rate/bases_per_second for now
    .drift_std=5.0
};
profile_t minion_r9_rna_prof = {
    .digitisation = 8192,
    .sample_rate = 3012,
    //.bases_per_second = 70,
    .range = 1126.47,
    .offset_mean=4.65491888,
    .offset_std=4.115262472,
    .median_before_mean=242.6584118,
    .median_before_std=10.60230888,
    .drift_mean=43.0,  //this must be sample_rate/bases_per_second for now
    .drift_std=35.0
};
profile_t prom_r9_rna_prof = {
    .digitisation = 2048,
    .sample_rate = 3000,
    //.bases_per_second = 70,
    .range = 548.788269,
    .offset_mean=-231.9440589,
    .offset_std=12.87185278,
    .median_before_mean=238.5286796,
    .median_before_std=21.1871794,
    .drift_mean=43.0,  //this must be sample_rate/bases_per_second for now
    .drift_std=35.0
};

const char* polya = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
const char *adaptor_dna = "GGCGTCTGCTTGGGTGTTTAACCTTTTTTTTTTAATGTACTTCGTTCAGTTACGTATTGCT";
const char *adaptor_rna = "TGATGATGAGGGATAGACGATGGTTGTTTCTGTTGGTGCTGATATTGCTTTTTTTTTTTTTATGATGCAAGATACGCAC";

typedef struct{
    int8_t ideal;
    int8_t full_contigs;
    int8_t ideal_time;
    int8_t ideal_amp;

    int32_t rlen;
    int64_t seed;

    int8_t rna;
    const char *model_file;
    int8_t prefix;

} opt_sim_t;

typedef struct {

    int64_t ref_pos;
    int64_t rand_strand;
    nrng_t *rand_time;
    grng_t *rand_rlen;
    nrng_t *rand_offset;
    nrng_t *rand_median_before;
    profile_t profile;
    model_t *model;
    nrng_t **kmer_gen;
    uint32_t kmer_size;
    uint32_t num_kmer;

    //opt
    opt_sim_t opt;

} core_sim_t;

static struct option long_options[] = {
    {"verbose", required_argument, 0, 'v'},        //0 verbosity level [1]
    {"help", no_argument, 0, 'h'},                 //1
    {"version", no_argument, 0, 'V'},              //2
    {"output",required_argument, 0, 'o'},          //3 output to a file [stdout]
    {"ideal", no_argument, 0, 0},                  //4 ideal signals with no noise
    {"full-contigs", no_argument, 0, 0},           //5 simulate full contigs without breaking
    {"nreads", required_argument, 0, 'n'},         //6 number of reads
    {"fasta", required_argument, 0, 'q'},          //7 fasta perfect
    {"rlen", required_argument, 0, 'r'},           //8 median read
    {"seed", required_argument, 0, 0 },            //9 seed
    {"ideal-time", no_argument, 0, 0 },            //10 no time domain noise
    {"ideal-amp", no_argument, 0, 0 },             //11 no amplitude domain noise
    {"drift-mean", required_argument, 0, 0 },      //12 drift mean
    {"profile", required_argument, 0, 'm' },       //13 parameter profile
    {"kmer-model", required_argument, 0, 0},       //14 custom nucleotide k-mer model file
    {"prefix", no_argument, 0, 0},                 //15 prefix such as adaptor
    {0, 0, 0, 0}};


profile_t set_profile(char *prof_name, opt_sim_t *opt){
    opt->rna = 0;
    if(strcmp(prof_name, "dna-r9-min") == 0){
        return minion_r9_dna_prof;
    }else if(strcmp(prof_name, "dna-r9-prom") == 0){
        return prom_r9_dna_prof;
    }else if(strcmp(prof_name, "rna-r9-min") == 0){
        opt->rna = 1;
        return minion_r9_rna_prof;
    }else if(strcmp(prof_name, "rna-r9-prom") == 0){
        opt->rna = 1;
        return prom_r9_rna_prof;
    }else{
        ERROR("Unknown profile: %s\n", prof_name);
        exit(EXIT_FAILURE);
    }
}

static nrng_t* init_nrng(int64_t seed,double mean, double std){
    nrng_t *rng = (nrng_t *)malloc(sizeof(nrng_t));
    rng->m = mean;
    rng->s = std;
    rng->x = seed;
    return rng;
}

static grng_t* init_grng(int64_t seed,double alpha, double beta){
    grng_t *rng = (grng_t *)malloc(sizeof(grng_t));
    rng->a = alpha;
    rng->b = beta;
    rng->x = seed;
    return rng;
}

static void free_nrng(nrng_t *rng){
    free(rng);
}

static void free_grng(grng_t *rng){
    free(rng);
}

static double rng(int64_t *xp){
    int64_t x = *xp;
    int64_t x_new = (16807 * (x % 127773)) - (2836 * (x / 127773));
    if (x_new > 0) x = x_new; else x = x_new + 2147483647;
    *xp = x_new;
    return (double)x/2147483647;
}

static double nrng(nrng_t *r){
    double u = 0.0;
    double t = 0.0;
    while (u == 0.0) u = rng(&r->x);
    while (t == 0.0) t = 2.0 * 3.14159265 * rng(&r->x);
    double x = sqrt(-2.0 * log(u)) * cos(t);
    return ((x * r->s) + r->m);
}

static double grng(grng_t *r){
    double s = 0.0;
    for(int i=0; i<r->a; i++){
        s += -log(1-rng(&r->x));
    }
    return s*(r->b);
}

static void init_opt_sim(opt_sim_t *opt){
    opt->ideal = 0;
    opt->ideal_time = 0;
    opt->ideal_amp = 0;
    opt->full_contigs = 0;
    opt->rlen = 10000;
    opt->seed = 0;
    opt->rna = 0;
    opt->model_file = NULL;
    opt->prefix = 0;
}

static core_sim_t *init_core_sim(opt_sim_t opt, profile_t p){
    core_sim_t *core = (core_sim_t *)malloc(sizeof(core_sim_t));
    core->opt = opt;

    core->profile = p;
    model_t *m = core->model = (model_t*)malloc(sizeof(model_t) * MAX_NUM_KMER);
    uint32_t k = 0;
    if (opt.model_file) {
        k=read_model(core->model, opt.model_file, MODEL_TYPE_NUCLEOTIDE);
    } else {
        if(opt.rna){
            INFO("%s","builtin RNA nucleotide model loaded");
            k=set_model(core->model, MODEL_ID_RNA_NUCLEOTIDE);
        }
        else{
            k=set_model(core->model, MODEL_ID_DNA_NUCLEOTIDE);
        }
    }

    core->kmer_size = k;
    uint32_t n = core->num_kmer = (uint32_t)(1 << 2*k);

    core->ref_pos = opt.seed;
    core->rand_strand = opt.seed+1;
    core->rand_time = init_nrng(opt.seed+2, p.drift_mean, p.drift_std);
    core->rand_rlen = init_grng(opt.seed+3, 2.0, opt.rlen);
    core->rand_offset = init_nrng(opt.seed+4, p.offset_mean, p.offset_std);
    core->rand_median_before = init_nrng(opt.seed+5, p.median_before_mean, p.median_before_std);

    core->kmer_gen = (nrng_t **)malloc(sizeof(nrng_t *) * n);
    for (uint32_t i = 0; i < n; i++){
        core->kmer_gen[i] = init_nrng(opt.seed+i, m[i].level_mean, m[i].level_stdv);
    }

    return core;
}

void free_core_sim(core_sim_t *core){
    for (uint32_t i = 0; i < core->num_kmer; i++){
        free_nrng(core->kmer_gen[i]);
    }
    free(core->kmer_gen);

    free(core->model);
    free_nrng(core->rand_time);
    free_grng(core->rand_rlen);
    free_nrng(core->rand_offset);
    free_nrng(core->rand_median_before);
    free(core);
}


//todo can be optimised for memory by better encoding if necessary
static ref_t *load_ref(const char *genome){

    gzFile fp;
    kseq_t *seq;
    int l;
    fp = gzopen(genome, "r");
    F_CHK(fp,genome);
    seq = kseq_init(fp);
    MALLOC_CHK(seq);


    ref_t *ref = (ref_t *) malloc(sizeof(ref_t));

    int c = 1;
    ref->ref_lengths = (int32_t *) malloc(sizeof(int32_t));
    ref->ref_names = (char **) malloc(sizeof(char *));
    ref->ref_seq = (char **) malloc(sizeof(char *));
    ref->sum = 0;

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
        strcpy(ref->ref_names[i], seq->name.s);
        ref->ref_seq[i] = (char *) malloc((l+1)*sizeof(char));
        strcpy(ref->ref_seq[i], seq->seq.s);

        ref->sum += l;
        i++;

    }

    ref->num_ref = i;

    kseq_destroy(seq);
    gzclose(fp);

    return ref;

}

static void free_ref_sim(ref_t *ref){

    for(int i=0;i<ref->num_ref;i++){
        free(ref->ref_names[i]);
        free(ref->ref_seq[i]);
    }

    free(ref->ref_lengths);
    free(ref->ref_names);
    free(ref->ref_seq);
    free(ref);
}



static void set_header_attributes(slow5_file_t *sp, int8_t rna){

    slow5_hdr_t *header=sp->header;

    //add a header group attribute called run_id
    if (slow5_hdr_add("run_id", header) < 0){
        ERROR("%s","Error adding run_id attribute\n");
        exit(EXIT_FAILURE);
    }
    //add another header group attribute called asic_id
    if (slow5_hdr_add("asic_id", header) < 0){
        ERROR("%s","Error adding asic_id attribute\n");
        exit(EXIT_FAILURE);
    }
    //add another header group attribute called asic_id
    if (slow5_hdr_add("exp_start_time", header) < 0){
        ERROR("%s","Error adding asic_id attribute\n");
        exit(EXIT_FAILURE);
    }
    //add another header group attribute called flow_cell_id
    if (slow5_hdr_add("flow_cell_id", header) < 0){
        ERROR("%s","Error adding flow_cell_id attribute\n");
        exit(EXIT_FAILURE);
    }
    //add another header group attribute called experiment_type
    if (slow5_hdr_add("experiment_type", header) < 0){
        ERROR("%s","Error adding experiment_type attribute\n");
        exit(EXIT_FAILURE);
    }


    //set the run_id attribute to "run_0" for read group 0
    if (slow5_hdr_set("run_id", "run_0", 0, header) < 0){
        ERROR("%s","Error setting run_id attribute in read group 0\n");
        exit(EXIT_FAILURE);
    }
    //set the asic_id attribute to "asic_0" for read group 0
    if (slow5_hdr_set("asic_id", "asic_id_0", 0, header) < 0){
        ERROR("%s","Error setting asic_id attribute in read group 0\n");
        exit(EXIT_FAILURE);
    }
    //set the exp_start_time attribute to "2022-07-20T00:00:00Z" for read group 0
    if (slow5_hdr_set("exp_start_time", "2022-07-20T00:00:00Z", 0, header) < 0){
        ERROR("%s","Error setting exp_start_time attribute in read group 0\n");
        exit(EXIT_FAILURE);
    }
    //set the flow_cell_id attribute to "FAN00000" for read group 0
    if (slow5_hdr_set("flow_cell_id", "FAN00000", 0, header) < 0){
        ERROR("%s","Error setting flow_cell_id attribute in read group 0\n");
        exit(EXIT_FAILURE);
    }
    //set the experiment_type attribute to genomic_dna or rna for read group 0
    const char* experiment_type = rna ? "rna" : "genomic_dna" ;
    if (slow5_hdr_set("experiment_type", experiment_type, 0, header) < 0){
        ERROR("%s","Error setting experiment_type attribute in read group 0\n");
        exit(EXIT_FAILURE);
    }

}

static void set_header_aux_fields(slow5_file_t *sp){

    //add auxilliary field: channel number
    if (slow5_aux_add("channel_number", SLOW5_STRING, sp->header) < 0){
        ERROR("%s","Error adding channel_number auxilliary field\n");
        exit(EXIT_FAILURE);
    }

    //add auxilliary field: median_before
    if (slow5_aux_add("median_before", SLOW5_DOUBLE, sp->header) < 0) {
        ERROR("%s","Error adding median_before auxilliary field\n");
        exit(EXIT_FAILURE);
    }

    //add auxilliary field: read_number
    if(slow5_aux_add("read_number", SLOW5_INT32_T, sp->header) < 0){
        ERROR("%s","Error adding read_number auxilliary field\n");
        exit(EXIT_FAILURE);
    }
    //add auxilliary field: start_mux
    if(slow5_aux_add("start_mux", SLOW5_UINT8_T, sp->header) < 0){
        ERROR("%s","Error adding start_mux auxilliary field\n");
        exit(EXIT_FAILURE);
    }
    //add auxilliary field: start_time
    if(slow5_aux_add("start_time", SLOW5_UINT64_T, sp->header) < 0){
        ERROR("%s","Error adding start_time auxilliary field\n");
        exit(EXIT_FAILURE);
    }

}

static void set_record_primary_fields(profile_t *profile, slow5_rec_t *slow5_record, char *read_id, double offset, int64_t len_raw_signal, int16_t *raw_signal){

    slow5_record -> read_id = read_id;
    slow5_record-> read_id_len = strlen(slow5_record -> read_id);
    slow5_record -> read_group = 0;
    slow5_record -> digitisation = profile->digitisation;
    slow5_record -> offset = offset;
    slow5_record -> range = profile->range;
    slow5_record -> sampling_rate = profile->sample_rate;
    slow5_record -> len_raw_signal = len_raw_signal;
    slow5_record -> raw_signal = raw_signal;

}

static void set_record_aux_fields(slow5_rec_t *slow5_record, slow5_file_t *sp, double median_before, int32_t read_number, uint64_t start_time){

    const char *channel_number = "0";
    uint8_t start_mux = 0;

    if(slow5_aux_set_string(slow5_record, "channel_number", channel_number, sp->header) < 0){
        ERROR("%s","Error setting channel_number auxilliary field\n");
        exit(EXIT_FAILURE);
    }

    if(slow5_aux_set(slow5_record, "median_before", &median_before, sp->header) < 0){
        ERROR("%s","Error setting median_before auxilliary field\n");
        exit(EXIT_FAILURE);
    }

    if(slow5_aux_set(slow5_record, "read_number", &read_number, sp->header) < 0){
        ERROR("%s","Error setting read_number auxilliary field\n");
        exit(EXIT_FAILURE);
    }

    if(slow5_aux_set(slow5_record, "start_mux", &start_mux, sp->header) < 0){
        ERROR("%s","Error setting start_mux auxilliary field\n");
        exit(EXIT_FAILURE);
    }

    if(slow5_aux_set(slow5_record, "start_time", &start_time, sp->header) < 0){
        ERROR("%s","Error setting start_time auxilliary field\n");
        exit(EXIT_FAILURE);
    }


}


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


int16_t * gen_sig_core_seq(core_sim_t *core, int16_t *raw_signal, int64_t* n, int64_t *c, double offset, const char *read, int32_t len){

    profile_t *profile = &core->profile;
    uint32_t kmer_size = core->kmer_size;
    model_t *pore_model = core->model;

    int8_t ideal = core->opt.ideal;
    int8_t ideal_time = core->opt.ideal_time;
    int8_t ideal_amp = core->opt.ideal_amp;
    int sps = (int)profile->drift_mean;

    int64_t n_kmers = len-kmer_size+1;


    for (int i=0; i< n_kmers; i++){
        uint32_t kmer_rank = get_kmer_rank(read+i, kmer_size);
        if(!(ideal || ideal_time)){
            sps = round(nrng(core->rand_time));
        }
        for(int j=0; j<sps; j++){
            if(*n==*c){
                *c *= 2;
                raw_signal = (int16_t *)realloc(raw_signal, *c*sizeof(int16_t));
            }
            float s = 0;
            if(ideal || ideal_amp){
                s=pore_model[kmer_rank].level_mean;
            } else {
                s=nrng(core->kmer_gen[kmer_rank]);
            }
            raw_signal[*n] = s*(profile->digitisation)/(profile->range)-(offset);
            *n = *n+1;
        }
    }
    return raw_signal;

}

int16_t *gen_prefix_rna(core_sim_t *core, int16_t *raw_signal, int64_t* n, int64_t *c, double offset){

    profile_t *profile = &core->profile;
    //float s;
    //polya
    // raw_signal=gen_sig_core_seq(core, raw_signal, n, c, offset, polya, strlen(polya));

    //adaptor
    int st = *n-(strlen(adaptor_rna)*(int)profile->drift_mean);
    // raw_signal=gen_sig_core_seq(core, raw_signal, n, c, offset, adaptor_rna, strlen(adaptor_rna));
    int end = *n;
    //fprintf(stderr, "adaptor_rna: %d %d\n", st, end);
    int16_t off = 30*(profile->digitisation)/(profile->range);
    for(int i=st; i<end; i++){
        raw_signal[i] = raw_signal[i]-off;
    }

    const char *stall = "AAAAAGAAAAAACCCCCCCCCCCCCCCCCC";
    raw_signal=gen_sig_core_seq(core, raw_signal, n, c, offset, stall, strlen(stall));


    return raw_signal;
}

char *attach_prefix(core_sim_t *core, const char *read, int32_t *len){

    char *s = (char *)read;
    if(core->opt.rna){
        int polya_len = strlen(polya);
        int adapt_len = strlen(adaptor_rna);
        char *seq = malloc((*len+polya_len+adapt_len+1)*sizeof(char));
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
        strncpy(seq, stall, stall_len);
        strncpy(seq+stall_len, adaptor_dna, adapt_len);
        strncpy(seq+stall_len+adapt_len, read, *len);
        *len = *len+stall_len+adapt_len;
        seq[*len] = '\0';
        s = seq;
    }
    return s;
}

int16_t *gen_sig_core(core_sim_t *core, const char *read, int32_t len, double *offset, double *median_before, int64_t *len_raw_signal){

    profile_t *profile = &core->profile;
    uint32_t kmer_size = core->kmer_size;
    //model_t *pore_model = core->model;

    int8_t ideal = core->opt.ideal;
    //int8_t ideal_time = core->opt.ideal_time;
    //int8_t ideal_amp = core->opt.ideal_amp;

    int64_t n_kmers = len-kmer_size+1;
    int64_t n=0;
    int sps = (int)profile->drift_mean;

    int64_t c = n_kmers * sps + 2000;
    int16_t *raw_signal = (int16_t *)malloc(c*sizeof(int16_t));

    if(ideal) {
        *offset = profile->offset_mean;
        *median_before = profile->median_before_mean;
    } else {
        *offset = nrng(core->rand_offset);
        *median_before = nrng(core->rand_median_before);
    }

    //todo adaptor sequence and prefix
    // if(core->opt.prefix && !core->opt.rna){
    //     gen_prefix_dna(raw_signal,&n,&c, profile, *offset);
    // }

    char *tmpread = NULL;
    if(core->opt.prefix){
        read = tmpread = attach_prefix(core, read, &len);
    }
    //fprintf(stderr, "read: %s\n", read);

    raw_signal = gen_sig_core_seq(core, raw_signal, &n, &c, *offset, read, len);

    if(core->opt.prefix && core->opt.rna){
        raw_signal=gen_prefix_rna(core, raw_signal,&n,&c, *offset);
    }
    if(tmpread){
        free(tmpread);
    }
    assert(n<=c);

    *len_raw_signal = n;
    return raw_signal;
}


static inline int16_t *gen_sig(core_sim_t *core, const char *read, int32_t len, double *offset, double *median_before, int64_t *len_raw_signal, int8_t rna){
    int16_t *sig = gen_sig_core(core, read, len, offset, median_before, len_raw_signal);
    if(rna){
        for(int i=0; i<*len_raw_signal/2; i++){
            int16_t tmp = sig[i];
            sig[i] = sig[*len_raw_signal-1-i];
            sig[*len_raw_signal-1-i] = tmp;
        }
    }
    return sig;
}


static inline int32_t is_bad_read(char *seq, int32_t len){
    if (len < 200){
        return -1;
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

char *gen_read_dna(core_sim_t *core, ref_t *ref, char **ref_id, int32_t *ref_pos, int32_t *rlen, char *c){

    char *seq = NULL;
    while(1){

        int len = grng(core->rand_rlen);
        int64_t ref_sum_pos = round(rng(&core->ref_pos)*ref->sum); //serialised pos
        assert(ref_sum_pos <= ref->sum);
        int64_t s = 0;
        int seq_i = 0;
        for(seq_i=0; seq_i<ref->num_ref; seq_i++){ //check manually if logic is right
            s+=ref->ref_lengths[seq_i];
            if(s>=ref_sum_pos){
                *ref_id = ref->ref_names[seq_i];
                *ref_pos = ref_sum_pos-s+ref->ref_lengths[seq_i];
                break;
            }
        }
        assert(s<=ref->sum);

        int64_t strand = round(rng(&core->rand_strand));
        *c = strand ? '+' : '-';

        seq= (char *)malloc((len+1)*sizeof(char));
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
                VERBOSE("Too short read: %d. %s:%d-%d. Trying again!",200,*ref_id,*ref_pos,*ref_pos+*rlen);
            } else{
                VERBOSE("Too many Ns in read: %d. %s:%d-%d. Trying again!",nc,*ref_id,*ref_pos,*ref_pos+*rlen);
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


char *gen_read_rna(core_sim_t *core, ref_t *ref, char **ref_id, int32_t *ref_pos, int32_t *rlen, char *c){


    char *seq = NULL;
    while(1){

        //int len = grng(core->rand_rlen); //for now the whole transcript is is simulated

        int seq_i = round(rng(&core->ref_pos)*(ref->num_ref-1)); //random transcript
        int len = ref->ref_lengths[seq_i];
        *ref_pos=0;
        *ref_id = ref->ref_names[seq_i];

        //int64_t strand = round(rng(&core->rand_strand));
        *c = '+'; //strand is alwats plus

        seq= (char *)malloc((len+1)*sizeof(char));
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
                VERBOSE("Too short read: %d. %s:%d-%d. Trying again!",200,*ref_id,*ref_pos,*ref_pos+*rlen);
            } else{
                VERBOSE("Too many Ns in read: %d. %s:%d-%d. Trying again!",nc,*ref_id,*ref_pos,*ref_pos+*rlen);
            }
            free(seq);
        }

    }

    return seq;

}



static inline char *gen_read(core_sim_t *core, ref_t *ref, char **ref_id, int32_t *ref_pos, int32_t *rlen, char *c, int8_t rna){
    return rna ? gen_read_rna(core, ref, ref_id, ref_pos, rlen, c): gen_read_dna(core, ref, ref_id, ref_pos, rlen, c);
}

int sim_main(int argc, char* argv[], double realtime0) {

    const char* optstring = "o:hVn:q:r:x:";

    int longindex = 0;
    int32_t c = -1;

    FILE *fp_help = stderr;
    char *output_file = NULL;
    char *fasta = NULL;

    opt_sim_t opt;
    init_opt_sim(&opt);

    profile_t p = prom_r9_dna_prof;

    int nreads = 4000;

    //parse the user args
    while ((c = getopt_long(argc, argv, optstring, long_options, &longindex)) >= 0) {
        if (c=='V'){
            fprintf(stdout,"sigsim %s\n",SIGSIM_VERSION);
            exit(EXIT_SUCCESS);
        } else if (c=='h'){
            fp_help = stdout;
        } else if (c=='x'){
            p = set_profile(optarg, &opt);
        } else if(c == 'o'){
            output_file=optarg;
        } else if (c == 0 && longindex == 4){   //generate ideal signal
            opt.ideal = 1;
        } else if (c == 0 && longindex == 5){  //generate signal for complete contigs
            opt.full_contigs =1 ;
        } else if (c == 'n'){
            nreads = atoi(optarg);
        } else if (c == 'q'){
            fasta = optarg;
        } else if (c == 'r'){
            opt.rlen = atoi(optarg);
        } else if (c == 0 && longindex == 9){  //seed
            opt.seed = atoi(optarg);
        } else if (c == 0 && longindex == 10){ //ideal-time
            opt.ideal_time = 1;
        } else if (c == 0 && longindex == 11){ //ideal-amp
            opt.ideal_amp = 1;
        } else if (c == 0 && longindex == 12){ //drift-mean
            p.drift_mean = atof(optarg);
        } else if (c == 0 && longindex == 14) { //custom nucleotide model file
            opt.model_file = optarg;
        } else if (c == 0 && longindex == 15) { //prefix
            opt.prefix = 1;
        }
    }

    if (argc-optind<1 || output_file==NULL ||  fp_help == stdout) {
        fprintf(fp_help,"Usage: sigsim [OPTIONS] ref.fa -o out_signal.blow5\n");
        fprintf(fp_help,"\nbasic options:\n");
        fprintf(fp_help,"   -o FILE                    SLOW5/BLOW5 file to write\n");
        fprintf(fp_help,"   -x STR                     parameter profile (always applied before other options) [dna-r9-prom]\n");
        fprintf(fp_help,"                              e.g., dna-r9-min, dna-r9-prom, rna-r9-min, rna-r9-prom\n");
        fprintf(fp_help,"   -n INT                     Number of reads to simulate (ignored if --full-contigs) [%d]\n", nreads);
        fprintf(fp_help,"   -q FILE                    FASTA file to write simulated reads with no errors\n");
        fprintf(fp_help,"   -r INT                     Median read length (estimate only, ignored if --full-contigs or dRNA) [%d]\n",opt.rlen);
        fprintf(fp_help,"   -h                         help\n");
        fprintf(fp_help,"   --full-contigs             Generate signals for complete contigs (ineffective for gDNA)\n");
        fprintf(fp_help,"   --ideal                    Generate ideal signals with no noise\n");
        fprintf(fp_help,"   --version                  print version\n");

        fprintf(fp_help,"\nadvanced options:\n");
        fprintf(fp_help,"   --seed INT                 Seed or random generators (if 0, will be autogenerated) [%ld]\n",opt.seed);
        fprintf(fp_help,"   --ideal-time               Generate signals with no time domain noise\n");
        fprintf(fp_help,"   --ideal-amp                Generate signals with no amplitiude domain noise\n");
        fprintf(fp_help,"   --drift-mean FLOAT         Mean of drift rate [%f]\n",p.drift_mean);
        fprintf(fp_help,"   --kmer-model FILE          custom nucleotide k-mer model file (format similar to https://github.com/hasindu2008/f5c/blob/master/test/r9-models/r9.4_450bps.nucleotide.6mer.template.model)\n");
        fprintf(fp_help,"   --prefix                   generate prefixes such as adaptor (and polya for RNA)\n");
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }

    //todo check args
    //-n incompatible with --full-contigs
    //-r incompatible with --full-contigs

    if (opt.seed == 0){
        opt.seed = realtime0;
        VERBOSE("Using random seed: %ld\n",opt.seed);
    }

    char *refname = argv[optind];
    ref_t *ref = load_ref(refname);

    slow5_file_t *sp = slow5_open(output_file, "w");
    if(sp==NULL){
        ERROR("Error opening file %s!\n",output_file);
        exit(EXIT_FAILURE);
    }

    set_header_attributes(sp, opt.rna);
    set_header_aux_fields(sp);

    if(slow5_hdr_write(sp) < 0){
        ERROR("%s","Error writing header!\n");
        exit(EXIT_FAILURE);
    }

    int n = nreads;
    double median_before = 0;
    int64_t n_samples = 0;
    double offset = 0;
    int64_t len_raw_signal =0;
    char *rid = NULL;
    char *seq = NULL;
    int32_t rlen = 0;
    char strand = '+';
    int32_t ref_pos_st = 0;
    int32_t ref_pos_end = 0;

    if(opt.full_contigs){
        n = ref->num_ref;
    }

    core_sim_t *core = init_core_sim(opt, p);

    FILE *fp_fasta = NULL;
    if(fasta){
        fp_fasta = fopen(fasta,"w");
        if(fp_fasta==NULL){
            fprintf(stderr, "Error opening file %s\n",fasta);
            exit(EXIT_FAILURE);
        }
    }

    for(int i=0;i<n;i++){

        slow5_rec_t *slow5_record = slow5_rec_init();
        if(slow5_record == NULL){
            fprintf(stderr,"Could not allocate space for a slow5 record.");
            exit(EXIT_FAILURE);
        }

        if(opt.full_contigs){
            rid = ref->ref_names[i];
            rlen = ref->ref_lengths[i];
            seq = ref->ref_seq[i];
            strand = '+';
            ref_pos_st = 0;
            ref_pos_end = rlen;
        } else {
            seq=gen_read(core, ref, &rid, &ref_pos_st, &rlen, &strand, opt.rna);
            ref_pos_end = ref_pos_st+rlen;
        }
        int16_t *raw_signal=gen_sig(core, seq, rlen, &offset, &median_before, &len_raw_signal, opt.rna);

        char *read_id= (char *)malloc(sizeof(char)*(1000));
        sprintf(read_id,"S1_%d!%s!%d!%d!%c",i+1, rid, ref_pos_st, ref_pos_end, strand);
        if(fasta){
            fprintf(fp_fasta,">%s\n",read_id);
            fprintf(fp_fasta,"%s\n",seq);
        }

        set_record_primary_fields(&core->profile, slow5_record, read_id, offset, len_raw_signal, raw_signal);
        set_record_aux_fields(slow5_record, sp, median_before, i, n_samples);
        n_samples+=len_raw_signal;

        if (slow5_write(slow5_record, sp) < 0){
            fprintf(stderr,"Error writing record!\n");
            exit(EXIT_FAILURE);
        }

        if(!opt.full_contigs){
            free(seq);
        }
        slow5_rec_free(slow5_record);

        if((i+1)%1000==0) {
            VERBOSE("%d reads done",i+1);
        }

    }

    if(fasta){
        fclose(fp_fasta);
    }

    free_core_sim(core);

    slow5_close(sp);
    free_ref_sim(ref);

    return 0;
}
