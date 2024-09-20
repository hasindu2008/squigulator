/* @file  sim.c
** nanopore signal simulator
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
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>

#include "version.h"
#include "sq.h"
#include "ref.h"
#include "format.h"
#include "seq.h"
#include "rand.h"
#include "error.h"
#include "misc.h"

void set_header_attributes(slow5_file_t *sp, int8_t rna, int8_t r10, double sample_frequency);
void set_header_aux_fields(slow5_file_t *sp, int8_t ont_friendly);
void set_record_primary_fields(profile_t *profile, slow5_rec_t *slow5_record, char *read_id, double offset, int64_t len_raw_signal, int16_t *raw_signal);
void set_record_aux_fields(slow5_rec_t *slow5_record, slow5_file_t *sp, double median_before, int32_t read_number, uint64_t start_time, int8_t ont_friendly);
uint32_t read_model(model_t* model, const char* file, uint32_t type);
uint32_t set_model(model_t* model, uint32_t model_id);


///////////////////////////////////////////// Profiles /////////////////////////////////////////////

profile_t minion_r9_dna_prof = {
    .digitisation = 8192,
    .sample_rate = 4000,
    .bps = 450,
    .range = 1443.030273,
    .offset_mean=13.7222605,
    .offset_std=10.25279688,
    .median_before_mean=200.815801,
    .median_before_std=20.48933762,
    .dwell_mean=9.0, //this must be sample_rate/bps for now
    .dwell_std=4.0
};
profile_t prom_r9_dna_prof = {
    .digitisation = 2048,
    .sample_rate = 4000,
    .bps = 450,
    .range = 748.5801,
    .offset_mean=-237.4102,
    .offset_std=14.1575,
    .median_before_mean=214.2890337,
    .median_before_std=18.0127916,
    .dwell_mean=9.0, //this must be sample_rate/bps for now
    .dwell_std=4.0
};
profile_t minion_r9_rna_prof = {
    .digitisation = 8192,
    .sample_rate = 3012,
    .bps = 70,
    .range = 1126.47,
    .offset_mean=4.65491888,
    .offset_std=4.115262472,
    .median_before_mean=242.6584118,
    .median_before_std=10.60230888,
    .dwell_mean=43.0,  //this must be sample_rate/bps for now
    .dwell_std=35.0
};
profile_t prom_r9_rna_prof = {
    .digitisation = 2048,
    .sample_rate = 3000,
    .bps = 70,
    .range = 548.788269,
    .offset_mean=-231.9440589,
    .offset_std=12.87185278,
    .median_before_mean=238.5286796,
    .median_before_std=21.1871794,
    .dwell_mean=43.0,  //this must be sample_rate/bps for now
    .dwell_std=35.0
};
profile_t prom_r10_dna_prof = {
    .digitisation = 2048,
    .sample_rate = 5000,
    .bps = 400,
    .range = 281.345551,
    .offset_mean=-127.5655735,
    .offset_std=19.377283387665,
    .median_before_mean=189.87607393756,
    .median_before_std=15.788097978713,
    .dwell_mean=13.0, //this must be sample_rate/bps for now
    .dwell_std=4.0
};
profile_t minion_r10_dna_prof = {
    .digitisation = 8192,
    .sample_rate = 5000,
    .bps = 400,
    .range = 1536.598389,
    .offset_mean=13.380569389019,
    .offset_std=16.311471649012,
    .median_before_mean=202.15407438804,
    .median_before_std=13.406139241768,
    .dwell_mean=13.0, //this must be sample_rate/bps for now
    .dwell_std=4.0
};
profile_t prom_rna004_rna_prof = {
    .digitisation = 2048,
    .sample_rate = 4000,
    .bps = 130,
    .range =  299.432068,
    .offset_mean=-259.421128,
    .offset_std=16.010841823643,
    .median_before_mean=205.63935594369,
    .median_before_std=8.3994882799157,
    .dwell_mean=31.0, //this must be sample_rate/bps for now
    .dwell_std=0.0
};
profile_t minion_rna004_rna_prof = {
    .digitisation = 8192,
    .sample_rate = 4000,
    .bps = 130,
    .range = 1437.976685,
    .offset_mean=12.47686423863,
    .offset_std=10.442126577137,
    .median_before_mean=205.08496731088,
    .median_before_std=8.6671292866233,
    .dwell_mean=31.0, //this must be sample_rate/bps for now
    .dwell_std=0.0
};

static inline profile_t set_profile(char *prof_name, opt_t *opt){
    //opt->rna = 0;
    if(strcmp(prof_name, "dna-r9-min") == 0){
        return minion_r9_dna_prof;
    }else if(strcmp(prof_name, "dna-r9-prom") == 0){
        return prom_r9_dna_prof;
    }else if(strcmp(prof_name, "rna-r9-min") == 0){
        opt->flag |= SQ_RNA;
        return minion_r9_rna_prof;
    }else if(strcmp(prof_name, "rna-r9-prom") == 0){
        opt->flag |= SQ_RNA;
        return prom_r9_rna_prof;
    }else if(strcmp(prof_name, "dna-r10-prom") == 0){
        INFO("%s", "dna-r10-prom is 5kHz from squigulator v0.3.0 onwards. Specify --sample-rate 4000 for old 4kHz.")
        WARNING("%s","Parameters and models for dna-r10-prom 5khz are still crude. If you have good modification-free data, please share!");
        opt->flag |= SQ_R10;
        return prom_r10_dna_prof;
    }else if(strcmp(prof_name, "dna-r10-min") == 0){
        INFO("%s", "dna-r10-min is 5kHz from squigulator v0.3.0 onwards. Specify --sample-rate 4000 for old 4kHz.")
        WARNING("%s","Parameters and models for dna-r10-min 5khz are still crude. If you have good modification-free data, please share!");
        opt->flag |= SQ_R10;
        return minion_r10_dna_prof;
    }else if(strcmp(prof_name, "rna004-min") == 0){
        WARNING("%s","Parameters and models for rna004-min are still crude. If you have good IVT data, please share!");
        opt->flag |= SQ_R10;
        opt->flag |= SQ_RNA;
        return minion_rna004_rna_prof;
    }else if(strcmp(prof_name, "rna004-prom") == 0){
        WARNING("%s","Parameters and models for rna004-prom are still crude. If you have good IVT data, please share!");
        opt->flag |= SQ_R10;
        opt->flag |= SQ_RNA;
        return prom_rna004_rna_prof;
    }else{
        ERROR("Unknown profile: %s\n", prof_name);
        exit(EXIT_FAILURE);
    }
}

static void print_model_stat(profile_t p){
    VERBOSE("digitisation: %.1f; sample_rate: %.1f; range: %.1f; offset_mean: %.1f; offset_std: %.1f; dwell_mean: %.1f; dwell_std: %.1f",
        p.digitisation,p.sample_rate,p.range,p.offset_mean,p.offset_std,p.dwell_mean,p.dwell_std);
}

///////////////////////////////////////////// CORE and DB /////////////////////////////////////////////

static void init_opt(opt_t *opt){
    //opt->ideal = 0;
    //opt->ideal_time = 0;
    //opt->ideal_amp = 0;
    //opt->full_contigs = 0;
    opt->rlen = 10000;
    opt->seed = 0;
    //opt->rna = 0;
    opt->model_file = NULL;
    //opt->prefix = 0;
    opt->flag = 0;
    opt->num_thread = 8;
    opt->batch_size = 1000;
    opt->amp_noise = 1;
    opt->meth_model_file = NULL;
    opt->meth_freq = NULL;
}

static void init_rand(core_t *core){

    uint32_t n = core->num_kmer;
    model_t *m = core->model;
    opt_t opt = core->opt;
    int t = opt.num_thread;
    profile_t p = core->profile;

    core->ref_pos  = (int64_t *) malloc(t*sizeof(int64_t));             MALLOC_CHK(core->ref_pos);
    core->rand_strand = (int64_t *) malloc(t*sizeof(int64_t));          MALLOC_CHK(core->rand_strand);
    core->rand_time = (nrng_t **) malloc(t*sizeof(nrng_t *));           MALLOC_CHK(core->rand_time);
    core->rand_rlen = (grng_t **) malloc(t*sizeof(grng_t *));           MALLOC_CHK(core->rand_rlen);
    core->rand_offset = (nrng_t **) malloc(t*sizeof(nrng_t *));         MALLOC_CHK(core->rand_offset);
    core->rand_median_before = (nrng_t **) malloc(t*sizeof(nrng_t *));  MALLOC_CHK(core->rand_median_before);
    core->kmer_gen = (nrng_t ***) malloc(t*sizeof(nrng_t **));          MALLOC_CHK(core->kmer_gen);

    if(core->opt.meth_freq){
        m = core->cpgmodel;
        core->rand_meth = (int64_t *) malloc(t*sizeof(int64_t));         MALLOC_CHK(core->rand_meth);
    } else {
        core->rand_meth = NULL;
    }

    int64_t seed = opt.seed;
    for(int i=0; i<t; i++){
        core->ref_pos[i] = seed;
        core->rand_strand[i] = seed+1;
        core->rand_time[i] = init_nrng(seed+2, p.dwell_mean, p.dwell_std);
        core->rand_rlen[i] = init_grng(seed+3, 2.0, opt.rlen/2);
        core->rand_offset[i] = init_nrng(seed+4, p.offset_mean, p.offset_std);
        core->rand_median_before[i] = init_nrng(seed+5, p.median_before_mean, p.median_before_std);

        core->kmer_gen[i] = (nrng_t **)malloc(sizeof(nrng_t *) * n); MALLOC_CHK(core->kmer_gen[i]);
        for (uint32_t j = 0; j < n; j++){
            core->kmer_gen[i][j] = init_nrng(seed+j, m[j].level_mean, m[j].level_stdv*opt.amp_noise);
        }

        if(core->opt.meth_freq){
            core->rand_meth[i] = seed+6;
        }

        seed += (n+10);
    }
}

static core_t *init_core(opt_t opt, profile_t p, char *refname, char *output_file, char *fasta, char *paf, char *sam, char *trans_count){
    core_t *core = (core_t *)malloc(sizeof(core_t)); MALLOC_CHK(core);
    core->opt = opt;

    core->profile = p;
    core->model = (model_t*)malloc(sizeof(model_t) * MAX_NUM_KMER); //todo not very memory efficient - do dynamically
    MALLOC_CHK(core->model);

    uint32_t k = 0;
    if (opt.model_file) {
        k=read_model(core->model, opt.model_file, MODEL_TYPE_NUCLEOTIDE);
    } else {
        if(opt.flag & SQ_R10){ //R10 or RNA004
            if(opt.flag & SQ_RNA){
                INFO("%s","builtin RNA004 nucleotide model loaded");
                k=set_model(core->model, MODEL_ID_RNA_RNA004_NUCLEOTIDE);
            } else {
                INFO("%s","builtin DNA R10 nucleotide model loaded");
                k=set_model(core->model, MODEL_ID_DNA_R10_NUCLEOTIDE);
            }
        }
        else { //R9
            if(opt.flag & SQ_RNA){
                INFO("%s","builtin RNA R9 nucleotide model loaded");
                k=set_model(core->model, MODEL_ID_RNA_R9_NUCLEOTIDE);
             }
            else{
                INFO("%s","builtin DNA R9 nucleotide model loaded");
                k=set_model(core->model, MODEL_ID_DNA_R9_NUCLEOTIDE);
            }
        }
    }

    core->kmer_size = k;
    core->num_kmer = (uint32_t)(1 << 2*k);

    uint32_t kmer_size_meth=0;
    if(opt.meth_freq){
        core->cpgmodel = (model_t*)malloc(sizeof(model_t) * MAX_NUM_KMER_METH);
        MALLOC_CHK(core->cpgmodel);

        if (opt.meth_model_file) {
            kmer_size_meth=read_model(core->cpgmodel, opt.meth_model_file, MODEL_TYPE_METH);
        } else {

            if(opt.flag & SQ_RNA){
                ERROR("%s","Methylation simulation is not supported for RNA data.");
                exit(EXIT_FAILURE);
            }

            if(opt.flag & SQ_R10){
                ERROR("%s","Methylation simulation is not yet supported for R10 data.");
                exit(EXIT_FAILURE);
                INFO("%s","builtin DNA R10 cpg model loaded");
                kmer_size_meth=set_model(core->cpgmodel, MODEL_ID_DNA_R10_CPG);
                WARNING("%s","Methylation simulation for R10 is still crude. If you have good control data, please share!");
            } else {
                INFO("%s","builtin DNA R9 cpg model loaded");
                kmer_size_meth=set_model(core->cpgmodel, MODEL_ID_DNA_R9_CPG);
            }
        }
        if(k != kmer_size_meth){
            ERROR("The k-mer size of the nucleotide model (%d) and the methylation model (%d) should be the same.",k,kmer_size_meth);
            exit(EXIT_FAILURE);
        }
        core->num_kmer = (uint32_t)pow(5,k);

    } else {
        core->cpgmodel = NULL;
    }


    init_rand(core);

    core->ref = load_ref(refname);
    if(trans_count!=NULL){
        load_trans_count(trans_count,core->ref);
        assert((opt.flag & SQ_RNA) || (opt.flag & SQ_CDNA));
    }
    if(opt.meth_freq!=NULL){
        load_meth_freq(opt.meth_freq, core->ref);
    }

    core->sp = slow5_open(output_file, "w");
    if(core->sp==NULL){
        ERROR("Error opening file %s!\n",output_file);
        exit(EXIT_FAILURE);
    }

    set_header_attributes(core->sp, opt.flag & SQ_RNA ? 1 : 0, opt.flag & SQ_R10 ? 1 : 0, p.sample_rate);
    set_header_aux_fields(core->sp, opt.flag & SQ_ONT ? 1 : 0);

    if(slow5_hdr_write(core->sp) < 0){
        ERROR("%s","Error writing header!\n");
        exit(EXIT_FAILURE);
    }

    core->fp_fasta = NULL;
    if(fasta){
        core->fp_fasta = fopen(fasta,"w");
        if(core->fp_fasta==NULL){
            fprintf(stderr, "Error opening file %s\n",fasta);
            exit(EXIT_FAILURE);
        }
    }
    core->fp_paf = NULL;
    if(paf){
        core->fp_paf = fopen(paf,"w");
        if(core->fp_paf==NULL){
            fprintf(stderr, "Error opening file %s\n",paf);
            exit(EXIT_FAILURE);
        }
    }
    core->fp_sam = NULL;
    if(sam){
        core->fp_sam = fopen(sam,"w");
        if(core->fp_sam==NULL){
            fprintf(stderr, "Error opening file %s\n",sam);
            exit(EXIT_FAILURE);
        }
        sam_hdr_wr(core->fp_sam, core->ref);
    }

    core->n_samples = 0;
    core->total_reads = 0;

    return core;
}

void free_core(core_t *core){

    for(int i=0; i<core->opt.num_thread; i++){
        free_nrng(core->rand_time[i]);
        free_grng(core->rand_rlen[i]);
        free_nrng(core->rand_offset[i]);
        free_nrng(core->rand_median_before[i]);
        for (uint32_t j = 0; j < core->num_kmer; j++){
            free_nrng(core->kmer_gen[i][j]);
        }
        free(core->kmer_gen[i]);
    }

    free(core->kmer_gen);
    free(core->rand_time);
    free(core->rand_rlen);
    free(core->rand_offset);
    free(core->rand_median_before);
    free(core->ref_pos);
    free(core->rand_strand);
    if(core->opt.meth_freq){
        free(core->rand_meth);
    }

    free(core->model);
    free(core->cpgmodel);

    if(core->fp_fasta){
        fclose(core->fp_fasta);
    }
    if(core->fp_paf){
        fclose(core->fp_paf);
    }
    if(core->fp_sam){
        fclose(core->fp_sam);
    }
    slow5_close(core->sp);
    free_ref_sim(core->ref);

    free(core);
}

/* initialise a data batch */
db_t* init_db(core_t* core, int32_t n_rec) {
    db_t* db = (db_t*)(malloc(sizeof(db_t)));
    MALLOC_CHK(db);

    db->capacity_rec = core->opt.batch_size;
    db->n_rec = n_rec;

    db->mem_records = (char**)(calloc(db->capacity_rec,sizeof(void*)));
    MALLOC_CHK(db->mem_records);
    db->mem_bytes = (size_t*)(calloc(db->capacity_rec,sizeof(size_t)));
    MALLOC_CHK(db->mem_bytes);
    if(core->fp_fasta){
        db->fasta = (char**)(calloc(db->capacity_rec,sizeof(char*)));
        MALLOC_CHK(db->fasta);
    } else {
        db->fasta = NULL;
    }

    if(core->fp_paf){
        db->paf = (char**)(calloc(db->capacity_rec,sizeof(char*)));
        MALLOC_CHK(db->paf);
    } else {
        db->paf = NULL;
    }
    if(core->fp_sam){
        db->sam = (char**)(calloc(db->capacity_rec,sizeof(char*)));
        MALLOC_CHK(db->sam);
    } else {
        db->sam = NULL;
    }

    return db;
}

/* completely free a data batch */
void free_db(db_t* db) {

    int32_t i = 0;
    for (i = 0; i < db->n_rec; ++i) {
        free(db->mem_records[i]);
        if(db->fasta){
            free(db->fasta[i]);
        }
        if(db->paf){
            free(db->paf[i]);
        }
        if(db->sam){
            free(db->sam[i]);
        }
    }

    free(db->mem_records);
    free(db->mem_bytes);
    if(db->fasta){
        free(db->fasta);
    }
    if(db->paf){
        free(db->paf);
    }
    if(db->sam){
        free(db->sam);
    }

    free(db);
}

void fake_uuid(char *read_id, int64_t num){
    if(num>999999999999){
        ERROR("Too many reads that my lazy fake uuid generator cannot handle %ld is too large. Open an issue and I will fix this",num);
        exit(EXIT_FAILURE);
    }
    sprintf(read_id,"00000000-0000-0000-0000-%012d",(int)num);
}


//void gen_prefix_dna(int16_t *raw_signal, int64_t* n, int64_t *c, profile_t *profile, double offset);

char *gen_read(core_t *core, char **ref_id, int32_t *ref_len, int32_t *ref_pos, int32_t *rlen, char *c, int8_t rna, int tid);
int16_t *gen_sig(core_t *core, const char *read, int32_t len, double *offset, double *median_before, int64_t *len_raw_signal, int8_t rna, int tid, aln_t *aln);


/* process the ith read in the batch db */
void work_per_single_read(core_t* core,db_t* db, int32_t i, int tid) {

    opt_t opt = core->opt;
    ref_t *ref = core->ref;
    slow5_file_t *sp = core->sp;

    int8_t rna = opt.flag & SQ_RNA ? 1 : 0;

    double median_before = 0;
    double offset = 0;
    int64_t len_raw_signal =0;
    char *rid = NULL;
    char *seq = NULL;
    int32_t rlen = 0;
    char strand = '+';
    int32_t ref_pos_st = 0;
    int32_t ref_pos_end = 0;
    int32_t ref_len = 0;

    slow5_rec_t *slow5_record = slow5_rec_init();

    if(slow5_record == NULL){
        fprintf(stderr,"Could not allocate space for a slow5 record.");
        exit(EXIT_FAILURE);
    }

    aln_t *aln = NULL;
    if(core->fp_paf || core->fp_sam) aln = init_aln();

    if(opt.flag & SQ_FULL_CONTIG){
        rid = ref->ref_names[core->total_reads+i];
        rlen = ref->ref_lengths[core->total_reads+i];
        seq = ref->ref_seq[core->total_reads+i];
        strand = '+';
        ref_pos_st = 0;
        ref_pos_end = rlen;
    } else {
        seq=gen_read(core, &rid, &ref_len, &ref_pos_st, &rlen, &strand, rna, tid);
        ref_pos_end = ref_pos_st+rlen;
    }
    if(rlen*core->profile.dwell_mean >= UINT32_MAX ){
        WARNING("Read %s:%d-%d length*dwell_mean is too large: %ld. Double check parameters. May go out of memory.",rid,ref_pos_st,ref_pos_end,(int64_t)(rlen*core->profile.dwell_mean));
    }
    int16_t *raw_signal=gen_sig(core, seq, rlen, &offset, &median_before, &len_raw_signal, rna, tid, aln);
    assert(raw_signal != NULL && len_raw_signal > 0);
    if(len_raw_signal >= UINT32_MAX){
        ERROR("Read %s:%d-%d has too many samples: %ld. Double check parameters.",rid,ref_pos_st,ref_pos_end,len_raw_signal);
        exit(EXIT_FAILURE);
    }

    char *read_id= (char *)malloc(sizeof(char)*(10000));
    MALLOC_CHK(read_id);
    if(opt.flag & SQ_ONT){
        fake_uuid(read_id, core->total_reads+i+1);
    } else {
        sprintf(read_id,"S1_%ld!%s!%d!%d!%c",core->total_reads+i+1, rid, ref_pos_st, ref_pos_end, strand);
    }
    if(core->fp_fasta){
        db->fasta[i] = (char *)malloc(sizeof(char)*(strlen(read_id)+strlen(seq)+10)); //+10 bit inefficent - for now
        MALLOC_CHK(db->fasta[i]);
        sprintf(db->fasta[i],">%s\n%s\n",read_id,seq);
    }
    if(core->fp_paf || core->fp_sam){
        int64_t n_kmer = rlen - core->kmer_size+1;
        assert(n_kmer > 0);
        aln->read_id = read_id;
        aln->len_raw_signal = len_raw_signal;

        aln->strand = strand;
        aln->si_st_ref = rna ? ref_pos_end - core->kmer_size+1 : ref_pos_st;
        aln->si_end_ref = rna ? ref_pos_st : ref_pos_end - core->kmer_size+1;
        if(core->opt.flag & SQ_PAF_REF){
            aln->tid = rid;
            aln->tlen = !(opt.flag & SQ_FULL_CONTIG) ?  ref_len - core->kmer_size+1: n_kmer;
            aln->t_st = aln->si_st_ref;
            aln->t_end = aln->si_end_ref;
        } else {
            aln->tid = read_id;
            aln->tlen = n_kmer;
            aln->t_st = rna ? n_kmer : 0;
            aln->t_end = rna ? 0 : n_kmer;
        }
        if (core->fp_paf) db->paf[i] = paf_str(aln);
        if (core->fp_sam) db->sam[i] = sam_str(aln,seq,rid,ref_pos_st);

        free_aln(aln);
    }

    int64_t n_samples = __sync_fetch_and_add(&core->n_samples, len_raw_signal);
    set_record_primary_fields(&core->profile, slow5_record, read_id, offset, len_raw_signal, raw_signal);
    set_record_aux_fields(slow5_record, sp, median_before, core->total_reads+i, n_samples, opt.flag & SQ_ONT ? 1 : 0);

    //encode to a buffer
    if (slow5_encode(&db->mem_records[i], &db->mem_bytes[i], slow5_record, sp) < 0){
        fprintf(stderr,"Error encoding record\n");
        exit(EXIT_FAILURE);
    }

    if(!(opt.flag & SQ_FULL_CONTIG)){
        free(seq);
    }
    slow5_rec_free(slow5_record);


}

void work_db(core_t* core, db_t* db, void (*func)(core_t*,db_t*,int,int));

void process_db(core_t* core,db_t* db){
    double proc_start = realtime();
    work_db(core, db, work_per_single_read);
    double proc_end = realtime();
    core->process_db_time += (proc_end-proc_start);
}

/* write the output for a processed data batch */
void output_db(core_t* core, db_t* db) {

    double output_start = realtime();

    int32_t i = 0;
    for (i = 0; i < db->n_rec; i++) {

        //write the buffer
        if (slow5_write_bytes(db->mem_records[i], db->mem_bytes[i], core->sp) < 0){
            fprintf(stderr,"Error writing record!\n");
            exit(EXIT_FAILURE);
        }
        if(core->fp_fasta){
            fprintf(core->fp_fasta,"%s",db->fasta[i]);
        }
        if(core->fp_paf){
            fprintf(core->fp_paf,"%s",db->paf[i]);
        }
        if(core->fp_sam){
            fprintf(core->fp_sam,"%s",db->sam[i]);
        }
    }

    core->total_reads += db->n_rec;

    //core->read_index = core->read_index + db->n_rec;
    double output_end = realtime();
    core->output_time += (output_end-output_start);

}


///////////////////////////////////////////// SIMMAIN /////////////////////////////////////////////

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
    {"dwell-mean", required_argument, 0, 0 },      //12 dwell mean
    {"profile", required_argument, 0, 'm' },       //13 parameter profile
    {"kmer-model", required_argument, 0, 0},       //14 custom nucleotide k-mer model file
    {"prefix", required_argument, 0, 0},           //15 prefix such as adaptor
    {"dwell-std", required_argument, 0, 0 },       //16 dwell std
    {"threads", required_argument, 0, 't'},        //17 number of threads [8]
    {"batchsize", required_argument, 0, 'K'},      //18 batchsize - number of reads processed at once [1000]
    {"paf", required_argument, 0, 'c'},            //19 output paf file with alognments
    {"amp-noise", required_argument, 0, 0 },       //20 amplitude noise factor
    {"paf-ref", no_argument, 0, 0 },               //21 in paf, use ref as target
    {"sam", required_argument, 0, 'a'},            //22 output paf file with alognments
    {"coverage", required_argument, 0, 'f'},       //23 coverage
    {"digitisation", required_argument, 0, 0 },    //24 digitisation
    {"sample-rate", required_argument, 0, 0 },     //25 sample_rate
    {"range", required_argument, 0, 0 },           //26 range
    {"offset-mean", required_argument, 0, 0 },     //27 offset_mean
    {"offset-std", required_argument, 0, 0 },      //28 offset-std
    {"bps", required_argument, 0, 0 },             //29 bases per second
    {"median-before-mean", required_argument, 0, 0 },             //30 bases per second
    {"median-before-std", required_argument, 0, 0 },             //31 bases per second
    {"trans-count", required_argument, 0, 0 },             //32 transcript count
    {"trans-trunc", required_argument, 0, 0 },             //33 transcript truncate
    {"cdna", no_argument, 0, 0 },                   //34 cdna
    {"ont-friendly", required_argument, 0, 0},             //35 ont-friendly
    {"meth-freq", required_argument, 0, 0 },                   //36 meth-freq
    {"meth-model", required_argument, 0, 0 },                   //37 meth-model
    {0, 0, 0, 0}};


typedef struct {
    int8_t n;
    int8_t r;
    int8_t full_contigs;
    int8_t f;
    int8_t dwell_mean;
    int8_t bps;
    int8_t sample_rate;
    int8_t trans_count;
    int8_t trans_trunc;
    int8_t cdna;
} opt_gvn_t;

static void print_help(FILE *fp_help, opt_t opt, profile_t p, int64_t nreads) {

    fprintf(fp_help,"Usage: squigulator [OPTIONS] ref.fa -o out_signal.blow5\n");
    fprintf(fp_help,"\nbasic options:\n");
    fprintf(fp_help,"   -o FILE                    SLOW5/BLOW5 file to write\n");
    fprintf(fp_help,"   -x STR                     parameter profile (always applied before other options) [dna-r9-prom]\n");
    fprintf(fp_help,"                              e.g., dna-r9-min, dna-r9-prom, rna-r9-min, rna-r9-prom, dna-r10-min, dna-r10-prom, rna004-min, rna004-prom\n");
    fprintf(fp_help,"   -n INT                     number of reads to simulate [%ld]\n", nreads);
    fprintf(fp_help,"   -r INT                     mean read length (estimate only, unused for direct RNA) [%d]\n",opt.rlen);
    fprintf(fp_help,"   -f INT                     fold coverage to simulate (incompatible with -n)\n");
    fprintf(fp_help,"   -t INT                     number of threads [%d]\n",opt.num_thread);

    fprintf(fp_help,"   -h                         help\n");
    fprintf(fp_help,"   --ideal                    generate ideal signals with no noise\n");
    fprintf(fp_help,"   --version                  print version\n");
    fprintf(fp_help,"   --verbose INT              verbosity level [%d]\n",(int)get_log_level());
    fprintf(fp_help,"   --full-contigs             generate signals for complete contigs/sequences in the input (incompatible with -n, -r & -n)\n");

    fprintf(fp_help,"\nadvanced options:\n");
    fprintf(fp_help,"   -K INT                     batch size (max number of reads created at once) [%d]\n",opt.batch_size);
    fprintf(fp_help,"   -q FILE                    FASTA file to write simulated reads with no errors\n");
    fprintf(fp_help,"   -c FILE                    PAF file to write the alignment of simulated reads\n");
    fprintf(fp_help,"   -a FILE                    SAM file to write the alignment of simulated reads\n");
    fprintf(fp_help,"   --ideal-amp                generate signals with no amplitiude domain noise\n");
    fprintf(fp_help,"   --ideal-time               generate signals with no time domain noise\n");
    fprintf(fp_help,"   --amp-noise FLOAT          amplitude domain noise factor [%.1f]\n",opt.amp_noise);
    fprintf(fp_help,"   --dwell-mean FLOAT         mean of number of signal samples per base [%.1f]\n",p.dwell_mean);
    fprintf(fp_help,"   --dwell-std FLOAT          standard deviation of number of signal samples per base [%.1f]\n",p.dwell_std);
    fprintf(fp_help,"   --bps INT                  translocation speed in bases per second (incompatible with --dwell-mean) [%ld]\n",(long)(p.sample_rate/p.dwell_mean));
    fprintf(fp_help,"   --prefix=yes|no            generate prefixes such as adaptor (and polya for RNA) [no]\n");
    fprintf(fp_help,"   --seed INT                 seed or random generators (if 0, will be autogenerated) [%ld]\n",opt.seed);
    fprintf(fp_help,"   --paf-ref                  in paf output, use the reference as the target instead of read (needs -c)\n");
    fprintf(fp_help,"   --cdna                     generate cDNA reads (only for dna profiles & the reference must a transcriptome)\n");
    fprintf(fp_help,"   --trans-count FILE         simulate relative abundance using tsv file [transcript name, count]  (for direct-rna & cDNA)\n");
    fprintf(fp_help,"   --trans-trunc=yes/no       simulate transcript truncation (only for direct-rna) [no]\n");
    fprintf(fp_help,"   --ont-friendly=yes/no      generate fake uuid for readids and add a dummy end_reason [no]\n");
    fprintf(fp_help,"   --meth-freq FILE           simulate CpG methylation using frequency tsv file [chr, 0-based pos, freq] (for DNA)\n");

    fprintf(fp_help,"\ndeveloper options (less tested):\n");
    fprintf(fp_help,"   --digitisation FLOAT       ADC digitisation [%.1f]\n",p.digitisation);
    fprintf(fp_help,"   --sample-rate FLOAT        ADC sampling rate [%.1f]\n",p.sample_rate);
    fprintf(fp_help,"   --range FLOAT              ADC range [%.1f]\n",p.range);
    fprintf(fp_help,"   --offset-mean FLOAT        ADC offset mean [%.1f]\n",p.offset_mean);
    fprintf(fp_help,"   --offset-std FLOAT         ADC offset standard deviation [%.1f]\n",p.offset_std);
    fprintf(fp_help,"   --median-before-mean FLOAT Median before mean [%.1f]\n",p.median_before_mean);
    fprintf(fp_help,"   --median-before-std FLOAT  Median before standard deviation [%.1f]\n",p.median_before_std);
    fprintf(fp_help,"   --kmer-model FILE          custom nucleotide k-mer model file (format similar to f5c models)\n");
    fprintf(fp_help,"   --meth-model FILE          custom methylation k-mer model file (format similar to f5c models)\n");
    fprintf(fp_help,"\n");
    fprintf(fp_help,"See the manual page on GitHub for more details, options and the format of input/output files.\n");

    if(fp_help == stdout){
        exit(EXIT_SUCCESS);
    }
    exit(EXIT_FAILURE);

}

static inline void check_noneg_farg(float arg, char *arg_name){
    if(arg < 0){
        ERROR("%s should be non negative. You entered %.1f",arg_name, arg);
        exit(EXIT_FAILURE);
    }
}
static inline void check_pos_farg(float arg, char *arg_name){
    if(arg < 1){
        ERROR("%s should be larger than 0.0. You entered %.1f.",arg_name, arg);
        exit(EXIT_FAILURE);
    }
}

static inline void check_pos_iarg(int64_t arg, char *arg_name){
    if(arg < 1){
        ERROR("%s should be larger than 0. You entered %ld.",arg_name, arg);
        exit(EXIT_FAILURE);
    }
}

static void check_args(opt_gvn_t opt_gvn, int8_t rna, opt_t opt, char *paf) {
    //-n incompatible with --full-contigs
    if(opt_gvn.n && opt_gvn.full_contigs){
        WARNING("%s","Option -n is ignored when --full-contigs is set");
    }
    //-r incompatible with --full-contigs
    if(opt_gvn.r && opt_gvn.full_contigs){
        WARNING("%s","Option -r is ignored when --full-contigs is set");
    }
    //-f incompatible with --full-contigs
    if(opt_gvn.f && opt_gvn.full_contigs){
        WARNING("%s","Option -f is ignored when --full-contigs is set");
    }
    if (rna && opt_gvn.r){
        WARNING("%s","Option -r is ignored for RNA. Complete transcripts are simulated.");
    }
    if (opt.flag & SQ_PAF_REF && paf==NULL){
        WARNING("%s","Option --paf-ref is ineffective without -c.");
    }
    if (opt_gvn.f && opt_gvn.n){
        WARNING("%s","Option -n is ignored when -f is set");
    }
    if (opt_gvn.dwell_mean && opt_gvn.bps){
        WARNING("%s","Option --dwell-mean is ignored when --bps is provided");
    }
    if (opt_gvn.dwell_mean && opt_gvn.sample_rate){
        WARNING("%s","Option --dwell-mean is ignored when --sample-rate is provided");
    }
    if (opt_gvn.trans_count && !(rna || opt_gvn.cdna)){
        ERROR("%s","Option --trans-count requires an RNA profile or for a DNA profile --cdna be set.");
        exit(EXIT_FAILURE);
    }
    if (opt_gvn.trans_trunc && !(rna || opt_gvn.cdna)){
        ERROR("%s","Option --trans-trunc requires an RNA profile and for DNA --cdna be set.");
        exit(EXIT_FAILURE);
    }
    if (opt_gvn.cdna && rna){
        ERROR("%s","Option --cdna is only valid with DNA profiles.");
        exit(EXIT_FAILURE);
    }

}

//parse yes or no arguments : taken from minimap2
static inline void yes_or_no(opt_t* opt, uint64_t flag, int long_idx, const char* arg, int yes_to_set)
{
    if (yes_to_set) {
        if (strcmp(arg, "yes") == 0 || strcmp(arg, "y") == 0) {
            opt->flag |= flag;
        } else if (strcmp(arg, "no") == 0 || strcmp(arg, "n") == 0) {
            opt->flag &= ~flag;
        } else {
            WARNING("option '--%s' only accepts 'yes' or 'no'.",
                    long_options[long_idx].name);
        }
    } else {
        if (strcmp(arg, "yes") == 0 || strcmp(arg, "y") == 0) {
            opt->flag &= ~flag;
        } else if (strcmp(arg, "no") == 0 || strcmp(arg, "n") == 0) {
            opt->flag |= flag;
        } else {
            WARNING("option '--%s' only accepts 'yes' or 'no'.",
                    long_options[long_idx].name);
        }
    }
}

int sim_main(int argc, char* argv[], double realtime0) {

    const char* optstring = "o:hVn:q:r:x:v:K:t:c:a:f:";

    int longindex = 0;
    int32_t c = -1;

    FILE *fp_help = stderr;
    char *output_file = NULL;
    char *fasta = NULL;
    char *paf = NULL;
    char *sam = NULL;
    char *trans_count = NULL;

    opt_t opt;
    init_opt(&opt);

    profile_t p = prom_r9_dna_prof;

    int64_t nreads = 4000;
    int64_t coverage = -1;

    opt_gvn_t opt_gvn = {0};
    int8_t x_gvn = 0;

    //parse the user args
    while ((c = getopt_long(argc, argv, optstring, long_options, &longindex)) >= 0) {
        if (c=='V'){
            fprintf(stdout,"squigulator %s\n",SQ_VERSION);
            exit(EXIT_SUCCESS);
        }else if (c=='v'){
            int v = atoi(optarg);
            set_log_level((enum sq_log_level_opt)v);
        } else if (c=='h'){
            fp_help = stdout;
        } else if (c=='x'){
            p = set_profile(optarg, &opt);
            if(x_gvn){
                WARNING("%s","Providing -x multiple times may lead to unspecified behaviour.")
            }
            x_gvn = 1;
        } else if(c == 'o'){
            output_file=optarg;
        } else if (c == 0 && longindex == 4){   //generate ideal signal
            opt.flag |= SQ_IDEAL;
        } else if (c == 0 && longindex == 5){  //generate signal for complete contigs
            opt.flag |= SQ_FULL_CONTIG;
            opt_gvn.full_contigs = 1;
        } else if (c == 'n'){
            nreads = atoi(optarg);
            opt_gvn.n = 1;
            check_pos_iarg(nreads,"Number of reads");
        } else if (c == 'q'){
            fasta = optarg;
        } else if (c == 'r'){
            opt.rlen = atoi(optarg);
            opt_gvn.r = 1;
            if(opt.rlen < 200){
                WARNING("Read length %d is too short. Set to minimum length 200. Open an issue if you want short reads.",opt.rlen);
                opt.rlen = 200;
            }
        } else if (c == 'K') {
            opt.batch_size = atoi(optarg);
            check_pos_iarg(opt.batch_size, "Batch size");
        } else if (c == 't') {
            opt.num_thread = atoi(optarg);
            check_pos_iarg(opt.num_thread, "Number of threads");
        } else if (c == 'c'){
            paf = optarg;
        } else if (c == 'a'){
            sam = optarg;
        } else if (c == 'f'){
            opt_gvn.f = 1;
            coverage = atoi(optarg);
            check_pos_iarg(coverage,"Coverage");
        } else if (c == 0 && longindex == 9){  //seed
            opt.seed = atoi(optarg);
        } else if (c == 0 && longindex == 10){ //ideal-time
            opt.flag |= SQ_IDEAL_TIME;
        } else if (c == 0 && longindex == 11){ //ideal-amp
            opt.flag |= SQ_IDEAL_AMP;
        } else if (c == 0 && longindex == 12){ //dwell-mean
            opt_gvn.dwell_mean = 1;
            p.dwell_mean = atof(optarg);
            check_pos_farg(p.dwell_mean, "--dwell-mean");
        } else if (c == 0 && longindex == 14) { //custom nucleotide model file
            opt.model_file = optarg;
        } else if (c == 0 && longindex == 15) { //prefix
            yes_or_no(&opt, SQ_PREFIX, longindex, optarg, 1);
        } else if (c == 0 && longindex == 16) { //dwell-std
            p.dwell_std = atof(optarg);
            check_noneg_farg(p.dwell_std, "--dwell-std");
        } else if (c == 0 && longindex == 20) { //amp-noise
            opt.amp_noise = atof(optarg);
            check_noneg_farg(opt.amp_noise, "--amp-noise");
        } else if (c == 0 && longindex == 21) { //paf-ref
            opt.flag |= SQ_PAF_REF;
        } else if (c == 0 && longindex == 24) { //digitisation
            p.digitisation = atof(optarg);
            check_pos_farg(p.digitisation, "--digitisation"); //must in theory check if power of two
        } else if (c == 0 && longindex == 25) { //sample_rate
            opt_gvn.sample_rate = 1;
            p.sample_rate = atof(optarg);
            check_pos_farg(p.sample_rate, "--sample-rate");
        } else if (c == 0 && longindex == 26) { //range
            p.range = atof(optarg);
        } else if (c == 0 && longindex == 27) { //offset_mean
            p.offset_mean = atof(optarg);
        } else if (c == 0 && longindex == 28) { //offset_std
            p.offset_std = atof(optarg);
        } else if (c == 0 && longindex == 29) { //bases per second
            p.bps = atof(optarg);
            opt_gvn.bps = 1;
            check_pos_iarg(p.bps,"Bases per second");
        } else if (c == 0 && longindex == 30) { //median before mean
            p.median_before_mean = atof(optarg);
        } else if (c == 0 && longindex == 31) { //median before std
            p.median_before_std = atof(optarg);
        } else if (c == 0 && longindex == 32){ //transcript count
            trans_count = optarg;
            opt_gvn.trans_count = 1;
        } else if (c == 0 && longindex == 33){ //transcript trunctate
            opt_gvn.trans_trunc = 1;
            yes_or_no(&opt, SQ_TRANS_TRUNC, longindex, optarg, 1);
            WARNING("%s","Option --trans-trunc is experimental. Please report any issues.")
        } else if (c == 0 && longindex == 34){ //cdna
            opt_gvn.cdna = 1;
            opt.flag |= SQ_CDNA;
        } else if (c == 0 && longindex == 35){ //ont-friendly
            yes_or_no(&opt, SQ_ONT, longindex, optarg, 1);
        } else if (c == 0 && longindex == 36){ //meth freq
            opt.meth_freq = optarg;
            WARNING("%s","Option --meth-freq is experimental. Please report any issues.")
        } else if (c == 0 && longindex == 37){ //meth model
            opt.meth_model_file = optarg;
            //WARNING("%s","Option --meth-model is experimental. Please report any issues.")
        } else if (c == '?'){
            exit(EXIT_FAILURE);
        } else {
            exit(EXIT_FAILURE);
        }
    }

    if (argc-optind<1 || argc-optind> 1 || output_file==NULL ||  fp_help == stdout) {
        print_help(fp_help, opt, p, nreads);
    }

    int8_t rna = opt.flag & SQ_RNA ? 1 : 0;

    //check args
    check_args(opt_gvn, rna, opt, paf);

    if((opt.flag & SQ_CDNA) && (opt.flag & SQ_TRANS_TRUNC)){
        ERROR("%s","Option --trans-trunc is not yet implemnted for --cdna.");
        exit(EXIT_FAILURE);
    }

    if((opt.flag & SQ_CDNA) && (opt.flag & SQ_PREFIX)){
        ERROR("%s","Option --prefix is not yet implemnted for --cdna.");
        exit(EXIT_FAILURE);
    }

    if (opt.seed == 0){
        opt.seed = realtime0;
        VERBOSE("Using random seed: %ld",opt.seed);
    }

    char *refname = argv[optind];

    int64_t n = nreads;

    if(opt_gvn.dwell_mean){
        p.bps = round(p.sample_rate / p.dwell_mean);
        VERBOSE("dwell_mean=%.1f,sample_rate=%.1f,set bps=sample_rate/dwell_mean=%.1f",p.dwell_mean, p.sample_rate, p.bps);
    }

    if(opt_gvn.sample_rate || opt_gvn.bps){
        p.dwell_mean = round(p.sample_rate / p.bps);
        VERBOSE("bps=%.1f,sample_rate=%.1f,set dwell_mean=sample_rate/bps=%.1f",p.bps,p.sample_rate,p.dwell_mean);
    } else {
        assert(round(p.sample_rate / p.bps) == (int)(p.dwell_mean));
    }

    core_t *core = init_core(opt, p, refname, output_file, fasta, paf, sam, trans_count);

    if(opt.flag & SQ_FULL_CONTIG){
        n = core->ref->num_ref;
    } else {
        if (coverage > 0) {
            if(rna){
                n = (int64_t) (core->ref->num_ref * coverage);
                VERBOSE("rna mode: %ld reads will be simulated for each transcript. ntranscripts=%d, total reads %ld",coverage,core->ref->num_ref,n);
            } else {
                n = (int64_t) (core->ref->sum * coverage / opt.rlen);
                VERBOSE("dna mode: %ld reads of mean length %d will be simulated to achive requested %ldX coverage for the %ld base genome",n,opt.rlen,coverage,core->ref->sum);
            }
        }
    }

    print_model_stat(core->profile);

    int64_t done = 0;
    while(done < n){
        int64_t batch = n-done < opt.batch_size  ? n-done : opt.batch_size;
        assert(batch > 0);
        db_t *db = init_db(core, batch);
        process_db(core, db);
        output_db(core, db);
        free_db(db);
        done += batch;
        VERBOSE("%ld/%ld reads done",done,n);
    }
    assert(done == n);

    free_core(core);

    return 0;
}
