/* @file  ref.c
** Stupidly simple signal simulator
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

#include <zlib.h>
#include "ref.h"
#include "error.h"
#include "str.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

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

    kseq_destroy(seq);
    gzclose(fp);

    VERBOSE("Loaded %d reference sequences with total length %f Mbases",ref->num_ref,ref->sum/1000000.0);

    return ref;

}

void free_ref_sim(ref_t *ref){

    for(int i=0;i<ref->num_ref;i++){
        free(ref->ref_names[i]);
        free(ref->ref_seq[i]);
    }

    free(ref->ref_lengths);
    free(ref->ref_names);
    free(ref->ref_seq);
    free(ref);
}