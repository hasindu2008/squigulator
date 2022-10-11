/* @file ref.h
**
******************************************************************************/

#ifndef SQ_REF_H
#define SQ_REF_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>


//todo : can make more efficient using bit encoding
static inline uint32_t get_rank(char base) {
    if (base == 'A' || base == 'a' || base == 'R' || base == 'W' || base == 'M' || base == 'D' || base == 'H' || base == 'V') {
        return 0;
    } else if (base == 'C' || base == 'c' || base == 'Y' || base == 'B') {
        return 1;
    } else if (base == 'G' || base == 'g' || base == 'S' || base == 'K') {
        return 2;
    } else if (base == 'T' || base == 't' || base == 'U') {
        return 3;
    } else {
        WARNING("A None ACGT base found : %c", base);
        return 0;
    }
}

// return the lexicographic rank of the kmer amongst all strings of
// length k for this alphabet
static inline uint32_t get_kmer_rank(const char* str, uint32_t k) {
    //uint32_t p = 1;
    uint32_t r = 0;

    // from last base to first
    for (uint32_t i = 0; i < k; ++i) {
        //r += rank(str[k - i - 1]) * p;
        //p *= size();
        r += get_rank(str[k - i - 1]) << (i << 1);
    }
    return r;
}

static inline char complement(char c){
    char r;
    switch (c){
        case 'A':
        case 'a':
            r='T';
            break;
        case 'C':
        case 'c':
            r='G';
            break;
        case 'G':
        case 'g':
            r='C';
            break;
        case 'T':
        case 't':
            r='A';
            break;
        default:
            r='T';
            break;
    }
    return r;
}

static inline char *reverse_complement(char *f){
    char *r = (char *)malloc(strlen(f) + 1);
    for(unsigned int i=0; i<strlen(f); i++){
        r[i] = complement(f[strlen(f)-i-1]);
    }
    //TODO: new line forgotten?
    return r;
}


#endif
