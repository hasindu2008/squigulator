
/* @file  rand.h
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

#ifndef SQ_RAND_H
#define SQ_RAND_H

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include "error.h"


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


static inline nrng_t* init_nrng(int64_t seed,double mean, double std){
    nrng_t *rng = (nrng_t *)malloc(sizeof(nrng_t));
    MALLOC_CHK(rng);
    rng->m = mean;
    rng->s = std;
    rng->x = seed;
    return rng;
}

static inline grng_t* init_grng(int64_t seed,double alpha, double beta){
    grng_t *rng = (grng_t *)malloc(sizeof(grng_t));
    MALLOC_CHK(rng);
    rng->a = alpha;
    rng->b = beta;
    rng->x = seed;
    return rng;
}

static inline void free_nrng(nrng_t *rng){
    free(rng);
}

static inline void free_grng(grng_t *rng){
    free(rng);
}

static inline double rng(int64_t *xp){
    int64_t x = *xp;
    int64_t x_new = (16807 * (x % 127773)) - (2836 * (x / 127773));
    if (x_new > 0) x = x_new; else x = x_new + 2147483647;
    *xp = x_new;
    return (double)x/2147483647;
}

static inline double nrng(nrng_t *r){
    double u = 0.0;
    double t = 0.0;
    while (u == 0.0) u = rng(&r->x);
    while (t == 0.0) t = 2.0 * 3.14159265 * rng(&r->x);
    double x = sqrt(-2.0 * log(u)) * cos(t);
    return ((x * r->s) + r->m);
}

static inline double grng(grng_t *r){
    double s = 0.0;
    for(int i=0; i<r->a; i++){
        s += -log(1-rng(&r->x));
    }
    return s*(r->b);
}

#endif