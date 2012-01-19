// Written by Vasily Volkov.
// Copyright (c) 2009, The Regents of the University of California.
// All rights reserved.

#include <time.h>
#include <math.h>

//
//  generates uniformly distributed entries in [-1,1]
//
class random_matrix_t
{
protected:
    int n;
    unsigned int seed;
    float *a;
public:
    random_matrix_t( int _n ) : n(_n)
    {
        seed = (unsigned int)time( NULL );
        a = new float[n];
    }
    ~random_matrix_t( ) { delete a; }
    virtual float copy( float *A, int lda );
    virtual void dgemv( char trans, double *x, double *y );
};

//
//  same as above + practical trick to enforce s.p.d.
//
class random_spd_matrix_t : public random_matrix_t
{
    void enforce_spd( float *x ) { x[0] += 2*sqrt((float)n); }
public:
    random_spd_matrix_t( int n ) : random_matrix_t( n ) {}
    float copy( float *A, int lda );
    void dgemv( char trans, double *x, double *y );
};
