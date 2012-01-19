// Written by Vasily Volkov.
// Copyright (c) 2009, The Regents of the University of California.
// All rights reserved.

#include "random_matrix_t.h"

//
//  random number generator helper
//
#include <mkl.h>
class random_t
{
    VSLStreamStatePtr stream;
public:
    random_t( unsigned int seed )          { vslNewStream( &stream, VSL_BRNG_MCG31, seed ); }
    ~random_t( )                           { vslDeleteStream( &stream ); };
    void generate( float *dst, int count ) { vsRngUniform( VSL_METHOD_SUNIFORM_STD, stream, count, dst, -1.0f, 1.0f ); }
};

//
//  a random matrix
//
float random_matrix_t::copy( float *A, int lda )
{
    int ione = 1;
    float col, norm = 0;
    random_t random( seed );
    for( int i = 0; i < n; i++ )
    {
        random.generate( A+i*lda, n );
        col = sasum( &n, A+i*lda, &ione );
        norm = norm > col ? norm : col;
    }
    return norm;
}

void random_matrix_t::dgemv( char trans, double *x, double *y )
{
    random_t random( seed );
    for( int i = 0; i < n; i++ )
    { 
        random.generate( a, n );
        if( trans == 'N' || trans == 'n' )
        {
            for( int j = 0; j < n; j++ )
                y[j] -= a[j]*x[i];
        }
        else
        {
            double sum = 0;
            for( int j = 0; j < n; j++ )
                sum += a[j]*x[j];
            y[i] -= sum;
        }
    }
}

//
//  a random s.p.d. matrix
//
float random_spd_matrix_t::copy( float *A, int lda )
{
    int ione = 1;
    float col, norm = 0;
    random_t random( seed );
    for( int i = 0; i < n; i++ )
    {
        random.generate( A+i*lda+i, n-i );
        enforce_spd( A+i*lda+i );

        int tail = n-i;
        col = sasum( &i, A+i, &lda ) + sasum( &tail, A+i+i*lda, &ione );
        norm = norm > col ? norm : col;
    }
    return norm;
}

void random_spd_matrix_t::dgemv( char trans, double *x, double *y )
{
    random_t random( seed );
    for( int i = 0; i < n; i++ )
    { 
        random.generate( a+i, n-i );
        enforce_spd( a+i );

        for( int j = i; j < n; j++ )
            y[j] -= a[j]*x[i];

        double sum = 0;
        for( int j = i+1; j < n; j++ )
            sum += a[j]*x[j];
        y[i] -= sum;
    }
}
