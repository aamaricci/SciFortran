// Written by Vasily Volkov.
// Copyright (c) 2009, The Regents of the University of California.
// All rights reserved.

#include <string.h>
#include <stdlib.h>
#include <mkl_blas.h>
#include "random_matrix_t.h"
#include "residual_norm.h"

//
//  subset of BLAS/LAPACK routines that accept input in single precision
//  but compute in double precision
//

void permute( char trans, int n, int *ipiv, double *x )
{
    bool forward = trans == 'N' || trans == 'n';
    for( int j = 0; j < n-1; j++ )
    {
        int i = forward ? j : n-2-j;
        int ip = ipiv[i]-1;
        double t = x[i];
        x[i] = x[ip];
        x[ip] = t;
    }
}

void dtrmv_U( int n, float *A, int lda, double *x )
{
    for( int j = 0; j < n; j++ )
    {
        double temp = x[j];
        for( int i = 0; i < j; i++ )
            x[i] += A[i+j*lda] * temp;
        x[j] = A[j+j*lda] * temp;
    }
}

void dtrmv_L( char diag, int n, float *A, int lda, double *x )
{
    bool nounit = diag == 'N' || diag == 'n';
    for( int j = n-1; j >= 0; j-- )
    {
        double temp = x[j];
        for( int i = n-1; i > j; i-- )
            x[i] += A[i+j*lda] * temp;
        if( nounit )
            x[j] = A[j+j*lda] * temp;
    }
}

void dtrmv_LT( char diag, int n, float *A, int lda, double *x )
{
    bool nounit = diag == 'N' || diag == 'n';
    for( int j = 0; j < n; j++ )
    {
        double temp = x[j];
        temp = nounit ? A[j+j*lda] * temp : temp;
        for( int i = j+1; i < n; i++ )
            temp += A[i+j*lda] * x[i];
        x[j] = temp;
    }
}

void dtrmv_UT( int n, float *A, int lda, double *x )
{
    for( int j = n-1; j >= 0; j-- )
    {
        double temp = 0;
        for( int i = 0; i <= j ; i++ )
            temp += A[i+j*lda] * x[i];
        x[j] = temp;
    }
}

void dorm2r_L( char trans, int n, float *A, int lda, float *tau, double *x )
{
    bool transposed = trans == 'T' || trans == 't';
    int i0 = transposed ? 0 : n-1;
    int id = transposed ? 1 : -1;
    int in = transposed ? n : -1;
    for( int i = i0; i != in; i += id )
    {
        double a = x[i];
        for( int j = i+1; j < n; j++ )
            a += A[i*lda+j]*x[j];
        a *= tau[i];

        x[i] -= a;
        for( int j = i+1; j < n; j++ )
            x[j] -= a*A[i*lda+j];
    }
}

//
//  given a factorization, such as A = Q*R, and vector x, computes A*x
//
void Cholesky_factors_t::dgemv( char trans, double *x )
{
    dtrmv_LT( 'N', n, A, lda, x );
    dtrmv_L( 'N', n, A, lda, x );
}

void LU_factors_t::dgemv( char trans, double *x )
{
    if( trans == 'N' || trans == 'n' )
    {
        dtrmv_U( n, A, lda, x );
        dtrmv_L( 'U', n, A, lda, x );
        permute( 'T', n, ipiv, x );
    }
    else
    {
        permute( 'N', n, ipiv, x );
        dtrmv_LT( 'U', n, A, lda, x );
        dtrmv_UT( n, A, lda, x );
    }
}

void QR_factors_t::dgemv( char trans, double *x )
{
    if( trans == 'N' || trans == 'n' )
    {
        dtrmv_U( n, A, lda, x );
        dorm2r_L( 'N', n, A, lda, tau, x );
    }
    else
    {
        dorm2r_L( 'T', n, A, lda, tau, x );
        dtrmv_UT( n, A, lda, x );
    }
}

//
//  estimates 1-norm of the residual using Hager's norm estimator
//  see Ch.2.4.3 in Demmel, J. W. 1997. Applied Numerical Linear Algebra, SIAM.
//
float estimate_residual( factors_t &factors, random_matrix_t &matrix )
{
    int n = factors.n;

    double *x = (double*)malloc( 4*n*sizeof(double) );
    double *y  = x + 1*n;
    double *z  = x + 2*n;
    double *xi = x + 3*n;
    
    for( int i = 0; i < n; i++ )
        x[i] = 1.0 / n;
    
    int ione = 1;
    double fone = 1.0, fzero = 0.0;
    for( int k = 0; k < 5; k++ )
    {
        dcopy( &n, x, &ione, y, &ione );
        factors.dgemv( 'N', y );
        matrix.dgemv( 'N', x, y );
        
        for( int i = 0; i < n; i++ )
            xi[i] = y[i] > 0.0 ? 1.0 : -1.0;
        
        dcopy( &n, xi, &ione, z, &ione );
        factors.dgemv( 'T', z );
        matrix.dgemv( 'T', xi, z );
        
        int j = idamax( &n, z, &ione ) - 1;
        if( fabs(z[j]) <= ddot( &n, z, &ione, x, &ione ) )
            break;
        memset( x, 0, n*sizeof(double) );
        x[j] = 1.0;
    }
    
    float y1 = (float)dasum( &n, y, &ione );
    
    free( x );
    
    return y1;
}
