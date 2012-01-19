// Written by Vasily Volkov.
// Copyright (c) 2009, The Regents of the University of California.
// All rights reserved.

//
//  helper classes
//
class factors_t
{
public:
    int n;
    float *A;
    int lda;
    factors_t( int _n, float *_A, int _lda ) : n(_n), A(_A), lda(_lda) { }
    virtual void dgemv( char trans, double *x ) = 0;
};

class LU_factors_t : public factors_t
{
public:
    int *ipiv;
    LU_factors_t( int _n, float *_A, int _lda, int *_ipiv ) : factors_t( _n, _A, _lda ), ipiv(_ipiv) { }
    void dgemv( char trans, double *x );
};

class Cholesky_factors_t : public factors_t
{
public:
    Cholesky_factors_t( int _n, float *_A, int _lda ) : factors_t( _n, _A, _lda ) { }
    void dgemv( char trans, double *x );
};

class QR_factors_t : public factors_t
{
public:
    float *tau;
    QR_factors_t( int _n, float *_A, int _lda, float *_tau ) : factors_t( _n, _A, _lda ), tau(_tau) { }
    void dgemv( char trans, double *x );
};

class random_matrix_t;

//
//  the exported function
//
float estimate_residual( factors_t &factors, random_matrix_t &A );
