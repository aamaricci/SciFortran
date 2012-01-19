// Written by Vasily Volkov.
// Copyright (c) 2009, The Regents of the University of California.
// All rights reserved.

#include <stdio.h>
#include <string.h>
#include <mkl.h>
#include <cuda_runtime.h>
#include "gpu_lapack.h"
#include "random_matrix_t.h"
#include "residual_norm.h"

const float eps = 1.1920929e-7f;

//
//  timer
//
#ifdef _WIN32
#include <windows.h>
double read_timer( )
{	
    static bool initialized = false;
    static double dfreq;
    static LARGE_INTEGER seconds0;
    LARGE_INTEGER temp;
    if( !initialized )
    {
        QueryPerformanceFrequency (&temp);
        dfreq = 1.0/temp.QuadPart;
        QueryPerformanceCounter (&seconds0);
        initialized = true;
    }
    QueryPerformanceCounter (&temp);
    return (temp.QuadPart - seconds0.QuadPart)*dfreq;
}	
#else
#include <sys/time.h>
double read_timer( )
{
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if( !initialized )
    {
        gettimeofday( &start, NULL );
        initialized = true;
    }

    gettimeofday( &end, NULL );

    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}
#endif

//
//  This code provides examples of using the library
//
int main(int argc, char** argv)
{   
    //
    //  process command line options
    //
    int idevice = 0;
    for( int i = 1; i < argc-1 ; i ++ )
        if( strcmp( argv[i], "-device" ) == 0 )
            idevice = atoi( argv[i+1] );
    
    bool use_gpu = true;
    for( int i = 1; i < argc; i ++ )
        if( strcmp( argv[i], "-cpu" ) == 0 )
            use_gpu = false;
    
    bool pinned = use_gpu;
    
    //
    //  initialize
    //
    if( cudaSetDevice( idevice ) != cudaSuccess )
    {
        printf( "Failed to set device %d\n", idevice );
        return 1;
    }
    
    struct cudaDeviceProp prop;
    cudaGetDeviceProperties( &prop, idevice );
    printf( "\nDevice: %s, %.0f MHz clock, %.0f MB memory.\n\n", prop.name, prop.clockRate/1000.f, prop.totalGlobalMem/1024.f/1024.f );
    
    if( !gpu_lapack_init( ) )
    {
        printf( "failed to initialize the library\n" );
        return 1;
    }
    
    //
    // solve a set of problems
    //
    printf( "Errors reported are 1-norms of the residual such as ||A-QR||_1.\n\n" );
    printf( "            Cholesky            LU               QR      \n" );
    printf( "         --------------   --------------   --------------\n" );
    printf( "   N     Gflop/s  error   Gflop/s  error   Gflop/s  error\n" );
    printf( " -----   --------------   --------------   --------------\n" );
    float *A = NULL;
    for( int n = 1000; ; n += 1000 )
    {
        int lda = ((n+15)&~15);

        //
        //  allocate matrix in pinned memory
        //
        if( pinned )
        {
            cudaError_t ret = cudaMallocHost( (void**)&A, n*lda*sizeof(float) );
            if( ret != cudaSuccess )
                break;
        }
        else
        {
            A = (float*)malloc( n*lda*sizeof(float) );
            if( A == NULL )
                break;
        }
        
        //
        //  Common variables
        //
        random_matrix_t matrix( n );
        random_spd_matrix_t spd_matrix( n );
        
        int info;
        double seconds;
        float residual_norm, matrix_norm;

        printf( " %5d", n );
        
        //
        //  Cholesky factorization
        //
        matrix_norm = spd_matrix.copy( A, lda );
        
        seconds = read_timer( );
        use_gpu ? gpu_spotrf(       n, A,  lda,  info ) 
                :     spotrf( "L", &n, A, &lda, &info );
        seconds = read_timer( ) - seconds;
        
        if( info != 0 )
            break;
        
        Cholesky_factors_t Cholesky_factors( n, A, lda );
        residual_norm = estimate_residual( Cholesky_factors, spd_matrix );
        printf( "   %7.2f %6.2f", 1e-9*n*n*n/3/seconds, residual_norm / matrix_norm / eps );
        
        //
        //  LU factorization
        //
        int *ipiv = (int*)malloc( n*sizeof(int) );
        matrix_norm = matrix.copy( A, lda );
        
        seconds = read_timer( );
        use_gpu ? gpu_sgetrf(  n,     A,  lda, ipiv,  info )
                :     sgetrf( &n, &n, A, &lda, ipiv, &info );
        seconds = read_timer( ) - seconds;
        
        if( info != 0 )
        {
            free( ipiv );
            break;
        }
        
        LU_factors_t LU_factors( n, A, lda, ipiv );
        residual_norm = estimate_residual( LU_factors, matrix );
        printf( "   %7.2f %6.2f", 2e-9*n*n*n/3/seconds, residual_norm / matrix_norm / eps );
        
        free( ipiv );
        
        //
        //  QR factorization
        //
        float requirement;
        int lwork = -1;
        use_gpu ? gpu_sgeqrf(  n,     A,  lda, A, &requirement,  lwork,  info )
                :     sgeqrf( &n, &n, A, &lda, A, &requirement, &lwork, &info );
        if( info != 0 )
            break;
        lwork = (int)requirement;
        float *work = (float*)malloc( lwork*sizeof(int) );
        float *tau = (float*)malloc( n*sizeof(float) );
        
        matrix_norm = matrix.copy( A, lda );
        
        seconds = read_timer( );
        use_gpu ? gpu_sgeqrf(  n,     A,  lda, tau, work,  lwork,  info ) 
                :     sgeqrf( &n, &n, A, &lda, tau, work, &lwork, &info );
        seconds = read_timer( ) - seconds;
        
        if( info < 0 )
        {
            free( work );
            free( tau );
            break;
        }
        
        QR_factors_t QR_factors( n, A, lda, tau );
        residual_norm = estimate_residual( QR_factors, matrix );
        printf( "   %7.2f %6.2f\n", 4e-9*n*n*n/3/seconds, residual_norm / matrix_norm / eps );
    	
        free( tau );
        free( work );
        
        if( pinned )
            cudaFreeHost( A );
        else
            free( A );
        
        A = NULL;
    }
    
    printf( "\n" );
    
    //
    //  release resources
    //
    if( A != NULL )
    {
        if( pinned )
            cudaFreeHost( A );
        else
            free( A );
    }
    
    gpu_lapack_release( );
    
    return 0;
}	
