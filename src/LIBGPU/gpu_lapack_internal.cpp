// Written by Vasily Volkov.
// Copyright (c) 2009, The Regents of the University of California.
// All rights reserved.

#include "gpu_lapack_internal.h"
#include <cuda.h>

const unsigned int MB = 1<<20;

static unsigned int reserved, allocated;
static char *pool = NULL;
float *cpu_buffer = NULL;

//
//  Release the allocated GPU resources
//
void gpu_lapack_release( )
{
    if( pool != NULL )
    {
    	Q( cudaFree( pool ) );
        pool = NULL;
    }
    if( cpu_buffer != NULL )
    {
    	free( cpu_buffer );
        cpu_buffer = NULL;
    }
    reserved = allocated = 0;
}

//
//  Allocate all resources that we might ever need
//
bool gpu_lapack_init( )
{
    gpu_lapack_release( );
    
    cpu_buffer = (float*)malloc( nb*nb*sizeof(float) );
    if( cpu_buffer == NULL )
        return false;
    
    if( cublasInit( ) != CUBLAS_STATUS_SUCCESS )
    {
        free( cpu_buffer );
        return false;
    }
    
    //
    //  pre-allocate as much GPU memory as possible
    //
    unsigned int total;
    Q( cuMemGetInfo( &reserved, &total ) );
    while( cudaMalloc( (void**)&pool, reserved ) != cudaSuccess )
    {
        reserved -= MB;
        if( reserved < MB )
        {
            free( cpu_buffer );
            Q( cublasShutdown( ) );
            return false;
        }
    }
    
    //
    //  reset the error states
    //
    cudaGetLastError( );
    
    return true;
}

//
//  Reset the partitioning of the GPU memory
//
void gpu_malloc_reset( )
{
    allocated = 0;
}

//
//  Get a piece of the allocated GPU memory
//
p2_t gpu_malloc2D( int m, int n )
{
    if( m <= 0 || n <= 0 )
        return p2_t( NULL, 0 );
    
    //
    //  pad to align columns and avoid bank contention
    //
    unsigned int m_padded = (m+63)&~63|64;
    
    //
    //  pad to avoid page failures in custom BLAS3 kernels
    //
    unsigned int n_padded = (n+31)&~31;
    
    //
    //  allocate memory
    //
    p2_t p2( NULL, m_padded );
    unsigned int size = sizeof(float) * m_padded * n_padded;
    if( allocated + size <= reserved )
    {
        p2.A = (float*)(pool + allocated);
        allocated = (allocated + size + 63) & ~63;
    }
    return p2;
}
