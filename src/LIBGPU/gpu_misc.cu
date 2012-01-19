// Written by Vasily Volkov.
// Copyright (c) 2009, The Regents of the University of California.
// All rights reserved.

#include "gpu_lapack_internal.h"

//
//  make matrix explicitly lower triangular, unit diagonal
//
static __global__ void enforceLU_device( float *matrix, int lda )
{ 
    int i = threadIdx.x;
    int j = blockIdx.x;
    if( i <= j )
        matrix[i + j*lda] = (i == j) ? 1 : 0;
}

extern "C" void gpu_enforceLU( int n, p2_t matrix )
{
    if( n > 0 )
    {
    	enforceLU_device<<< n, n >>>( matrix.A, matrix.lda );
        Q( cudaGetLastError( ) );
    }
}
