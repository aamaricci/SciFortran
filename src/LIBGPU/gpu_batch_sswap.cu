// Written by Vasily Volkov.
// Copyright (c) 2009, The Regents of the University of California.
// All rights reserved.

#include "gpu_lapack_internal.h"

#define SWAPS_PER_RUN 64
#define VL 64

typedef struct
{   
    int nswaps, n;
	float *A;
	int lda;
	short ipiv[SWAPS_PER_RUN];
} run_t;

//
//  Do many BLAS's sswap-s in one kernel
//
static __global__ void batch_sswap( run_t run )
{   
	unsigned int tid = threadIdx.x + VL * blockIdx.x;
	if( tid >= run.n )
        return;
    
	float *A = run.A + tid;
	for( int i = 0; i < run.nswaps; i++ )
	{
		int j = run.ipiv[i];
		float temp   = A[i*run.lda];
		A[i*run.lda] = A[j*run.lda];
		A[j*run.lda] = temp;
	}
}	

extern "C" void gpu_batch_sswap( int batchsize, int n, int *ipiv, p2_t matrix )
{	
	if( n <= 0 || batchsize <= 0 )
		return;
    
	for( int i = 0; i < batchsize; i += SWAPS_PER_RUN )
	{
        int nb = imin( SWAPS_PER_RUN, batchsize - i );
        run_t run = { nb, n, matrix(0,i).A, matrix.lda };
		for( int j = 0; j < nb; j++ )
			run.ipiv[j] = ipiv[i+j] - i - 1;
        
        batch_sswap<<< (n+VL-1)/VL, VL >>>( run );
        Q( cudaGetLastError( ) );
	}
}	
