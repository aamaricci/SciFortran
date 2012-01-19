// Written by Vasily Volkov.
// Copyright (c) 2009, The Regents of the University of California.
// All rights reserved.

#include "gpu_lapack_internal.h"

#define BLOCK_SIZE 8 // hand-tuned parameter

static __global__ void transpose_device( float *dst, int ldd, float *src, int lds )
{ 
	src += blockIdx.x*32 + threadIdx.x + ( blockIdx.y*32 + threadIdx.y ) * lds;
	dst += blockIdx.y*32 + threadIdx.x + ( blockIdx.x*32 + threadIdx.y ) * ldd;

	__shared__ float a[32][33];

    //
    //  load 32x32 block
    //
	for( int i = 0; i < 32; i += BLOCK_SIZE )
		 a[i+threadIdx.y][threadIdx.x] = src[i*lds];

	__syncthreads();

    //
    //  store transposed block
    //
	for( int i = 0; i < 32; i += BLOCK_SIZE )
        dst[i*ldd] = a[threadIdx.x][i+threadIdx.y];
}

static __global__ void transpose_inplace_device( float *matrix, int lda, int half, int parity )
{			
    bool bottom = blockIdx.x + parity > blockIdx.y;
	int ibx = bottom ? (blockIdx.x + parity - 1) : (blockIdx.y + half - parity);
	int iby = bottom ? blockIdx.y       : (blockIdx.x + half);
	
	ibx *= 32;
	iby *= 32;
	int inx = threadIdx.x;
	int iny = threadIdx.y;
	
	__shared__ float a[32][33], b[32][33];
    
    //
    //  load 32x32 block
    //
	float *A = matrix + ibx + inx + ( iby + iny ) * lda;
	for( int i = 0; i < 32; i += BLOCK_SIZE )
		 a[i+threadIdx.y][threadIdx.x] = A[i*lda];
	
	if( ibx == iby )
	{
        //
        //  this is a diagonal block
        //
		__syncthreads();
        
        //
        //  store transposed block
        //
        for( int i = 0; i < 32; i += BLOCK_SIZE )
	    	 A[i*lda] = a[threadIdx.x][i+threadIdx.y];
	}
	else
	{
        //
        //  load the opposite 32x32 block
        //
		float *B = matrix + iby + inx + ( ibx + iny ) * lda;
    	for( int i = 0; i < 32; i += BLOCK_SIZE )
	    	 b[i+threadIdx.y][threadIdx.x] = B[i*lda];
        
        __syncthreads();
        
        //
        //  store transposed blocks in reverse order
        //
        for( int i = 0; i < 32; i += BLOCK_SIZE )
	    	 A[i*lda] = b[threadIdx.x][i+threadIdx.y];
    	for( int i = 0; i < 32; i += BLOCK_SIZE )
	    	 B[i*lda] = a[threadIdx.x][i+threadIdx.y];
	}
}	

extern "C" void gpu_transpose( int m, int n, p2_t dst, p2_t src )
{
    if( m <= 0 || n <= 0 )
        return;

	dim3 threads( 32, BLOCK_SIZE, 1 );
	dim3 grid( (m+31)/32, (n+31)/32, 1 );
	transpose_device<<< grid, threads >>>( dst.A, dst.lda, src.A, src.lda );
    Q( cudaGetLastError( ) );
}
	
extern "C" void gpu_transpose_inplace( int n, p2_t matrix )
{
    if( n <= 0 )
        return;

	int in = (n+31) / 32;
	dim3 threads( 32, BLOCK_SIZE );
	dim3 grid( in|1, in/2+(in&1) );
	transpose_inplace_device<<< grid, threads >>>( matrix.A, matrix.lda, grid.y, in&1 );
    Q( cudaGetLastError( ) );
}
