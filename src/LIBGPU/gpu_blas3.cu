// Written by Vasily Volkov.
// Copyright (c) 2009, The Regents of the University of California.
// All rights reserved.

#include "gpu_lapack_internal.h"
//#define USE_CUBLAS

//
//  have to unroll some of the loops manually
//
__device__ void rank1_update( float a, const float *b, float *c )
{
	c[0] += a*b[0];
	c[1] += a*b[1];
	c[2] += a*b[2];
	c[3] += a*b[3];
	c[4] += a*b[4];
	c[5] += a*b[5];
	c[6] += a*b[6];
	c[7] += a*b[7];
	c[8] += a*b[8];
	c[9] += a*b[9];
	c[10] += a*b[10];
	c[11] += a*b[11];
	c[12] += a*b[12];
	c[13] += a*b[13];
	c[14] += a*b[14];
	c[15] += a*b[15];
}

__device__ void rankk_update( int k, const float *A, int lda, const float *b, int ldb, float *c )
{
    if( k <= 0 ) return;

    int i = 0;
    rank1_update( A[0], &b[i*ldb], c ); if( ++i >= k ) return; A += lda;
    rank1_update( A[0], &b[i*ldb], c ); if( ++i >= k ) return; A += lda;
    rank1_update( A[0], &b[i*ldb], c ); if( ++i >= k ) return; A += lda;
    rank1_update( A[0], &b[i*ldb], c ); if( ++i >= k ) return; A += lda;
    
    rank1_update( A[0], &b[i*ldb], c ); if( ++i >= k ) return; A += lda;
    rank1_update( A[0], &b[i*ldb], c ); if( ++i >= k ) return; A += lda;
    rank1_update( A[0], &b[i*ldb], c ); if( ++i >= k ) return; A += lda;
    rank1_update( A[0], &b[i*ldb], c ); if( ++i >= k ) return; A += lda;
    
    rank1_update( A[0], &b[i*ldb], c ); if( ++i >= k ) return; A += lda;
    rank1_update( A[0], &b[i*ldb], c ); if( ++i >= k ) return; A += lda;
    rank1_update( A[0], &b[i*ldb], c ); if( ++i >= k ) return; A += lda;
    rank1_update( A[0], &b[i*ldb], c ); if( ++i >= k ) return; A += lda;
    
    rank1_update( A[0], &b[i*ldb], c ); if( ++i >= k ) return; A += lda;
    rank1_update( A[0], &b[i*ldb], c ); if( ++i >= k ) return; A += lda;
    rank1_update( A[0], &b[i*ldb], c );
}

__device__ void store_block( int num, float alpha, float *c, float beta, float *C, int ldc )
{
    if( num <= 0 ) return;

    if( beta == 0 )
    {
        //
        //  for the case when C is initialized with inf or NaN
        //
        int i = 0; 
        C[0] = alpha*c[i++]; if( i >= num ) return; C += ldc;  
        C[0] = alpha*c[i++]; if( i >= num ) return; C += ldc;  
        C[0] = alpha*c[i++]; if( i >= num ) return; C += ldc;  
        C[0] = alpha*c[i++]; if( i >= num ) return; C += ldc;  
        
        C[0] = alpha*c[i++]; if( i >= num ) return; C += ldc;  
        C[0] = alpha*c[i++]; if( i >= num ) return; C += ldc;  
        C[0] = alpha*c[i++]; if( i >= num ) return; C += ldc;  
        C[0] = alpha*c[i++]; if( i >= num ) return; C += ldc;  

        C[0] = alpha*c[i++]; if( i >= num ) return; C += ldc;  
        C[0] = alpha*c[i++]; if( i >= num ) return; C += ldc;  
        C[0] = alpha*c[i++]; if( i >= num ) return; C += ldc;  
        C[0] = alpha*c[i++]; if( i >= num ) return; C += ldc;  

        C[0] = alpha*c[i++]; if( i >= num ) return; C += ldc;  
        C[0] = alpha*c[i++]; if( i >= num ) return; C += ldc;  
        C[0] = alpha*c[i++]; if( i >= num ) return; C += ldc;  
        C[0] = alpha*c[i++];
    }
    else
    {
        int i = 0; 
        C[0] = alpha*c[i++] + beta*C[0]; if( i >= num ) return; C += ldc;  
        C[0] = alpha*c[i++] + beta*C[0]; if( i >= num ) return; C += ldc;  
        C[0] = alpha*c[i++] + beta*C[0]; if( i >= num ) return; C += ldc;  
        C[0] = alpha*c[i++] + beta*C[0]; if( i >= num ) return; C += ldc;  
        
        C[0] = alpha*c[i++] + beta*C[0]; if( i >= num ) return; C += ldc;  
        C[0] = alpha*c[i++] + beta*C[0]; if( i >= num ) return; C += ldc;  
        C[0] = alpha*c[i++] + beta*C[0]; if( i >= num ) return; C += ldc;  
        C[0] = alpha*c[i++] + beta*C[0]; if( i >= num ) return; C += ldc;  

        C[0] = alpha*c[i++] + beta*C[0]; if( i >= num ) return; C += ldc;  
        C[0] = alpha*c[i++] + beta*C[0]; if( i >= num ) return; C += ldc;  
        C[0] = alpha*c[i++] + beta*C[0]; if( i >= num ) return; C += ldc;  
        C[0] = alpha*c[i++] + beta*C[0]; if( i >= num ) return; C += ldc;  

        C[0] = alpha*c[i++] + beta*C[0]; if( i >= num ) return; C += ldc;  
        C[0] = alpha*c[i++] + beta*C[0]; if( i >= num ) return; C += ldc;  
        C[0] = alpha*c[i++] + beta*C[0]; if( i >= num ) return; C += ldc;  
        C[0] = alpha*c[i++] + beta*C[0];
    }
}

//
//  C = alpha*A*B + beta*C
//
static __global__ void sgemmNN_device( int m, int n, const float *A, int lda, const float *B, int ldb, float* C, int ldc, int k, float alpha, float beta )
{
	const int inx = threadIdx.x;
	const int iny = threadIdx.y;
	const int ibx = blockIdx.x * 64;
	const int iby = blockIdx.y * 16;
	const int row = ibx + inx + iny*16;
	
	A += row;
	B += inx + ( iby + iny ) * ldb;
	C += row  + iby * ldc;
	
	float c[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    
	__shared__ float b[16][17];
	for( ; k > 0; k -= 16 )
	{
#pragma unroll
		for( int i = 0; i < 16; i += 4 )
			b[inx][iny+i]  = B[i*ldb];
		__syncthreads();

        if( k < 16 )
            break;

#pragma unroll
	    for( int i = 0; i < 16; i++, A += lda )
		    rank1_update( A[0], &b[i][0], c ); 
	    __syncthreads();
		
		B += 16;
	};

    rankk_update( k, A, lda, &b[0][0], 17, c );

    if( row >= m ) 
        return;
    
    store_block( n - iby, alpha, c, beta, C, ldc);
}	

//
//  C = alpha*A*B^T + beta*C
//
static __global__ void sgemmNT_device( int m, int n, const float *A, int lda, const float *B, int ldb, float* C, int ldc, int k, float alpha, float beta )
{
	const int inx = threadIdx.x;
	const int iny = threadIdx.y;
	const int ibx = blockIdx.x * 64;
	const int iby = blockIdx.y * 16;
	const int row = ibx + inx + iny*16;
	
	A += row;
	B += iby + inx + iny * ldb;
	C += row  + iby * ldc;
	
	float c[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    
	__shared__ float b[16][16];
	for( ; k > 0; k -= 16 )
	{
#pragma unroll
		for( int i = 0; i < 16; i += 4 )
			b[iny+i][inx]  = B[i*ldb];
		__syncthreads();
        
		if( k < 16 ) 
            break;
        
#pragma unroll
	    for( int i = 0; i < 16; i++, A += lda )
		    rank1_update( A[0], &b[i][0], c ); 
        __syncthreads();

        B += 16*ldb;
    }
    
    rankk_update( k, A, lda, &b[0][0], 16, c );
    
    if( row >= m ) 
        return;
    
    store_block( n - iby, alpha, c, beta, C, ldc);
}	

//
//  C = alpha*A*A^T + beta*C for lower triangular C
//
extern "C" __global__ void ssyrkLN_device( const float *A, int lda, float* C, int ldc, int k, float alpha, float beta, int parity, int n )
{	
	const int by = blockIdx.y/4;
	const bool bottom = ( blockIdx.x + parity > by );
	int ibx = bottom ? blockIdx.x + parity - 1 : by + gridDim.y/4 - parity;
	int iby = bottom ? blockIdx.y              : blockIdx.x*4 + (blockIdx.y%4) + gridDim.y;
	ibx *= 64;
	iby *= 16;
	
	const int inx = threadIdx.x;
	const int iny = threadIdx.y;
	const int row  = ibx + inx + iny*16;
	
	const float *B = A + iby + inx + iny * lda;
	A += row;
	C += row  + iby * ldc;
	
	float c[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    
	__shared__ float b[16][16];
	for( ; k > 0; k -= 16 )
	{
#pragma unroll
		for( int i = 0; i < 16; i += 4 )
			b[iny+i][inx]  = B[i*lda];
		__syncthreads();

		if( k < 16 ) 
            break;

#pragma unroll
	    for( int i = 0; i < 16; i++, A += lda )
		    rank1_update( A[0], &b[i][0], c ); 
        __syncthreads();

        B += 16*lda;
    }
    
    rankk_update( k, A, lda, &b[0][0], 16, c );
    
    if( row >= n ) 
        return;
    
    store_block( row + 1 - iby, alpha, c, beta, C, ldc);
}

//
//  Symmetric rank k update
//  See http://www.netlib.org/blas/ssyrk.f
//
extern "C" void gpu_ssyrkLN( int n, int k, float alpha, const p2_t A, float beta, p2_t C )
{	
	if( n <= 0 || k <= 0 )
		return;
#ifdef USE_CUBLAS
    cublasSsyrk( 'L', 'N', n, k, alpha, A.A, A.lda, beta, C.A, C.lda );
    Q( cublasGetError( ) );
#else
	int in = (n+63) / 64;
	int parity = in&1;
	dim3 grid( in|1, 2*(in+parity) );
	dim3 threads( 16, 4 );
	ssyrkLN_device<<< grid, threads >>>( A.A, A.lda, C.A, C.lda, k, alpha, beta, parity, n );
    Q( cudaGetLastError( ) );
#endif
}	

//
//  Matrix-matrix multiplications
//  See http://www.netlib.org/blas/sgemm.f
//
extern "C" void gpu_sgemmNN( int m, int n, int k, float alpha, const p2_t A, const p2_t B, float beta, p2_t C )
{	
	if( m <= 0 || n <= 0 || k <= 0 )
		return;
#ifdef USE_CUBLAS
    cublasSgemm( 'N', 'N', m, n, k, alpha, A.A, A.lda, B.A, B.lda, beta, C.A, C.lda );
    Q( cublasGetError( ) );
#else
	dim3 grid( (m+63)/64, (n+15)/16 ), threads( 16, 4 );
	sgemmNN_device<<<grid, threads>>>( m, n, A.A, A.lda, B.A, B.lda, C.A, C.lda, k, alpha, beta );
    Q( cudaGetLastError( ) );
#endif
}	

extern "C" void gpu_sgemmNT( int m, int n, int k, float alpha, const p2_t A, const p2_t B, float beta, p2_t C )
{	
	if( m <= 0 || n <= 0 || k <= 0 )
		return;
#ifdef USE_CUBLAS
    cublasSgemm( 'N', 'T', m, n, k, alpha, A.A, A.lda, B.A, B.lda, beta, C.A, C.lda );
    Q( cublasGetError( ) );
#else
	dim3 grid( (m+63)/64, (n+15)/16 ), threads( 16, 4 );
	sgemmNT_device<<<grid, threads>>>( m, n, A.A, A.lda, B.A, B.lda, C.A, C.lda, k, alpha, beta );
    Q( cudaGetLastError( ) );
#endif
}	
