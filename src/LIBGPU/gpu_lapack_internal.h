// Written by Vasily Volkov.
// Copyright (c) 2009, The Regents of the University of California.
// All rights reserved.

#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <cublas.h>

inline int imin( int a, int b ) { return a < b ? a : b; }
#define Q( condition ) {if( (condition) != 0 ) { fprintf( stderr, "\n FAILURE in %s, line %d\n", __FILE__, __LINE__ );exit( 1 );}}

const int nb = 64;
extern float *cpu_buffer;

//
//  Pointer to matrix
//
struct p2_t
{
    float *A;
	int lda;

    p2_t() {}
	p2_t( float *_A, int _lda ) : A(_A),lda(_lda) {}
	p2_t operator() (int i, int j) { return p2_t(A+i+j*lda,lda); };
	float &at (int i, int j) { return A[i+j*lda]; };
};

//
//  Memory allocation
//
void gpu_malloc_reset( );
p2_t gpu_malloc2D( int m, int n );

//
//  Shortcuts to CUDA and CUBLAS calls
//
inline void upload( int m, int n, p2_t dst, p2_t src )
{	
	if( m > 0 && n > 0 )
		Q( cudaMemcpy2D( dst.A, dst.lda*sizeof(float), src.A, src.lda*sizeof(float), m*sizeof(float), n, cudaMemcpyHostToDevice ) );
}	

inline void uploadL( int n, p2_t dst, p2_t src )
{	
    const int nb = 128;
	for( int i = 0; i < n; i += nb )
        upload( n-i, imin(nb,n-i), dst(i,i), src(i,i) );
}	

inline void download( int m, int n, p2_t dst, p2_t src )
{	
	if( m > 0 && n > 0 )
		Q( cudaMemcpy2D( dst.A, dst.lda*sizeof(float), src.A, src.lda*sizeof(float), m*sizeof(float), n, cudaMemcpyDeviceToHost ) );
}	

inline void gpu_copy( int m, int n, p2_t dst, p2_t src )
{	
	if( m > 0 && n > 0 )
		Q( cudaMemcpy2D( dst.A, dst.lda*sizeof(float), src.A, src.lda*sizeof(float), m*sizeof(float), n, cudaMemcpyDeviceToDevice ) );
}	

inline void gpu_strsmRUNU( int m, int n, float alpha, p2_t A, p2_t B )
{	
	if( m > 0 && n > 0 )
    {
        cublasStrsm( 'R', 'U', 'N', 'U', m, n, alpha, A.A, A.lda, B.A, B.lda );
        Q( cublasGetError( ) );
    }
}

//
//  Custom GPU kernels
//
extern "C" void gpu_sgemmNN( int m, int n, int k, float alpha, const p2_t A, const p2_t B, float beta, p2_t C );
extern "C" void gpu_sgemmNT( int m, int n, int k, float alpha, const p2_t A, const p2_t B, float beta, p2_t C );
extern "C" void gpu_ssyrkLN( int n, int k, float alpha, const p2_t A, float beta, p2_t C );
extern "C" void gpu_batch_sswap( int batchsize, int n, int *ipiv, p2_t matrix );
extern "C" void gpu_enforceLU( int n, p2_t matrix );
extern "C" void gpu_transpose( int m, int n, p2_t dst, p2_t src );
extern "C" void gpu_transpose_inplace( int n, p2_t matrix );
