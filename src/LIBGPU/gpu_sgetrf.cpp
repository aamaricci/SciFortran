// Written by Vasily Volkov.
// Copyright (c) 2009, The Regents of the University of California.
// All rights reserved.

#include "gpu_lapack_internal.h"
#include <string.h>
#include <mkl_blas.h> // Intel platform is assumed
#include <mkl_lapack.h>

//
//  LU factorization
//  See http://www.netlib.org/lapack/single/sgetrf.f
//
//  set fast_trsm = true  to do TRSM via inv+GEMM
//  set fast_trsm = false for the standard approach
//
void gpu_sgetrf( int n, float *cpu_A, int lda, int *cpu_ipiv, int &info, bool fast_trsm )
{	
    //
    //  partition preallocated GPU memory
    //
    gpu_malloc_reset( );
    p2_t gpu_matrix = gpu_malloc2D( n,  n );
    p2_t gpu_buffer = gpu_malloc2D( n,  nb );
    p2_t gpu_L      = gpu_malloc2D( nb, nb );
    if( gpu_matrix.A == NULL || gpu_buffer.A == NULL || gpu_L.A == NULL)
    {
        info = -1;
        return;
    }
    
    //
    //  CPU pointer helpers
    //
    p2_t cpu_matrix( cpu_A, lda );
    p2_t cpu_L( cpu_buffer, nb );
    
    //
    //  upload matrix to GPU
    //
    upload( n, n-nb, gpu_matrix(0,nb), cpu_matrix(0,nb) );
    gpu_transpose_inplace( n, gpu_matrix );
    
    //
    //  iterate through block columns
    //
    info = 0;
    for( int i = 0; i < n; i += nb )
    {
        int h = n - i;
        int w = h < nb ? h : nb;
        
        if( i > 0 )
        {
            //
            //  download panel to CPU
            //
            gpu_transpose( w, h, gpu_buffer, gpu_matrix(i,i) );
            download(  h, w, cpu_matrix(i,i), gpu_buffer );
            
            //
            //  enqueue right update on GPU
            //
            if( fast_trsm )
            {
                gpu_copy( h-nb, nb, gpu_buffer, gpu_matrix(i+nb,i-nb) );
                gpu_sgemmNT( h-nb, nb, nb, 1, gpu_buffer, gpu_L, 0, gpu_matrix(i+nb,i-nb) );
            }
            else
                gpu_strsmRUNU( h-nb, nb, 1, gpu_matrix(i-nb,i-nb), gpu_matrix(i+nb,i-nb) );
            gpu_sgemmNN( h-nb, h, nb, -1, gpu_matrix(i+nb,i-nb), gpu_matrix(i-nb,i), 1, gpu_matrix(i+nb,i) );
        }
        
        //
        //  factorize panel on CPU
        //
        int iinfo;
        sgetrf( &h, &w, cpu_matrix(i,i).A, &cpu_matrix.lda, cpu_ipiv+i, &iinfo );
        if( info == 0 && iinfo > 0 )
            info = i + iinfo;
        
        //
        //  pivot matrix on GPU
        //
        gpu_batch_sswap( w, n, cpu_ipiv+i, gpu_matrix(0,i) );
        for( int j = 0; j < w; j++ )
            cpu_ipiv[i+j] += i;
        
        //
        //  invert L's diagonal block on CPU
        //
        if( fast_trsm && h > nb )
        {
            //
            //  make identity matrix I
            //
            memset( cpu_L.A, 0, sizeof(float)*nb*nb );
            for( int j = 0; j < nb; j++ )
                cpu_L.at(j,j) = 1;
            
            //
            //  solve LX = I on CPU
            //
            float fone = 1;
            strsm( "L", "L", "N", "U", &nb, &nb, &fone, cpu_matrix(i,i).A, &cpu_matrix.lda, cpu_L.A, &cpu_L.lda );
            upload( nb, nb, gpu_L, cpu_L );
        }
        
        //
        //  upload panel to GPU
        //
        upload( h, w, gpu_buffer, cpu_matrix(i,i) );
        gpu_transpose( h, w, gpu_matrix(i,i), gpu_buffer );
        
        //
        //  look-ahead on GPU
        //
        if( h > nb )
        {
            if( fast_trsm )
            {
                gpu_copy( nb, nb, gpu_buffer, gpu_matrix(i+nb,i) );
                gpu_sgemmNT( nb, nb, nb, 1, gpu_buffer, gpu_L, 0, gpu_matrix(i+nb,i) );
            }
            else
                gpu_strsmRUNU( nb, nb, 1, gpu_matrix(i,i), gpu_matrix(i+nb,i) );
            gpu_sgemmNN( nb, h-nb, nb, -1, gpu_matrix(i+nb,i), gpu_matrix(i,i+nb), 1, gpu_matrix(i+nb,i+nb) );
        }
    }
	
    //
    //  download matrix back to CPU
    //
    gpu_transpose_inplace( n, gpu_matrix );
    download( n, n, cpu_matrix, gpu_matrix );
}	
