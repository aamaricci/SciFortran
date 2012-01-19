// Written by Vasily Volkov.
// Copyright (c) 2009, The Regents of the University of California.
// All rights reserved.

#include "gpu_lapack_internal.h"
#include <mkl_blas.h> // Intel platform is assumed
#include <mkl_lapack.h>

//
//  Cholesky factorization
//  See http://www.netlib.org/lapack/single/spotrf.f
//
void gpu_spotrf( int n, float *cpu_A, int lda, int &info )
{	
    //
    //  partition preallocated GPU memory
    //
    gpu_malloc_reset( );
    p2_t gpu_matrix = gpu_malloc2D( n, n );
    if( gpu_matrix.A == NULL )
    {
        info = -1;
        return;
    }
    
    //
    //  CPU pointer helpers
    //
    p2_t cpu_matrix(  cpu_A, lda );
    
    //
    //  upload matrix
    //
    uploadL( n-nb, gpu_matrix(nb,nb), cpu_matrix(nb,nb) );
    
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
            //  look ahead on GPU
            //
      	    gpu_ssyrkLN( w, nb, -1, gpu_matrix(i,i-nb), 1, gpu_matrix(i,i) );
	        gpu_sgemmNT( h-w, w, nb, -1, gpu_matrix(i+nb,i-nb), gpu_matrix(i,i-nb), 1, gpu_matrix(i+nb,i) );

            //
            //  download panel to CPU
            //
            download( h, w, cpu_matrix(i,i), gpu_matrix(i,i) );
            
            //
            //  enqueue right update on GPU
            //
            if( h > nb )
                gpu_ssyrkLN( h-nb, nb, -1, gpu_matrix(i+nb,i-nb), 1, gpu_matrix(i+nb,i+nb) );
        }
        
        //
        //  factorize diagonal block on CPU
        //
        spotrf( "L", &w, cpu_matrix(i,i).A, &cpu_matrix.lda, &info );
        if( info != 0 )
        {
            info += i;
            return;
        }
        
        if( h > nb )
        {
            //
            //  compute the rest of panel on CPU
            //
            float fone = 1;
            int m = h - nb;
            strsm( "R", "L", "T", "N", &m, &nb, &fone, cpu_matrix(i,i).A, &cpu_matrix.lda, cpu_matrix(i+nb,i).A, &cpu_matrix.lda );
            
            //
            //  upload panel to GPU
            //
            upload( h, w, gpu_matrix(i,i), cpu_matrix(i,i) );
        }
    }
}	
