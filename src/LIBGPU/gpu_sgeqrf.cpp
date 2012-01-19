// Written by Vasily Volkov.
// Copyright (c) 2009, The Regents of the University of California.
// All rights reserved.

#include "gpu_lapack_internal.h"
#include <mkl_lapack.h> // Intel platform is assumed

//
//  QR factorization
//  See http://www.netlib.org/lapack/single/sgetrf.f
//
void gpu_sgeqrf( int n, float *cpu_A, int lda, float *cpu_tau, float *cpu_work, int lwork, int &info )
{	
    //
    //  respond on workspace query
    //
    if( lwork == -1 )
    {
        int w = nb;
        sgeqrf( &n, &w, cpu_A, &lda, cpu_tau, cpu_work, &lwork, &info );
        return;
    }
    
    //
    //  partition preallocated GPU memory
    //
    gpu_malloc_reset( );
    p2_t gpu_matrix = gpu_malloc2D( n,  n );
    p2_t gpu_TV     = gpu_malloc2D( nb, n );
    p2_t gpu_TVA    = gpu_malloc2D( nb, n );
    if( gpu_matrix.A == NULL || gpu_TV.A == NULL || gpu_TVA.A == NULL )
    {
        info = -1;
        return;
    }
    p2_t gpu_T = gpu_TVA;
    
    //
    //  CPU pointer helpers
    //
    p2_t cpu_matrix( cpu_A, lda );
    p2_t cpu_T = p2_t( cpu_buffer, nb );
    
    //
    //  upload matrix to GPU
    //
    upload( n, n-nb, gpu_matrix( 0, nb ), cpu_matrix( 0, nb ) );
    
    //
    //  iterate through block columns
    //
    for( int i = 0; i < n; i += nb )
    {
        int h = n - i;
        int w = h < nb ? h : nb;
        
        if( i > 0 )
        {
            //
            //  look ahead on GPU
            //
            gpu_sgemmNN( nb, w, h+nb, 1, gpu_TV, gpu_matrix(i-nb,i), 0, gpu_TVA );
            gpu_sgemmNN( h+nb, w, nb, -1, gpu_matrix(i-nb,i-nb), gpu_TVA, 1, gpu_matrix(i-nb,i) );
            
            //
            //  download panel to CPU
            //
            download( n, w, cpu_matrix(0,i), gpu_matrix(0,i) );
			
            //
            //  enqueue right update on GPU
            //
            gpu_sgemmNN( nb, h-nb, h+nb, 1, gpu_TV, gpu_matrix(i-nb,i+nb), 0, gpu_TVA );
            gpu_sgemmNN( h+nb, h-nb, nb, -1, gpu_matrix(i-nb,i-nb), gpu_TVA, 1, gpu_matrix(i-nb,i+nb) );
        }
		
        //
        //  factorize panel on CPU
        //
        sgeqrf( &h, &w, cpu_matrix(i,i).A, &cpu_matrix.lda, cpu_tau+i, cpu_work, &lwork, &info );
        
        if( h > nb )
        {
            //
            //  generate the triangular factor T on CPU; transpose to avoid sgemmTT/sgemmTN
            //
            slarft( "F", "C", &h, &w, cpu_matrix(i,i).A, &cpu_matrix.lda, cpu_tau+i, cpu_T.A, &cpu_T.lda );
            for( int j = 0; j < nb; j++ )
                for( int k = 0; k < j; k++ )
                {
                    cpu_T.at(j,k) = cpu_T.at(k,j);
                    cpu_T.at(k,j) = 0;
                }
            
            //
            //  upload T and V to GPU
            //
            upload( nb, nb, gpu_T, cpu_T );
            upload( h, nb, gpu_matrix(i,i), cpu_matrix(i,i) );
    		
            //
            //  compute T^T * V^T on GPU
            //
            gpu_enforceLU( nb, gpu_matrix(i,i) );
            gpu_sgemmNT( nb, h, nb, 1, gpu_T, gpu_matrix(i,i), 0, gpu_TV );
        }
    }
    
    info = 0;
}
