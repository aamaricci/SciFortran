// Written by Vasily Volkov.
// Copyright (c) 2009, The Regents of the University of California.
// All rights reserved.

#pragma once

bool gpu_lapack_init( );

void gpu_spotrf( int n, float *A, int lda, int &info );
void gpu_sgetrf( int n, float *A, int lda, int *ipiv, int &info, bool fast_trsm = true );
void gpu_sgeqrf( int n, float *A, int lda, float *tau, float *work, int lwork, int &info );

void gpu_lapack_release( );
