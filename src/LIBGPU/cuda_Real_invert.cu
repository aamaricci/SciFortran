
///********************************************************************
//*  File: cuda_invert.cu
//*
//*  Description:
//*     
//*   Mainfunction to compute an inverse Matrix from a positive definite 
//*   Matrix on the GPU. The Routine is using the Gaus Seidel Matrix invertion 
//*   algorithm. 
//*   
//*   1   2   1   |  1  0   0               1   0   0  |  -2.5   1.5   0.5
//*   2   3   1   |  0  1   0       =>      0   1   0  |   1.5  -0.5  -0.5 
//*   1   1   2   |  0  0   1               0   0   1  |   0.5  -0.5   0.5 
//*   Inputmatrix       E                       E          Inverse Matrix
//*
//*  Arguments: 
//*   - double *A            Input Matrix 1D, no data changes
//*   - double *invA         Output Matrix 1D, the inverse datamatrix  
//*   - int size            Matrix dimension in size, size = height = size
//*     
//*  Used custom kernels rutines:
//*   - GPUsetIdentity          
//*   - eliminateBlockL   
//*   - adjustColL         
//*   - eliminateColL      
//*   - eliminateRestL     
//*
//*   - eliminateBlockU    
//*   - adjustColU         
//*   - eliminateColU      
//*   - eliminateRestU     
//*
//*********************************************************************


//********************************************
//********************************************
//********************************************
//********************************************


#include <stdio.h>
#include <stdlib.h>
#include <cutil.h>
#include <cuda_runtime.h>


#define MAXSIZE 16



//************************************************
// Kernel Mul
//************************************************

extern "C" void Mul(int size, const double* A, const double* B, double* C)
{

 double aa,bb;
// c_ij = a_ik * b_kj
for(int i =0; i<size; i++) { for(int j =0; j<size; j++) {
 C[j*size + i]=0.;
 for(int k =0; k<size; k++) 
 {
 aa = A[ k * size + i];
 bb = B[ j * size + k];
 C[ j*size + i ] += aa*bb;
 };};};

}

//************************************************
// Kernel GPUsetIdentity 
//************************************************

__global__ void GPUsetIdentity(int BLOCKSIZE,double* matrix, int size) 
{
    int tx       = threadIdx.x;
    int bx       = blockIdx.x;
    int offset   = bx * BLOCKSIZE + tx;
    matrix[offset * size + offset] = 1 ;
}

//************************************************
// Kernel eliminateBlockL
//************************************************


__global__ void eliminateBlockL(int BLOCKSIZE, double *dInData, int size)
{
    int BLOCKSIZEMINUS1 = BLOCKSIZE-1;
    int tx = threadIdx.x;
    int ty = threadIdx.y;

    __shared__ double triangleBlock[MAXSIZE][MAXSIZE];
    
    triangleBlock[ty][tx] = dInData[ ty * size + tx];
             
    __syncthreads ();
    

    //i equals the current row
    
    for (int i = 0; i < BLOCKSIZEMINUS1; i++)
    {
        
        // calculate the pivot element to get the current row i to zero
        
        
        double pivotEl = triangleBlock[ty][i] / triangleBlock[i][i];

        __syncthreads ();       // Each pivotEl have to be calculated and store in the registers

        if (ty > i)             // If all cols (ty) are below the current row (i)?
        {
             if (tx > i)         // The element is right to the current row, subtract the element
            {
                 triangleBlock[ty][tx] -= pivotEl * triangleBlock[i][tx];
            }
                
            if (tx == i)        // Store the pivot element in the current row
            {
                triangleBlock[ty][tx] = pivotEl;
            }
                
        }
            
        __syncthreads ();       // Wait for each thread
          
      }
      
           dInData[ty * size + tx] = triangleBlock[ty][tx]; // Write the result back to memory
}

//************************************************
// Kernel eliminateBlockU
//************************************************

__global__ void eliminateBlockU(int BLOCKSIZE,double *dInData, int size)
{
     int BLOCKSIZEMINUS1 = BLOCKSIZE-1;
          
    int tx = threadIdx.x;
    int ty = threadIdx.y;

    __shared__ double triangleBlock[MAXSIZE][MAXSIZE];

    triangleBlock[ty][tx] = dInData[ty * size + tx];
    
    __syncthreads ();


    //i equals the current row

    for (int i = BLOCKSIZEMINUS1; i > 0; i--)
    {
        // calculate the pivot element to get the current row i to zero

        double pivotEl = triangleBlock[ty][i] / triangleBlock[i][i];

        __syncthreads ();       // Each pivotEl have to be calculated and store in the registers

        if (ty < i)             // If all rows (ty) are above the current row (i)?
        {

            if (tx < i)         // The element is left to the current row, subtract the element
            {
                triangleBlock[ty][tx] -= pivotEl * triangleBlock[i][tx];
            }

            if (tx == i)        // Store the pivot element in the current row
            {
                triangleBlock[ty][tx] = pivotEl;
            }

        }

        __syncthreads ();        // Wait for each thread

    }

         dInData[ty * size + tx] = triangleBlock[ty][tx];       //Write the result back to device memory
    
}

//************************************************
// Kernel adjustRowL
//************************************************

                                                   
__global__ void adjustRowL(int BLOCKSIZE, double *dMatrixIn, double *dMatrixInDiag, double *dMatrixInv, int size, int diagEl)
{
     int BLOCKSIZEMINUS1 = BLOCKSIZE-1;
     
    int tx = threadIdx.x;
    int ty = threadIdx.y;
    int bx = blockIdx.x;

    __shared__ double pivotBlock[MAXSIZE][MAXSIZE];
    __shared__ double inBlock[MAXSIZE][MAXSIZE];
    __shared__ double invBlock[MAXSIZE][MAXSIZE];


    //*
    //* Adjust the rest blocks which are right from the prepared block of step 1
    //* and adjust the inverse blocks
    //*
    
    if (bx * BLOCKSIZE > (diagEl + 1))
    {

       pivotBlock[ty][tx] = dMatrixInDiag[ty * size + tx];
          inBlock[ty][tx] = dMatrixIn[ty  * size + bx * BLOCKSIZE +tx];
         invBlock[ty][tx] = dMatrixInv[ty * size + bx * BLOCKSIZE +tx];

        __syncthreads ();

        //i equals the current row where the pivot elements are stored

        for (int i = 0; i < BLOCKSIZEMINUS1; i++)
        {
            // if the cols are below  
            
            if (ty > i)
            {
                double pivot = pivotBlock[ty][i];
                
                //Subtract the row
                
                inBlock[ty][tx]  -= inBlock[i][tx]  * pivot;
                invBlock[ty][tx] -= invBlock[i][tx] * pivot;
            }

            __syncthreads ();
        }
        
        // Store the results back in device memory

               dMatrixIn[ty * size + bx * BLOCKSIZE +tx] = inBlock[ty][tx];
              dMatrixInv[ty * size + bx * BLOCKSIZE +tx] = invBlock[ty][tx];
         
    }
    
    //* Adjust the last blocks from the indentity matrix which are left 
    
    else
    {
          pivotBlock[ty][tx] = dMatrixInDiag[ty * size + tx];
          invBlock[ty][tx]   = dMatrixInv[ty * size + bx * BLOCKSIZE +tx];
         
        __syncthreads ();
        
        for (int i = 0; i < BLOCKSIZEMINUS1; i++) //last changed
        {
            if (ty > i)
            { invBlock[ty][tx] -= invBlock[i][tx] * pivotBlock[ty][i];}

            __syncthreads ();
        }
             dMatrixInv[ty * size + bx * BLOCKSIZE + tx] = invBlock[ty][tx];
    }
    
}

//************************************************
// Kernel adjustRowU
//************************************************

__global__ void adjustRowU(int BLOCKSIZE,double *dMatrixIn, double *dMatrixInv, int size, int diagEl)
{
    int BLOCKSIZEMINUS1 = BLOCKSIZE-1;
          
    int tx = threadIdx.x;
    int ty = threadIdx.y;
    int bx = blockIdx.x;

    __shared__ double pivotBlock[MAXSIZE][MAXSIZE];
    __shared__ double invBlock[MAXSIZE][MAXSIZE];

     pivotBlock[ty][tx] =  dMatrixIn[ty * size + tx];
       invBlock[ty][tx] = dMatrixInv[ty * size + bx * BLOCKSIZE +tx];
    
    __syncthreads ();

    for (int i = BLOCKSIZEMINUS1; i > 0; i--)
    {
        if (ty < i)
        { invBlock[ty][tx] -= invBlock[i][tx] * pivotBlock[ty][i];}

        __syncthreads ();
    }

         dMatrixInv[ty * size + bx * BLOCKSIZE +tx] = invBlock[ty][tx];

        __syncthreads ();

}

//************************************************
// Kernel eliminateColL
//************************************************

__global__ void eliminateColL(int BLOCKSIZE, double *dMatrixIn, int size, int diagEl)
{
     
    int tx = threadIdx.x; int ty = threadIdx.y;

    // bx is used to adress the Blocks above the precalculated block from step 1
    
    int bx = blockIdx.x;    

    //only the blocks can enter which are above the precalculated block from step 1
    
    if (bx * BLOCKSIZE > (diagEl + 1))
    {
        int offset = diagEl * size;
        int blockOffset = bx * BLOCKSIZE *size;

        __shared__ double pivotBlock[MAXSIZE][MAXSIZE];
        __shared__ double inBlock[MAXSIZE][MAXSIZE];

        pivotBlock[ty][tx] = dMatrixIn[offset + ty * size + tx];       // The Block from step 1
        inBlock[ty][tx]    = dMatrixIn[blockOffset + ty * size + tx];  // each Block which is above the pivotBlock
        
        __syncthreads ();
    
        //iterate through the block und calculate the pivot elements
    
        for (int i = 0; i < BLOCKSIZE; i++)
        {
            double pivotEl = inBlock[ty][i] / pivotBlock[i][i];

            __syncthreads ();

            //adjust all values right to the current interation step
           
            if (tx > i)
            {
                //substract the row
                inBlock[ty][tx] -= pivotBlock[i][tx] * pivotEl;
            }
           
            //store the pivot element in the col
           
            else
            {
                inBlock[ty][i] = pivotEl;
            }

            __syncthreads ();
        }

          dMatrixIn[blockOffset + ty * size + tx] = inBlock[ty][tx];
        
    }
}


//************************************************
// Kernel eliminateColU
//************************************************

__global__ void eliminateColU(int BLOCKSIZE, double *dMatrixIn, int size, int diagEl)
{
    int BLOCKSIZEMINUS1 = BLOCKSIZE-1;
          
    int tx = threadIdx.x; int ty = threadIdx.y; int bx = blockIdx.x;

    if (bx * BLOCKSIZE <diagEl)
    {
        int offset = diagEl * size;
        int blockOffset = bx * BLOCKSIZE *size;

        __shared__ double pivotBlock[MAXSIZE][MAXSIZE];
        __shared__ double inBlock[MAXSIZE][MAXSIZE];

        pivotBlock[ty][tx] = dMatrixIn[offset + ty * size + tx];
        inBlock[ty][tx]    = dMatrixIn[blockOffset + ty * size + tx];
        
        __syncthreads ();


        for (int i = BLOCKSIZEMINUS1; i >= 0; i--)
        {
            double pivotEl = inBlock[ty][i] / pivotBlock[i][i];

            __syncthreads ();

            if (tx < i)
            {
                inBlock[ty][tx] -= pivotBlock[i][tx] * pivotEl;
            }
            else //* if (tx == i)
            {
                inBlock[ty][i] = pivotEl;
            }

            __syncthreads ();
        }

             dMatrixIn[blockOffset + ty * size + tx] = inBlock[ty][tx];
        
    }
}


//************************************************
// Kernel eliminateRestL
//************************************************

__global__ void eliminateRestL(int BLOCKSIZE, double *dMatrixIn, double *dMatrixInv, int size,int diagEl)
{
     
    int tx = threadIdx.x; int ty = threadIdx.y;
    int bx = blockIdx.x;  int by = blockIdx.y;

    __shared__ double pivEl[MAXSIZE][MAXSIZE];
    __shared__ double pivBlock[MAXSIZE][MAXSIZE];
    __shared__ double inBlock[MAXSIZE][MAXSIZE];

    //rest of the unadjusted Matrix which is right above the diagEl

    if (bx * BLOCKSIZE > (diagEl + 1) && by * BLOCKSIZE > (diagEl + 1))
    {

        int blockOffset = by * BLOCKSIZE * size + bx * BLOCKSIZE;
        int blockPivElOffset = by * BLOCKSIZE * size + diagEl;
        int blockPivOffset = diagEl * size + bx * BLOCKSIZE;

        inBlock[ty][tx] = dMatrixIn[blockOffset + ty * size + tx];
        pivEl[ty][tx] = dMatrixIn[blockPivElOffset + ty * size + tx];
        pivBlock[ty][tx] = dMatrixIn[blockPivOffset + ty * size + tx];
        
        __syncthreads ();

        //Subtract each row from the input Matrix =>dMatrixIn
        
        for (int i = 0; i < BLOCKSIZE; i++)
        {
            inBlock[ty][tx] -= pivEl[ty][i] * pivBlock[i][tx];
        }
        
        __syncthreads ();
        
        if( (blockOffset + ty*size + tx) < (size*size) )
        {
             dMatrixIn[blockOffset + ty * size + tx] = inBlock[ty][tx];
        }
        
        __syncthreads ();

             inBlock[ty][tx] = dMatrixInv[blockOffset + ty * size + tx];
             pivBlock[ty][tx] = dMatrixInv[blockPivOffset + ty * size + tx];
        
        __syncthreads ();

        //Subtract each row from the invers Matrix =>dMatrixInv
        
        for (int i = 0; i < BLOCKSIZE; i++)
        {
            inBlock[ty][tx] -= pivEl[ty][i] * pivBlock[i][tx];
        }

        __syncthreads ();
        
        dMatrixInv[blockOffset + ty * size + tx] = inBlock[ty][tx];
        
    }
    
    //Adjust the left Blocks from the invers matrix which are left from the diagEl
    
    else if (by * BLOCKSIZE > (diagEl + 1))
    {
        int blockOffset = by * BLOCKSIZE * size + bx * BLOCKSIZE;
        int blockPivElOffset = by * BLOCKSIZE * size + diagEl;
        int blockPivOffset = diagEl * size + bx * BLOCKSIZE;

        pivEl[ty][tx] = dMatrixIn[blockPivElOffset + ty * size + tx];
        inBlock[ty][tx] = dMatrixInv[blockOffset + ty * size + tx];
        pivBlock[ty][tx] = dMatrixInv[blockPivOffset + ty * size + tx];
        
        __syncthreads ();

        for (int i = 0; i < BLOCKSIZE; i++)
        {
            inBlock[ty][tx] -= pivEl[ty][i] * pivBlock[i][tx];
        }

        __syncthreads ();

        dMatrixInv[blockOffset + ty * size + tx] = inBlock[ty][tx];
        
    }
}

//************************************************
// Kernel eliminateRestU
//************************************************

__global__ void eliminateRestU(int BLOCKSIZE, double *dMatrixIn, double *dMatrixInv, int size,int diagEl)
{
    int BLOCKSIZEMINUS1 = BLOCKSIZE-1;
          
    int tx = threadIdx.x ; int ty = threadIdx.y;
    int bx = blockIdx.x  ; int by = blockIdx.y;

    __shared__ double pivEl[MAXSIZE][MAXSIZE];
    __shared__ double pivBlock[MAXSIZE][MAXSIZE];
    __shared__ double inBlock[MAXSIZE][MAXSIZE];

    if ((bx * BLOCKSIZE + 1) <diagEl && (by * BLOCKSIZE +1) <diagEl)   
    {

        int blockOffset = by * BLOCKSIZE * size + bx * BLOCKSIZE;
        int blockPivElOffset = by * BLOCKSIZE * size + diagEl;
        int blockPivOffset = diagEl * size + bx * BLOCKSIZE;

        inBlock[ty][tx] = dMatrixIn[blockOffset + ty * size + tx];
        pivEl[ty][tx] = dMatrixIn[blockPivElOffset + ty * size + tx];
        pivBlock[ty][tx] = dMatrixIn[blockPivOffset + ty * size + tx];
        
        __syncthreads ();

        for (int i = BLOCKSIZEMINUS1; i >= 0; i--)
        {
            inBlock[ty][tx] -= pivEl[ty][i] * pivBlock[i][tx];
        }
        
        __syncthreads ();
        
        if( (blockOffset + ty*size + tx) < (size*size) )
        { 
             dMatrixIn[blockOffset + ty * size + tx] = inBlock[ty][tx];
        }
        
        __syncthreads ();

         inBlock[ty][tx] = dMatrixInv[blockOffset + ty * size + tx];
         pivBlock[ty][tx] = dMatrixInv[blockPivOffset + ty * size + tx];
        
        __syncthreads ();

        for (int i = BLOCKSIZEMINUS1; i >= 0; i--)
        {
            inBlock[ty][tx] -= pivEl[ty][i] * pivBlock[i][tx];
        }

        __syncthreads ();
        
        dMatrixInv[blockOffset + ty * size + tx] = inBlock[ty][tx];
        
    }
    
    
    else if (by * BLOCKSIZE < (diagEl))
    {
        int blockOffset = by * BLOCKSIZE *size + bx * BLOCKSIZE;
        int blockPivElOffset = by * BLOCKSIZE *size + diagEl;
        int blockPivOffset = diagEl * size + bx * BLOCKSIZE;

        pivEl[ty][tx] = dMatrixIn[blockPivElOffset + ty * size + tx];
        inBlock[ty][tx] = dMatrixInv[blockOffset + ty * size + tx];
        pivBlock[ty][tx] = dMatrixInv[blockPivOffset + ty * size + tx];
        
        __syncthreads ();

        for (int i = BLOCKSIZEMINUS1; i >= 0; i--)
        {
            inBlock[ty][tx] -= pivEl[ty][i] * pivBlock[i][tx];
        }

        __syncthreads ();

        dMatrixInv[blockOffset + ty * size + tx] = inBlock[ty][tx];
          
    }
}

//************************************************
// Kernel normalizeDiag
//************************************************

__global__ void normalizeDiag(int BLOCKSIZE,double *diagMatrix, double *invMatrix, int size,int row)
{
          
    int tx = threadIdx.x; int ty = threadIdx.y; int bx = blockIdx.x;

    int blockOffset = bx * BLOCKSIZE;

    __shared__ double diagEl[MAXSIZE];

    if (tx == ty ) { diagEl[ty] = diagMatrix[row + ty * size + tx]; }

    __syncthreads ();

    invMatrix[blockOffset + ty * size + tx] =
    invMatrix[blockOffset + ty * size + tx] / diagEl[ty];
    
}

//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
// Wrapper Fortran / C
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************

 extern "C" void cuda_invert__(int *pBlockSize, double** A_, double** invA_, int *pndim)
{

    int i; int size = *pndim; int BlockSize = *pBlockSize;

    if(BlockSize>MAXSIZE) {printf("error Block size cuda invert too big");return;}
    dim3 idyThreads  ( BlockSize );    
    dim3 idyBlocks   ( size / BlockSize );
    dim3 nThreads    ( BlockSize, BlockSize );   
    dim3 nBlocks     ( size / BlockSize );        
    dim3 nBlocksRest ( size / BlockSize, size / BlockSize );
 
    int pinned=1 ; if(size<7200) pinned=0; 

//  printf(" BlockSize=%d , size=%d \n \n \n ", BlockSize , size); 
//  printf(" number of blocks : %d %d \n ", size/BlockSize, size/BlockSize  ); 

    cudaThreadExit();

    int mat_size=size*size*sizeof(double);
    double *invA ; double *A  ; double *invAd ; double *Ad;

//---------------------------------------------------------------//
    if(pinned==0){
//    printf(" copy memory host to gpu  \n ");
      cudaMalloc((void**) &invA, mat_size);
      cudaMalloc((void**) &A   , mat_size);
      cudaMemcpy(A   , A_   , mat_size, cudaMemcpyHostToDevice);
      cudaMemcpy(invA, invA_, mat_size, cudaMemcpyHostToDevice);
    }else{
      cudaSetDevice(0); cudaSetDeviceFlags( cudaDeviceMapHost ); 
      cudaHostAlloc( (void**) &Ad,    mat_size, cudaHostAllocMapped | cudaHostAllocPortable );
      cudaHostAlloc( (void**) &invAd, mat_size, cudaHostAllocMapped | cudaHostAllocPortable );
      cudaHostGetDevicePointer((void**) &A   , Ad,    0 );
      cudaHostGetDevicePointer((void**) &invA, invAd, 0 );
      cudaMemcpy(Ad,    A_,    mat_size, cudaMemcpyHostToHost);
      cudaMemcpy(invAd, invA_, mat_size, cudaMemcpyHostToHost);
      cudaEventSynchronize(0);
    }
//---------------------------------------------------------------//

//  Calculate the Identitymatrix 
    GPUsetIdentity<<<idyBlocks,idyThreads>>>(BlockSize,invA,size);

    cudaThreadSynchronize();

  //calculate the right diagonal Matrix (L)

    for (i = 0; i < size; i += BlockSize)
    {
        int offset = i * size + i;
        
        // * step 1:
        // *  calculate the triangle matrix
        // *  store the pivot elements to left part of the triangel

        eliminateBlockL<<< 1, nThreads >>> (BlockSize,A + offset, size);
        cudaThreadSynchronize ();
    
        // * step 2:
        // *  calculate the rest of the rows with the pivot elements from step 1
        
        adjustRowL<<< nBlocks, nThreads >>> (BlockSize,A + i * size, A + offset, invA + i * size, size, i);
        cudaThreadSynchronize ();
    
        //* step 3:
        //*    Fill the colls below the block with the pivot elements they are used
        //*    to get the colls to zero and multiply with the row
        
        eliminateColL<<< nBlocks, nThreads >>> (BlockSize,A + i, size, i);
        cudaThreadSynchronize ();

        //* step 4:
        //*  Adjust the rest of the Matrix with the calculated pivot Elements
        //*  El_new_0 -= (p0+p1+p2..+p15) * El_piv_0
        
        eliminateRestL<<< nBlocksRest, nThreads >>> (BlockSize, A, invA, size, i);
        cudaThreadSynchronize ();

    }

    cudaEventSynchronize(0);


    //Set the left lower diagonalmatrix to zero (async?)
    
    for (i = 1; i < size; i++)
    { 
      int offset = i * size;
      if(pinned==0){
         cudaMemset ((void *) (A  + offset), 0, i*sizeof(double));}
      else{
         memset ((void *) (Ad  + offset), 0, i*sizeof(double));
      };
      cudaEventSynchronize(0);
    }
    cudaThreadSynchronize ();


    //calculate the right diagonal Matrix (U)
    
    for (i = (size - BlockSize); i >= 0; i -= BlockSize)
    {
        int offset = i * size + i;

        //* step 1:
        //*  calculate the triangle matrix
        //*  store the pivot elements to left part of the triangel
        
        eliminateBlockU<<< 1, nThreads >>> (BlockSize,A + offset, size);
        cudaThreadSynchronize ();

        //* step 2:
        //*  calculate the rest of the rows with the pivot elements from step 1

        int rowOffset = i * size;
        adjustRowU<<< nBlocks, nThreads >>> (BlockSize,A + offset,invA + rowOffset, size, i);
        cudaThreadSynchronize ();

        //* step 3:
        //*  Fill the colls below the block with the pivot elements they are used
        //*      to get the colls to zero and multiply with the row
  
        eliminateColU<<< nBlocks, nThreads >>> (BlockSize,A + i, size, i);
        cudaThreadSynchronize ();

        //* step 4:
        //*  Adjust the rest of the Matrix with the calculated pivot Elements
        //*  El_new_0 -= (p0+p1+p2..+p15) * El_piv_0
        
        eliminateRestU<<< nBlocksRest, nThreads >>> (BlockSize,A, invA, size, i);
        cudaThreadSynchronize ();

    }
 
     cudaEventSynchronize(0);
 
    //* force the diagonal entries to 1

    for (i = 0; i < size; i += BlockSize)
    {
        int rowOffset = i * size;
        normalizeDiag<<< nBlocks, nThreads >>> (BlockSize,A + rowOffset, invA + rowOffset, size, i);
        cudaThreadSynchronize ();
    }

//---------------------------------------------------------------//
  if(pinned==0){
     cudaMemcpy(invA_,invA,mat_size,cudaMemcpyDeviceToHost);
     cudaFree(invA); 
     cudaFree(A);  
   }else{
     cudaMemcpy(invA_,invAd, mat_size,cudaMemcpyHostToHost);
     cudaFreeHost(Ad); 
     cudaFreeHost(invAd); 
     cudaEventSynchronize(0);
   };
//---------------------------------------------------------------//


//  printf(" SUCCESFULL DIAG. \n \n \n ");

    return;
}

//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
