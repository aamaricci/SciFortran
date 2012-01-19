
///********************************************************************
//*  File: cuda_invert.cu
//*
//*  Description:
//*     
//*     Mainfunction to compute an inverse Matrix from a positive definite 
//*   Matrix on the GPU. The Routine is using the Gaus Seidel Matrix invertion 
//*   algorithm. 
//*   
//*   1   2   1   |  1  0   0               1   0   0  |  -2.5   1.5   0.5
//*   2   3   1   |  0  1   0       =>      0   1   0  |   1.5  -0.5  -0.5 
//*   1   1   2   |  0  0   1               0   0   1  |   0.5  -0.5   0.5 
//*   Inputmatrix       E                       E          Inverse Matrix
//*
//*  Arguments: 
//*       - float *A            Input Matrix 1D, no data changes
//*   - float *invA         Output Matrix 1D, the inverse datamatrix  
//*   - int size            Matrix dimension in size, width = height = size
//*     
//*  Used custom kernels rutines:
//*       - GPUsetIdentity          
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


#include <stdio.h>
#include <stdlib.h>
#include <cutil.h>
#include <cuComplex.h>
#include <cuda_runtime.h>

#define MAXSIZE 16

//************************************************
// Kernel GPUsetIdentity 
//************************************************

__global__ void GPUsetIdentity(int BLOCKSIZE,cuDoubleComplex* matrix,int width)
{

    int tx = threadIdx.x; int bx = blockIdx.x;

    int offset = bx * BLOCKSIZE + tx;
    
    matrix[offset * width + offset] = make_cuDoubleComplex(1.0,0.0);

}

//************************************************
// Kernel eliminateBlockL
//************************************************

__global__ void eliminateBlock(char S, int BLOCKSIZE, cuDoubleComplex *dInData, int size)
{
     int BLOCKSIZEMINUS1 = BLOCKSIZE-1;
     
    int tx = threadIdx.x; int ty = threadIdx.y;

    __shared__ cuDoubleComplex triangleBlock[MAXSIZE][MAXSIZE];

    if(S == 'L')
    {
         triangleBlock[ty][tx] = dInData[ ty * size + tx];
         
         __syncthreads ();

        //i equals the current row
    
        for (int i = 0; i < BLOCKSIZEMINUS1; i++)
        {
        
        // calculate the pivot element to get the current row i to zero
        
        
            cuDoubleComplex pivotEl = cuCdiv( triangleBlock[ty][i] , triangleBlock[i][i] );

            __syncthreads ();       // Each pivotEl have to be calculated and store in the registers

            if (ty > i)             // If all cols (ty) are below the current row (i)?
            {
                 if (tx > i)         // The element is right to the current row, subtract the element
                {
                     triangleBlock[ty][tx] = cuCsub(triangleBlock[ty][tx], cuCmul(pivotEl,triangleBlock[i][tx]));
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
    
    if(S == 'U')
    {

         triangleBlock[ty][tx] = dInData[ty * size + tx];
    
    __syncthreads ();

    //i equals the current row

    for (int i = BLOCKSIZEMINUS1; i > 0; i--)
    {
        // calculate the pivot element to get the current row i to zero

        cuDoubleComplex pivotEl = cuCdiv(triangleBlock[ty][i] , triangleBlock[i][i]);

        __syncthreads ();       // Each pivotEl have to be calculated and store in the registers

        if (ty < i)             // If all rows (ty) are above the current row (i)?
        {

            if (tx < i)         // The element is left to the current row, subtract the element
            {
                triangleBlock[ty][tx] = cuCsub(triangleBlock[ty][tx],cuCmul(pivotEl,triangleBlock[i][tx]));
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
    
}

//************************************************
// Kernel adjustRowL
//************************************************

                                                   
__global__ void adjustRowL(int BLOCKSIZE, cuDoubleComplex *dMatrixIn, cuDoubleComplex *dMatrixInDiag, cuDoubleComplex *dMatrixInv, int width, int diagEl)
{
    int BLOCKSIZEMINUS1 = BLOCKSIZE-1;
     
    int tx = threadIdx.x; int ty = threadIdx.y; int bx = blockIdx.x;

    __shared__ cuDoubleComplex pivotBlock[MAXSIZE][MAXSIZE];
    __shared__ cuDoubleComplex inBlock[MAXSIZE][MAXSIZE];
    __shared__ cuDoubleComplex invBlock[MAXSIZE][MAXSIZE];

    //* Adjust the rest blocks which are right from the prepared block of step 1 and adjust the inverse blocks
    
    if (bx * BLOCKSIZE > (diagEl + 1))
    {

        pivotBlock[ty][tx] = dMatrixInDiag[ty * width + tx];
        inBlock[ty][tx]    = dMatrixIn[ty * width + bx * BLOCKSIZE +tx];
        invBlock[ty][tx]   = dMatrixInv[ty * width + bx * BLOCKSIZE +tx];

        __syncthreads ();

        //i equals the current row where the pivot elements are stored

        for (int i = 0; i < BLOCKSIZEMINUS1; i++)
        {
            // if the cols are below  
            
            if (ty > i)
            {
                cuDoubleComplex pivot = pivotBlock[ty][i];
                
                //Subtract the row
                
                inBlock[ty][tx]  = cuCsub(inBlock[ty][tx], cuCmul(inBlock[i][tx],pivot));
                invBlock[ty][tx] = cuCsub(invBlock[ty][tx], cuCmul(invBlock[i][tx],pivot));
            }

            __syncthreads ();
        }
        
        // Store the results back in device memory

        dMatrixIn[ty * width + bx * BLOCKSIZE +tx] = inBlock[ty][tx];
        dMatrixInv[ty * width + bx * BLOCKSIZE +tx] = invBlock[ty][tx];
         
    }
    
    //* Adjust the last blocks from the indentity matrix which are left 
    
    else
    {
      pivotBlock[ty][tx] = dMatrixInDiag[ty * width + tx];
        invBlock[ty][tx] = dMatrixInv[ty * width + bx * BLOCKSIZE +tx];
         
        __syncthreads ();
        
        for (int i = 0; i < BLOCKSIZEMINUS1; i++) //last changed
        {
            if (ty > i)
            {
                cuDoubleComplex pivot = pivotBlock[ty][i];
                invBlock[ty][tx] = cuCsub(invBlock[ty][tx] ,cuCmul(invBlock[i][tx],pivot));
            }

            __syncthreads ();
        }

             dMatrixInv[ty * width + bx * BLOCKSIZE + tx] = invBlock[ty][tx];
                
    }
    
}

//************************************************
// Kernel adjustRowU
//************************************************


__global__ void adjustRowU(int BLOCKSIZE,cuDoubleComplex *dMatrixIn, cuDoubleComplex *dMatrixInv, int width,int diagEl)
{
    int BLOCKSIZEMINUS1 = BLOCKSIZE-1;
          
    int tx = threadIdx.x; int ty = threadIdx.y; int bx = blockIdx.x;

    __shared__ cuDoubleComplex pivotBlock[MAXSIZE][MAXSIZE];
    __shared__ cuDoubleComplex invBlock[MAXSIZE][MAXSIZE];

     pivotBlock[ty][tx] = dMatrixIn[ty * width + tx];
     invBlock[ty][tx] = dMatrixInv[ty * width + bx * BLOCKSIZE +tx];
    
    __syncthreads ();

    for (int i = BLOCKSIZEMINUS1; i > 0; i--)
    {
        if (ty < i)
        {
            cuDoubleComplex pivot = pivotBlock[ty][i];

            invBlock[ty][tx] = cuCsub(invBlock[ty][tx],cuCmul(invBlock[i][tx],pivot));
        }

        __syncthreads ();
    }

     dMatrixInv[ty * width + bx * BLOCKSIZE +tx] = invBlock[ty][tx];

}

//************************************************
// Kernel eliminateColL
//************************************************

__global__ void eliminateColL(int BLOCKSIZE, cuDoubleComplex *dMatrixIn, int size, int diagEl)
{
     
    int tx = threadIdx.x; int ty = threadIdx.y;

    // bx is used to adress the Blocks above the precalculated block from step 1
    
    int bx = blockIdx.x;    

    //only the blocks can enter which are above the precalculated block from step 1
    
    if (bx * BLOCKSIZE > (diagEl + 1))
    {
        int offset = diagEl * size;
        int blockOffset = bx * BLOCKSIZE *size;

        __shared__ cuDoubleComplex pivotBlock[MAXSIZE][MAXSIZE];
        __shared__ cuDoubleComplex inBlock[MAXSIZE][MAXSIZE];

        pivotBlock[ty][tx] = dMatrixIn[offset + ty * size + tx];   // The Block from step 1
           inBlock[ty][tx] = dMatrixIn[blockOffset + ty * size + tx]; // each Block which is above the pivotBlock

        __syncthreads ();
    
        //iterate through the block und calculate the pivot elements
    
        for (int i = 0; i < BLOCKSIZE; i++)
        {
            cuDoubleComplex pivotEl = cuCdiv(inBlock[ty][i],pivotBlock[i][i]);

            __syncthreads ();

            //adjust all values right to the current interation step
           
            if (tx > i)
            {
                //substract the row
                inBlock[ty][tx] = cuCsub(inBlock[ty][tx], cuCmul(pivotBlock[i][tx],pivotEl));
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

__global__ void eliminateColU(int BLOCKSIZE, cuDoubleComplex *dMatrixIn, int size, int diagEl)
{
    int BLOCKSIZEMINUS1 = BLOCKSIZE-1;
          
    int tx = threadIdx.x; int ty = threadIdx.y; int bx = blockIdx.x;

    if (bx * BLOCKSIZE <diagEl)
    {
        int offset = diagEl * size;
        int blockOffset = bx * BLOCKSIZE *size;

        __shared__ cuDoubleComplex pivotBlock[MAXSIZE][MAXSIZE];
        __shared__ cuDoubleComplex inBlock[MAXSIZE][MAXSIZE];

       pivotBlock[ty][tx] = dMatrixIn[offset + ty * size + tx];
          inBlock[ty][tx] = dMatrixIn[blockOffset + ty * size + tx];
        
        __syncthreads ();

        for (int i = BLOCKSIZEMINUS1; i >= 0; i--)
        {
            cuDoubleComplex pivotEl = cuCdiv(inBlock[ty][i],pivotBlock[i][i]);

            __syncthreads ();

            if (tx < i)
            {
                inBlock[ty][tx] = cuCsub(inBlock[ty][tx],cuCmul(pivotBlock[i][tx],pivotEl));
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

__global__ void eliminateRestL(int BLOCKSIZE, cuDoubleComplex *dMatrixIn, cuDoubleComplex *dMatrixInv, int size,int diagEl)
{
     
    int tx = threadIdx.x; int ty = threadIdx.y;
    int bx = blockIdx.x;  int by = blockIdx.y;

    __shared__ cuDoubleComplex pivEl[MAXSIZE][MAXSIZE];
    __shared__ cuDoubleComplex pivBlock[MAXSIZE][MAXSIZE];
    __shared__ cuDoubleComplex inBlock[MAXSIZE][MAXSIZE];

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
            inBlock[ty][tx] = cuCsub(inBlock[ty][tx], cuCmul(pivEl[ty][i],pivBlock[i][tx]));
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
            inBlock[ty][tx] = cuCsub(inBlock[ty][tx], cuCmul(pivEl[ty][i],pivBlock[i][tx]));
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
            inBlock[ty][tx] = cuCsub(inBlock[ty][tx], cuCmul(pivEl[ty][i],pivBlock[i][tx]));
        }

        __syncthreads ();

       dMatrixInv[blockOffset + ty * size + tx] = inBlock[ty][tx];
        
    }
}

//************************************************
// Kernel eliminateRestU
//************************************************

__global__ void eliminateRestU(int BLOCKSIZE, cuDoubleComplex *dMatrixIn, cuDoubleComplex *dMatrixInv, int size,int diagEl)
{
    int BLOCKSIZEMINUS1 = BLOCKSIZE-1;
          
    int tx = threadIdx.x; int ty = threadIdx.y;
    int bx = blockIdx.x;  int by = blockIdx.y;

    __shared__ cuDoubleComplex pivEl[MAXSIZE][MAXSIZE];
    __shared__ cuDoubleComplex pivBlock[MAXSIZE][MAXSIZE];
    __shared__ cuDoubleComplex inBlock[MAXSIZE][MAXSIZE];

    //rest der unbearbeiteten Matrix bearbeiten

    if ((bx * BLOCKSIZE + 1) <diagEl && (by * BLOCKSIZE +1) <diagEl)     //linke seite von in; 0-pivblock
    {
        int blockOffset      = by * BLOCKSIZE * size + bx * BLOCKSIZE;
        int blockPivElOffset = by * BLOCKSIZE * size + diagEl;
        int blockPivOffset   = diagEl * size + bx * BLOCKSIZE;

         inBlock[ty][tx] = dMatrixIn[blockOffset + ty * size + tx];
           pivEl[ty][tx] = dMatrixIn[blockPivElOffset + ty * size + tx];
        pivBlock[ty][tx] = dMatrixIn[blockPivOffset + ty * size + tx];
        
        __syncthreads ();

        for (int i = BLOCKSIZEMINUS1; i >= 0; i--)
        {
            inBlock[ty][tx] = cuCsub(inBlock[ty][tx], cuCmul(pivEl[ty][i],pivBlock[i][tx]));
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
            inBlock[ty][tx] = cuCsub(inBlock[ty][tx], cuCmul(pivEl[ty][i],pivBlock[i][tx]));
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
            inBlock[ty][tx] = cuCsub(inBlock[ty][tx], cuCmul(pivEl[ty][i],pivBlock[i][tx]));
        }

        __syncthreads ();

        dMatrixInv[blockOffset + ty * size + tx] = inBlock[ty][tx];
          
    }
}


//************************************************
// Kernel normalizeDiag
//************************************************

__global__ void normalizeDiag(int BLOCKSIZE,cuDoubleComplex *diagMatrix, cuDoubleComplex *invMatrix, int size,int row)
{
          
    int tx = threadIdx.x; int ty = threadIdx.y; int bx = blockIdx.x;

    int blockOffset = bx * BLOCKSIZE;

    __shared__ cuDoubleComplex diagEl[MAXSIZE];

    if (tx == ty )
    {
        diagEl[ty] = diagMatrix[row + ty * size + tx];
    }
    __syncthreads ();

    invMatrix[blockOffset + ty * size + tx] = cuCdiv(
    invMatrix[blockOffset + ty * size + tx],diagEl[ty]);
    
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
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************


extern "C" void cuda_complex_invert_(int *pBlockSize, cuDoubleComplex** A_, cuDoubleComplex** invA_, int *pndim)
{
    int i; int size = *pndim; int BlockSize = *pBlockSize;

    dim3 idyThreads (BlockSize);    
    dim3 idyBlocks ( size / BlockSize );
    dim3 nThreads (BlockSize, BlockSize);   
    dim3 nBlocks ( size / BlockSize );        
    dim3 nBlocksRest ( size / BlockSize, size / BlockSize);

    unsigned int mat_size=size*size*sizeof(cuDoubleComplex);

    int pinned = 1 ; if(size<3600) pinned=0;

//  printf(" copy memory host to gpu \n ");

    cuDoubleComplex *invA ; cuDoubleComplex *A  ; cuDoubleComplex *invAd ; cuDoubleComplex *Ad;


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



    GPUsetIdentity <<< idyBlocks, idyThreads >>> (BlockSize,invA, size);


    cudaThreadSynchronize ();

    //calculate the right diagonal Matrix (L)

    for (i = 0; i < size; i += BlockSize)
    {
        int offset = i * size + i;

        // *  step 1:
        // *  calculate the triangle matrix
        // *  store the pivot elements to left part of the triangel

        eliminateBlock<<< 1, nThreads >>> ('L',BlockSize,A + offset, size);
        
        cudaThreadSynchronize ();
    
        // * step 2:
        // *  calculate the rest of the rows with the pivot elements from step 1
        
        adjustRowL<<< nBlocks, nThreads >>> (BlockSize,A + i * size, A + offset, invA + i * size, size, i);
        
        cudaThreadSynchronize ();
    
        //* step 3:
        //* Fill the colls below the block with the pivot elements they are used
        //*    to get the colls to zero and multiply with the row
        
        eliminateColL<<< nBlocks, nThreads >>> (BlockSize,A + i, size, i);
        
        cudaThreadSynchronize ();

        //* step 4:
        //*  Adjust the rest of the Matrix with the calculated pivot Elements
        //*  El_new_0 -= (p0+p1+p2..+p15) * El_piv_0
        
        eliminateRestL<<< nBlocksRest, nThreads >>> (BlockSize,A, invA, size, i);
        
        cudaThreadSynchronize ();
    }
    


   //Set the left lower diagonalmatrix to zero (async?)

    for (i = 1; i < size; i++)
    {
      int offset = i * size;
      if(pinned==0){
         cudaMemset ((void *) (A  + offset), 0, i*sizeof(cuDoubleComplex));}
      else{
         memset ((void *) (Ad  + offset), 0, i*sizeof(cuDoubleComplex));
      };
      cudaEventSynchronize(0);
    }
    cudaThreadSynchronize ();




    //calculate the right diagonal Matrix (U)
    
    for (i = (size - BlockSize); i >= 0; i -= BlockSize)
    {
        int offset = i * size + i;

        //*  step 1:
        //*  calculate the triangle matrix
        //*  store the pivot elements to left part of the triangel
        
        eliminateBlock<<< 1, nThreads >>> ('U',BlockSize,A + offset, size);
        
        cudaThreadSynchronize ();

        //* step 2:
        //*  calculate the rest of the rows with the pivot elements from step 1
        
        adjustRowU<<< nBlocks, nThreads >>> (BlockSize,A + offset,invA + i*size, size, i);
        
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

extern "C" void cuda_c_invert_(int *pBlockSize, double** A_, double** A_i, double** invA_, double** invA_i, int *pndim)
{

  int i; int nlines=*pndim; int ncolumns=*pndim; int ntot=nlines*ncolumns ;
  cuDoubleComplex **aa; cuDoubleComplex **inva;
  aa   = (cuDoubleComplex**) malloc(sizeof(cuDoubleComplex)*ntot);
  inva = (cuDoubleComplex**) malloc(sizeof(cuDoubleComplex)*ntot);

  for(i=0; i<ntot; i++)
   {
  *aa[i]=make_cuDoubleComplex(*A_[i],*A_i[i])
  ;};

  int k = nlines/2  ;
  cuda_complex_invert_(&k, aa, inva, &nlines);

  for(i=0; i<ntot; i++)
   {
  *invA_[i]  =cuCreal(*inva[i]);
  *invA_i[i] =cuCimag(*inva[i]);
  };

  free(inva);free(aa);

}

//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************

