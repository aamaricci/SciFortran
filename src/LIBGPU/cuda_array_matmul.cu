 
 #include <stdio.h>
 #include <stdlib.h>
 #include <cuda_runtime.h>
 #include <cublas.h>
 #include <cuComplex.h>

#define MAX_BLOCK 22

cuDoubleComplex *collect_array;


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

texture<int2,1,cudaReadModeElementType> tex;
inline void   bind_x(cuDoubleComplex *x) {    cudaBindTexture(0,tex,x); };
inline void unbind_x()                   {  cudaUnbindTexture(  tex  ); };
__inline__  __device__ cuDoubleComplex fetch_x(const int& i)
  {
         int  jj = 2*(i-1);
         int2 v  = tex1Dfetch(tex,jj);
      double rr  = __hiloint2double(v.y, v.x);
              v  = tex1Dfetch(tex,jj+1);
      double im  = __hiloint2double(v.y, v.x);
      return make_cuDoubleComplex(rr,im);
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
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************

__global__ void matmul_array_cuda( int nfrequ, int nnn, cuDoubleComplex *collect )


{

    cuDoubleComplex zero; zero=make_cuDoubleComplex(0.0,0.0); 
   
    __shared__ cuDoubleComplex shared_[256];
               cuDoubleComplex E[MAX_BLOCK*MAX_BLOCK];
    int ifrequ;

    for(int iorb=0;iorb<nnn*nnn;iorb++) shared_[iorb]=zero; ifrequ=threadIdx.x;

    if(ifrequ<nfrequ)
    { 
       for(int iorb=0;iorb<nnn*nnn;iorb++) E[iorb] = collect[iorb+nnn*nnn*ifrequ]; 
    };

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
//********************************************
//********************************************

extern "C" void complex_array_of_matrices_matmul(int* nnn_, int* nfrequ_, cuDoubleComplex* collect_ )

{

  int nnn ; int nfrequ; nnn=*nnn_; nfrequ=*nfrequ_;

  if(nnn>MAX_BLOCK)   { printf( " sum of inverse complex cuda : matrices are too big!!!!! \n"); };
  if(nfrequ>512)      { printf( " too many blocks in inverse array cuda \n ");};

  cudaMalloc (  (void**)&collect_array  , nfrequ*nnn*nnn*sizeof(cuDoubleComplex) ); 
  cudaMemcpy ( collect_array ,collect_  , nnn*nnn*nfrequ*sizeof(cuDoubleComplex) , cudaMemcpyHostToDevice);

  matmul_array_cuda <<<1,512>>> ( nfrequ,nnn,collect_array);

  cudaEventSynchronize(0); cudaThreadSynchronize();
  cudaMemcpy (collect_ ,collect_array  , nfrequ*nnn*nnn*sizeof(cuDoubleComplex) , cudaMemcpyDeviceToHost);
  cudaFree(collect_array);


}

//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
