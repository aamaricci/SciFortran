#define MAX_SIZE 16

#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <cutil.h>

//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************

__global__ void Muld_gpu_fortran(int BLOCKSIZE, double* A, double* B, double* C, int ha, int wA, int wB)
{

// A: ha x wA
// B: wA x wB
// C: ha x wB

int bx = blockIdx.x; int by = blockIdx.y; int tx = threadIdx.x; int ty = threadIdx.y;


int aBegin  =   BLOCKSIZE * by                 ;
int aEnd    =   BLOCKSIZE * by  + (wA - 1)*ha  ;
int bBegin  =   BLOCKSIZE * bx *   wA          ;

int aStep   =   BLOCKSIZE * ha                 ;
int bStep   =   BLOCKSIZE                      ;

int c       =   BLOCKSIZE * by + BLOCKSIZE * bx * ha  ;

double Csub = 0;

/*******************************************************************************/
for (int a = aBegin, b = bBegin; a <= aEnd; a += aStep, b += bStep) 
{
__shared__ double As[MAX_SIZE][MAX_SIZE]; __shared__ double Bs[MAX_SIZE][MAX_SIZE];
   As[ty][tx] = A[a    +      ty + tx*ha];   Bs[ty][tx] = B[b   +      ty + tx*wA ];
__syncthreads();
  for (int k = 0; k < BLOCKSIZE; ++k) Csub += As[ty][k] * Bs[k][tx];
__syncthreads();
}
   C[c  +      ty + tx*ha ] = Csub ;
/*******************************************************************************/

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

__global__ void Muld_gpu_cxx(int BLOCKSIZE, double* A, double* B, double* C, int ha, int wA, int wB)
{

// A: ha x wA
// B: wA x wB
// C: ha x wB

int bx = blockIdx.x; int by = blockIdx.y; int tx = threadIdx.x; int ty = threadIdx.y;


int aBegin  =  wA * BLOCKSIZE * by           ;
int aEnd    =  wA * BLOCKSIZE * by  + wA - 1 ;
int bBegin  =       BLOCKSIZE * bx           ;

int aStep   =       BLOCKSIZE                ;
int bStep   =       BLOCKSIZE * wB           ;

int c       =  wB * BLOCKSIZE * by + BLOCKSIZE * bx ;

double Csub = 0;

/*******************************************************************************/
for (int a = aBegin, b = bBegin; a <= aEnd; a += aStep, b += bStep) 
{
__shared__ double As[MAX_SIZE][MAX_SIZE]; __shared__ double Bs[MAX_SIZE][MAX_SIZE];
 As[ty][tx] = A[a    + wA * ty + tx   ];   Bs[ty][tx] = B[b   + wB * ty + tx    ];
__syncthreads();
for (int k = 0; k < BLOCKSIZE; ++k) Csub += As[ty][k] * Bs[k][tx];
__syncthreads();
}
 C[c  + wB * ty + tx    ] = Csub ;
/*******************************************************************************/

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

extern "C" void matmul_gpu_fortran_(int* pBLOCKSIZE, double** A, double** B, double** C,  int* hA, int* wA, int* wB)
{

// A: ha x wA // B: wA x wB // C: ha x wB

int verbose   =  0 ;
int pinned    =  1 ;
if(*hA<5000 && *wA<5000 && *wB<5000) pinned=0;

int BLOCKSIZE = *pBLOCKSIZE; 

int sizea; int sizeb; int sizec; 
double *Ad  ; double *Bd  ; double *Cd  ;
double *Add ; double *Bdd ; double *Cdd ;
sizea = *hA * *wA * sizeof(double);
sizeb = *wA * *wB * sizeof(double);
sizec = *hA * *wB * sizeof(double);

if(verbose){printf(" allocate \n ");};
if(pinned==0){
 if(verbose) printf(" allocate memory on device \n");
  cudaMalloc((void**) &Ad, sizea); 
  cudaMalloc((void**) &Bd, sizeb); 
  cudaMalloc((void**) &Cd, sizec);
}else{

  if(verbose) printf( " allocate pinned memory \n ");
  cudaSetDevice(0); cudaSetDeviceFlags( cudaDeviceMapHost );
  cudaHostAlloc((void**)&Ad, sizea, cudaHostAllocMapped | cudaHostAllocPortable );
  cudaHostAlloc((void**)&Bd, sizeb, cudaHostAllocMapped | cudaHostAllocPortable );
  cudaHostAlloc((void**)&Cd, sizec, cudaHostAllocMapped | cudaHostAllocPortable );
  cudaHostGetDevicePointer((void**)&Add, Ad, 0 );
  cudaHostGetDevicePointer((void**)&Bdd, Bd, 0 );
  cudaHostGetDevicePointer((void**)&Cdd, Cd, 0 );
}

 if(verbose){printf(" copy inputs \n ");};
 if(pinned==0){
   cudaMemcpy(Bd, B, sizeb, cudaMemcpyHostToDevice);
   cudaMemcpy(Ad, A, sizea, cudaMemcpyHostToDevice);
 }else{
   cudaMemcpy(Bd, B, sizeb, cudaMemcpyHostToHost);
   cudaMemcpy(Ad, A, sizea, cudaMemcpyHostToHost);
 };

 dim3 dimBlock(BLOCKSIZE, BLOCKSIZE);
 dim3 dimGrid(*wB / dimBlock.x, *hA / dimBlock.y);
 if(verbose) printf(" start device kernel \n ");

 if(pinned==0){
 Muld_gpu_fortran<<<dimGrid, dimBlock>>>(BLOCKSIZE, Ad, Bd, Cd, *hA, *wA, *wB);
}else{
 Muld_gpu_fortran<<<dimGrid, dimBlock>>>(BLOCKSIZE, Add , Bdd , Cdd ,*hA, *wA, *wB );
}

 if(verbose){printf(" copy results \n ");};
 if(pinned==0){
 cudaMemcpy(C,Cd,sizec,cudaMemcpyDeviceToHost); 
 }else{
 cudaEventSynchronize(0);
 cudaMemcpy(C,Cd,sizec,cudaMemcpyHostToHost);
 };

 if(verbose){printf(" free memory \n ");};
 if(pinned==1){
  cudaFreeHost(Ad); cudaFreeHost(Bd); cudaFreeHost(Cd);
 }else{
  cudaFree(Ad); cudaFree(Bd); cudaFree(Cd);
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

extern "C" void matmul_gpu_cxx(int* pBLOCKSIZE, double** A, double** B, double** C,  int* hA, int* wA, int* wB)
{

// A: ha x wA
// B: wA x wB
// C: ha x wB

int verbose   =  0 ;
int pinned    =  1 ;
if(*hA<5000 && *wA<5000 && *wB<5000) pinned=0;
int BLOCKSIZE=*pBLOCKSIZE;

int sizea; int sizeb; int sizec;
double *Ad  ; double *Bd  ; double *Cd  ;
double *Add ; double *Bdd ; double *Cdd ;
sizea = *hA * *wA * sizeof(double);
sizeb = *wA * *wB * sizeof(double);
sizec = *hA * *wB * sizeof(double);

if(verbose){printf(" allocate \n ");};
if(pinned==0){
 if(verbose) printf(" allocate memory on device \n");
  cudaMalloc((void**) &Ad, sizea);
  cudaMalloc((void**) &Bd, sizeb);
  cudaMalloc((void**) &Cd, sizec);
}else{

 if(verbose) printf( " allocate pinned memory \n ");
  cudaSetDevice(0); cudaSetDeviceFlags( cudaDeviceMapHost );
  cudaHostAlloc((void**)&Ad, sizea, cudaHostAllocMapped | cudaHostAllocPortable );
  cudaHostAlloc((void**)&Bd, sizeb, cudaHostAllocMapped | cudaHostAllocPortable );
  cudaHostAlloc((void**)&Cd, sizec, cudaHostAllocMapped | cudaHostAllocPortable );
  cudaHostGetDevicePointer((void**)&Add, Ad, 0 );
  cudaHostGetDevicePointer((void**)&Bdd, Bd, 0 );
  cudaHostGetDevicePointer((void**)&Cdd, Cd, 0 );
}

 if(verbose){printf(" copy inputs \n ");};
 if(pinned==0){
   cudaMemcpy(Bd, B, sizeb, cudaMemcpyHostToDevice);
   cudaMemcpy(Ad, A, sizea, cudaMemcpyHostToDevice);
 }else{
   cudaMemcpy(Bd, B, sizeb, cudaMemcpyHostToHost);
   cudaMemcpy(Ad, A, sizea, cudaMemcpyHostToHost);
 };

 dim3 dimBlock(BLOCKSIZE, BLOCKSIZE);
 dim3 dimGrid(*wB / dimBlock.x, *hA / dimBlock.y);

 if(pinned==0){
 Muld_gpu_cxx<<<dimGrid, dimBlock>>>(BLOCKSIZE, Ad, Bd, Cd, *hA, *wA, *wB);
}else{
 Muld_gpu_cxx<<<dimGrid, dimBlock>>>(BLOCKSIZE, Add , Bdd , Cdd ,*hA, *wA, *wB );
}

 if(verbose){printf(" copy results \n ");};
 if(pinned==0){
 cudaMemcpy(C,Cd,sizec,cudaMemcpyDeviceToHost);
 }else{
 cudaEventSynchronize(0);
 cudaMemcpy(C,Cd,sizec,cudaMemcpyHostToHost);
 };

 if(verbose){printf(" free memory \n ");};
 if(pinned==1){
  cudaFreeHost(Ad); cudaFreeHost(Bd); cudaFreeHost(Cd);
 }else{
  cudaFree(Ad); cudaFree(Bd); cudaFree(Cd);
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

