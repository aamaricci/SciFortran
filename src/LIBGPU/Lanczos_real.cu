
 #include <stdio.h>
 #include <stdlib.h>
 #include <cuda_runtime.h>
 #include <cublas.h>

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

const int use_texture_lr = 1 ;

#define WARP_SIZE 32
#define MAX_BLOCK 65500

texture<int2,1,cudaReadModeElementType> tex;

inline void   bind_x(double *x, int N) {    cudaBindTexture(0,tex,x,N*sizeof(double)); };
inline void unbind_x()                 {  cudaUnbindTexture(  tex  ); };

__inline__    __device__ double fetch_x(const int& i)
  {   int2 v = tex1Dfetch(tex,i-1); return __hiloint2double(v.y, v.x); }

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

  //************************************************//
  //               Kernel Hmult                     //
  //************************************************//

     //---------------------------------------//
     //---------------------------------------//
     //---------------------------------------//
     //---------------------------------------//

  __global__  void Hmult_ker(int ngrid,int BLOCK_SIZE, int num_rows, double *y, const double *x, const double *QUART, const double *diagsz, 
                             const int *noffsz, const int *Aj, const double *Ax, const int *offdia, const int use_texture_lr)
{
   __shared__ double sdata[16][WARP_SIZE]; 
   __shared__ int     ptrs[32][2];
   __shared__ double  temp[32][3];

          int row ;  
    const int warp_lane   = threadIdx.y; 
    const int thread_lane = threadIdx.x;

    int row_start; int row_end; int jj;

   for(int iii=0; iii<=ngrid; iii++){

   row = BLOCK_SIZE * (blockIdx.y+iii*MAX_BLOCK) + threadIdx.y ;

    if(row<num_rows)
   {
        if(thread_lane==0) ptrs[warp_lane][0]=offdia[row];
        if(thread_lane==1) ptrs[warp_lane][1]=noffsz[row];
        if(thread_lane==2) temp[warp_lane][0]=QUART[row];
        if(thread_lane==3) temp[warp_lane][1]=diagsz[row];
        
         if(use_texture_lr==0)
         {
           y[row]=(temp[warp_lane][0]+temp[warp_lane][1]) * x[row];
         }else{
           y[row]=(temp[warp_lane][0]+temp[warp_lane][1]) * fetch_x(row+1);
         };

        row_start = ptrs[warp_lane][0] ; row_end = ptrs[warp_lane][1]+row_start ; 

        sdata[threadIdx.y][threadIdx.x]=0.0;
        
       if(use_texture_lr==1)
       { 
        for(jj=row_start+thread_lane;jj<row_end;jj+=WARP_SIZE) sdata[threadIdx.y][threadIdx.x]+=Ax[jj] * fetch_x(Aj[jj]); 
       }else{
        for(jj=row_start+thread_lane;jj<row_end;jj+=WARP_SIZE) sdata[threadIdx.y][threadIdx.x]+=Ax[jj] *     x[Aj[jj]-1]; 
       }

        if (thread_lane < 16) { sdata[threadIdx.y][threadIdx.x] += sdata[threadIdx.y][threadIdx.x + 16]; };
        if (thread_lane <  8) { sdata[threadIdx.y][threadIdx.x] += sdata[threadIdx.y][threadIdx.x +  8]; };
        if (thread_lane <  4) { sdata[threadIdx.y][threadIdx.x] += sdata[threadIdx.y][threadIdx.x +  4]; };
        if (thread_lane <  2) { sdata[threadIdx.y][threadIdx.x] += sdata[threadIdx.y][threadIdx.x +  2]; };
        if (thread_lane <  1) { sdata[threadIdx.y][threadIdx.x] += sdata[threadIdx.y][threadIdx.x +  1]; };

        if (thread_lane == 0) y[row] += sdata[threadIdx.y][threadIdx.x];
   };
   };
}
     //---------------------------------------//
     //---------------------------------------//
     //---------------------------------------//
     //---------------------------------------//
     //---------------------------------------//
     //---------------------------------------//
     //---------------------------------------//

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

extern "C" void hmult_sz_real_cuda_rout_(int *pblocksize, int *offdiasize, int *roff, 
int* ntot, double* QUART, double* diagsz, double* vec_in, double* vec_out, int* noffsz, 
int* rankoffsz, double* offdiagsz )
{

 //---------------------------------------------------------------------------------------//
     int verbose=0;
 //---------------------------------------------------------------------------------------//

     if(verbose==1) printf(" start Hmult GPU \n" );

     int blocksize= *pblocksize; int size = *ntot;
 
    int nb=(size-size % blocksize)/blocksize+1; int ngrid=nb/MAX_BLOCK; if(ngrid>0) nb=MAX_BLOCK;
    dim3 bl(1,nb),th(WARP_SIZE,blocksize);

    if(verbose==1) printf( " --------------- \n  Nblock=%d Ngrid=%d \n ----------------- \n ",nb,ngrid);

     double  *vec_in_gpu,*vec_out_gpu,*QUART_gpu,*diagsz_gpu,*offdiagsz_gpu;
     int     *noffsz_gpu,*rankoffsz_gpu;
     double  *vec_in_gpu_pointer,*vec_out_gpu_pointer;
     double  *QUART_gpu_pointer,*diagsz_gpu_pointer,*offdiagsz_gpu_pointer;
     int     *noffsz_gpu_pointer,*rankoffsz_gpu_pointer;
     int     *offdia_gpu, *offdia_gpu_pointer;

     if(verbose==1) printf(" GPU , size of Lanczos vector = %d ", size);

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

       cudaSetDevice(0); cudaSetDeviceFlags( cudaDeviceMapHost ); cudaEventSynchronize(0);

       if(use_texture_lr==1){
         cudaMalloc((void**)&vec_in_gpu_pointer,sizeof(double)*size);
         bind_x(vec_in_gpu_pointer,size);
       }else{
         cudaHostAlloc((void**)&vec_in_gpu    , size*sizeof(double),        cudaHostAllocMapped | cudaHostAllocPortable );
         cudaHostGetDevicePointer((void**)  &vec_in_gpu_pointer  ,  vec_in_gpu    , 0 );
       }

       cudaHostAlloc((void**)&vec_out_gpu   , size*sizeof(double),        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&QUART_gpu     , size*sizeof(double),        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&diagsz_gpu    , size*sizeof(double),        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&offdiagsz_gpu , *offdiasize*sizeof(double), cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&noffsz_gpu    , size*sizeof(int)   ,        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&rankoffsz_gpu , *roff*sizeof(int)  ,        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&offdia_gpu    , size*sizeof(int)   ,        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaEventSynchronize(0);

       cudaHostGetDevicePointer((void**)   &vec_out_gpu_pointer ,  vec_out_gpu    , 0 );
       cudaHostGetDevicePointer((void**)   &QUART_gpu_pointer ,    QUART_gpu      , 0 );
       cudaHostGetDevicePointer((void**)   &diagsz_gpu_pointer ,   diagsz_gpu     , 0 );
       cudaHostGetDevicePointer((void**)   &offdiagsz_gpu_pointer, offdiagsz_gpu  , 0 );
       cudaHostGetDevicePointer((void**)   &noffsz_gpu_pointer,    noffsz_gpu     , 0 );
       cudaHostGetDevicePointer((void**)   &rankoffsz_gpu_pointer, rankoffsz_gpu  , 0 );
       cudaHostGetDevicePointer((void**)   &offdia_gpu_pointer,    offdia_gpu     , 0 );
       cudaEventSynchronize(0);

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//
       if(use_texture_lr==0) {
        cudaMemcpy(vec_in_gpu,    vec_in,      size*sizeof(double),        cudaMemcpyHostToHost);
       }else{
        cudaMemcpy(vec_in_gpu_pointer,    vec_in,      size*sizeof(double),        cudaMemcpyHostToDevice);
       }
       cudaMemcpy(vec_out_gpu,   vec_out,     size*sizeof(double),        cudaMemcpyHostToHost);
       cudaMemcpy(QUART_gpu,     QUART,       size*sizeof(double),        cudaMemcpyHostToHost);
       cudaMemcpy(diagsz_gpu,    diagsz,      size*sizeof(double),        cudaMemcpyHostToHost);
       cudaMemcpy(offdiagsz_gpu, offdiagsz,   *offdiasize*sizeof(double), cudaMemcpyHostToHost);
       cudaMemcpy(noffsz_gpu,    noffsz ,     size*sizeof(int),           cudaMemcpyHostToHost);
       cudaMemcpy(rankoffsz_gpu, rankoffsz,   *roff*sizeof(int),          cudaMemcpyHostToHost);
       cudaEventSynchronize(0);
       offdia_gpu[0]=0; for(int istate=1; istate<size; istate++) { offdia_gpu[istate]=offdia_gpu[istate-1]+noffsz[istate-1]; };
       cudaEventSynchronize(0);
 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

  if(verbose==1) printf(" call kernel  \n ");
  Hmult_ker<<<bl,th>>>(ngrid,blocksize,size,vec_out_gpu_pointer,vec_in_gpu_pointer,QUART_gpu_pointer,
  diagsz_gpu_pointer,noffsz_gpu_pointer,rankoffsz_gpu_pointer,offdiagsz_gpu_pointer,offdia_gpu_pointer,use_texture_lr); 
  cudaEventSynchronize(0); cudaThreadSynchronize();
  if(verbose==1) printf(" .....done.....  \n ");

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

     cudaMemcpy(vec_out,vec_out_gpu,size*sizeof(double),cudaMemcpyHostToHost);
     cudaEventSynchronize(0);
     cudaFreeHost(vec_out_gpu);
     if(use_texture_lr==0){
      cudaFreeHost(vec_in_gpu);}
     else{
      unbind_x();
      cudaFree(vec_in_gpu_pointer);
     }
     cudaFreeHost(offdia_gpu);
     cudaFreeHost(QUART_gpu);
     cudaFreeHost(diagsz_gpu);
     cudaFreeHost(rankoffsz_gpu);
     cudaFreeHost(offdiagsz_gpu);
     cudaFreeHost(noffsz_gpu);
     cudaEventSynchronize(0);

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

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
//********************************************

 void one_step_lanczos_cuda(int blocksize, int Niter, int size, int iter, double *diag, 
 double *subdiag, double *vec_tmp_gpu, double *vec_in_gpu, double *vec_out_gpu, double *vec_tmp_gpu_pointer, 
 double *vec_in_gpu_pointer, double *vec_out_gpu_pointer, double *QUART_gpu_pointer, double *diagsz_gpu_pointer, 
 int *noffsz_gpu_pointer, int *rankoffsz_gpu_pointer, double *offdiagsz_gpu_pointer, int *offdia_gpu_pointer, 
 double *QUART_gpu, double *diagsz_gpu, int *noffsz_gpu, int *rankoffsz_gpu, double *offdiagsz_gpu, int *offdia_gpu)

{
    int verbose=0 ; int psize=size; 

    int nb=(size-size % blocksize)/blocksize+1 ; int ngrid=nb/MAX_BLOCK; if(ngrid>0)nb=MAX_BLOCK; dim3 bl(1,nb),th(WARP_SIZE,blocksize);

    if(verbose==1) printf( " \n --------------- \n  Nblock=%d \n Ngrid=%d \n blocksize=%d \n ----------------- \n ",nb,ngrid,blocksize);

    if(verbose==1) printf ( " Sdot, vec in norm \n ");
    double normv = sqrt(cublasDdot(size,vec_in_gpu_pointer,1,vec_in_gpu_pointer,1)); cudaEventSynchronize(0); cudaThreadSynchronize();
    if(verbose==1) printf( " norm=%f \n ", normv);

    cublasDscal(size,1.0/normv,vec_in_gpu_pointer,1); cudaEventSynchronize(0); cudaThreadSynchronize(); 

    if(verbose==1) printf( " call kernel ... \n ");
    Hmult_ker<<<bl,th>>>(ngrid,blocksize,psize,vec_out_gpu_pointer,vec_in_gpu_pointer,QUART_gpu_pointer,diagsz_gpu_pointer,
                         noffsz_gpu_pointer,rankoffsz_gpu_pointer,offdiagsz_gpu_pointer,offdia_gpu_pointer,use_texture_lr); 

   if(verbose==1) printf( " done.... \n " );
   cudaEventSynchronize(0); cudaThreadSynchronize();

   if(iter>0){cublasDaxpy(size,-subdiag[iter],vec_tmp_gpu_pointer,1,vec_out_gpu_pointer,1);}; cudaEventSynchronize(0); cudaThreadSynchronize();

  if(use_texture_lr==0){
   cudaMemcpy(vec_tmp_gpu,vec_in_gpu,size*sizeof(double),cudaMemcpyHostToHost);cudaEventSynchronize(0);cudaThreadSynchronize();
  }else{
   cudaMemcpy(vec_tmp_gpu,vec_in_gpu_pointer,size*sizeof(double),cudaMemcpyDeviceToHost);cudaEventSynchronize(0);cudaThreadSynchronize();
  }

   diag[iter]=cublasDdot(size, vec_out_gpu_pointer,1,vec_in_gpu_pointer,1);cudaEventSynchronize(0);cudaThreadSynchronize();

   cublasDaxpy(size,-diag[iter],vec_tmp_gpu_pointer,1,vec_out_gpu_pointer,1); cudaEventSynchronize(0); cudaThreadSynchronize();
   normv = sqrt(cublasDdot(size, vec_out_gpu_pointer,1,vec_out_gpu_pointer,1)); cudaEventSynchronize(0); cudaThreadSynchronize();

   if(iter<Niter-1) subdiag[iter+1]=normv;

  if(use_texture_lr==0){
   cudaMemcpy( vec_in_gpu, vec_out_gpu, size*sizeof(double), cudaMemcpyHostToHost); cudaEventSynchronize(0); cudaThreadSynchronize();
  }else{
   cudaMemcpy( vec_in_gpu_pointer, vec_out_gpu, size*sizeof(double), cudaMemcpyHostToDevice); cudaEventSynchronize(0); cudaThreadSynchronize();
  }
 
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

extern "C" void lanczos_real_dynamic_cuda_(int *pblocksize, int *Niter_lanczos_,int *offdiasize, int *roff, int *ntot, 
double *QUART, double *diagsz, int *noffsz, int *rankoffsz, double *offdiagsz, double *diag, double *subdiag , 
double *vecinit)
{

 //---------------------------------------------------------------------------------------//
     int verbose=0;
 //---------------------------------------------------------------------------------------//

     int Niter_lanczos= *Niter_lanczos_; 

     if(verbose==1) printf(" start Lanczos Real on GPU \n" );

     int size = *ntot; int blocksize= *pblocksize; 

     double   *vec_in_gpu,*vec_out_gpu,*QUART_gpu,*diagsz_gpu,*offdiagsz_gpu,*vec_tmp_gpu;
     int      *noffsz_gpu,*rankoffsz_gpu;
     double   *vec_in_gpu_pointer,*vec_out_gpu_pointer,*vec_tmp_gpu_pointer;
     double   *QUART_gpu_pointer,*diagsz_gpu_pointer,*offdiagsz_gpu_pointer;
     int      *noffsz_gpu_pointer,*rankoffsz_gpu_pointer;
     int      *offdia_gpu, *offdia_gpu_pointer;

     if(verbose==1) printf(" GPU get eigenvalues, size of Lanczos vectors = %d ", size);

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

       cudaSetDevice(0); cudaSetDeviceFlags( cudaDeviceMapHost ); cudaEventSynchronize(0);

       cudaHostAlloc((void**)&vec_tmp_gpu   , size*sizeof(double),        cudaHostAllocMapped | cudaHostAllocPortable );
  
       if(use_texture_lr==1){
         cudaMalloc((void**)&vec_in_gpu_pointer,sizeof(double)*size);
         bind_x(vec_in_gpu_pointer,size);
       }else{
         cudaHostAlloc((void**)&vec_in_gpu    , size*sizeof(double),        cudaHostAllocMapped | cudaHostAllocPortable );
         cudaHostGetDevicePointer((void**)  &vec_in_gpu_pointer  ,  vec_in_gpu    , 0 );
       }

       cudaHostAlloc((void**)&vec_out_gpu   , size*sizeof(double),        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&QUART_gpu     , size*sizeof(double),        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&diagsz_gpu    , size*sizeof(double),        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&offdiagsz_gpu , *offdiasize*sizeof(double), cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&noffsz_gpu    , size*sizeof(int)   ,        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&rankoffsz_gpu , *roff*sizeof(int)  ,        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&offdia_gpu    , size*sizeof(int)   ,        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaEventSynchronize(0);

       cudaHostGetDevicePointer((void**)  &vec_tmp_gpu_pointer ,  vec_tmp_gpu   , 0 );
       cudaHostGetDevicePointer((void**)  &vec_out_gpu_pointer ,  vec_out_gpu   , 0 );
       cudaHostGetDevicePointer((void**)  &QUART_gpu_pointer ,    QUART_gpu     , 0 );
       cudaHostGetDevicePointer((void**)  &diagsz_gpu_pointer ,   diagsz_gpu    , 0 );
       cudaHostGetDevicePointer((void**)  &offdiagsz_gpu_pointer, offdiagsz_gpu , 0 );
       cudaHostGetDevicePointer((void**)  &noffsz_gpu_pointer,    noffsz_gpu    , 0 );
       cudaHostGetDevicePointer((void**)  &rankoffsz_gpu_pointer, rankoffsz_gpu , 0 );
       cudaHostGetDevicePointer((void**)  &offdia_gpu_pointer,    offdia_gpu    , 0 );
       cudaEventSynchronize(0);

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//
       cudaMemcpy(QUART_gpu,     QUART,       size*sizeof(double),        cudaMemcpyHostToHost);
       cudaMemcpy(diagsz_gpu,    diagsz,      size*sizeof(double),        cudaMemcpyHostToHost);
       cudaMemcpy(offdiagsz_gpu, offdiagsz,   *offdiasize*sizeof(double), cudaMemcpyHostToHost);
       cudaMemcpy(noffsz_gpu,    noffsz ,     size*sizeof(int),           cudaMemcpyHostToHost);
       cudaMemcpy(rankoffsz_gpu, rankoffsz,   *roff*sizeof(int),          cudaMemcpyHostToHost);
       cudaEventSynchronize(0);
       offdia_gpu[0]=0; for(int istate=1; istate<size; istate++) { offdia_gpu[istate]=offdia_gpu[istate-1]+noffsz[istate-1]; };
       cudaEventSynchronize(0);
 //---------------------------------------------------------------------------------------//

  if(verbose==1) printf( " initialized, now run actual Lanczos \n " );

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

  cudaThreadSynchronize(); cudaEventSynchronize(0); 

 if(use_texture_lr==0){
   cudaMemcpy(vec_in_gpu,vecinit,size*sizeof(double),cudaMemcpyHostToHost);
 }else{
   cudaMemcpy(vec_out_gpu,vecinit,size*sizeof(double),cudaMemcpyHostToHost);
   cudaMemcpy(vec_in_gpu_pointer,vec_out_gpu,size*sizeof(double),cudaMemcpyHostToDevice);
 }

  cudaThreadSynchronize(); cudaEventSynchronize(0);

  for(int iter=0;iter<Niter_lanczos;iter++){
   if(verbose==1) printf( " Lanczos iteration %d / %d \n", iter,Niter_lanczos);
   one_step_lanczos_cuda(blocksize,Niter_lanczos,size,iter,diag,subdiag,vec_tmp_gpu,vec_in_gpu,vec_out_gpu,vec_tmp_gpu_pointer,vec_in_gpu_pointer,vec_out_gpu_pointer,QUART_gpu_pointer,diagsz_gpu_pointer,noffsz_gpu_pointer,rankoffsz_gpu_pointer,offdiagsz_gpu_pointer,offdia_gpu_pointer,QUART_gpu,diagsz_gpu,noffsz_gpu,rankoffsz_gpu, offdiagsz_gpu, offdia_gpu);

  };

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

     cudaEventSynchronize(0);
     cudaFreeHost(vec_tmp_gpu);
     if(use_texture_lr==0){
      cudaFreeHost(vec_in_gpu);}
     else{
      unbind_x();
      cudaFree(vec_in_gpu_pointer);
     }
     cudaFreeHost(vec_out_gpu);
     cudaFreeHost(offdia_gpu);
     cudaFreeHost(QUART_gpu);
     cudaFreeHost(diagsz_gpu);
     cudaFreeHost(rankoffsz_gpu);
     cudaFreeHost(offdiagsz_gpu);
     cudaFreeHost(noffsz_gpu);
     cudaEventSynchronize(0);

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

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
//********************************************

extern "C" void lanczos_real_cuda_(int *pblocksize,  int *Niter_lanczos_,int *offdiasize, int *roff, int *ntot, 
double *QUART, double *diagsz, int *noffsz, int *rankoffsz, double *offdiagsz, double *diag, double *subdiag )
{

 //---------------------------------------------------------------------------------------//
     int verbose=0;
 //---------------------------------------------------------------------------------------//

     int Niter_lanczos= *Niter_lanczos_; 
     if(verbose==1) printf(" start Lanczos Real on GPU \n" );
     int size = *ntot; int blocksize = *pblocksize;

     double   *vec_in_gpu,*vec_out_gpu,*QUART_gpu,*diagsz_gpu,*offdiagsz_gpu,*vec_tmp_gpu;
     int      *noffsz_gpu,*rankoffsz_gpu;
     double   *vec_in_gpu_pointer,*vec_out_gpu_pointer,*vec_tmp_gpu_pointer;
     double   *QUART_gpu_pointer,*diagsz_gpu_pointer,*offdiagsz_gpu_pointer;
     int      *noffsz_gpu_pointer,*rankoffsz_gpu_pointer;
     int      *offdia_gpu, *offdia_gpu_pointer;


     if(verbose==1) printf(" GPU get eigenvalues, size of Lanczos vectors = %d \n", size);

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

       cudaSetDevice(0); cudaSetDeviceFlags( cudaDeviceMapHost ); cudaEventSynchronize(0);

       if(use_texture_lr==1){
         if(verbose==1) printf(" GPU go for texture \n");
         cudaMalloc((void**)&vec_in_gpu_pointer,sizeof(double)*size);
         bind_x(vec_in_gpu_pointer,size);
       }else{
         if(verbose==1) printf(" GPU HostAlloc ");
         cudaHostAlloc((void**)&vec_in_gpu    , size*sizeof(double),        cudaHostAllocMapped | cudaHostAllocPortable );
         cudaHostGetDevicePointer((void**) &vec_in_gpu_pointer ,vec_in_gpu , 0 );
       }

       if(verbose==1) printf(" GPU Host Alloc arrays \n");
       cudaHostAlloc((void**)&vec_tmp_gpu   , size*sizeof(double),        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&vec_out_gpu   , size*sizeof(double),        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&QUART_gpu     , size*sizeof(double),        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&diagsz_gpu    , size*sizeof(double),        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&offdiagsz_gpu , *offdiasize*sizeof(double), cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&noffsz_gpu    , size*sizeof(int)   ,        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&rankoffsz_gpu , *roff*sizeof(int)  ,        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&offdia_gpu    , size*sizeof(int)   ,        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaEventSynchronize(0);

       if(verbose==1) printf(" GPU allocate Device Pointer \n");
       cudaHostGetDevicePointer((void**)   &vec_tmp_gpu_pointer ,  vec_tmp_gpu    , 0 );
       cudaHostGetDevicePointer((void**)   &vec_out_gpu_pointer ,  vec_out_gpu    , 0 );
       cudaHostGetDevicePointer((void**)   &QUART_gpu_pointer ,    QUART_gpu      , 0 );
       cudaHostGetDevicePointer((void**)   &diagsz_gpu_pointer ,   diagsz_gpu     , 0 );
       cudaHostGetDevicePointer((void**)   &offdiagsz_gpu_pointer, offdiagsz_gpu  , 0 );
       cudaHostGetDevicePointer((void**)   &noffsz_gpu_pointer,    noffsz_gpu     , 0 );
       cudaHostGetDevicePointer((void**)   &rankoffsz_gpu_pointer, rankoffsz_gpu  , 0 );
       cudaHostGetDevicePointer((void**)   &offdia_gpu_pointer,    offdia_gpu     , 0 );
       cudaEventSynchronize(0);

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//
       if(verbose==1) printf(" GPU MemCpy arrays \n");
       cudaMemcpy(QUART_gpu,     QUART,       size*sizeof(double),        cudaMemcpyHostToHost);
       if(verbose==1) printf(" GPU MemCpy arrays diagsz \n");
       cudaMemcpy(diagsz_gpu,    diagsz,      size*sizeof(double),        cudaMemcpyHostToHost);
       if(verbose==1) printf(" GPU MemCpy arrays offdiagsz \n");
       cudaMemcpy(offdiagsz_gpu, offdiagsz,   *offdiasize*sizeof(double), cudaMemcpyHostToHost);
       if(verbose==1) printf(" GPU MemCpy arrays noffsz \n");
       cudaMemcpy(noffsz_gpu,    noffsz ,     size*sizeof(int),           cudaMemcpyHostToHost);
       if(verbose==1) printf(" GPU MemCpy arrays rankoffsz \n");
       cudaMemcpy(rankoffsz_gpu, rankoffsz,   *roff*sizeof(int),          cudaMemcpyHostToHost);
       if(verbose==1) printf(" Build up offdiag GPU \n");
       cudaEventSynchronize(0);
       offdia_gpu[0]=0; for(int istate=1; istate<size; istate++) { offdia_gpu[istate]=offdia_gpu[istate-1]+noffsz[istate-1]; };
       cudaEventSynchronize(0);
 //---------------------------------------------------------------------------------------//

  if(verbose==1) printf( " initialized, now run actual Lanczos \n " );

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

 if(use_texture_lr==0){
  for(int i=0;i<size;i++) vec_in_gpu[i]=1.0; cudaThreadSynchronize();
 }else{
   for(int i=0;i<size;i++) vec_out_gpu[i]=1.0; cudaThreadSynchronize();
   cudaMemcpy(vec_in_gpu_pointer,vec_out_gpu,size*sizeof(double),cudaMemcpyHostToDevice);
 }

  for(int iter=0;iter<Niter_lanczos;iter++){
    if(verbose==1) printf( " Lanczos iteration %d / %d \n", iter,Niter_lanczos);
    one_step_lanczos_cuda(blocksize,Niter_lanczos,size,iter,diag,subdiag,vec_tmp_gpu,vec_in_gpu,vec_out_gpu,vec_tmp_gpu_pointer,vec_in_gpu_pointer,vec_out_gpu_pointer,QUART_gpu_pointer,diagsz_gpu_pointer,noffsz_gpu_pointer,rankoffsz_gpu_pointer,offdiagsz_gpu_pointer,offdia_gpu_pointer,QUART_gpu,diagsz_gpu,noffsz_gpu,rankoffsz_gpu, offdiagsz_gpu, offdia_gpu);
  };

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

     cudaEventSynchronize(0);
     cudaFreeHost(vec_tmp_gpu);
     if(use_texture_lr==0){
      cudaFreeHost(vec_in_gpu);}
     else{
      unbind_x();
      cudaFree(vec_in_gpu_pointer);
     }
     cudaFreeHost(vec_out_gpu);
     cudaFreeHost(offdia_gpu);
     cudaFreeHost(QUART_gpu);
     cudaFreeHost(diagsz_gpu);
     cudaFreeHost(rankoffsz_gpu);
     cudaFreeHost(offdiagsz_gpu);
     cudaFreeHost(noffsz_gpu);
     cudaEventSynchronize(0);

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

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

extern "C" void lanczos_real_get_gs_cuda_(int *pblocksize, int *Niter_lanczos_,int *offdiasize, int *roff, 
int *ntot, double *QUART, double *diagsz, int *noffsz, int *rankoffsz, double *offdiagsz, double *vecp, double *GS)
{

 //---------------------------------------------------------------------------------------//
     int verbose=0;
 //---------------------------------------------------------------------------------------//

     if(verbose==1) printf(" start Lanczos Real on GPU \n" );
     int Niter_lanczos= *Niter_lanczos_;
     double diag[Niter_lanczos], subdiag[Niter_lanczos];
     int size = *ntot; int blocksize = *pblocksize;
     double *vec_in_gpu,*vec_out_gpu,*QUART_gpu,*diagsz_gpu,*offdiagsz_gpu,*vec_tmp_gpu;
     int *noffsz_gpu,*rankoffsz_gpu;
     double *vec_in_gpu_pointer,*vec_out_gpu_pointer,*vec_tmp_gpu_pointer;
     double *QUART_gpu_pointer,*diagsz_gpu_pointer,*offdiagsz_gpu_pointer;
     int *noffsz_gpu_pointer,*rankoffsz_gpu_pointer;
     int *offdia_gpu, *offdia_gpu_pointer;

     double *GS_gpu, *GS_gpu_pointer;

     if(verbose==1) printf(" GPU get Ground State, size of Lanczos vectors = %d ", size);

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

       cudaSetDevice(0); cudaSetDeviceFlags( cudaDeviceMapHost ); cudaEventSynchronize(0);

       cudaHostAlloc((void**)&GS_gpu        , size*sizeof(double),        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&vec_tmp_gpu   , size*sizeof(double),        cudaHostAllocMapped | cudaHostAllocPortable );

       if(use_texture_lr==1){
        cudaMalloc((void**)&vec_in_gpu_pointer,sizeof(double)*size);
        bind_x(vec_in_gpu_pointer,size);
      }else{
        cudaHostAlloc((void**)&vec_in_gpu    , size*sizeof(double),        cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostGetDevicePointer((void**)   &vec_in_gpu_pointer  ,  vec_in_gpu     , 0 );
      }

       cudaHostAlloc((void**)&vec_out_gpu   , size*sizeof(double),        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&QUART_gpu     , size*sizeof(double),        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&diagsz_gpu    , size*sizeof(double),        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&offdiagsz_gpu , *offdiasize*sizeof(double), cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&noffsz_gpu    , size*sizeof(int)   ,        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&rankoffsz_gpu , *roff*sizeof(int)  ,        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&offdia_gpu    , size*sizeof(int)   ,        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaEventSynchronize(0);

       cudaHostGetDevicePointer((void**)   &GS_gpu_pointer ,       GS_gpu         , 0 );
       cudaHostGetDevicePointer((void**)   &vec_tmp_gpu_pointer ,  vec_tmp_gpu    , 0 );
       cudaHostGetDevicePointer((void**)   &vec_out_gpu_pointer ,  vec_out_gpu    , 0 );
       cudaHostGetDevicePointer((void**)   &QUART_gpu_pointer ,    QUART_gpu      , 0 );
       cudaHostGetDevicePointer((void**)   &diagsz_gpu_pointer ,   diagsz_gpu     , 0 );
       cudaHostGetDevicePointer((void**)   &offdiagsz_gpu_pointer, offdiagsz_gpu  , 0 );
       cudaHostGetDevicePointer((void**)   &noffsz_gpu_pointer,    noffsz_gpu     , 0 );
       cudaHostGetDevicePointer((void**)   &rankoffsz_gpu_pointer, rankoffsz_gpu  , 0 );
       cudaHostGetDevicePointer((void**)   &offdia_gpu_pointer,    offdia_gpu     , 0 );
       cudaEventSynchronize(0);

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//
       cudaMemcpy(QUART_gpu,     QUART,       size*sizeof(double),        cudaMemcpyHostToHost);
       cudaMemcpy(diagsz_gpu,    diagsz,      size*sizeof(double),        cudaMemcpyHostToHost);
       cudaMemcpy(offdiagsz_gpu, offdiagsz,   *offdiasize*sizeof(double), cudaMemcpyHostToHost);
       cudaMemcpy(noffsz_gpu,    noffsz ,     size*sizeof(int),           cudaMemcpyHostToHost);
       cudaMemcpy(rankoffsz_gpu, rankoffsz,   *roff*sizeof(int),          cudaMemcpyHostToHost);
       cudaEventSynchronize(0);
       offdia_gpu[0]=0; for(int istate=1; istate<size; istate++) { offdia_gpu[istate]=offdia_gpu[istate-1]+noffsz[istate-1]; };
       cudaEventSynchronize(0);
 //---------------------------------------------------------------------------------------//

  if(verbose==1) printf( " initialized, now run actual Lanczos \n " );

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

  memset ((void **)GS_gpu, 0, size*sizeof(double));

 if(use_texture_lr==0){
  for(int i=0;i<size;i++) vec_in_gpu[i]=1.0; cudaThreadSynchronize();
 }else{
   for(int i=0;i<size;i++) vec_out_gpu[i]=1.0; cudaThreadSynchronize();
   cudaMemcpy(vec_in_gpu_pointer,vec_out_gpu,size*sizeof(double),cudaMemcpyHostToDevice);
 }

  double coef= 1.0/sqrt(cublasDdot(size,vec_in_gpu_pointer,1,vec_in_gpu_pointer,1))*vecp[0]; cudaThreadSynchronize();  cudaEventSynchronize(0);
  cublasDaxpy(size,coef,vec_in_gpu_pointer,1,GS_gpu_pointer,1); cudaThreadSynchronize(); cudaEventSynchronize(0);cudaThreadSynchronize();

   one_step_lanczos_cuda(blocksize, Niter_lanczos,size,0,diag,subdiag,vec_tmp_gpu,vec_in_gpu,vec_out_gpu,vec_tmp_gpu_pointer,vec_in_gpu_pointer,vec_out_gpu_pointer,QUART_gpu_pointer,diagsz_gpu_pointer,noffsz_gpu_pointer,rankoffsz_gpu_pointer,offdiagsz_gpu_pointer,offdia_gpu_pointer,QUART_gpu,diagsz_gpu,noffsz_gpu,rankoffsz_gpu, offdiagsz_gpu, offdia_gpu);

  for(int iter=1;iter<Niter_lanczos-1;iter++){

   coef=1.0/sqrt(cublasDdot(size,vec_out_gpu_pointer,1,vec_out_gpu_pointer,1)); cublasDscal (size,coef,vec_out_gpu_pointer,1);  cudaEventSynchronize(0);
   cublasDaxpy(size,vecp[iter],vec_out_gpu_pointer,1,GS_gpu_pointer,1); cudaThreadSynchronize(); cudaEventSynchronize(0);

   one_step_lanczos_cuda(blocksize,Niter_lanczos,size,iter,diag,subdiag,vec_tmp_gpu,vec_in_gpu,vec_out_gpu,vec_tmp_gpu_pointer,vec_in_gpu_pointer,vec_out_gpu_pointer,QUART_gpu_pointer,diagsz_gpu_pointer,noffsz_gpu_pointer,rankoffsz_gpu_pointer,offdiagsz_gpu_pointer,offdia_gpu_pointer,QUART_gpu,diagsz_gpu,noffsz_gpu,rankoffsz_gpu, offdiagsz_gpu, offdia_gpu);  cudaEventSynchronize(0);cudaThreadSynchronize();
  };

   coef=1.0/sqrt(cublasDdot(size,vec_out_gpu_pointer,1,vec_out_gpu_pointer,1)); cublasDscal(size,coef,vec_out_gpu_pointer,1);cudaEventSynchronize(0);
   cublasDaxpy(size,vecp[Niter_lanczos-1],vec_out_gpu_pointer,1,GS_gpu_pointer,1); cudaThreadSynchronize(); cudaEventSynchronize(0); 

   coef=1.0/sqrt(cublasDdot(size,GS_gpu_pointer,1,GS_gpu_pointer,1)); cudaThreadSynchronize();  cudaEventSynchronize(0);
   cublasDscal(size,coef,GS_gpu_pointer,1); cudaEventSynchronize(0); cudaThreadSynchronize();

   cudaMemcpy(GS,GS_gpu,size*sizeof(double),cudaMemcpyHostToHost);

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

     cudaEventSynchronize(0);
     cudaFreeHost(vec_tmp_gpu);
     if(use_texture_lr==0){
      cudaFreeHost(vec_in_gpu);}
     else{
      unbind_x();
      cudaFree(vec_in_gpu_pointer);
     }
     cudaFreeHost(vec_out_gpu);
     cudaFreeHost(offdia_gpu);
     cudaFreeHost(QUART_gpu);
     cudaFreeHost(diagsz_gpu);
     cudaFreeHost(rankoffsz_gpu);
     cudaFreeHost(offdiagsz_gpu);
     cudaFreeHost(noffsz_gpu);
     cudaEventSynchronize(0);

     cudaFreeHost(GS_gpu); cudaEventSynchronize(0);

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

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
