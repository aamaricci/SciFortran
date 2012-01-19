
 #include <stdio.h>
 #include <stdlib.h>
 #include <cuda_runtime.h>
 #include <cublas.h>
 #include <cuComplex.h>

//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************

const int use_texture_lc = 1 ;

#define WARP_SIZE 32
#define MAX_BLOCK 65500

texture<int2,1,cudaReadModeElementType> tex;

//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************

inline void cublasZscal_no_device (int n, cuDoubleComplex alpha, cuDoubleComplex *y, int incy)
 { for (int i=0; i<n; i++) { y[i] = cuCmul( alpha, y[i] ); } }

//********************************************
//********************************************
//********************************************
//********************************************
//********************************************

inline void cublasZaxpy_no_device (int n, cuDoubleComplex alpha, cuDoubleComplex *x, int incx,cuDoubleComplex *y, int incy)
 { for (int i=0; i<n; i++) { y[i] = cuCadd( y[i],  cuCmul( alpha, x[i] ) ); } }

//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************

inline cuDoubleComplex cublasZdotu_no_device (int n, cuDoubleComplex *x, int incx, cuDoubleComplex *y, int incy)
 { 
  cuDoubleComplex dot_; dot_=make_cuDoubleComplex(0.0,0.0) ;
  for (int i=0; i<n; i++) { dot_ = cuCadd(dot_,cuCmul( cuConj(x[i]), y[i] ) ); } 
  return dot_;
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

inline void   bind_x(cuDoubleComplex *x, int N) {    cudaBindTexture(0,tex,x,N*sizeof(cuDoubleComplex)); };
inline void unbind_x()                          {  cudaUnbindTexture(  tex  ); };

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

  //************************************************//
  //               Kernel Hmult                     //
  //************************************************//

  __global__  void Hmult_ker_complex_(int ngrid,int blocksize, int psize, cuDoubleComplex *vec_out, 
                const cuDoubleComplex *vec_in,  const double *QUART, const double *diagsz, const int *noffsz, const int *rankoffsz, 
		const cuDoubleComplex *offdiagsz, const int *offdia, const int use_texture_lc)
{
   int nhoffdiag,jstate,noff; 
   cuDoubleComplex hoffdiag,tmp; 
   int istate;

   for(int iii=0;iii<=ngrid;iii++){

    istate = (blockIdx.x+iii*MAX_BLOCK)*blocksize + threadIdx.x ;

    if(istate < psize ){
    tmp                =  make_cuDoubleComplex(QUART[istate]+diagsz[istate],0.0);
    vec_out[istate]    =  cuCmul(tmp,vec_in[istate]); 
            noff       =  noffsz[istate] ;
            nhoffdiag  =  offdia[istate] ;
 
   for(int irank=0;irank<noff;irank++)
    { jstate           =  rankoffsz[nhoffdiag+irank]  ;
      hoffdiag         =  cuConj(offdiagsz[nhoffdiag+irank]);
      vec_out[istate]  =  cuCadd(  vec_out[istate],  cuCmul(hoffdiag,vec_in[jstate-1]) ) ;
    };
  };
  };
}

  //************************************************//
  //               Kernel Hmult                     //
  //************************************************//

     //---------------------------------------//
     //---------------------------------------//
     //---------------------------------------//
     //---------------------------------------//

  __global__  void Hmult_ker_complex(int ngrid,int BLOCK_SIZE, int num_rows, cuDoubleComplex *y,
         const cuDoubleComplex *x, const double *QUART, const double *diagsz, const int *noffsz, const int *Aj,
         const cuDoubleComplex *Ax, const int *offdia, const int use_texture_lc)
{
   __shared__ cuDoubleComplex   sdata[16][WARP_SIZE];
   __shared__ int                ptrs[32][2];
   __shared__ double            temp2[32][2];

   const int warp_lane   = threadIdx.y ; const int thread_lane = threadIdx.x ; 

   int row_start; int row_end; int jj, row;

   for(int iii=0;iii<=ngrid;iii++){

    row = BLOCK_SIZE * (blockIdx.y + iii*MAX_BLOCK) + threadIdx.y;

    if(row<num_rows)
   {
        if(thread_lane==0)   ptrs[warp_lane][0] = offdia[row];
        if(thread_lane==1)   ptrs[warp_lane][1] = noffsz[row];
        if(thread_lane==2)  temp2[warp_lane][0] = QUART[row];
        if(thread_lane==3)  temp2[warp_lane][1] = diagsz[row];

        row_start = ptrs[warp_lane][0] ; row_end = ptrs[warp_lane][1]+row_start ;

        if(use_texture_lc==0)  {
          y[row] = cuCmul( make_cuDoubleComplex( temp2[warp_lane][0]+temp2[warp_lane][1],0.0 )  , x[row]);
        }else{
          y[row] = cuCmul( make_cuDoubleComplex( temp2[warp_lane][0]+temp2[warp_lane][1],0.0 )  , fetch_x(row+1) );
        };

        sdata[threadIdx.y][threadIdx.x]=make_cuDoubleComplex(0.0,0.0);

        if(use_texture_lc==1)
        {
        for(jj=row_start+thread_lane;jj<row_end;jj+=WARP_SIZE) 
           sdata[threadIdx.y][threadIdx.x]=cuCadd(sdata[threadIdx.y][threadIdx.x], cuCmul(cuConj(Ax[jj]),fetch_x(Aj[jj])));
        }else{
        for(jj=row_start+thread_lane;jj<row_end;jj+=WARP_SIZE) 
           sdata[threadIdx.y][threadIdx.x]=cuCadd(sdata[threadIdx.y][threadIdx.x], cuCmul(cuConj(Ax[jj]),x[Aj[jj]-1]));
        };

        if (thread_lane < 16) { sdata[threadIdx.y][threadIdx.x] = cuCadd(sdata[threadIdx.y][threadIdx.x], sdata[threadIdx.y][threadIdx.x + 16] ); };
        if (thread_lane <  8) { sdata[threadIdx.y][threadIdx.x] = cuCadd(sdata[threadIdx.y][threadIdx.x], sdata[threadIdx.y][threadIdx.x +  8] ); };
        if (thread_lane <  4) { sdata[threadIdx.y][threadIdx.x] = cuCadd(sdata[threadIdx.y][threadIdx.x], sdata[threadIdx.y][threadIdx.x +  4] ); };
        if (thread_lane <  2) { sdata[threadIdx.y][threadIdx.x] = cuCadd(sdata[threadIdx.y][threadIdx.x], sdata[threadIdx.y][threadIdx.x +  2] ); };
        if (thread_lane <  1) { sdata[threadIdx.y][threadIdx.x] = cuCadd(sdata[threadIdx.y][threadIdx.x], sdata[threadIdx.y][threadIdx.x +  1] ); };
        if (thread_lane == 0)   y[row] = cuCadd(y[row],sdata[threadIdx.y][threadIdx.x]) ;
   }; }; }

//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
__global__ void norm_vec_ker(int size , cuDoubleComplex *x, double *normv)
{
    int blocksize = blockDim.x ;
   __shared__ double  temp[512]; __shared__ double temp2[512]; int row; int  ii;

    temp[threadIdx.x]=0.;
    for(row=threadIdx.x; row<size; row+=blocksize)
   {
    temp2[threadIdx.x] = cuCabs(x[row]); temp2[threadIdx.x]*= temp2[threadIdx.x];
    temp[threadIdx.x] += temp2[threadIdx.x];
   };
   __syncthreads ();
    *normv=0.0; for(ii=0;ii<blocksize;ii++){ *normv+=temp[ii];}; *normv=sqrt(*normv); 
}
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
__global__ void real_scal_vec_ker( int size , cuDoubleComplex *x, cuDoubleComplex *y, double *normv)
{
   int blocksize = blockDim.x ; __shared__ double temp[512]; int row , ii; 

   temp[threadIdx.x]=0.0;
   for(row=threadIdx.x; row<size; row+=blocksize) { temp[threadIdx.x] += cuCreal( cuCmul(  cuConj( x[row] ) ,y[row]  ) ); };

  __syncthreads ();

   *normv=0.0; for(ii=0;ii<blocksize ;ii++){ *normv+=temp[ii];};
}
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************

__global__ void normalize_vec_ker(int ngrid,  int size ,cuDoubleComplex *x, double *normv)
{
    int blocksize = blockDim.x ;
    int row ;
    for(int iii=0;iii<=ngrid;iii++){
    row = blocksize * (blockIdx.x+iii*MAX_BLOCK) + threadIdx.x ;
    if(row<size) { x[row]=cuCmul( make_cuDoubleComplex(1.0 / *normv ,0.0), x[row] ); };};
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

extern "C" void hmult_sz_complex_cuda_rout_(int *pblocksize, int *offdiasize, int *roff, int* ntot, 
double* QUART, double* diagsz, cuDoubleComplex* vec_in, cuDoubleComplex* vec_out, int* noffsz, int* rankoffsz, 
cuDoubleComplex* offdiagsz )
{

 //---------------------------------------------------------------------------------------//
     int verbose=0;
 //---------------------------------------------------------------------------------------//

     if(verbose==1) printf(" start Hmult GPU \n" );

     int blocksize=*pblocksize; int size = *ntot; 

    int nb=(size-size % blocksize)/blocksize+1; int ngrid=nb/MAX_BLOCK; if(ngrid>0) nb=MAX_BLOCK;
    dim3 bl(1,nb),th(WARP_SIZE,blocksize);

    if(verbose==1) printf( " --------------- \n  Nblock=%d Ngrid=%d \n ----------------- \n ",nb,ngrid);

     cuDoubleComplex  *vec_in_gpu,*vec_out_gpu;
     double           *QUART_gpu,*diagsz_gpu;
     cuDoubleComplex  *offdiagsz_gpu;
     int              *noffsz_gpu,*rankoffsz_gpu;
     cuDoubleComplex  *vec_in_gpu_pointer,*vec_out_gpu_pointer;
     double           *QUART_gpu_pointer,*diagsz_gpu_pointer;
     cuDoubleComplex  *offdiagsz_gpu_pointer;
     int              *noffsz_gpu_pointer,*rankoffsz_gpu_pointer;
     int              *offdia_gpu, *offdia_gpu_pointer;

     if(verbose==1) printf(" GPU , size of Lanczos vector = %d ", size);

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

       cudaSetDevice(0); cudaSetDeviceFlags( cudaDeviceMapHost ); cudaEventSynchronize(0);

       if(use_texture_lc==1){
         cudaMalloc((void**)&vec_in_gpu_pointer,sizeof(cuDoubleComplex)*size);
         bind_x(vec_in_gpu_pointer,size);
       }else{
         cudaHostAlloc((void**)&vec_in_gpu    , size*sizeof(cuDoubleComplex),  cudaHostAllocMapped | cudaHostAllocPortable );
         cudaHostGetDevicePointer((void**)  &vec_in_gpu_pointer  ,  vec_in_gpu    , 0 );
       }

       cudaHostAlloc((void**)&vec_out_gpu   , size*sizeof(cuDoubleComplex),        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&QUART_gpu     , size*sizeof(double),                 cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&diagsz_gpu    , size*sizeof(double),                 cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&offdiagsz_gpu , *offdiasize*sizeof(cuDoubleComplex), cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&noffsz_gpu    , size*sizeof(int)   ,                 cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&rankoffsz_gpu , *roff*sizeof(int)  ,                 cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&offdia_gpu    , size*sizeof(int)   ,                 cudaHostAllocMapped | cudaHostAllocPortable );
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

      if(use_texture_lc==0){
       cudaMemcpy(vec_in_gpu, vec_in, size*sizeof(cuDoubleComplex), cudaMemcpyHostToHost);
      }else{
        cudaMemcpy(vec_in_gpu_pointer, vec_in, size*sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
      }
       cudaMemcpy(vec_out_gpu,   vec_out,     size*sizeof(cuDoubleComplex),        cudaMemcpyHostToHost);
       cudaMemcpy(QUART_gpu,     QUART,       size*sizeof(double),        cudaMemcpyHostToHost);
       cudaMemcpy(diagsz_gpu,    diagsz,      size*sizeof(double),        cudaMemcpyHostToHost);
       cudaMemcpy(offdiagsz_gpu, offdiagsz,   *offdiasize*sizeof(cuDoubleComplex), cudaMemcpyHostToHost);
       cudaMemcpy(noffsz_gpu,    noffsz ,     size*sizeof(int),           cudaMemcpyHostToHost);
       cudaMemcpy(rankoffsz_gpu, rankoffsz,   *roff*sizeof(int),          cudaMemcpyHostToHost);
       cudaEventSynchronize(0);
       offdia_gpu[0]=0; for(int istate=1; istate<size; istate++) { offdia_gpu[istate]=offdia_gpu[istate-1]+noffsz[istate-1]; };
       cudaEventSynchronize(0);

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

  if(verbose==1) printf(" call kernel  \n ");
  Hmult_ker_complex<<<bl,th>>>(ngrid,blocksize,size,vec_out_gpu_pointer,vec_in_gpu_pointer,QUART_gpu_pointer,diagsz_gpu_pointer,
                              noffsz_gpu_pointer,rankoffsz_gpu_pointer,offdiagsz_gpu_pointer,offdia_gpu_pointer,use_texture_lc); 
                              cudaEventSynchronize(0); cudaThreadSynchronize();
  if(verbose==1) printf(" .....done.....  \n ");

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

     cudaMemcpy(vec_out,vec_out_gpu,size*sizeof(cuDoubleComplex),cudaMemcpyHostToHost);
     cudaEventSynchronize(0);
     cudaFreeHost(vec_out_gpu);
    if(use_texture_lc==0){
     cudaFreeHost(vec_in_gpu);
    }else{
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

 void one_step_lanczos_cuda_complex(int blocksize, int Niter, int size, int iter, double *diag, double *subdiag, 
 cuDoubleComplex *vec_tmp_gpu, cuDoubleComplex *vec_in_gpu, cuDoubleComplex *vec_out_gpu, 
 cuDoubleComplex *vec_tmp_gpu_pointer, cuDoubleComplex *vec_in_gpu_pointer, cuDoubleComplex *vec_out_gpu_pointer, 
 double *QUART_gpu_pointer, double *diagsz_gpu_pointer, int *noffsz_gpu_pointer, int *rankoffsz_gpu_pointer, 
 cuDoubleComplex *offdiagsz_gpu_pointer, int *offdia_gpu_pointer, double *QUART_gpu, double *diagsz_gpu, int *noffsz_gpu, 
 int *rankoffsz_gpu, cuDoubleComplex *offdiagsz_gpu, int *offdia_gpu)
{

   int psize=size; int verbose=0;

   double normv; double *normv_ker; cudaMalloc((void**)&normv_ker,sizeof(double)); double *normv_loc; normv_loc = &normv;
   cuDoubleComplex coef;

   int nb  =(size -size % blocksize)/blocksize + 1; int ngrid=nb/MAX_BLOCK; if(ngrid>0) nb=MAX_BLOCK;
   dim3 bl(1,nb),th(WARP_SIZE,blocksize);
   if(verbose==1) printf( " --------------- \n  Nblock=%d Ngrid=%d \n ----------------- \n ",nb,ngrid);

   int nb2 =(size-size % 256)/256 + 1; 
   int ngrid2=nb2/MAX_BLOCK; if(ngrid2>0) nb2=MAX_BLOCK;

   norm_vec_ker<<<1,512>>>(size,vec_in_gpu_pointer,normv_ker); cudaEventSynchronize(0); cudaThreadSynchronize();
   normalize_vec_ker<<<nb2,256>>>(ngrid2,size,vec_in_gpu_pointer,normv_ker);

   Hmult_ker_complex<<<bl,th>>>(ngrid,blocksize,psize,vec_out_gpu_pointer,vec_in_gpu_pointer,QUART_gpu_pointer,diagsz_gpu_pointer,noffsz_gpu_pointer,rankoffsz_gpu_pointer,offdiagsz_gpu_pointer,offdia_gpu_pointer,use_texture_lc); cudaEventSynchronize(0); cudaThreadSynchronize();

  if(iter>0){ coef=make_cuDoubleComplex(-subdiag[iter],0.0); cublasZaxpy_no_device(size,coef,vec_tmp_gpu,1,vec_out_gpu,1); }; cudaEventSynchronize(0); cudaThreadSynchronize();

  if(use_texture_lc==0){
   cudaMemcpy(vec_tmp_gpu,vec_in_gpu,size*sizeof(cuDoubleComplex),cudaMemcpyHostToHost);cudaEventSynchronize(0);cudaThreadSynchronize();
  }else{
   cudaMemcpy(vec_tmp_gpu,vec_in_gpu_pointer,size*sizeof(cuDoubleComplex),cudaMemcpyDeviceToHost);cudaEventSynchronize(0);cudaThreadSynchronize();
  }

   real_scal_vec_ker<<<1,512>>> (size,vec_out_gpu_pointer,vec_in_gpu_pointer,normv_ker);
   cudaMemcpy(normv_loc,normv_ker,sizeof(double),cudaMemcpyDeviceToHost);
   diag[iter]=*normv_loc ; 

   cudaEventSynchronize(0);cudaThreadSynchronize();
   coef=make_cuDoubleComplex(-diag[iter],0.0);
   cublasZaxpy_no_device(size,coef,vec_tmp_gpu,1,vec_out_gpu,1); cudaEventSynchronize(0); cudaThreadSynchronize();

   normv = sqrt(cuCabs(cublasZdotu_no_device(size,vec_out_gpu,1,vec_out_gpu,1))); cudaEventSynchronize(0); cudaThreadSynchronize();

   if(iter<Niter-1) subdiag[iter+1]=normv;

  if(use_texture_lc==0){
   cudaMemcpy( vec_in_gpu, vec_out_gpu, size*sizeof(cuDoubleComplex), cudaMemcpyHostToHost); cudaEventSynchronize(0); cudaThreadSynchronize();
  }else{
   cudaMemcpy( vec_in_gpu_pointer, vec_out_gpu, size*sizeof(cuDoubleComplex), cudaMemcpyHostToDevice); cudaEventSynchronize(0); cudaThreadSynchronize();
  }

  cudaFree(normv_ker);
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

extern "C" void lanczos_dynamic_cuda_complex_(int *pblocksize, int *Niter_lanczos_,int *offdiasize, int *roff, 
int *ntot, double *QUART, double *diagsz, int *noffsz, int *rankoffsz, cuDoubleComplex *offdiagsz, double *diag, 
double *subdiag , cuDoubleComplex *vecinit)
{

 //---------------------------------------------------------------------------------------//
     int verbose=0;
 //---------------------------------------------------------------------------------------//

     int Niter_lanczos=*Niter_lanczos_; 

     if(verbose==1) printf(" start Lanczos Real on GPU \n" );

     int size = *ntot; int blocksize=*pblocksize; 

     cuDoubleComplex   *vec_in_gpu,*vec_out_gpu;
     double            *QUART_gpu,*diagsz_gpu;
     cuDoubleComplex   *offdiagsz_gpu,*vec_tmp_gpu;

     int               *noffsz_gpu,*rankoffsz_gpu;
     cuDoubleComplex   *vec_in_gpu_pointer,*vec_out_gpu_pointer,*vec_tmp_gpu_pointer;
     double            *QUART_gpu_pointer,*diagsz_gpu_pointer;
     cuDoubleComplex   *offdiagsz_gpu_pointer;
     int               *noffsz_gpu_pointer,*rankoffsz_gpu_pointer;
     int               *offdia_gpu, *offdia_gpu_pointer;

     if(verbose==1) printf(" GPU get eigenvalues, size of Lanczos vectors = %d ", size);

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

       cudaSetDevice(0); cudaSetDeviceFlags( cudaDeviceMapHost ); cudaEventSynchronize(0);

       if(use_texture_lc==1){
         cudaMalloc((void**)&vec_in_gpu_pointer,sizeof(cuDoubleComplex)*size);
         bind_x(vec_in_gpu_pointer,size);
       }else{
         cudaHostAlloc((void**)&vec_in_gpu    , size*sizeof(cuDoubleComplex),        cudaHostAllocMapped | cudaHostAllocPortable );
         cudaHostGetDevicePointer((void**)  &vec_in_gpu_pointer  ,  vec_in_gpu    , 0 );
       }

       cudaHostAlloc((void**)&vec_tmp_gpu   , size*sizeof(cuDoubleComplex),        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&vec_out_gpu   , size*sizeof(cuDoubleComplex),        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&QUART_gpu     , size*sizeof(double),                 cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&diagsz_gpu    , size*sizeof(double),                 cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&offdiagsz_gpu , *offdiasize*sizeof(cuDoubleComplex), cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&noffsz_gpu    , size*sizeof(int)   ,                 cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&rankoffsz_gpu , *roff*sizeof(int)  ,                 cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&offdia_gpu    , size*sizeof(int)   ,                 cudaHostAllocMapped | cudaHostAllocPortable );
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

       cudaMemcpy(QUART_gpu,     QUART,       size*sizeof(double),        cudaMemcpyHostToHost);
       cudaMemcpy(diagsz_gpu,    diagsz,      size*sizeof(double),        cudaMemcpyHostToHost);
       cudaMemcpy(offdiagsz_gpu, offdiagsz,   *offdiasize*sizeof(cuDoubleComplex), cudaMemcpyHostToHost);
       cudaMemcpy(noffsz_gpu,    noffsz ,     size*sizeof(int),           cudaMemcpyHostToHost);
       cudaMemcpy(rankoffsz_gpu, rankoffsz,   *roff*sizeof(int),          cudaMemcpyHostToHost);
       cudaEventSynchronize(0);
       offdia_gpu[0]=0; for(int istate=1; istate<size; istate++) { offdia_gpu[istate]=offdia_gpu[istate-1]+noffsz[istate-1]; };
       cudaEventSynchronize(0);

  if(verbose==1) printf( " initialized, now run actual Lanczos \n " );

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

  cudaThreadSynchronize(); cudaEventSynchronize(0); 

 if(use_texture_lc==0){
   cudaMemcpy(vec_in_gpu,vecinit,size*sizeof(cuDoubleComplex),cudaMemcpyHostToHost);
 }else{
   cudaMemcpy(vec_out_gpu,vecinit,size*sizeof(cuDoubleComplex),cudaMemcpyHostToHost);
   cudaMemcpy(vec_in_gpu_pointer,vec_out_gpu,size*sizeof(cuDoubleComplex),cudaMemcpyHostToDevice);
 }

  for(int iter=0;iter<Niter_lanczos;iter++){
   if(verbose==1) printf( " Lanczos iteration %d / %d \n", iter,Niter_lanczos);
   one_step_lanczos_cuda_complex(blocksize,Niter_lanczos,size,iter,diag,subdiag,vec_tmp_gpu,vec_in_gpu,vec_out_gpu,vec_tmp_gpu_pointer,vec_in_gpu_pointer,vec_out_gpu_pointer,QUART_gpu_pointer,diagsz_gpu_pointer,noffsz_gpu_pointer,rankoffsz_gpu_pointer,offdiagsz_gpu_pointer,offdia_gpu_pointer,QUART_gpu,diagsz_gpu,noffsz_gpu,rankoffsz_gpu, offdiagsz_gpu, offdia_gpu);

  };

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

     cudaEventSynchronize(0);
     cudaFreeHost(vec_tmp_gpu);
    if(use_texture_lc==0){
     cudaFreeHost(vec_in_gpu);
    }else{
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

extern "C" void lanczos_cuda_complex_(int *pblocksize,  int *Niter_lanczos_,int *offdiasize, int *roff, int *ntot, 
double *QUART, double *diagsz, int *noffsz, int *rankoffsz, cuDoubleComplex *offdiagsz, double *diag, double *subdiag )
{

 //---------------------------------------------------------------------------------------//
     int verbose=0;
 //---------------------------------------------------------------------------------------//

     int Niter_lanczos=*Niter_lanczos_; 
     if(verbose==1) printf(" start Lanczos Real on GPU \n" );
     int size = *ntot; int blocksize = *pblocksize;

     cuDoubleComplex  *vec_in_gpu,*vec_out_gpu;
     double           *QUART_gpu,*diagsz_gpu;
     cuDoubleComplex  *offdiagsz_gpu,*vec_tmp_gpu;
     int *noffsz_gpu, *rankoffsz_gpu;
     cuDoubleComplex  *vec_in_gpu_pointer,*vec_out_gpu_pointer,*vec_tmp_gpu_pointer;
     double           *QUART_gpu_pointer,*diagsz_gpu_pointer;
     cuDoubleComplex  *offdiagsz_gpu_pointer;
     int              *noffsz_gpu_pointer,*rankoffsz_gpu_pointer;
     int              *offdia_gpu, *offdia_gpu_pointer;

     if(verbose==1) printf(" GPU get eigenvalues, size of Lanczos vectors = %d ", size);

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

       cudaSetDevice(0); cudaSetDeviceFlags( cudaDeviceMapHost ); cudaEventSynchronize(0);

       if(use_texture_lc==1){
         cudaMalloc((void**)&vec_in_gpu_pointer,sizeof(cuDoubleComplex)*size);
         bind_x(vec_in_gpu_pointer,size);
       }else{
         cudaHostAlloc((void**)&vec_in_gpu    , size*sizeof(cuDoubleComplex),        cudaHostAllocMapped | cudaHostAllocPortable );
         cudaHostGetDevicePointer((void**) &vec_in_gpu_pointer ,vec_in_gpu , 0 );
       }

       cudaHostAlloc((void**)&vec_tmp_gpu   , size*sizeof(cuDoubleComplex),        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&vec_out_gpu   , size*sizeof(cuDoubleComplex),        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&QUART_gpu     , size*sizeof(double),                 cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&diagsz_gpu    , size*sizeof(double),                 cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&offdiagsz_gpu , *offdiasize*sizeof(cuDoubleComplex), cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&noffsz_gpu    , size*sizeof(int)   ,                 cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&rankoffsz_gpu , *roff*sizeof(int)  ,                 cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&offdia_gpu    , size*sizeof(int)   ,                 cudaHostAllocMapped | cudaHostAllocPortable );
       cudaEventSynchronize(0);

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

       cudaMemcpy(QUART_gpu,     QUART,       size*sizeof(double),        cudaMemcpyHostToHost);
       cudaMemcpy(diagsz_gpu,    diagsz,      size*sizeof(double),        cudaMemcpyHostToHost);
       cudaMemcpy(offdiagsz_gpu, offdiagsz,   *offdiasize*sizeof(cuDoubleComplex), cudaMemcpyHostToHost);
       cudaMemcpy(noffsz_gpu,    noffsz ,     size*sizeof(int),           cudaMemcpyHostToHost);
       cudaMemcpy(rankoffsz_gpu, rankoffsz,   *roff*sizeof(int),          cudaMemcpyHostToHost);
       cudaEventSynchronize(0);
       offdia_gpu[0]=0; for(int istate=1; istate<size; istate++) { offdia_gpu[istate]=offdia_gpu[istate-1]+noffsz[istate-1]; };
       cudaEventSynchronize(0);

  if(verbose==1) printf( " initialized, now run actual Lanczos \n " );

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

 if(use_texture_lc==0){
  for(int i=0;i<size;i++) vec_in_gpu[i]=make_cuDoubleComplex(1.0,0.0); cudaThreadSynchronize();
 }else{
   for(int i=0;i<size;i++) vec_out_gpu[i]=make_cuDoubleComplex(1.0,0.0); cudaThreadSynchronize();
   cudaMemcpy(vec_in_gpu_pointer,vec_out_gpu,size*sizeof(cuDoubleComplex),cudaMemcpyHostToDevice);
 }

  for(int iter=0;iter<Niter_lanczos;iter++){
    if(verbose==1) printf( " Lanczos iteration %d / %d \n", iter,Niter_lanczos);
    one_step_lanczos_cuda_complex(blocksize,Niter_lanczos,size,iter,diag,subdiag,vec_tmp_gpu,vec_in_gpu,vec_out_gpu,vec_tmp_gpu_pointer,vec_in_gpu_pointer,vec_out_gpu_pointer,QUART_gpu_pointer,diagsz_gpu_pointer,noffsz_gpu_pointer,rankoffsz_gpu_pointer,offdiagsz_gpu_pointer,offdia_gpu_pointer,QUART_gpu,diagsz_gpu,noffsz_gpu,rankoffsz_gpu, offdiagsz_gpu, offdia_gpu);
  };

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

     cudaEventSynchronize(0);
     cudaFreeHost(vec_tmp_gpu);
    if(use_texture_lc==0){
     cudaFreeHost(vec_in_gpu);
    }else{
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

extern "C" void lanczos_get_gs_cuda_complex_(int *pblocksize, int *Niter_lanczos_,int *offdiasize, int *roff, 
int *ntot, double *QUART, double *diagsz, int *noffsz, int *rankoffsz, cuDoubleComplex *offdiagsz, double *vecp, 
cuDoubleComplex *GS)
{

 //---------------------------------------------------------------------------------------//
     int verbose=0;
 //---------------------------------------------------------------------------------------//

     if(verbose==1) printf(" start Lanczos Real on GPU \n" );

     int Niter_lanczos=*Niter_lanczos_;
     double diag[Niter_lanczos], subdiag[Niter_lanczos];
     int size = *ntot; int blocksize =*pblocksize;

     cuDoubleComplex  *vec_in_gpu,*vec_out_gpu;
     double           *QUART_gpu,*diagsz_gpu;
     cuDoubleComplex  *offdiagsz_gpu,*vec_tmp_gpu;
     int              *noffsz_gpu,*rankoffsz_gpu;
     cuDoubleComplex  *vec_in_gpu_pointer,*vec_out_gpu_pointer,*vec_tmp_gpu_pointer;
     double           *QUART_gpu_pointer,*diagsz_gpu_pointer;
     cuDoubleComplex  *offdiagsz_gpu_pointer;
     int              *noffsz_gpu_pointer,*rankoffsz_gpu_pointer;
     int              *offdia_gpu, *offdia_gpu_pointer;
     cuDoubleComplex  *GS_gpu, *GS_gpu_pointer;

     if(verbose==1) printf(" \n GPU get Ground State, size of Lanczos vectors = %d \n ", size);

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

       cudaSetDevice(0); cudaSetDeviceFlags( cudaDeviceMapHost ); cudaEventSynchronize(0);

       if(use_texture_lc==1){
        cudaMalloc((void**)&vec_in_gpu_pointer,sizeof(cuDoubleComplex)*size);
        bind_x(vec_in_gpu_pointer,size);
      }else{
        cudaHostAlloc((void**)&vec_in_gpu    , size*sizeof(cuDoubleComplex),        cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostGetDevicePointer((void**)   &vec_in_gpu_pointer  ,  vec_in_gpu     , 0 );
      }

       cudaHostAlloc((void**)&GS_gpu        , size*sizeof(cuDoubleComplex),        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&vec_tmp_gpu   , size*sizeof(cuDoubleComplex),        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&vec_out_gpu   , size*sizeof(cuDoubleComplex),        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&QUART_gpu     , size*sizeof(double),                 cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&diagsz_gpu    , size*sizeof(double),                 cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&offdiagsz_gpu , *offdiasize*sizeof(cuDoubleComplex), cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&noffsz_gpu    , size*sizeof(int)   ,                 cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&rankoffsz_gpu , *roff*sizeof(int)  ,                 cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&offdia_gpu    , size*sizeof(int)   ,                 cudaHostAllocMapped | cudaHostAllocPortable );
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

       cudaMemcpy(QUART_gpu,     QUART,       size*sizeof(double),        cudaMemcpyHostToHost);
       cudaMemcpy(diagsz_gpu,    diagsz,      size*sizeof(double),        cudaMemcpyHostToHost);
       cudaMemcpy(offdiagsz_gpu, offdiagsz,   *offdiasize*sizeof(cuDoubleComplex), cudaMemcpyHostToHost);
       cudaMemcpy(noffsz_gpu,    noffsz ,     size*sizeof(int),           cudaMemcpyHostToHost);
       cudaMemcpy(rankoffsz_gpu, rankoffsz,   *roff*sizeof(int),          cudaMemcpyHostToHost);
       cudaEventSynchronize(0);
       offdia_gpu[0]=0; for(int istate=1; istate<size; istate++) { offdia_gpu[istate]=offdia_gpu[istate-1]+noffsz[istate-1]; };
       cudaEventSynchronize(0);

  if(verbose==1) printf( " initialized, now run actual Lanczos \n " );

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

  memset ((void **)GS_gpu, 0, size*sizeof(cuDoubleComplex)); cudaThreadSynchronize();  cudaEventSynchronize(0);

 if(use_texture_lc==0){
  for(int i=0;i<size;i++) vec_in_gpu[i]=make_cuDoubleComplex(1.0,0.0); cudaThreadSynchronize();
 }else{
   for(int i=0;i<size;i++) vec_out_gpu[i]=make_cuDoubleComplex(1.0,0.0); cudaThreadSynchronize();
   cudaMemcpy(vec_in_gpu_pointer,vec_out_gpu,size*sizeof(cuDoubleComplex),cudaMemcpyHostToDevice);
 }

 if(verbose==1) printf(" start lanczos iterations \n ");

 double *normv_ker; double normv; double *normv_loc; cudaMalloc((void**)&normv_ker,sizeof(double)); cuDoubleComplex coef; 

 if(use_texture_lc==0){
   coef = make_cuDoubleComplex(1.0/sqrt(cuCabs(cublasZdotu_no_device(size,vec_in_gpu,1,vec_in_gpu,1)))*vecp[0],0); cudaThreadSynchronize();  cudaEventSynchronize(0);
   cublasZaxpy_no_device(size,coef,vec_in_gpu,1,GS_gpu,1); cudaThreadSynchronize(); cudaEventSynchronize(0);cudaThreadSynchronize();
 }else{
   norm_vec_ker<<<1,512>>>(size,vec_in_gpu_pointer,normv_ker); cudaEventSynchronize(0); cudaThreadSynchronize();
   normv_loc=&normv; cudaMemcpy(normv_loc,normv_ker,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(vec_out_gpu,vec_in_gpu,size*sizeof(cuDoubleComplex),cudaMemcpyDeviceToHost); 
   coef=make_cuDoubleComplex(vecp[0]/normv,0.0);
   cublasZaxpy_no_device(size,coef,vec_out_gpu,1,GS_gpu,1); cudaThreadSynchronize(); cudaEventSynchronize(0);cudaThreadSynchronize();
 };

  if(verbose==1) printf( " first step ... \n ");

  one_step_lanczos_cuda_complex(blocksize, Niter_lanczos,size,0,diag,subdiag,vec_tmp_gpu,vec_in_gpu,vec_out_gpu,vec_tmp_gpu_pointer,vec_in_gpu_pointer,vec_out_gpu_pointer,QUART_gpu_pointer,diagsz_gpu_pointer,noffsz_gpu_pointer,rankoffsz_gpu_pointer,offdiagsz_gpu_pointer,offdia_gpu_pointer,QUART_gpu,diagsz_gpu,noffsz_gpu,rankoffsz_gpu, offdiagsz_gpu, offdia_gpu);

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//
 
  for(int iter=1;iter<Niter_lanczos-1;iter++){

   if(verbose==1) printf( " iterations = %d \n " , iter);

   coef = make_cuDoubleComplex(1.0/sqrt(cuCabs(cublasZdotu_no_device(size,vec_out_gpu,1,vec_out_gpu,1))),0.0); 
   cublasZscal_no_device (size,coef, vec_out_gpu,1);  cudaEventSynchronize(0);

   coef = make_cuDoubleComplex(vecp[iter],0);
   cublasZaxpy_no_device(size,coef,vec_out_gpu,1,GS_gpu,1); cudaThreadSynchronize(); cudaEventSynchronize(0);

   one_step_lanczos_cuda_complex(blocksize,Niter_lanczos,size,iter,diag,subdiag,vec_tmp_gpu,vec_in_gpu,vec_out_gpu,vec_tmp_gpu_pointer,vec_in_gpu_pointer,vec_out_gpu_pointer,QUART_gpu_pointer,diagsz_gpu_pointer,noffsz_gpu_pointer,rankoffsz_gpu_pointer,offdiagsz_gpu_pointer,offdia_gpu_pointer,QUART_gpu,diagsz_gpu,noffsz_gpu,rankoffsz_gpu, offdiagsz_gpu, offdia_gpu);  cudaEventSynchronize(0);cudaThreadSynchronize();
  };

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

   if(verbose==1) printf("done...\n");

   coef=make_cuDoubleComplex(1.0/sqrt(cuCabs(cublasZdotu_no_device(size,vec_out_gpu,1,vec_out_gpu,1))),0.0); 
   cublasZscal_no_device(size,coef,vec_out_gpu,1);cudaEventSynchronize(0);

   coef=make_cuDoubleComplex(vecp[Niter_lanczos-1],0.0);
   cublasZaxpy_no_device(size,coef,vec_out_gpu,1,GS_gpu,1); cudaThreadSynchronize(); cudaEventSynchronize(0); 

   coef=make_cuDoubleComplex(1.0/sqrt(cuCabs(cublasZdotu_no_device(size,GS_gpu,1,GS_gpu,1))),0.0);  cudaThreadSynchronize();  cudaEventSynchronize(0);
   cublasZscal_no_device(size,coef,GS_gpu,1); cudaEventSynchronize(0); cudaThreadSynchronize();

   cudaMemcpy(GS,GS_gpu,size*sizeof(cuDoubleComplex),cudaMemcpyHostToHost); cudaEventSynchronize(0); cudaThreadSynchronize();

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

     cudaEventSynchronize(0);
     cudaFreeHost(vec_tmp_gpu);
    if(use_texture_lc==0){
     cudaFreeHost(vec_in_gpu);
    }else{
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

     cudaFree(normv_ker); cudaEventSynchronize(0); cudaFreeHost(GS_gpu);cudaFree(GS_gpu_pointer); cudaEventSynchronize(0);

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
