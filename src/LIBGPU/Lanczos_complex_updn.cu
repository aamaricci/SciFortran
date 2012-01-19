
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
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************

#define IBSET(a,b) ((a) |= (1<<(b)))
#define IBCLR(a,b) ((a) &= ~(1<<(b)))
#define BTEST(a,b) ((a) & (1<<(b)))

#define WARP_SIZE 32
#define MAX_BLOCK 65500
#define BLOCK_DIM 16
#define USE_TRANSPOSE 1
const int use_texture_updo_comp = 0 ;


texture<int2,1,cudaReadModeElementType> tex, texdn;

inline void   bind_x(cuDoubleComplex *x, int size) { cudaBindTexture(0,tex,x,size*sizeof(cuDoubleComplex)); };
inline void unbind_x()                    { cudaUnbindTexture(tex); };

inline void   bind_x_(cuDoubleComplex *y, int size) { cudaBindTexture(0,texdn,y,size*sizeof(cuDoubleComplex)); };
inline void unbind_x_()                    { cudaUnbindTexture(texdn);    };

__inline__  __device__ cuDoubleComplex fetch_x(const int& i)
  { int  jj = 2*(i); int2 v  = tex1Dfetch(tex,jj); double rr  = __hiloint2double(v.y, v.x); v  = tex1Dfetch(tex,jj+1);
    double im  = __hiloint2double(v.y, v.x); return make_cuDoubleComplex(rr,im); }
__inline__  __device__ cuDoubleComplex fetch_x_(const int& i)
  { int  jj = 2*(i); int2 v  = tex1Dfetch(texdn,jj); double rr  = __hiloint2double(v.y, v.x); v  = tex1Dfetch(texdn,jj+1);
    double im  = __hiloint2double(v.y, v.x); return make_cuDoubleComplex(rr,im); }

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

inline void cublasZscal_no_device (int n, cuDoubleComplex alpha, cuDoubleComplex *y, int incy)
 { for (int i=0; i<n; i++) { y[i] = cuCmul( alpha, y[i] ); } }

//********************************************
//********************************************

inline void cublasZaxpy_no_device (int n, cuDoubleComplex alpha, cuDoubleComplex *x, int incx,cuDoubleComplex *y, int incy)
 { for (int i=0; i<n; i++) { y[i] = cuCadd( y[i],  cuCmul( alpha, x[i] ) ); } }

//********************************************
//********************************************

inline cuDoubleComplex cublasZdotu_no_device (int n, cuDoubleComplex *x, int incx, cuDoubleComplex *y, int incy)
 {
  cuDoubleComplex dot_; dot_=make_cuDoubleComplex(0.0,0.0) ;
  for (int i=0; i<n; i++) { dot_ = cuCadd(dot_,cuCmul( cuConj(x[i]), y[i] ) ); } return dot_;
 }

//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************

__global__ void norm_vec_ker_complex_updo(int size , cuDoubleComplex *x, double *normv)

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
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************

__global__ void real_scal_vec_ker_complex_updo( int size , cuDoubleComplex *x, cuDoubleComplex *y, double *normv)

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
//********************************************

__global__ void normalize_vec_ker_complex_updo(int ngrid, int size ,cuDoubleComplex *x, double *normv)
{
    int blocksize = blockDim.x ;
    int row ;
    for(int iii=0;iii<=ngrid;iii++)
    {
     row= blocksize * (blockIdx.x + iii*MAX_BLOCK)+ threadIdx.x ;
     if(row<size)
      { x[row]=cuCmul( make_cuDoubleComplex(1.0 / *normv ,0.0), x[row] ); };
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

  //************************************************//
  //               Kernel Hmult                     //
  //************************************************//

     //---------------------------------------//
     //---------------------------------------//
     //---------------------------------------//
     //---------------------------------------//
     //---------------------------------------//
     //---------------------------------------//
     //---------------------------------------//
     //---------------------------------------//

__global__ void transpose_ker_complex(int nb1, int nb2, int ngrid, cuDoubleComplex *output, cuDoubleComplex *data, int width, int height)
{

       __shared__ cuDoubleComplex block[BLOCK_DIM][BLOCK_DIM+1];

       // nb1 horizontal  nb2 vertical ; 
       // bl1 index block horizontal bl2 index block vertical

      int iv,index,bl1,bl2;

      for(iv=0; iv<=ngrid; iv++) 
       {

        index=blockIdx.x+iv*MAX_BLOCK;
 
        bl1= index % nb1 ;  
        bl2= (index - bl1)/nb1 ;

       // read
        unsigned int xIndex1 = bl1*BLOCK_DIM + threadIdx.x;
        unsigned int yIndex1 = bl2*BLOCK_DIM + threadIdx.y;
        if((xIndex1 < width) && (yIndex1 < height))
        { unsigned int index_in = yIndex1 * width + xIndex1;
          block[threadIdx.y][threadIdx.x] = data[index_in];
        }
        __syncthreads();

       // write
        xIndex1 =bl2* BLOCK_DIM + threadIdx.x;
        yIndex1 =bl1* BLOCK_DIM + threadIdx.y;
        if((xIndex1 < height) && (yIndex1 < width))
        {      unsigned int index_out = yIndex1 * height + xIndex1;
               output[index_out] = block[threadIdx.x][threadIdx.y];
        };

        };
}

     //---------------------------------------//
     //---------------------------------------//
     //---------------------------------------//
     //---------------------------------------//
     //---------------------------------------//

  __global__  void Hmult_ker_updo_comp0(int nup, int ndn, int norb, int ngrid, int BLOCK_SIZE, int num_rows, cuDoubleComplex *y, const cuDoubleComplex *x,
                             const double *QUART, const double *diagup, const double *diagdn,
                             const int *noffup, const int *noffdn, const int *rankoffup, const int *rankoffdn,
                             const cuDoubleComplex *offdiagup, const cuDoubleComplex *offdiagdn, const int *offdiaup, const int *offdiadn,
                             const int *statesup, const int *statesdn, const int *UMASK,
                             const int use_texture_updo_comp, const int *iorbup, const int *iorbdn)
{

   int istate,iup,idn; int iii; double diagdn_ ; const int thread_lane = threadIdx.x; const int interval = WARP_SIZE;

// DIAGONAL TERMS, COALESCENT READING

   ////////////////////////////////////////////////////////////
   for(iii=0; iii<=ngrid; iii++){
    idn = BLOCK_SIZE * (blockIdx.y+iii*MAX_BLOCK) + threadIdx.y ;
    if(idn<ndn)
    {
     diagdn_=diagdn[idn]; 
     for(iup=thread_lane; iup<nup; iup+=interval ) 
      { istate       =   iup + idn * nup;
       if(use_texture_updo_comp==0)
       {
        y[istate]   = cuCadd( y[istate], cuCmul( make_cuDoubleComplex( diagup[iup]+diagdn_,0.0 ) , x[istate]));
       }else{
        y[istate]   = cuCadd( y[istate], cuCmul( make_cuDoubleComplex( diagup[iup]+diagdn_,0.0 ) , fetch_x(istate)));
       };
      };
    };
   };
   ////////////////////////////////////////////////////////////

}

     //---------------------------------------//
     //---------------------------------------//
     //---------------------------------------//
     //---------------------------------------//
     //---------------------------------------//
     //---------------------------------------//
     //---------------------------------------//

  __global__  void Hmult_ker_updo_comp1(int nup, int ndn, int norb, int ngrid, int BLOCK_SIZE, int num_rows, 
                             cuDoubleComplex *y, const cuDoubleComplex *x, 
                             const double *QUART, const double *diagup, const double *diagdn, const int *noffup, const int *noffdn,  
                             const int *rankoffup, const int *rankoffdn, const cuDoubleComplex *offdiagup, cuDoubleComplex *offdiagdn, 
                             const int *offdiaup, const int *offdiadn, const int *statesup, const int *statesdn, const int *UMASK,
                             const int use_texture_updo_comp, const int *iorbup, const int *iorbdn)
{
   __shared__ double QUART_[32*32];
   __shared__ short int iorbup_[32],iorbdn_[32],UMASK_[32*32];

   int istate,iup,idn ;  int ii,iii,iorb,statesup_; double quart_;

   for(iorb=threadIdx.x;iorb<norb*norb;iorb+=BLOCK_SIZE)  UMASK_[iorb]= UMASK[iorb];
   for(iorb=threadIdx.x;iorb<norb*norb;iorb+=BLOCK_SIZE)  QUART_[iorb]= QUART[iorb];
   for(iorb=threadIdx.x;iorb<norb;iorb+=BLOCK_SIZE)      iorbup_[iorb]=iorbup[iorb];
   for(iorb=threadIdx.x;iorb<norb;iorb+=BLOCK_SIZE)      iorbdn_[iorb]=iorbdn[iorb];

   __syncthreads();

// DIAGONAL TERMS, QUARTIC TERMS, OPTIMIZED FOR SHARED MEMORY

   ////////////////////////////////////////////////////////////
   for(iii=0; iii<=ngrid; iii++){
    iup = BLOCK_SIZE * (blockIdx.x+iii*MAX_BLOCK) + threadIdx.x ;
    if(iup<nup)
    {
      statesup_=statesup[iup];
      for(ii=0;ii<norb;ii++) 
      {
       if(UMASK_[ii*norb+ii]!=0 && BTEST(statesup_,iorbup_[ii]-1)>0)
        {
         quart_ = QUART_[ii*norb+ii];
         for(idn=0; idn<ndn; idn+=1)
          {
           if(BTEST(statesdn[idn],iorbdn_[ii]-1)>0){
            istate       = iup + idn * nup;
            if(use_texture_updo_comp==1){
             y[istate] = cuCadd(y[istate], cuCmul( make_cuDoubleComplex(quart_,0.0) , fetch_x(istate))); 
            }else{
             y[istate] = cuCadd(y[istate], cuCmul( make_cuDoubleComplex(quart_,0.0) , x[istate]));
            };
           };
          };
        };
      };
   };}; 
   ////////////////////////////////////////////////////////////

 }

     //---------------------------------------//
     //---------------------------------------//
     //---------------------------------------//
     //---------------------------------------//

  __global__  void Hmult_ker_updo_comp2(int nup, int ndn, int norb, int ngrid, int BLOCK_SIZE, int num_rows,  
                             cuDoubleComplex *y, const cuDoubleComplex *x, 
                             const double *QUART, const double *diagup, const double *diagdn,
                             const int *noffup, const int *noffdn, const int *rankoffup, const int *rankoffdn,
                             const cuDoubleComplex *offdiagup, const cuDoubleComplex *offdiagdn,
                             const int *offdiaup, const int *offdiadn, const int *statesup, const int *statesdn,
                             const int *UMASK, const int use_texture_updo_comp, const int *iorbup, const int *iorbdn)
{
   int jstate,iup ; int istatemin,istatemax,noff,nhoffdiag,ii,iii,irank,jup; cuDoubleComplex hoffdiag ;

// HOPPING UP

   for(iii=0; iii<=ngrid; iii++){
     iup = BLOCK_SIZE * (blockIdx.x+iii*MAX_BLOCK) + threadIdx.x ;
     if(iup<nup)
    {
   ////////////////////////////////////////////////////////////
      istatemin  = iup; istatemax  = iup + (ndn-1) * nup;
      noff       = noffup[iup];
      nhoffdiag  = offdiaup[iup];
      for(irank=0;irank<noff;irank++)
       {
       hoffdiag    =  cuConj(offdiagup[nhoffdiag+irank]);
       jup         =  rankoffup[nhoffdiag+irank]-1;
       jstate      = jup  ;
       for(ii=istatemin;ii<=istatemax;ii+=nup)
        {
        if(use_texture_updo_comp==1){
         y[ii] =cuCadd(y[ii],cuCmul( hoffdiag , fetch_x(jstate)));
        }else{
         y[ii] =cuCadd(y[ii],cuCmul( hoffdiag , x[jstate]));
        };
         jstate += nup;
        };
      };
   ////////////////////////////////////////////////////////////

   }; };

 }

     //---------------------------------------//
     //---------------------------------------//
     //---------------------------------------//
     //---------------------------------------//

  __global__  void Hmult_ker_updo_comp3(int nup, int ndn, int norb, int ngrid, int BLOCK_SIZE, int num_rows, 
                             cuDoubleComplex *y, const cuDoubleComplex *x, 
                             const double *QUART, const double *diagup, const double *diagdn,
                             const int *noffup, const int *noffdn, const int *rankoffup, const int *rankoffdn,
                             const cuDoubleComplex *offdiagup, const cuDoubleComplex *offdiagdn, const int *offdiaup, const int *offdiadn,
                             const int *statesup, const int *statesdn, const int *UMASK, const int use_texture_updo_comp, const int *iorbup, const int *iorbdn)
{

   int jstate,idn ; int istatemin,istatemax,noff,nhoffdiag,ii,iii,irank,jdn;
   cuDoubleComplex hoffdiag ; const int thread_lane = threadIdx.x;

// HOPPING DN

   ////////////////////////////////////////////////////////////
   for(iii=0; iii<=ngrid; iii++){
     idn = BLOCK_SIZE * (blockIdx.y+iii*MAX_BLOCK) + threadIdx.y ;
     if(idn<ndn)
    {
      istatemin   = idn*nup;
      istatemax   = idn*nup+nup-1;
      noff        = noffdn[idn];
      nhoffdiag   = offdiadn[idn];
      for(irank=0;irank<noff;irank++)
       {
       hoffdiag    =  cuConj(offdiagdn[nhoffdiag+irank]);
       jdn         =  rankoffdn[nhoffdiag+irank]-1;
       jstate = jdn*nup+thread_lane ;
       for(ii=istatemin+thread_lane;ii<=istatemax;ii+=WARP_SIZE)
        { 
         if(use_texture_updo_comp==0){
            y[ii]  = cuCadd(y[ii],cuCmul(hoffdiag , x[jstate])); 
         }else{
          if(use_texture_updo_comp==1){
            y[ii]  = cuCadd(y[ii],cuCmul(hoffdiag , fetch_x(jstate)));
           }else{
            y[ii]  = cuCadd(y[ii],cuCmul(hoffdiag , fetch_x_(jstate)));
           };
         };
          jstate += WARP_SIZE; };
       };
    };
   };
   ////////////////////////////////////////////////////////////

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
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************

 void one_step_lanczos_cuda(int nup, int ndn, int norb, int blocksize, int Niter, int size, int iter, double *diag, 
 double *subdiag, cuDoubleComplex *vec_tmp_gpu, cuDoubleComplex *vec_in_gpu, cuDoubleComplex *vec_out_gpu, cuDoubleComplex *vec_tmp_gpu_pointer, 
 cuDoubleComplex *vec_in_gpu_pointer, cuDoubleComplex *vec_out_gpu_pointer, double *QUART_gpu_pointer, double *diagup_gpu_pointer, double *diagdn_gpu_pointer,
 int *noffup_gpu_pointer, int *noffdn_gpu_pointer, int *rankoffup_gpu_pointer, int *rankoffdn_gpu_pointer, 
 cuDoubleComplex *offdiagup_gpu_pointer, cuDoubleComplex *offdiagdn_gpu_pointer, int *offdiaup_gpu_pointer, int *offdiadn_gpu_pointer,
 int *statesup_gpu_pointer,int *statesdn_gpu_pointer, int *UMASK_gpu_pointer,int *iorbup_gpu_pointer,int *iorbdn_gpu_pointer,
 cuDoubleComplex *tmp_transp, cuDoubleComplex *tmp_transp_pointer)

{
   int verbose=0 ; int psize=size; int nb, ngrid,nb2,nb1;
   int nb3 =(size-size % 256 )/256 + 1; int ngrid2=nb3/MAX_BLOCK; if(ngrid2>0) nb3=MAX_BLOCK;
   double normv; double *normv_ker; cudaMalloc((void**)&normv_ker,sizeof(double)); double *normv_loc; normv_loc = &normv; cuDoubleComplex coef;

   norm_vec_ker_complex_updo<<<1,512>>>(size,vec_in_gpu_pointer,normv_ker); cudaEventSynchronize(0); cudaThreadSynchronize();
   normalize_vec_ker_complex_updo<<<nb3,256>>>(ngrid2, size, vec_in_gpu_pointer, normv_ker);

   //------------------------------------------------------------------------------------------------------------------//

   for(int i=0;i<size;i++) vec_out_gpu[i]=make_cuDoubleComplex(0.0,0.0); cudaThreadSynchronize();  cudaEventSynchronize(0);

   //------------------------------------------------------------------------------------------------------------------//
   //------------------------------------------------------------------------------------------------------------------//
  if(USE_TRANSPOSE==1)
   {

   nb1=(nup-nup % BLOCK_DIM)/BLOCK_DIM+1; nb2=(ndn-ndn % BLOCK_DIM)/BLOCK_DIM+1; nb=nb1*nb2; ngrid=nb/MAX_BLOCK; if(ngrid>0)nb=MAX_BLOCK;
   dim3 threads(BLOCK_DIM, BLOCK_DIM, 1);
   transpose_ker_complex<<<nb,threads>>>(nb1,nb2,ngrid,tmp_transp_pointer,vec_in_gpu_pointer,nup,ndn);
   cudaEventSynchronize(0); cudaThreadSynchronize();

   int bb=blocksize/WARP_SIZE; if(bb<1) bb=1; nb=(nup-nup % bb)/bb+1; ngrid=nb/MAX_BLOCK; if(ngrid>0)nb=MAX_BLOCK; dim3 t1(1,nb),t2(WARP_SIZE,bb);
   Hmult_ker_updo_comp3<<<t1,t2>>>(ndn,nup,norb,ngrid,bb,psize,vec_out_gpu_pointer,tmp_transp_pointer,QUART_gpu_pointer,
                         diagdn_gpu_pointer,diagup_gpu_pointer,noffdn_gpu_pointer,noffup_gpu_pointer,rankoffdn_gpu_pointer,rankoffup_gpu_pointer,
                         offdiagdn_gpu_pointer,offdiagup_gpu_pointer,offdiadn_gpu_pointer,offdiaup_gpu_pointer,
                         statesdn_gpu_pointer,statesup_gpu_pointer,UMASK_gpu_pointer,-use_texture_updo_comp,iorbdn_gpu_pointer,iorbup_gpu_pointer);
   if(verbose==1) printf( " done.... \n " );
   cudaEventSynchronize(0); cudaThreadSynchronize();

   nb2=(nup-nup % BLOCK_DIM)/BLOCK_DIM+1; nb1=(ndn-ndn % BLOCK_DIM)/BLOCK_DIM+1; nb=nb1*nb2; ngrid=nb/MAX_BLOCK; if(ngrid>0)nb=MAX_BLOCK;
   transpose_ker_complex<<<nb,threads>>>(nb1,nb2,ngrid,vec_in_gpu_pointer,tmp_transp_pointer,ndn,nup);
   cudaEventSynchronize(0); cudaThreadSynchronize();

   transpose_ker_complex<<<nb,threads>>>(nb1,nb2,ngrid,tmp_transp_pointer,vec_out_gpu_pointer,ndn,nup);
   cudaEventSynchronize(0);cudaThreadSynchronize();
 
   if(use_texture_updo_comp==0){
    cudaMemcpy(vec_out_gpu,tmp_transp,size*sizeof(cuDoubleComplex),cudaMemcpyHostToHost);cudaEventSynchronize(0);cudaThreadSynchronize();
   }else{
    cudaMemcpy(vec_out_gpu,tmp_transp_pointer,size*sizeof(cuDoubleComplex),cudaMemcpyDeviceToHost);cudaEventSynchronize(0);cudaThreadSynchronize();
   }

   cudaEventSynchronize(0);cudaThreadSynchronize();
   if(verbose==1) printf( "transposition applied....\n");

   //------------------------------------------------------------------------------------------------------------------//
  }else{
   nb=(nup-nup % blocksize)/blocksize+1; ngrid=nb/MAX_BLOCK; if(ngrid>0)nb=MAX_BLOCK;
   Hmult_ker_updo_comp2<<<nb,blocksize>>>(nup,ndn,norb,ngrid,blocksize,psize,vec_out_gpu_pointer,vec_in_gpu_pointer,QUART_gpu_pointer,
                         diagup_gpu_pointer,diagdn_gpu_pointer,noffup_gpu_pointer,noffdn_gpu_pointer,rankoffup_gpu_pointer,rankoffdn_gpu_pointer,
                         offdiagup_gpu_pointer,offdiagdn_gpu_pointer,offdiaup_gpu_pointer,offdiadn_gpu_pointer,
                         statesup_gpu_pointer,statesdn_gpu_pointer,UMASK_gpu_pointer,use_texture_updo_comp
                        ,iorbup_gpu_pointer,iorbdn_gpu_pointer);
   if(verbose==1) printf( " done.... \n " );
   cudaEventSynchronize(0); cudaThreadSynchronize();
  };
   //-------------------------------------------------------------------------------------------------------------------//
   //-------------------------------------------------------------------------------------------------------------------//

   int bb=blocksize/WARP_SIZE; if(bb<1) bb=1; nb=(ndn-ndn % bb)/bb+1; ngrid=nb/MAX_BLOCK; if(ngrid>0)nb=MAX_BLOCK; dim3 t1(1,nb), t2(WARP_SIZE,bb);
   Hmult_ker_updo_comp3<<<t1,t2>>>(nup,ndn,norb,ngrid,bb,psize,vec_out_gpu_pointer,vec_in_gpu_pointer,QUART_gpu_pointer,
                         diagup_gpu_pointer,diagdn_gpu_pointer,noffup_gpu_pointer,noffdn_gpu_pointer,rankoffup_gpu_pointer,rankoffdn_gpu_pointer,
                         offdiagup_gpu_pointer,offdiagdn_gpu_pointer,offdiaup_gpu_pointer,offdiadn_gpu_pointer,
                         statesup_gpu_pointer,statesdn_gpu_pointer,UMASK_gpu_pointer,use_texture_updo_comp,iorbup_gpu_pointer,iorbdn_gpu_pointer);
   if(verbose==1) printf( " done.... \n " );
   cudaEventSynchronize(0); cudaThreadSynchronize();
   //------------------------------------------------------------------------------------------------------------------//
   bb=blocksize/WARP_SIZE; if(bb<1) bb=1; nb=(ndn- ndn % bb)/bb+1; ngrid=nb/MAX_BLOCK; if(ngrid>0)nb=MAX_BLOCK; dim3 tt1(1,nb), tt2(WARP_SIZE,bb); 
   Hmult_ker_updo_comp0<<<tt1,tt2>>>(nup,ndn,norb,ngrid,bb,psize,
                         vec_out_gpu_pointer,vec_in_gpu_pointer,QUART_gpu_pointer, diagup_gpu_pointer,diagdn_gpu_pointer,
                         noffup_gpu_pointer,noffdn_gpu_pointer,rankoffup_gpu_pointer,rankoffdn_gpu_pointer,
                         offdiagup_gpu_pointer,offdiagdn_gpu_pointer,offdiaup_gpu_pointer,offdiadn_gpu_pointer,
                         statesup_gpu_pointer,statesdn_gpu_pointer,UMASK_gpu_pointer,use_texture_updo_comp
                        ,iorbup_gpu_pointer,iorbdn_gpu_pointer);
   if(verbose==1) printf( " done.... \n " );
   cudaEventSynchronize(0); cudaThreadSynchronize();
   //------------------------------------------------------------------------------------------------------------------//
   nb=(nup-nup % blocksize)/blocksize+1; ngrid=nb/MAX_BLOCK; if(ngrid>0)nb=MAX_BLOCK; 
   if(verbose==1) printf( " call kernel ... \n ");
   Hmult_ker_updo_comp1<<<nb,blocksize>>>(nup,ndn,norb,ngrid,blocksize,psize,vec_out_gpu_pointer,vec_in_gpu_pointer,QUART_gpu_pointer, 
                         diagup_gpu_pointer,diagdn_gpu_pointer, noffup_gpu_pointer,noffdn_gpu_pointer,rankoffup_gpu_pointer,rankoffdn_gpu_pointer,
                         offdiagup_gpu_pointer,offdiagdn_gpu_pointer,offdiaup_gpu_pointer,offdiadn_gpu_pointer,
                         statesup_gpu_pointer,statesdn_gpu_pointer,UMASK_gpu_pointer,use_texture_updo_comp
                        ,iorbup_gpu_pointer,iorbdn_gpu_pointer); 
   if(verbose==1) printf( " done.... \n " );
   cudaEventSynchronize(0); cudaThreadSynchronize();
   //------------------------------------------------------------------------------------------------------------------//

  if(iter>0){ coef=make_cuDoubleComplex(-subdiag[iter],0.0); cublasZaxpy_no_device(size,coef,vec_tmp_gpu,1,vec_out_gpu,1); }; cudaEventSynchronize(0); cudaThreadSynchronize();

  if(use_texture_updo_comp==0){
   cudaMemcpy(vec_tmp_gpu,vec_in_gpu,size*sizeof(cuDoubleComplex),cudaMemcpyHostToHost);cudaEventSynchronize(0);cudaThreadSynchronize();
  }else{
   cudaMemcpy(vec_tmp_gpu,vec_in_gpu_pointer,size*sizeof(cuDoubleComplex),cudaMemcpyDeviceToHost);cudaEventSynchronize(0);cudaThreadSynchronize();
  }

   real_scal_vec_ker_complex_updo<<<1,512>>> (size,vec_out_gpu_pointer,vec_in_gpu_pointer,normv_ker);
   cudaMemcpy(normv_loc,normv_ker,sizeof(double),cudaMemcpyDeviceToHost);
   diag[iter]=*normv_loc ;

   cudaEventSynchronize(0);cudaThreadSynchronize();
   coef=make_cuDoubleComplex(-diag[iter],0.0);
   cublasZaxpy_no_device(size,coef,vec_tmp_gpu,1,vec_out_gpu,1); cudaEventSynchronize(0); cudaThreadSynchronize();
   normv = sqrt(cuCabs(cublasZdotu_no_device(size,vec_out_gpu,1,vec_out_gpu,1))); cudaEventSynchronize(0); cudaThreadSynchronize();
   if(iter<Niter-1) subdiag[iter+1]=normv;

   if(use_texture_updo_comp==0){
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
//********************************************
//********************************************
//********************************************

extern "C" void lanczos_complex_updo_cuda_(int *norbs, int *pblocksize, int *Niter_lanczos_,
   int *offdiasizeup, int* offdiasizedn, int *roffup, int *roffdn, int *ntot,  int *sizup, int *sizdn, 
   double *QUART, double *diagup, double *diagdn, int *noffup, int *noffdn, int *rankoffup, int *rankoffdn, 
   cuDoubleComplex *offdiagup, cuDoubleComplex *offdiagdn, double *diag, double *subdiag, int *UMASK, int *statesup, int *statesdn , int *iorbup, int *iorbdn)
{

 //---------------------------------------------------------------------------------------//
     int verbose=0;
 //---------------------------------------------------------------------------------------//

     int Niter_lanczos= *Niter_lanczos_; int nup=*sizup; int ndn=*sizdn; int norb=*norbs;

     if(verbose==1) printf(" start Lanczos Complex on GPU \n" );

     int size = *ntot; int blocksize= *pblocksize * WARP_SIZE ; 

     cuDoubleComplex   *vec_in_gpu,*vec_out_gpu;
     double            *QUART_gpu,*diagup_gpu,*diagdn_gpu;
     cuDoubleComplex   *offdiagup_gpu,*offdiagdn_gpu,*vec_tmp_gpu;
     int               *noffup_gpu,*noffdn_gpu,*rankoffup_gpu,*rankoffdn_gpu;
     cuDoubleComplex   *vec_in_gpu_pointer,*vec_out_gpu_pointer,*vec_tmp_gpu_pointer;
     double            *QUART_gpu_pointer,*diagup_gpu_pointer,*diagdn_gpu_pointer;
     cuDoubleComplex   *offdiagup_gpu_pointer,*offdiagdn_gpu_pointer;
     int               *noffup_gpu_pointer,*noffdn_gpu_pointer,*rankoffup_gpu_pointer,*rankoffdn_gpu_pointer;
     int               *offdiaup_gpu,*offdiadn_gpu, *offdiaup_gpu_pointer,*offdiadn_gpu_pointer;
     int               *statesup_gpu,*statesdn_gpu,*statesup_gpu_pointer,*statesdn_gpu_pointer;
     int               *UMASK_gpu,*UMASK_gpu_pointer;
     int               *iorbup_gpu,*iorbdn_gpu,*iorbup_gpu_pointer,*iorbdn_gpu_pointer;
     cuDoubleComplex   *tmp_transp,*tmp_transp_pointer;

     if(verbose==1) printf(" GPU get eigenvalues, size of Lanczos vectors = %d ", size);

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

       cudaSetDevice(0); cudaSetDeviceFlags( cudaDeviceMapHost ); cudaEventSynchronize(0);

       cudaHostAlloc((void**)&vec_tmp_gpu   , size*sizeof(cuDoubleComplex),  cudaHostAllocMapped | cudaHostAllocPortable );
 
       if(USE_TRANSPOSE==1){ 
        if(use_texture_updo_comp==0){
          cudaHostAlloc((void**)&tmp_transp, size*sizeof(cuDoubleComplex),   cudaHostAllocMapped | cudaHostAllocPortable );
          cudaHostGetDevicePointer((void**)  &tmp_transp_pointer,tmp_transp,0);
         }else{
          cudaMalloc((void**)&tmp_transp_pointer,size*sizeof(cuDoubleComplex) );
          bind_x_(tmp_transp_pointer,size);
         };
       };
 
       if(use_texture_updo_comp==1){
         cudaMalloc((void**)&vec_in_gpu_pointer,sizeof(cuDoubleComplex)*size);
         bind_x(vec_in_gpu_pointer,size);
       }else{
         cudaHostAlloc((void**)&vec_in_gpu , size*sizeof(cuDoubleComplex),   cudaHostAllocMapped | cudaHostAllocPortable );
         cudaHostGetDevicePointer((void**)  &vec_in_gpu_pointer  ,  vec_in_gpu    , 0 );
       }

       //-----------------------------------------------------------------------------------------------//
       //-----------------------------------------------------------------------------------------------//
        cudaHostAlloc((void**)&statesup_gpu    , nup *sizeof(int),           cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&statesdn_gpu    , ndn *sizeof(int),           cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&diagdn_gpu    , ndn*sizeof(double),           cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&offdiagdn_gpu , *offdiasizedn*sizeof(cuDoubleComplex), cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&noffdn_gpu    , ndn*sizeof(int)   ,           cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&rankoffdn_gpu , *roffdn*sizeof(int)  ,        cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&diagup_gpu    , nup*sizeof(double),           cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&offdiagup_gpu , *offdiasizeup*sizeof(cuDoubleComplex), cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&noffup_gpu    , nup*sizeof(int)   ,           cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&rankoffup_gpu , *roffup*sizeof(int)  ,        cudaHostAllocMapped | cudaHostAllocPortable );

        cudaMemcpy(statesup_gpu,  statesup,     nup * sizeof(int),           cudaMemcpyHostToHost);
        cudaMemcpy(statesdn_gpu,  statesdn,     ndn * sizeof(int),           cudaMemcpyHostToHost);
        cudaMemcpy(diagup_gpu,    diagup,       nup*sizeof(double),          cudaMemcpyHostToHost);
        cudaMemcpy(offdiagup_gpu, offdiagup,   *offdiasizeup*sizeof(cuDoubleComplex), cudaMemcpyHostToHost);
        cudaMemcpy(noffup_gpu,    noffup ,      nup*sizeof(int),             cudaMemcpyHostToHost);
        cudaMemcpy(rankoffup_gpu, rankoffup,   *roffup*sizeof(int),          cudaMemcpyHostToHost);
        cudaMemcpy(diagdn_gpu,    diagdn,       ndn*sizeof(double),          cudaMemcpyHostToHost);
        cudaMemcpy(offdiagdn_gpu, offdiagdn,   *offdiasizedn*sizeof(cuDoubleComplex), cudaMemcpyHostToHost);
        cudaMemcpy(noffdn_gpu,    noffdn ,      ndn*sizeof(int),             cudaMemcpyHostToHost);
        cudaMemcpy(rankoffdn_gpu, rankoffdn,   *roffdn*sizeof(int),          cudaMemcpyHostToHost);
        cudaEventSynchronize(0);

        cudaHostGetDevicePointer((void**)  &statesup_gpu_pointer , statesup_gpu  , 0 );
        cudaHostGetDevicePointer((void**)  &statesdn_gpu_pointer , statesdn_gpu  , 0 );
        cudaHostGetDevicePointer((void**)  &diagdn_gpu_pointer ,   diagdn_gpu    , 0 );
        cudaHostGetDevicePointer((void**)  &offdiagdn_gpu_pointer, offdiagdn_gpu , 0 );
        cudaHostGetDevicePointer((void**)  &noffdn_gpu_pointer,    noffdn_gpu    , 0 );
        cudaHostGetDevicePointer((void**)  &rankoffdn_gpu_pointer, rankoffdn_gpu , 0 );
        cudaHostGetDevicePointer((void**)  &diagup_gpu_pointer ,   diagup_gpu    , 0 );
        cudaHostGetDevicePointer((void**)  &offdiagup_gpu_pointer, offdiagup_gpu , 0 );
        cudaHostGetDevicePointer((void**)  &noffup_gpu_pointer,    noffup_gpu    , 0 );
        cudaHostGetDevicePointer((void**)  &rankoffup_gpu_pointer, rankoffup_gpu , 0 );
        cudaEventSynchronize(0);
       //-----------------------------------------------------------------------------------------------//
       //-----------------------------------------------------------------------------------------------//

       cudaHostAlloc((void**)&offdiaup_gpu    , nup*sizeof(int)   ,   cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&offdiadn_gpu    , ndn*sizeof(int)   ,   cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostGetDevicePointer((void**)  &offdiadn_gpu_pointer,    offdiadn_gpu, 0 );
       cudaHostGetDevicePointer((void**)  &offdiaup_gpu_pointer,    offdiaup_gpu, 0 );
       cudaEventSynchronize(0);
       offdiaup_gpu[0]=0; for(int istate=1; istate<nup; istate++) { offdiaup_gpu[istate]=offdiaup_gpu[istate-1]+noffup[istate-1]; };
       offdiadn_gpu[0]=0; for(int istate=1; istate<ndn; istate++) { offdiadn_gpu[istate]=offdiadn_gpu[istate-1]+noffdn[istate-1]; };
       cudaEventSynchronize(0);
       cudaEventSynchronize(0);

       cudaHostAlloc((void**)&iorbup_gpu     , norb*sizeof(int),            cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&iorbdn_gpu     , norb*sizeof(int),            cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&QUART_gpu     , norb*norb*sizeof(double),     cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&UMASK_gpu,      norb*norb*sizeof(int),        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaMemcpy(iorbup_gpu,     iorbup,      norb*sizeof(int),            cudaMemcpyHostToHost);
       cudaMemcpy(iorbdn_gpu,     iorbdn,      norb*sizeof(int),            cudaMemcpyHostToHost);
       cudaMemcpy(QUART_gpu,     QUART,      norb*norb*sizeof(double),      cudaMemcpyHostToHost);
       cudaMemcpy(UMASK_gpu,     UMASK,      norb*norb*sizeof(int),         cudaMemcpyHostToHost);
       cudaHostGetDevicePointer((void**)  &iorbup_gpu_pointer ,    iorbup_gpu   , 0 );
       cudaHostGetDevicePointer((void**)  &iorbdn_gpu_pointer ,    iorbdn_gpu   , 0 );
       cudaHostGetDevicePointer((void**)  &QUART_gpu_pointer ,    QUART_gpu     , 0 );
       cudaHostGetDevicePointer((void**)  &UMASK_gpu_pointer ,   UMASK_gpu      , 0 );

       cudaHostAlloc((void**)&vec_out_gpu   , size*sizeof(cuDoubleComplex),        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaEventSynchronize(0);
       cudaHostGetDevicePointer((void**)  &vec_tmp_gpu_pointer ,  vec_tmp_gpu   , 0 );
       cudaHostGetDevicePointer((void**)  &vec_out_gpu_pointer ,  vec_out_gpu   , 0 );
       cudaEventSynchronize(0);

 //---------------------------------------------------------------------------------------//

  if(verbose==1) printf( " initialized, now run actual Lanczos \n " );

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

  cudaThreadSynchronize(); cudaEventSynchronize(0); 
 if(use_texture_updo_comp==0){
  for(int i=0;i<size;i++) vec_in_gpu[i]=make_cuDoubleComplex(1.0,0.0); cudaThreadSynchronize();
 }else{
   for(int i=0;i<size;i++) vec_out_gpu[i]=make_cuDoubleComplex(1.0,0.0); cudaThreadSynchronize();
   cudaMemcpy(vec_in_gpu_pointer,vec_out_gpu,size*sizeof(cuDoubleComplex),cudaMemcpyHostToDevice);
 }
   cudaThreadSynchronize(); cudaEventSynchronize(0);

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//
  for(int iter=0;iter<Niter_lanczos;iter++){
   if(verbose==1) printf( " Lanczos iteration %d / %d \n", iter,Niter_lanczos);
   one_step_lanczos_cuda(nup,ndn,norb,blocksize,Niter_lanczos,size,iter,diag,subdiag,vec_tmp_gpu,vec_in_gpu,vec_out_gpu,
               vec_tmp_gpu_pointer,vec_in_gpu_pointer,vec_out_gpu_pointer,QUART_gpu_pointer,
               diagup_gpu_pointer,diagdn_gpu_pointer,noffup_gpu_pointer,noffdn_gpu_pointer,
               rankoffup_gpu_pointer,rankoffdn_gpu_pointer,offdiagup_gpu_pointer,offdiagdn_gpu_pointer,
               offdiaup_gpu_pointer,offdiadn_gpu_pointer,statesup_gpu_pointer,statesdn_gpu_pointer,UMASK_gpu_pointer
              ,iorbup_gpu_pointer,iorbdn_gpu_pointer,tmp_transp,tmp_transp_pointer);
  };
 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

     cudaEventSynchronize(0);
     cudaFreeHost(vec_tmp_gpu);

     if(use_texture_updo_comp==0){
       cudaFreeHost(vec_in_gpu);
       if(USE_TRANSPOSE==1){ cudaFreeHost(tmp_transp);}
     }else{
       cudaFree(vec_in_gpu_pointer); 
       if(USE_TRANSPOSE==1){
         cudaFree(tmp_transp_pointer); unbind_x_();
       };
     }

     cudaFreeHost(vec_out_gpu);  cudaFreeHost(iorbup_gpu);   cudaFreeHost(iorbdn_gpu); 
     cudaFreeHost(QUART_gpu);    cudaFreeHost(UMASK_gpu);
     cudaFreeHost(offdiaup_gpu); cudaFreeHost(offdiadn_gpu);
     cudaFreeHost(statesup_gpu);  cudaFreeHost(statesdn_gpu);  cudaFreeHost(diagup_gpu);  
     cudaFreeHost(rankoffup_gpu);
     cudaFreeHost(offdiagup_gpu); cudaFreeHost(noffup_gpu);    cudaFreeHost(diagdn_gpu);
     cudaFreeHost(rankoffdn_gpu); cudaFreeHost(offdiagdn_gpu); cudaFreeHost(noffdn_gpu);

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
//********************************************
//********************************************

extern "C" void lanczos_complex_updo_dynamic_cuda_(int *norbs, int *pblocksize, 
   int *Niter_lanczos_,int *offdiasizeup, int* offdiasizedn, int *roffup, int *roffdn, int *ntot,  int *sizup, int *sizdn, 
   double *QUART, double *diagup, double *diagdn, int *noffup, int *noffdn, int *rankoffup, int *rankoffdn, 
   cuDoubleComplex *offdiagup, cuDoubleComplex *offdiagdn, double *diag, double *subdiag, int *UMASK, int *statesup, int *statesdn , 
   int *iorbup, int *iorbdn, cuDoubleComplex *vecinit)
{

 //---------------------------------------------------------------------------------------//
     int verbose=0;
 //---------------------------------------------------------------------------------------//

     int Niter_lanczos= *Niter_lanczos_; int nup=*sizup; int ndn=*sizdn; int norb=*norbs;

     if(verbose==1) printf(" start Lanczos Complex on GPU \n" );

     int size = *ntot; int blocksize= *pblocksize * WARP_SIZE ; 

     cuDoubleComplex   *vec_in_gpu,*vec_out_gpu;
     double            *QUART_gpu,*diagup_gpu,*diagdn_gpu;
     cuDoubleComplex   *offdiagup_gpu,*offdiagdn_gpu,*vec_tmp_gpu;
     int               *noffup_gpu,*noffdn_gpu,*rankoffup_gpu,*rankoffdn_gpu;
     cuDoubleComplex   *vec_in_gpu_pointer,*vec_out_gpu_pointer,*vec_tmp_gpu_pointer;
     double            *QUART_gpu_pointer,*diagup_gpu_pointer,*diagdn_gpu_pointer;
     cuDoubleComplex   *offdiagup_gpu_pointer,*offdiagdn_gpu_pointer;
     int               *noffup_gpu_pointer,*noffdn_gpu_pointer,*rankoffup_gpu_pointer,*rankoffdn_gpu_pointer;
     int               *offdiaup_gpu,*offdiadn_gpu, *offdiaup_gpu_pointer,*offdiadn_gpu_pointer;
     int               *statesup_gpu,*statesdn_gpu,*statesup_gpu_pointer,*statesdn_gpu_pointer;
     int               *UMASK_gpu,*UMASK_gpu_pointer;
     int               *iorbup_gpu,*iorbdn_gpu,*iorbup_gpu_pointer,*iorbdn_gpu_pointer;
     cuDoubleComplex   *tmp_transp,*tmp_transp_pointer;

     if(verbose==1) printf(" GPU get eigenvalues, size of Lanczos vectors = %d ", size);

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

       cudaSetDevice(0); cudaSetDeviceFlags( cudaDeviceMapHost ); cudaEventSynchronize(0);

       cudaHostAlloc((void**)&vec_tmp_gpu   , size*sizeof(cuDoubleComplex),  cudaHostAllocMapped | cudaHostAllocPortable );
 
       if(USE_TRANSPOSE==1){ 
        if(use_texture_updo_comp==0){
          cudaHostAlloc((void**)&tmp_transp, size*sizeof(cuDoubleComplex),   cudaHostAllocMapped | cudaHostAllocPortable );
          cudaHostGetDevicePointer((void**)  &tmp_transp_pointer,tmp_transp,0);
         }else{
          cudaMalloc((void**)&tmp_transp_pointer,size*sizeof(cuDoubleComplex) );
          bind_x_(tmp_transp_pointer,size);
         };
       };
 
       if(use_texture_updo_comp==1){
         cudaMalloc((void**) &vec_in_gpu_pointer,sizeof(cuDoubleComplex)*size);
         bind_x(vec_in_gpu_pointer,size);
       }else{
         cudaHostAlloc((void**)&vec_in_gpu , size*sizeof(cuDoubleComplex),   cudaHostAllocMapped | cudaHostAllocPortable );
         cudaHostGetDevicePointer((void**)  &vec_in_gpu_pointer  ,  vec_in_gpu    , 0 );
       }

       //-----------------------------------------------------------------------------------------------//
       //-----------------------------------------------------------------------------------------------//
        cudaHostAlloc((void**)&statesup_gpu    , nup *sizeof(int),           cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&statesdn_gpu    , ndn *sizeof(int),           cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&diagdn_gpu    , ndn*sizeof(double),           cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&offdiagdn_gpu , *offdiasizedn*sizeof(cuDoubleComplex), cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&noffdn_gpu    , ndn*sizeof(int)   ,           cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&rankoffdn_gpu , *roffdn*sizeof(int)  ,        cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&diagup_gpu    , nup*sizeof(double),           cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&offdiagup_gpu , *offdiasizeup*sizeof(cuDoubleComplex), cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&noffup_gpu    , nup*sizeof(int)   ,           cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&rankoffup_gpu , *roffup*sizeof(int)  ,        cudaHostAllocMapped | cudaHostAllocPortable );

        cudaMemcpy(statesup_gpu,  statesup,     nup * sizeof(int),           cudaMemcpyHostToHost);
        cudaMemcpy(statesdn_gpu,  statesdn,     ndn * sizeof(int),           cudaMemcpyHostToHost);
        cudaMemcpy(diagup_gpu,    diagup,       nup*sizeof(double),          cudaMemcpyHostToHost);
        cudaMemcpy(offdiagup_gpu, offdiagup,   *offdiasizeup*sizeof(cuDoubleComplex), cudaMemcpyHostToHost);
        cudaMemcpy(noffup_gpu,    noffup ,      nup*sizeof(int),             cudaMemcpyHostToHost);
        cudaMemcpy(rankoffup_gpu, rankoffup,   *roffup*sizeof(int),          cudaMemcpyHostToHost);
        cudaMemcpy(diagdn_gpu,    diagdn,       ndn*sizeof(double),          cudaMemcpyHostToHost);
        cudaMemcpy(offdiagdn_gpu, offdiagdn,   *offdiasizedn*sizeof(cuDoubleComplex), cudaMemcpyHostToHost);
        cudaMemcpy(noffdn_gpu,    noffdn ,      ndn*sizeof(int),             cudaMemcpyHostToHost);
        cudaMemcpy(rankoffdn_gpu, rankoffdn,   *roffdn*sizeof(int),          cudaMemcpyHostToHost);
        cudaEventSynchronize(0);

        cudaHostGetDevicePointer((void**)  &statesup_gpu_pointer , statesup_gpu  , 0 );
        cudaHostGetDevicePointer((void**)  &statesdn_gpu_pointer , statesdn_gpu  , 0 );
        cudaHostGetDevicePointer((void**)  &diagdn_gpu_pointer ,   diagdn_gpu    , 0 );
        cudaHostGetDevicePointer((void**)  &offdiagdn_gpu_pointer, offdiagdn_gpu , 0 );
        cudaHostGetDevicePointer((void**)  &noffdn_gpu_pointer,    noffdn_gpu    , 0 );
        cudaHostGetDevicePointer((void**)  &rankoffdn_gpu_pointer, rankoffdn_gpu , 0 );
        cudaHostGetDevicePointer((void**)  &diagup_gpu_pointer ,   diagup_gpu    , 0 );
        cudaHostGetDevicePointer((void**)  &offdiagup_gpu_pointer, offdiagup_gpu , 0 );
        cudaHostGetDevicePointer((void**)  &noffup_gpu_pointer,    noffup_gpu    , 0 );
        cudaHostGetDevicePointer((void**)  &rankoffup_gpu_pointer, rankoffup_gpu , 0 );
        cudaEventSynchronize(0);
       //-----------------------------------------------------------------------------------------------//
       //-----------------------------------------------------------------------------------------------//

       cudaHostAlloc((void**)&offdiaup_gpu  , nup*sizeof(int)   ,   cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&offdiadn_gpu  , ndn*sizeof(int)   ,   cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostGetDevicePointer((void**)  &offdiadn_gpu_pointer,    offdiadn_gpu, 0 );
       cudaHostGetDevicePointer((void**)  &offdiaup_gpu_pointer,    offdiaup_gpu, 0 );
       cudaEventSynchronize(0);
       offdiaup_gpu[0]=0; for(int istate=1; istate<nup; istate++) { offdiaup_gpu[istate]=offdiaup_gpu[istate-1]+noffup[istate-1]; };
       offdiadn_gpu[0]=0; for(int istate=1; istate<ndn; istate++) { offdiadn_gpu[istate]=offdiadn_gpu[istate-1]+noffdn[istate-1]; };
       cudaEventSynchronize(0);
       cudaEventSynchronize(0);

       cudaHostAlloc((void**)&iorbup_gpu     , norb*sizeof(int),            cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&iorbdn_gpu     , norb*sizeof(int),            cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&QUART_gpu     , norb*norb*sizeof(double),     cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&UMASK_gpu,      norb*norb*sizeof(int),        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaMemcpy(iorbup_gpu,     iorbup,      norb*sizeof(int),            cudaMemcpyHostToHost);
       cudaMemcpy(iorbdn_gpu,     iorbdn,      norb*sizeof(int),            cudaMemcpyHostToHost);
       cudaMemcpy(QUART_gpu,     QUART,      norb*norb*sizeof(double),      cudaMemcpyHostToHost);
       cudaMemcpy(UMASK_gpu,     UMASK,      norb*norb*sizeof(int),         cudaMemcpyHostToHost);
       cudaHostGetDevicePointer((void**)  &iorbup_gpu_pointer ,    iorbup_gpu   , 0 );
       cudaHostGetDevicePointer((void**)  &iorbdn_gpu_pointer ,    iorbdn_gpu   , 0 );
       cudaHostGetDevicePointer((void**)  &QUART_gpu_pointer ,    QUART_gpu     , 0 );
       cudaHostGetDevicePointer((void**)  &UMASK_gpu_pointer ,   UMASK_gpu      , 0 );

       cudaHostAlloc((void**)&vec_out_gpu , size*sizeof(cuDoubleComplex),  cudaHostAllocMapped | cudaHostAllocPortable );
       cudaEventSynchronize(0);
       cudaHostGetDevicePointer((void**)  &vec_tmp_gpu_pointer ,  vec_tmp_gpu   , 0 );
       cudaHostGetDevicePointer((void**)  &vec_out_gpu_pointer ,  vec_out_gpu   , 0 );
       cudaEventSynchronize(0);

 //---------------------------------------------------------------------------------------//

  if(verbose==1) printf( " initialized, now run actual Lanczos \n " );

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

  cudaThreadSynchronize(); cudaEventSynchronize(0); 
  if(use_texture_updo_comp==0){
   cudaMemcpy(vec_in_gpu,vecinit,size*sizeof(cuDoubleComplex),cudaMemcpyHostToHost);
  }else{
   cudaMemcpy(vec_out_gpu,vecinit,size*sizeof(cuDoubleComplex),cudaMemcpyHostToHost);
   cudaMemcpy(vec_in_gpu_pointer,vec_out_gpu,size*sizeof(cuDoubleComplex),cudaMemcpyHostToDevice);
  }
  cudaThreadSynchronize(); cudaEventSynchronize(0);

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//
  for(int iter=0;iter<Niter_lanczos;iter++){
   if(verbose==1) printf( " Lanczos iteration %d / %d \n", iter,Niter_lanczos);
   one_step_lanczos_cuda(nup,ndn,norb,blocksize,Niter_lanczos,size,iter,diag,subdiag,vec_tmp_gpu,vec_in_gpu,vec_out_gpu,
               vec_tmp_gpu_pointer,vec_in_gpu_pointer,vec_out_gpu_pointer,QUART_gpu_pointer,
               diagup_gpu_pointer,diagdn_gpu_pointer,noffup_gpu_pointer,noffdn_gpu_pointer,
               rankoffup_gpu_pointer,rankoffdn_gpu_pointer,offdiagup_gpu_pointer,offdiagdn_gpu_pointer,
               offdiaup_gpu_pointer,offdiadn_gpu_pointer,statesup_gpu_pointer,statesdn_gpu_pointer,UMASK_gpu_pointer
              ,iorbup_gpu_pointer,iorbdn_gpu_pointer,tmp_transp,tmp_transp_pointer);
  };
 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

     cudaEventSynchronize(0);
     cudaFreeHost(vec_tmp_gpu);

     if(use_texture_updo_comp==0){
       cudaFreeHost(vec_in_gpu);
       if(USE_TRANSPOSE==1){ cudaFreeHost(tmp_transp);}
     }else{
       cudaFree(vec_in_gpu_pointer); 
       if(USE_TRANSPOSE==1){
         cudaFree(tmp_transp_pointer); unbind_x_();
       };
     }

     cudaFreeHost(vec_out_gpu);  cudaFreeHost(iorbup_gpu);   cudaFreeHost(iorbdn_gpu); 
     cudaFreeHost(QUART_gpu);    cudaFreeHost(UMASK_gpu);
     cudaFreeHost(offdiaup_gpu); cudaFreeHost(offdiadn_gpu);
     cudaFreeHost(statesup_gpu);  cudaFreeHost(statesdn_gpu);  cudaFreeHost(diagup_gpu);  
     cudaFreeHost(rankoffup_gpu);
     cudaFreeHost(offdiagup_gpu); cudaFreeHost(noffup_gpu);    cudaFreeHost(diagdn_gpu);
     cudaFreeHost(rankoffdn_gpu); cudaFreeHost(offdiagdn_gpu); cudaFreeHost(noffdn_gpu);

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
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************

extern "C" void lanczos_complex_updo_gs_cuda_(int *norbs, int *pblocksize, 
   int *Niter_lanczos_,int *offdiasizeup, int* offdiasizedn, int *roffup, int *roffdn, int *ntot,  int *sizup, int *sizdn, 
   double *QUART, double *diagup, double *diagdn, int *noffup, int *noffdn, int *rankoffup, int *rankoffdn, 
   cuDoubleComplex *offdiagup, cuDoubleComplex *offdiagdn,  int *UMASK, int *statesup, int *statesdn, 
   int *iorbup, int *iorbdn, double *vecp, double *GS)
{

     double diag[*Niter_lanczos_],subdiag[*Niter_lanczos_];

 //---------------------------------------------------------------------------------------//
     int verbose=0;
 //---------------------------------------------------------------------------------------//

     int Niter_lanczos= *Niter_lanczos_; int nup=*sizup; int ndn=*sizdn; int norb=*norbs;

     if(verbose==1) printf(" start Lanczos Complex on GPU \n" );

     int size = *ntot; int blocksize= *pblocksize * WARP_SIZE ; 

     cuDoubleComplex   *vec_in_gpu,*vec_out_gpu;
     double            *QUART_gpu,*diagup_gpu,*diagdn_gpu;
     cuDoubleComplex   *offdiagup_gpu,*offdiagdn_gpu,*vec_tmp_gpu;
     int               *noffup_gpu,*noffdn_gpu,*rankoffup_gpu,*rankoffdn_gpu;
     cuDoubleComplex   *vec_in_gpu_pointer,*vec_out_gpu_pointer,*vec_tmp_gpu_pointer;
     double            *QUART_gpu_pointer,*diagup_gpu_pointer,*diagdn_gpu_pointer;
     cuDoubleComplex   *offdiagup_gpu_pointer,*offdiagdn_gpu_pointer;
     int               *noffup_gpu_pointer,*noffdn_gpu_pointer,*rankoffup_gpu_pointer,*rankoffdn_gpu_pointer;
     int               *offdiaup_gpu,*offdiadn_gpu, *offdiaup_gpu_pointer,*offdiadn_gpu_pointer;
     int               *statesup_gpu,*statesdn_gpu,*statesup_gpu_pointer,*statesdn_gpu_pointer;
     int               *UMASK_gpu,*UMASK_gpu_pointer;
     int               *iorbup_gpu,*iorbdn_gpu,*iorbup_gpu_pointer,*iorbdn_gpu_pointer;
     cuDoubleComplex   *tmp_transp,*tmp_transp_pointer;

     cuDoubleComplex   *GS_gpu, *GS_gpu_pointer;

     if(verbose==1) printf(" GPU get Ground State, size of Lanczos vectors = %d ", size);

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

       cudaSetDevice(0); cudaSetDeviceFlags( cudaDeviceMapHost ); cudaEventSynchronize(0);

       cudaHostAlloc((void**)&vec_tmp_gpu   , size*sizeof(cuDoubleComplex),  cudaHostAllocMapped | cudaHostAllocPortable );

       cudaHostAlloc((void**)&GS_gpu        , size*sizeof(cuDoubleComplex),        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostGetDevicePointer((void**)   &GS_gpu_pointer ,       GS_gpu         , 0 );
 
       if(USE_TRANSPOSE==1){ 
        if(use_texture_updo_comp==0){
          cudaHostAlloc((void**)&tmp_transp, size*sizeof(cuDoubleComplex),   cudaHostAllocMapped | cudaHostAllocPortable );
          cudaHostGetDevicePointer((void**)  &tmp_transp_pointer,tmp_transp,0);
         }else{
          cudaMalloc((void**)&tmp_transp_pointer,size*sizeof(cuDoubleComplex) );
          bind_x_(tmp_transp_pointer,size);
         };
       };
 
       if(use_texture_updo_comp==1){
         cudaMalloc((void**)&vec_in_gpu_pointer,sizeof(cuDoubleComplex)*size);
         bind_x(vec_in_gpu_pointer,size);
       }else{
         cudaHostAlloc((void**)&vec_in_gpu , size*sizeof(cuDoubleComplex),   cudaHostAllocMapped | cudaHostAllocPortable );
         cudaHostGetDevicePointer((void**)  &vec_in_gpu_pointer  ,  vec_in_gpu    , 0 );
       }

       //-----------------------------------------------------------------------------------------------//
       //-----------------------------------------------------------------------------------------------//
        cudaHostAlloc((void**)&statesup_gpu  , nup *sizeof(int),           cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&statesdn_gpu  , ndn *sizeof(int),           cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&diagdn_gpu    , ndn*sizeof(double),           cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&offdiagdn_gpu , *offdiasizedn*sizeof(cuDoubleComplex), cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&noffdn_gpu    , ndn*sizeof(int)   ,           cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&rankoffdn_gpu , *roffdn*sizeof(int)  ,        cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&diagup_gpu    , nup*sizeof(double),           cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&offdiagup_gpu , *offdiasizeup*sizeof(cuDoubleComplex), cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&noffup_gpu    , nup*sizeof(int)   ,           cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&rankoffup_gpu , *roffup*sizeof(int)  ,        cudaHostAllocMapped | cudaHostAllocPortable );

        cudaMemcpy(statesup_gpu,  statesup,     nup * sizeof(int),           cudaMemcpyHostToHost);
        cudaMemcpy(statesdn_gpu,  statesdn,     ndn * sizeof(int),           cudaMemcpyHostToHost);
        cudaMemcpy(diagup_gpu,    diagup,       nup*sizeof(double),          cudaMemcpyHostToHost);
        cudaMemcpy(offdiagup_gpu, offdiagup,   *offdiasizeup*sizeof(cuDoubleComplex), cudaMemcpyHostToHost);
        cudaMemcpy(noffup_gpu,    noffup ,      nup*sizeof(int),             cudaMemcpyHostToHost);
        cudaMemcpy(rankoffup_gpu, rankoffup,   *roffup*sizeof(int),          cudaMemcpyHostToHost);
        cudaMemcpy(diagdn_gpu,    diagdn,       ndn*sizeof(double),          cudaMemcpyHostToHost);
        cudaMemcpy(offdiagdn_gpu, offdiagdn,   *offdiasizedn*sizeof(cuDoubleComplex), cudaMemcpyHostToHost);
        cudaMemcpy(noffdn_gpu,    noffdn ,      ndn*sizeof(int),             cudaMemcpyHostToHost);
        cudaMemcpy(rankoffdn_gpu, rankoffdn,   *roffdn*sizeof(int),          cudaMemcpyHostToHost);
        cudaEventSynchronize(0);

        cudaHostGetDevicePointer((void**)  &statesup_gpu_pointer , statesup_gpu  , 0 );
        cudaHostGetDevicePointer((void**)  &statesdn_gpu_pointer , statesdn_gpu  , 0 );
        cudaHostGetDevicePointer((void**)  &diagdn_gpu_pointer ,   diagdn_gpu    , 0 );
        cudaHostGetDevicePointer((void**)  &offdiagdn_gpu_pointer, offdiagdn_gpu , 0 );
        cudaHostGetDevicePointer((void**)  &noffdn_gpu_pointer,    noffdn_gpu    , 0 );
        cudaHostGetDevicePointer((void**)  &rankoffdn_gpu_pointer, rankoffdn_gpu , 0 );
        cudaHostGetDevicePointer((void**)  &diagup_gpu_pointer ,   diagup_gpu    , 0 );
        cudaHostGetDevicePointer((void**)  &offdiagup_gpu_pointer, offdiagup_gpu , 0 );
        cudaHostGetDevicePointer((void**)  &noffup_gpu_pointer,    noffup_gpu    , 0 );
        cudaHostGetDevicePointer((void**)  &rankoffup_gpu_pointer, rankoffup_gpu , 0 );
        cudaEventSynchronize(0);
       //-----------------------------------------------------------------------------------------------//
       //-----------------------------------------------------------------------------------------------//

       cudaHostAlloc((void**)&offdiaup_gpu    , nup*sizeof(int)   ,   cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&offdiadn_gpu    , ndn*sizeof(int)   ,   cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostGetDevicePointer((void**)  &offdiadn_gpu_pointer,    offdiadn_gpu, 0 );
       cudaHostGetDevicePointer((void**)  &offdiaup_gpu_pointer,    offdiaup_gpu, 0 );
       cudaEventSynchronize(0);
       offdiaup_gpu[0]=0; for(int istate=1; istate<nup; istate++) { offdiaup_gpu[istate]=offdiaup_gpu[istate-1]+noffup[istate-1]; };
       offdiadn_gpu[0]=0; for(int istate=1; istate<ndn; istate++) { offdiadn_gpu[istate]=offdiadn_gpu[istate-1]+noffdn[istate-1]; };
       cudaEventSynchronize(0);
       cudaEventSynchronize(0);

       cudaHostAlloc((void**)&iorbup_gpu     , norb*sizeof(int),            cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&iorbdn_gpu     , norb*sizeof(int),            cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&QUART_gpu     , norb*norb*sizeof(double),     cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&UMASK_gpu,      norb*norb*sizeof(int),        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaMemcpy(iorbup_gpu,     iorbup,      norb*sizeof(int),            cudaMemcpyHostToHost);
       cudaMemcpy(iorbdn_gpu,     iorbdn,      norb*sizeof(int),            cudaMemcpyHostToHost);
       cudaMemcpy(QUART_gpu,     QUART,      norb*norb*sizeof(double),      cudaMemcpyHostToHost);
       cudaMemcpy(UMASK_gpu,     UMASK,      norb*norb*sizeof(int),         cudaMemcpyHostToHost);
       cudaHostGetDevicePointer((void**)  &iorbup_gpu_pointer ,    iorbup_gpu   , 0 );
       cudaHostGetDevicePointer((void**)  &iorbdn_gpu_pointer ,    iorbdn_gpu   , 0 );
       cudaHostGetDevicePointer((void**)  &QUART_gpu_pointer ,    QUART_gpu     , 0 );
       cudaHostGetDevicePointer((void**)  &UMASK_gpu_pointer ,   UMASK_gpu      , 0 );

       cudaHostAlloc((void**)&vec_out_gpu   , size*sizeof(cuDoubleComplex),  cudaHostAllocMapped | cudaHostAllocPortable );
       cudaEventSynchronize(0);
       cudaHostGetDevicePointer((void**)  &vec_tmp_gpu_pointer ,  vec_tmp_gpu   , 0 );
       cudaHostGetDevicePointer((void**)  &vec_out_gpu_pointer ,  vec_out_gpu   , 0 );
       cudaEventSynchronize(0);

 //---------------------------------------------------------------------------------------//

  if(verbose==1) printf( " initialized, now run actual Lanczos \n " );

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

  cudaThreadSynchronize(); cudaEventSynchronize(0); 
 if(use_texture_updo_comp==0){
  for(int i=0;i<size;i++) vec_in_gpu[i]=make_cuDoubleComplex(1.0,0.0); cudaThreadSynchronize();
 }else{
   for(int i=0;i<size;i++) vec_out_gpu[i]=make_cuDoubleComplex(1.0,0.0); cudaThreadSynchronize();
   cudaMemcpy(vec_in_gpu_pointer,vec_out_gpu,size*sizeof(cuDoubleComplex),cudaMemcpyHostToDevice);
 }
   cudaThreadSynchronize(); cudaEventSynchronize(0);

  memset ((void **)GS_gpu, 0, size*sizeof(cuDoubleComplex));
  double *normv_ker; double normv; double *normv_loc; cudaMalloc((void**)&normv_ker,sizeof(double)); cuDoubleComplex coef;

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

 if(use_texture_updo_comp==0){
   coef = make_cuDoubleComplex(1.0/sqrt(cuCabs(cublasZdotu_no_device(size,vec_in_gpu,1,vec_in_gpu,1)))*vecp[0],0); cudaThreadSynchronize();  cudaEventSynchronize(0);
   cublasZaxpy_no_device(size,coef,vec_in_gpu,1,GS_gpu,1); cudaThreadSynchronize(); cudaEventSynchronize(0);cudaThreadSynchronize();
 }else{
   norm_vec_ker_complex_updo<<<1,512>>>(size,vec_in_gpu_pointer,normv_ker); cudaEventSynchronize(0); cudaThreadSynchronize();
   normv_loc=&normv; cudaMemcpy(normv_loc,normv_ker,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(vec_out_gpu,vec_in_gpu,size*sizeof(cuDoubleComplex),cudaMemcpyDeviceToHost);
   coef=make_cuDoubleComplex(vecp[0]/normv,0.0);
   cublasZaxpy_no_device(size,coef,vec_out_gpu,1,GS_gpu,1); cudaThreadSynchronize(); cudaEventSynchronize(0);cudaThreadSynchronize();
 };

   one_step_lanczos_cuda(nup,ndn,norb,blocksize,Niter_lanczos,size,0,diag,subdiag,vec_tmp_gpu,vec_in_gpu,vec_out_gpu,
               vec_tmp_gpu_pointer,vec_in_gpu_pointer,vec_out_gpu_pointer,QUART_gpu_pointer,
               diagup_gpu_pointer,diagdn_gpu_pointer,noffup_gpu_pointer,noffdn_gpu_pointer,
               rankoffup_gpu_pointer,rankoffdn_gpu_pointer,offdiagup_gpu_pointer,offdiagdn_gpu_pointer,
               offdiaup_gpu_pointer,offdiadn_gpu_pointer,statesup_gpu_pointer,statesdn_gpu_pointer,UMASK_gpu_pointer
              ,iorbup_gpu_pointer,iorbdn_gpu_pointer,tmp_transp,tmp_transp_pointer);

  for(int iter=1;iter<Niter_lanczos-1;iter++){

   coef = make_cuDoubleComplex(1.0/sqrt(cuCabs(cublasZdotu_no_device(size,vec_out_gpu,1,vec_out_gpu,1))),0.0);
   cublasZscal_no_device (size,coef, vec_out_gpu,1);  cudaEventSynchronize(0);

   coef = make_cuDoubleComplex(vecp[iter],0);
   cublasZaxpy_no_device(size,coef,vec_out_gpu,1,GS_gpu,1); cudaThreadSynchronize(); cudaEventSynchronize(0);

   one_step_lanczos_cuda(nup,ndn,norb,blocksize,Niter_lanczos,size,iter,diag,subdiag,vec_tmp_gpu,vec_in_gpu,vec_out_gpu,
               vec_tmp_gpu_pointer,vec_in_gpu_pointer,vec_out_gpu_pointer,QUART_gpu_pointer,
               diagup_gpu_pointer,diagdn_gpu_pointer,noffup_gpu_pointer,noffdn_gpu_pointer,
               rankoffup_gpu_pointer,rankoffdn_gpu_pointer, offdiagup_gpu_pointer,offdiagdn_gpu_pointer,
               offdiaup_gpu_pointer,offdiadn_gpu_pointer,statesup_gpu_pointer,statesdn_gpu_pointer,UMASK_gpu_pointer
              ,iorbup_gpu_pointer,iorbdn_gpu_pointer,tmp_transp,tmp_transp_pointer);

  };

   if(verbose==1) printf("done...\n");

   coef=make_cuDoubleComplex(1.0/sqrt(cuCabs(cublasZdotu_no_device(size,vec_out_gpu,1,vec_out_gpu,1))),0.0);
   cublasZscal_no_device(size,coef,vec_out_gpu,1);cudaEventSynchronize(0);

   coef=make_cuDoubleComplex(vecp[Niter_lanczos-1],0.0);
   cublasZaxpy_no_device(size,coef,vec_out_gpu,1,GS_gpu,1); cudaThreadSynchronize(); cudaEventSynchronize(0);

   coef=make_cuDoubleComplex(1.0/sqrt(cuCabs(cublasZdotu_no_device(size,GS_gpu,1,GS_gpu,1))),0.0);  cudaThreadSynchronize();  cudaEventSynchronize(0);
   cublasZscal_no_device(size,coef,GS_gpu,1); cudaEventSynchronize(0); cudaThreadSynchronize();

   cudaMemcpy(GS,GS_gpu,size*sizeof(cuDoubleComplex),cudaMemcpyHostToHost);

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

     cudaFree(normv_ker);

     cudaEventSynchronize(0);
     cudaFreeHost(vec_tmp_gpu);

     if(use_texture_updo_comp==0){
       cudaFreeHost(vec_in_gpu);
       if(USE_TRANSPOSE==1){ cudaFreeHost(tmp_transp);}
     }else{
       cudaFree(vec_in_gpu_pointer); 
       if(USE_TRANSPOSE==1){
         cudaFree(tmp_transp_pointer); unbind_x_();
       };
     }

     cudaFreeHost(vec_out_gpu);  cudaFreeHost(iorbup_gpu);   cudaFreeHost(iorbdn_gpu); 
     cudaFreeHost(QUART_gpu);    cudaFreeHost(UMASK_gpu);
     cudaFreeHost(offdiaup_gpu); cudaFreeHost(offdiadn_gpu);
     cudaFreeHost(statesup_gpu);  cudaFreeHost(statesdn_gpu);  cudaFreeHost(diagup_gpu);  
     cudaFreeHost(rankoffup_gpu);
     cudaFreeHost(offdiagup_gpu); cudaFreeHost(noffup_gpu);    cudaFreeHost(diagdn_gpu);
     cudaFreeHost(rankoffdn_gpu); cudaFreeHost(offdiagdn_gpu); cudaFreeHost(noffdn_gpu);

     cudaFreeHost(GS_gpu);

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
//********************************************
//********************************************
//********************************************
//********************************************

