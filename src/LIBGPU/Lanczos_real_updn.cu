
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

#define IBSET(a,b) ((a) |= (1<<(b)))
#define IBCLR(a,b) ((a) &= ~(1<<(b)))
#define BTEST(a,b) ((a) & (1<<(b)))

#define WARP_SIZE 32
#define MAX_BLOCK 65500 
#define BLOCK_DIM 16
#define USE_TRANSPOSE 1
const int use_texture_updo = 0 ;

texture<int2,1,cudaReadModeElementType> tex, texdn;


inline void   bind_x(double *x, int size) { cudaBindTexture(0,tex,x,size*sizeof(double)); };
inline void unbind_x()                    { cudaUnbindTexture(tex); };

inline void   bind_x_(double *y, int size) { cudaBindTexture(0,texdn,y,size*sizeof(double)); };
inline void unbind_x_()                    { cudaUnbindTexture(texdn);    };


__inline__    __device__ double fetch_x(const int& i)
  {   int2 v = tex1Dfetch(tex,i); return __hiloint2double(v.y, v.x); }
__inline__    __device__ double fetch_x_(const int& i)
  {   int2 v = tex1Dfetch(texdn,i); return __hiloint2double(v.y, v.x); }

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
     //---------------------------------------//
     //---------------------------------------//
     //---------------------------------------//
     //---------------------------------------//

__global__ void transpose_ker(int nb1, int nb2, int ngrid, double *output, double *data, int width, int height)
{

       __shared__ double block[BLOCK_DIM][BLOCK_DIM+1];

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

  __global__  void Hmult_ker_updo0(int nup, int ndn, int norb, int ngrid, int BLOCK_SIZE, int num_rows, double *y, const double *x,
                             const double *QUART, const double *diagup, const double *diagdn,
                             const int *noffup, const int *noffdn, const int *rankoffup, const int *rankoffdn,
                             const double *offdiagup, const double *offdiagdn, const int *offdiaup, const int *offdiadn,
                             const int *statesup, const int *statesdn, const int *UMASK,
                             const int use_texture_updo, const int *iorbup, const int *iorbdn)
{

   int istate,iup,idn;
   int iii;
   double diagdn_ ;
   const int thread_lane = threadIdx.x; const int interval = WARP_SIZE;

// DIAGONAL TERMS, COALESCENT READING

   ////////////////////////////////////////////////////////////
   for(iii=0; iii<=ngrid; iii++){
    idn = BLOCK_SIZE * (blockIdx.y+iii*MAX_BLOCK) + threadIdx.y ;
    if(idn<ndn)
    {
     diagdn_=diagdn[idn]; 
     for(iup=thread_lane; iup<nup; iup+=interval ) 
      { istate       =   iup + idn * nup;
       if(use_texture_updo==0){
        y[istate]   += ( diagup[iup] + diagdn_ ) * x[istate];
       }else{
        y[istate]   += ( diagup[iup] + diagdn_ ) * fetch_x(istate);
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

  __global__  void Hmult_ker_updo1(int nup, int ndn, int norb, int ngrid, int BLOCK_SIZE, int num_rows, double *y, const double *x, 
                             const double *QUART, const double *diagup, const double *diagdn, const int *noffup, const int *noffdn,  
                             const int *rankoffup, const int *rankoffdn, const double *offdiagup, const double *offdiagdn, 
                             const int *offdiaup, const int *offdiadn, const int *statesup, const int *statesdn, const int *UMASK,
                             const int use_texture_updo, const int *iorbup, const int *iorbdn)
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
            if(use_texture_updo==1){
             y[istate] += quart_ * fetch_x(istate); 
            }else{
             y[istate] += quart_ * x[istate];
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

  __global__  void Hmult_ker_updo2(int nup, int ndn, int norb, int ngrid, int BLOCK_SIZE, int num_rows, double *y, const double *x, 
                             const double *QUART, const double *diagup, const double *diagdn,
                             const int *noffup, const int *noffdn,
                             const int *rankoffup, const int *rankoffdn,
                             const double *offdiagup, const double *offdiagdn,
                             const int *offdiaup, const int *offdiadn,
                             const int *statesup, const int *statesdn,
                             const int *UMASK, const int use_texture_updo, const int *iorbup, const int *iorbdn)
{
   int jstate,iup ;
   int istatemin,istatemax,noff,nhoffdiag,ii,iii,irank,jup;
   double hoffdiag ;

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
       hoffdiag    =  offdiagup[nhoffdiag+irank];
       jup         =  rankoffup[nhoffdiag+irank]-1;
       jstate      =  jup  ;
       for(ii=istatemin;ii<=istatemax;ii+=nup)
        {
        if(use_texture_updo==1){
         y[ii] += hoffdiag * fetch_x(jstate);
        }else{
         y[ii] += hoffdiag * x[jstate];
        };
         jstate  +=   nup;
        };
      };
   ////////////////////////////////////////////////////////////

   }; };

 }

     //---------------------------------------//
     //---------------------------------------//
     //---------------------------------------//
     //---------------------------------------//

  __global__  void Hmult_ker_updo3(int nup, int ndn, int norb, int ngrid, int BLOCK_SIZE, int num_rows, double *y, const double *x, 
                             const double *QUART, const double *diagup, const double *diagdn,
                             const int *noffup, const int *noffdn, const int *rankoffup, const int *rankoffdn,
                             const double *offdiagup, const double *offdiagdn, const int *offdiaup, const int *offdiadn,
                             const int *statesup, const int *statesdn, const int *UMASK, const int use_texture_updo, const int *iorbup, const int *iorbdn)
{

   int jstate,idn ; 
   int istatemin,istatemax,noff,nhoffdiag,ii,iii,irank;
   const int thread_lane = threadIdx.x;
   double hoffdiag;
   int jdn; 

// HOPPING DN

   ////////////////////////////////////////////////////////////
   for(iii=0; iii<=ngrid; iii++){
     idn = BLOCK_SIZE * (blockIdx.y+iii*MAX_BLOCK) + threadIdx.y ;
     if(idn<ndn)
    {
      istatemin = idn*nup;
      istatemax = idn*nup+nup-1;
      noff      = noffdn[idn];
      nhoffdiag = offdiadn[idn];
      for(irank=0;irank<noff;irank++)
       {
       hoffdiag =  offdiagdn[nhoffdiag+irank];
       jdn      =  rankoffdn[nhoffdiag+irank]-1;
       jstate   =  jdn*nup+thread_lane ;
       for(ii=istatemin+thread_lane;ii<=istatemax;ii+=WARP_SIZE)
        { 
         if(use_texture_updo==0){
          y[ii]  += hoffdiag * x[jstate]; 
         }else{
          if(use_texture_updo==1){
            y[ii]  += hoffdiag * fetch_x(jstate);
           }else{
            y[ii]  += hoffdiag * fetch_x_(jstate);
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
 double *subdiag, double *vec_tmp_gpu, double *vec_in_gpu, double *vec_out_gpu, double *vec_tmp_gpu_pointer, 
 double *vec_in_gpu_pointer, double *vec_out_gpu_pointer, double *QUART_gpu_pointer, double *diagup_gpu_pointer, double *diagdn_gpu_pointer,
 int *noffup_gpu_pointer, int *noffdn_gpu_pointer, int *rankoffup_gpu_pointer, int *rankoffdn_gpu_pointer, 
 double *offdiagup_gpu_pointer,double *offdiagdn_gpu_pointer, int *offdiaup_gpu_pointer, int *offdiadn_gpu_pointer,
 int *statesup_gpu_pointer,int *statesdn_gpu_pointer, int *UMASK_gpu_pointer,int *iorbup_gpu_pointer,int *iorbdn_gpu_pointer,double *tmp_transp,
 double *tmp_transp_pointer)

{
   int verbose=0 ; int psize=size; int nb, ngrid,nb2,nb1;

   if(verbose==1) printf ( " BLOCKSIZE = %d \n ", blocksize );
   if(verbose==1) printf ( " Sdot, vec in norm \n ");
   double normv = sqrt(cublasDdot(size,vec_in_gpu_pointer,1,vec_in_gpu_pointer,1)); cudaEventSynchronize(0); cudaThreadSynchronize();
   if(verbose==1) printf( " norm=%f \n ", normv);

   cublasDscal(size,1.0/normv,vec_in_gpu_pointer,1); cudaEventSynchronize(0); cudaThreadSynchronize(); 


   //------------------------------------------------------------------------------------------------------------------//

   if(verbose==1) printf( " reset memory \n ");
   memset ((void **)vec_out_gpu, 0,size*sizeof(double));
   cudaThreadSynchronize();  cudaEventSynchronize(0);

   //------------------------------------------------------------------------------------------------------------------//
   //------------------------------------------------------------------------------------------------------------------//
  if(USE_TRANSPOSE==1)
   {

   nb1=(nup-nup % BLOCK_DIM)/BLOCK_DIM+1; nb2=(ndn-ndn % BLOCK_DIM)/BLOCK_DIM+1; nb=nb1*nb2; ngrid=nb/MAX_BLOCK; if(ngrid>0)nb=MAX_BLOCK;
   dim3 threads(BLOCK_DIM, BLOCK_DIM, 1);
   transpose_ker<<<nb,threads>>>(nb1,nb2,ngrid,tmp_transp_pointer,vec_in_gpu_pointer,nup,ndn);
   cudaEventSynchronize(0); cudaThreadSynchronize();

   int bb=blocksize/WARP_SIZE; if(bb<1) bb=1; nb=(nup-nup % bb)/bb+1; ngrid=nb/MAX_BLOCK; if(ngrid>0)nb=MAX_BLOCK; dim3 t1(1,nb),t2(WARP_SIZE,bb);
   Hmult_ker_updo3<<<t1,t2>>>(ndn,nup,norb,ngrid,bb,psize,vec_out_gpu_pointer,tmp_transp_pointer,QUART_gpu_pointer,
                         diagdn_gpu_pointer,diagup_gpu_pointer,noffdn_gpu_pointer,noffup_gpu_pointer,rankoffdn_gpu_pointer,rankoffup_gpu_pointer,
                         offdiagdn_gpu_pointer,offdiagup_gpu_pointer,offdiadn_gpu_pointer,offdiaup_gpu_pointer,
                         statesdn_gpu_pointer,statesup_gpu_pointer,UMASK_gpu_pointer,-use_texture_updo,iorbdn_gpu_pointer,iorbup_gpu_pointer);
   if(verbose==1) printf( " done.... \n " );
   cudaEventSynchronize(0); cudaThreadSynchronize();

   nb2=(nup-nup % BLOCK_DIM)/BLOCK_DIM+1; nb1=(ndn-ndn % BLOCK_DIM)/BLOCK_DIM+1; nb=nb1*nb2; ngrid=nb/MAX_BLOCK; if(ngrid>0)nb=MAX_BLOCK;
   transpose_ker<<<nb,threads>>>(nb1,nb2,ngrid,vec_in_gpu_pointer,tmp_transp_pointer,ndn,nup);
   cudaEventSynchronize(0); cudaThreadSynchronize();

   transpose_ker<<<nb,threads>>>(nb1,nb2,ngrid,tmp_transp_pointer,vec_out_gpu_pointer,ndn,nup);
   cudaEventSynchronize(0);cudaThreadSynchronize();
 
   if(use_texture_updo==0){
    cudaMemcpy(vec_out_gpu,tmp_transp,size*sizeof(double),cudaMemcpyHostToHost);cudaEventSynchronize(0);cudaThreadSynchronize();
   }else{
    cudaMemcpy(vec_out_gpu,tmp_transp_pointer,size*sizeof(double),cudaMemcpyDeviceToHost);cudaEventSynchronize(0);cudaThreadSynchronize();
   }

   cudaEventSynchronize(0);cudaThreadSynchronize();
   if(verbose==1) printf( "transposition applied....\n");

   //------------------------------------------------------------------------------------------------------------------//
  }else{
   nb=(nup-nup % blocksize)/blocksize+1; ngrid=nb/MAX_BLOCK; if(ngrid>0)nb=MAX_BLOCK;
   Hmult_ker_updo2<<<nb,blocksize>>>(nup,ndn,norb,ngrid,blocksize,psize,vec_out_gpu_pointer,vec_in_gpu_pointer,QUART_gpu_pointer,
                         diagup_gpu_pointer,diagdn_gpu_pointer,noffup_gpu_pointer,noffdn_gpu_pointer,rankoffup_gpu_pointer,rankoffdn_gpu_pointer,
                         offdiagup_gpu_pointer,offdiagdn_gpu_pointer,offdiaup_gpu_pointer,offdiadn_gpu_pointer,
                         statesup_gpu_pointer,statesdn_gpu_pointer,UMASK_gpu_pointer,use_texture_updo
                        ,iorbup_gpu_pointer,iorbdn_gpu_pointer);
   if(verbose==1) printf( " done.... \n " );
   cudaEventSynchronize(0); cudaThreadSynchronize();
  };
   //------------------------------------------------------------------------------------------------------------------//
   //-------------------------------------------------------------------------------------------------------------------//

   int bb=blocksize/WARP_SIZE; if(bb<1) bb=1; nb=(ndn-ndn % bb)/bb+1; ngrid=nb/MAX_BLOCK; if(ngrid>0)nb=MAX_BLOCK; dim3 t1(1,nb), t2(WARP_SIZE,bb);
   Hmult_ker_updo3<<<t1,t2>>>(nup,ndn,norb,ngrid,bb,psize,vec_out_gpu_pointer,vec_in_gpu_pointer,QUART_gpu_pointer,
                         diagup_gpu_pointer,diagdn_gpu_pointer,noffup_gpu_pointer,noffdn_gpu_pointer,rankoffup_gpu_pointer,rankoffdn_gpu_pointer,
                         offdiagup_gpu_pointer,offdiagdn_gpu_pointer,offdiaup_gpu_pointer,offdiadn_gpu_pointer,
                         statesup_gpu_pointer,statesdn_gpu_pointer,UMASK_gpu_pointer,use_texture_updo,iorbup_gpu_pointer,iorbdn_gpu_pointer);
   if(verbose==1) printf( " done.... \n " );
   cudaEventSynchronize(0); cudaThreadSynchronize();
   //------------------------------------------------------------------------------------------------------------------//
   if(verbose==1) printf(" start diag terms...\n");
   bb=blocksize/WARP_SIZE; if(bb<1) bb=1; nb=(ndn- ndn % bb)/bb+1; ngrid=nb/MAX_BLOCK; if(ngrid>0)nb=MAX_BLOCK; dim3 tt1(1,nb), tt2(WARP_SIZE,bb); 
   if(verbose==1) printf(" now kernel....\n");
   Hmult_ker_updo0<<<tt1,tt2>>>(nup,ndn,norb,ngrid,bb,psize,
                         vec_out_gpu_pointer,vec_in_gpu_pointer,QUART_gpu_pointer, diagup_gpu_pointer,diagdn_gpu_pointer,
                         noffup_gpu_pointer,noffdn_gpu_pointer,rankoffup_gpu_pointer,rankoffdn_gpu_pointer,
                         offdiagup_gpu_pointer,offdiagdn_gpu_pointer,offdiaup_gpu_pointer,offdiadn_gpu_pointer,
                         statesup_gpu_pointer,statesdn_gpu_pointer,UMASK_gpu_pointer,use_texture_updo
                        ,iorbup_gpu_pointer,iorbdn_gpu_pointer);
   if(verbose==1) printf( " done.... \n " );
   cudaEventSynchronize(0); cudaThreadSynchronize();
   //------------------------------------------------------------------------------------------------------------------//
   nb=(nup-nup % blocksize)/blocksize+1; ngrid=nb/MAX_BLOCK; if(ngrid>0)nb=MAX_BLOCK; 
   if(verbose==1) printf( " call kernel ... nblock = %d ngrid=%d norb=%d blocksize=%d \n ",nb, ngrid, norb, blocksize );
   Hmult_ker_updo1<<<nb,blocksize>>>(nup,ndn,norb,ngrid,blocksize,psize,vec_out_gpu_pointer,vec_in_gpu_pointer,QUART_gpu_pointer, 
                         diagup_gpu_pointer,diagdn_gpu_pointer, noffup_gpu_pointer,noffdn_gpu_pointer,rankoffup_gpu_pointer,rankoffdn_gpu_pointer,
                         offdiagup_gpu_pointer,offdiagdn_gpu_pointer,offdiaup_gpu_pointer,offdiadn_gpu_pointer,
                         statesup_gpu_pointer,statesdn_gpu_pointer,UMASK_gpu_pointer,use_texture_updo
                        ,iorbup_gpu_pointer,iorbdn_gpu_pointer); 
   if(verbose==1) printf( " done.... \n " );
   cudaEventSynchronize(0); cudaThreadSynchronize();
   //------------------------------------------------------------------------------------------------------------------//

  if(verbose==1) printf( " scalar product, get subdiag \n " );
 
  if(iter>0){cublasDaxpy(size,-subdiag[iter],vec_tmp_gpu_pointer,1,vec_out_gpu_pointer,1);}; cudaEventSynchronize(0); cudaThreadSynchronize();

  if(verbose==1) printf( "get vec tmp, subdiag[%d]=%F  \n ", iter, subdiag[iter]);

  if(use_texture_updo==0){
   cudaMemcpy(vec_tmp_gpu,vec_in_gpu,size*sizeof(double),cudaMemcpyHostToHost);cudaEventSynchronize(0);cudaThreadSynchronize();
  }else{
   cudaMemcpy(vec_tmp_gpu,vec_in_gpu_pointer,size*sizeof(double),cudaMemcpyDeviceToHost);cudaEventSynchronize(0);cudaThreadSynchronize();
  }

   diag[iter]=cublasDdot(size, vec_out_gpu_pointer,1,vec_in_gpu_pointer,1);cudaEventSynchronize(0);cudaThreadSynchronize();
   if(verbose==1) printf( " get diag iter , diag[%d]=%f \n ", iter,diag[iter]);
   cublasDaxpy(size,-diag[iter],vec_tmp_gpu_pointer,1,vec_out_gpu_pointer,1); cudaEventSynchronize(0); cudaThreadSynchronize();
   normv = sqrt(cublasDdot(size, vec_out_gpu_pointer,1,vec_out_gpu_pointer,1)); cudaEventSynchronize(0); cudaThreadSynchronize();

   if(iter<Niter-1) subdiag[iter+1]=normv;

   if(verbose==1) printf( " copy back vec_out_gpu to vec_int_gpu , normv=%f \n ",normv );

   if(use_texture_updo==0){
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

extern "C" void lanczos_real_updo_cuda_(int *norbs, int *pblocksize, int *Niter_lanczos_,
   int *offdiasizeup, int* offdiasizedn, 
   int *roffup, int *roffdn, int *ntot,  int *sizup, int *sizdn, 
   double *QUART, double *diagup, double *diagdn, int *noffup, int *noffdn, int *rankoffup, int *rankoffdn, 
   double *offdiagup, double *offdiagdn, double *diag, double *subdiag, int *UMASK, int *statesup, int *statesdn , int *iorbup, int *iorbdn)
{

 //---------------------------------------------------------------------------------------//
     int verbose=0;
 //---------------------------------------------------------------------------------------//

     int Niter_lanczos= *Niter_lanczos_; int nup=*sizup; int ndn=*sizdn; int norb=*norbs;

     if(verbose==1) printf(" start Lanczos Real on GPU \n" );
     if(verbose==1) printf(" norbitals=%d \n ", norb );

     int size = *ntot; int blocksize= *pblocksize * WARP_SIZE ; 

     double   *vec_in_gpu,*vec_out_gpu,*QUART_gpu,*diagup_gpu,*diagdn_gpu,*offdiagup_gpu,*offdiagdn_gpu,*vec_tmp_gpu;
     int      *noffup_gpu,*noffdn_gpu,*rankoffup_gpu,*rankoffdn_gpu;
     double   *vec_in_gpu_pointer,*vec_out_gpu_pointer,*vec_tmp_gpu_pointer;
     double   *QUART_gpu_pointer,*diagup_gpu_pointer,*diagdn_gpu_pointer,*offdiagup_gpu_pointer,*offdiagdn_gpu_pointer;
     int      *noffup_gpu_pointer,*noffdn_gpu_pointer,*rankoffup_gpu_pointer,*rankoffdn_gpu_pointer;
     int      *offdiaup_gpu,*offdiadn_gpu, *offdiaup_gpu_pointer,*offdiadn_gpu_pointer;
     int      *statesup_gpu,*statesdn_gpu,*statesup_gpu_pointer,*statesdn_gpu_pointer;
     int      *UMASK_gpu,*UMASK_gpu_pointer;
     int      *iorbup_gpu,*iorbdn_gpu,*iorbup_gpu_pointer,*iorbdn_gpu_pointer;
     double   *tmp_transp,*tmp_transp_pointer;

     if(verbose==1) printf(" GPU get eigenvalues, size of Lanczos vectors = %d ", size);
     if(verbose==1) printf(" use texture = %d \n ", use_texture_updo);
 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

       cudaSetDevice(0); cudaSetDeviceFlags( cudaDeviceMapHost ); cudaEventSynchronize(0);

       cudaHostAlloc((void**)&vec_tmp_gpu   , size*sizeof(double),  cudaHostAllocMapped | cudaHostAllocPortable );
 
       if(USE_TRANSPOSE==1){ 
        if(use_texture_updo==0){
          cudaHostAlloc((void**)&tmp_transp, size*sizeof(double),   cudaHostAllocMapped | cudaHostAllocPortable );
          cudaHostGetDevicePointer((void**)  &tmp_transp_pointer,tmp_transp,0);
         }else{
          cudaMalloc((void**)&tmp_transp_pointer,size*sizeof(double) );
          bind_x_(tmp_transp_pointer,size);
         };
       };
 
       if(use_texture_updo==1){
         cudaMalloc((void**)&vec_in_gpu_pointer,sizeof(double)*size);
         bind_x(vec_in_gpu_pointer,size);
       }else{
         cudaHostAlloc((void**)&vec_in_gpu , size*sizeof(double),   cudaHostAllocMapped | cudaHostAllocPortable );
         cudaHostGetDevicePointer((void**)  &vec_in_gpu_pointer  ,  vec_in_gpu    , 0 );
       }

       //-----------------------------------------------------------------------------------------------//
       //-----------------------------------------------------------------------------------------------//
        cudaHostAlloc((void**)&statesup_gpu    , nup *sizeof(int),           cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&statesdn_gpu    , ndn *sizeof(int),           cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&diagdn_gpu    , ndn*sizeof(double),           cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&offdiagdn_gpu , *offdiasizedn*sizeof(double), cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&noffdn_gpu    , ndn*sizeof(int)   ,           cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&rankoffdn_gpu , *roffdn*sizeof(int)  ,        cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&diagup_gpu    , nup*sizeof(double),           cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&offdiagup_gpu , *offdiasizeup*sizeof(double), cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&noffup_gpu    , nup*sizeof(int)   ,           cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&rankoffup_gpu , *roffup*sizeof(int)  ,        cudaHostAllocMapped | cudaHostAllocPortable );

        cudaMemcpy(statesup_gpu,  statesup,     nup * sizeof(int),           cudaMemcpyHostToHost);
        cudaMemcpy(statesdn_gpu,  statesdn,     ndn * sizeof(int),           cudaMemcpyHostToHost);
        cudaMemcpy(diagup_gpu,    diagup,       nup*sizeof(double),          cudaMemcpyHostToHost);
        cudaMemcpy(offdiagup_gpu, offdiagup,   *offdiasizeup*sizeof(double), cudaMemcpyHostToHost);
        cudaMemcpy(noffup_gpu,    noffup ,      nup*sizeof(int),             cudaMemcpyHostToHost);
        cudaMemcpy(rankoffup_gpu, rankoffup,   *roffup*sizeof(int),          cudaMemcpyHostToHost);
        cudaMemcpy(diagdn_gpu,    diagdn,       ndn*sizeof(double),          cudaMemcpyHostToHost);
        cudaMemcpy(offdiagdn_gpu, offdiagdn,   *offdiasizedn*sizeof(double), cudaMemcpyHostToHost);
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
       if(verbose==1){
       printf( " \n number of up states = %d \n ", nup);
       printf( " \n number of dn states = %d \n ", ndn);
       printf( " \n total number of off-diag terms up=%d dn=%d \n " , offdiaup_gpu[nup-1],offdiadn_gpu[ndn-1]);
       };
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

       cudaHostAlloc((void**)&vec_out_gpu   , size*sizeof(double),        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaEventSynchronize(0);
       cudaHostGetDevicePointer((void**)  &vec_tmp_gpu_pointer ,  vec_tmp_gpu   , 0 );
       cudaHostGetDevicePointer((void**)  &vec_out_gpu_pointer ,  vec_out_gpu   , 0 );
       cudaEventSynchronize(0);

 //---------------------------------------------------------------------------------------//

  if(verbose==1) printf( " initialized, now run actual Lanczos \n " );

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

  cudaThreadSynchronize(); cudaEventSynchronize(0); 
 if(use_texture_updo==0){
  for(int i=0;i<size;i++) vec_in_gpu[i]=1.0; cudaThreadSynchronize();
 }else{
   for(int i=0;i<size;i++) vec_out_gpu[i]=1.0; cudaThreadSynchronize();
   cudaMemcpy(vec_in_gpu_pointer,vec_out_gpu,size*sizeof(double),cudaMemcpyHostToDevice);
 }
   cudaThreadSynchronize(); cudaEventSynchronize(0);

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//
  for(int iter=0;iter<Niter_lanczos;iter++){
   if(verbose==1) printf( " Lanczos iteration %d / %d \n", iter,Niter_lanczos);
   one_step_lanczos_cuda(nup,ndn,norb,blocksize,Niter_lanczos,size,iter,diag,subdiag,vec_tmp_gpu,vec_in_gpu,vec_out_gpu,
               vec_tmp_gpu_pointer,vec_in_gpu_pointer,vec_out_gpu_pointer,QUART_gpu_pointer,
               diagup_gpu_pointer,diagdn_gpu_pointer,noffup_gpu_pointer,noffdn_gpu_pointer,
               rankoffup_gpu_pointer,rankoffdn_gpu_pointer, offdiagup_gpu_pointer,offdiagdn_gpu_pointer,
               offdiaup_gpu_pointer,offdiadn_gpu_pointer,statesup_gpu_pointer,statesdn_gpu_pointer,UMASK_gpu_pointer
              ,iorbup_gpu_pointer,iorbdn_gpu_pointer,tmp_transp,tmp_transp_pointer);
  };
 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

     cudaEventSynchronize(0);
     cudaFreeHost(vec_tmp_gpu);

     if(use_texture_updo==0){
       cudaFreeHost(vec_in_gpu);
       if(USE_TRANSPOSE==1){ cudaFreeHost(tmp_transp);}
     }else{
       cudaFree(vec_in_gpu_pointer); 
       if(USE_TRANSPOSE==1){
         cudaFree(tmp_transp_pointer);
         unbind_x_();
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

extern "C" void lanczos_real_updo_dynamic_cuda_(int *norbs, int *pblocksize, 
   int *Niter_lanczos_,int *offdiasizeup, int* offdiasizedn, 
   int *roffup, int *roffdn, int *ntot,  int *sizup, int *sizdn, 
   double *QUART, double *diagup, double *diagdn, int *noffup, int *noffdn, int *rankoffup, int *rankoffdn, 
   double *offdiagup, double *offdiagdn, double *diag, double *subdiag, int *UMASK, int *statesup, int *statesdn , 
   int *iorbup, int *iorbdn, double *vecinit)
{

 //---------------------------------------------------------------------------------------//
     int verbose=0;
 //---------------------------------------------------------------------------------------//

     int Niter_lanczos= *Niter_lanczos_; int nup=*sizup; int ndn=*sizdn; int norb=*norbs;

     if(verbose==1) printf(" start Lanczos Real on GPU \n" );
     if(verbose==1) printf(" norbitals=%d \n ", norb );

     int size = *ntot; int blocksize= *pblocksize * WARP_SIZE ; 

     double   *vec_in_gpu,*vec_out_gpu,*QUART_gpu,*diagup_gpu,*diagdn_gpu,*offdiagup_gpu,*offdiagdn_gpu,*vec_tmp_gpu;
     int      *noffup_gpu,*noffdn_gpu,*rankoffup_gpu,*rankoffdn_gpu;
     double   *vec_in_gpu_pointer,*vec_out_gpu_pointer,*vec_tmp_gpu_pointer;
     double   *QUART_gpu_pointer,*diagup_gpu_pointer,*diagdn_gpu_pointer,*offdiagup_gpu_pointer,*offdiagdn_gpu_pointer;
     int      *noffup_gpu_pointer,*noffdn_gpu_pointer,*rankoffup_gpu_pointer,*rankoffdn_gpu_pointer;
     int      *offdiaup_gpu,*offdiadn_gpu, *offdiaup_gpu_pointer,*offdiadn_gpu_pointer;
     int      *statesup_gpu,*statesdn_gpu,*statesup_gpu_pointer,*statesdn_gpu_pointer;
     int      *UMASK_gpu,*UMASK_gpu_pointer;
     int      *iorbup_gpu,*iorbdn_gpu,*iorbup_gpu_pointer,*iorbdn_gpu_pointer;
     double   *tmp_transp,*tmp_transp_pointer;

     if(verbose==1) printf(" GPU get eigenvalues, size of Lanczos vectors = %d ", size);

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

       cudaSetDevice(0); cudaSetDeviceFlags( cudaDeviceMapHost ); cudaEventSynchronize(0);

       cudaHostAlloc((void**)&vec_tmp_gpu   , size*sizeof(double),  cudaHostAllocMapped | cudaHostAllocPortable );
 
       if(USE_TRANSPOSE==1){ 
        if(use_texture_updo==0){
          cudaHostAlloc((void**)&tmp_transp, size*sizeof(double),   cudaHostAllocMapped | cudaHostAllocPortable );
          cudaHostGetDevicePointer((void**)  &tmp_transp_pointer,tmp_transp,0);
         }else{
          cudaMalloc((void**)&tmp_transp_pointer,size*sizeof(double) );
          bind_x_(tmp_transp_pointer,size);
         };
       };
 
       if(use_texture_updo==1){
         cudaMalloc((void**)&vec_in_gpu_pointer,sizeof(double)*size);
         bind_x(vec_in_gpu_pointer,size);
       }else{
         cudaHostAlloc((void**)&vec_in_gpu , size*sizeof(double),   cudaHostAllocMapped | cudaHostAllocPortable );
         cudaHostGetDevicePointer((void**)  &vec_in_gpu_pointer  ,  vec_in_gpu    , 0 );
       }

       //-----------------------------------------------------------------------------------------------//
       //-----------------------------------------------------------------------------------------------//
        cudaHostAlloc((void**)&statesup_gpu    , nup *sizeof(int),           cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&statesdn_gpu    , ndn *sizeof(int),           cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&diagdn_gpu    , ndn*sizeof(double),           cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&offdiagdn_gpu , *offdiasizedn*sizeof(double), cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&noffdn_gpu    , ndn*sizeof(int)   ,           cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&rankoffdn_gpu , *roffdn*sizeof(int)  ,        cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&diagup_gpu    , nup*sizeof(double),           cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&offdiagup_gpu , *offdiasizeup*sizeof(double), cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&noffup_gpu    , nup*sizeof(int)   ,           cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&rankoffup_gpu , *roffup*sizeof(int)  ,        cudaHostAllocMapped | cudaHostAllocPortable );

        cudaMemcpy(statesup_gpu,  statesup,     nup * sizeof(int),           cudaMemcpyHostToHost);
        cudaMemcpy(statesdn_gpu,  statesdn,     ndn * sizeof(int),           cudaMemcpyHostToHost);
        cudaMemcpy(diagup_gpu,    diagup,       nup*sizeof(double),          cudaMemcpyHostToHost);
        cudaMemcpy(offdiagup_gpu, offdiagup,   *offdiasizeup*sizeof(double), cudaMemcpyHostToHost);
        cudaMemcpy(noffup_gpu,    noffup ,      nup*sizeof(int),             cudaMemcpyHostToHost);
        cudaMemcpy(rankoffup_gpu, rankoffup,   *roffup*sizeof(int),          cudaMemcpyHostToHost);
        cudaMemcpy(diagdn_gpu,    diagdn,       ndn*sizeof(double),          cudaMemcpyHostToHost);
        cudaMemcpy(offdiagdn_gpu, offdiagdn,   *offdiasizedn*sizeof(double), cudaMemcpyHostToHost);
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

       cudaHostAlloc((void**)&vec_out_gpu   , size*sizeof(double),        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaEventSynchronize(0);
       cudaHostGetDevicePointer((void**)  &vec_tmp_gpu_pointer ,  vec_tmp_gpu   , 0 );
       cudaHostGetDevicePointer((void**)  &vec_out_gpu_pointer ,  vec_out_gpu   , 0 );
       cudaEventSynchronize(0);

 //---------------------------------------------------------------------------------------//

  if(verbose==1) printf( " initialized, now run actual Lanczos \n " );

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

  cudaThreadSynchronize(); cudaEventSynchronize(0); 
  if(use_texture_updo==0){
   cudaMemcpy(vec_in_gpu,vecinit,size*sizeof(double),cudaMemcpyHostToHost);
  }else{
   cudaMemcpy(vec_out_gpu,vecinit,size*sizeof(double),cudaMemcpyHostToHost);
   cudaMemcpy(vec_in_gpu_pointer,vec_out_gpu,size*sizeof(double),cudaMemcpyHostToDevice);
  }
  cudaThreadSynchronize(); cudaEventSynchronize(0);

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//
  for(int iter=0;iter<Niter_lanczos;iter++){
   if(verbose==1) printf( " Lanczos iteration %d / %d \n", iter,Niter_lanczos);
   one_step_lanczos_cuda(nup,ndn,norb,blocksize,Niter_lanczos,size,iter,diag,subdiag,vec_tmp_gpu,vec_in_gpu,vec_out_gpu,
               vec_tmp_gpu_pointer,vec_in_gpu_pointer,vec_out_gpu_pointer,QUART_gpu_pointer,
               diagup_gpu_pointer,diagdn_gpu_pointer,noffup_gpu_pointer,noffdn_gpu_pointer,
               rankoffup_gpu_pointer,rankoffdn_gpu_pointer, offdiagup_gpu_pointer,offdiagdn_gpu_pointer,
               offdiaup_gpu_pointer,offdiadn_gpu_pointer,statesup_gpu_pointer,statesdn_gpu_pointer,UMASK_gpu_pointer
              ,iorbup_gpu_pointer,iorbdn_gpu_pointer,tmp_transp,tmp_transp_pointer);
  };
 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

     cudaEventSynchronize(0);
     cudaFreeHost(vec_tmp_gpu);

     if(use_texture_updo==0){
       cudaFreeHost(vec_in_gpu);
       if(USE_TRANSPOSE==1){ cudaFreeHost(tmp_transp);}
     }else{
       cudaFree(vec_in_gpu_pointer); 
       if(USE_TRANSPOSE==1){
         cudaFree(tmp_transp_pointer);
         unbind_x_();
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

extern "C" void lanczos_real_updo_gs_cuda_(int *norbs, int *pblocksize, 
   int *Niter_lanczos_,int *offdiasizeup, int* offdiasizedn, 
   int *roffup, int *roffdn, int *ntot,  int *sizup, int *sizdn, 
   double *QUART, double *diagup, double *diagdn, int *noffup, int *noffdn, int *rankoffup, int *rankoffdn, 
   double *offdiagup, double *offdiagdn,  int *UMASK, int *statesup, int *statesdn , 
   int *iorbup, int *iorbdn, double *vecp, double *GS)
{

     double diag[*Niter_lanczos_],subdiag[*Niter_lanczos_];

 //---------------------------------------------------------------------------------------//
     int verbose=0;
 //---------------------------------------------------------------------------------------//

     int Niter_lanczos= *Niter_lanczos_; int nup=*sizup; int ndn=*sizdn; int norb=*norbs;

     if(verbose==1) printf(" start Lanczos Real on GPU \n" );

     int size = *ntot; int blocksize= *pblocksize * WARP_SIZE ; 

     double   *vec_in_gpu,*vec_out_gpu,*QUART_gpu,*diagup_gpu,*diagdn_gpu,*offdiagup_gpu,*offdiagdn_gpu,*vec_tmp_gpu;
     int      *noffup_gpu,*noffdn_gpu,*rankoffup_gpu,*rankoffdn_gpu;
     double   *vec_in_gpu_pointer,*vec_out_gpu_pointer,*vec_tmp_gpu_pointer;
     double   *QUART_gpu_pointer,*diagup_gpu_pointer,*diagdn_gpu_pointer,*offdiagup_gpu_pointer,*offdiagdn_gpu_pointer;
     int      *noffup_gpu_pointer,*noffdn_gpu_pointer,*rankoffup_gpu_pointer,*rankoffdn_gpu_pointer;
     int      *offdiaup_gpu,*offdiadn_gpu, *offdiaup_gpu_pointer,*offdiadn_gpu_pointer;
     int      *statesup_gpu,*statesdn_gpu,*statesup_gpu_pointer,*statesdn_gpu_pointer;
     int      *UMASK_gpu,*UMASK_gpu_pointer;
     int      *iorbup_gpu,*iorbdn_gpu,*iorbup_gpu_pointer,*iorbdn_gpu_pointer;
     double   *tmp_transp,*tmp_transp_pointer;

     double *GS_gpu, *GS_gpu_pointer;

     if(verbose==1) printf(" GPU get Ground State, size of Lanczos vectors = %d ", size);

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

       cudaSetDevice(0); cudaSetDeviceFlags( cudaDeviceMapHost ); cudaEventSynchronize(0);

       cudaHostAlloc((void**)&vec_tmp_gpu   , size*sizeof(double),  cudaHostAllocMapped | cudaHostAllocPortable );

       cudaHostAlloc((void**)&GS_gpu        , size*sizeof(double),        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostGetDevicePointer((void**)   &GS_gpu_pointer ,       GS_gpu         , 0 );
 
       if(USE_TRANSPOSE==1){ 
        if(use_texture_updo==0){
          cudaHostAlloc((void**)&tmp_transp, size*sizeof(double),   cudaHostAllocMapped | cudaHostAllocPortable );
          cudaHostGetDevicePointer((void**)  &tmp_transp_pointer,tmp_transp,0);
         }else{
          cudaMalloc((void**)&tmp_transp_pointer,size*sizeof(double) );
          bind_x_(tmp_transp_pointer,size);
         };
       };
 
       if(use_texture_updo==1){
         cudaMalloc((void**)&vec_in_gpu_pointer,sizeof(double)*size);
         bind_x(vec_in_gpu_pointer,size);
       }else{
         cudaHostAlloc((void**)&vec_in_gpu , size*sizeof(double),   cudaHostAllocMapped | cudaHostAllocPortable );
         cudaHostGetDevicePointer((void**)  &vec_in_gpu_pointer  ,  vec_in_gpu    , 0 );
       }

       //-----------------------------------------------------------------------------------------------//
       //-----------------------------------------------------------------------------------------------//
        cudaHostAlloc((void**)&statesup_gpu    , nup *sizeof(int),           cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&statesdn_gpu    , ndn *sizeof(int),           cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&diagdn_gpu    , ndn*sizeof(double),           cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&offdiagdn_gpu , *offdiasizedn*sizeof(double), cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&noffdn_gpu    , ndn*sizeof(int)   ,           cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&rankoffdn_gpu , *roffdn*sizeof(int)  ,        cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&diagup_gpu    , nup*sizeof(double),           cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&offdiagup_gpu , *offdiasizeup*sizeof(double), cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&noffup_gpu    , nup*sizeof(int)   ,           cudaHostAllocMapped | cudaHostAllocPortable );
        cudaHostAlloc((void**)&rankoffup_gpu , *roffup*sizeof(int)  ,        cudaHostAllocMapped | cudaHostAllocPortable );

        cudaMemcpy(statesup_gpu,  statesup,     nup * sizeof(int),           cudaMemcpyHostToHost);
        cudaMemcpy(statesdn_gpu,  statesdn,     ndn * sizeof(int),           cudaMemcpyHostToHost);
        cudaMemcpy(diagup_gpu,    diagup,       nup*sizeof(double),          cudaMemcpyHostToHost);
        cudaMemcpy(offdiagup_gpu, offdiagup,   *offdiasizeup*sizeof(double), cudaMemcpyHostToHost);
        cudaMemcpy(noffup_gpu,    noffup ,      nup*sizeof(int),             cudaMemcpyHostToHost);
        cudaMemcpy(rankoffup_gpu, rankoffup,   *roffup*sizeof(int),          cudaMemcpyHostToHost);
        cudaMemcpy(diagdn_gpu,    diagdn,       ndn*sizeof(double),          cudaMemcpyHostToHost);
        cudaMemcpy(offdiagdn_gpu, offdiagdn,   *offdiasizedn*sizeof(double), cudaMemcpyHostToHost);
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

       cudaHostAlloc((void**)&vec_out_gpu   , size*sizeof(double),        cudaHostAllocMapped | cudaHostAllocPortable );
       cudaEventSynchronize(0);
       cudaHostGetDevicePointer((void**)  &vec_tmp_gpu_pointer ,  vec_tmp_gpu   , 0 );
       cudaHostGetDevicePointer((void**)  &vec_out_gpu_pointer ,  vec_out_gpu   , 0 );
       cudaEventSynchronize(0);

 //---------------------------------------------------------------------------------------//

  if(verbose==1) printf( " initialized, now run actual Lanczos \n " );

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

  cudaThreadSynchronize(); cudaEventSynchronize(0); 
 if(use_texture_updo==0){
  for(int i=0;i<size;i++) vec_in_gpu[i]=1.0; cudaThreadSynchronize();
 }else{
   for(int i=0;i<size;i++) vec_out_gpu[i]=1.0; cudaThreadSynchronize();
   cudaMemcpy(vec_in_gpu_pointer,vec_out_gpu,size*sizeof(double),cudaMemcpyHostToDevice);
 }
   cudaThreadSynchronize(); cudaEventSynchronize(0);

  memset ((void **)GS_gpu, 0, size*sizeof(double));

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

  double coef= 1.0/sqrt(cublasDdot(size,vec_in_gpu_pointer,1,vec_in_gpu_pointer,1))*vecp[0]; cudaThreadSynchronize();  cudaEventSynchronize(0);
  cublasDaxpy(size,coef,vec_in_gpu_pointer,1,GS_gpu_pointer,1); cudaThreadSynchronize(); cudaEventSynchronize(0);cudaThreadSynchronize();

   one_step_lanczos_cuda(nup,ndn,norb,blocksize,Niter_lanczos,size,0,diag,subdiag,vec_tmp_gpu,vec_in_gpu,vec_out_gpu,
               vec_tmp_gpu_pointer,vec_in_gpu_pointer,vec_out_gpu_pointer,QUART_gpu_pointer,
               diagup_gpu_pointer,diagdn_gpu_pointer,noffup_gpu_pointer,noffdn_gpu_pointer,
               rankoffup_gpu_pointer,rankoffdn_gpu_pointer, offdiagup_gpu_pointer,offdiagdn_gpu_pointer,
               offdiaup_gpu_pointer,offdiadn_gpu_pointer,statesup_gpu_pointer,statesdn_gpu_pointer,UMASK_gpu_pointer
              ,iorbup_gpu_pointer,iorbdn_gpu_pointer,tmp_transp,tmp_transp_pointer);

  for(int iter=1;iter<Niter_lanczos-1;iter++){

   coef=1.0/sqrt(cublasDdot(size,vec_out_gpu_pointer,1,vec_out_gpu_pointer,1)); cublasDscal (size,coef,vec_out_gpu_pointer,1);  cudaEventSynchronize(0);
   cublasDaxpy(size,vecp[iter],vec_out_gpu_pointer,1,GS_gpu_pointer,1); cudaThreadSynchronize(); cudaEventSynchronize(0);

   one_step_lanczos_cuda(nup,ndn,norb,blocksize,Niter_lanczos,size,iter,diag,subdiag,vec_tmp_gpu,vec_in_gpu,vec_out_gpu,
               vec_tmp_gpu_pointer,vec_in_gpu_pointer,vec_out_gpu_pointer,QUART_gpu_pointer,
               diagup_gpu_pointer,diagdn_gpu_pointer,noffup_gpu_pointer,noffdn_gpu_pointer,
               rankoffup_gpu_pointer,rankoffdn_gpu_pointer, offdiagup_gpu_pointer,offdiagdn_gpu_pointer,
               offdiaup_gpu_pointer,offdiadn_gpu_pointer,statesup_gpu_pointer,statesdn_gpu_pointer,UMASK_gpu_pointer
              ,iorbup_gpu_pointer,iorbdn_gpu_pointer,tmp_transp,tmp_transp_pointer);

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

     if(use_texture_updo==0){
       cudaFreeHost(vec_in_gpu);
       if(USE_TRANSPOSE==1){ cudaFreeHost(tmp_transp);}
     }else{
       cudaFree(vec_in_gpu_pointer); 
       if(USE_TRANSPOSE==1){
         cudaFree(tmp_transp_pointer);
         unbind_x_();
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

