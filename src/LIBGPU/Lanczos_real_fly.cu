#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <cublas.h>

#define IBSET(a,b) ((a) |= (1<<(b)))
#define IBCLR(a,b) ((a) &= ~(1<<(b)))
#define BTEST(a,b) ((a) & (1<<(b)))

const int use_texture_lrf = 1 ;
      int use_texture_rank = 1 ;

#define WARP_SIZE 32
#define MAX_BLOCK 65500 
#define max_rank 225000000

texture<int2,1,cudaReadModeElementType> tex,texb,texc,texbc;
texture<int,1,cudaReadModeElementType>  texr;

inline void   bind_x_r(int *x)         {cudaBindTexture(0,texr,x);};
inline void unbind_x_r()               {cudaUnbindTexture(texr);};

inline void   bind_x(double *x, int N) {cudaBindTexture(0,tex,x,N*sizeof(double));};
inline void unbind_x()                 {cudaUnbindTexture(tex);};
inline void   bind_x_(double *x, double *y, double *z ) 
                                       {cudaBindTexture(0,texb,x); cudaBindTexture(0,texc,y); cudaBindTexture(0,texbc,z);};
inline void unbind_x_()                {cudaUnbindTexture(texb);cudaUnbindTexture(texc);cudaUnbindTexture(texbc);};

__inline__  __device__ int    fetch_r(const int& i) {  int  v = tex1Dfetch(texr,i);  return v;                          }
__inline__  __device__ double fetch_x(const int& i) {  int2 v = tex1Dfetch(tex,i);   return __hiloint2double(v.y, v.x); }
__inline__  __device__ double fetchb(const int& i)  {  int2 v = tex1Dfetch(texb,i);  return __hiloint2double(v.y, v.x); }
__inline__  __device__ double fetchc(const int& i)  {  int2 v = tex1Dfetch(texc,i);  return __hiloint2double(v.y, v.x); }
__inline__  __device__ double fetchbc(const int& i) {  int2 v = tex1Dfetch(texbc,i); return __hiloint2double(v.y, v.x); }

//********************************************
//********************************************
//********************************************
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

  __global__  void Hmult_ker_fly(int ngrid,int blocksize, int psize, double *vec_out, const double *vec_in, const double *quart, const double* Eb, 
    const double* Ec, const double* Vbc, const int* sector_states, const int* sector_rank, const int norbs, const int bathnorbs, const int impnorbs, 
    const int* imporbs_, const int* bathorbs_, const int* maskEb_, const int* maskEc_, const int* maskVbc_, const int use_texture_lrf,const int use_texture_rank)

{
   __shared__ double sdata[16][WARP_SIZE];
   __shared__ short int maskEc[32*32],maskEb[32*32],maskVbc[32*32],imporbs[32],bathorbs[32];

          int istate ;  
    const int warp_lane   = threadIdx.y; 
    const int thread_lane = threadIdx.x; 

    int        jstate,kets_out,kets_in,jj,iorb,jorb,n1,n2; 
    double     hoffdiag,vecin;
    short int  fermion_sign;

    if(warp_lane==0) for(iorb=thread_lane;iorb<impnorbs*impnorbs;iorb+=WARP_SIZE)   maskEc[iorb]=maskEc_[iorb];
    if(warp_lane==1) for(iorb=thread_lane;iorb<bathnorbs*bathnorbs;iorb+=WARP_SIZE) maskEb[iorb]=maskEb_[iorb];
    if(warp_lane==2) for(iorb=thread_lane;iorb<bathnorbs*impnorbs;iorb+=WARP_SIZE)  maskVbc[iorb]=maskVbc_[iorb];
    if(warp_lane==3) for(iorb=thread_lane;iorb<bathnorbs;iorb+=WARP_SIZE)           bathorbs[iorb]=bathorbs_[iorb];
    if(warp_lane==4) for(iorb=thread_lane;iorb<impnorbs;iorb+=WARP_SIZE)            imporbs[iorb]=imporbs_[iorb];
    __syncthreads();

   for(int iii=0;iii<=ngrid;iii++){
    istate=blocksize * ( blockIdx.y + iii*MAX_BLOCK )+ threadIdx.y ;
    if(istate < psize)
  {

    if(use_texture_lrf==0){ vecin =  vec_in[istate];
                     }else{ vecin = fetch_x(istate); };

    kets_in = sector_states[istate];
    if(thread_lane==0){
      sdata[threadIdx.y][threadIdx.x] = quart[istate] * vecin;
    }else{
      sdata[threadIdx.y][threadIdx.x] = 0.0;
    };

    //////////////////////////////////////////////////////////////////////////

   if(use_texture_lrf==0){
    for(iorb=thread_lane;iorb<impnorbs ;iorb+=WARP_SIZE)
      {if(BTEST(kets_in,imporbs[iorb]-1 )>0) sdata[threadIdx.y][threadIdx.x] += Ec[iorb* impnorbs+iorb] * vecin;};
    for(iorb=thread_lane;iorb<bathnorbs;iorb+=WARP_SIZE)
      {if(BTEST(kets_in,bathorbs[iorb]-1)>0) sdata[threadIdx.y][threadIdx.x] += Eb[iorb*bathnorbs+iorb] * vecin;};
   }else{
    for(iorb=thread_lane;iorb<impnorbs ;iorb+=WARP_SIZE)
      {if(BTEST(kets_in,imporbs[iorb]-1 )>0) sdata[threadIdx.y][threadIdx.x] += fetchc(iorb* impnorbs+iorb) * vecin;};
    for(iorb=thread_lane;iorb<bathnorbs;iorb+=WARP_SIZE)
      {if(BTEST(kets_in,bathorbs[iorb]-1)>0) sdata[threadIdx.y][threadIdx.x] += fetchb(iorb*bathnorbs+iorb) * vecin;};
   }

    //////////////////////////////////////////////////////////////////////////
     for(jorb=thread_lane;jorb<impnorbs;jorb+=WARP_SIZE){
      if(BTEST(kets_in,imporbs[jorb]-1)>0){

       for(iorb=0;iorb<impnorbs;iorb++){
       if(iorb!=jorb){
       if(maskEc[jorb*impnorbs+iorb]!=0){

        kets_out=kets_in; IBCLR(kets_out,imporbs[jorb]-1);
        if(BTEST(kets_out,imporbs[iorb]-1)==0)
        {
         IBSET(kets_out,imporbs[iorb]-1);
         n1=imporbs[iorb]-1;n2=imporbs[jorb]-1; fermion_sign=1;
         if(n1<n2){ for(jj=n1+1; jj<=n2-1; jj++) { if(BTEST(kets_out,jj)>0) fermion_sign=-fermion_sign; };};
         if(n2<n1){ for(jj=n2+1; jj<=n1-1; jj++) { if(BTEST(kets_out,jj)>0) fermion_sign=-fermion_sign; };};
           if(use_texture_rank==0){
            jstate                            = sector_rank[kets_out]-1;}else{
            jstate                            = fetch_r(kets_out)-1;};
           if(use_texture_lrf==0){
            hoffdiag                          = Ec[jorb*impnorbs+iorb] * fermion_sign;
            sdata[threadIdx.y][threadIdx.x]  += hoffdiag*vec_in[jstate]; 
           }else{
            sdata[threadIdx.y][threadIdx.x]  += fetch_x(jstate)*fetchc(jorb*impnorbs+iorb)*fermion_sign;
           }
        }; }; }; }; }; };
    //////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////
     for(jorb=thread_lane;jorb<bathnorbs;jorb+=WARP_SIZE){
      if(BTEST(kets_in,bathorbs[jorb]-1)>0){

       for(iorb=0;iorb<bathnorbs;iorb++){
       if(iorb!=jorb){
       if(maskEb[jorb*bathnorbs+iorb]!=0){

        kets_out=kets_in; IBCLR(kets_out,bathorbs[jorb]-1);
        if(BTEST(kets_out,bathorbs[iorb]-1)==0)
        {
         IBSET(kets_out,bathorbs[iorb]-1);
         n1=bathorbs[iorb]-1;n2=bathorbs[jorb]-1; fermion_sign=1;
         if(n1<n2){ for(int jj=n1+1; jj<=n2-1; jj++) { if(BTEST(kets_out,jj)>0) fermion_sign=-fermion_sign; };};
         if(n2<n1){ for(int jj=n2+1; jj<=n1-1; jj++) { if(BTEST(kets_out,jj)>0) fermion_sign=-fermion_sign; };};
         if(use_texture_rank==0){
            jstate                            = sector_rank[kets_out]-1;}else{
            jstate                            = fetch_r(kets_out)-1;};
          if(use_texture_lrf==0){
           hoffdiag                         = Eb[jorb*bathnorbs+iorb]*fermion_sign;
           sdata[threadIdx.y][threadIdx.x] += hoffdiag*vec_in[jstate];
          }else{
           sdata[threadIdx.y][threadIdx.x] += fetch_x(jstate)*fetchb(jorb*bathnorbs+iorb)*fermion_sign;
          }
        }; }; }; }; }; };
    //////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////
      for(jorb=thread_lane;jorb<impnorbs;jorb+=WARP_SIZE){ 
      if(BTEST(kets_in,imporbs[jorb]-1)>0){

      for(iorb=0;iorb<bathnorbs;iorb++){if(maskVbc[jorb*bathnorbs+iorb]!=0){
        kets_out=kets_in; IBCLR(kets_out,imporbs[jorb]-1);
        if(BTEST(kets_out,bathorbs[iorb]-1)==0)
        {
         IBSET(kets_out,bathorbs[iorb]-1);
         n1=imporbs[jorb]-1;n2=bathorbs[iorb]-1; fermion_sign=1;
         if(n1<n2){ for(int jj=n1+1; jj<=n2-1; jj++) { if(BTEST(kets_out,jj)>0) fermion_sign=-fermion_sign; };};
         if(n2<n1){ for(int jj=n2+1; jj<=n1-1; jj++) { if(BTEST(kets_out,jj)>0) fermion_sign=-fermion_sign; };};
          if(use_texture_rank==0){
            jstate                            = sector_rank[kets_out]-1;}else{
            jstate                            = fetch_r(kets_out)-1;};
           if(use_texture_lrf==0){
             hoffdiag                         =  Vbc[jorb*bathnorbs+iorb]*fermion_sign;
             sdata[threadIdx.y][threadIdx.x]  += hoffdiag*vec_in[jstate];
           }else{
             sdata[threadIdx.y][threadIdx.x]  += fetch_x(jstate)*fetchbc(jorb*bathnorbs+iorb)*fermion_sign;
           }
        }; }; }; }; };
    //////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////
    for(iorb=thread_lane;iorb<bathnorbs;iorb+=WARP_SIZE){ 
      if(BTEST(kets_in,bathorbs[iorb]-1)>0){

        for(jorb=0;jorb<impnorbs;jorb++){ if(maskVbc[jorb*bathnorbs+iorb]!=0){
        kets_out=kets_in; IBCLR(kets_out,bathorbs[iorb]-1);
        if(BTEST(kets_out,imporbs[jorb]-1)==0)
        {
         IBSET(kets_out,imporbs[jorb]-1);
         n1=imporbs[jorb]-1;n2=bathorbs[iorb]-1; fermion_sign=1;
         if(n1<n2){ for(int jj=n1+1; jj<=n2-1; jj++) { if(BTEST(kets_out,jj)>0) fermion_sign=-fermion_sign; };};
         if(n2<n1){ for(int jj=n2+1; jj<=n1-1; jj++) { if(BTEST(kets_out,jj)>0) fermion_sign=-fermion_sign; };};
          if(use_texture_rank==0){
            jstate                            = sector_rank[kets_out]-1;}else{
            jstate                            = fetch_r(kets_out)-1;};
           if(use_texture_lrf==0){
             hoffdiag                         =  Vbc[jorb*bathnorbs+iorb]*fermion_sign;
             sdata[threadIdx.y][threadIdx.x]  += hoffdiag*vec_in[jstate];
           }else{
             sdata[threadIdx.y][threadIdx.x]  += fetch_x(jstate)*fetchbc(jorb*bathnorbs+iorb)*fermion_sign;
           }
        }; }; }; }; };
    //////////////////////////////////////////////////////////////////////////

        if (thread_lane < 16) { sdata[threadIdx.y][threadIdx.x] += sdata[threadIdx.y][threadIdx.x + 16]; };
        if (thread_lane <  8) { sdata[threadIdx.y][threadIdx.x] += sdata[threadIdx.y][threadIdx.x +  8]; };
        if (thread_lane <  4) { sdata[threadIdx.y][threadIdx.x] += sdata[threadIdx.y][threadIdx.x +  4]; };
        if (thread_lane <  2) { sdata[threadIdx.y][threadIdx.x] += sdata[threadIdx.y][threadIdx.x +  2]; };
        if (thread_lane <  1) { sdata[threadIdx.y][threadIdx.x] += sdata[threadIdx.y][threadIdx.x +  1]; };
        if (thread_lane == 0)   vec_out[istate] = sdata[threadIdx.y][threadIdx.x];

    //////////////////////////////////////////////////////////////////////////

  };
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
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************

 void one_step_lanczos_fly_cuda(int blocksize, int Niter, int size, int iter,  double *vec_tmp_gpu, double *vec_in_gpu, double *vec_out_gpu,
     double *vec_tmp_gpu_pointer, double *vec_in_gpu_pointer, double *vec_out_gpu_pointer, double *quart_gpu_pointer , double *diag, double *subdiag
   , double* Eb, double* Ec, double* Vbc, int* sector_states, int* sector_ranks,int norbs, int bathnorbs,int impnorbs,int* imporbs, int* bathorbs, 
     int* maskEb, int* maskEc, int* maskVbc )

{

    int verbose=0; int psize=size; if(verbose==1) printf( " define blocksize \n ");
 
    int nb=(size - size % blocksize) / blocksize + 1 ; int ngrid=nb/MAX_BLOCK; if(ngrid>0) nb=MAX_BLOCK;
    if(verbose==1) printf( " --------------- \n  Nblock=%d Ngrid=%d \n ----------------- \n ",nb,ngrid);
    dim3 bl(1,nb),th(WARP_SIZE,blocksize);

    if(verbose==1) printf( " Sdot, calculate vecin norm \n ");

    double normv = sqrt(cublasDdot(size,vec_in_gpu_pointer,1,vec_in_gpu_pointer,1)); cudaEventSynchronize(0); cudaThreadSynchronize();

    if(verbose==1) printf( " done norm=%f \n ", normv);
 
    cublasDscal(size,1.0/normv,vec_in_gpu_pointer,1); cudaEventSynchronize(0); cudaThreadSynchronize(); 

    cudaEventSynchronize(0); cudaThreadSynchronize();
    Hmult_ker_fly<<<bl,th>>>(ngrid,blocksize,psize,vec_out_gpu_pointer,vec_in_gpu_pointer,quart_gpu_pointer,Eb,Ec,Vbc,sector_states,sector_ranks,norbs,
                             bathnorbs,impnorbs,imporbs,bathorbs,maskEb,maskEc,maskVbc,use_texture_lrf,use_texture_rank); 
    cudaEventSynchronize(0); cudaThreadSynchronize();

    if(iter>0){cublasDaxpy(size,-subdiag[iter],vec_tmp_gpu_pointer,1,vec_out_gpu_pointer,1);}; cudaEventSynchronize(0); cudaThreadSynchronize();

   if(verbose==1) printf(" iteration carried on, now copy vecin to vectmp \n ");
 
   if(use_texture_lrf==0){
    cudaMemcpy(vec_tmp_gpu,vec_in_gpu,size*sizeof(double),cudaMemcpyHostToHost);cudaEventSynchronize(0);cudaThreadSynchronize();
   }else{
    cudaMemcpy(vec_tmp_gpu,vec_in_gpu_pointer,size*sizeof(double),cudaMemcpyDeviceToHost);cudaEventSynchronize(0);cudaThreadSynchronize();
   }

   if(verbose==1) printf(" now scalar product .... \n " );
    
    diag[iter]=cublasDdot(size, vec_out_gpu_pointer,1,vec_in_gpu_pointer,1);cudaEventSynchronize(0);cudaThreadSynchronize();

    if(verbose==1) printf(" now calculate scalprod \n " );
    cublasDaxpy(size,-diag[iter],vec_tmp_gpu_pointer,1,vec_out_gpu_pointer,1); cudaEventSynchronize(0); cudaThreadSynchronize();
    normv = sqrt(cublasDdot(size, vec_out_gpu_pointer,1,vec_out_gpu_pointer,1)); cudaEventSynchronize(0); cudaThreadSynchronize();

    if(iter<Niter-1) subdiag[iter+1]=normv;

  if(verbose==1) printf( " subdiag(%d)= %f \n " , iter+1,subdiag[iter+1]);
  if(verbose==1) printf( " now copy vecout to vecin \n " );

  if(use_texture_lrf==0){
    cudaMemcpy( vec_in_gpu, vec_out_gpu, size*sizeof(double), cudaMemcpyHostToHost); cudaEventSynchronize(0); cudaThreadSynchronize();
  }else{
    cudaMemcpy( vec_in_gpu_pointer, vec_out_gpu, size*sizeof(double), cudaMemcpyHostToDevice); cudaEventSynchronize(0); cudaThreadSynchronize();
  }
  if(verbose==1) printf( "done now return \n ");

  cudaEventSynchronize(0); cudaThreadSynchronize();

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

extern "C" void lanczos_real_fly_cuda_(int *dimen_in, int *pblocksize, int* norbs_in, int *Niter_lanczos_, int* ntot, double* quart, 
      double* diag, double* subdiag, double* Eb_in, double* Ec_in, double* Vbc_in, int* sector_states_in, int* sector_ranks_in 
     ,int* bathnorbs_in,int* impnorbs_in,int* imporbs_in, int* bathorbs_in, int* maskEb_in, int* maskEc_in, int* maskVbc_in)
{

  if(*dimen_in<max_rank){ use_texture_rank = 1;}else{use_texture_rank = 0;};

 //---------------------------------------------------------------------------------------//
     int verbose=0;
 //---------------------------------------------------------------------------------------//

     int Niter_lanczos=*Niter_lanczos_; int norbs=*norbs_in; int dimen=*dimen_in;
     
     if(verbose==1) printf(" start Lanczos Real on GPU \n" );

     int size      = *ntot; int blocksize = *pblocksize;
     int bathnorbs = *bathnorbs_in; int impnorbs=*impnorbs_in;

     double  *vec_in_gpu,*vec_out_gpu,*vec_tmp_gpu,*quart_gpu;
     double  *vec_in_gpu_pointer,*vec_out_gpu_pointer,*vec_tmp_gpu_pointer,*quart_gpu_pointer;
     double  *Eb,*Ec,*Vbc,*Eb_p,*Ec_p,*Vbc_p;

     int *sector_states_p,*sector_ranks_p,*imporbs_p,*bathorbs_p,*maskEb_p,*maskEc_p,*maskVbc_p;
     int *sector_states,*sector_ranks,*imporbs,*bathorbs,*maskEb,*maskEc,*maskVbc;

     if(verbose==1) printf(" GPU get eigenvalues, size of Lanczos vectors = %d ", size);

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

       cudaSetDevice(0); cudaSetDeviceFlags( cudaDeviceMapHost ); cudaEventSynchronize(0);
       if(verbose==1) printf( " using texture = %d \n " , use_texture_lrf);

       if(use_texture_lrf==1){
         cudaMalloc((void**)&vec_in_gpu_pointer,sizeof(double)*size);
         cudaMalloc((void**)&Eb_p  , bathnorbs*bathnorbs*sizeof(double) );
         cudaMalloc((void**)&Ec_p  ,  impnorbs*impnorbs*sizeof(double) );
         cudaMalloc((void**)&Vbc_p , bathnorbs*impnorbs*sizeof(double) );
         cudaMemcpy(Eb_p  ,  Eb_in,  bathnorbs*bathnorbs*sizeof(double), cudaMemcpyHostToDevice);
         cudaMemcpy(Ec_p  ,  Ec_in,   impnorbs*impnorbs*sizeof(double),  cudaMemcpyHostToDevice);
         cudaMemcpy(Vbc_p , Vbc_in,  bathnorbs*impnorbs*sizeof(double),  cudaMemcpyHostToDevice);
         bind_x(vec_in_gpu_pointer,size);
         bind_x_(Eb_p,Ec_p,Vbc_p);
       }else{
         cudaHostAlloc((void**)&vec_in_gpu , size*sizeof(double)  , cudaHostAllocMapped | cudaHostAllocPortable );
         cudaHostGetDevicePointer((void**)   &vec_in_gpu_pointer  , vec_in_gpu , 0 );
         cudaHostAlloc((void**)&Eb         , bathnorbs*bathnorbs*sizeof(double), cudaHostAllocMapped | cudaHostAllocPortable );
         cudaHostAlloc((void**)&Ec         , impnorbs*impnorbs*sizeof(double)  , cudaHostAllocMapped | cudaHostAllocPortable );
         cudaHostAlloc((void**)&Vbc        , bathnorbs*impnorbs*sizeof(double) , cudaHostAllocMapped | cudaHostAllocPortable );
         cudaMemcpy(Eb,  Eb_in,  bathnorbs*bathnorbs*sizeof(double), cudaMemcpyHostToHost);
         cudaMemcpy(Ec,  Ec_in,  impnorbs*impnorbs*sizeof(double),   cudaMemcpyHostToHost);
         cudaMemcpy(Vbc, Vbc_in, bathnorbs*impnorbs*sizeof(double),  cudaMemcpyHostToHost);
         cudaHostGetDevicePointer((void**)   &Eb_p  , Eb  , 0 );
         cudaHostGetDevicePointer((void**)   &Ec_p  , Ec  , 0 );
         cudaHostGetDevicePointer((void**)   &Vbc_p , Vbc , 0 );
       }
      cudaEventSynchronize(0); cudaThreadSynchronize();

      if(verbose==1) printf( " bind sector_ranks to texture if necessary = %d \n ", use_texture_rank);
      if(verbose==1) printf( " dimension of Hilbert space= %d \n ", dimen);

      if(use_texture_rank==0){
       cudaHostAlloc((void**)&sector_ranks  , dimen*sizeof(int), cudaHostAllocMapped | cudaHostAllocPortable );
       cudaMemcpy(sector_ranks , sector_ranks_in , dimen*sizeof(int), cudaMemcpyHostToHost);
       cudaHostGetDevicePointer((void**)   &sector_ranks_p      , sector_ranks   , 0 );
      }else{
       cudaMalloc((void**)&sector_ranks_p  , dimen*sizeof(int) );
       cudaMemcpy(         sector_ranks_p , sector_ranks_in , dimen*sizeof(int), cudaMemcpyHostToDevice);
       bind_x_r(sector_ranks_p);
      };
       cudaEventSynchronize(0); cudaThreadSynchronize();

       cudaHostAlloc((void**)&vec_tmp_gpu   , size*sizeof(double),                cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&vec_out_gpu   , size*sizeof(double),                cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&quart_gpu     , size*sizeof(double),                cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&sector_states , size*sizeof(int),                   cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&imporbs       , impnorbs*sizeof(int),               cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&bathorbs      , bathnorbs*sizeof(int),              cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&maskEb        , bathnorbs*bathnorbs*sizeof(int),    cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&maskEc        , impnorbs*impnorbs*sizeof(int),      cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&maskVbc       , bathnorbs*impnorbs*sizeof(int),     cudaHostAllocMapped | cudaHostAllocPortable );

       cudaEventSynchronize(0);
       cudaMemcpy(quart_gpu, quart , size*sizeof(double), cudaMemcpyHostToHost);
       cudaMemcpy(sector_states, sector_states_in, size*sizeof(int), cudaMemcpyHostToHost);
       cudaMemcpy(imporbs, imporbs_in , impnorbs*sizeof(int),  cudaMemcpyHostToHost);
       cudaMemcpy(bathorbs,bathorbs_in, bathnorbs*sizeof(int), cudaMemcpyHostToHost);
       cudaMemcpy(maskEb  ,maskEb_in , bathnorbs *bathnorbs  *sizeof(int), cudaMemcpyHostToHost);
       cudaMemcpy(maskEc  ,maskEc_in , impnorbs  *impnorbs   *sizeof(int), cudaMemcpyHostToHost);
       cudaMemcpy(maskVbc ,maskVbc_in, impnorbs  *bathnorbs  *sizeof(int), cudaMemcpyHostToHost);

       cudaEventSynchronize(0);
       cudaHostGetDevicePointer((void**)   &vec_tmp_gpu_pointer , vec_tmp_gpu    , 0 );
       cudaHostGetDevicePointer((void**)   &vec_out_gpu_pointer , vec_out_gpu    , 0 );
       cudaHostGetDevicePointer((void**)   &quart_gpu_pointer   , quart_gpu      , 0 );
       cudaHostGetDevicePointer((void**)   &sector_states_p     , sector_states  , 0 );
       cudaHostGetDevicePointer((void**)   &imporbs_p           , imporbs        , 0 );
       cudaHostGetDevicePointer((void**)   &bathorbs_p          , bathorbs       , 0 );
       cudaHostGetDevicePointer((void**)   &maskEb_p            , maskEb         , 0 );
       cudaHostGetDevicePointer((void**)   &maskEc_p            , maskEc         , 0 );
       cudaHostGetDevicePointer((void**)   &maskVbc_p           , maskVbc        , 0 );
       cudaEventSynchronize(0);

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

  if(verbose==1) printf( " initialized, now run actual Lanczos \n " );

 if(use_texture_lrf==0){
  for(int i=0;i<size;i++) vec_in_gpu[i]=1.0; cudaThreadSynchronize();
 }else{
   for(int i=0;i<size;i++) vec_out_gpu[i]=1.0; cudaThreadSynchronize();
   cudaMemcpy(vec_in_gpu_pointer,vec_out_gpu,size*sizeof(double),cudaMemcpyHostToDevice);
 }

  cudaEventSynchronize(0); cudaThreadSynchronize();
  for(int iter=0;iter<Niter_lanczos;iter++){

  if(verbose==1) printf( " Lanczos iteration %d / %d \n", iter,Niter_lanczos);

  one_step_lanczos_fly_cuda(blocksize,Niter_lanczos,size,iter,vec_tmp_gpu,vec_in_gpu, vec_out_gpu,vec_tmp_gpu_pointer,
                            vec_in_gpu_pointer,vec_out_gpu_pointer,quart_gpu_pointer,diag,subdiag,Eb_p,Ec_p,Vbc_p,sector_states_p,sector_ranks_p,norbs,
                            bathnorbs,impnorbs,imporbs_p,bathorbs_p,maskEb_p,maskEc_p,maskVbc_p);
  };

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

   cudaEventSynchronize(0); 
   if(use_texture_lrf==0)
   {
    cudaFreeHost(vec_in_gpu);cudaFreeHost(Eb);cudaFreeHost(Ec);cudaFreeHost(Vbc);
   }else{
    unbind_x(); unbind_x_(); cudaFree(vec_in_gpu_pointer); cudaFree(Eb_p);cudaFree(Ec_p);cudaFree(Vbc_p);
   }
   cudaEventSynchronize(0); cudaThreadSynchronize();

   if(use_texture_rank==0){ cudaFreeHost(sector_ranks);}else{unbind_x_r(); cudaFree(sector_ranks_p);};

   cudaFreeHost(vec_tmp_gpu);    cudaFreeHost(vec_out_gpu); 
   cudaFreeHost(sector_states);  cudaFreeHost(quart_gpu);
   cudaFreeHost(bathorbs);       cudaFreeHost(imporbs);        cudaFreeHost(maskVbc);
   cudaFreeHost(maskEb);         cudaFreeHost(maskEc);
   cudaEventSynchronize(0);      cudaThreadSynchronize();

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

extern "C" void lanczos_real_fly_dynamic_cuda_(int *dimen_in, int *pblocksize, int* norbs_in, int *Niter_lanczos_, int* ntot, double* quart, 
        double* diag, double* subdiag, double* Eb_in, double* Ec_in, double* Vbc_in, int* sector_states_in, int* sector_ranks_in 
       ,int* bathnorbs_in,int* impnorbs_in,int* imporbs_in, int* bathorbs_in, int* maskEb_in, int* maskEc_in, int* maskVbc_in, double* vecinit)
{

  if(*dimen_in<max_rank){ use_texture_rank = 1;}else{use_texture_rank = 0;};

 //---------------------------------------------------------------------------------------//
     int verbose=0;
 //---------------------------------------------------------------------------------------//

     int Niter_lanczos=*Niter_lanczos_; int norbs=*norbs_in; int dimen=*dimen_in;
     
     if(verbose==1) printf(" start Lanczos Real on GPU \n" );

     int size = *ntot; int blocksize = *pblocksize;
     int bathnorbs=*bathnorbs_in; int impnorbs=*impnorbs_in;

     double *vec_in_gpu,*vec_out_gpu,*vec_tmp_gpu,*quart_gpu;
     double *vec_in_gpu_pointer,*vec_out_gpu_pointer,*vec_tmp_gpu_pointer,*quart_gpu_pointer;
     double *Eb,*Ec,*Vbc,*Eb_p,*Ec_p,*Vbc_p;

     int *sector_states_p,*sector_ranks_p,*imporbs_p,*bathorbs_p,*maskEb_p,*maskEc_p,*maskVbc_p;
     int *sector_states,*sector_ranks,*imporbs,*bathorbs,*maskEb,*maskEc,*maskVbc;

     if(verbose==1) printf(" GPU get eigenvalues, size of Lanczos vectors = %d ", size);

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

      cudaSetDevice(0); cudaSetDeviceFlags( cudaDeviceMapHost ); cudaEventSynchronize(0);

      if(use_texture_lrf==1){
         cudaMalloc((void**)&vec_in_gpu_pointer,sizeof(double)*size);
         cudaMalloc((void**)&Eb_p  , bathnorbs*bathnorbs*sizeof(double) );
         cudaMalloc((void**)&Ec_p  , impnorbs*impnorbs*sizeof(double) );
         cudaMalloc((void**)&Vbc_p , bathnorbs*impnorbs*sizeof(double) );
         cudaMemcpy(Eb_p,  Eb_in,  bathnorbs*bathnorbs*sizeof(double), cudaMemcpyHostToDevice);
         cudaMemcpy(Ec_p,  Ec_in,  impnorbs*impnorbs*sizeof(double),   cudaMemcpyHostToDevice);
         cudaMemcpy(Vbc_p, Vbc_in, bathnorbs*impnorbs*sizeof(double),  cudaMemcpyHostToDevice);
         bind_x(vec_in_gpu_pointer,size); bind_x_(Eb_p,Ec_p,Vbc_p);
       }else{
         cudaHostAlloc((void**)&vec_in_gpu    , size*sizeof(double),        cudaHostAllocMapped | cudaHostAllocPortable );
         cudaHostGetDevicePointer((void**)      &vec_in_gpu_pointer  ,  vec_in_gpu    , 0 );
         cudaHostAlloc((void**)&Eb            , bathnorbs*bathnorbs*sizeof(double), cudaHostAllocMapped | cudaHostAllocPortable );
         cudaHostAlloc((void**)&Ec            , impnorbs*impnorbs*sizeof(double),   cudaHostAllocMapped | cudaHostAllocPortable );
         cudaHostAlloc((void**)&Vbc           , bathnorbs*impnorbs*sizeof(double),  cudaHostAllocMapped | cudaHostAllocPortable );
         cudaMemcpy(Eb,  Eb_in,  bathnorbs*bathnorbs*sizeof(double), cudaMemcpyHostToHost);
         cudaMemcpy(Ec,  Ec_in,  impnorbs*impnorbs*sizeof(double),   cudaMemcpyHostToHost);
         cudaMemcpy(Vbc, Vbc_in, bathnorbs*impnorbs*sizeof(double),  cudaMemcpyHostToHost);
         cudaHostGetDevicePointer((void**) &Eb_p  , Eb  , 0 );
         cudaHostGetDevicePointer((void**) &Ec_p  , Ec  , 0 );
         cudaHostGetDevicePointer((void**) &Vbc_p , Vbc , 0 );
      }

      if(use_texture_rank==0){
       cudaHostAlloc((void**)&sector_ranks  , dimen*sizeof(int), cudaHostAllocMapped | cudaHostAllocPortable );
       cudaMemcpy(sector_ranks , sector_ranks_in , dimen*sizeof(int), cudaMemcpyHostToHost);
       cudaHostGetDevicePointer((void**)   &sector_ranks_p      , sector_ranks   , 0 );
      }else{
       cudaMalloc((void**)&sector_ranks_p  , dimen*sizeof(int) );
       cudaMemcpy(         sector_ranks_p , sector_ranks_in , dimen*sizeof(int), cudaMemcpyHostToDevice);
       bind_x_r(sector_ranks_p);
      };

       cudaHostAlloc((void**)&vec_tmp_gpu   , size*sizeof(double),               cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&vec_out_gpu   , size*sizeof(double),               cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&quart_gpu     , size*sizeof(double),               cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&sector_states , size*sizeof(int),                  cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&imporbs       , impnorbs*sizeof(int),              cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&bathorbs      , bathnorbs*sizeof(int),             cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&maskEb        , bathnorbs*bathnorbs*sizeof(int),   cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&maskEc        , impnorbs*impnorbs*sizeof(int),     cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&maskVbc       , bathnorbs*impnorbs*sizeof(int),    cudaHostAllocMapped | cudaHostAllocPortable );

       cudaEventSynchronize(0);
       cudaMemcpy(quart_gpu,   quart        , size*sizeof(double), cudaMemcpyHostToHost);
       cudaMemcpy(sector_states, sector_states_in, size*sizeof(int), cudaMemcpyHostToHost);
       cudaMemcpy(imporbs, imporbs_in , impnorbs*sizeof(int),  cudaMemcpyHostToHost);
       cudaMemcpy(bathorbs,bathorbs_in, bathnorbs*sizeof(int), cudaMemcpyHostToHost);
       cudaMemcpy(maskEb  ,maskEb_in , bathnorbs *bathnorbs  *sizeof(int), cudaMemcpyHostToHost);
       cudaMemcpy(maskEc  ,maskEc_in , impnorbs  *impnorbs   *sizeof(int), cudaMemcpyHostToHost);
       cudaMemcpy(maskVbc ,maskVbc_in, impnorbs  *bathnorbs  *sizeof(int), cudaMemcpyHostToHost);
       cudaEventSynchronize(0);

       cudaHostGetDevicePointer((void**)   &vec_tmp_gpu_pointer , vec_tmp_gpu    , 0 );
       cudaHostGetDevicePointer((void**)   &vec_out_gpu_pointer , vec_out_gpu    , 0 );
       cudaHostGetDevicePointer((void**)   &quart_gpu_pointer   , quart_gpu      , 0 );
       cudaHostGetDevicePointer((void**)   &sector_states_p     , sector_states  , 0 );
       cudaHostGetDevicePointer((void**)   &imporbs_p           , imporbs        , 0 );
       cudaHostGetDevicePointer((void**)   &bathorbs_p          , bathorbs       , 0 );
       cudaHostGetDevicePointer((void**)   &maskEb_p            , maskEb         , 0 );
       cudaHostGetDevicePointer((void**)   &maskEc_p            , maskEc         , 0 );
       cudaHostGetDevicePointer((void**)   &maskVbc_p           , maskVbc        , 0 );
       cudaEventSynchronize(0);

       cudaThreadSynchronize(); cudaEventSynchronize(0); 
       if(use_texture_lrf==0){
        cudaMemcpy(vec_in_gpu,vecinit,size*sizeof(double),cudaMemcpyHostToHost);
       }else{
        cudaMemcpy(vec_out_gpu,vecinit,size*sizeof(double),cudaMemcpyHostToHost);
        cudaMemcpy(vec_in_gpu_pointer,vec_out_gpu,size*sizeof(double),cudaMemcpyHostToDevice);
       }
       cudaThreadSynchronize(); cudaEventSynchronize(0);

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

  if(verbose==1) printf( " initialized, now run actual Lanczos \n " );

  for(int iter=0;iter<Niter_lanczos;iter++){

    if(verbose==1) printf( " Lanczos iteration %d / %d \n", iter,Niter_lanczos);

    one_step_lanczos_fly_cuda(blocksize,Niter_lanczos,size,iter,vec_tmp_gpu,vec_in_gpu,
                    vec_out_gpu,vec_tmp_gpu_pointer,vec_in_gpu_pointer,vec_out_gpu_pointer,
                    quart_gpu_pointer,diag,subdiag,Eb_p,Ec_p,Vbc_p,sector_states_p,sector_ranks_p,norbs,
                    bathnorbs,impnorbs,imporbs_p,bathorbs_p,maskEb_p,maskEc_p,maskVbc_p);
  };

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//
   cudaEventSynchronize(0); 
   if(use_texture_lrf==0)
   {
     cudaFreeHost(vec_in_gpu);cudaFreeHost(Eb);cudaFreeHost(Ec);cudaFreeHost(Vbc);
   }else{
     unbind_x(); unbind_x_();
     cudaFree(vec_in_gpu_pointer); cudaFree(Eb_p);cudaFree(Ec_p);cudaFree(Vbc_p);
   }
   if(use_texture_rank==0){ cudaFreeHost(sector_ranks);}else{unbind_x_r(); cudaFree(sector_ranks_p);};

   cudaFreeHost(vec_tmp_gpu);   cudaFreeHost(vec_out_gpu);   cudaFreeHost(quart_gpu);
   cudaFreeHost(sector_states); cudaFreeHost(bathorbs);      cudaFreeHost(imporbs);
   cudaFreeHost(maskEb);        cudaFreeHost(maskEc);        cudaFreeHost(maskVbc);
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

extern "C" void lanczos_real_fly_gs_cuda_(int *dimen_in, int *pblocksize, int* norbs_in, int *Niter_lanczos_, int* ntot, double* quart, 
                                       double* Eb_in, double* Ec_in, double* Vbc_in, int* sector_states_in, int* sector_ranks_in 
                                      ,int* bathnorbs_in,int* impnorbs_in,int* imporbs_in, int* bathorbs_in, int* maskEb_in, 
                                       int* maskEc_in, int* maskVbc_in,double *vecp,double *GS)
{

  if(*dimen_in<max_rank){ use_texture_rank = 1;}else{use_texture_rank = 0;};

 //---------------------------------------------------------------------------------------//
     int verbose=0;
 //---------------------------------------------------------------------------------------//

     int Niter_lanczos=*Niter_lanczos_; int norbs=*norbs_in; int dimen=*dimen_in;
     
     if(verbose==1) printf(" start Lanczos Real on GPU \n" );

     int size = *ntot; int blocksize = *pblocksize;
     int bathnorbs=*bathnorbs_in; int impnorbs=*impnorbs_in;

     double    *vec_in_gpu,*vec_out_gpu,*vec_tmp_gpu,*quart_gpu;
     double    *vec_in_gpu_pointer,*vec_out_gpu_pointer,*vec_tmp_gpu_pointer,*quart_gpu_pointer;
     double    *Eb,*Ec,*Vbc,*Eb_p,*Ec_p,*Vbc_p;
     int       *sector_states_p,*sector_ranks_p,*imporbs_p,*bathorbs_p;
     int       *maskEb_p,*maskEc_p,*maskVbc_p,*sector_states,*sector_ranks,*imporbs,*bathorbs;
     int       *maskEb,*maskEc,*maskVbc;
     double    *GS_gpu,*GS_gpu_pointer;

     double diag[Niter_lanczos], subdiag[Niter_lanczos];

     if(verbose==1) printf(" GPU get eigenvalues, size of Lanczos vectors = %d ", size);

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

       cudaSetDevice(0); cudaSetDeviceFlags( cudaDeviceMapHost ); cudaEventSynchronize(0);
       if(use_texture_lrf==1){
         cudaMalloc((void**)&vec_in_gpu_pointer,sizeof(double)*size);
         cudaMalloc((void**)&Eb_p  , bathnorbs*bathnorbs*sizeof(double) );
         cudaMalloc((void**)&Ec_p  , impnorbs*impnorbs*sizeof(double) );
         cudaMalloc((void**)&Vbc_p , bathnorbs*impnorbs*sizeof(double) );
         cudaMemcpy(Eb_p,  Eb_in,  bathnorbs*bathnorbs*sizeof(double), cudaMemcpyHostToDevice);
         cudaMemcpy(Ec_p,  Ec_in,  impnorbs*impnorbs*sizeof(double),   cudaMemcpyHostToDevice);
         cudaMemcpy(Vbc_p, Vbc_in, bathnorbs*impnorbs*sizeof(double),  cudaMemcpyHostToDevice);
         bind_x(vec_in_gpu_pointer,size);
         bind_x_(Eb_p,Ec_p,Vbc_p);
       }else{
         cudaHostAlloc((void**)&vec_in_gpu    , size*sizeof(double),        cudaHostAllocMapped | cudaHostAllocPortable );
         cudaHostGetDevicePointer((void**)      &vec_in_gpu_pointer  ,  vec_in_gpu    , 0 );
         cudaHostAlloc((void**)&Eb            , bathnorbs*bathnorbs*sizeof(double), cudaHostAllocMapped | cudaHostAllocPortable );
         cudaHostAlloc((void**)&Ec            , impnorbs*impnorbs*sizeof(double),   cudaHostAllocMapped | cudaHostAllocPortable );
         cudaHostAlloc((void**)&Vbc           , bathnorbs*impnorbs*sizeof(double),  cudaHostAllocMapped | cudaHostAllocPortable );
         cudaMemcpy(Eb,  Eb_in,  bathnorbs*bathnorbs*sizeof(double), cudaMemcpyHostToHost);
         cudaMemcpy(Ec,  Ec_in,  impnorbs*impnorbs*sizeof(double),   cudaMemcpyHostToHost);
         cudaMemcpy(Vbc, Vbc_in, bathnorbs*impnorbs*sizeof(double),  cudaMemcpyHostToHost);
         cudaHostGetDevicePointer((void**)  &Eb_p   , Eb  , 0 );
         cudaHostGetDevicePointer((void**)  &Ec_p   , Ec  , 0 );
         cudaHostGetDevicePointer((void**)  &Vbc_p  , Vbc , 0 );
       }

       if(use_texture_rank==0){
         cudaHostAlloc((void**)&sector_ranks  , dimen*sizeof(int), cudaHostAllocMapped | cudaHostAllocPortable );
         cudaMemcpy(sector_ranks , sector_ranks_in , dimen*sizeof(int), cudaMemcpyHostToHost);
         cudaHostGetDevicePointer((void**)   &sector_ranks_p      , sector_ranks   , 0 );
       }else{
         cudaMalloc((void**)&sector_ranks_p  , dimen*sizeof(int) );
         cudaMemcpy(         sector_ranks_p , sector_ranks_in , dimen*sizeof(int), cudaMemcpyHostToDevice);
         bind_x_r(sector_ranks_p);
       };

       cudaHostAlloc((void**)&vec_tmp_gpu   , size*sizeof(double),              cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&vec_out_gpu   , size*sizeof(double),              cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&quart_gpu     , size*sizeof(double),              cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&sector_states , size*sizeof(int),                 cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&imporbs       , impnorbs*sizeof(int),             cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&bathorbs      , bathnorbs*sizeof(int),            cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&maskEb        , bathnorbs*bathnorbs*sizeof(int),  cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&maskEc        , impnorbs*impnorbs*sizeof(int),    cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&maskVbc       , bathnorbs*impnorbs*sizeof(int),   cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&GS_gpu        , size*sizeof(double),              cudaHostAllocMapped | cudaHostAllocPortable );

       cudaEventSynchronize(0);
       cudaMemcpy(quart_gpu,   quart        , size*sizeof(double), cudaMemcpyHostToHost);
       cudaMemcpy(sector_states, sector_states_in, size*sizeof(int), cudaMemcpyHostToHost);
       cudaMemcpy(imporbs, imporbs_in , impnorbs*sizeof(int),  cudaMemcpyHostToHost);
       cudaMemcpy(bathorbs,bathorbs_in, bathnorbs*sizeof(int), cudaMemcpyHostToHost);
       cudaMemcpy(maskEb  ,maskEb_in , bathnorbs *bathnorbs  *sizeof(int), cudaMemcpyHostToHost);
       cudaMemcpy(maskEc  ,maskEc_in , impnorbs  *impnorbs   *sizeof(int), cudaMemcpyHostToHost);
       cudaMemcpy(maskVbc ,maskVbc_in, impnorbs  *bathnorbs  *sizeof(int), cudaMemcpyHostToHost);

       cudaEventSynchronize(0);
       cudaHostGetDevicePointer((void**)   &GS_gpu_pointer ,      GS_gpu        , 0 );
       cudaHostGetDevicePointer((void**)   &vec_tmp_gpu_pointer , vec_tmp_gpu   , 0 );
       cudaHostGetDevicePointer((void**)   &vec_out_gpu_pointer , vec_out_gpu   , 0 );
       cudaHostGetDevicePointer((void**)   &quart_gpu_pointer   , quart_gpu     , 0 );
       cudaHostGetDevicePointer((void**)   &sector_states_p     , sector_states , 0 );
       cudaHostGetDevicePointer((void**)   &imporbs_p           , imporbs       , 0 );
       cudaHostGetDevicePointer((void**)   &bathorbs_p          , bathorbs      , 0 );
       cudaHostGetDevicePointer((void**)   &maskEb_p            , maskEb        , 0 );
       cudaHostGetDevicePointer((void**)   &maskEc_p            , maskEc        , 0 );
       cudaHostGetDevicePointer((void**)   &maskVbc_p           , maskVbc       , 0 );
       cudaEventSynchronize(0);

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

   memset ((void **)GS_gpu, 0, size*sizeof(double));

    if(use_texture_lrf==0){
     for(int i=0;i<size;i++) vec_in_gpu[i]=1.0; cudaThreadSynchronize();
    }else{
     for(int i=0;i<size;i++) vec_out_gpu[i]=1.0; cudaThreadSynchronize();
     cudaMemcpy(vec_in_gpu_pointer,vec_out_gpu,size*sizeof(double),cudaMemcpyHostToDevice);
    }

   double coef= 1.0/sqrt(cublasDdot(size,vec_in_gpu_pointer,1,vec_in_gpu_pointer,1))*vecp[0]; cudaThreadSynchronize();  cudaEventSynchronize(0);
   cublasDaxpy(size,coef,vec_in_gpu_pointer,1,GS_gpu_pointer,1); cudaThreadSynchronize(); cudaEventSynchronize(0);cudaThreadSynchronize();

    one_step_lanczos_fly_cuda(blocksize,Niter_lanczos,size,0,vec_tmp_gpu,vec_in_gpu,
                    vec_out_gpu,vec_tmp_gpu_pointer,vec_in_gpu_pointer,vec_out_gpu_pointer,
                    quart_gpu_pointer,diag,subdiag,Eb_p,Ec_p,Vbc_p,sector_states_p,sector_ranks_p,norbs,
                    bathnorbs,impnorbs,imporbs_p,bathorbs_p,maskEb_p,maskEc_p,maskVbc_p);

  for(int iter=1;iter<Niter_lanczos-1;iter++){

   coef=1.0/sqrt(cublasDdot(size,vec_out_gpu_pointer,1,vec_out_gpu_pointer,1)); cublasDscal (size,coef,vec_out_gpu_pointer,1);  cudaEventSynchronize(0);
   cublasDaxpy(size,vecp[iter],vec_out_gpu_pointer,1,GS_gpu_pointer,1); cudaThreadSynchronize(); cudaEventSynchronize(0);

    one_step_lanczos_fly_cuda(blocksize,Niter_lanczos,size,iter,vec_tmp_gpu,vec_in_gpu,
                    vec_out_gpu,vec_tmp_gpu_pointer,vec_in_gpu_pointer,vec_out_gpu_pointer,
                    quart_gpu_pointer,diag,subdiag,Eb_p,Ec_p,Vbc_p,sector_states_p,sector_ranks_p,norbs,
                    bathnorbs,impnorbs,imporbs_p,bathorbs_p,maskEb_p,maskEc_p,maskVbc_p);

  };

   coef=1.0/sqrt(cublasDdot(size,vec_out_gpu_pointer,1,vec_out_gpu_pointer,1)); cublasDscal(size,coef,vec_out_gpu_pointer,1);cudaEventSynchronize(0);
   cublasDaxpy(size,vecp[Niter_lanczos-1],vec_out_gpu_pointer,1,GS_gpu_pointer,1); cudaThreadSynchronize(); cudaEventSynchronize(0);

   coef=1.0/sqrt(cublasDdot(size,GS_gpu_pointer,1,GS_gpu_pointer,1)); cudaThreadSynchronize();  cudaEventSynchronize(0);
   cublasDscal(size,coef,GS_gpu_pointer,1); cudaEventSynchronize(0); cudaThreadSynchronize();

   cudaMemcpy(GS,GS_gpu,size*sizeof(double),cudaMemcpyHostToHost);

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//
   cudaEventSynchronize(0); 
   if(use_texture_lrf==0)
   {
    cudaFreeHost(vec_in_gpu);cudaFreeHost(Eb);cudaFreeHost(Ec);cudaFreeHost(Vbc);
   }else{
    unbind_x(); unbind_x_();
    cudaFree(vec_in_gpu_pointer); cudaFree(Eb_p);cudaFree(Ec_p);cudaFree(Vbc_p);
   }
   if(use_texture_rank==0){ cudaFreeHost(sector_ranks);}else{unbind_x_r(); cudaFree(sector_ranks_p);};

   cudaFreeHost(GS_gpu);        cudaFreeHost(vec_tmp_gpu);  cudaFreeHost(vec_out_gpu);   cudaFreeHost(quart_gpu);
   cudaFreeHost(sector_states); 
   cudaFreeHost(bathorbs);      cudaFreeHost(imporbs);
   cudaFreeHost(maskEb);        cudaFreeHost(maskEc);       cudaFreeHost(maskVbc);
 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//
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
