#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <cublas.h>
#include <cuComplex.h>

//********************************************
//********************************************
//********************************************

#define IBSET(a,b) ((a) |= (1<<(b)))
#define IBCLR(a,b) ((a) &= ~(1<<(b)))
#define BTEST(a,b) ((a) & (1<<(b)))

//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************

const int use_texture_lcf = 1 ;
      int use_texture_rank_c = 1 ;

#define WARP_SIZE 32
#define MAX_BLOCK 65500
#define max_rank 225000000

texture<int2,1,cudaReadModeElementType> tex,texb,texc,texbc;
texture<int,1,cudaReadModeElementType>  texr;

inline void   bind_x_r(int *x)         {cudaBindTexture(0,texr,x);};
inline void unbind_x_r()               {cudaUnbindTexture(texr);};

inline void   bind_x(cuDoubleComplex *x, int N) {cudaBindTexture(0,tex,x,N*sizeof(cuDoubleComplex));};
inline void unbind_x()                          {cudaUnbindTexture(tex);};

inline void   bind_x_(cuDoubleComplex *x, cuDoubleComplex *y, cuDoubleComplex *z ) 
                                       {cudaBindTexture(0,texb,x); cudaBindTexture(0,texc,y); cudaBindTexture(0,texbc,z);};
inline void unbind_x_()                {cudaUnbindTexture(texb);cudaUnbindTexture(texc);cudaUnbindTexture(texbc);};

__inline__  __device__ int    fetch_r(const int& i) {  int  v = tex1Dfetch(texr,i);  return v; }

__inline__  __device__ cuDoubleComplex fetch_x(const int& i)
  { int  jj = 2*(i); int2 v  = tex1Dfetch(tex,jj); double rr  = __hiloint2double(v.y, v.x); v  = tex1Dfetch(tex,jj+1);
    double im  = __hiloint2double(v.y, v.x); return make_cuDoubleComplex(rr,im); }
__inline__  __device__ cuDoubleComplex fetchb(const int& i)
  { int  jj = 2*(i); int2 v  = tex1Dfetch(texb,jj); double rr  = __hiloint2double(v.y, v.x); v  = tex1Dfetch(texb,jj+1);
    double im  = __hiloint2double(v.y, v.x); return make_cuDoubleComplex(rr,im); }
__inline__  __device__ cuDoubleComplex fetchc(const int& i)
  { int  jj = 2*(i); int2 v  = tex1Dfetch(texc,jj); double rr  = __hiloint2double(v.y, v.x); v  = tex1Dfetch(texc,jj+1);
    double im  = __hiloint2double(v.y, v.x); return make_cuDoubleComplex(rr,im); }
__inline__  __device__ cuDoubleComplex fetchbc(const int& i)
  { int  jj = 2*(i); int2 v  = tex1Dfetch(texbc,jj); double rr  = __hiloint2double(v.y, v.x); v  = tex1Dfetch(texbc,jj+1);
    double im  = __hiloint2double(v.y, v.x); return make_cuDoubleComplex(rr,im); }

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
//********************************************
//********************************************
//********************************************
//********************************************

__global__ void norm_vec_ker_complex(int size , cuDoubleComplex *x, double *normv)

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

__global__ void real_scal_vec_ker_complex( int size , cuDoubleComplex *x, cuDoubleComplex *y, double *normv)

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

__global__ void normalize_vec_ker_complex(int ngrid, int size ,cuDoubleComplex *x, double *normv)
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

  //************************************************//
  //               Kernel Hmult                     //
  //************************************************//

  __global__  void Hmult_ker_complex_fly(int ngrid, int blocksize, int psize, cuDoubleComplex *vec_out, const cuDoubleComplex *vec_in, 
    const double *quart, const cuDoubleComplex* Eb, const cuDoubleComplex* Ec, const cuDoubleComplex* Vbc, 
    const int* sector_states, const int* sector_rank, const int norbs, const int bathnorbs, const int impnorbs, 
    const int* imporbs_, const int* bathorbs_, const int* maskEb_, const int* maskEc_, const int* maskVbc_, const int use_texture_lcf,
    const int use_texture_rank_c)

{
   __shared__ cuDoubleComplex sdata[16][WARP_SIZE];
   __shared__ short int maskEc[32*32],maskEb[32*32],maskVbc[32*32];
   __shared__ short int imporbs[32],bathorbs[32];

    int istate ;  const int warp_lane   = threadIdx.y; const int thread_lane = threadIdx.x; 

    int jstate,kets_out,kets_in,jj,iorb,jorb,n1,n2; 
    cuDoubleComplex hoffdiag,vecin;
    short int fermion_sign;

    if(warp_lane==0) for(iorb=thread_lane;iorb<impnorbs*impnorbs;iorb+=WARP_SIZE)   maskEc[iorb]=maskEc_[iorb];
    if(warp_lane==1) for(iorb=thread_lane;iorb<bathnorbs*bathnorbs;iorb+=WARP_SIZE) maskEb[iorb]=maskEb_[iorb];
    if(warp_lane==2) for(iorb=thread_lane;iorb<bathnorbs*impnorbs;iorb+=WARP_SIZE)  maskVbc[iorb]=maskVbc_[iorb];
    if(warp_lane==3) for(iorb=thread_lane;iorb<bathnorbs;iorb+=WARP_SIZE)           bathorbs[iorb]=bathorbs_[iorb];
    if(warp_lane==4) for(iorb=thread_lane;iorb<impnorbs;iorb+=WARP_SIZE)            imporbs[iorb]=imporbs_[iorb];

    __syncthreads();


  for(int iii=0;iii<=ngrid;iii++){

    istate = blocksize * ( blockIdx.y + iii*MAX_BLOCK )  + threadIdx.y ;

    if(istate < psize)
  {

    if(use_texture_lcf==0){ vecin =  vec_in[istate];
                     }else{ vecin = fetch_x(istate); };

    kets_in = sector_states[istate];
    if(thread_lane==0){
      sdata[threadIdx.y][threadIdx.x] = cuCmul(make_cuDoubleComplex(quart[istate],0.0) , vecin);
    }else{
      sdata[threadIdx.y][threadIdx.x] = make_cuDoubleComplex(0.0,0.0) ;
    };

    //////////////////////////////////////////////////////////////////////////

   if(use_texture_lcf==0){
    for(iorb=thread_lane;iorb<impnorbs ;iorb+=WARP_SIZE)
      {if(BTEST(kets_in,imporbs[iorb]-1 )>0) sdata[threadIdx.y][threadIdx.x]=cuCadd(sdata[threadIdx.y][threadIdx.x],cuCmul(Ec[iorb* impnorbs+iorb],vecin));};
    for(iorb=thread_lane;iorb<bathnorbs;iorb+=WARP_SIZE)
      {if(BTEST(kets_in,bathorbs[iorb]-1)>0) sdata[threadIdx.y][threadIdx.x]=cuCadd(sdata[threadIdx.y][threadIdx.x],cuCmul(Eb[iorb*bathnorbs+iorb],vecin));};
   }else{
    for(iorb=thread_lane;iorb<impnorbs ;iorb+=WARP_SIZE)
      {if(BTEST(kets_in,imporbs[iorb]-1 )>0) sdata[threadIdx.y][threadIdx.x]=cuCadd(sdata[threadIdx.y][threadIdx.x],cuCmul(fetchc(iorb* impnorbs+iorb),vecin));};
    for(iorb=thread_lane;iorb<bathnorbs;iorb+=WARP_SIZE)
      {if(BTEST(kets_in,bathorbs[iorb]-1)>0) sdata[threadIdx.y][threadIdx.x]=cuCadd(sdata[threadIdx.y][threadIdx.x],cuCmul(fetchb(iorb*bathnorbs+iorb),vecin));};
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
           if(use_texture_rank_c==0){
             jstate                            = sector_rank[kets_out]-1;}else{
             jstate                            = fetch_r(kets_out)-1;};
           if(use_texture_lcf==0){
             hoffdiag                          = Ec[jorb*impnorbs+iorb];
             hoffdiag                          = make_cuDoubleComplex(cuCreal(hoffdiag)*fermion_sign,-cuCimag(hoffdiag)*fermion_sign );
             sdata[threadIdx.y][threadIdx.x]   =cuCadd(sdata[threadIdx.y][threadIdx.x], cuCmul(hoffdiag,vec_in[jstate]));
           }else{
             hoffdiag                          = fetchc(jorb*impnorbs+iorb);
             hoffdiag                          = make_cuDoubleComplex(cuCreal(hoffdiag)*fermion_sign,-cuCimag(hoffdiag)*fermion_sign );
             sdata[threadIdx.y][threadIdx.x]   =cuCadd(sdata[threadIdx.y][threadIdx.x], cuCmul(fetch_x(jstate),hoffdiag));
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
         if(use_texture_rank_c==0){
            jstate                            = sector_rank[kets_out]-1;}else{
            jstate                            = fetch_r(kets_out)-1;};
            if(use_texture_lcf==0){
             hoffdiag                          = Eb[jorb*bathnorbs+iorb];
             hoffdiag                          = make_cuDoubleComplex(cuCreal(hoffdiag)*fermion_sign,-cuCimag(hoffdiag)*fermion_sign );
             sdata[threadIdx.y][threadIdx.x]   = cuCadd(sdata[threadIdx.y][threadIdx.x],cuCmul(hoffdiag,vec_in[jstate]));
           }else{
             hoffdiag                          = fetchb(jorb*bathnorbs+iorb);
             hoffdiag                          = make_cuDoubleComplex(cuCreal(hoffdiag)*fermion_sign,-cuCimag(hoffdiag)*fermion_sign );
             sdata[threadIdx.y][threadIdx.x]   =cuCadd(sdata[threadIdx.y][threadIdx.x], cuCmul(fetch_x(jstate),hoffdiag));
           };
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
          if(use_texture_rank_c==0){
            jstate                            = sector_rank[kets_out]-1;}else{
            jstate                            = fetch_r(kets_out)-1;};
            if(use_texture_lcf==0){
             hoffdiag                          = Vbc[jorb*bathnorbs+iorb];
             hoffdiag                          = make_cuDoubleComplex(cuCreal(hoffdiag)*fermion_sign,-cuCimag(hoffdiag)*fermion_sign );
             sdata[threadIdx.y][threadIdx.x]   =cuCadd(sdata[threadIdx.y][threadIdx.x], cuCmul(hoffdiag,vec_in[jstate]));
           }else{
             hoffdiag                          = fetchbc(jorb*bathnorbs+iorb);
             hoffdiag                          = make_cuDoubleComplex(cuCreal(hoffdiag)*fermion_sign,-cuCimag(hoffdiag)*fermion_sign );
             sdata[threadIdx.y][threadIdx.x]   =cuCadd(sdata[threadIdx.y][threadIdx.x], cuCmul(fetch_x(jstate),hoffdiag));
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
          if(use_texture_rank_c==0){
            jstate                            = sector_rank[kets_out]-1;}else{
            jstate                            = fetch_r(kets_out)-1;};
           if(use_texture_lcf==0){
             hoffdiag                          = Vbc[jorb*bathnorbs+iorb];
             hoffdiag                          = make_cuDoubleComplex(cuCreal(hoffdiag)*fermion_sign,-cuCimag(hoffdiag)*fermion_sign );
             sdata[threadIdx.y][threadIdx.x]   =cuCadd(sdata[threadIdx.y][threadIdx.x], cuCmul(hoffdiag,vec_in[jstate]));
           }else{
             hoffdiag                          = fetchbc(jorb*bathnorbs+iorb);
             hoffdiag                          = make_cuDoubleComplex(cuCreal(hoffdiag)*fermion_sign,-cuCimag(hoffdiag)*fermion_sign );
             sdata[threadIdx.y][threadIdx.x]   =cuCadd(sdata[threadIdx.y][threadIdx.x], cuCmul(fetch_x(jstate),hoffdiag));
           }
        }; }; }; }; };
    //////////////////////////////////////////////////////////////////////////

        if (thread_lane < 16) { sdata[threadIdx.y][threadIdx.x] = cuCadd(sdata[threadIdx.y][threadIdx.x], sdata[threadIdx.y][threadIdx.x + 16] ); };
        if (thread_lane <  8) { sdata[threadIdx.y][threadIdx.x] = cuCadd(sdata[threadIdx.y][threadIdx.x], sdata[threadIdx.y][threadIdx.x +  8] ); };
        if (thread_lane <  4) { sdata[threadIdx.y][threadIdx.x] = cuCadd(sdata[threadIdx.y][threadIdx.x], sdata[threadIdx.y][threadIdx.x +  4] ); };
        if (thread_lane <  2) { sdata[threadIdx.y][threadIdx.x] = cuCadd(sdata[threadIdx.y][threadIdx.x], sdata[threadIdx.y][threadIdx.x +  2] ); };
        if (thread_lane <  1) { sdata[threadIdx.y][threadIdx.x] = cuCadd(sdata[threadIdx.y][threadIdx.x], sdata[threadIdx.y][threadIdx.x +  1] ); };
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

 void one_step_lanczos_complex_fly_cuda(int blocksize, int Niter, int size, int iter, cuDoubleComplex *vec_tmp_gpu, cuDoubleComplex *vec_in_gpu, 
     cuDoubleComplex *vec_out_gpu, cuDoubleComplex *vec_tmp_gpu_pointer, cuDoubleComplex *vec_in_gpu_pointer, cuDoubleComplex *vec_out_gpu_pointer,
     double *quart_gpu_pointer , double *diag, double *subdiag , cuDoubleComplex* Eb, cuDoubleComplex* Ec, cuDoubleComplex* Vbc,
     int* sector_states, int* sector_ranks,int norbs, int bathnorbs,int impnorbs,int* imporbs, int* bathorbs, int* maskEb, int* maskEc, int* maskVbc )

{
   int verbose=0; int psize=size; if(verbose==1) printf( " define blocksize \n ");
   double normv; double *normv_ker; cudaMalloc((void**)&normv_ker,sizeof(double)); double *normv_loc; normv_loc = &normv; cuDoubleComplex coef;

   int nb=(size -size % blocksize)/blocksize + 1; int ngrid=nb/MAX_BLOCK; if(ngrid>0) nb=MAX_BLOCK;
   dim3 bl(1,nb),th(WARP_SIZE,blocksize);
   if(verbose==1) printf( " --------------- \n Nblock=%d Ngrid=%d \n ----------------- \n ",nb,ngrid);

   int nb2 =(size-size % 256 )/256 + 1; 
   int ngrid2=nb2/MAX_BLOCK; if(ngrid2>0) nb2=MAX_BLOCK;

   norm_vec_ker_complex<<<1,512>>>(size,vec_in_gpu_pointer,normv_ker); cudaEventSynchronize(0); cudaThreadSynchronize();
   normalize_vec_ker_complex<<<nb2,256>>>(ngrid2, size, vec_in_gpu_pointer, normv_ker);
 
   if(verbose==1) printf( " one step Lanczos, norm vecin=%f, size of problem=%d \n " , normv, size); 
   if(verbose==1) printf( " use texture optimization = %d " , use_texture_lcf );

   Hmult_ker_complex_fly<<<bl,th>>>(ngrid,blocksize,psize,vec_out_gpu_pointer,vec_in_gpu_pointer,quart_gpu_pointer,Eb,Ec,Vbc,sector_states,
                                   sector_ranks,norbs,bathnorbs,impnorbs,imporbs,bathorbs,maskEb,maskEc,maskVbc,use_texture_lcf,
                                   use_texture_rank_c); 
   cudaEventSynchronize(0); cudaThreadSynchronize();

  if(iter>0){ coef=make_cuDoubleComplex(-subdiag[iter],0.0); cublasZaxpy_no_device(size,coef,vec_tmp_gpu,1,vec_out_gpu,1); }; cudaEventSynchronize(0); cudaThreadSynchronize();

  if(use_texture_lcf==0){
   cudaMemcpy(vec_tmp_gpu,vec_in_gpu,size*sizeof(cuDoubleComplex),cudaMemcpyHostToHost);cudaEventSynchronize(0);cudaThreadSynchronize();
  }else{
   cudaMemcpy(vec_tmp_gpu,vec_in_gpu_pointer,size*sizeof(cuDoubleComplex),cudaMemcpyDeviceToHost);cudaEventSynchronize(0);cudaThreadSynchronize();
  }

   real_scal_vec_ker_complex<<<1,512>>> (size,vec_out_gpu_pointer,vec_in_gpu_pointer,normv_ker);
   cudaMemcpy(normv_loc,normv_ker,sizeof(double),cudaMemcpyDeviceToHost);
   diag[iter]=*normv_loc ;

   cudaEventSynchronize(0);cudaThreadSynchronize();
   coef=make_cuDoubleComplex(-diag[iter],0.0);
   cublasZaxpy_no_device(size,coef,vec_tmp_gpu,1,vec_out_gpu,1); cudaEventSynchronize(0); cudaThreadSynchronize();

   normv = sqrt(cuCabs(cublasZdotu_no_device(size,vec_out_gpu,1,vec_out_gpu,1))); cudaEventSynchronize(0); cudaThreadSynchronize();

   if(iter<Niter-1) subdiag[iter+1]=normv;

  if(use_texture_lcf==0){
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

extern "C" void lanczos_complex_fly_cuda_(int *dimen_in, int *pblocksize, int* norbs_in, int *Niter_lanczos_, int* ntot, double* quart, 
      double* diag, double* subdiag, cuDoubleComplex* Eb_in, cuDoubleComplex* Ec_in, cuDoubleComplex* Vbc_in, int* sector_states_in, int* sector_ranks_in 
     ,int* bathnorbs_in,int* impnorbs_in,int* imporbs_in, int* bathorbs_in, int* maskEb_in, int* maskEc_in, int* maskVbc_in)
{

  if(*dimen_in<max_rank){ use_texture_rank_c = 1;}else{use_texture_rank_c = 0;};

 //---------------------------------------------------------------------------------------//
     int verbose=0;
 //---------------------------------------------------------------------------------------//

     int Niter_lanczos=*Niter_lanczos_; int norbs=*norbs_in; int dimen=*dimen_in;
     
     if(verbose==1) printf(" start Lanczos Real on GPU \n" );

     int size      = *ntot; int blocksize = *pblocksize;
     int bathnorbs = *bathnorbs_in; int impnorbs=*impnorbs_in;

     cuDoubleComplex  *vec_in_gpu,*vec_out_gpu,*vec_tmp_gpu;
     double           *quart_gpu;
     cuDoubleComplex  *vec_in_gpu_pointer,*vec_out_gpu_pointer,*vec_tmp_gpu_pointer;
     double           *quart_gpu_pointer;
     cuDoubleComplex  *Eb,*Ec,*Vbc,*Eb_p,*Ec_p,*Vbc_p;

     int *sector_states_p,*sector_ranks_p,*imporbs_p,*bathorbs_p,*maskEb_p,*maskEc_p,*maskVbc_p;
     int *sector_states,*sector_ranks,*imporbs,*bathorbs,*maskEb,*maskEc,*maskVbc;

     if(verbose==1) printf(" GPU get eigenvalues, size of Lanczos vectors = %d ", size);

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

       cudaSetDevice(0); cudaSetDeviceFlags( cudaDeviceMapHost ); cudaEventSynchronize(0);

       if(use_texture_lcf==1){
         cudaMalloc((void**)&vec_in_gpu_pointer,sizeof(cuDoubleComplex)*size);
         cudaMalloc((void**)&Eb_p  , bathnorbs*bathnorbs*sizeof(cuDoubleComplex) );
         cudaMalloc((void**)&Ec_p  ,  impnorbs*impnorbs*sizeof(cuDoubleComplex) );
         cudaMalloc((void**)&Vbc_p , bathnorbs*impnorbs*sizeof(cuDoubleComplex) );
         cudaMemcpy(Eb_p  ,  Eb_in,  bathnorbs*bathnorbs*sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
         cudaMemcpy(Ec_p  ,  Ec_in,   impnorbs*impnorbs*sizeof(cuDoubleComplex),  cudaMemcpyHostToDevice);
         cudaMemcpy(Vbc_p , Vbc_in,  bathnorbs*impnorbs*sizeof(cuDoubleComplex),  cudaMemcpyHostToDevice);
         bind_x(vec_in_gpu_pointer,size);
         bind_x_(Eb_p,Ec_p,Vbc_p);
       }else{
         cudaHostAlloc((void**)&vec_in_gpu , size*sizeof(cuDoubleComplex)  , cudaHostAllocMapped | cudaHostAllocPortable );
         cudaHostGetDevicePointer((void**)   &vec_in_gpu_pointer  , vec_in_gpu , 0 );
         cudaHostAlloc((void**)&Eb         , bathnorbs*bathnorbs*sizeof(cuDoubleComplex), cudaHostAllocMapped | cudaHostAllocPortable );
         cudaHostAlloc((void**)&Ec         , impnorbs*impnorbs*sizeof(cuDoubleComplex)  , cudaHostAllocMapped | cudaHostAllocPortable );
         cudaHostAlloc((void**)&Vbc        , bathnorbs*impnorbs*sizeof(cuDoubleComplex) , cudaHostAllocMapped | cudaHostAllocPortable );
         cudaMemcpy(Eb,  Eb_in,  bathnorbs*bathnorbs*sizeof(cuDoubleComplex), cudaMemcpyHostToHost);
         cudaMemcpy(Ec,  Ec_in,  impnorbs*impnorbs*sizeof(cuDoubleComplex),   cudaMemcpyHostToHost);
         cudaMemcpy(Vbc, Vbc_in, bathnorbs*impnorbs*sizeof(cuDoubleComplex),  cudaMemcpyHostToHost);
         cudaHostGetDevicePointer((void**)   &Eb_p  , Eb  , 0 );
         cudaHostGetDevicePointer((void**)   &Ec_p  , Ec  , 0 );
         cudaHostGetDevicePointer((void**)   &Vbc_p , Vbc , 0 );
       }

      if(verbose==1) printf( " bind sector_ranks to texture if necessary \n ");

      if(use_texture_rank_c==0){
       cudaHostAlloc((void**)&sector_ranks , dimen*sizeof(int), cudaHostAllocMapped | cudaHostAllocPortable );
       cudaMemcpy(sector_ranks , sector_ranks_in , dimen*sizeof(int), cudaMemcpyHostToHost);
       cudaHostGetDevicePointer((void**)  &sector_ranks_p      , sector_ranks   , 0 );
      }else{
       cudaMalloc((void**)&sector_ranks_p , dimen*sizeof(int) );
       cudaMemcpy(         sector_ranks_p , sector_ranks_in , dimen*sizeof(int), cudaMemcpyHostToDevice);
       bind_x_r(sector_ranks_p);
      };

       cudaHostAlloc((void**)&vec_tmp_gpu   , size*sizeof(cuDoubleComplex),       cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&vec_out_gpu   , size*sizeof(cuDoubleComplex),       cudaHostAllocMapped | cudaHostAllocPortable );
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

 if(use_texture_lcf==0){
  for(int i=0;i<size;i++) vec_in_gpu[i]=make_cuDoubleComplex(1.0,0.0); cudaThreadSynchronize();
 }else{
   for(int i=0;i<size;i++) vec_out_gpu[i]=make_cuDoubleComplex(1.0,0.0); cudaThreadSynchronize();
   cudaMemcpy(vec_in_gpu_pointer,vec_out_gpu,size*sizeof(cuDoubleComplex),cudaMemcpyHostToDevice);
 }

  for(int iter=0;iter<Niter_lanczos;iter++){

    if(verbose==1) printf( " Lanczos iteration %d / %d \n", iter,Niter_lanczos);

    one_step_lanczos_complex_fly_cuda(blocksize,Niter_lanczos,size,iter,vec_tmp_gpu,vec_in_gpu,
                    vec_out_gpu,vec_tmp_gpu_pointer,vec_in_gpu_pointer,vec_out_gpu_pointer,
                    quart_gpu_pointer,diag,subdiag,Eb_p,Ec_p,Vbc_p,sector_states_p,sector_ranks_p,norbs,
                    bathnorbs,impnorbs,imporbs_p,bathorbs_p,maskEb_p,maskEc_p,maskVbc_p);
  };

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

   cudaEventSynchronize(0); 
   if(use_texture_lcf==0)
   {
    cudaFreeHost(vec_in_gpu);cudaFreeHost(Eb);cudaFreeHost(Ec);cudaFreeHost(Vbc);
   }else{
    unbind_x(); unbind_x_(); cudaFree(vec_in_gpu_pointer); cudaFree(Eb_p);cudaFree(Ec_p);cudaFree(Vbc_p);
   }

   if(use_texture_rank_c==0){ cudaFreeHost(sector_ranks);}else{unbind_x_r(); cudaFree(sector_ranks_p);};

   cudaFreeHost(vec_tmp_gpu);    cudaFreeHost(vec_out_gpu); 
   cudaFreeHost(sector_states);  cudaFreeHost(quart_gpu);
   cudaFreeHost(bathorbs);       cudaFreeHost(imporbs);        cudaFreeHost(maskVbc);
   cudaFreeHost(maskEb);         cudaFreeHost(maskEc);
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

extern "C" void lanczos_complex_fly_dynamic_cuda_(int *dimen_in, int *pblocksize, int* norbs_in, int *Niter_lanczos_, int* ntot, double* quart, 
        double* diag, double* subdiag, cuDoubleComplex* Eb_in, cuDoubleComplex* Ec_in, cuDoubleComplex* Vbc_in, int* sector_states_in, int* sector_ranks_in 
       ,int* bathnorbs_in,int* impnorbs_in,int* imporbs_in, int* bathorbs_in, int* maskEb_in, int* maskEc_in, int* maskVbc_in, cuDoubleComplex* vecinit)
{

  if(*dimen_in<max_rank){ use_texture_rank_c = 1;}else{use_texture_rank_c = 0;};

 //---------------------------------------------------------------------------------------//
     int verbose=0;
 //---------------------------------------------------------------------------------------//

     int Niter_lanczos=*Niter_lanczos_; int norbs=*norbs_in; int dimen=*dimen_in;
     
     if(verbose==1) printf(" start Lanczos Real on GPU \n" );

     int size = *ntot; int blocksize = *pblocksize; int bathnorbs=*bathnorbs_in; int impnorbs=*impnorbs_in;

     cuDoubleComplex *vec_in_gpu,*vec_out_gpu,*vec_tmp_gpu;
     double          *quart_gpu;
     cuDoubleComplex *vec_in_gpu_pointer,*vec_out_gpu_pointer,*vec_tmp_gpu_pointer;
     double          *quart_gpu_pointer;
     cuDoubleComplex *Eb,*Ec,*Vbc,*Eb_p,*Ec_p,*Vbc_p;

     int *sector_states_p,*sector_ranks_p,*imporbs_p,*bathorbs_p,*maskEb_p,*maskEc_p,*maskVbc_p;
     int *sector_states,*sector_ranks,*imporbs,*bathorbs,*maskEb,*maskEc,*maskVbc;

     if(verbose==1) printf(" GPU get eigenvalues, size of Lanczos vectors = %d ", size);

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

      cudaSetDevice(0); cudaSetDeviceFlags( cudaDeviceMapHost ); cudaEventSynchronize(0);

      if(use_texture_lcf==1){
         cudaMalloc((void**)&vec_in_gpu_pointer,sizeof(cuDoubleComplex)*size);
         cudaMalloc((void**)&Eb_p  , bathnorbs*bathnorbs*sizeof(cuDoubleComplex) );
         cudaMalloc((void**)&Ec_p  , impnorbs*impnorbs*sizeof(cuDoubleComplex) );
         cudaMalloc((void**)&Vbc_p , bathnorbs*impnorbs*sizeof(cuDoubleComplex) );
         cudaMemcpy(Eb_p,  Eb_in,  bathnorbs*bathnorbs*sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
         cudaMemcpy(Ec_p,  Ec_in,  impnorbs*impnorbs*sizeof(cuDoubleComplex),   cudaMemcpyHostToDevice);
         cudaMemcpy(Vbc_p, Vbc_in, bathnorbs*impnorbs*sizeof(cuDoubleComplex),  cudaMemcpyHostToDevice);
         bind_x(vec_in_gpu_pointer,size); bind_x_(Eb_p,Ec_p,Vbc_p);
       }else{
         cudaHostAlloc((void**)&vec_in_gpu    , size*sizeof(cuDoubleComplex),        cudaHostAllocMapped | cudaHostAllocPortable );
         cudaHostGetDevicePointer((void**)      &vec_in_gpu_pointer  ,  vec_in_gpu    , 0 );
         cudaHostAlloc((void**)&Eb            , bathnorbs*bathnorbs*sizeof(cuDoubleComplex), cudaHostAllocMapped | cudaHostAllocPortable );
         cudaHostAlloc((void**)&Ec            , impnorbs*impnorbs*sizeof(cuDoubleComplex),   cudaHostAllocMapped | cudaHostAllocPortable );
         cudaHostAlloc((void**)&Vbc           , bathnorbs*impnorbs*sizeof(cuDoubleComplex),  cudaHostAllocMapped | cudaHostAllocPortable );
         cudaMemcpy(Eb,  Eb_in,  bathnorbs*bathnorbs*sizeof(cuDoubleComplex), cudaMemcpyHostToHost);
         cudaMemcpy(Ec,  Ec_in,  impnorbs*impnorbs*sizeof(cuDoubleComplex),   cudaMemcpyHostToHost);
         cudaMemcpy(Vbc, Vbc_in, bathnorbs*impnorbs*sizeof(cuDoubleComplex),  cudaMemcpyHostToHost);
         cudaHostGetDevicePointer((void**) &Eb_p  , Eb  , 0 );
         cudaHostGetDevicePointer((void**) &Ec_p  , Ec  , 0 );
         cudaHostGetDevicePointer((void**) &Vbc_p , Vbc , 0 );
      }

      if(use_texture_rank_c==0){
       cudaHostAlloc((void**)&sector_ranks  , dimen*sizeof(int), cudaHostAllocMapped | cudaHostAllocPortable );
       cudaMemcpy(sector_ranks , sector_ranks_in , dimen*sizeof(int), cudaMemcpyHostToHost);
       cudaHostGetDevicePointer((void**)   &sector_ranks_p      , sector_ranks   , 0 );
      }else{
       cudaMalloc((void**)&sector_ranks_p  , dimen*sizeof(int) );
       cudaMemcpy(         sector_ranks_p , sector_ranks_in , dimen*sizeof(int), cudaMemcpyHostToDevice);
       bind_x_r(sector_ranks_p);
      };

       cudaHostAlloc((void**)&vec_tmp_gpu   , size*sizeof(cuDoubleComplex),       cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&vec_out_gpu   , size*sizeof(cuDoubleComplex),       cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&quart_gpu     , size*sizeof(double),                cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&sector_states , size*sizeof(int),                   cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&imporbs       , impnorbs*sizeof(int),               cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&bathorbs      , bathnorbs*sizeof(int),              cudaHostAllocMapped | cudaHostAllocPortable );
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
       if(use_texture_lcf==0){
        cudaMemcpy(vec_in_gpu,vecinit,size*sizeof(cuDoubleComplex),cudaMemcpyHostToHost);
       }else{
        cudaMemcpy(vec_out_gpu,vecinit,size*sizeof(cuDoubleComplex),cudaMemcpyHostToHost);
        cudaMemcpy(vec_in_gpu_pointer,vec_out_gpu,size*sizeof(cuDoubleComplex),cudaMemcpyHostToDevice);
       }
       cudaThreadSynchronize(); cudaEventSynchronize(0);

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

  if(verbose==1) printf( " initialized, now run actual Lanczos \n " );

  for(int iter=0;iter<Niter_lanczos;iter++){

    if(verbose==1) printf( " Lanczos iteration %d / %d \n", iter,Niter_lanczos);

    one_step_lanczos_complex_fly_cuda(blocksize,Niter_lanczos,size,iter,vec_tmp_gpu,vec_in_gpu,
                    vec_out_gpu,vec_tmp_gpu_pointer,vec_in_gpu_pointer,vec_out_gpu_pointer,
                    quart_gpu_pointer,diag,subdiag,Eb_p,Ec_p,Vbc_p,sector_states_p,sector_ranks_p,norbs,
                    bathnorbs,impnorbs,imporbs_p,bathorbs_p,maskEb_p,maskEc_p,maskVbc_p);
  };

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//
   cudaEventSynchronize(0); 
   if(use_texture_lcf==0)
   {
     cudaFreeHost(vec_in_gpu);cudaFreeHost(Eb);cudaFreeHost(Ec);cudaFreeHost(Vbc);
   }else{
     unbind_x(); unbind_x_();
     cudaFree(vec_in_gpu_pointer); cudaFree(Eb_p);cudaFree(Ec_p);cudaFree(Vbc_p);
   }
   if(use_texture_rank_c==0){ cudaFreeHost(sector_ranks);}else{unbind_x_r(); cudaFree(sector_ranks_p);};

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

extern "C" void lanczos_complex_fly_gs_cuda_(int *dimen_in, int *pblocksize, int* norbs_in, int *Niter_lanczos_, int* ntot, double* quart, 
  cuDoubleComplex* Eb_in, cuDoubleComplex* Ec_in, cuDoubleComplex* Vbc_in, int* sector_states_in, int* sector_ranks_in 
 ,int* bathnorbs_in,int* impnorbs_in,int* imporbs_in, int* bathorbs_in, int* maskEb_in, int* maskEc_in, int* maskVbc_in,double *vecp,cuDoubleComplex *GS)
{

  if(*dimen_in<max_rank){ use_texture_rank_c = 1;}else{use_texture_rank_c = 0;};

 //---------------------------------------------------------------------------------------//
     int verbose=0;
 //---------------------------------------------------------------------------------------//

     int Niter_lanczos=*Niter_lanczos_; int norbs=*norbs_in; int dimen=*dimen_in;
     
     if(verbose==1) printf(" start Lanczos Real on GPU \n" );

     int size = *ntot; int blocksize = *pblocksize;
     int bathnorbs=*bathnorbs_in; int impnorbs=*impnorbs_in;

     cuDoubleComplex    *vec_in_gpu,*vec_out_gpu,*vec_tmp_gpu;
     double             *quart_gpu;
     cuDoubleComplex    *vec_in_gpu_pointer,*vec_out_gpu_pointer,*vec_tmp_gpu_pointer;
     double             *quart_gpu_pointer;
     cuDoubleComplex    *Eb,*Ec,*Vbc,*Eb_p,*Ec_p,*Vbc_p;
     int                *sector_states_p,*sector_ranks_p,*imporbs_p,*bathorbs_p;
     int                *maskEb_p,*maskEc_p,*maskVbc_p,*sector_states,*sector_ranks,*imporbs,*bathorbs;
     int                *maskEb,*maskEc,*maskVbc;
     cuDoubleComplex    *GS_gpu,*GS_gpu_pointer;

     double diag[Niter_lanczos], subdiag[Niter_lanczos];

     if(verbose==1) printf(" GPU get eigenvalues, size of Lanczos vectors = %d ", size);

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

       cudaSetDevice(0); cudaSetDeviceFlags( cudaDeviceMapHost ); cudaEventSynchronize(0);
       if(use_texture_lcf==1){
         cudaMalloc((void**)&vec_in_gpu_pointer,sizeof(cuDoubleComplex)*size);
         cudaMalloc((void**)&Eb_p  , bathnorbs*bathnorbs*sizeof(cuDoubleComplex) );
         cudaMalloc((void**)&Ec_p  , impnorbs*impnorbs*sizeof(cuDoubleComplex) );
         cudaMalloc((void**)&Vbc_p , bathnorbs*impnorbs*sizeof(cuDoubleComplex) );
         cudaMemcpy(Eb_p,  Eb_in,  bathnorbs*bathnorbs*sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
         cudaMemcpy(Ec_p,  Ec_in,  impnorbs*impnorbs*sizeof(cuDoubleComplex),   cudaMemcpyHostToDevice);
         cudaMemcpy(Vbc_p, Vbc_in, bathnorbs*impnorbs*sizeof(cuDoubleComplex),  cudaMemcpyHostToDevice);
         bind_x(vec_in_gpu_pointer,size);
         bind_x_(Eb_p,Ec_p,Vbc_p);
       }else{
         cudaHostAlloc((void**)&vec_in_gpu    , size*sizeof(cuDoubleComplex),        cudaHostAllocMapped | cudaHostAllocPortable );
         cudaHostGetDevicePointer((void**)      &vec_in_gpu_pointer  ,  vec_in_gpu    , 0 );
         cudaHostAlloc((void**)&Eb            , bathnorbs*bathnorbs*sizeof(cuDoubleComplex), cudaHostAllocMapped | cudaHostAllocPortable );
         cudaHostAlloc((void**)&Ec            , impnorbs*impnorbs*sizeof(cuDoubleComplex),   cudaHostAllocMapped | cudaHostAllocPortable );
         cudaHostAlloc((void**)&Vbc           , bathnorbs*impnorbs*sizeof(cuDoubleComplex),  cudaHostAllocMapped | cudaHostAllocPortable );
         cudaMemcpy(Eb,  Eb_in,  bathnorbs*bathnorbs*sizeof(cuDoubleComplex), cudaMemcpyHostToHost);
         cudaMemcpy(Ec,  Ec_in,  impnorbs*impnorbs*sizeof(cuDoubleComplex),   cudaMemcpyHostToHost);
         cudaMemcpy(Vbc, Vbc_in, bathnorbs*impnorbs*sizeof(cuDoubleComplex),  cudaMemcpyHostToHost);
         cudaHostGetDevicePointer((void**)  &Eb_p   , Eb  , 0 );
         cudaHostGetDevicePointer((void**)  &Ec_p   , Ec  , 0 );
         cudaHostGetDevicePointer((void**)  &Vbc_p  , Vbc , 0 );
       }

       if(use_texture_rank_c==0){
         cudaHostAlloc((void**)&sector_ranks  , dimen*sizeof(int), cudaHostAllocMapped | cudaHostAllocPortable );
         cudaMemcpy(sector_ranks , sector_ranks_in , dimen*sizeof(int), cudaMemcpyHostToHost);
         cudaHostGetDevicePointer((void**)   &sector_ranks_p      , sector_ranks   , 0 );
       }else{
         cudaMalloc((void**)&sector_ranks_p  , dimen*sizeof(int) );
         cudaMemcpy(         sector_ranks_p , sector_ranks_in , dimen*sizeof(int), cudaMemcpyHostToDevice);
         bind_x_r(sector_ranks_p);
       };

       cudaHostAlloc((void**)&vec_tmp_gpu   , size*sizeof(cuDoubleComplex),     cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&vec_out_gpu   , size*sizeof(cuDoubleComplex),     cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&quart_gpu     , size*sizeof(double),              cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&sector_states , size*sizeof(int),                 cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&imporbs       , impnorbs*sizeof(int),             cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&bathorbs      , bathnorbs*sizeof(int),            cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&maskEb        , bathnorbs*bathnorbs*sizeof(int),  cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&maskEc        , impnorbs*impnorbs*sizeof(int),    cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&maskVbc       , bathnorbs*impnorbs*sizeof(int),   cudaHostAllocMapped | cudaHostAllocPortable );
       cudaHostAlloc((void**)&GS_gpu        , size*sizeof(cuDoubleComplex),     cudaHostAllocMapped | cudaHostAllocPortable );

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

   memset ((void **)GS_gpu, 0, size*sizeof(cuDoubleComplex));

    if(use_texture_lcf==0){
     for(int i=0;i<size;i++) vec_in_gpu[i]=make_cuDoubleComplex(1.0,0.0); cudaThreadSynchronize();
    }else{
     for(int i=0;i<size;i++) vec_out_gpu[i]=make_cuDoubleComplex(1.0,0.0); cudaThreadSynchronize();
     cudaMemcpy(vec_in_gpu_pointer,vec_out_gpu,size*sizeof(cuDoubleComplex),cudaMemcpyHostToDevice);
    }

    double *normv_ker; double normv; double *normv_loc; cudaMalloc((void**)&normv_ker,sizeof(double)); cuDoubleComplex coef;

 if(use_texture_lcf==0){
   coef = make_cuDoubleComplex(1.0/sqrt(cuCabs(cublasZdotu_no_device(size,vec_in_gpu,1,vec_in_gpu,1)))*vecp[0],0); cudaThreadSynchronize();  cudaEventSynchronize(0);
   cublasZaxpy_no_device(size,coef,vec_in_gpu,1,GS_gpu,1); cudaThreadSynchronize(); cudaEventSynchronize(0);cudaThreadSynchronize();
 }else{
   norm_vec_ker_complex<<<1,512>>>(size,vec_in_gpu_pointer,normv_ker); cudaEventSynchronize(0); cudaThreadSynchronize();
   normv_loc=&normv; cudaMemcpy(normv_loc,normv_ker,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(vec_out_gpu,vec_in_gpu,size*sizeof(cuDoubleComplex),cudaMemcpyDeviceToHost);
   coef=make_cuDoubleComplex(vecp[0]/normv,0.0);
   cublasZaxpy_no_device(size,coef,vec_out_gpu,1,GS_gpu,1); cudaThreadSynchronize(); cudaEventSynchronize(0);cudaThreadSynchronize();
 };

 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

  if(verbose==1) printf( " first step ... \n ");

    one_step_lanczos_complex_fly_cuda(blocksize,Niter_lanczos,size,0,vec_tmp_gpu,vec_in_gpu,
                    vec_out_gpu,vec_tmp_gpu_pointer,vec_in_gpu_pointer,vec_out_gpu_pointer,
                    quart_gpu_pointer,diag,subdiag,Eb_p,Ec_p,Vbc_p,sector_states_p,sector_ranks_p,norbs,
                    bathnorbs,impnorbs,imporbs_p,bathorbs_p,maskEb_p,maskEc_p,maskVbc_p);

  for(int iter=1;iter<Niter_lanczos-1;iter++){

   coef = make_cuDoubleComplex(1.0/sqrt(cuCabs(cublasZdotu_no_device(size,vec_out_gpu,1,vec_out_gpu,1))),0.0);
   cublasZscal_no_device (size,coef, vec_out_gpu,1);  cudaEventSynchronize(0);

   coef = make_cuDoubleComplex(vecp[iter],0);
   cublasZaxpy_no_device(size,coef,vec_out_gpu,1,GS_gpu,1); cudaThreadSynchronize(); cudaEventSynchronize(0);

   one_step_lanczos_complex_fly_cuda(blocksize,Niter_lanczos,size,iter,vec_tmp_gpu,vec_in_gpu,
                    vec_out_gpu,vec_tmp_gpu_pointer,vec_in_gpu_pointer,vec_out_gpu_pointer,
                    quart_gpu_pointer,diag,subdiag,Eb_p,Ec_p,Vbc_p,sector_states_p,sector_ranks_p,norbs,
                    bathnorbs,impnorbs,imporbs_p,bathorbs_p,maskEb_p,maskEc_p,maskVbc_p);

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
 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//
 //---------------------------------------------------------------------------------------//

   cudaFree(normv_ker);
 
   cudaEventSynchronize(0); 
   if(use_texture_lcf==0)
   {
    cudaFreeHost(vec_in_gpu);cudaFreeHost(Eb);cudaFreeHost(Ec);cudaFreeHost(Vbc);
   }else{
    unbind_x(); unbind_x_();
    cudaFree(vec_in_gpu_pointer); cudaFree(Eb_p);cudaFree(Ec_p);cudaFree(Vbc_p);
   }
   if(use_texture_rank_c==0){ cudaFreeHost(sector_ranks);}else{unbind_x_r(); cudaFree(sector_ranks_p);};

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
