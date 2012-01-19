 
 #include <stdio.h>
 #include <stdlib.h>
 #include <cuda_runtime.h>
 #include <cublas.h>
 #include <cuComplex.h>


#define MAX_BLOCK 22

cuDoubleComplex *Eb , *frequ , *collect_array;


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

  __global__ void build_array_inverse_array (int nthreads, int nfrequ, int nnn, double *collect )

{

    double sum,x;
    double Eb_[32*32];

    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////

    for(int ifrequ=threadIdx.x;ifrequ<nfrequ;ifrequ+=nthreads) {

    for(int iorb=0;iorb<nnn*nnn;iorb++)   Eb_[iorb]  = collect[iorb+nnn*nnn*ifrequ];

    for (int i=1; i < nnn; i++) Eb_[i] /= Eb_[0]; 

    for (int i=1; i < nnn; i++)  
      {
      for (int j=i; j < nnn; j++)  { 
        sum = 0.0;
        for (int k = 0; k < i; k++)
            sum += Eb_[j*nnn+k] * Eb_[k*nnn+i];
        Eb_[j*nnn+i] -= sum;
        }
      if (i == nnn-1) continue;

    for (int j=i+1; j < nnn; j++)  
       {  
        sum = 0.0;
        for (int k = 0; k < i; k++) sum += Eb_[i*nnn+k]*Eb_[k*nnn+j];
        Eb_[i*nnn+j] = (Eb_[i*nnn+j]-sum) / Eb_[i*nnn+i];
        }
      }

    for ( int i = 0; i < nnn; i++ )  
      for ( int j = i; j < nnn; j++ )  {
        x = 1.0;
        if ( i != j ) {
          x = 0.0;
          for ( int k = i; k < j; k++ )
              x -= Eb_[j*nnn+k]*Eb_[k*nnn+i];
          }
        Eb_[j*nnn+i] = x / Eb_[j*nnn+j];
        }

    for ( int i = 0; i < nnn; i++ )  
      for ( int j = i; j < nnn; j++ )  {
        if ( i == j ) continue;
        sum = 0.0;
        for ( int k = i; k < j; k++ )
            sum += Eb_[k*nnn+j]*( (i==k) ? 1.0 : Eb_[i*nnn+k] );
        Eb_[i*nnn+j] = -sum;
        }

    for ( int i = 0; i < nnn; i++ )   
      for ( int j = 0; j < nnn; j++ )  {
        sum = 0.0;
        for ( int k = ((i>j)?i:j); k < nnn; k++ )
            sum += ((j==k)?1.0:Eb_[j*nnn+k])*Eb_[k*nnn+i];
        Eb_[j*nnn+i] = sum;
        }

    for(int iorb=0;iorb<nnn*nnn;iorb++) collect[iorb+ifrequ*nnn*nnn]=Eb_[iorb];

   };

    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////

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

  __global__ void build_array_inverse_collect (int nthreads, int nfrequ, int nnn, double* Eb,double *collect, double *vec )

{

    double sum,x;

    __shared__ double shared_[32*32];
               double Eb_[32*32];
    
    for(int iorb=0;iorb<nnn*nnn;iorb++) shared_[iorb]=Eb[iorb];

    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////

    for(int ifrequ=threadIdx.x;ifrequ<nfrequ;ifrequ+=nthreads) {

    for(int iorb=0;iorb<nnn*nnn;iorb++)   Eb_[iorb         ]  = -shared_[iorb];
    for(int iorb=0;iorb<nnn    ;iorb++)   Eb_[iorb*nnn+iorb] +=  vec[ifrequ];

    for (int i=1; i < nnn; i++) Eb_[i] /= Eb_[0]; 

    for (int i=1; i < nnn; i++)  
      {
      for (int j=i; j < nnn; j++)  { 
        sum = 0.0;
        for (int k = 0; k < i; k++)
            sum += Eb_[j*nnn+k] * Eb_[k*nnn+i];
        Eb_[j*nnn+i] -= sum;
        }
      if (i == nnn-1) continue;

    for (int j=i+1; j < nnn; j++)  
       {  
        sum = 0.0;
        for (int k = 0; k < i; k++) sum += Eb_[i*nnn+k]*Eb_[k*nnn+j];
        Eb_[i*nnn+j] = (Eb_[i*nnn+j]-sum) / Eb_[i*nnn+i];
        }
      }

    for ( int i = 0; i < nnn; i++ )  
      for ( int j = i; j < nnn; j++ )  {
        x = 1.0;
        if ( i != j ) {
          x = 0.0;
          for ( int k = i; k < j; k++ )
              x -= Eb_[j*nnn+k]*Eb_[k*nnn+i];
          }
        Eb_[j*nnn+i] = x / Eb_[j*nnn+j];
        }

    for ( int i = 0; i < nnn; i++ )  
      for ( int j = i; j < nnn; j++ )  {
        if ( i == j ) continue;
        sum = 0.0;
        for ( int k = i; k < j; k++ )
            sum += Eb_[k*nnn+j]*( (i==k) ? 1.0 : Eb_[i*nnn+k] );
        Eb_[i*nnn+j] = -sum;
        }

    for ( int i = 0; i < nnn; i++ )   
      for ( int j = 0; j < nnn; j++ )  {
        sum = 0.0;
        for ( int k = ((i>j)?i:j); k < nnn; k++ )
            sum += ((j==k)?1.0:Eb_[j*nnn+k])*Eb_[k*nnn+i];
        Eb_[j*nnn+i] = sum;
        }

    for(int iorb=0;iorb<nnn*nnn;iorb++) collect[iorb+ifrequ*nnn*nnn]=Eb_[iorb];

   };

    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////

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

  __global__ void build_array_inverse (int nthreads, int nfrequ, int nnn, double* Eb, double *vec )

{

    double sum,x;

    __shared__ double shared_[32*32];
               double Eb_[32*32],tot[32*32];
    
    for(int iorb=0;iorb<nnn*nnn;iorb++) shared_[iorb]=Eb[iorb];

    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////

    for(int iorb=0;iorb<nnn*nnn;iorb++) tot[iorb]=0.0;

    for(int ifrequ=threadIdx.x;ifrequ<nfrequ;ifrequ+=nthreads) {

    for(int iorb=0;iorb<nnn*nnn;iorb++)   Eb_[iorb         ]  = -shared_[iorb];
    for(int iorb=0;iorb<nnn    ;iorb++)   Eb_[iorb*nnn+iorb] +=  vec[ifrequ];

    for (int i=1; i < nnn; i++) Eb_[i] /= Eb_[0]; 

    for (int i=1; i < nnn; i++)  
      {
      for (int j=i; j < nnn; j++)  { 
        sum = 0.0;
        for (int k = 0; k < i; k++)
            sum += Eb_[j*nnn+k] * Eb_[k*nnn+i];
        Eb_[j*nnn+i] -= sum;
        }
      if (i == nnn-1) continue;

    for (int j=i+1; j < nnn; j++)  
       {  
        sum = 0.0;
        for (int k = 0; k < i; k++) sum += Eb_[i*nnn+k]*Eb_[k*nnn+j];
        Eb_[i*nnn+j] = (Eb_[i*nnn+j]-sum) / Eb_[i*nnn+i];
        }
      }

    for ( int i = 0; i < nnn; i++ )  
      for ( int j = i; j < nnn; j++ )  {
        x = 1.0;
        if ( i != j ) {
          x = 0.0;
          for ( int k = i; k < j; k++ )
              x -= Eb_[j*nnn+k]*Eb_[k*nnn+i];
          }
        Eb_[j*nnn+i] = x / Eb_[j*nnn+j];
        }

    for ( int i = 0; i < nnn; i++ )  
      for ( int j = i; j < nnn; j++ )  {
        if ( i == j ) continue;
        sum = 0.0;
        for ( int k = i; k < j; k++ )
            sum += Eb_[k*nnn+j]*( (i==k) ? 1.0 : Eb_[i*nnn+k] );
        Eb_[i*nnn+j] = -sum;
        }

    for ( int i = 0; i < nnn; i++ )   
      for ( int j = 0; j < nnn; j++ )  {
        sum = 0.0;
        for ( int k = ((i>j)?i:j); k < nnn; k++ )
            sum += ((j==k)?1.0:Eb_[j*nnn+k])*Eb_[k*nnn+i];
        Eb_[j*nnn+i] = sum;
        }

    for(int iorb=0;iorb<nnn*nnn;iorb++) tot[iorb]+=Eb_[iorb];

   };

    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////

   for(int iorb=0;iorb<nnn*nnn;iorb++) Eb_[iorb]=0.0;
   for(int ii=0;ii<nthreads;ii++) 
   {
    if(threadIdx.x==ii){ for(int iorb=0;iorb<nnn*nnn;iorb++) shared_[iorb] = tot[iorb];};
                         __syncthreads();
                         for(int iorb=0;iorb<nnn*nnn;iorb++) Eb_[iorb]+= shared_[iorb];  
   };
                         __syncthreads();
   for(int iorb=0;iorb<nnn*nnn;iorb++) Eb[iorb] = Eb_[iorb];

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

extern "C" void sum_of_inverse_frequ_(int* nnn_, int* nfrequ_,  double* Eb_, double* totsum_ , double *frequ_)

{

  int nnn ; int nfrequ; nnn=*nnn_; 
  nfrequ=*nfrequ_;

  double *Eb , *frequ ;

  cudaMalloc (  (void**)&Eb     , nnn*nnn   *sizeof(double) );
  cudaMalloc (  (void**)&frequ  , nfrequ *sizeof(double) );
  cudaMemcpy ( Eb      ,Eb_     , nnn*nnn*sizeof(double)   , cudaMemcpyHostToDevice);
  cudaMemcpy ( frequ   ,frequ_  , nfrequ*sizeof(double) , cudaMemcpyHostToDevice);

  if(nnn>32) { printf( " sum of inverse real cuda (1) : matrices are too big!!!!! \n");};

  build_array_inverse <<<1,512>>> ( 512, nfrequ, nnn, Eb, frequ );

  cudaEventSynchronize(0); cudaThreadSynchronize();

  cudaMemcpy(totsum_ , Eb , nnn*nnn*sizeof(double) , cudaMemcpyDeviceToHost);
  cudaFree(frequ); cudaFree(Eb);

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

extern "C" void sum_of_inverse_frequ_collect_(int* nnn_, int* nfrequ_,  double* Eb_, double* totsum_ , double *frequ_)

{

  int nnn ; int nfrequ; nnn=*nnn_; 
  nfrequ=*nfrequ_;

  double *Eb , *frequ , *collect;

  cudaMalloc (  (void**)&Eb     , nnn*nnn*sizeof(double) );
  cudaMalloc (  (void**)&frequ  , nfrequ *sizeof(double) );
  cudaMalloc (  (void**)&collect, nnn*nnn*nfrequ *sizeof(double) );

  cudaMemcpy ( Eb      ,Eb_     , nnn*nnn*sizeof(double)   , cudaMemcpyHostToDevice);
  cudaMemcpy ( frequ   ,frequ_  , nfrequ*sizeof(double) , cudaMemcpyHostToDevice);

  if(nnn>32) { printf( " sum of inverse real cuda (2) : matrices are too big!!!!! \n" );};

  build_array_inverse_collect <<<1,512>>> ( 512, nfrequ, nnn, Eb, collect, frequ );

  cudaEventSynchronize(0); cudaThreadSynchronize();

  cudaMemcpy(totsum_ , collect , nnn*nnn*nfrequ*sizeof(double) , cudaMemcpyDeviceToHost);
  cudaFree(frequ); cudaFree(Eb);

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

extern "C" void sum_of_inverse_frequ_array_(int* nnn_, int* nfrequ_, double* totsum_)

{

  int nnn ; int nfrequ; nnn=*nnn_; 
   nfrequ=*nfrequ_;

  double *collect;

  cudaMalloc (  (void**)&collect, nnn*nnn*nfrequ *sizeof(double) );
  cudaMemcpy ( collect      ,totsum_     , nnn*nnn*nfrequ*sizeof(double)   , cudaMemcpyHostToDevice);

  if(nnn>32) { printf( " sum of inverse real cuda (3) : matrices are too big!!!!! \n" );};

  build_array_inverse_array <<<1,512>>> ( 512, nfrequ, nnn, collect );

  cudaEventSynchronize(0); cudaThreadSynchronize();

  cudaMemcpy(totsum_ , collect , nnn*nnn*nfrequ*sizeof(double) , cudaMemcpyDeviceToHost);
  cudaFree(collect); 

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

  __global__ void build_array_inverse_complex (int blocksize, int nfrequ, int nnn, cuDoubleComplex* Eb, cuDoubleComplex *vec )

{
    cuDoubleComplex sum,x,moins1,zero,uno;
    zero=make_cuDoubleComplex(0.0,0.0); moins1=make_cuDoubleComplex(-1.0,0.0); uno=make_cuDoubleComplex(1.0,0.0);

    __shared__ cuDoubleComplex sum_[512];
    __shared__ cuDoubleComplex shared_[MAX_BLOCK*MAX_BLOCK]; 
               cuDoubleComplex Eb_[MAX_BLOCK*MAX_BLOCK]; 
               cuDoubleComplex tot_[MAX_BLOCK*MAX_BLOCK];

    for(int iorb=0;iorb<nnn*nnn;iorb++) shared_[iorb]=cuCmul(Eb[iorb],moins1);
    for(int iorb=0;iorb<nnn*nnn;iorb++) tot_[iorb]=zero;

    for(int ifrequ=threadIdx.x;ifrequ<nfrequ;ifrequ+=blocksize) {

    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////

    for(int iorb=0;iorb<nnn*nnn;iorb++)   Eb_[iorb         ]      = shared_[iorb];
    for(int iorb=0;iorb<nnn    ;iorb++)   Eb_[iorb*nnn+iorb]      = cuCadd(Eb_[iorb*nnn+iorb],vec[ifrequ]);

    for (int i=1; i < nnn; i++) Eb_[i] = cuCdiv(Eb_[i],Eb_[0]); 

    for (int i=1; i < nnn; i++)  
      {
      for (int j=i; j < nnn; j++)  { 
        sum = zero;
        for (int k = 0; k < i; k++)
            sum = cuCadd(sum, cuCmul(Eb_[j*nnn+k] , Eb_[k*nnn+i]));
        Eb_[j*nnn+i] = cuCsub(Eb_[j*nnn+i] ,sum);
      }
      if (i == nnn-1) continue;

    for (int j=i+1; j < nnn; j++)  
       {  
        sum = zero;
        for (int k = 0; k < i; k++) sum = cuCadd(sum,cuCmul(Eb_[i*nnn+k],Eb_[k*nnn+j]));
        Eb_[i*nnn+j] = cuCdiv( cuCsub( Eb_[i*nnn+j],sum),Eb_[i*nnn+i]);
        }
      }

    for ( int i = 0; i < nnn; i++ )  
      for ( int j = i; j < nnn; j++ )  {
        x = uno ;
        if ( i != j ) {
          x = zero;
          for ( int k = i; k < j; k++ )
              x = cuCsub(x, cuCmul(Eb_[j*nnn+k],Eb_[k*nnn+i]));
          }
        Eb_[j*nnn+i] = cuCdiv(x , Eb_[j*nnn+j]);
        }

    for ( int i = 0; i < nnn; i++ )  
      for ( int j = i; j < nnn; j++ )  {
        if ( i == j ) continue;
        sum = zero;
        for ( int k = i; k < j; k++ )
           if(i==k){
            sum = cuCadd(sum,Eb_[k*nnn+j]);
           }else{
            sum = cuCadd(sum,cuCmul(Eb_[k*nnn+j],Eb_[i*nnn+k] ));
           };
        Eb_[i*nnn+j] = cuCmul(sum,moins1);
        }

    for ( int i = 0; i < nnn; i++ )   
      for ( int j = 0; j < nnn; j++ )  {
        sum = zero;
        for ( int k = ((i>j)?i:j); k < nnn; k++ )
          if(j==k){
            sum = cuCadd(sum,Eb_[k*nnn+i]);
          }else{
            sum = cuCadd(sum,cuCmul(Eb_[j*nnn+k],Eb_[k*nnn+i]));
          }
        Eb_[j*nnn+i] = sum;
        }

      for(int iorb=0;iorb<nnn*nnn;iorb++) 
      { tot_[iorb]=cuCadd(tot_[iorb],Eb_[iorb]); };

   };

    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////
  
   // parallel reduction
   for(int iorb=0;iorb<nnn*nnn;iorb++) 
    {
     int i =threadIdx.x;
     sum_[i]=tot_[iorb]; 
     __syncthreads();
     for(int bit=blocksize/2; bit>0; bit/=2)
      {
      cuDoubleComplex t = cuCadd(sum_[i],sum_[i^bit]); __syncthreads();
      sum_[i]=t ;__syncthreads();
      }
      Eb[iorb]=sum_[i];
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

extern "C" void sum_of_inverse_frequ_complex_(int* nnn_, int* nfrequ_,  cuDoubleComplex* Eb_, cuDoubleComplex* totsum_ , 
                                             cuDoubleComplex *frequ_, int* firstlast)

{

  int nnn ; int nfrequ; nnn=*nnn_; 
  nfrequ=*nfrequ_;

  const int nthreads = 512 ;

  if(nnn>MAX_BLOCK)  { printf( " sum of inverse complex cuda (1) : matrices are too big!!!!! \n"); };

if(*firstlast==1){
  cudaMalloc (  (void**)&Eb     , nnn*nnn   *sizeof(cuDoubleComplex) );
  cudaMalloc (  (void**)&frequ  , nfrequ *sizeof(cuDoubleComplex) );
};

  cudaMemcpy ( Eb      ,Eb_     , nnn*nnn*sizeof(cuDoubleComplex)   , cudaMemcpyHostToDevice);
if(*firstlast==1){
  cudaMemcpy ( frequ   ,frequ_  , nfrequ*sizeof(cuDoubleComplex) , cudaMemcpyHostToDevice);
}

  build_array_inverse_complex <<<1,nthreads>>> ( nthreads, nfrequ, nnn, Eb, frequ);
  cudaEventSynchronize(0); cudaThreadSynchronize();

  cudaMemcpy(totsum_ , Eb , nnn*nnn*sizeof(cuDoubleComplex) , cudaMemcpyDeviceToHost);
  if(*firstlast==2){ cudaFree(frequ); cudaFree(Eb);  };

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
//********************************************
//********************************************

  __global__ void build_array_inverse_complex_collect (int blocksize, int nfrequ, int nnn, cuDoubleComplex* Eb, 
                                 cuDoubleComplex *tot,cuDoubleComplex *vec )

{
    cuDoubleComplex sum,x,moins1,zero,uno;
    zero=make_cuDoubleComplex(0.0,0.0); moins1=make_cuDoubleComplex(-1.0,0.0); uno=make_cuDoubleComplex(1.0,0.0);

    __shared__ cuDoubleComplex shared_[MAX_BLOCK*MAX_BLOCK]; 
               cuDoubleComplex Eb_[MAX_BLOCK*MAX_BLOCK]; 

    for(int iorb=0;iorb<nnn*nnn;iorb++) shared_[iorb]=cuCmul(Eb[iorb],moins1);

    for(int ifrequ=threadIdx.x;ifrequ<nfrequ;ifrequ+=blocksize) {

    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////

    for(int iorb=0;iorb<nnn*nnn;iorb++)   Eb_[iorb         ]      = shared_[iorb];
    for(int iorb=0;iorb<nnn    ;iorb++)   Eb_[iorb*nnn+iorb]      = cuCadd(Eb_[iorb*nnn+iorb],vec[ifrequ]);

    for (int i=1; i < nnn; i++) Eb_[i] = cuCdiv(Eb_[i],Eb_[0]); 

    for (int i=1; i < nnn; i++)  
      {
      for (int j=i; j < nnn; j++)  { 
        sum = zero;
        for (int k = 0; k < i; k++)
            sum = cuCadd(sum, cuCmul(Eb_[j*nnn+k] , Eb_[k*nnn+i]));
        Eb_[j*nnn+i] = cuCsub(Eb_[j*nnn+i] ,sum);
      }
      if (i == nnn-1) continue;

    for (int j=i+1; j < nnn; j++)  
       {  
        sum = zero;
        for (int k = 0; k < i; k++) sum = cuCadd(sum,cuCmul(Eb_[i*nnn+k],Eb_[k*nnn+j]));
        Eb_[i*nnn+j] = cuCdiv( cuCsub( Eb_[i*nnn+j],sum),Eb_[i*nnn+i]);
        }
      }

    for ( int i = 0; i < nnn; i++ )  
      for ( int j = i; j < nnn; j++ )  {
        x = uno ;
        if ( i != j ) {
          x = zero;
          for ( int k = i; k < j; k++ )
              x = cuCsub(x, cuCmul(Eb_[j*nnn+k],Eb_[k*nnn+i]));
          }
        Eb_[j*nnn+i] = cuCdiv(x , Eb_[j*nnn+j]);
        }

    for ( int i = 0; i < nnn; i++ )  
      for ( int j = i; j < nnn; j++ )  {
        if ( i == j ) continue;
        sum = zero;
        for ( int k = i; k < j; k++ )
           if(i==k){
            sum = cuCadd(sum,Eb_[k*nnn+j]);
           }else{
            sum = cuCadd(sum,cuCmul(Eb_[k*nnn+j],Eb_[i*nnn+k] ));
           };
        Eb_[i*nnn+j] = cuCmul(sum,moins1);
        }

    for ( int i = 0; i < nnn; i++ )   
      for ( int j = 0; j < nnn; j++ )  {
        sum = zero;
        for ( int k = ((i>j)?i:j); k < nnn; k++ )
          if(j==k){
            sum = cuCadd(sum,Eb_[k*nnn+i]);
          }else{
            sum = cuCadd(sum,cuCmul(Eb_[j*nnn+k],Eb_[k*nnn+i]));
          }
        Eb_[j*nnn+i] = sum;
        }

      for(int iorb=0;iorb<nnn*nnn;iorb++) { tot[iorb+ifrequ*nnn*nnn]=Eb_[iorb];};

   };

    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////
  
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

  __global__ void build_array_inverse_complex_array (int blocksize, int nfrequ, int nnn, cuDoubleComplex *tot )

{
    cuDoubleComplex sum,x,moins1,zero,uno;
    zero=make_cuDoubleComplex(0.0,0.0); moins1=make_cuDoubleComplex(-1.0,0.0); uno=make_cuDoubleComplex(1.0,0.0);
    cuDoubleComplex Eb_[MAX_BLOCK*MAX_BLOCK]; 
    int ifrequ=blocksize*blockIdx.x+threadIdx.x;

    if(ifrequ<nfrequ) {

    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////

    for(int iorb=0;iorb<nnn*nnn;iorb++) { Eb_[iorb] = tot[iorb+nnn*nnn*ifrequ];}

    for (int i=1; i < nnn; i++) Eb_[i] = cuCdiv(Eb_[i],Eb_[0]); 

    for (int i=1; i < nnn; i++)  
      {
      for (int j=i; j < nnn; j++)  { 
        sum = zero;
        for (int k = 0; k < i; k++)
            sum = cuCadd(sum, cuCmul(Eb_[j*nnn+k] , Eb_[k*nnn+i]));
        Eb_[j*nnn+i] = cuCsub(Eb_[j*nnn+i] ,sum);
      }
      if (i == nnn-1) continue;

    for (int j=i+1; j < nnn; j++)  
       {  
        sum = zero;
        for (int k = 0; k < i; k++) sum = cuCadd(sum,cuCmul(Eb_[i*nnn+k],Eb_[k*nnn+j]));
        Eb_[i*nnn+j] = cuCdiv( cuCsub( Eb_[i*nnn+j],sum),Eb_[i*nnn+i]);
        }
      }

    for ( int i = 0; i < nnn; i++ )  
      for ( int j = i; j < nnn; j++ )  {
        x = uno ;
        if ( i != j ) {
          x = zero;
          for ( int k = i; k < j; k++ )
              x = cuCsub(x, cuCmul(Eb_[j*nnn+k],Eb_[k*nnn+i]));
          }
        Eb_[j*nnn+i] = cuCdiv(x , Eb_[j*nnn+j]);
        }

    for ( int i = 0; i < nnn; i++ )  
      for ( int j = i; j < nnn; j++ )  {
        if ( i == j ) continue;
        sum = zero;
        for ( int k = i; k < j; k++ )
           if(i==k){
            sum = cuCadd(sum,Eb_[k*nnn+j]);
           }else{
            sum = cuCadd(sum,cuCmul(Eb_[k*nnn+j],Eb_[i*nnn+k] ));
           };
        Eb_[i*nnn+j] = cuCmul(sum,moins1);
        }

    for ( int i = 0; i < nnn; i++ )   
      for ( int j = 0; j < nnn; j++ )  {
        sum = zero;
        for ( int k = ((i>j)?i:j); k < nnn; k++ )
          if(j==k){
            sum = cuCadd(sum,Eb_[k*nnn+i]);
          }else{
            sum = cuCadd(sum,cuCmul(Eb_[j*nnn+k],Eb_[k*nnn+i]));
          }
        Eb_[j*nnn+i] = sum;
        }

      for(int iorb=0;iorb<nnn*nnn;iorb++) { tot[iorb+ifrequ*nnn*nnn]=Eb_[iorb];};

   };

    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////
  
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
//********************************************

extern "C" void sum_of_inverse_frequ_complex_collect_(int* nnn_, int* nfrequ_,  cuDoubleComplex* Eb_, 
                                             cuDoubleComplex* collect_ , 
                                             cuDoubleComplex *frequ_, int* firstlast)

{

  int nnn ; int nfrequ; nnn=*nnn_; 
  nfrequ=*nfrequ_;

  const int nthreads = 512 ;

  if(nnn>MAX_BLOCK)           { printf( " sum of inverse complex cuda (2) : matrices are too big!!!!! \n"); };

if(*firstlast==1){
  cudaMalloc (  (void**)&Eb     , nnn*nnn   *sizeof(cuDoubleComplex) );
  cudaMalloc (  (void**)&frequ  , nfrequ *sizeof(cuDoubleComplex) );
  cudaMalloc (  (void**)&collect_array  , nfrequ*nnn*nnn*sizeof(cuDoubleComplex) );
};

  cudaMemcpy ( Eb      ,Eb_     , nnn*nnn*sizeof(cuDoubleComplex)   , cudaMemcpyHostToDevice);
if(*firstlast==1){
  cudaMemcpy ( frequ   ,frequ_  , nfrequ*sizeof(cuDoubleComplex) , cudaMemcpyHostToDevice);
}

  build_array_inverse_complex_collect <<<1,nthreads>>> ( nthreads,nfrequ,nnn,Eb,collect_array,frequ);
  cudaEventSynchronize(0); cudaThreadSynchronize();
  cudaMemcpy (collect_ ,collect_array  , nfrequ*nnn*nnn*sizeof(cuDoubleComplex) , cudaMemcpyDeviceToHost);
  cudaEventSynchronize(0); cudaThreadSynchronize();

  if(*firstlast==2){ cudaFree(frequ); cudaFree(Eb); cudaFree(collect_array);};

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

extern "C" void sum_of_inverse_frequ_complex_array_(int* nnn_, int* nfrequ_, cuDoubleComplex* collect_ , int* firstlast)

{

  int nnn ; int nfrequ; nnn=*nnn_; nfrequ=*nfrequ_;
  const int nthreads = 256 ;

  if(nnn>MAX_BLOCK)           { printf( " sum of inverse complex cuda (3) : matrices are too big!!!!! \n"); };
  if(nfrequ/nthreads+1>65500) { printf( " too many blocks in inverse array cuda \n ");};

if(*firstlast==1){
  cudaMalloc (  (void**)&collect_array  , nfrequ*nnn*nnn*sizeof(cuDoubleComplex) );
};
  cudaMemcpy ( collect_array ,collect_  , nnn*nnn*nfrequ*sizeof(cuDoubleComplex) , cudaMemcpyHostToDevice);

  build_array_inverse_complex_array <<<nfrequ/nthreads+1,nthreads>>> ( nthreads,nfrequ,nnn,collect_array);
  cudaEventSynchronize(0); cudaThreadSynchronize();
  
  cudaMemcpy (collect_ ,collect_array  , nfrequ*nnn*nnn*sizeof(cuDoubleComplex) , cudaMemcpyDeviceToHost);

  if(*firstlast==2){cudaFree(collect_array);};

}

//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
//********************************************
