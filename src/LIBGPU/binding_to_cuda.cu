#include </opt/cuda_2.3/src/fortran.c>
#include </opt/cuda_2.3/include/cufft.h>

extern "C" void cufftplan1d_(cufftHandle* plan, int* nx, cufftType* cucu, int* batch)
{ int i=cufftPlan1d(plan,*nx,*cucu,*batch); }

extern "C" void cufftexecc2c_(int* plan, cufftComplex** data, cufftComplex** data2, int* cucu, int* nx, int* batch)
{ int i=cufftExecC2C(*plan,*data,*data2,*cucu); }

extern "C" void cufftdestroy_(int* plan)
{ cufftdestroy_(plan); }

