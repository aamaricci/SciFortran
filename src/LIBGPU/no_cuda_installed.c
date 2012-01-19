#include <stdio.h>
#include <stdlib.h>

 ///////////
 // TOOLS //
 ///////////
 
extern "C" void cuda_complex_invert_(int *pBlockSize, double** A_, double** invA_, int *pndim){}
extern "C" void cuda_c_invert_(int *pBlockSize, double** A_, double** A_i, double** invA_, double** invA_i, int *pndim){}
extern "C" void matmul_gpu_fortran_(int* pBLOCKSIZE, double** A, double** B, double** C,  int* hA, int* wA, int* wB){}
extern "C" void matmul_gpu_cxx(int* pBLOCKSIZE, double** A, double** B, double** C,  int* hA, int* wA, int* wB){}
extern "C" void matmul_gpu_fortran_c_(int* pBLOCKSIZE, double** A, double** B, double** C,  int* hA, int* wA, int* wB){}
extern "C" void matmul_gpu_cxx_c(int* pBLOCKSIZE, double** A, double** B, double** C,  int* hA, int* wA, int* wB){}
extern "C" void Mul(int size, const double* A, const double* B, double* C){}
extern "C" void cuda_invert__(int *pBlockSize, double** A_, double** invA_, int *pndim){}
extern "C" void cufftplan1d_(int* plan, int* nx, int* cucu, int* batch) { }
extern "C" void cufftexecc2c_(int* plan, float** data, float** data2, int* cucu, int* nx, int* batch) { }
extern "C" void cufftdestroy_(int* plan) { }


 //////////////////  
 // LANCZOS REAL //
 //////////////////

 
extern "C" void hmult_sz_real_cuda_rout_(int *pblocksize, int *offdiasize, int *roff,
int* ntot, double* QUART, double* diagsz, double* vec_in, double* vec_out, int* noffsz, int* rankoffsz, double* offdiagsz ){};

extern "C" void lanczos_test_gpu_(int *blocksize, int *ntot){};

extern "C" void lanczos_real_dynamic_cuda_(int *pblocksize, int *Niter_lanczos_,int *offdiasize, int *roff, int *ntot,
double *QUART, double *diagsz, int *noffsz, int *rankoffsz, double *offdiagsz, double *diag, double *subdiag , double *vecinit){};

extern "C" void lanczos_real_cuda_(int *pblocksize,  int *Niter_lanczos_,int *offdiasize, int *roff, int *ntot,
double *QUART, double *diagsz, int *noffsz, int *rankoffsz, double *offdiagsz, double *diag, double *subdiag ){};

extern "C" void lanczos_real_get_gs_cuda_(int *pblocksize, int *Niter_lanczos_,int *offdiasize, int *roff,
int *ntot, double *QUART, double *diagsz, int *noffsz, int *rankoffsz, double *offdiagsz, double *vecp, double *GS){}

extern "C" void lanczos_real_fly_cuda_(int *dimen_in, int *pblocksize, int* norbs_in, int *Niter_lanczos_, int* ntot, double* quart,
double* diag, double* subdiag, double* Eb_in, double* Ec_in, double* Vbc_in, int* sector_states_in, int* sector_ranks_in
,int* bathnorbs_in,int* impnorbs_in,int* imporbs_in, int* bathorbs_in, int* maskEb_in, int* maskEc_in, int* maskVbc_in){}

extern "C" void lanczos_real_fly_dynamic_cuda_(int *dimen_in, int *pblocksize, int* norbs_in, int *Niter_lanczos_, int* ntot, double* quart,
double* diag, double* subdiag, double* Eb_in, double* Ec_in, double* Vbc_in, int* sector_states_in, int* sector_ranks_in
,int* bathnorbs_in,int* impnorbs_in,int* imporbs_in, int* bathorbs_in, int* maskEb_in, int* maskEc_in, int* maskVbc_in, double* vecinit){}

extern "C" void lanczos_real_fly_gs_cuda_(int *dimen_in, int *pblocksize, int* norbs_in, int *Niter_lanczos_, int* ntot, double* quart,
double* Eb_in, double* Ec_in, double* Vbc_in, int* sector_states_in, int* sector_ranks_in ,int* bathnorbs_in,int* impnorbs_in,
int* imporbs_in, int* bathorbs_in, int* maskEb_in, int* maskEc_in, int* maskVbc_in,double *vecp,double *GS){}

extern "C" void lanczos_real_updo_cuda_(int *norbs, int *pblocksize, int *Niter_lanczos_,int *offdiasizeup, int* offdiasizedn,
int *roffup, int *roffdn, int *ntot,  int *sizup, int *sizdn,
double *QUART, double *diagup, double *diagdn, int *noffup, int *noffdn, int *rankoffup, int *rankoffdn,
double *offdiagup, double *offdiagdn, double *diag, double *subdiag, int *UMASK, int *statesup, int *statesdn , int *iorbup, int *iorbdn){};

extern "C" void lanczos_real_updo_dynamic_cuda_(int *norbs, int *pblocksize, int *Niter_lanczos_,int *offdiasizeup, int* offdiasizedn,
int *roffup, int *roffdn, int *ntot,  int *sizup, int *sizdn,
double *QUART, double *diagup, double *diagdn, int *noffup, int *noffdn, int *rankoffup, int *rankoffdn,
double *offdiagup, double *offdiagdn, double *diag, double *subdiag, int *UMASK, int *statesup, int *statesdn ,
int *iorbup, int *iorbdn, double *vecinit){};

extern "C" void lanczos_real_updo_gs_cuda_(int *norbs, int *pblocksize, int *Niter_lanczos_,int *offdiasizeup, int* offdiasizedn,
   int *roffup, int *roffdn, int *ntot,  int *sizup, int *sizdn,
   double *QUART, double *diagup, double *diagdn, int *noffup, int *noffdn, int *rankoffup, int *rankoffdn,
   double *offdiagup, double *offdiagdn, double *diag, double *subdiag, int *UMASK, int *statesup, int *statesdn ,
   int *iorbup, int *iorbdn, double *vecp, double *GS){};




 /////////////////////
 // LANCZOS COMPLEX //
 /////////////////////


extern "C" void hmult_sz_complex_cuda_rout_(int *pblocksize, int *offdiasize, int *roff, int* ntot,
double* QUART, double* diagsz, double* vec_in, double* vec_out, int* noffsz, int* rankoffsz, double* offdiagsz ){};

extern "C" void lanczos_dynamic_cuda_complex_(int *pblocksize, int *Niter_lanczos_,int *offdiasize, int *roff,
int *ntot, double *QUART, double *diagsz, int *noffsz, int *rankoffsz, double *offdiagsz, double *diag,
double *subdiag , double *vecinit){};

extern "C" void lanczos_cuda_complex_(int *pblocksize,  int *Niter_lanczos_,int *offdiasize, int *roff, int *ntot,
double *QUART, double *diagsz, int *noffsz, int *rankoffsz, double *offdiagsz, double *diag, double *subdiag ){};

extern "C" void lanczos_get_gs_cuda_complex_(int *pblocksize, int *Niter_lanczos_,int *offdiasize, int *roff,
int *ntot, double *QUART, double *diagsz, int *noffsz, int *rankoffsz, double *offdiagsz, double *vecp, double *GS){};

extern "C" void lanczos_complex_fly_gs_cuda_(int *dimen_in, int *pblocksize, int* norbs_in, int *Niter_lanczos_, int* ntot, double* quart,
  double* Eb_in, double* Ec_in, double* Vbc_in, int* sector_states_in, int* sector_ranks_in
 ,int* bathnorbs_in,int* impnorbs_in,int* imporbs_in, int* bathorbs_in, int* maskEb_in, int* maskEc_in, int* maskVbc_in,double *vecp,double *GS){};

extern "C" void lanczos_complex_fly_dynamic_cuda_(int *dimen_in, int *pblocksize, int* norbs_in, int *Niter_lanczos_, int* ntot, 
double* quart, double* diag, double* subdiag, double* Eb_in, double* Ec_in, double* Vbc_in, 
int* sector_states_in, int* sector_ranks_in ,int* bathnorbs_in,int* impnorbs_in,int* imporbs_in, int* bathorbs_in, int* maskEb_in, 
int* maskEc_in, int* maskVbc_in, double* vecinit){};

extern "C" void lanczos_complex_fly_cuda_(int *dimen_in, int *pblocksize, int* norbs_in, int *Niter_lanczos_, int* ntot, double* quart,
double* diag, double* subdiag, double* Eb_in, double* Ec_in, double* Vbc_in, int* sector_states_in, 
int* sector_ranks_in ,int* bathnorbs_in,int* impnorbs_in,int* imporbs_in, int* bathorbs_in, int* maskEb_in, int* maskEc_in, 
int* maskVbc_in) {};

extern "C" void lanczos_complex_updo_gs_cuda_(int *norbs, int *pblocksize,
   int *Niter_lanczos_,int *offdiasizeup, int* offdiasizedn, int *roffup, int *roffdn, int *ntot,  int *sizup, int *sizdn,
   double *QUART, double *diagup, double *diagdn, int *noffup, int *noffdn, int *rankoffup, int *rankoffdn,
   double *offdiagup, double *offdiagdn,  int *UMASK, int *statesup, int *statesdn,
   int *iorbup, int *iorbdn, double *vecp, double *GS){}

extern "C" void lanczos_complex_updo_dynamic_cuda_(int *norbs, int *pblocksize,
   int *Niter_lanczos_,int *offdiasizeup, int* offdiasizedn, int *roffup, int *roffdn, int *ntot,  int *sizup, int *sizdn,
   double *QUART, double *diagup, double *diagdn, int *noffup, int *noffdn, int *rankoffup, int *rankoffdn,
   double *offdiagup, double *offdiagdn, double *diag, double *subdiag, int *UMASK, int *statesup, int *statesdn ,
   int *iorbup, int *iorbdn, double *vecinit){}

extern "C" void lanczos_complex_updo_cuda_(int *norbs, int *pblocksize, int *Niter_lanczos_,
   int *offdiasizeup, int* offdiasizedn, int *roffup, int *roffdn, int *ntot,  int *sizup, int *sizdn,
   double *QUART, double *diagup, double *diagdn, int *noffup, int *noffdn, int *rankoffup, int *rankoffdn,
   double *offdiagup, double *offdiagdn, double *diag, double *subdiag, int *UMASK, int *statesup, int *statesdn , int *iorbup, int *iorbdn){}


 /////////////////////
 // SUM OF INVERSE  //
 /////////////////////

extern "C" void sum_of_inverse_frequ_(int* nnn_, int* nfrequ_,  double* Eb_, double* totsum_ , double *frequ_){}
extern "C" void sum_of_inverse_frequ_collect_(int* nnn_, int* nfrequ_,  double* Eb_, double* totsum_ , double *frequ_){}
extern "C" void sum_of_inverse_frequ_array_(int* nnn_, int* nfrequ_, double* totsum_){}
extern "C" void sum_of_inverse_frequ_complex_(int* nnn_, int* nfrequ_,  double* Eb_, double* totsum_ , double *frequ_, int* firstlast){}
extern "C" void sum_of_inverse_frequ_complex_collect_(int* nnn_, int* nfrequ_,  double* Eb_, double* collect_ , double *frequ_, int* firstlast){}
extern "C" void sum_of_inverse_frequ_complex_array_(int* nnn_, int* nfrequ_, double* collect_ , int* firstlast){}



