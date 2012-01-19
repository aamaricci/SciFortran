module fortran_cuda

 use fortran_cuda_routines

!-----------------------------------------------------------------------------------!
!                                  --------------                                   !
!                                      WARNING                                      !
!                                  --------------                                   !
!  in this module, the conversion from Fortran (array storage in columns) to C      !
!  (array storage in lines) is taken into account, but can be tricky. For           !
!  matrices diagonalization, I use the fact :                                       !
!  (A-1)A=Id, A_T (A-1)_T =1, (A-1)_T=(A_T)-1                                       !
!  so C receives A_T, diagonalize it and get (A_T)-1 inside C, and fortran          !
!  gets ((A_T)-1)_T, which is A_-1                                                  ! 
!  in MATMUL_cuda the call is done to a C routine which takes into account          !
!  the fortran column memory format.                                                !
!-----------------------------------------------------------------------------------!


 INTERFACE cuda_array_of_inverse
   MODULE PROCEDURE cuda_array_of_inverse_r,cuda_array_of_inverse_c
 END INTERFACE

 INTERFACE cuda_array_of_inverse_collect
   MODULE PROCEDURE cuda_array_of_inverse_collect_c,cuda_array_of_inverse_collect_r
 END INTERFACE

 INTERFACE cuda_array_of_inverse_array
   MODULE PROCEDURE cuda_array_of_inverse_array_r,cuda_array_of_inverse_array_c
 END INTERFACE


  contains

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

 subroutine cuda_array_of_inverse_r(nnn,nfrequ,Eb,totsum,frequ )
 implicit none
 integer :: nnn,nfrequ
 real(8) :: Eb(nnn,nnn),totsum(nnn,nnn),frequ(nfrequ)

 interface 
  subroutine sum_of_inverse_frequ(nnn,nfrequ,Eb,totsum,frequ)
   implicit none
   integer :: nnn,nfrequ
   real(8) :: Eb(nnn*nnn),totsum(nnn*nnn),frequ(nfrequ)  
  end subroutine
 end interface
  
 call sum_of_inverse_frequ(nnn,nfrequ,Eb,totsum,frequ)

 end subroutine

           !------------------------------!
           !------------------------------!
           !------------------------------!
           !------------------------------!

 subroutine cuda_array_of_inverse_c(nnn,nfrequ,Eb,totsum,frequ,firstlast )
 implicit none
 integer    :: nnn,nfrequ,firstlast
 complex(8) :: Eb(nnn,nnn),totsum(nnn,nnn),frequ(nfrequ)

 interface
  subroutine sum_of_inverse_frequ_complex(nnn,nfrequ,Eb,totsum,frequ,firstlast)
   implicit none
   integer    :: nnn,nfrequ,firstlast
   complex(8) :: Eb(nnn*nnn),totsum(nnn*nnn),frequ(nfrequ)
  end subroutine
 end interface

 call sum_of_inverse_frequ_complex(nnn,nfrequ,Eb,totsum,frequ,firstlast)

 end subroutine

           !------------------------------!
           !------------------------------!
           !------------------------------!
           !------------------------------!

 subroutine cuda_array_of_inverse_collect_c(nnn,nfrequ,Eb,cc1,cc2,frequ,firstlast)
 implicit none
 integer    :: nnn,nfrequ,firstlast
 complex(8) :: Eb(nnn,nnn),frequ(nfrequ),cc(nnn,nnn,nfrequ)
 real(8)    :: cc1(nnn,nnn,nfrequ),cc2(nnn,nnn,nfrequ)
 interface
  subroutine sum_of_inverse_frequ_complex_collect(nnn,nfrequ,Eb,cc,frequ,firstlast)
   implicit none
   integer    :: nnn,nfrequ,firstlast
   complex(8) :: Eb(nnn*nnn),cc(nnn*nnn*nfrequ),frequ(nfrequ)
  end subroutine
 end interface
 call sum_of_inverse_frequ_complex_collect(nnn,nfrequ,Eb,cc,frequ,firstlast)
 cc1=real(cc)
 cc2=aimag(cc)
 end subroutine

           !------------------------------!
           !------------------------------!
           !------------------------------!
           !------------------------------!

 subroutine cuda_array_of_inverse_array_c(nnn,nfrequ,cc,firstlast)
 implicit none
 integer    :: nnn,nfrequ,firstlast
 complex(8) :: cc(nnn,nnn,nfrequ)
 interface
  subroutine sum_of_inverse_frequ_complex_array(nnn,nfrequ,cc,firstlast)
   implicit none
   integer    :: nnn,nfrequ,firstlast
   complex(8) :: cc(nnn*nnn*nfrequ)
  end subroutine
 end interface
 call sum_of_inverse_frequ_complex_array(nnn,nfrequ,cc,firstlast)
 end subroutine

           !------------------------------!
           !------------------------------!
           !------------------------------!
           !------------------------------!

 subroutine cuda_array_of_inverse_collect_r(nnn,nfrequ,Eb,cc,frequ,firstlast)
 implicit none
 integer    :: nnn,nfrequ,firstlast
 real(8)    :: Eb(nnn,nnn),frequ(nfrequ),cc(nnn,nnn,nfrequ)
 interface
  subroutine sum_of_inverse_frequ_collect(nnn,nfrequ,Eb,cc,frequ,firstlast)
   implicit none
   integer    :: nnn,nfrequ,firstlast
   real(8)    :: Eb(nnn*nnn),cc(nnn*nnn*nfrequ),frequ(nfrequ)
  end subroutine
 end interface
 call sum_of_inverse_frequ_collect(nnn,nfrequ,Eb,cc,frequ,firstlast)
 end subroutine

           !------------------------------!
           !------------------------------!
           !------------------------------!
           !------------------------------!

 subroutine cuda_array_of_inverse_array_r(nnn,nfrequ,cc,firstlast)
 implicit none
 integer    :: nnn,nfrequ,firstlast
 real(8)    :: cc(nnn,nnn,nfrequ)
 interface
  subroutine sum_of_inverse_frequ_array(nnn,nfrequ,cc,firstlast)
   implicit none
   integer    :: nnn,nfrequ,firstlast
   real(8)    :: cc(nnn*nnn*nfrequ)
  end subroutine
 end interface
 call sum_of_inverse_frequ_array(nnn,nfrequ,cc,firstlast)
 end subroutine

           !------------------------------!
           !------------------------------!
           !------------------------------!
           !------------------------------!


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

      subroutine FFT_CUDA_real_4_(NX,BATCH,data)
      implicit none
    !----------------------------------------! 
    ! TYPICAL VALUES ARE : NX=256,BATCH=10   !
    !----------------------------------------!
      integer    :: CUFFT_FORWARD, CUFFT_INVERSE
      parameter(CUFFT_FORWARD=-1, CUFFT_INVERSE=1)
      integer    :: CUFFT_R2C, CUFFT_C2R, CUFFT_C2C
      parameter(CUFFT_R2C=X"2a", CUFFT_C2R=X"2c", CUFFT_C2C=X"29")
      integer    :: cudaMemcpyHostToDevice, cudaMemcpyDeviceToHost
      parameter(cudaMemcpyHostToDevice=1, cudaMemcpyDeviceToHost=2)
      integer    :: NX,BATCH
      real       :: PI
      parameter (PI=3.14159)
      integer    :: plan
      integer    :: i,err
      complex(4) :: data(:)

      interface
       subroutine cufftplan1d(plan,nx,cucu,batch)
        implicit none
        integer :: batch,nx,plan,cucu
       end subroutine
       subroutine cufftexecc2c(plan,data,data2,cucu,nx,batch)
        implicit none
        integer    :: batch,nx,plan,cucu
        complex(4) :: data(nx*batch),data2(nx*batch)
       end subroutine
       subroutine cufftdestroy(plan)
        implicit none
        integer :: plan
       end subroutine
      end interface

       if(size(data,1)/=NX*BATCH) stop 'FFT_CUDA_real_4_ size of data is different from NX*BATCH, critical stop'

 !     Create a 1D FFT plan. 
       call cufftplan1d(plan, NX, CUFFT_C2C, BATCH)
 !     Use the CUFFT plan to transform the signal in place.
       call cufftexecc2c(plan, data, data, CUFFT_FORWARD,NX,BATCH)
 !     Inverse transform the signal in place.
       call cufftexecc2c(plan, data, data, CUFFT_INVERSE,NX,BATCH)
 !     Destroy the CUFFT plan.
       call cufftdestroy(plan)

      return
      end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

subroutine diago_cuda_it_c(lin,mat)
implicit none
 integer                :: lin,block
 complex(8)             :: mat(lin**2),inv(lin**2)
 complex(8),allocatable :: mat_(:),inv_(:)
 integer                :: i,j,k,l,lin_

 interface
  subroutine cuda_complex_invert(k,aa,inva,nlines)
   complex(8) :: aa(2*nlines**2),inva(2*nlines**2)
   integer    :: k,nlines
  end subroutine
 end interface

 block=0
 do i=16,2,-1
  if(mod(lin,i)==0) then
   block=i
  !write(*,*) 'using block size : ', block
   exit
  endif
 enddo

 if(block==0)then
  if(allocated(mat_)) deallocate(mat_); if(allocated(inv_)) deallocate(inv_)
  lin_= lin-mod(lin,16)
  lin_= (lin_/16+1)*16

  allocate(mat_(lin_**2),inv_(lin_**2))

  do i=1,lin
   do j=1,lin
    mat_( (i-1)*lin_+j ) = mat( (i-1)*lin+j )
   enddo
  enddo
  do i=1,lin_
   do j=lin+1,lin_
    mat_( (i-1)*lin_+j ) = 0.
   enddo
  enddo
  do j=1,lin_
   do i=lin+1,lin_
    mat_( (i-1)*lin_+j ) = 0.
   enddo
  enddo
  do i=lin+1,lin_
   do j=lin+1,lin_
    mat_( (i-1)*lin_+j ) = 0.
   enddo
  enddo
  do i=lin+1,lin_
   mat_( (i-1)*lin_ + i ) = 1.d0
  enddo
  inv_=0.d0
  call cuda_complex_invert(16,mat_,inv_,lin_)
  do i=1,lin
   do j=1,lin
    mat( (i-1)*lin+j ) = inv_( (i-1)*lin_+j )
   enddo
  enddo
  if(allocated(mat_)) deallocate(mat_); if(allocated(inv_)) deallocate(inv_)

 else

  inv=0.d0
  call cuda_complex_invert(block,mat,inv,lin)
  mat=inv

 endif

return
end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

subroutine diago_cuda_it_r(lin,mat)
implicit none
 integer             :: lin,block
 real(8)             :: mat(lin**2),inv(lin**2)
 real(8),allocatable :: mat_(:),inv_(:)
 integer             :: i,j,k,l,lin_

 interface
  subroutine cuda_invert_(k,aa,inva,nlines)
   real(8) :: aa(nlines**2),inva(nlines**2)
   integer :: k,nlines
  end subroutine
 end interface

 block=0
 do i=16,2,-1
  if(mod(lin,i)==0) then
   block=i
   !write(*,*) 'using block size : ', block
   exit
  endif
 enddo

 if(block==0)then
  if(allocated(mat_)) deallocate(mat_); if(allocated(inv_)) deallocate(inv_)
  lin_= lin-mod(lin,16); lin_= (lin_/16+1)*16
  allocate(mat_(lin_**2),inv_(lin_**2))
  do i=1,lin
   do j=1,lin
    mat_( (i-1)*lin_+j ) = mat( (i-1)*lin+j )
   enddo
  enddo
  do i=1,lin_
   do j=lin+1,lin_
    mat_( (i-1)*lin_+j ) = 0.
   enddo
  enddo
  do j=1,lin_
   do i=lin+1,lin_
    mat_( (i-1)*lin_+j ) = 0.
   enddo
  enddo
  do i=lin+1,lin_
   do j=lin+1,lin_
    mat_( (i-1)*lin_+j ) = 0.
   enddo
  enddo
  do i=lin+1,lin_
   mat_( (i-1)*lin_ + i ) = 1.d0
  enddo
  inv_=0.d0
  call cuda_invert_(16,mat_,inv_,lin_)
  do i=1,lin
   do j=1,lin
    mat( (i-1)*lin+j ) = inv_( (i-1)*lin_+j )
   enddo
  enddo
  if(allocated(mat_)) deallocate(mat_); if(allocated(inv_)) deallocate(inv_)
 else 
  inv=0.d0
  call cuda_invert_(block,mat,inv,lin)
  mat=inv
 endif

end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

 subroutine test_cuda
 implicit none

 integer,parameter :: lin=31
 real(8)           :: mat(lin,lin),inv(lin,lin)
 real(4)           :: dd
 integer           :: i,j,k,l

 do i=1,lin; do j=1,lin;  call random_number(dd); mat(i,j)=dd ; enddo; enddo ; 

 if(lin<20)then
 write(*,*) '============================'
 write(*,*) 'matrix'
 do i=1,lin
  do j=1,lin
    write(*,*) i,j,mat(i,j)
  enddo
 enddo
 write(*,*) '============================'
 else
  write(*,*) 'maxval matrix input : ', maxval(abs(mat))
 endif

 inv=mat
 call diago_cuda_it_r(lin,inv)

 if(lin<20)then
 write(*,*) '============================'
 write(*,*) 'inverse matrix'
 do i=1,lin
  do j=1,lin
    write(*,*) i,j,inv(i,j)
  enddo
 enddo
 write(*,*) '============================'
 else
  write(*,*) 'maxval inv output : ', maxval(abs(inv))
 endif

 write(*,*) 'check...........'
 write(*,*) maxval(abs(MATMUL(mat,inv))),minval(abs(MATMUL(mat,inv)))
 write(*,*) 'done............'

 stop
 end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

 subroutine test_cuda_c
 implicit none

 integer,parameter :: lin=4
 complex(8)        :: mat(lin,lin),inv(lin,lin)
 real(4)           :: dd,dd1,dd2
 integer           :: i,j,k,l

 do i=1,lin; do j=1,lin; call random_number(dd1);call random_number(dd2);  mat(i,j)=cmplx(dd1,dd2,8) ; enddo; enddo ;

 write(*,*) '============================'
 write(*,*) 'matrix'
 do i=1,lin
  do j=1,lin
    write(*,*) i,j,mat(i,j)
  enddo
 enddo
 write(*,*) '============================'

 inv=mat
 call diago_cuda_it_c(lin,mat)

 write(*,*) '============================'
 write(*,*) 'inverse matrix'
 do i=1,lin
  do j=1,lin
    write(*,*) i,j,mat(i,j)
  enddo
 enddo
 write(*,*) '============================'
                                                   
 write(*,*) 'check...........'
 write(*,*) maxval(abs(MATMUL(mat,inv))),minval(abs(MATMUL(mat,inv)))
 write(*,*) 'done............'

 stop
 end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

subroutine test_cuda_MATRIXMUL
IMPLICIT NONE
INTEGER,PARAMETER :: DIMA=4,DIMB=6,DIMC=8
real(8)           :: M(DIMA,DIMB),N(DIMB,DIMC),P(DIMA,DIMC),Q(DIMA,DIMC)
REAL(4)           :: dd
INTEGER           :: I,J,K

M=0.
DO i=1,dima
do j=1,dimb
  call random_number(dd)
  M(i,j)=dd
enddo
ENDDO
N=0.
DO i=1,dimb
do j=1,dimc
  call random_number(dd)
  N(i,j)=dd
enddo
ENDDO

 P=MATMUL(M,N)
 CALL matmul_cuda(2,M,N,Q,DIMA,DIMB,DIMC)

 write(*,*) '============================='
 write(*,*) 'error : ', maxval(abs(P-Q))
 write(*,*) '============================='

STOP
END subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

subroutine MATMUL_cuda(block,M,N,Q,dima,dimb,dimc)
implicit none
INTEGER           :: DIMA,DIMB,DIMC
REAL(8)           :: M(DIMA*DIMB),N(DIMB*DIMC),Q(DIMA*DIMC)
INTEGER           :: I,J,K,block
interface
  subroutine matmul_gpu_fortran(block,M,N,Q,a,b,c)
  implicit none
    integer :: a,b,c,block
    real(8) :: M(a*b),N(b*c),Q(a*c)
   end subroutine
  end interface
  CALL matmul_gpu_fortran(block,M,N,Q,DIMA,DIMB,DIMC)
end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

subroutine MATMUL_cuda_c(block,M,N,Q,dima,dimb,dimc)
implicit none
INTEGER           :: DIMA,DIMB,DIMC
complex(8)        :: M(DIMA*DIMB),N(DIMB*DIMC),Q(DIMA*DIMC)
INTEGER           :: I,J,K,block
interface
  subroutine matmul_gpu_fortran_c(block,M,N,Q,a,b,c)
  implicit none
    integer    :: a,b,c,block
    complex(8) :: M(a*b),N(b*c),Q(a*c)
   end subroutine
  end interface
  CALL matmul_gpu_fortran_c(block,M,N,Q,DIMA,DIMB,DIMC)
end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

subroutine matmulcuda(M,N,Q,dima,dimb,dimc)
implicit none
INTEGER             :: DIMA,DIMB,DIMC,block
real(8)             :: M(DIMA,DIMB),N(DIMB,DIMC),Q(DIMA,DIMC)
real(8),allocatable :: M_(:,:),N_(:,:),Q_(:,:)
INTEGER             :: I,J,K,i_,j_

  if(allocated(M_)) deallocate(M_); if(allocated(N_)) deallocate(N_); if(allocated(Q_)) deallocate(Q_)
  i = size(M,1); i_= i-mod(i,16); i_= (i_/16+1)*16
  j = size(M,2); j_= j-mod(j,16); j_= (j_/16+1)*16
  allocate(M_(i_,j_));M_=0.; do k=1,i; M_(k,k)=1.d0; enddo; M_(1:i,1:j)=M(1:i,1:j)
  i = size(N,1); i_= i-mod(i,16); i_= (i_/16+1)*16
  j = size(N,2); j_= j-mod(j,16); j_= (j_/16+1)*16
  allocate(N_(i_,j_));N_=0.; do k=1,i; N_(k,k)=1.d0; enddo; N_(1:i,1:j)=N(1:i,1:j)
  i=size(M_,1);j=size(M_,2);k=size(N_,2)
  allocate(Q_(i,j))  
  CALL matmul_cuda(16,M_,N_,Q_,i,j,k)
  Q(1:dima,1:dimc)=Q_(1:dima,1:dimc)
  if(allocated(M_)) deallocate(M_); if(allocated(N_)) deallocate(N_); if(allocated(Q_)) deallocate(Q_)

return
end subroutine

  !-------------------------------------------------------------------------------!

subroutine matmulcuda_c(M,N,Q,dima,dimb,dimc)
implicit none
INTEGER                :: DIMA,DIMB,DIMC,block
complex(8)             :: M(DIMA,DIMB),N(DIMB,DIMC),Q(DIMA,DIMC)
complex(8),allocatable :: M_(:,:),N_(:,:),Q_(:,:)
INTEGER                :: I,J,K,i_,j_

  if(allocated(M_)) deallocate(M_); if(allocated(N_)) deallocate(N_); if(allocated(Q_)) deallocate(Q_)
  i = size(M,1); i_= i-mod(i,16); i_= (i_/16+1)*16
  j = size(M,2); j_= j-mod(j,16); j_= (j_/16+1)*16
  allocate(M_(i_,j_));M_=0.; do k=1,i; M_(k,k)=1.d0; enddo; M_(1:i,1:j)=M(1:i,1:j)
  i = size(N,1); i_= i-mod(i,16); i_= (i_/16+1)*16
  j = size(N,2); j_= j-mod(j,16); j_= (j_/16+1)*16
  allocate(N_(i_,j_));N_=0.; do k=1,i; N_(k,k)=1.d0; enddo; N_(1:i,1:j)=N(1:i,1:j)
  i=size(M_,1);j=size(M_,2);k=size(N_,2)
  allocate(Q_(i,j))  
  CALL matmul_cuda_c(16,M_,N_,Q_,i,j,k)
  Q(1:dima,1:dimc)=Q_(1:dima,1:dimc)
  if(allocated(M_)) deallocate(M_); if(allocated(N_)) deallocate(N_); if(allocated(Q_)) deallocate(Q_)

return
end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

end module

