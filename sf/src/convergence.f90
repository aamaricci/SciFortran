!+-------------------------------------------------------------------+
!PROGRAM  : 
!TYPE     : function
!PURPOSE  : 
!+-------------------------------------------------------------------+
function dv_check_convergence(Xnew,eps,N1,N2,id,index,total) result(convergence)
  integer,optional         :: id,index,total
  integer                  :: id_,index_,total_
  integer                  :: i,Msize
  logical                  :: convergence  
  integer                  :: N1,N2
  real(8)                  :: eps
  real(8)                  :: error,err
  real(8)                  :: M,S
  real(8)                  :: Xnew(:)
  real(8),save,allocatable :: Xold(:,:)
  integer,save             :: success=0,check=1
  id_=0;if(present(id))id_=id
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  if(mpiID==id_)then
     Msize=size(Xnew)
     if(.not.allocated(Xold))then
        allocate(Xold(total_,Msize))
        Xold=0.d0
     endif
     S=0.d0 ; M=0.d0
     do i=1,Msize
        M=M + abs(Xnew(i)-Xold(index_,i))
        S=S + abs(Xnew(i))
     enddo
     err= M/S
     Xold(index_,:)=Xnew
     include "convergence_write_error_file_V.f90"
     if(err < eps)success=success+1
     if(err > eps)success=0
     convergence=.false.
     include "convergence_check.f90"
     include "convergence_print_error_msg_V.f90"
     check=check+1
  endif
end function dv_check_convergence

function zv_check_convergence(Xnew,eps,N1,N2,id,index,total) result(convergence)
  integer,optional         :: id,index,total
  integer                  :: id_,index_,total_
  integer                     :: i,Msize
  logical                     :: convergence  
  integer                     :: N1,N2
  real(8)                     :: eps
  real(8)                     :: error,err
  real(8)                     :: M,S
  complex(8)                  :: Xnew(:)
  complex(8),save,allocatable :: Xold(:,:)
  integer,save                :: success=0,check=1
  character(len=2)            :: label
  id_=0;if(present(id))id_=id
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  if(mpiID==id_)then
     Msize=size(Xnew)
     if(.not.allocated(Xold))then
        allocate(Xold(total_,Msize))
        Xold=(0.d0,0.d0)
     endif
     S=0.d0 ; M=0.d0
     do i=1,Msize
        M=M + abs(Xnew(i)-Xold(index_,i))
        S=S + abs(Xnew(i))
     enddo
     err= M/S
     Xold(index_,:)=Xnew
     include "convergence_write_error_file_V.f90"
     if(err < eps)success=success+1
     if(err > eps)success=0
     convergence=.false.
     include "convergence_check.f90"
     include "convergence_print_error_msg_V.f90"
     check=check+1
  endif
end function zv_check_convergence

function dm_check_convergence(Xnew,eps,N1,N2,id,index,total,tight) result(convergence)
  integer,optional         :: id,index,total
  integer                  :: id_,index_,total_
  integer                  :: i,j,Msize1,Msize2
  logical,optional         :: tight
  logical                  :: convergence,strict
  integer                  :: N1,N2
  real(8)                  :: eps
  real(8)                  :: error(2),err
  real(8),allocatable      :: M(:),S(:),Verror(:)
  real(8)                  :: Xnew(:,:)
  real(8),save,allocatable :: Xold(:,:,:)
  integer,save             :: success=0,check=1
  character(len=2)            :: label
  strict=.false.;if(present(tight))strict=tight
  id_=0;if(present(id))id_=id
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  if(mpiID==id_)then
     Msize1=size(Xnew,1);Msize2=size(Xnew,2)
     if(.not.allocated(Xold))then
        allocate(Xold(total_,Msize1,Msize2))
        Xold=0.d0
     endif
     if(.not.allocated(M))allocate(M(Msize1))
     if(.not.allocated(S))allocate(S(Msize1))
     if(.not.allocated(Verror))allocate(Verror(Msize1))
     S=0.d0 ; M=0.d0
     do i=1,Msize2
        M=M + abs(Xnew(:,i)-Xold(index_,:,i))
        S=S + abs(Xnew(:,i))
     enddo
     Verror= M/S
     error(1)=maxval(Verror)
     error(2)=minval(Verror)
     err=sum(Verror)/dble(size(Verror))
     Xold(index_,:,:)=Xnew
     include "convergence_write_error_file_M.f90"
     if(strict)err=error(1)
     if(err < eps)success=success+1
     if(err > eps)success=0
     convergence=.false.
     include "convergence_check.f90"
     include "convergence_print_error_msg_M.f90"
     check=check+1
  endif
end function dm_check_convergence

function zm_check_convergence(Xnew,eps,N1,N2,id,index,total,tight) result(convergence)
  integer,optional            :: id,index,total
  integer                     :: id_,index_,total_
  integer                     :: i,j,Msize1,Msize2
  logical,optional            :: tight
  logical                     :: convergence,strict
  integer                     :: N1,N2
  real(8)                     :: eps
  real(8)                     :: error(2),err
  real(8),allocatable         :: M(:),S(:),Verror(:)
  complex(8)                  :: Xnew(:,:)
  complex(8),save,allocatable :: Xold(:,:,:)
  integer,save                :: success=0,check=1
  character(len=2)            :: label
  strict=.false.;if(present(tight))strict=tight
  id_=0;if(present(id))id_=id
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  if(mpiID==id_)then
     Msize1=size(Xnew,1);Msize2=size(Xnew,2)
     if(.not.allocated(Xold))then
        allocate(Xold(total_,Msize1,Msize2))
        Xold=0.d0
     endif
     if(.not.allocated(M))allocate(M(Msize1))
     if(.not.allocated(S))allocate(S(Msize1))
     if(.not.allocated(Verror))allocate(Verror(Msize1))
     S=0.d0 ; M=0.d0
     do i=1,Msize2
        M=M + abs(Xnew(:,i)-Xold(index_,:,i))
        S=S + abs(Xnew(:,i))
     enddo
     Verror= M/S
     error(1)=maxval(Verror)
     error(2)=minval(Verror)
     err=sum(Verror)/dble(size(Verror))
     Xold(index_,:,:)=Xnew
     include "convergence_write_error_file_M.f90"
     if(strict)err=error(1)
     if(err < eps)success=success+1
     if(err > eps)success=0
     convergence=.false.
     include "convergence_check.f90"
     include "convergence_print_error_msg_M.f90"
     check=check+1
  endif
  !call MPI_BCAST(convergence,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpiERR)
end function zm_check_convergence

!******************************************************************
!******************************************************************
!******************************************************************
