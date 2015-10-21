!+-------------------------------------------------------------------+
!PURPOSE  : 
!+-------------------------------------------------------------------+
function i0_check_convergence_global(Xnew,eps,N1,N2,id,file,index,total,oerr,reset,extend) result(convergence)
  integer,intent(in)                           :: Xnew(:)
  real(8),intent(inout)                        :: eps
  real(8),optional                             :: oerr
  logical,optional                             :: reset, extend
  logical                                      :: reset_, extend_
  integer,intent(inout)                        :: N1,N2
  integer,optional                             :: id,index,total
  integer                                      :: id_,index_,total_
  integer                                      :: i,j,Msum
  logical                                      :: convergence  
  real(8)                                      :: error,err
  real(8)                                      :: M
  integer,save,allocatable                     :: Xold(:,:)
  integer,save                                 :: success=0,check=1
  character(len=2)                             :: label
  character(len=*),optional                    :: file
  character(len=100)                           :: file_
  file_='error.err';if(present(file))file_=reg(file)
  reset_=.true.;if(present(reset))reset_=reset
  extend_=.false.;if(present(extend))extend_=extend
  id_=0;if(present(id))id_=id
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  Msum=size(Xnew)
  if(.not.allocated(Xold))then
     allocate(Xold(total_,Msum))
     Xold=0.d0
  endif
  M=0.d0
  do i=1,Msum
     M=M + abs(Xnew(i)-Xold(index_,i))/abs(Xnew(i))
  enddo
  err= M/real(Msum,8)
  Xold(index_,:)=Xnew
  include "error_write_file_dim0.f90"
  if(err < eps)then
     success=success+1
  else
     success=0
  endif
  convergence=.false.
  include "error_test_convergence.f90"
  include "error_msg_dim0.f90"
  if(present(oerr))oerr=err
  check=check+1
end function i0_check_convergence_global

function i1_check_convergence_global(Xnew,eps,N1,N2,id,file,index,total,oerr,reset,extend) result(convergence)
  integer,intent(in)                           :: Xnew(:,:)
  real(8),intent(inout)                        :: eps
  real(8),optional                             :: oerr
  logical,optional                             :: reset, extend
  logical                                      :: reset_, extend_
  integer,intent(inout)                        :: N1,N2
  integer,optional                             :: id,index,total
  integer                                      :: id_,index_,total_
  integer                                      :: i,j,Msize1,Msum
  logical                                      :: convergence  
  real(8)                                      :: error(2),err
  real(8),dimension(size(Xnew,1))              :: M,Verror
  integer,save,allocatable                     :: Xold(:,:,:)
  integer,save                                 :: success=0,check=1
  character(len=2)                             :: label
  character(len=*),optional                    :: file
  character(len=100)                           :: file_
  file_='error.err';if(present(file))file_=reg(file)
  reset_=.true.;if(present(reset))reset_=reset
  extend_=.false.;if(present(extend))extend_=extend
  id_=0;if(present(id))id_=id
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  Msize1=size(Xnew,1) ; Msum=size(Xnew,2)
  if(.not.allocated(Xold))then
     allocate(Xold(total_,Msize1,Msum))
     Xold=0.d0
  endif
  M=0.d0 ; Verror=0.d0
  do i=1,Msum
     M=M + abs(Xnew(:,i)-Xold(index_,:,i))/abs(Xnew(:,i))
  enddo
  Verror= M/real(Msum,8)
  error(1)=maxval(Verror)
  error(2)=minval(Verror)
  err=sum(Verror)/dble(size(Verror))
  Xold(index_,:,:)=Xnew
  include "error_write_file_dim1.f90"
  if(err < eps)then
     success=success+1
  else
     success=0
  endif
  convergence=.false.
  include "error_test_convergence.f90"
  include "error_msg_dim1.f90"
  if(present(oerr))oerr=err
  check=check+1
end function i1_check_convergence_global

function i2_check_convergence_global(Xnew,eps,N1,N2,id,file,index,total,oerr,reset,extend) result(convergence)
  integer,intent(in)                           :: Xnew(:,:,:)
  real(8),intent(inout)                        :: eps
  real(8),optional                             :: oerr
  logical,optional                             :: reset, extend
  logical                                      :: reset_, extend_
  integer,intent(inout)                        :: N1,N2
  integer,optional                             :: id,index,total
  integer                                      :: id_,index_,total_
  integer                                      :: i,j,Msize1,Msize2,Msum
  logical                                      :: convergence  
  real(8)                                      :: error(2),err
  real(8),dimension(size(Xnew,1),size(Xnew,2)) :: M,Verror
  integer,save,allocatable                     :: Xold(:,:,:,:)
  integer,save                                 :: success=0,check=1
  character(len=2)                             :: label
  character(len=*),optional                    :: file
  character(len=100)                           :: file_
  file_='error.err';if(present(file))file_=reg(file)
  reset_=.true.;if(present(reset))reset_=reset
  extend_=.false.;if(present(extend))extend_=extend
  id_=0;if(present(id))id_=id
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  Msize1=size(Xnew,1);Msize2=size(Xnew,2);Msum=size(Xnew,3)
  if(.not.allocated(Xold))then
     allocate(Xold(total_,Msize1,Msize2,Msum))
     Xold=0.d0
  endif
  M=0.d0 ; Verror=0.d0
  do i=1,Msum
     M=M + abs(Xnew(:,:,i)-Xold(index_,:,:,i))/abs(Xnew(:,:,i))
  enddo
  Verror= M/real(Msum,8)
  error(1)=maxval(Verror)
  error(2)=minval(Verror)
  err=sum(Verror)/dble(size(Verror))
  Xold(index_,:,:,:)=Xnew
  include "error_write_file_dim2.f90"
  if(err < eps)then
     success=success+1
  else
     success=0
  endif
  convergence=.false.
  include "error_test_convergence.f90"
  include "error_msg_dim2.f90"
  if(present(oerr))oerr=err
  check=check+1
end function i2_check_convergence_global


!----------------------------------------------------------------------


function d0_check_convergence_global(Xnew,eps,N1,N2,id,file,index,total,oerr,reset,extend) result(convergence)
  real(8),intent(in)                           :: Xnew(:)
  real(8),intent(inout)                        :: eps
  real(8),optional                             :: oerr
  logical,optional                             :: reset, extend
  logical                                      :: reset_, extend_
  integer,intent(inout)                        :: N1,N2
  integer,optional                             :: id,index,total
  integer                                      :: id_,index_,total_
  integer                                      :: i,j,Msum
  logical                                      :: convergence  
  real(8)                                      :: error,err
  real(8)                                      :: M
  real(8),save,allocatable                     :: Xold(:,:)
  integer,save                                 :: success=0,check=1
  character(len=2)                             :: label
  character(len=*),optional                    :: file
  character(len=100)                           :: file_
  file_='error.err';if(present(file))file_=reg(file)
  reset_=.true.;if(present(reset))reset_=reset
  extend_=.false.;if(present(extend))extend_=extend
  id_=0;if(present(id))id_=id
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  Msum=size(Xnew)
  if(.not.allocated(Xold))then
     allocate(Xold(total_,Msum))
     Xold=0.d0
  endif
  M=0.d0
  do i=1,Msum
     M=M + abs(Xnew(i)-Xold(index_,i))/abs(Xnew(i))
  enddo
  err= M/real(Msum,8)
  Xold(index_,:)=Xnew
  include "error_write_file_dim0.f90"
  if(err < eps)then
     success=success+1
  else
     success=0
  endif
  convergence=.false.
  include "error_test_convergence.f90"
  include "error_msg_dim0.f90"
  if(present(oerr))oerr=err
  check=check+1
end function d0_check_convergence_global

function d1_check_convergence_global(Xnew,eps,N1,N2,id,file,index,total,oerr,reset,extend) result(convergence)
  real(8),intent(in)                           :: Xnew(:,:)
  real(8),intent(inout)                        :: eps
  real(8),optional                             :: oerr
  logical,optional                             :: reset, extend
  logical                                      :: reset_, extend_
  integer,intent(inout)                        :: N1,N2
  integer,optional                             :: id,index,total
  integer                                      :: id_,index_,total_
  integer                                      :: i,j,Msize1,Msum
  logical                                      :: convergence  
  real(8)                                      :: error(2),err
  real(8),dimension(size(Xnew,1))              :: M,Verror
  real(8),save,allocatable                     :: Xold(:,:,:)
  integer,save                                 :: success=0,check=1
  character(len=2)                             :: label
  character(len=*),optional                    :: file
  character(len=100)                           :: file_
  file_='error.err';if(present(file))file_=reg(file)
  reset_=.true.;if(present(reset))reset_=reset
  extend_=.false.;if(present(extend))extend_=extend
  id_=0;if(present(id))id_=id
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  Msize1=size(Xnew,1) ; Msum=size(Xnew,2)
  if(.not.allocated(Xold))then
     allocate(Xold(total_,Msize1,Msum))
     Xold=0.d0
  endif
  M=0.d0 ; Verror=0.d0
  do i=1,Msum
     M=M + abs(Xnew(:,i)-Xold(index_,:,i))/abs(Xnew(:,i))
  enddo
  Verror= M/real(Msum,8)
  error(1)=maxval(Verror)
  error(2)=minval(Verror)
  err=sum(Verror)/dble(size(Verror))
  Xold(index_,:,:)=Xnew
  include "error_write_file_dim1.f90"
  if(err < eps)then
     success=success+1
  else
     success=0
  endif
  convergence=.false.
  include "error_test_convergence.f90"
  include "error_msg_dim1.f90"
  if(present(oerr))oerr=err
  check=check+1
end function d1_check_convergence_global

function d2_check_convergence_global(Xnew,eps,N1,N2,id,file,index,total,oerr,reset,extend) result(convergence)
  real(8),intent(in)                           :: Xnew(:,:,:)
  real(8),intent(inout)                        :: eps
  real(8),optional                             :: oerr
  logical,optional                             :: reset, extend
  logical                                      :: reset_, extend_
  integer,intent(inout)                        :: N1,N2
  integer,optional                             :: id,index,total
  integer                                      :: id_,index_,total_
  integer                                      :: i,j,Msize1,Msize2,Msum
  logical                                      :: convergence  
  real(8)                                      :: error(2),err
  real(8),dimension(size(Xnew,1),size(Xnew,2)) :: M,Verror
  real(8),save,allocatable                     :: Xold(:,:,:,:)
  integer,save                                 :: success=0,check=1
  character(len=2)                             :: label
  character(len=*),optional                    :: file
  character(len=100)                           :: file_
  file_='error.err';if(present(file))file_=reg(file)
  reset_=.true.;if(present(reset))reset_=reset
  extend_=.false.;if(present(extend))extend_=extend
  id_=0;if(present(id))id_=id
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  Msize1=size(Xnew,1);Msize2=size(Xnew,2);Msum=size(Xnew,3)
  if(.not.allocated(Xold))then
     allocate(Xold(total_,Msize1,Msize2,Msum))
     Xold=0.d0
  endif
  M=0.d0 ; Verror=0.d0
  do i=1,Msum
     M=M + abs(Xnew(:,:,i)-Xold(index_,:,:,i))/abs(Xnew(:,:,i))
  enddo
  Verror= M/real(Msum,8)
  error(1)=maxval(Verror)
  error(2)=minval(Verror)
  err=sum(Verror)/dble(size(Verror))
  Xold(index_,:,:,:)=Xnew
  include "error_write_file_dim2.f90"
  if(err < eps)then
     success=success+1
  else
     success=0
  endif
  convergence=.false.
  include "error_test_convergence.f90"
  include "error_msg_dim2.f90"
  if(present(oerr))oerr=err
  check=check+1
end function d2_check_convergence_global


!----------------------------------------------------------------------


function z0_check_convergence_global(Xnew,eps,N1,N2,id,file,index,total,oerr,reset,extend) result(convergence)
  complex(8),intent(in)                        :: Xnew(:)
  real(8),intent(inout)                        :: eps
  real(8),optional                             :: oerr
  logical,optional                             :: reset, extend
  logical                                      :: reset_, extend_
  integer,intent(inout)                        :: N1,N2
  integer,optional                             :: id,index,total
  integer                                      :: id_,index_,total_
  integer                                      :: i,j,Msum
  logical                                      :: convergence  
  real(8)                                      :: error,err
  real(8)                                      :: M
  complex(8),save,allocatable                  :: Xold(:,:)
  integer,save                                 :: success=0,check=1
  character(len=2)                             :: label
  character(len=*),optional                    :: file
  character(len=100)                           :: file_
  file_='error.err';if(present(file))file_=reg(file)
  reset_=.true.;if(present(reset))reset_=reset
  extend_=.false.;if(present(extend))extend_=extend
  id_=0;if(present(id))id_=id
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  Msum=size(Xnew)
  if(.not.allocated(Xold))then
     allocate(Xold(total_,Msum))
     Xold=0.d0
  endif
  M=0.d0
  do i=1,Msum
     M=M + abs(Xnew(i)-Xold(index_,i))/abs(Xnew(i))
  enddo
  err= M/real(Msum,8)
  Xold(index_,:)=Xnew
  include "error_write_file_dim0.f90"
  if(err < eps)then
     success=success+1
  else
     success=0
  endif
  convergence=.false.
  include "error_test_convergence.f90"
  include "error_msg_dim0.f90"
  if(present(oerr))oerr=err
  check=check+1
end function z0_check_convergence_global

function z1_check_convergence_global(Xnew,eps,N1,N2,id,file,index,total,oerr,reset,extend) result(convergence)
  complex(8),intent(in)                        :: Xnew(:,:)
  real(8),intent(inout)                        :: eps
  real(8),optional                             :: oerr
  logical,optional                             :: reset, extend
  logical                                      :: reset_, extend_
  integer,intent(inout)                        :: N1,N2
  integer,optional                             :: id,index,total
  integer                                      :: id_,index_,total_
  integer                                      :: i,j,Msize1,Msum
  logical                                      :: convergence  
  real(8)                                      :: error(2),err
  real(8),dimension(size(Xnew,1))              :: M,Verror
  complex(8),save,allocatable                  :: Xold(:,:,:)
  integer,save                                 :: success=0,check=1
  character(len=2)                             :: label
  character(len=*),optional                    :: file
  character(len=100)                           :: file_
  file_='error.err';if(present(file))file_=reg(file)
  reset_=.true.;if(present(reset))reset_=reset
  extend_=.false.;if(present(extend))extend_=extend
  id_=0;if(present(id))id_=id
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  Msize1=size(Xnew,1) ; Msum=size(Xnew,2)
  if(.not.allocated(Xold))then
     allocate(Xold(total_,Msize1,Msum))
     Xold=0.d0
  endif
  M=0.d0 ; Verror=0.d0
  do i=1,Msum
     M=M + abs(Xnew(:,i)-Xold(index_,:,i))/abs(Xnew(:,i))
  enddo
  Verror= M/real(Msum,8)
  error(1)=maxval(Verror)
  error(2)=minval(Verror)
  err=sum(Verror)/dble(size(Verror))
  Xold(index_,:,:)=Xnew
  include "error_write_file_dim1.f90"
  if(err < eps)then
     success=success+1
  else
     success=0
  endif
  convergence=.false.
  include "error_test_convergence.f90"
  include "error_msg_dim1.f90"
  if(present(oerr))oerr=err
  check=check+1
end function z1_check_convergence_global

function z2_check_convergence_global(Xnew,eps,N1,N2,id,file,index,total,oerr,reset,extend) result(convergence)
  complex(8),intent(in)                        :: Xnew(:,:,:)
  real(8),intent(inout)                        :: eps
  real(8),optional                             :: oerr
  logical,optional                             :: reset, extend
  logical                                      :: reset_, extend_
  integer,intent(inout)                        :: N1,N2
  integer,optional                             :: id,index,total
  integer                                      :: id_,index_,total_
  integer                                      :: i,j,Msize1,Msize2,Msum
  logical                                      :: convergence  
  real(8)                                      :: error(2),err
  real(8),dimension(size(Xnew,1),size(Xnew,2)) :: M,Verror
  complex(8),save,allocatable                  :: Xold(:,:,:,:)
  integer,save                                 :: success=0,check=1
  character(len=2)                             :: label
  character(len=*),optional                    :: file
  character(len=100)                           :: file_
  file_='error.err';if(present(file))file_=reg(file)
  reset_=.true.;if(present(reset))reset_=reset
  extend_=.false.;if(present(extend))extend_=extend
  id_=0;if(present(id))id_=id
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index

  Msize1=size(Xnew,1);Msize2=size(Xnew,2);Msum=size(Xnew,3)
  if(.not.allocated(Xold))then
     allocate(Xold(total_,Msize1,Msize2,Msum))
     Xold=0.d0
  endif
  M=0.d0 ; Verror=0.d0
  do i=1,Msum
     M=M + abs(Xnew(:,:,i)-Xold(index_,:,:,i))/abs(Xnew(:,:,i))
  enddo
  Verror= M/real(Msum,8)
  error(1)=maxval(Verror)
  error(2)=minval(Verror)
  err=sum(Verror)/dble(size(Verror))
  Xold(index_,:,:,:)=Xnew
  include "error_write_file_dim2.f90"
  if(err < eps)then
     success=success+1
  else
     success=0
  endif
  convergence=.false.
  include "error_test_convergence.f90"
  include "error_msg_dim2.f90"
  if(present(oerr))oerr=err
  check=check+1
end function z2_check_convergence_global
