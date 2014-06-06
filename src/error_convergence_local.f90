!+-------------------------------------------------------------------+
!PURPOSE  : 
!+-------------------------------------------------------------------+
function i0_check_convergence_local(Xnew,eps,N1,N2,id,file,index,total,oerr,reset,extend) result(convergence)
  integer,intent(in)                           :: Xnew
  real(8),intent(inout)                        :: eps
  real(8),optional                             :: oerr
  logical,optional                             :: reset, extend
  logical                                      :: reset_, extend_
  integer,intent(inout)                        :: N1,N2
  integer,optional                             :: id,index,total
  integer                                      :: id_,index_,total_
  integer                                      :: i,j,Msize
  logical                                      :: convergence  
  real(8)                                      :: error,err
  integer,save,allocatable                     :: Xold(:)
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
  if(.not.allocated(Xold))then
     allocate(Xold(total_))
     Xold=0.d0
     include "error_init_file_dim0.f90"
  endif
  err=abs(Xnew-Xold(index_))
  if(check==1)err=1.d0
  Xold(index_)=Xnew
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
end function i0_check_convergence_local

function i1_check_convergence_local(Xnew,eps,N1,N2,id,file,index,total,strict,oerr,reset,extend) result(convergence)
  integer,intent(in)                           :: Xnew(:)
  real(8),intent(inout)                        :: eps
  real(8),optional                             :: oerr
  logical,optional                             :: reset, extend
  logical                                      :: reset_, extend_
  integer,intent(inout)                        :: N1,N2
  integer,optional                             :: id,index,total
  integer                                      :: id_,index_,total_
  integer                                      :: i,j,Msize1
  logical                                      :: convergence  
  real(8)                                      :: error(2),err
  real(8),dimension(size(Xnew))                :: M,Verror
  integer,save,allocatable                     :: Xold(:,:)
  integer,save                                 :: success=0,check=1
  character(len=2)                             :: label
  logical,optional                             :: strict
  logical                                      :: strict_
  character(len=*),optional                    :: file
  character(len=100)                           :: file_
  file_='error.err';if(present(file))file_=reg(file)
  reset_=.true.;if(present(reset))reset_=reset
  extend_=.false.;if(present(extend))extend_=extend
  strict_=.false.;if(present(strict))strict_=strict
  id_=0;if(present(id))id_=id
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  Msize1=size(Xnew)
  if(.not.allocated(Xold))then
     allocate(Xold(total_,Msize1))
     Xold=0.d0
     include "error_init_file_dim1.f90"
  endif
  Verror=abs(Xnew-Xold(index_,:))
  if(check==1)Verror=1.d0
  error(1)=maxval(Verror)
  error(2)=minval(Verror)
  err=sum(Verror)/dble(size(Verror))
  if(strict_)err=error(1)
  Xold(index_,:)=Xnew
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
end function i1_check_convergence_local

function i2_check_convergence_local(Xnew,eps,N1,N2,id,file,index,total,strict,oerr,reset,extend) result(convergence)
  integer,intent(in)                           :: Xnew(:,:)
  real(8),intent(inout)                        :: eps
  real(8),optional                             :: oerr
  logical,optional                             :: reset, extend
  logical                                      :: reset_, extend_
  integer,intent(inout)                        :: N1,N2
  integer,optional                             :: id,index,total
  integer                                      :: id_,index_,total_
  integer                                      :: i,j,Msize1,Msize2
  logical                                      :: convergence  
  real(8)                                      :: error(2),err
  real(8),dimension(size(Xnew,1),size(Xnew,2)) :: Verror
  integer,save,allocatable                     :: Xold(:,:,:)
  integer,save                                 :: success=0,check=1
  character(len=2)                             :: label
  logical,optional                             :: strict
  logical                                      :: strict_
  character(len=*),optional                    :: file
  character(len=100)                           :: file_
  file_='error.err';if(present(file))file_=reg(file)
  reset_=.true.;if(present(reset))reset_=reset
  extend_=.false.;if(present(extend))extend_=extend
  strict_=.false.;if(present(strict))strict_=strict
  id_=0;if(present(id))id_=id
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  Msize1=size(Xnew,1);Msize2=size(Xnew,2)
  if(.not.allocated(Xold))then
     allocate(Xold(total_,Msize1,Msize2))
     Xold=0.d0
     include "error_init_file_dim2.f90"
  endif
  Verror=abs(Xnew-Xold(index_,:,:))
  if(check==1)Verror=1.d0
  error(1)=maxval(Verror)
  error(2)=minval(Verror)
  err=sum(Verror)/dble(size(Verror))
  if(strict_)err=error(1)
  Xold(index_,:,:)=Xnew
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
end function i2_check_convergence_local


!----------------------------------------------------------------------


function d0_check_convergence_local(Xnew,eps,N1,N2,id,file,index,total,oerr,reset,extend) result(convergence)
  real(8),intent(in)                           :: Xnew
  real(8),intent(inout)                        :: eps
  real(8),optional                             :: oerr
  logical,optional                             :: reset, extend
  logical                                      :: reset_, extend_
  integer,intent(inout)                        :: N1,N2
  integer,optional                             :: id,index,total
  integer                                      :: id_,index_,total_
  integer                                      :: i,j,Msize
  logical                                      :: convergence  
  real(8)                                      :: error,err
  real(8),save,allocatable                     :: Xold(:)
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
  if(.not.allocated(Xold))then
     allocate(Xold(total_))
     Xold=0.d0
     include "error_init_file_dim0.f90"
  endif
  err=abs(Xnew-Xold(index_))
  if(check==1)err=1.d0
  Xold(index_)=Xnew
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
end function d0_check_convergence_local

function d1_check_convergence_local(Xnew,eps,N1,N2,id,file,index,total,strict,oerr,reset,extend) result(convergence)
  real(8),intent(in)                           :: Xnew(:)
  real(8),intent(inout)                        :: eps
  real(8),optional                             :: oerr
  logical,optional                             :: reset, extend
  logical                                      :: reset_, extend_
  integer,intent(inout)                        :: N1,N2
  integer,optional                             :: id,index,total
  integer                                      :: id_,index_,total_
  integer                                      :: i,j,Msize1
  logical                                      :: convergence  
  real(8)                                      :: error(2),err
  real(8),dimension(size(Xnew))                :: Verror
  real(8),save,allocatable                     :: Xold(:,:)
  integer,save                                 :: success=0,check=1
  character(len=2)                             :: label
  logical,optional                             :: strict
  logical                                      :: strict_
  character(len=*),optional                    :: file
  character(len=100)                           :: file_
  file_='error.err';if(present(file))file_=reg(file)
  reset_=.true.;if(present(reset))reset_=reset
  extend_=.false.;if(present(extend))extend_=extend
  strict_=.false.;if(present(strict))strict_=strict
  id_=0;if(present(id))id_=id
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  Msize1=size(Xnew)
  if(.not.allocated(Xold))then
     allocate(Xold(total_,Msize1))
     Xold=0.d0
     include "error_init_file_dim1.f90"
  endif
  Verror=abs(Xnew-Xold(index_,:))
  if(check==1)Verror=1.d0
  error(1)=maxval(Verror)
  error(2)=minval(Verror)
  err=sum(Verror)/dble(size(Verror))
  if(strict_)err=error(1)
  Xold(index_,:)=Xnew
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
end function d1_check_convergence_local

function d2_check_convergence_local(Xnew,eps,N1,N2,id,file,index,total,strict,oerr,reset,extend) result(convergence)
  real(8),intent(in)                           :: Xnew(:,:)
  real(8),intent(inout)                        :: eps
  real(8),optional                             :: oerr
  logical,optional                             :: reset, extend
  logical                                      :: reset_, extend_
  integer,intent(inout)                        :: N1,N2
  integer,optional                             :: id,index,total
  integer                                      :: id_,index_,total_
  integer                                      :: i,j,Msize1,Msize2
  logical                                      :: convergence  
  real(8)                                      :: error(2),err
  real(8),dimension(size(Xnew,1),size(Xnew,2)) :: Verror
  real(8),save,allocatable                     :: Xold(:,:,:)
  integer,save                                 :: success=0,check=1
  character(len=2)                             :: label
  logical,optional                             :: strict
  logical                                      :: strict_
  character(len=*),optional                    :: file
  character(len=100)                           :: file_
  file_='error.err';if(present(file))file_=reg(file)
  reset_=.true.;if(present(reset))reset_=reset
  extend_=.false.;if(present(extend))extend_=extend
  strict_=.false.;if(present(strict))strict_=strict
  id_=0;if(present(id))id_=id
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  Msize1=size(Xnew,1);Msize2=size(Xnew,2)
  if(.not.allocated(Xold))then
     allocate(Xold(total_,Msize1,Msize2))
     Xold=0.d0
     include "error_init_file_dim2.f90"
  endif
  Verror=abs(Xnew-Xold(index_,:,:))
  if(check==1)Verror=1.d0
  error(1)=maxval(Verror)
  error(2)=minval(Verror)
  err=sum(Verror)/dble(size(Verror))
  if(strict_)err=error(1)
  Xold(index_,:,:)=Xnew
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
end function d2_check_convergence_local


!----------------------------------------------------------------------


function z0_check_convergence_local(Xnew,eps,N1,N2,id,file,index,total,oerr,reset,extend) result(convergence)
  complex(8),intent(in)                        :: Xnew
  real(8),intent(inout)                        :: eps
  real(8),optional                             :: oerr
  logical,optional                             :: reset, extend
  logical                                      :: reset_, extend_
  integer,intent(inout)                        :: N1,N2
  integer,optional                             :: id,index,total
  integer                                      :: id_,index_,total_
  integer                                      :: i,j,Msize
  logical                                      :: convergence  
  real(8)                                      :: error,err
  complex(8),save,allocatable                  :: Xold(:)
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
  if(.not.allocated(Xold))then
     allocate(Xold(total_))
     Xold=0.d0
     include "error_init_file_dim0.f90"
  endif
  err=abs(Xnew-Xold(index_))
  if(check==1)err=1.d0
  Xold(index_)=Xnew
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
end function z0_check_convergence_local

function z1_check_convergence_local(Xnew,eps,N1,N2,id,file,index,total,strict,oerr,reset,extend) result(convergence)
  complex(8),intent(in)                        :: Xnew(:)
  real(8),intent(inout)                        :: eps
  real(8),optional                             :: oerr
  logical,optional                             :: reset, extend
  logical                                      :: reset_, extend_
  integer,intent(inout)                        :: N1,N2
  integer,optional                             :: id,index,total
  integer                                      :: id_,index_,total_
  integer                                      :: i,j,Msize1
  logical                                      :: convergence  
  real(8)                                      :: error(2),err
  real(8),dimension(size(Xnew))                :: Verror
  complex(8),save,allocatable                  :: Xold(:,:)
  integer,save                                 :: success=0,check=1
  character(len=2)                             :: label
  logical,optional                             :: strict
  logical                                      :: strict_
  character(len=*),optional                    :: file
  character(len=100)                           :: file_
  file_='error.err';if(present(file))file_=reg(file)
  reset_=.true.;if(present(reset))reset_=reset
  extend_=.false.;if(present(extend))extend_=extend
  strict_=.false.;if(present(strict))strict_=strict
  id_=0;if(present(id))id_=id
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  Msize1=size(Xnew)
  if(.not.allocated(Xold))then
     allocate(Xold(total_,Msize1))
     Xold=0.d0
     include "error_init_file_dim1.f90"
  endif
  Verror=abs(Xnew-Xold(index_,:))
  if(check==1)Verror=1.d0
  error(1)=maxval(Verror)
  error(2)=minval(Verror)
  err=sum(Verror)/dble(size(Verror))
  if(strict_)err=error(1)
  Xold(index_,:)=Xnew
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
end function z1_check_convergence_local

function z2_check_convergence_local(Xnew,eps,N1,N2,id,file,index,total,strict,oerr,reset,extend) result(convergence)
  complex(8),intent(in)                        :: Xnew(:,:)
  real(8),intent(inout)                        :: eps
  real(8),optional                             :: oerr
  logical,optional                             :: reset, extend
  logical                                      :: reset_, extend_
  integer,intent(inout)                        :: N1,N2
  integer,optional                             :: id,index,total
  integer                                      :: id_,index_,total_
  integer                                      :: i,j,Msize1,Msize2
  logical                                      :: convergence  
  real(8)                                      :: error(2),err
  real(8),dimension(size(Xnew,1),size(Xnew,2)) :: Verror
  complex(8),save,allocatable                  :: Xold(:,:,:)
  integer,save                                 :: success=0,check=1
  character(len=2)                             :: label
  logical,optional                             :: strict
  logical                                      :: strict_
  character(len=*),optional                    :: file
  character(len=100)                           :: file_
  file_='error.err';if(present(file))file_=reg(file)
  reset_=.true.;if(present(reset))reset_=reset
  extend_=.false.;if(present(extend))extend_=extend
  strict_=.false.;if(present(strict))strict_=strict
  id_=0;if(present(id))id_=id
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  Msize1=size(Xnew,1);Msize2=size(Xnew,2)
  if(.not.allocated(Xold))then
     allocate(Xold(total_,Msize1,Msize2))
     Xold=0.d0
     include "error_init_file_dim2.f90"
  endif
  Verror=abs(Xnew-Xold(index_,:,:))
  if(check==1)Verror=1.d0
  error(1)=maxval(Verror)
  error(2)=minval(Verror)
  err=sum(Verror)/dble(size(Verror))
  if(strict_)err=error(1)
  Xold(index_,:,:)=Xnew
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
end function z2_check_convergence_local
