!+-------------------------------------------------------------------+
!PURPOSE  : 
!+-------------------------------------------------------------------+
function i0_check_convergence_scalar(Xnew,eps,N1,N2,id,file,index,total) result(convergence)
  integer,intent(in)       :: Xnew
  real(8),intent(in)       :: eps
  integer,intent(in)       :: N1,N2
  integer,optional         :: id,index,total
  integer                  :: id_,index_,total_
  integer                  :: i,j,Msize
  logical                  :: convergence  
  real(8)                  :: error,err
  integer,save,allocatable :: Xold(:)
  integer,save             :: success=0,check=1
  character(len=2)         :: label
  character(len=*),optional:: file
  character(len=100)       :: file_
  file_='error.err';if(present(file))file_=reg(file)
  id_=0;if(present(id))id_=id
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  if(mpiID==id_)then
     if(.not.allocated(Xold))then
        allocate(Xold(total_))
        Xold=0.d0
     endif
     err=abs(Xnew-Xold(index_))
     if(check==1)err=1.d0
     Xold(index_)=Xnew
     include "tools_convergence/write_error_file_dim0.f90"
     if(err < eps)then
        success=success+1
     else
        success=0
     endif
     convergence=.false.
     include "tools_convergence/convergence_test.f90"
     include "tools_convergence/print_error_msg_dim0.f90"
     check=check+1
  endif
end function i0_check_convergence_scalar

function i1_check_convergence_scalar(Xnew,eps,N1,N2,id,file,index,total,strict) result(convergence)
  integer,intent(in)            :: Xnew(:)
  real(8),intent(in)            :: eps
  integer,intent(in)            :: N1,N2
  integer,optional              :: id,index,total
  integer                       :: id_,index_,total_
  integer                       :: i,j,Msize1
  logical                       :: convergence  
  real(8)                       :: error(2),err
  real(8),dimension(size(Xnew)) :: M,Verror
  integer,save,allocatable      :: Xold(:,:)
  integer,save                  :: success=0,check=1
  character(len=2)              :: label
  logical,optional              :: strict
  logical                       :: strict_
  character(len=*),optional:: file
  character(len=100)       :: file_
  file_='error.err';if(present(file))file_=reg(file)
  strict_=.false.;if(present(strict))strict_=strict
  id_=0;if(present(id))id_=id
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  if(mpiID==id_)then
     Msize1=size(Xnew)
     if(.not.allocated(Xold))then
        allocate(Xold(total_,Msize1))
        Xold=0.d0
     endif
     Verror=abs(Xnew-Xold(index_,:))
     if(check==1)Verror=1.d0
     error(1)=maxval(Verror)
     error(2)=minval(Verror)
     err=sum(Verror)/dble(size(Verror))
     if(strict_)err=error(1)
     Xold(index_,:)=Xnew
     include "tools_convergence/write_error_file_dim1.f90"
     if(err < eps)then
        success=success+1
     else
        success=0
     endif
     convergence=.false.
     include "tools_convergence/convergence_test.f90"
     include "tools_convergence/print_error_msg_dim1.f90"
     check=check+1
  endif
end function i1_check_convergence_scalar

function i2_check_convergence_scalar(Xnew,eps,N1,N2,id,file,index,total,strict) result(convergence)
  integer,intent(in)                           :: Xnew(:,:)
  real(8),intent(in)                           :: eps
  integer,intent(in)                           :: N1,N2
  integer,optional                             :: id,index,total
  integer                                      :: id_,index_,total_
  integer                                      :: i,j,Msize1,Msize2
  logical                                      :: convergence  
  real(8)                                      :: error(2),err
  real(8),dimension(size(Xnew,1),size(Xnew,2)) :: Verror
  integer,save,allocatable                     :: Xold(:,:,:)
  integer,save                                 :: success=0,check=1
  character(len=2)         :: label
  logical,optional              :: strict
  logical                       :: strict_
  character(len=*),optional:: file
  character(len=100)       :: file_
  file_='error.err';if(present(file))file_=reg(file)
  strict_=.false.;if(present(strict))strict_=strict
  id_=0;if(present(id))id_=id
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  if(mpiID==id_)then
     Msize1=size(Xnew,1);Msize2=size(Xnew,2)
     if(.not.allocated(Xold))then
        allocate(Xold(total_,Msize1,Msize2))
        Xold=0.d0
     endif
     Verror=abs(Xnew-Xold(index_,:,:))
     if(check==1)Verror=1.d0
     error(1)=maxval(Verror)
     error(2)=minval(Verror)
     err=sum(Verror)/dble(size(Verror))
     if(strict_)err=error(1)
     Xold(index_,:,:)=Xnew
     include "tools_convergence/write_error_file_dim2.f90"
     if(err < eps)then
        success=success+1
     else
        success=0
     endif
     convergence=.false.
     include "tools_convergence/convergence_test.f90"
     include "tools_convergence/print_error_msg_dim2.f90"
     check=check+1
  endif
end function i2_check_convergence_scalar


!----------------------------------------------------------------------


function d0_check_convergence_scalar(Xnew,eps,N1,N2,id,file,index,total) result(convergence)
  real(8),intent(in)       :: Xnew
  real(8),intent(in)       :: eps
  integer,intent(in)       :: N1,N2
  integer,optional         :: id,index,total
  integer                  :: id_,index_,total_
  integer                  :: i,j,Msize
  logical                  :: convergence  
  real(8)                  :: error,err
  real(8),save,allocatable :: Xold(:)
  integer,save             :: success=0,check=1
  character(len=2)         :: label
  character(len=*),optional:: file
  character(len=100)       :: file_
  file_='error.err';if(present(file))file_=reg(file)
  id_=0;if(present(id))id_=id
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  if(mpiID==id_)then
     if(.not.allocated(Xold))then
        allocate(Xold(total_))
        Xold=0.d0
     endif
     err=abs(Xnew-Xold(index_))
     if(check==1)err=1.d0
     Xold(index_)=Xnew
     include "tools_convergence/write_error_file_dim0.f90"
     if(err < eps)then
        success=success+1
     else
        success=0
     endif
     convergence=.false.
     include "tools_convergence/convergence_test.f90"
     include "tools_convergence/print_error_msg_dim0.f90"
     check=check+1
  endif
end function d0_check_convergence_scalar

function d1_check_convergence_scalar(Xnew,eps,N1,N2,id,file,index,total,strict) result(convergence)
  real(8),intent(in)            :: Xnew(:)
  real(8),intent(in)            :: eps
  integer,intent(in)            :: N1,N2
  integer,optional              :: id,index,total
  integer                       :: id_,index_,total_
  integer                       :: i,j,Msize1
  logical                       :: convergence  
  real(8)                       :: error(2),err
  real(8),dimension(size(Xnew)) :: Verror
  real(8),save,allocatable      :: Xold(:,:)
  integer,save                  :: success=0,check=1
  character(len=2)         :: label
  logical,optional              :: strict
  logical                       :: strict_
  character(len=*),optional:: file
  character(len=100)       :: file_
  file_='error.err';if(present(file))file_=reg(file)
  strict_=.false.;if(present(strict))strict_=strict
  id_=0;if(present(id))id_=id
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  if(mpiID==id_)then
     Msize1=size(Xnew)
     if(.not.allocated(Xold))then
        allocate(Xold(total_,Msize1))
        Xold=0.d0
     endif
     Verror=abs(Xnew-Xold(index_,:))
     if(check==1)Verror=1.d0
     error(1)=maxval(Verror)
     error(2)=minval(Verror)
     err=sum(Verror)/dble(size(Verror))
     if(strict_)err=error(1)
     Xold(index_,:)=Xnew
     include "tools_convergence/write_error_file_dim1.f90"
     if(err < eps)then
        success=success+1
     else
        success=0
     endif
     convergence=.false.
     include "tools_convergence/convergence_test.f90"
     include "tools_convergence/print_error_msg_dim1.f90"
     check=check+1
  endif
end function d1_check_convergence_scalar

function d2_check_convergence_scalar(Xnew,eps,N1,N2,id,file,index,total,strict) result(convergence)
  real(8),intent(in)                           :: Xnew(:,:)
  real(8),intent(in)                           :: eps
  integer,intent(in)                           :: N1,N2
  integer,optional                             :: id,index,total
  integer                                      :: id_,index_,total_
  integer                                      :: i,j,Msize1,Msize2
  logical                                      :: convergence  
  real(8)                                      :: error(2),err
  real(8),dimension(size(Xnew,1),size(Xnew,2)) :: Verror
  real(8),save,allocatable                     :: Xold(:,:,:)
  integer,save                                 :: success=0,check=1
  character(len=2)         :: label
  logical,optional              :: strict
  logical                       :: strict_
  character(len=*),optional:: file
  character(len=100)       :: file_
  file_='error.err';if(present(file))file_=reg(file)
  strict_=.false.;if(present(strict))strict_=strict
  id_=0;if(present(id))id_=id
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  if(mpiID==id_)then
     Msize1=size(Xnew,1);Msize2=size(Xnew,2)
     if(.not.allocated(Xold))then
        allocate(Xold(total_,Msize1,Msize2))
        Xold=0.d0
     endif
     Verror=abs(Xnew-Xold(index_,:,:))
     if(check==1)Verror=1.d0
     error(1)=maxval(Verror)
     error(2)=minval(Verror)
     err=sum(Verror)/dble(size(Verror))
     if(strict_)err=error(1)
     Xold(index_,:,:)=Xnew
     include "tools_convergence/write_error_file_dim2.f90"
     if(err < eps)then
        success=success+1
     else
        success=0
     endif
     convergence=.false.
     include "tools_convergence/convergence_test.f90"
     include "tools_convergence/print_error_msg_dim2.f90"
     check=check+1
  endif
end function d2_check_convergence_scalar


!----------------------------------------------------------------------


function z0_check_convergence_scalar(Xnew,eps,N1,N2,id,file,index,total) result(convergence)
  complex(8),intent(in)       :: Xnew
  real(8),intent(in)          :: eps
  integer,intent(in)          :: N1,N2
  integer,optional            :: id,index,total
  integer                     :: id_,index_,total_
  integer                     :: i,j,Msize
  logical                     :: convergence  
  real(8)                     :: error,err
  complex(8),save,allocatable :: Xold(:)
  integer,save                :: success=0,check=1
  character(len=2)         :: label
  character(len=*),optional:: file
  character(len=100)       :: file_
  file_='error.err';if(present(file))file_=reg(file)
  id_=0;if(present(id))id_=id
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  if(mpiID==id_)then
     if(.not.allocated(Xold))then
        allocate(Xold(total_))
        Xold=0.d0
     endif
     err=abs(Xnew-Xold(index_))
     if(check==1)err=1.d0
     Xold(index_)=Xnew
     include "tools_convergence/write_error_file_dim0.f90"
     if(err < eps)then
        success=success+1
     else
        success=0
     endif
     convergence=.false.
     include "tools_convergence/convergence_test.f90"
     include "tools_convergence/print_error_msg_dim0.f90"
     check=check+1
  endif
end function z0_check_convergence_scalar

function z1_check_convergence_scalar(Xnew,eps,N1,N2,id,file,index,total,strict) result(convergence)
  complex(8),intent(in)         :: Xnew(:)
  real(8),intent(in)            :: eps
  integer,intent(in)            :: N1,N2
  integer,optional              :: id,index,total
  integer                       :: id_,index_,total_
  integer                       :: i,j,Msize1
  logical                       :: convergence  
  real(8)                       :: error(2),err
  real(8),dimension(size(Xnew)) :: Verror
  complex(8),save,allocatable   :: Xold(:,:)
  integer,save                  :: success=0,check=1
  character(len=2)         :: label
  logical,optional              :: strict
  logical                       :: strict_
  character(len=*),optional:: file
  character(len=100)       :: file_
  file_='error.err';if(present(file))file_=reg(file)
  strict_=.false.;if(present(strict))strict_=strict
  id_=0;if(present(id))id_=id
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  if(mpiID==id_)then
     Msize1=size(Xnew)
     if(.not.allocated(Xold))then
        allocate(Xold(total_,Msize1))
        Xold=0.d0
     endif
     Verror=abs(Xnew-Xold(index_,:))
     if(check==1)Verror=1.d0
     error(1)=maxval(Verror)
     error(2)=minval(Verror)
     err=sum(Verror)/dble(size(Verror))
     if(strict_)err=error(1)
     Xold(index_,:)=Xnew
     include "tools_convergence/write_error_file_dim1.f90"
     if(err < eps)then
        success=success+1
     else
        success=0
     endif
     convergence=.false.
     include "tools_convergence/convergence_test.f90"
     include "tools_convergence/print_error_msg_dim1.f90"
     check=check+1
  endif
end function z1_check_convergence_scalar

function z2_check_convergence_scalar(Xnew,eps,N1,N2,id,file,index,total,strict) result(convergence)
  complex(8),intent(in)                        :: Xnew(:,:)
  real(8),intent(in)                           :: eps
  integer,intent(in)                           :: N1,N2
  integer,optional                             :: id,index,total
  integer                                      :: id_,index_,total_
  integer                                      :: i,j,Msize1,Msize2
  logical                                      :: convergence  
  real(8)                                      :: error(2),err
  real(8),dimension(size(Xnew,1),size(Xnew,2)) :: Verror
  complex(8),save,allocatable                  :: Xold(:,:,:)
  integer,save                                 :: success=0,check=1
  character(len=2)         :: label
  logical,optional              :: strict
  logical                       :: strict_
  character(len=*),optional:: file
  character(len=100)       :: file_
  file_='error.err';if(present(file))file_=reg(file)
  strict_=.false.;if(present(strict))strict_=strict
  id_=0;if(present(id))id_=id
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  if(mpiID==id_)then
     Msize1=size(Xnew,1);Msize2=size(Xnew,2)
     if(.not.allocated(Xold))then
        allocate(Xold(total_,Msize1,Msize2))
        Xold=0.d0
     endif
     Verror=abs(Xnew-Xold(index_,:,:))
     if(check==1)Verror=1.d0
     error(1)=maxval(Verror)
     error(2)=minval(Verror)
     err=sum(Verror)/dble(size(Verror))
     if(strict_)err=error(1)
     Xold(index_,:,:)=Xnew
     include "tools_convergence/write_error_file_dim2.f90"
     if(err < eps)then
        success=success+1
     else
        success=0
     endif
     convergence=.false.
     include "tools_convergence/convergence_test.f90"
     include "tools_convergence/print_error_msg_dim2.f90"
     check=check+1
  endif
end function z2_check_convergence_scalar
