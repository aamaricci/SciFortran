!+-------------------------------------------------------------------+
!PROGRAM  : 
!TYPE     : function
!PURPOSE  : 
!+-------------------------------------------------------------------+
function i0_check_convergence_function_local(Xnew,eps,N1,N2,id,index,total) result(convergence)
  integer,intent(in)       :: Xnew(:)
  real(8),intent(in)       :: eps
  integer,intent(in)       :: N1,N2
  integer,optional         :: id,index,total
  integer                  :: id_,index_,total_
  integer                  :: i,j,Msum
  logical                  :: convergence  
  real(8)                  :: error,err
  real(8)                  :: M
  integer,save,allocatable :: Xold(:,:)
  integer,save             :: success=0,check=1
  character(len=2)         :: label
  id_=0;if(present(id))id_=id
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  if(mpiID==id_)then
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
     include "convergence_write_error_file_dim0.f90"
     if(err < eps)then
        success=success+1
     else
        success=0
     endif
     convergence=.false.
     include "convergence_test.f90"
     include "convergence_print_error_msg_dim0.f90"
     check=check+1
  endif
end function i0_check_convergence_function_local

function i1_check_convergence_function_local(Xnew,eps,N1,N2,id,index,total) result(convergence)
  integer,intent(in)              :: Xnew(:,:)
  real(8),intent(in)              :: eps
  integer,intent(in)              :: N1,N2
  integer,optional                :: id,index,total
  integer                         :: id_,index_,total_
  integer                         :: i,j,Msize1,Msum
  logical                         :: convergence  
  real(8)                         :: error(2),err
  real(8),dimension(size(Xnew,1)) :: M,Verror
  integer,save,allocatable        :: Xold(:,:,:)
  integer,save                    :: success=0,check=1
  character(len=2)         :: label
  id_=0;if(present(id))id_=id
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  if(mpiID==id_)then
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
     include "convergence_write_error_file_dim1.f90"
     if(err < eps)then
        success=success+1
     else
        success=0
     endif
     convergence=.false.
     include "convergence_test.f90"
     include "convergence_print_error_msg_dim1.f90"
     check=check+1
  endif
end function i1_check_convergence_function_local

function i2_check_convergence_function_local(Xnew,eps,N1,N2,id,index,total) result(convergence)
  integer,intent(in)                                        :: Xnew(:,:,:)
  real(8),intent(in)                                        :: eps
  integer,intent(in)                                        :: N1,N2
  integer,optional                                          :: id,index,total
  integer                                                   :: id_,index_,total_
  integer                                                   :: i,j,Msize1,Msize2,Msum
  logical                                                   :: convergence  
  real(8)                                                   :: error(2),err
  real(8),dimension(size(Xnew,1),size(Xnew,2)) :: M,Verror
  integer,save,allocatable                                  :: Xold(:,:,:,:)
  integer,save                                              :: success=0,check=1
  character(len=2)         :: label
  id_=0;if(present(id))id_=id
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  if(mpiID==id_)then
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
     include "convergence_write_error_file_dim2.f90"
     if(err < eps)then
        success=success+1
     else
        success=0
     endif
     convergence=.false.
     include "convergence_test.f90"
     include "convergence_print_error_msg_dim2.f90"
     check=check+1
  endif
end function i2_check_convergence_function_local


!----------------------------------------------------------------------


function d0_check_convergence_function_local(Xnew,eps,N1,N2,id,index,total) result(convergence)
  real(8),intent(in)       :: Xnew(:)
  real(8),intent(in)       :: eps
  integer,intent(in)       :: N1,N2
  integer,optional         :: id,index,total
  integer                  :: id_,index_,total_
  integer                  :: i,j,Msum
  logical                  :: convergence  
  real(8)                  :: error,err
  real(8)                  :: M
  real(8),save,allocatable :: Xold(:,:)
  integer,save             :: success=0,check=1
  character(len=2)         :: label
  id_=0;if(present(id))id_=id
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  if(mpiID==id_)then
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
     include "convergence_write_error_file_dim0.f90"
     if(err < eps)then
        success=success+1
     else
        success=0
     endif
     convergence=.false.
     include "convergence_test.f90"
     include "convergence_print_error_msg_dim0.f90"
     check=check+1
  endif
end function d0_check_convergence_function_local

function d1_check_convergence_function_local(Xnew,eps,N1,N2,id,index,total) result(convergence)
  real(8),intent(in)              :: Xnew(:,:)
  real(8),intent(in)              :: eps
  integer,intent(in)              :: N1,N2
  integer,optional                :: id,index,total
  integer                         :: id_,index_,total_
  integer                         :: i,j,Msize1,Msum
  logical                         :: convergence  
  real(8)                         :: error(2),err
  real(8),dimension(size(Xnew,1)) :: M,Verror
  real(8),save,allocatable        :: Xold(:,:,:)
  integer,save                    :: success=0,check=1
  character(len=2)         :: label
  id_=0;if(present(id))id_=id
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  if(mpiID==id_)then
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
     include "convergence_write_error_file_dim1.f90"
     if(err < eps)then
        success=success+1
     else
        success=0
     endif
     convergence=.false.
     include "convergence_test.f90"
     include "convergence_print_error_msg_dim1.f90"
     check=check+1
  endif
end function d1_check_convergence_function_local

function d2_check_convergence_function_local(Xnew,eps,N1,N2,id,index,total) result(convergence)
  real(8),intent(in)                                        :: Xnew(:,:,:)
  real(8),intent(in)                                        :: eps
  integer,intent(in)                                        :: N1,N2
  integer,optional                                          :: id,index,total
  integer                                                   :: id_,index_,total_
  integer                                                   :: i,j,Msize1,Msize2,Msum
  logical                                                   :: convergence  
  real(8)                                                   :: error(2),err
  real(8),dimension(size(Xnew,1),size(Xnew,2)) :: M,Verror
  real(8),save,allocatable                                  :: Xold(:,:,:,:)
  integer,save                                              :: success=0,check=1
  character(len=2)         :: label
  id_=0;if(present(id))id_=id
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  if(mpiID==id_)then
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
     include "convergence_write_error_file_dim2.f90"
     if(err < eps)then
        success=success+1
     else
        success=0
     endif
     convergence=.false.
     include "convergence_test.f90"
     include "convergence_print_error_msg_dim2.f90"
     check=check+1
  endif
end function d2_check_convergence_function_local


!----------------------------------------------------------------------


function z0_check_convergence_function_local(Xnew,eps,N1,N2,id,index,total) result(convergence)
  complex(8),intent(in)       :: Xnew(:)
  real(8),intent(in)          :: eps
  integer,intent(in)          :: N1,N2
  integer,optional            :: id,index,total
  integer                     :: id_,index_,total_
  integer                     :: i,j,Msum
  logical                     :: convergence  
  real(8)                     :: error,err
  real(8)                     :: M
  complex(8),save,allocatable :: Xold(:,:)
  integer,save                :: success=0,check=1
  character(len=2)         :: label
  id_=0;if(present(id))id_=id
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  if(mpiID==id_)then
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
     include "convergence_write_error_file_dim0.f90"
     if(err < eps)then
        success=success+1
     else
        success=0
     endif
     convergence=.false.
     include "convergence_test.f90"
     include "convergence_print_error_msg_dim0.f90"
     check=check+1
  endif
end function z0_check_convergence_function_local

function z1_check_convergence_function_local(Xnew,eps,N1,N2,id,index,total) result(convergence)
  complex(8),intent(in)           :: Xnew(:,:)
  real(8),intent(in)              :: eps
  integer,intent(in)              :: N1,N2
  integer,optional                :: id,index,total
  integer                         :: id_,index_,total_
  integer                         :: i,j,Msize1,Msum
  logical                         :: convergence  
  real(8)                         :: error(2),err
  real(8),dimension(size(Xnew,1)) :: M,Verror
  complex(8),save,allocatable     :: Xold(:,:,:)
  integer,save                    :: success=0,check=1
  character(len=2)         :: label
  id_=0;if(present(id))id_=id
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  if(mpiID==id_)then
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
     include "convergence_write_error_file_dim1.f90"
     if(err < eps)then
        success=success+1
     else
        success=0
     endif
     convergence=.false.
     include "convergence_test.f90"
     include "convergence_print_error_msg_dim1.f90"
     check=check+1
  endif
end function z1_check_convergence_function_local

function z2_check_convergence_function_local(Xnew,eps,N1,N2,id,index,total) result(convergence)
  complex(8),intent(in)                                     :: Xnew(:,:,:)
  real(8),intent(in)                                        :: eps
  integer,intent(in)                                        :: N1,N2
  integer,optional                                          :: id,index,total
  integer                                                   :: id_,index_,total_
  integer                                                   :: i,j,Msize1,Msize2,Msum
  logical                                                   :: convergence  
  real(8)                                                   :: error(2),err
  real(8),dimension(size(Xnew,1),size(Xnew,2)) :: M,Verror
  complex(8),save,allocatable                               :: Xold(:,:,:,:)
  integer,save                                              :: success=0,check=1
  character(len=2)         :: label
  id_=0;if(present(id))id_=id
  total_=1;if(present(total))total_=total
  index_=1;if(present(index))index_=index
  if(mpiID==id_)then
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
     include "convergence_write_error_file_dim2.f90"
     if(err < eps)then
        success=success+1
     else
        success=0
     endif
     convergence=.false.
     include "convergence_test.f90"
     include "convergence_print_error_msg_dim2.f90"
     check=check+1
  endif
end function z2_check_convergence_function_local
