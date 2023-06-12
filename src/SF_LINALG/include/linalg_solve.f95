subroutine Dsolve_1rhs(A,b,trans)
  real(8),dimension(:,:),intent(in)  :: A
  real(8),dimension(:),intent(inout) :: b
  real(8),dimension(:,:),allocatable :: b_
  character(len=1),optional          :: trans
  character(len=1)                   :: trans_
  integer                            :: m,n,nrhs,lda,ldb
  integer                            :: info
  integer,dimension(:),allocatable   :: ipvt
  trans_="N";if(present(trans))trans_=trans
  lda = max(1,size(A,1))
  ldb = max(1,size(B))
  m   = size(A,1)
  n   = size(A,2)
  nrhs= 1
  allocate(ipvt(min(m,n)))
  call dgetrf(m,n,A,lda,ipvt,info)
  if(info/=0)stop "Error MATRIX/d_mat_solve_linear_system: dgetrf"    
  allocate(b_(ldb,nrhs))
  b_(:,1)=b
  call dgetrs(trans_,n,nrhs,A,lda,ipvt,b_,ldb,info)
  if(info/=0)stop "Error MATRIX/d_mat_solve_linear_system: dgetrs"
  b=b_(:,1)
  deallocate(ipvt,b_)
end subroutine dsolve_1rhs

subroutine Zsolve_1rhs(A,b,trans)
  complex(8),dimension(:,:),intent(in)  :: A
  complex(8),dimension(:),intent(inout) :: b
  complex(8),dimension(:,:),allocatable :: b_
  character(len=1),optional          :: trans
  character(len=1)                   :: trans_
  integer                            :: m,n,nrhs,lda,ldb
  integer                            :: info
  integer,dimension(:),allocatable   :: ipvt
  trans_="N";if(present(trans))trans_=trans
  lda = max(1,size(A,1))
  ldb = max(1,size(B))
  m   = size(A,1)
  n   = size(A,2)
  nrhs= 1
  allocate(ipvt(min(m,n)))
  call zgetrf(m,n,A,lda,ipvt,info)
  if(info/=0)stop "Error MATRIX/d_mat_solve_linear_system: dgetrf"    
  allocate(b_(ldb,nrhs))
  b_(:,1)=b
  call zgetrs(trans_,n,nrhs,A,lda,ipvt,b_,ldb,info)
  if(info/=0)stop "Error MATRIX/d_mat_solve_linear_system: dgetrs"
  b=b_(:,1)
  deallocate(ipvt,b_)
end subroutine Zsolve_1rhs

subroutine Dsolve_Mrhs(A,b,trans)
  real(8),dimension(:,:),intent(in)    :: A
  real(8),dimension(:,:),intent(inout) :: b
  character(len=1),optional          :: trans
  character(len=1)                   :: trans_
  integer                            :: m,n,nrhs,lda,ldb
  integer                            :: info
  integer,dimension(:),allocatable   :: ipvt
  trans_="N";if(present(trans))trans_=trans
  lda = max(1,size(A,1))
  ldb = max(1,size(B,1))
  m   = size(A,1)
  n   = size(A,2)
  nrhs= size(B,2)
  allocate(ipvt(min(m,n)))
  call dgetrf(m,n,A,lda,ipvt,info)
  if(info/=0)stop "Error MATRIX/d_mat_solve_linear_system: dgetrf"    
  call dgetrs(trans_,n,nrhs,A,lda,ipvt,b,ldb,info)
  if(info/=0)stop "Error MATRIX/d_mat_solve_linear_system: dgetrs"
  deallocate(ipvt)
end subroutine dsolve_Mrhs

subroutine Zsolve_Mrhs(A,b,trans)
  complex(8),dimension(:,:),intent(in)    :: A
  complex(8),dimension(:,:),intent(inout) :: b
  character(len=1),optional          :: trans
  character(len=1)                   :: trans_
  integer                            :: m,n,nrhs,lda,ldb
  integer                            :: info
  integer,dimension(:),allocatable   :: ipvt
  trans_="N";if(present(trans))trans_=trans
  lda = max(1,size(A,1))
  ldb = max(1,size(B,1))
  m   = size(A,1)
  n   = size(A,2)
  nrhs= size(B,2)   
  allocate(ipvt(n))
  call zgetrf(m,n,A,lda,ipvt,info)
  if(info/=0)stop "Error MATRIX/d_mat_solve_linear_system: dgetrf"    
  lda=n ; ldb=n ; nrhs=size(b,2)
  call zgetrs(trans_,n,nrhs,A,lda,ipvt,b,ldb,info)
  if(info/=0)stop "Error MATRIX/d_mat_solve_linear_system: dgetrs"
  deallocate(ipvt)
end subroutine Zsolve_Mrhs
