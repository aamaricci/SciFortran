MODULE SF_SPARSE_ARRAY_COO
  USE SF_SPARSE_COMMON
#ifdef _MPI
  USE MPI
#endif
  implicit none
  public :: sf_sparse_dmatrix_coo,sf_sparse_cmatrix_coo, sf_init_matrix_coo, sf_delete_matrix_coo, sf_insert_element_coo, sf_dump_matrix_coo

!COO MATRIX
  type sf_sparse_dmatrix_coo
     integer,dimension(:),allocatable :: rows
     integer,dimension(:),allocatable :: cols
     real(8),dimension(:),allocatable :: vals
     integer                                   :: Size
     integer                                   :: Nrow
     integer                                   :: Ncol
     logical                                   :: status=.false.
! #ifdef _MPI
!      type(sf_sparse_drow_coo),dimension(:),pointer :: loc
!      integer                                   :: istart=0 !global start index for MPI storage
!      integer                                   :: iend=0
!      integer                                   :: ishift=0
!      logical                                   :: mpi=.false.
! #endif
  end type sf_sparse_dmatrix_coo
  
  type sf_sparse_cmatrix_coo
     integer,dimension(:),allocatable     :: rows
     integer,dimension(:),allocatable     :: cols
     complex(8),dimension(:),allocatable :: vals
     integer                                   :: Size
     integer                                   :: Nrow
     integer                                   :: Ncol
     logical                                   :: status=.false.
! #ifdef _MPI
!      type(sf_sparse_crow_coo),dimension(:),pointer :: loc
!      integer                                   :: istart=0 !global start index for MPI storage
!      integer                                   :: iend=0
!      integer                                   :: ishift=0
!      logical                                   :: mpi=.false.
! #endif
  end type sf_sparse_cmatrix_coo
  
  !INIT SPARSE MATRICES 
  interface sf_init_matrix_coo
     module procedure :: sf_init_dmatrix_coo
     module procedure :: sf_init_cmatrix_coo
! #ifdef _MPI
!      module procedure :: mpi_sf_init_dmatrix_coo
!      module procedure :: mpi_sf_init_cmatrix_coo
! #endif
  end interface sf_init_matrix_coo
  
  !DELETE SPARSE MATRIX 
  interface sf_delete_matrix_coo
     module procedure :: sf_delete_dmatrix_coo
     module procedure :: sf_delete_cmatrix_coo
! #ifdef _MPI
!      module procedure :: mpi_sf_delete_dmatrix_coo
!      module procedure :: mpi_sf_delete_cmatrix_coo
! #endif
  end interface sf_delete_matrix_coo


  !INSERT ELEMENTS
  interface sf_insert_element_coo
     module procedure :: sf_insert_delement_coo
     module procedure :: sf_insert_celement_coo
! #ifdef _MPI
!      module procedure :: mpi_sf_insert_delement_coo
!      module procedure :: mpi_sf_insert_celement_coo
! #endif
  end interface sf_insert_element_coo

  
  !DUMP SPARSE MATRIX INTO STANDARD MATRIX
  interface sf_dump_matrix_coo
     module procedure :: sf_dump_dmatrix_coo
     module procedure :: sf_dump_cmatrix_coo
! #ifdef _MPI
!      module procedure :: mpi_sf_dump_matrix_coo
!      module procedure :: mpi_sf_dump_matrix_coo
! #endif
  end interface sf_dump_matrix_coo

contains  
  !+------------------------------------------------------------------+
  !PURPOSE:  initialize the sparse matrix list
  !+------------------------------------------------------------------+
  subroutine sf_init_dmatrix_coo(sparse,N,N1)
    type(sf_sparse_dmatrix_coo),intent(inout) :: sparse
    integer                               :: N
    integer,optional                      :: N1
    integer                               :: i
    !
#ifdef _DEBUG
    write(*,"(A)")"DEBUG sf_init_dmatrix_coo: allocate sparse"
#endif
    !put here a delete statement to avoid problems
    if(sparse%status)stop "sf_init_dmatrix_coo: already allocated can not init"
    !
    sparse%Nrow=N
    sparse%Ncol=N
    sparse%Size=0
    if(present(N1))sparse%Ncol=N1
    !
    allocate(sparse%rows(0))
    allocate(sparse%cols(0))
    allocate(sparse%vals(0))
    !
    sparse%status=.true.
    !
  end subroutine sf_init_dmatrix_coo
  
subroutine sf_init_cmatrix_coo(sparse,N,N1)
    type(sf_sparse_cmatrix_coo),intent(inout) :: sparse
    integer                               :: N
    integer,optional                      :: N1
    integer                               :: i
    !
#ifdef _DEBUG
    write(*,"(A)")"DEBUG sf_init_cmatrix_coo: allocate sparse"
#endif
    !put here a delete statement to avoid problems
    if(sparse%status)stop "sf_init_cmatrix_coo: already allocated can not init"
    !
    sparse%Nrow=N
    sparse%Ncol=N
    sparse%Size=0
    if(present(N1))sparse%Ncol=N1
    !
    allocate(sparse%rows(0))
    allocate(sparse%cols(0))
    allocate(sparse%vals(0))
    !
    sparse%status=.true.
    !
  end subroutine sf_init_cmatrix_coo



! #ifdef _MPI
!   subroutine mpi_sf_init_dmatrix_coo(MpiComm,sparse,N,N1)
!     integer                               :: MpiComm
!     type(sf_sparse_dmatrix_coo),intent(inout) :: sparse
!     integer                               :: N
!     integer,optional                      :: N1
!     integer                               :: i,Ncol,Nloc
!     !
! #ifdef _DEBUG
!     write(*,"(A)")"DEBUG MPI_sf_init_dmatrix_coo: allocate sparse"
! #endif
!     if(MpiComm==Mpi_Comm_Null)return
!     !
!     call sf_test_matrix_mpi(MpiComm,sparse,"mpi_sf_init_dmatrix_coo")
!     !
!     Ncol = N
!     if(present(N1))Ncol=N1
!     !
!     Nloc = sparse%iend-sparse%istart+1
!     !
!     call sf_init_matrix_coo(sparse,Nloc,Ncol)
!     !
!     allocate(sparse%loc(Nloc))
!     do i=1,Nloc
!        sparse%loc(i)%size=0
!        allocate(sparse%loc(i)%vals(0)) !empty array
!        allocate(sparse%loc(i)%cols(0)) !empty array
!     end do
!     !
!   end subroutine mpi_sf_init_dmatrix_coo
  
!   subroutine mpi_sf_init_cmatrix_coo(MpiComm,sparse,N,N1)
!     integer                               :: MpiComm
!     type(sf_sparse_cmatrix_coo),intent(inout) :: sparse
!     integer                               :: N
!     integer,optional                      :: N1
!     integer                               :: i,Ncol,Nloc
!     !
! #ifdef _DEBUG
!     write(*,"(A)")"DEBUG MPI_sf_init_cmatrix_coo: allocate sparse"
! #endif
!     if(MpiComm==Mpi_Comm_Null)return
!     !
!     call sf_test_matrix_mpi(MpiComm,sparse,"mpi_sf_init_cmatrix_coo")
!     !
!     Ncol = N
!     if(present(N1))Ncol=N1
!     !
!     Nloc = sparse%iend-sparse%istart+1
!     !
!     call sf_init_matrix_coo(sparse,Nloc,Ncol)
!     !
!     allocate(sparse%loc(Nloc))
!     do i=1,Nloc
!        sparse%loc(i)%size=0
!        allocate(sparse%loc(i)%vals(0)) !empty array
!        allocate(sparse%loc(i)%cols(0)) !empty array
!     end do
!     !
!   end subroutine mpi_sf_init_cmatrix_coo
! #endif

  
  !+------------------------------------------------------------------+
  !PURPOSE: delete an entire sparse matrix
  !+------------------------------------------------------------------+
  subroutine sf_delete_dmatrix_coo(sparse)
    type(sf_sparse_dmatrix_coo),intent(inout) :: sparse
    integer                           :: i
    !
#ifdef _DEBUG
    write(*,"(A)")"DEBUG sf_delete_dmatrix_coo: delete sparse"
#endif
    if(.not.sparse%status)return !stop "Error SPARSE/sf_delete_matrix: sparse is not allocated."
    !
    deallocate(sparse%rows)
    deallocate(sparse%cols)
    deallocate(sparse%vals)
    !
    sparse%Nrow=0
    sparse%Ncol=0
    sparse%Size=0
    sparse%status=.false.
  end subroutine sf_delete_dmatrix_coo

    subroutine sf_delete_cmatrix_coo(sparse)
    type(sf_sparse_cmatrix_coo),intent(inout) :: sparse
    integer                           :: i
    !
#ifdef _DEBUG
    write(*,"(A)")"DEBUG sf_delete_cmatrix_coo: delete sparse"
#endif
    if(.not.sparse%status)return !stop "Error SPARSE/sf_delete_matrix: sparse is not allocated."
    !
    deallocate(sparse%rows)
    deallocate(sparse%cols)
    deallocate(sparse%vals)
    !
    sparse%Nrow=0
    sparse%Ncol=0
    sparse%Size=0
    sparse%status=.false.
  end subroutine sf_delete_cmatrix_coo



! #ifdef _MPI
!   subroutine mpi_sf_delete_dmatrix_coo(MpiComm,sparse)
!     integer                              :: MpiComm
!     type(sf_sparse_dmatrix_coo),intent(inout) :: sparse
!     integer                              :: i
!     type(sf_sparse_drow_coo),pointer          :: row
!     !
! #ifdef _DEBUG
!     write(*,"(A)")"DEBUG MPI_sf_delete_dmatrix_coo: delete sparse"
! #endif
!     if(.not.sparse%status)return !stop "Error SPARSE/mpi_sp_delete_matrix: sparse is not allocated."
!     !
!     do i=1,sparse%Nrow
!        deallocate(sparse%row(i)%vals)
!        deallocate(sparse%row(i)%cols)
!        sparse%row(i)%Size  = 0
!        !
!        deallocate(sparse%loc(i)%vals)
!        deallocate(sparse%loc(i)%cols)
!        sparse%loc(i)%Size  = 0
!     enddo
!     deallocate(sparse%row)
!     deallocate(sparse%loc)
!     !
!     sparse%Nrow=0
!     sparse%Ncol=0
!     sparse%status=.false.
!     !
!     sparse%istart=0
!     sparse%iend=0
!     sparse%ishift=0
!     sparse%mpi=.false.
!     !
!   end subroutine mpi_sf_delete_dmatrix_coo
  
!   subroutine mpi_sf_delete_cmatrix_coo(MpiComm,sparse)
!     integer                              :: MpiComm
!     type(sf_sparse_cmatrix_coo),intent(inout) :: sparse
!     integer                              :: i
!     type(sf_sparse_crow_coo),pointer          :: row
!     !
! #ifdef _DEBUG
!     write(*,"(A)")"DEBUG MPI_sf_delete_cmatrix_coo: delete sparse"
! #endif
!     if(.not.sparse%status)return !stop "Error SPARSE/mpi_sp_delete_matrix: sparse is not allocated."
!     !
!     do i=1,sparse%Nrow
!        deallocate(sparse%row(i)%vals)
!        deallocate(sparse%row(i)%cols)
!        sparse%row(i)%Size  = 0
!        !
!        deallocate(sparse%loc(i)%vals)
!        deallocate(sparse%loc(i)%cols)
!        sparse%loc(i)%Size  = 0
!     enddo
!     deallocate(sparse%row)
!     deallocate(sparse%loc)
!     !
!     sparse%Nrow=0
!     sparse%Ncol=0
!     sparse%status=.false.
!     !
!     sparse%istart=0
!     sparse%iend=0
!     sparse%ishift=0
!     sparse%mpi=.false.
!     !
!   end subroutine mpi_sf_delete_cmatrix_coo
! #endif    


  !+------------------------------------------------------------------+
  !PURPOSE: insert an element value at position (i,j) in the sparse matrix
  !+------------------------------------------------------------------+
  subroutine sf_insert_delement_coo(sparse,value,i,j)
    type(sf_sparse_dmatrix_coo),intent(inout) :: sparse
    real(8),intent(in)                    :: value
    integer,intent(in)                    :: i,j
    integer                               :: k
    logical                               :: present
    !
#ifdef _DEBUG
    write(*,"(A,2I8)")"DEBUG sf_insert_delement_coo: insert element in sparse @",i,j
#endif
    !
    present=.false.
    do k=1,sparse%Size !Find position if present
       if( (i==sparse%rows(k)).and.(j==sparse%cols(k)))then
          present=.true.
          exit
       end if
    end do
    !
    if(present)then                            ! Add if present
       sparse%vals(k) = sparse%vals(k) + value !
    else
       call add_to(sparse%rows,i)
       call add_to(sparse%cols,j)
       call add_to(sparse%vals,value)
       sparse%Size = sparse%Size +1
    endif
    !
    if(sparse%Size > sparse%Ncol*sparse%Nrow)stop "sf_insert_delement_coo ERROR: sparse%Size > sparse%Ncol*sparse%Nrow"
    !
  end subroutine sf_insert_delement_coo

  subroutine sf_insert_celement_coo(sparse,value,i,j)
    type(sf_sparse_cmatrix_coo),intent(inout) :: sparse
    complex(8),intent(in)                 :: value
    integer,intent(in)                    :: i,j
    integer                               :: k
    logical                               :: present
    !
#ifdef _DEBUG
    write(*,"(A,2I8)")"DEBUG sf_insert_celement_coo: insert element in sparse @",i,j
#endif
    !
    present=.false.
    do k=1,sparse%Size !Find position if present
       if( (i==sparse%rows(k)).and.(j==sparse%cols(k)))then
          present=.true.
          exit
       end if
    end do
    !
    if(present)then                            ! Add if present
       sparse%vals(k) = sparse%vals(k) + value !
    else
       call add_to(sparse%rows,i)
       call add_to(sparse%cols,j)
       call add_to(sparse%vals,value)
       sparse%Size = sparse%Size +1
    endif
    !
    if(sparse%Size > sparse%Ncol*sparse%Nrow)stop "sf_insert_celement_coo ERROR: sparse%Size > sparse%Ncol*sparse%Nrow"
    !
  end subroutine sf_insert_celement_coo

! #ifdef _MPI
!   subroutine mpi_sf_insert_delement_coo(MpiComm,sparse,value,i,j)
!     integer                               :: MpiComm
!     type(sf_sparse_dmatrix_coo),intent(inout) :: sparse
!     real(8),intent(in)                    :: value
!     integer,intent(in)                    :: i,j
!     type(sf_sparse_drow_coo),pointer          :: row
!     integer                               :: column,pos
!     logical                               :: iadd
!     !
! #ifdef _DEBUG
!     write(*,"(A,2I8)")"DEBUG MPI_sf_insert_delement_coo: insert element in sparse @",i,j
! #endif
!     !
!     if(MpiComm==Mpi_Comm_Null)return
!     !
!     call sp_test_matrix_mpi(MpiComm,sparse," mpi_sf_insert_delement_coo")
!     !
!     column = j
!     !
!     row => sparse%row(i-sparse%Ishift)
!     !
!     iadd = .false.                          !check if column already exist
!     if(any(row%cols == column))then         !
!        pos = binary_search(row%cols,column) !find the position  column in %cols        
!        iadd=.true.                          !set Iadd to true
!     endif
!     !
!     if(iadd)then                            !this column exists so just sum it up       
!        row%vals(pos)=row%vals(pos) + value  !add up value to the current one in %vals
!     else                                    !this column is new. increase counter and store it 
!        call add_to(row%vals,value)
!        call add_to(row%cols,column)
!        row%Size = row%Size + 1
!     endif
!     !
!     if(row%Size > sparse%Ncol)stop "mpi_sp_insert_element_coo ERROR: row%Size > sparse%Ncol"
!     !
!   end subroutine mpi_sf_insert_delement_coo

!   subroutine mpi_sf_insert_celement_coo(MpiComm,sparse,value,i,j)
!     integer                               :: MpiComm
!     type(sf_sparse_cmatrix_coo),intent(inout) :: sparse
!     complex(8),intent(in)                 :: value
!     integer,intent(in)                    :: i,j
!     type(sf_sparse_crow_coo),pointer          :: row
!     integer                               :: column,pos
!     logical                               :: iadd
!     !
! #ifdef _DEBUG
!     write(*,"(A,2I8)")"DEBUG MPI_sf_insert_celement_coo: insert element in sparse @",i,j
! #endif
!     !
!     call sp_test_matrix_mpi(MpiComm,sparse," mpi_sf_insert_celement_coo")
!     !
!     column = j
!     !
!     if(column>=sparse%Istart.AND.column<=sparse%Iend)then
!        row => sparse%loc(i-sparse%Ishift)
!     else
!        row => sparse%row(i-sparse%Ishift)
!     endif
!     !
!     iadd = .false.                          !check if column already exist
!     if(any(row%cols == column))then         !
!        pos = binary_search(row%cols,column) !find the position  column in %cols        
!        iadd=.true.                          !set Iadd to true
!     endif
!     !
!     if(iadd)then                            !this column exists so just sum it up       
!        row%vals(pos)=row%vals(pos) + value  !add up value to the current one in %vals
!     else                                    !this column is new. increase counter and store it 
!        call add_to(row%vals,value)
!        call add_to(row%cols,column)
!        row%Size = row%Size + 1
!     endif
!     !
!     if(row%Size > sparse%Ncol)stop "mpi_sp_insert_element_coo ERROR: row%Size > sparse%Ncol"
!     !
!   end subroutine mpi_sf_insert_celement_coo
! #endif
  
  !+------------------------------------------------------------------+
  !PURPOSE: dump a sparse matrix into a regular 2dim array
  !+------------------------------------------------------------------+
  subroutine sf_dump_dmatrix_coo(sparse,matrix)
    type(sf_sparse_dmatrix_coo),intent(in)   :: sparse
    real(8),dimension(:,:),intent(inout) :: matrix
    real(8)                              :: val
    integer                              :: i,col,row,Ndim1,Ndim2
    !
#ifdef _DEBUG
    write(*,"(A)")"DEBUG sf_dump_dmatrix_coo: dump sparse"
#endif
    !
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)
    !
    if(sparse%Nrow/=Ndim1 .OR. sparse%Ncol/=Ndim2)stop "Warning SPARSE/dump_matrix: dimensions error"
    !
    do i=1,sparse%Size
       row=sparse%rows(i);  col=sparse%cols(i)
       matrix(row,col) = matrix(row,col) + sparse%vals(i)
    enddo
  end subroutine sf_dump_dmatrix_coo

  subroutine sf_dump_cmatrix_coo(sparse,matrix)
    type(sf_sparse_cmatrix_coo),intent(in)      :: sparse
    complex(8),dimension(:,:),intent(inout) :: matrix
    complex(8)                              :: vals
    integer                                 :: i,col,row,Ndim1,Ndim2
    !
#ifdef _DEBUG
    write(*,"(A)")"DEBUG sf_dump_cmatrix_coo: dump sparse"
#endif
    !
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)
    !
    if(sparse%Nrow/=Ndim1 .OR. sparse%Ncol/=Ndim2)stop "Warning SPARSE/dump_matrix: dimensions error"
    !
    do i=1,sparse%Size
       row=sparse%rows(i);  col=sparse%cols(i)
       matrix(row,col) = matrix(row,col) + sparse%vals(i)
    enddo
  end subroutine sf_dump_cmatrix_coo

! #ifdef _MPI
!   subroutine mpi_sf_dump_dmatrix_coo(MpiComm,sparse,matrix)
!     integer                              :: MpiComm
!     type(sf_sparse_dmatrix_coo),intent(in)   :: sparse
!     real(8),dimension(:,:),intent(inout) :: matrix
!     real(8),dimension(:,:),allocatable   :: matrix_tmp
!     integer                              :: i,impi,j,N1_,N2_,Ndim1,Ndim2,Nrow,Ncol
!     !
! #ifdef _DEBUG
!     write(*,"(A)")"DEBUG MPI_sf_dump_dmatrix_coo: dump sparse"
! #endif
!     !
!     call sp_test_matrix_mpi(MpiComm,sparse," mpi_sf_dump_dmatrix_coo")
!     !
!     Ndim1=size(matrix,1)
!     Ndim2=size(matrix,2)
!     !
!     N1_  = sparse%Nrow
!     N2_  = sparse%Ncol
!     Nrow = 0
!     Ncol = 0
!     call MPI_AllReduce(N1_,Nrow,1,MPI_Integer,MPI_SUM,MpiComm,MpiIerr)
!     call MPI_AllReduce(N2_,Ncol,1,MPI_Integer,MPI_MAX,MpiComm,MpiIerr)
!     !
!     if(Nrow>Ndim1 .OR. Ncol>Ndim2)stop "Warning SPARSE/mpi_dump_matrix: dimensions error"
!     !
!     allocate(matrix_tmp(Ndim1,Ndim2)) ; matrix_tmp=0d0
!     do i=sparse%Istart,sparse%Iend
!        impi = i - sparse%Ishift
!        do j=1,sparse%row(impi)%Size
!           matrix_tmp(i,sparse%row(impi)%cols(j))=matrix_tmp(i,sparse%row(impi)%cols(j))+sparse%row(impi)%vals(j)
!        enddo
!     enddo
!     !
!     call AllReduce_MPI(MpiCOmm,Matrix_tmp,Matrix)
!     !
!   end subroutine mpi_sf_dump_dmatrix_coo

!   subroutine mpi_sf_dump_cmatrix_coo(MpiComm,sparse,matrix)
!     integer                                 :: MpiComm
!     type(sf_sparse_cmatrix_coo),intent(in)      :: sparse
!     complex(8),dimension(:,:),intent(inout) :: matrix
!     complex(8),dimension(:,:),allocatable   :: matrix_tmp
!     integer                                 :: i,impi,j,N1_,N2_,Ndim1,Ndim2,Nrow,Ncol
!     !
! #ifdef _DEBUG
!     write(*,"(A)")"DEBUG MPI_sf_dump_cmatrix_coo: dump sparse"
! #endif
!     !
!     call sp_test_matrix_mpi(MpiComm,sparse," mpi_sf_dump_cmatrix_coo")
!     !
!     Ndim1=size(matrix,1)
!     Ndim2=size(matrix,2)
!     !
!     N1_  = sparse%Nrow
!     N2_  = sparse%Ncol
!     Nrow = 0
!     Ncol = 0
!     call MPI_AllReduce(N1_,Nrow,1,MPI_Integer,MPI_SUM,MpiComm,MpiIerr)
!     call MPI_AllReduce(N2_,Ncol,1,MPI_Integer,MPI_MAX,MpiComm,MpiIerr)
!     !
!     if(Nrow>Ndim1 .OR. Ncol>Ndim2)stop "Warning SPARSE/mpi_dump_matrix: dimensions error"
!     !
!     allocate(matrix_tmp(Ndim1,Ndim2)) ; matrix_tmp=zero
!     do i=sparse%Istart,sparse%Iend
!        impi = i - sparse%Ishift
!        !Local part:
!        do j=1,sparse%loc(impi)%Size
!           matrix_tmp(i,sparse%loc(impi)%cols(j))=matrix_tmp(i,sparse%loc(impi)%cols(j))+sparse%loc(impi)%vals(j)
!        enddo
!        !
!        !Non-local part:
!        do j=1,sparse%row(impi)%Size
!           matrix_tmp(i,sparse%row(impi)%cols(j))=matrix_tmp(i,sparse%row(impi)%cols(j))+sparse%row(impi)%cvals(j)
!        enddo
!     enddo
!     !
!     call AllReduce_MPI(MpiCOmm,Matrix_tmp,Matrix)
!     !
!   end subroutine mpi_sf_dump_cmatrix_coo
! #endif 

END MODULE SF_SPARSE_ARRAY_COO
