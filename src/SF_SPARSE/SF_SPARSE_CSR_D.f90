MODULE SF_SPARSE_CSR_D  
!:synopsis: Double real sparse matrices stored as dynamic rows
!Concrete implementation of :code:`sparse_dmatrix_csr`, extending the common
!:code:`sparse_matrix` parent. Import :code:`SF_SPARSE` for the combined public API.
!
!.. warning::
!
!   NEW column indices MUST increase strictly within each stored row.
!   :code:`insert` accumulates existing coordinates; :code:`fast_insert` requires
!   unique coordinates. Neither operation sorts indices or reorders values.
!
!Dimensions and coordinates are physical row/column counts and one-based indices.
!Matrices are local objects. Optional MPI collectives transfer complete objects.
!
  USE SF_SPARSE_COMMON, only: sparse_matrix, shape, append, binary_search
  USE SF_IOTOOLS, only: str, free_unit
  USE SF_MISC, only: assert_shape
#ifdef _MPI
  USE SF_MPI
  USE MPI
#endif
  implicit none
  private


  real(8),parameter :: zero=0d0
  ! All loop indices are local; matrices are local objects on each MPI rank.


  ! Dynamic sparse row; indices and values have matching lengths.
  type sparse_row
    !One dynamic sparse row. Indices are strictly increasing and values have matching length.
    !
     sequence
     integer                             :: size=0 !Stored-entry count, including exact zeros.
     integer,dimension(:),allocatable    :: cols !One-based column indices, strictly increasing.
     real(8),dimension(:),allocatable    :: vals !Values in the same order as the indices.
  end type sparse_row


  ! Storage contract: indices in each row MUST be strictly increasing.
  ! insert accumulates existing coordinates; NEW coordinates must be appended in order.
  ! fast_insert appends unique coordinates. Neither routine sorts indices or values.
  ! nnz counts stored entries, including exact zeros; algebra applies no tolerance.
  !SPARSE MATRIX STRUCTURE
  type, extends(sparse_matrix) :: sparse_dmatrix_csr
    !Double real local sparse matrix with inherited metadata and lifecycle interface.
    !Each stored row holds increasing column indices and matching values.
    !The inherited :code:`Nelements` cache is maintained by library operations;
    !direct storage edits must also maintain sizes, ordering, and this cache.
    !Use :code:`class(sparse_matrix)` for value-independent polymorphic operations.
    !
     type(sparse_row),dimension(:),allocatable :: row !Concrete sparse row storage; never distributed across MPI ranks.
   contains
     procedure,pass :: init       => sp_init_matrix_csr_d
     procedure,pass :: free       => sp_free_matrix_csr_d
     procedure,pass :: load       => sp_load_matrix_csr_d
     procedure,pass :: copy       => sp_copy_matrix_csr_d
     procedure,pass :: dump       => sp_dump_matrix_csr_d
     procedure,pass :: fast_insert=> sp_fast_insert_element_csr_d
     procedure,pass :: insert     => sp_insert_element_csr_d
     procedure,pass :: get        => sp_get_element_csr_d
     procedure,pass :: show       => sp_show_matrix_csr_d
     procedure,pass :: display    => sp_display_matrix_csr_d
     procedure,pass :: spy        => sp_spy_matrix_csr_d
     procedure,pass :: as_matrix  => sp_as_matrix_csr_d
     procedure,pass :: dgr        => sp_hconjg_matrix_csr_d
     procedure,pass :: t          => sp_transpose_matrix_csr_d
     procedure,pass :: transpose  => sp_transpose_matrix_csr_d
     procedure,pass :: nnz        => sp_nnz_matrix_csr_d
     procedure,pass :: dot        => sp_matmul_vector_csr_d
     procedure,pass :: write      => sp_write_matrix_csr_d
     procedure,pass :: read       => sp_read_matrix_csr_d
#ifdef _MPI
     procedure,pass :: bcast      => sp_bcast_matrix_csr_d
#endif
     procedure,pass :: bytes      => sp_get_bytes_matrix_csr_d
  end type sparse_dmatrix_csr


  interface sparse_dmatrix_csr
    !Dense-array constructor for :code:`sparse_dmatrix_csr`.
    !
     module procedure :: sp_construct_matrix_csr_d
  end interface sparse_dmatrix_csr

  interface as_csr
    !Construct this CSR numerical type from a dense array with :code:`as_csr(matrix)`.
    !
     module procedure :: sp_construct_matrix_csr_d
  end interface as_csr

  interface as_sparse
    !Dense-to-sparse constructor alias available when importing one storage orientation.
    !
     module procedure :: sp_construct_matrix_csr_d
  end interface as_sparse

  interface as_matrix
    !Return the physical dense representation of a sparse object.
    !
     module procedure :: sp_as_matrix_csr_d
  end interface as_matrix

  interface sparse
    !Dense-to-sparse constructor alias available when importing one storage orientation.
    !
     module procedure :: sp_construct_matrix_csr_d
  end interface sparse

  !EQUALITY with scalar and function (A=B, A=cmplx)
  interface assignment(=)
    !Scalar assignment changes stored values only. Matrix-to-matrix copying uses intrinsic assignment.
    !
     module procedure :: sp_matrix_equal_scalar_csr_d
  end interface assignment(=)

  !ADDITION
  interface operator (+)
    !Add equal-shape sparse matrices by merging sorted indices.
    !
     module procedure :: sp_plus_matrix_csr_d
  end interface operator (+)

  !SUBTRACTION
  interface operator (-)
    !Subtract equal-shape sparse matrices without numerical pruning.
    !
     module procedure :: sp_minus_matrix_csr_d
  end interface operator (-)

  !SCALAR PRODUCT
  interface operator(*)
    !Multiply stored values by an integer, real, or supported complex scalar.
    !
     module procedure :: sp_left_product_matrix_i_csr_d
     module procedure :: sp_left_product_matrix_d_csr_d
     !
     module procedure :: sp_right_product_matrix_i_csr_d
     module procedure :: sp_right_product_matrix_d_csr_d
  end interface operator(*)

  !SCALAR DIVISION
  interface operator(/)
    !Divide stored values by a nonzero integer, real, or supported complex scalar.
    !
     module procedure :: sp_right_division_matrix_i_csr_d
     module procedure :: sp_right_division_matrix_d_csr_d
  end interface operator(/)


  !Matrix-Matrix PRODUCT
  interface operator(.m.)
    !Ordinary sparse matrix multiplication, equivalent to :code:`matmul(A,B)`.
    !
     module procedure :: sp_matmul_matrix_csr_d
  end interface operator(.m.)

  !KRONECKER PRODUCT
  interface operator(.x.)
    !Kronecker product, equivalent to :code:`kron(A,B)`.
    !
     module procedure :: sp_kron_matrix_csr_d
  end interface operator(.x.)

  interface kron
    !Kronecker product preserving orientation and numerical type.
    !
     module procedure :: sp_kron_matrix_csr_d
  end interface kron


  !RETURN SHAPE OF THE SPARSE MATRIX

  !EXTEND TRANSPOSE TO SPARSE MATRIX
  intrinsic :: transpose
  interface transpose
    !Transpose without conjugating values; the output keeps the input storage orientation.
    !
     module procedure :: sp_transpose_matrix_csr_d
  end interface transpose

  ! Hermitian conjugate (ordinary transpose for real values).
  interface hconjg
    !Hermitian conjugate; equivalent to ordinary transpose for real values.
    !
     module procedure :: sp_hconjg_matrix_csr_d
  end interface hconjg

  intrinsic :: matmul
  interface matmul
    !Ordinary sparse matrix-matrix, matrix-vector, and vector-matrix products. Complex inputs are not conjugated.
    !
     module procedure :: sp_matmul_matrix_csr_d
     module procedure :: sp_matmul_vector_csr_d
     module procedure :: sp_left_matmul_vector_csr_d
  end interface matmul

#ifdef _MPI
  interface  AllGather_MPI
    !Collectively replicate rank-one or rank-two arrays of whole matrices from unique owners.
    !
     module procedure ::  sp_array_all_gather_1_csr_d
     module procedure ::  sp_array_all_gather_2_csr_d
  end interface AllGather_MPI
#endif


  public :: sparse_dmatrix_csr
  public :: as_csr
  public :: as_sparse
  public :: as_matrix
  public :: sparse
  public :: assignment(=)
  public :: operator(+)
  public :: operator(-)
  public :: operator(*)
  public :: operator(/)
  public :: operator(.m.)
  public :: operator(.x.)
  public :: kron
  public :: shape
  public :: transpose
  public :: hconjg
  public :: matmul
  public :: sp_eye_csr_d
#ifdef _MPI
  public :: AllGather_MPI
#endif


contains       


  !+------------------------------------------------------------------+
  !PURPOSE:  initialize the sparse matrix list
  !+------------------------------------------------------------------+
  subroutine sp_init_matrix_csr_d(sparse,Nrow,Ncol)
    !Initialize an empty local matrix through :code:`A%init(Nrow,Ncol)`.
    !The optional second dimension defaults to :f:var:`Nrow`. Existing storage is
    !released. Dimensions must be nonnegative; all row/column arrays start empty.
    !
    class(sparse_dmatrix_csr),intent(inout) :: sparse !Sparse object; concrete storage belongs to this MPI rank.
    integer,intent(in)                 :: Nrow !Physical number of rows; must be nonnegative.
    integer,intent(in),optional         :: Ncol !Physical number of columns; defaults to Nrow.
    integer                            :: i

    !
    if(sparse%status)call sparse%free()
    ! 
    if(Nrow<0)error stop "sparse init: negative dimension"
    sparse%Nrow=Nrow
    sparse%Ncol=Nrow
    if(present(Ncol))sparse%Ncol=Ncol
    if(sparse%Ncol<0)error stop "sparse init: negative dimension"
    !
    if(allocated(sparse%row))deallocate(sparse%row)
    allocate(sparse%row(Nrow))
    do i=1,Nrow
       sparse%row(i)%size=0
       allocate(sparse%row(i)%vals(0)) !empty array
       allocate(sparse%row(i)%cols(0)) !empty array
    end do
    sparse%Nelements=0
    sparse%status=.true.
  end subroutine sp_init_matrix_csr_d


  function sp_construct_matrix_csr_d(matrix) result(self)
    !Construct a sparse matrix from a rank-two dense array, preserving its shape.
    !Exact zero entries are omitted. Traversal creates sorted indices automatically.
    !Available through the type-name constructor, the orientation-specific
    !:code:`as_csr` or :code:`as_csc` interface, and local :code:`as_sparse` / :code:`sparse`.
    !
    real(8),dimension(:,:),intent(in)    :: matrix !Dense array with physical row/column shape, or returned dense representation.
    type(sparse_dmatrix_csr)                  :: self !Constructed sparse object, including dimensions, storage, and entry count.

    call self%load(matrix)

    self%Nelements=self%nnz()
  end function sp_construct_matrix_csr_d


  subroutine sp_copy_matrix_csr_d(sparse,b)
    !Deep-copy :f:var:`b` into :f:var:`sparse`, including dimensions and status.
    !An uninitialized source leaves the destination uninitialized. Numerical arrays
    !are copied independently. Ordinary matrix assignment also performs a deep copy.
    !
    class(sparse_dmatrix_csr),intent(inout) :: sparse !Sparse object; concrete storage belongs to this MPI rank.
    type(sparse_dmatrix_csr),intent(in) :: b !Second sparse operand or source matrix.
    if(allocated(b%row))then
      sparse%row=b%row
    else
      if(allocated(sparse%row))deallocate(sparse%row)
    endif
    sparse%Nrow=b%Nrow
    sparse%Ncol=b%Ncol
    sparse%Nelements=b%nnz()
    sparse%status=b%status
  end subroutine sp_copy_matrix_csr_d


  !+------------------------------------------------------------------+
  !PURPOSE: free an entire sparse matrix
  !+------------------------------------------------------------------+
  subroutine sp_free_matrix_csr_d(sparse)
    !Release all concrete sparse storage and reset dimensions, entry count, and status.
    !Safe for fresh objects and repeated calls. Overrides the common parent binding.
    !
    class(sparse_dmatrix_csr),intent(inout) :: sparse !Sparse object; concrete storage belongs to this MPI rank.
    if(allocated(sparse%row))deallocate(sparse%row)
    sparse%Nrow=0
    sparse%Ncol=0
    sparse%Nelements=0
    sparse%status=.false.
  end subroutine sp_free_matrix_csr_d


  !+------------------------------------------------------------------+
  !PURPOSE: load a dense matrix into a sparse one
  !+------------------------------------------------------------------+
  subroutine sp_load_matrix_csr_d(sparse,matrix)
    !Replace this object by the sparse representation of :f:var:`matrix`.
    !The physical shape is preserved and exact zeros are omitted. The construction
    !order satisfies the sorted-index contract without requiring a caller-side sort.
    !
    class(sparse_dmatrix_csr),intent(inout)   :: sparse !Sparse object; concrete storage belongs to this MPI rank.
    real(8),dimension(:,:),intent(in)    :: matrix !Dense array with physical row/column ordering.
    integer                              :: Ndim1,Ndim2
    integer :: i,j
    !
    call sparse%free()
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)
    !
    call sparse%init(Ndim1,Ndim2) !set status=T
    do i=1,Ndim1
       do j=1,Ndim2
          if(matrix(i,j)/=zero)call sp_insert_element_csr_d(sparse,matrix(i,j),i,j)
       enddo
    enddo
  end subroutine sp_load_matrix_csr_d


  !+------------------------------------------------------------------+
  !PURPOSE: dump a sparse matrix into a regular 2dim array
  !+------------------------------------------------------------------+
  subroutine sp_dump_matrix_csr_d(sparse,matrix)
    !Write the dense representation into an existing array of matching shape.
    !The destination is cleared first: this operation replaces its contents rather
    !than accumulating into them. The sparse object is unchanged.
    !
    class(sparse_dmatrix_csr),intent(in)         :: sparse !Sparse object; concrete storage belongs to this MPI rank.
    real(8),dimension(:,:),intent(inout)    :: matrix !Dense array with physical row/column ordering.
    integer                                 :: Ndim1,Ndim2
    integer :: i,j
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)
    call assert_shape(matrix,[sparse%Nrow,sparse%Ncol],"sp_dump_matrix_csr_d","Matrix")
    Matrix = zero
    do i=1,Ndim1
       do j=1,sparse%row(i)%Size
          matrix(i,sparse%row(i)%cols(j)) = matrix(i,sparse%row(i)%cols(j)) + sparse%row(i)%vals(j)
       enddo
    enddo
  end subroutine sp_dump_matrix_csr_d


  !+------------------------------------------------------------------+
  !PURPOSE: dump a sparse matrix into a regular 2dim array
  !+------------------------------------------------------------------+
  function sp_as_matrix_csr_d(sparse) result(matrix)
    !Return a dense rank-two array with physical shape :code:`[Nrow,Ncol]`.
    !Absent coordinates are zero. Available as :code:`A%as_matrix()` and
    !:code:`as_matrix(A)`. An uninitialized object has the default zero-by-zero shape.
    !
    class(sparse_dmatrix_csr),intent(in)               :: sparse
    real(8),dimension(sparse%Nrow,sparse%Ncol)    :: matrix !Dense array with physical row/column shape, or returned dense representation.
    integer :: i,j
    matrix = zero
    if(.not.sparse%status)return
    do i=1,sparse%Nrow
       do j=1,sparse%row(i)%Size
          matrix(i,sparse%row(i)%cols(j)) = sparse%row(i)%vals(j)
       enddo
    enddo
  end function sp_as_matrix_csr_d


  !+------------------------------------------------------------------+
  !PURPOSE: insert an element value at position (i,j) in the sparse matrix
  !+------------------------------------------------------------------+
  subroutine sp_insert_element_csr_d(sparse,value,i,j)
    !Add :f:var:`value` to a physical matrix coordinate, called as
    !:code:`A%insert(value,row,column)`. Existing coordinates are accumulated.
    !
    !.. warning::
    !
    !   NEW indices MUST be appended in strictly increasing order within each
    !   stored row (CSR) or column (CSC). This routine does not sort indices or
    !   reorder values. Revisiting an existing coordinate preserves the ordering.
    !
    !The matrix must be initialized and the physical coordinate must be in bounds.
    !Exact zero values remain stored and contribute to :code:`nnz()` / :code:`Nelements`.
    !
    class(sparse_dmatrix_csr),intent(inout) :: sparse !Sparse object; concrete storage belongs to this MPI rank.
    real(8),intent(in)                 :: value !Value to accumulate (insert) or append (fast_insert).
    integer,intent(in)                 :: i,j !Physical coordinate arguments; call with row, column.
    integer                            :: column,pos
    logical                            :: iadd
    !
    if(.not.sparse%status)error stop "sp_insert_element_csr_d: sparse.status=F"
    if(.not.sparse%status)error stop "sparse insert: uninitialized matrix"
    if(i<1.or.i>sparse%Nrow.or.j<1.or.j>sparse%Ncol)error stop "sparse insert: index out of bounds"
    column = j
    !
    iadd = .false.                          !check if column already exist
    if(any(sparse%row(i)%cols == column))then         !
       pos = binary_search(sparse%row(i)%cols,column) !find the position  column in %cols
       iadd=.true.                          !set Iadd to true
    endif
    !
    if(iadd)then                            !this column exists so just sum it up       
       sparse%row(i)%vals(pos)=sparse%row(i)%vals(pos) + value  !add up value to the current one in %vals
    else                                    !this column is new. increase counter and store it 
       call append(sparse%row(i)%vals,value)
       call append(sparse%row(i)%cols,column)
       sparse%row(i)%Size = sparse%row(i)%Size + 1
       sparse%Nelements=sparse%Nelements+1
    endif
    !
    if(sparse%row(i)%Size > sparse%Ncol)error stop "sp_insert_element_csr_d ERROR: row%Size > sparse%Ncol"
    !
  end subroutine sp_insert_element_csr_d

  ! Insertion contract: NEW column indices MUST increase within each CSR row.
  ! Existing coordinates may be accumulated by insert. No sorting is performed.
  ! fast_insert requires unique coordinates in strictly increasing order. 
  subroutine sp_fast_insert_element_csr_d(sparse,value,i,j)
    !Append a value at a physical coordinate through
    !:code:`A%fast_insert(value,row,column)`, without searching for an existing entry.
    !
    !.. warning::
    !
    !   The caller MUST guarantee unique coordinates and strictly increasing
    !   indices within each stored CSR row or CSC column. No sorting or duplicate
    !   accumulation is performed.
    !
    !The object must be initialized and the coordinate must be in bounds.
    !
    class(sparse_dmatrix_csr),intent(inout) :: sparse !Sparse object; concrete storage belongs to this MPI rank.
    real(8),intent(in)                 :: value !Value to accumulate (insert) or append (fast_insert).
    integer,intent(in)                 :: i,j !Physical coordinate arguments; call with row, column.
    integer                            :: column
    !
    if(.not.sparse%status)error stop "sparse insert: uninitialized matrix"
    if(i<1.or.i>sparse%Nrow.or.j<1.or.j>sparse%Ncol)error stop "sparse insert: index out of bounds"
    column = j
    !
    call append(sparse%row(i)%vals,value)
    call append(sparse%row(i)%cols,column)
    sparse%row(i)%Size = sparse%row(i)%Size + 1
    sparse%Nelements=sparse%Nelements+1
    !
    if(sparse%row(i)%Size > sparse%Ncol)error stop "sp_insert_element_csr_d ERROR: row%Size > sparse%Ncol"
    !
  end subroutine sp_fast_insert_element_csr_d


  !+------------------------------------------------------------------+
  !PURPOSE: get the element value at position (i,j) in the sparse matrix
  !+------------------------------------------------------------------+
  function sp_get_element_csr_d(sparse,i,j) result(val)
    !Return the value at a physical coordinate through :code:`A%get(row,column)`.
    !An absent coordinate returns exact zero. The object must be initialized, the
    !coordinate must be in bounds, and stored indices must satisfy the ordering contract.
    !
    class(sparse_dmatrix_csr),intent(in) :: sparse     !Sparse object; concrete storage belongs to this MPI rank.
    integer,intent(in)              :: i,j
    real(8)                         :: val !Value at the requested physical coordinate; zero when absent.
    integer                         :: pos
    if(.not.sparse%status)error stop "sp_get_element_csr_d: sparse.status=F"
    if(i<1.or.i>sparse%Nrow.or.j<1.or.j>sparse%Ncol)error stop "sparse get: index out of bounds"
    val=zero
    if(.not.any(sparse%row(i)%cols==j))return
    pos=binary_search(sparse%row(i)%cols,j)
    val=sparse%row(i)%vals(pos)
  end function sp_get_element_csr_d


  !+------------------------------------------------------------------+
  !PURPOSE: return the number of non-zero elements
  !+------------------------------------------------------------------+
  function sp_nnz_matrix_csr_d(sparse) result(nnz)
    !Count stored entries from the concrete row/column sizes.
    !Exact zeros are counted. Uninitialized objects return zero. This storage-based
    !count is independent of the inherited cache :code:`Nelements`.
    !
    class(sparse_dmatrix_csr),intent(in) :: sparse    
    integer                         :: nnz !Number of stored entries, including exact zeros.
    integer                         :: i

    nnz=0
    if(.not.sparse%status)return
    do i=1,sparse%Nrow
       nnz=nnz + sparse%row(i)%size
    enddo
  end function sp_nnz_matrix_csr_d


  !##################################################################
  !##################################################################
  !              SHOW  / SPY
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE: show matrix
  !+------------------------------------------------------------------+
  subroutine sp_show_matrix_csr_d(sparse,fmt,file)
    !Print the dense matrix in physical row/column order.
    !The default edit descriptor is :code:`ES16.8`; :f:var:`fmt` may supply another
    !real edit descriptor, also used for each component of complex values.
    !Output goes to unit 6 unless :f:var:`file` is supplied. An uninitialized object
    !produces no output. The object is unchanged.
    !
    class(sparse_dmatrix_csr),intent(in) :: sparse !Sparse object; concrete storage belongs to this MPI rank.
    character(*),intent(in),optional :: fmt,file
    character(:),allocatable :: format
    integer :: unit_,i,j
    if(.not.sparse%status)return
    format="ES16.8"
    if(present(fmt))format=fmt
    unit_=6
    if(present(file))open(newunit=unit_,file=file,status="replace",action="write")
    do i=1,sparse%Nrow
      do j=1,sparse%Ncol
        write(unit_,"("//format//",1X)",advance="no")sparse%get(i,j)
      enddo
      write(unit_,*)
    enddo
    if(present(file))close(unit_)
  end subroutine sp_show_matrix_csr_d

  subroutine sp_display_matrix_csr_d(sparse)
    !Print stored physical row indices, column indices, and values to standard output.
    !Absent coordinates are not printed. An uninitialized object produces no output.
    !
    class(sparse_dmatrix_csr),intent(in) :: sparse !Sparse object; concrete storage belongs to this MPI rank.
    integer :: i,j
    if(.not.sparse%status)return
    do i=1,sparse%Nrow
      do j=1,sparse%row(i)%size
        write(*,*)i,sparse%row(i)%cols(j),sparse%row(i)%vals(j)
      enddo
    enddo
  end subroutine sp_display_matrix_csr_d


  subroutine sp_write_matrix_csr_d(sparse,file,unit)
    !Serialize an initialized matrix as formatted dimensions and index/value records.
    !Supply exactly one of :f:var:`file` or :f:var:`unit`. A named file is replaced and
    !closed by this routine; a caller-owned unit remains open. The payload is
    !orientation-specific and must be read with the same storage and numerical type.
    !
    class(sparse_dmatrix_csr),intent(in) :: sparse !Sparse object; concrete storage belongs to this MPI rank.
    character(len=*),optional :: file !File name; owned and closed by this routine.
    integer,optional          :: unit !Caller-owned formatted unit, left open on return.
    integer                   :: unit_

    integer :: i,j
    !
    if(present(file).eqv.present(unit))error stop "sparse I/O: supply exactly one of file or unit"
    if(present(file))then
      open(newunit=unit_,file=file,status="replace",action="write")
    else
      unit_=unit
    endif
    !
    if(.not.sparse%status)error stop "sparse write: uninitialized matrix"
    write(unit_,*)sparse%Nrow,sparse%Ncol
    do i=1,sparse%Nrow
       write(unit_,*)sparse%row(i)%size
       do j=1,sparse%row(i)%size
          write(unit_,*)sparse%row(i)%cols(j),sparse%row(i)%vals(j)
       enddo
    enddo
    if(present(file))close(unit_)
  end subroutine sp_write_matrix_csr_d


  subroutine sp_read_matrix_csr_d(sparse,file,unit)
    !Replace this object by a matrix read from formatted dimensions and index/value records.
    !Supply exactly one of :f:var:`file` or :f:var:`unit`. A named file must exist and is
    !closed after reading; a caller-owned unit remains open. Stored sizes and index
    !bounds are checked, and indices must be strictly increasing within each record.
    !Use the same orientation and numerical type as the writer.
    !
    class(sparse_dmatrix_csr),intent(inout) :: sparse !Sparse object; concrete storage belongs to this MPI rank.
    character(len=*),optional          :: file !File name; owned and closed by this routine.
    integer,optional                   :: unit !Caller-owned formatted unit, left open on return.
    integer                            :: unit_
    integer                            ::Nrow,Ncol
    integer                            :: irow,i,nsize

    integer,allocatable                :: cols(:)
    real(8),allocatable                :: vals(:)

    !
    if(present(file).eqv.present(unit))error stop "sparse I/O: supply exactly one of file or unit"
    if(present(file))then
      open(newunit=unit_,file=file,status="old",action="read")
    else
      unit_=unit
    endif
    !
    if(sparse%status)call sparse%free()
    read(unit_,*)Nrow,Ncol
    call sparse%init(Nrow,Ncol)
    allocate(vals(Ncol));vals=0
    allocate(cols(Ncol));cols=0
    do irow=1,Nrow
       read(unit_,*)nsize
       if(nsize<0.or.nsize>Ncol)error stop "sparse read: invalid row size"
       if(nsize==0)cycle
       do i=1,nsize
          read(unit_,*)cols(i),vals(i)
       enddo
       if(any(cols(1:nsize)<1).or.any(cols(1:nsize)>Ncol))error stop "sparse read: invalid index"
       do i=2,nsize
         if(cols(i)<=cols(i-1))error stop "sparse read: indices must be strictly increasing"
       enddo
       sparse%row(irow)%size=nsize
       call append(sparse%row(irow)%cols,cols(1:nsize))
       call append(sparse%row(irow)%vals,vals(1:nsize))
    enddo
    if(present(file))close(unit_)
    deallocate(vals,cols)
    sparse%Nelements=sparse%nnz()
  end subroutine sp_read_matrix_csr_d


  !+------------------------------------------------------------------+
  !PURPOSE: pretty print a sparse matrix on a given unit using format fmt
  !+------------------------------------------------------------------+  
  subroutine sp_spy_matrix_csr_d(sparse,header)
    !Write stored coordinates to :code:`header_data.dat` and plotting instructions to
    !:code:`plot_header_commands.gp`. Coordinates are written as column, row for plotting.
    !The script requests PDF output through Gnuplot/PostScript utilities, but this
    !routine does not execute those programs. An uninitialized object produces no files.
    !
    class(sparse_dmatrix_csr)    :: sparse !Sparse object; concrete storage belongs to this MPI rank.
    character ( len = * )   :: header !Prefix for coordinate data and Gnuplot command files.
    integer                 :: N1,N2
    character ( len = 255 ) :: command_filename
    integer                 :: command_unit
    character ( len = 255 ) :: data_filename
    integer                 :: data_unit
    integer                 :: i, j
    integer                 :: nz_num
    !
    !  Create data file.
    !
    if(.not.sparse%status)return
    !
    N1 = sparse%Nrow
    N2 = sparse%Ncol
    data_filename = trim ( header ) // '_data.dat'
    open (unit=free_unit(data_unit), file = data_filename, status = 'replace' )
    nz_num = 0
    do i=1,N1
       do j=1,sparse%row(i)%size
          write(data_unit,'(2x,i6,2x,i6)') sparse%row(i)%cols(j),i
          nz_num = nz_num + 1
       enddo
    enddo
    close(data_unit)
    !
    !  Create command file.
    !
    command_filename = "plot_"//str(header)//'_commands.gp'
    open(unit = free_unit(command_unit), file = command_filename, status = 'replace' )
    write(command_unit,'(a)') '#unset key'
    write(command_unit,'(a)') 'set terminal postscript eps enhanced color font "Times-Roman,16"'
    write(command_unit,'(a)') 'set output "|ps2pdf -sEPSCrop - '//str(header)//".pdf"//'"'
    write(command_unit,'(a)') 'set size ratio -1'
    write(command_unit,'(a)') 'set xlabel "<--- J --->"'
    write(command_unit,'(a)') 'set ylabel "<--- I --->"'
    write(command_unit,'(a,i6,a)')'set title "',nz_num,' nonzeros for '//str(header)//'"'
    write(command_unit,'(a)') 'set timestamp'
    write(command_unit,'(a)' )'plot [x=1:'//str(N1)//'] [y='//str(N2)//':1] "'//&
         str(data_filename)//'" w p pt 5 ps 0.4 lc rgb "red"'
    close ( unit = command_unit )
    return
  end subroutine sp_spy_matrix_csr_d


  !##################################################################
  !##################################################################
  !               SPARSE MATRIX KRON PRODUCT 
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE:  Perform simple Kroenecker product of two sparse matrices
  !AxB(1+rowB*(i-1):rowB*i,1+colB*(j-1):colB*j)  =  A(i,j)*B(:,:)
  !+------------------------------------------------------------------+
  ! Each output row has exactly nnz(A_row)*nnz(B_row) entries.
  function sp_kron_matrix_csr_d(A,B) result(C)
    !Return the Kronecker product, preserving the storage orientation and numerical type.
    !The physical indexing convention is
    !:code:`C((i-1)*B%Nrow+k,(j-1)*B%Ncol+l)=A(i,j)*B(k,l)`.
    !Each output row/column is allocated to its exact structural entry count.
    !Both inputs must be initialized. Result dimensions must fit default integers.
    !Available as :code:`kron(A,B)` and :code:`A .x. B`.
    !
    type(sparse_dmatrix_csr),intent(in) :: A,B
    type(sparse_dmatrix_csr) :: C !Returned result; its physical shape and semantics are specified above.
    integer :: i,k,p,q,r,n,major,minor
    if(.not.A%status.or..not.B%status)error stop "sparse kron: uninitialized matrix"
    if(A%Nrow>huge(1)/max(1,B%Nrow).or.A%Ncol>huge(1)/max(1,B%Ncol)) &
      error stop "sparse kron: result dimensions overflow integer range"
    call C%init(A%Nrow*B%Nrow,A%Ncol*B%Ncol)
    do i=1,A%Nrow
      do k=1,B%Nrow
        major=(i-1)*B%Nrow+k
        n=A%row(i)%size*B%row(k)%size
        deallocate(C%row(major)%cols,C%row(major)%vals)
        allocate(C%row(major)%cols(n),C%row(major)%vals(n))
        C%row(major)%size=n
        r=0
        do p=1,A%row(i)%size
          do q=1,B%row(k)%size
            minor=(A%row(i)%cols(p)-1)*B%Ncol+B%row(k)%cols(q)
            r=r+1
            C%row(major)%cols(r)=minor
            C%row(major)%vals(r)=A%row(i)%vals(p)*B%row(k)%vals(q)
          enddo
        enddo
      enddo
    enddo

    C%Nelements=C%nnz()
  end function sp_kron_matrix_csr_d


  function sp_matmul_matrix_csr_d(A,B) result(AxB)
    !Return the ordinary sparse matrix product :math:`A B`.
    !Both operands must be initialized, have the same concrete type, and satisfy
    !:code:`A%Ncol == B%Nrow`. The result has shape :code:`[A%Nrow,B%Ncol]`.
    !The kernel intersects sorted index streams against a transposed operand;
    !no complex conjugation is applied. Exact-zero results are omitted.
    !Available as :code:`matmul(A,B)` and :code:`A .m. B`.
    !
    type(sparse_dmatrix_csr), intent(in)     :: A,B
    type(sparse_dmatrix_csr)                 :: Bt,AxB !Returned result; its physical shape and semantics are specified above.
    integer                             :: i,icol,j,jcol,k,ii
    integer,dimension(:),allocatable    :: cols
    real(8),dimension(:),allocatable    :: vals
    real(8)                             :: value
    integer, dimension(:), allocatable  :: indx_A, indx_B
    integer                             :: count
    !
    if(.not.A%status)error stop "sp_matmul_matrix_csr_d: A.status=F"
    if(.not.B%status)error stop "sp_matmul_matrix_csr_d: B.status=F"
    !
    ! Both inputs are ordinary local matrices.
    if(A%Ncol/=B%Nrow)error stop "sparse matmul: incompatible dimensions"
    call AxB%free
    call AxB%init(a%Nrow,b%Ncol)
    Bt = B%t()
    !
    allocate(cols(b%Ncol))
    allocate(vals(b%Ncol))
    !
    do i=1,AxB%Nrow
       ii  = 0
       cols= 0
       vals= 0
       do j=1,AxB%Ncol
          ! if(.NOT.check_intersection(A%row(i)%cols, Bt%row(j)%cols))cycle
          call get_intersection(A%row(i)%cols, Bt%row(j)%cols) !count = 0 cycle 
          value = zero
          do k=1,count
             icol = indx_A(k)
             jcol = indx_B(k)
             value    = value + A%row(i)%vals(icol)*Bt%row(j)%vals(jcol)
          enddo
          if(value==zero)cycle
          ii = ii + 1
          vals(ii) = value
          cols(ii) = j
       enddo
       AxB%row(i)%Size = ii
       call append(AxB%row(i)%vals,vals(1:ii))
       call append(AxB%row(i)%cols,cols(1:ii))
    enddo
    call Bt%free
    deallocate(cols,vals)
    AxB%Nelements=AxB%nnz()

  contains
    !
    !
    !Order O(size(A)+size(B)) << O(size(A)*size(B))
    subroutine get_intersection(A, B)
      integer, dimension(:), intent(in) :: A, B
      integer                           :: max_intersections
      integer                           :: i, j !local copies no interference
      !
      max_intersections = min(size(A), size(B))
      if(allocated(indx_A))then
        if(size(indx_A)<max_intersections)deallocate(indx_A,indx_B)
      endif
      if(.not.allocated(indx_A))allocate(indx_A(max_intersections),indx_B(max_intersections))
      !
      i = 1
      j = 1
      count = 0
      do while (i <= size(A) .and. j <= size(B))
         if (A(i) < B(j)) then
            i = i + 1
         else if (A(i) > B(j)) then
            j = j + 1
         else ! A(i) == B(j)
            count = count + 1
            indx_A(count) = i
            indx_B(count) = j
            i = i + 1
            j = j + 1
         end if
      end do
    end subroutine get_intersection
    !
  end function sp_matmul_matrix_csr_d


  function sp_matmul_vector_csr_d(H,v) result(Hv)
    !Return :math:`H v` as an allocated vector of length :code:`H%Nrow`.
    !The input length must be :code:`H%Ncol` and :f:var:`H` must be initialized.
    !No complex conjugation is applied. Available as :code:`matmul(H,v)` and
    !:code:`H%dot(v)`.
    !
    class(sparse_dmatrix_csr),intent(in) :: H !Initialized sparse matrix.
    real(8),intent(in) :: v(:)
    real(8),allocatable :: Hv(:) !Returned product vector with the length specified above.
    integer :: i,k
    if(.not.H%status)error stop "sparse matmul: uninitialized matrix"
    if(size(v)/=H%Ncol)error stop "sparse matmul: vector size mismatch"
    allocate(Hv(H%Nrow))
    Hv=zero
    do i=1,H%Nrow
      do k=1,H%row(i)%size
        Hv(i)=Hv(i)+H%row(i)%vals(k)*v(H%row(i)%cols(k))
      enddo
    enddo
  end function sp_matmul_vector_csr_d


  !##################################################################
  !##################################################################
  !               SPARSE MATRIX BASIC ALGEBRA 
  !##################################################################
  !##################################################################


  function sp_transpose_matrix_csr_d(a) result(c)
    !Return the transpose with shape :code:`[A%Ncol,A%Nrow]` and the same storage type.
    !Values are not conjugated. A count-allocate-fill pass constructs sorted output
    !indices. Available as :code:`transpose(A)`, :code:`A%transpose()`, and :code:`A%t()`.
    !
    class(sparse_dmatrix_csr), intent(in) :: a
    type(sparse_dmatrix_csr)              :: c !Returned result; its physical shape and semantics are specified above.
    integer                          :: i, j, k, col, M, pos
    integer, allocatable             :: nnz_per_row(:), idx(:)
    real(8)                          :: val
    if(.not.a%status)error stop "sp_transpose_matrix_csr_d: A.status=F"
    call c%init(a%Ncol,a%Nrow)       !tranpose
    !
    !Count how many elements are per row with a given col-index
    allocate(nnz_per_row(c%Nrow)) ! c%Nrow == a%Ncol
    nnz_per_row = 0
    do i=1,a%Nrow
       do j=1,a%row(i)%size
          col = a%row(i)%cols(j)
          nnz_per_row(col) = nnz_per_row(col) + 1
       enddo
    enddo
    !Allocate C manually to correct size
    do k = 1, c%Nrow
       M = nnz_per_row(k)
       c%row(k)%size = M
       if(allocated(c%row(k)%vals)) deallocate(c%row(k)%vals)
       if(allocated(c%row(k)%cols)) deallocate(c%row(k)%cols)
       allocate(c%row(k)%vals(M))
       allocate(c%row(k)%cols(M))
    enddo
    !Build C using a counter to get the incremental col position:
    allocate(idx(c%Nrow))
    idx = 0
    do i=1,a%Nrow
       do j=1,a%row(i)%size
          k      = a%row(i)%cols(j)  ! A-col => new C-row index
          val    = a%row(i)%vals(j)
          idx(k) = idx(k) + 1      !advance counter 
          pos    = idx(k)
          c%row(k)%cols(pos) = i
          c%row(k)%vals(pos) = val
       enddo
    enddo
    deallocate(nnz_per_row, idx)

    c%Nelements=c%nnz()
  end function sp_transpose_matrix_csr_d


  function sp_hconjg_matrix_csr_d(a) result(c)
    !Return the Hermitian conjugate with shape :code:`[A%Ncol,A%Nrow]`.
    !Complex values are conjugated; for real matrices this is the ordinary transpose.
    !The count-allocate-fill construction preserves sorted output indices.
    !Available as :code:`hconjg(A)` and :code:`A%dgr()`.
    !
    class(sparse_dmatrix_csr), intent(in) :: a
    type(sparse_dmatrix_csr)              :: c !Returned result; its physical shape and semantics are specified above.
    integer                          :: i, j, k, col, M, pos
    integer, allocatable             :: nnz_per_row(:), idx(:)
    real(8)                          :: val
    if(.not.a%status)error stop "sp_transpose_matrix_csr_d: A.status=F"
    call c%init(a%Ncol,a%Nrow)       !tranpose
    !
    !Count how many elements are per row with a given col-index
    allocate(nnz_per_row(c%Nrow)) ! c%Nrow == a%Ncol
    nnz_per_row = 0
    do i=1,a%Nrow
       do j=1,a%row(i)%size
          col = a%row(i)%cols(j)
          nnz_per_row(col) = nnz_per_row(col) + 1
       enddo
    enddo
    !Allocate C manually to correct size
    do k = 1, c%Nrow
       M = nnz_per_row(k)
       c%row(k)%size = M
       if(allocated(c%row(k)%vals)) deallocate(c%row(k)%vals)
       if(allocated(c%row(k)%cols)) deallocate(c%row(k)%cols)
       allocate(c%row(k)%vals(M))
       allocate(c%row(k)%cols(M))
    enddo
    !Build C using a counter to get the incremental col position:
    allocate(idx(c%Nrow))
    idx = 0
    do i=1,a%Nrow
       do j=1,a%row(i)%size
          k      = a%row(i)%cols(j)  ! A-col => new C-row index
          val    = a%row(i)%vals(j)
          idx(k) = idx(k) + 1      !advance counter 
          pos    = idx(k)
          c%row(k)%cols(pos) = i
          c%row(k)%vals(pos) = val
       enddo
    enddo
    deallocate(nnz_per_row, idx)

    c%Nelements=c%nnz()
  end function sp_hconjg_matrix_csr_d


  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix equality spA = spB. Deep copy
  !+------------------------------------------------------------------+


  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix scalar equality spA = const. 
  !+------------------------------------------------------------------+
  subroutine sp_matrix_equal_scalar_csr_d(a,c)
    !Assign a scalar to every stored value without changing the sparsity pattern.
    !The object must be initialized. Absent coordinates remain absent, even for a
    !nonzero scalar; assigning zero retains stored entries and their count.
    !
    type(sparse_dmatrix_csr),intent(inout) :: a !First sparse operand or source matrix.
    real(8),intent(in)                :: c !Scalar value; a divisor must be nonzero.
    integer :: i
    if(.not.a%status)error stop "sp_matrix_equal_matrix: A.status=F"
    do i=1,a%Nrow
       a%row(i)%vals = c
    enddo
  end subroutine sp_matrix_equal_scalar_csr_d


  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix addition spA + spB = spC
  !+------------------------------------------------------------------+
  function sp_plus_matrix_csr_d(a,b) result(c)
    !Return the sum of two initialized matrices of equal shape and concrete type.
    !A linear merge of sorted index streams constructs the result in the same
    !orientation. Exact cancellations remain stored; no numerical cutoff is applied.
    !
    type(sparse_dmatrix_csr), intent(in) :: a,b
    type(sparse_dmatrix_csr)             :: c !Returned result; its physical shape and semantics are specified above.
    integer                         :: i, na, nb, p, q, nn, col_a, col_b
    integer,allocatable             :: cols(:)
    real(8), allocatable            :: vals(:)
    real(8)                         :: va, vb

    if(.not.a%status) error stop "sp_plus_matrix_merge error: a.status=F"
    if(.not.b%status) error stop "sp_plus_matrix_merge error: b.status=F"
    if(a%Nrow /= b%Nrow) error stop "sp_plus_matrix_merge error: Nrow mismatch"
    if(a%Ncol /= b%Ncol) error stop "sp_plus_matrix_merge error: Ncol mismatch"
    !
    call c%free()
    call c%init(a%Nrow,a%Ncol)
    !
    do i=1,a%Nrow
       na = a%row(i)%size
       nb = b%row(i)%size
       if(na+nb == 0) cycle
       ! allocate worst-case storage
       allocate(cols(na+nb))
       allocate(vals(na+nb))
       p = 1
       q = 1
       nn = 0
       do while(p <= na .or. q <= nb)
          if(p > na) then
             ! copy remaining from b
             nn = nn + 1
             cols(nn) = b%row(i)%cols(q)
             vals(nn) = b%row(i)%vals(q)
             q = q + 1
          elseif(q > nb) then
             nn = nn + 1
             cols(nn) = a%row(i)%cols(p)
             vals(nn) = a%row(i)%vals(p)
             p = p + 1
          else
             col_a = a%row(i)%cols(p)
             col_b = b%row(i)%cols(q)
             if(col_a == col_b) then
                va = a%row(i)%vals(p)
                vb = b%row(i)%vals(q)
                nn = nn + 1
                cols(nn) = col_a
                vals(nn) = va + vb
                p = p + 1; q = q + 1
             elseif(col_a < col_b) then
                nn = nn + 1
                cols(nn) = col_a
                vals(nn) = a%row(i)%vals(p)
                p = p + 1
             else
                nn = nn + 1
                cols(nn) = col_b
                vals(nn) = b%row(i)%vals(q)
                q = q + 1
             endif
          endif
       end do
       !
       c%row(i)%size = nn
       call append(c%row(i)%cols,cols(1:nn))
       call append(c%row(i)%vals,vals(1:nn))
       deallocate(cols,vals)
    end do

    c%Nelements=c%nnz()
  end function sp_plus_matrix_csr_d


  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix difference spA - spB = spC
  !+------------------------------------------------------------------+
  function sp_minus_matrix_csr_d(a,b) result(c)
    !Return the difference of two initialized matrices of equal shape and concrete type.
    !The result uses the sorted merge of the first operand with the negated second
    !operand. Exact cancellations remain stored; no numerical cutoff is applied.
    !
    type(sparse_dmatrix_csr),intent(in) :: a,b
    type(sparse_dmatrix_csr) :: c !Returned result; its physical shape and semantics are specified above.
    c=sp_plus_matrix_csr_d(a,sp_left_product_matrix_i_csr_d(-1,b))

    c%Nelements=c%nnz()
  end function sp_minus_matrix_csr_d


  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix left scalar product const*spA = spC.
  ! type[Const]=integer(4),real(8),cmplx(8)
  !+------------------------------------------------------------------+
  function sp_left_product_matrix_i_csr_d(C,A) result(B)
    !Return the scalar times the sparse matrix, preserving shape, orientation, and stored entries.
    !The matrix must be initialized. Integer and real scalars are supported; complex
    !scalars are also supported by the complex implementation.
    !
    integer,intent(in)             :: C !Scalar value; a divisor must be nonzero.
    type(sparse_dmatrix_csr),intent(in) :: A
    type(sparse_dmatrix_csr)            :: B !Returned result; its physical shape and semantics are specified above.


    integer :: i
    if(.not.A%status)error stop "sp_left_product_matrix error: A.status=F"
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do i=1,a%Nrow
       b%row(i)%size=a%row(i)%size
       call append(b%row(i)%cols,a%row(i)%cols)
       call append(b%row(i)%vals,a%row(i)%vals*C)
    enddo

    B%Nelements=B%nnz()
  end function sp_left_product_matrix_i_csr_d

  function sp_left_product_matrix_d_csr_d(C,A) result(B)
    !Return the scalar times the sparse matrix, preserving shape, orientation, and stored entries.
    !The matrix must be initialized. Integer and real scalars are supported; complex
    !scalars are also supported by the complex implementation.
    !
    real(8),intent(in)             :: C !Scalar value; a divisor must be nonzero.
    type(sparse_dmatrix_csr),intent(in) :: A
    type(sparse_dmatrix_csr)            :: B !Returned result; its physical shape and semantics are specified above.


    integer :: i
    if(.not.A%status)error stop "sp_left_product_matrix error: A.status=F"
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do i=1,a%Nrow
       b%row(i)%size=a%row(i)%size
       call append(b%row(i)%cols,a%row(i)%cols)
       call append(b%row(i)%vals,a%row(i)%vals*C)
    enddo

    B%Nelements=B%nnz()
  end function sp_left_product_matrix_d_csr_d


  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix right scalar product spA*const = spC.
  ! type[Const]=integer(4),real(8),cmplx(8)
  !+------------------------------------------------------------------+
  function sp_right_product_matrix_i_csr_d(A,C) result(B)
    !Return the sparse matrix times a scalar, preserving shape, orientation, and stored entries.
    !The matrix must be initialized. This scaling does not remove exact zeros.
    !
    integer,intent(in)             :: C
    type(sparse_dmatrix_csr),intent(in) :: A !First sparse operand or source matrix.
    type(sparse_dmatrix_csr)            :: B !Returned result; its physical shape and semantics are specified above.


    integer :: i
    if(.not.A%status)error stop "sp_right_product_matrix error: A.status=F"
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do i=1,a%Nrow
       b%row(i)%size=a%row(i)%size
       call append(b%row(i)%cols,a%row(i)%cols)
       call append(b%row(i)%vals,a%row(i)%vals*C)
    enddo

    B%Nelements=B%nnz()
  end function sp_right_product_matrix_i_csr_d

  function sp_right_product_matrix_d_csr_d(A,C) result(B)
    !Return the sparse matrix times a scalar, preserving shape, orientation, and stored entries.
    !The matrix must be initialized. This scaling does not remove exact zeros.
    !
    real(8),intent(in)             :: C
    type(sparse_dmatrix_csr),intent(in) :: A !First sparse operand or source matrix.
    type(sparse_dmatrix_csr)            :: B !Returned result; its physical shape and semantics are specified above.


    integer :: i
    if(.not.A%status)error stop "sp_right_product_matrix error: A.status=F"
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do i=1,a%Nrow
       b%row(i)%size=a%row(i)%size
       call append(b%row(i)%cols,a%row(i)%cols)
       call append(b%row(i)%vals,a%row(i)%vals*C)
    enddo

    B%Nelements=B%nnz()
  end function sp_right_product_matrix_d_csr_d


  !+------------------------------------------------------------------+
  !PURPOSE:  Sparse matrix right scalar division spA/const = spC.
  ! type[Const]=integer(4),real(8),cmplx(8)
  !+------------------------------------------------------------------+
  function sp_right_division_matrix_i_csr_d(A,C) result(B)
    !Return the sparse matrix divided by a nonzero scalar.
    !The matrix must be initialized. Shape, orientation, and the stored sparsity
    !pattern are preserved; zero divisors must be excluded by the caller.
    !
    integer,intent(in)             :: C
    type(sparse_dmatrix_csr),intent(in) :: A !First sparse operand or source matrix.
    type(sparse_dmatrix_csr)            :: B !Returned result; its physical shape and semantics are specified above.


    integer :: i
    if(.not.A%status)error stop "sp_right_division_matrix error: A.status=F"
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do i=1,a%Nrow
       b%row(i)%size=a%row(i)%size
       call append(b%row(i)%cols,a%row(i)%cols)
       call append(b%row(i)%vals,a%row(i)%vals/C)
    enddo

    B%Nelements=B%nnz()
  end function sp_right_division_matrix_i_csr_d

  function sp_right_division_matrix_d_csr_d(A,C) result(B)
    !Return the sparse matrix divided by a nonzero scalar.
    !The matrix must be initialized. Shape, orientation, and the stored sparsity
    !pattern are preserved; zero divisors must be excluded by the caller.
    !
    real(8),intent(in)             :: C
    type(sparse_dmatrix_csr),intent(in) :: A !First sparse operand or source matrix.
    type(sparse_dmatrix_csr)            :: B !Returned result; its physical shape and semantics are specified above.


    integer :: i
    if(.not.A%status)error stop "sp_right_division_matrix error: A.status=F"
    call b%free()
    call b%init(a%Nrow,a%Ncol)
    do i=1,a%Nrow
       b%row(i)%size=a%row(i)%size
       call append(b%row(i)%cols,a%row(i)%cols)
       call append(b%row(i)%vals,a%row(i)%vals/C)
    enddo

    B%Nelements=B%nnz()
  end function sp_right_division_matrix_d_csr_d


  !+------------------------------------------------------------------+
  !PURPOSE:  Return the shape of a sparse matrix
  !+------------------------------------------------------------------+


  !+------------------------------------------------------------------+
  !PURPOSE:  Return the identiy sparse matrix of given dimension
  !+------------------------------------------------------------------+
  function sp_eye_csr_d(ndim) result(self)
    !Construct an initialized sparse identity of size :f:var:`ndim` in this concrete type.
    !The diagonal is assembled directly without a dense temporary. The dimension
    !must be nonnegative. Use the explicitly named identity constructor to select
    !orientation and numerical type.
    !
    type(sparse_dmatrix_csr) :: self !Constructed sparse object, including dimensions, storage, and entry count.
    integer,intent(in) :: ndim !Nonnegative number of identity rows and columns.
    integer :: i
    call self%init(ndim,ndim)
    do i=1,ndim
      call self%fast_insert(1d0,i,i)
    enddo

    self%Nelements=self%nnz()
  end function sp_eye_csr_d


#ifdef _MPI
  subroutine sp_bcast_matrix_csr_d(self,comm,root)
    !Collectively broadcast a complete local sparse matrix from :f:var:`root`.
    !Defaults are :code:`MPI_COMM_WORLD` and root zero. The root is a rank in the
    !supplied communicator. Dimensions, stored entries, cached count, and initialization
    !status are copied; existing receiver contents are replaced. An initialized empty
    !matrix remains distinct from an uninitialized object.
    !
    !Every participant must hold the same concrete dynamic type, including when
    !calling through :code:`class(sparse_matrix)`. This binding does not allocate or
    !change the dynamic type of the receiver. MPI must already be initialized.
    !
    class(sparse_dmatrix_csr), intent(inout) :: self !Sparse object or constructed result, including inherited metadata.
    integer,intent(in),optional         :: comm !MPI communicator; broadcast defaults to MPI_COMM_WORLD.
    integer,intent(in),optional         :: root !Source rank within comm; default zero.
    integer                             :: comm_
    integer                             :: root_
    integer                             :: rank, ierr, i, ncpu
    integer                             :: Nrow, Ncol, Nsize
    logical                             :: master, initialized

    !
    if(.not.check_MPI())error stop "sp_bcast error: check_MPI=F"
    comm_  = MPI_COMM_WORLD;if(present(comm))comm_ = comm
    root_  = 0 ;if(present(root))root_=root
    rank   = get_Rank_MPI(comm_)
    ncpu=get_Size_MPI(comm_)
    if(root_<0.or.root_>=ncpu)call MPI_Abort(comm_,1,ierr)
    master = (rank == root_)
    !
    if(master)initialized=self%status
    call Bcast_MPI(comm_,initialized,root=root_)
    if(.not.initialized)then
      call self%free()
      return
    endif
    if(master)Nrow = self%Nrow ; call Bcast_MPI(comm_,Nrow, root=root_)
    if(master)Ncol = self%Ncol ; call Bcast_MPI(comm_,Ncol, root=root_)
    !
    if(.not.master)call self%init(Nrow, Ncol)
    !
    do i=1,Nrow
       if(master)Nsize = self%row(i)%size
       call Bcast_MPI(comm_,Nsize, root=root_)
       !
       self%row(i)%size = Nsize !tautology for master
       if(Nsize==0)cycle
       if (.NOT.master) then
          if(allocated(self%row(i)%cols)) deallocate(self%row(i)%cols)
          if(allocated(self%row(i)%vals)) deallocate(self%row(i)%vals)
          allocate(self%row(i)%cols(Nsize))
          allocate(self%row(i)%vals(Nsize))
       endif
       !
       call Bcast_MPI(comm_,self%row(i)%cols, root=root_)       
       call Bcast_MPI(comm_,self%row(i)%vals, root=root_)
    end do
    self%Nelements=self%nnz()
  end subroutine sp_bcast_matrix_csr_d


  !+------------------------------------------------------------------+
  ! PURPOSE:  Performs an "All-Gather" operation on a distributed
  !           array of sparse_dmatrix_csr objects. After the call, every
  !           process will have the complete, assembled array.
  ! USAGE:    Assumes the array is distributed, where for each index `i`,
  !           the matrix A(i) is initialized (i.e., A(i)%status == .true.)
  !           on exactly ONE process.
  !+------------------------------------------------------------------+
  subroutine sp_array_all_gather_1_csr_d(comm,A)
    !Collectively assemble a rank-one array of whole local sparse matrices.
    !Array extents must agree across :f:var:`comm`. Initialization status identifies
    !ownership: each entry may have at most one owner. Unowned entries remain
    !uninitialized. Multiple owners or inconsistent extents cause :code:`MPI_Abort`.
    !On return all participants hold replicas; free replicas before another gather
    !that uses initialization status to express ownership. Matrix rows/columns are
    !not distributed objects.
    !
    class(sparse_dmatrix_csr),dimension(:),intent(inout) :: A !Rank-one or rank-two array of matrices; initialization marks unique ownership.
    integer,intent(in)                              :: comm !MPI communicator; broadcast defaults to MPI_COMM_WORLD.
    integer,allocatable,dimension(:)                :: local_ownership, global_ownership
    integer                                         :: i, N, rank, owner_rank, ierr
    logical                                         :: master

    integer :: dims(1), dmin(1), dmax(1)
    integer,allocatable :: owners(:),counts(:)
    !
    if (.not. check_MPI()) error stop "sp_array_all_gather: check_MPI=F"
    !
    dims=shape(A)
    call MPI_Allreduce(dims,dmin,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(dims,dmax,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(any(dmin/=dmax))call MPI_Abort(comm,2,ierr)
    N = size(A)
    rank   = get_Rank_MPI(comm)
    master = get_Master_MPI(comm)
    !
    ! Step 1: Determine the owner of each array element.
    allocate(local_ownership(N))
    allocate(global_ownership(N))
    !
    ! Each process creates a local map. If it owns A(i), it marks it with its rank.
    ! Otherwise, it marks it with -1 (or any value lower than a valid rank).
    do i=1,N
       if (A(i)%status) then
          local_ownership(i) = rank
       else
          local_ownership(i) = -1
       endif
    enddo
    !
    ! Use MPI_Allreduce with MPI_MAX. For each index `i`, the process with the
    ! valid rank will "win", and every process will get the correct owner rank.
    call MPI_Allreduce(local_ownership, global_ownership, N, MPI_INTEGER, MPI_MAX, comm, ierr)
    !
    allocate(owners(N),counts(N))
    owners=merge(1,0,local_ownership>=0)
    call MPI_Allreduce(owners,counts,size(A),MPI_INTEGER,MPI_SUM,comm,ierr)
    if(any(counts>1))call MPI_Abort(comm,3,ierr)
    ! Step 2: Iterate through the array and broadcast each element from its owner.
    do i = 1, N
       owner_rank = global_ownership(i)
       !
       ! If owner_rank is -1, it means NO process initialized this element.
       ! Skip the broadcast, leave A(i) uninitialized state on all processes.
       if (owner_rank == -1) then
          cycle
       else
          !Owner exist, bcast to all nodes
          call A(i)%bcast(comm,root=owner_rank)
       endif
       !
    enddo
    !
    deallocate(global_ownership)
    deallocate(local_ownership)
    call Barrier_MPI(comm)
  end subroutine sp_array_all_gather_1_csr_d


  subroutine sp_array_all_gather_2_csr_d(comm,A)
    !Collectively assemble a rank-two array of whole local sparse matrices.
    !Both array extents must agree across :f:var:`comm`. Each entry may be initialized
    !on at most one rank; unowned entries remain uninitialized. Multiple owners or
    !inconsistent extents cause :code:`MPI_Abort`. On return every rank holds the
    !complete array of replicas. The matrices themselves use ordinary local storage.
    !
    class(sparse_dmatrix_csr),dimension(:,:),intent(inout) :: A !Rank-one or rank-two array of matrices; initialization marks unique ownership.
    integer,intent(in)                                :: comm !MPI communicator; broadcast defaults to MPI_COMM_WORLD.
    integer,allocatable,dimension(:,:)                :: local_ownership,global_ownership
    integer                                           :: i,j,N,M,rank,owner_rank,ierr
    logical                                           :: master
    integer :: dims(2), dmin(2), dmax(2)
    integer,allocatable :: owners(:,:),counts(:,:)
    !
    if (.not. check_MPI()) error stop "sp_array_all_gather: check_MPI=F"
    !
    dims=shape(A)
    call MPI_Allreduce(dims,dmin,2,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(dims,dmax,2,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(any(dmin/=dmax))call MPI_Abort(comm,2,ierr)
    N = size(A,1)
    M = size(A,2)
    !
    rank   = get_Rank_MPI(comm)
    master = get_Master_MPI(comm)
    !
    ! Step 1: Determine the owner of each array element.
    allocate(local_ownership(N,M))
    allocate(global_ownership(N,M))
    global_ownership=0
    ! Each process creates a local map. If it owns A(i), it marks it with its rank.
    ! Otherwise, it marks it with -1 (or any value lower than a valid rank).
    do i=1,N
       do j=1,M
          if (A(i,j)%status) then
             local_ownership(i,j) = rank
          else
             local_ownership(i,j) = -1
          endif
       enddo
    enddo
    !
    ! Use MPI_Allreduce with MPI_MAX. For each index `i`, the process with the
    ! valid rank will "win", and every process will get the correct owner rank.
    call MPI_Allreduce(local_ownership, global_ownership, size(local_ownership), MPI_INTEGER, MPI_MAX, comm, ierr)
    !
    allocate(owners(N,M),counts(N,M))
    owners=merge(1,0,local_ownership>=0)
    call MPI_Allreduce(owners,counts,size(A),MPI_INTEGER,MPI_SUM,comm,ierr)
    if(any(counts>1))call MPI_Abort(comm,3,ierr)
    ! Step 2: Iterate through the array and broadcast each element from its owner.
    do i=1,N
       do j=1,M
          owner_rank = global_ownership(i,j)
          ! If owner_rank is -1, it means NO process initialized this element.
          ! Skip the broadcast, leave A(i) uninitialized state on all processes.
          if (owner_rank == -1) then
             cycle
          else
             !Owner exist, bcast to all nodes
             call A(i,j)%bcast(comm,root=owner_rank)
          endif
       enddo
    enddo
    !
    deallocate(global_ownership)
    deallocate(local_ownership)
    call Barrier_MPI(comm)
  end subroutine sp_array_all_gather_2_csr_d


#endif
  ! Decimal kB of stored payload, excluding allocatable descriptors and padding.
  function sp_get_bytes_matrix_csr_d(self) result(kbytes)
    !Estimate the local stored payload in decimal kilobytes (:math:`1\,kB=1000` bytes).
    !Includes parent metadata, row/column size integers, stored indices, and numerical
    !values. Excludes allocatable descriptors, allocation overhead, and padding.
    !Available as :code:`A%bytes()` with or without MPI; no collective is performed.
    !
    class(sparse_dmatrix_csr),intent(in) :: self !Constructed sparse object, including dimensions, storage, and entry count.
    real(8) :: kbytes !Estimated local payload in decimal kilobytes.
    kbytes=(3d0*storage_size(self%Nrow)/8d0 + storage_size(self%status)/8d0 + &
      self%Nrow*storage_size(self%Nrow)/8d0 + &
      dble(self%nnz())*(storage_size(self%Nrow)/8d0+storage_size(zero)/8d0))/1000d0
  end function sp_get_bytes_matrix_csr_d


  function sp_left_matmul_vector_csr_d(vector,matrix) result(Hv)
    !Return the ordinary left product :math:`v A`, with output length :code:`matrix%Ncol`.
    !The input vector length must be :code:`matrix%Nrow` and the matrix must be
    !initialized. Complex vector values are NOT conjugated. Available as
    !:code:`matmul(vector,matrix)`.
    !
    class(sparse_dmatrix_csr),intent(in) :: matrix !Dense array with physical row/column shape, or returned dense representation.
    real(8),intent(in) :: vector(:) !Input left vector; its length must match the matrix row count.
    real(8),allocatable :: Hv(:) !Returned product vector with the length specified above.
    integer :: i,k
    if(.not.matrix%status)error stop "sparse matmul: uninitialized matrix"
    if(size(vector)/=matrix%Nrow)error stop "sparse matmul: vector size mismatch"
    allocate(Hv(matrix%Ncol))
    Hv=zero
    do i=1,matrix%Nrow
      do k=1,matrix%row(i)%size
        Hv(matrix%row(i)%cols(k))=Hv(matrix%row(i)%cols(k))+vector(i)*matrix%row(i)%vals(k)
      enddo
    enddo
  end function sp_left_matmul_vector_csr_d

end module SF_SPARSE_CSR_D


