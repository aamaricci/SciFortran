MODULE SF_SPARSE_COMMON
!:synopsis: Common sparse-matrix parent, metadata, and storage helpers
!The concrete CSR/CSC numerical types extend :code:`sparse_matrix`.
!The parent defines the value-independent lifecycle and communication interface.
!A plain parent object contains metadata only; an allocatable
!:code:`class(sparse_matrix)` object can hold any concrete sparse child.
!
  
#ifdef _MPI
  USE SF_MPI, only: Bcast_MPI, check_MPI
  USE MPI
#endif
  implicit none
  
  interface append
    !Append an integer, double real, or double complex scalar/array without sorting.
    !
     module procedure :: append_I
     module procedure :: append_D
     module procedure :: append_Z
     module procedure :: append_IA, append_DA, append_ZA
  end interface append


  
  ! Common metadata and value-independent interface for all sparse types.
  ! Concrete children override lifecycle/count/communication for their storage.
  ! A plain parent object represents metadata only (no numerical entries).
  !SPARSE MATRIX
  type sparse_matrix
    !Common parent for the four concrete sparse types.
    !Use :code:`class(sparse_matrix)` for runtime selection of a concrete type.
    !Bindings :code:`init`, :code:`free`, :code:`nnz`, :code:`bytes`, and (with MPI)
    !:code:`bcast` are overridden by the children; :code:`shape` is inherited.
    !Numerical insertion and algebra require a concrete declaration or
    !:code:`select type`. A plain parent has no numerical row/column storage.
    !
     integer                                   :: Nrow=0 !Physical row count; zero for a fresh or freed object.
     integer                                   :: Ncol=0 !Physical column count; zero for a fresh or freed object.
     integer                                   :: Nelements=0 !Cached stored-entry count, including exact zeros; maintained by library operations.
     logical                                   :: status=.false. !True for initialized objects, including initialized empty matrices.
   contains
     procedure,pass :: shape     => shape_matrix
     procedure,pass :: init      => init_matrix
     procedure,pass :: free      => free_matrix
     procedure,pass :: nnz       => nnz_matrix
     procedure,pass :: bytes     => bytes_matrix
#ifdef _MPI
     procedure,pass :: bcast     => bcast_matrix
#endif
  end type sparse_matrix

  !RETURN SHAPE OF THE SPARSE MATRIX [Nrow,Ncol]
  interface shape
    !Return physical row/column counts for the parent or any concrete extension.
    !
     module procedure :: shape_matrix
  end interface shape
  
  public :: shape
  
contains

  subroutine init_matrix(sparse,Nrow,Ncol)
    !Initialize parent metadata with nonnegative dimensions and zero stored entries.
    !The optional column count defaults to the row count. Concrete children override
    !this binding to allocate their row/column storage.
    !
    class(sparse_matrix),intent(inout) :: sparse !Parent object; overridden bindings dispatch on its concrete dynamic type.
    integer,intent(in) :: Nrow !Nonnegative physical row count.
    integer,intent(in),optional :: Ncol !Optional physical column count; defaults to Nrow.
    call sparse%free()
    sparse%Nrow=Nrow
    sparse%Ncol=Nrow
    if(present(Ncol))sparse%Ncol=Ncol
    if(min(sparse%Nrow,sparse%Ncol)<0)error stop "sparse init: negative dimension"
    sparse%status=.true.
  end subroutine init_matrix

  subroutine free_matrix(sparse)
    !Reset parent dimensions, cached count, and initialization status.
    !Safe for fresh or already freed objects; children override it to release storage.
    !
    class(sparse_matrix),intent(inout) :: sparse !Parent object; overridden bindings dispatch on its concrete dynamic type.
    sparse%Nrow=0
    sparse%Ncol=0
    sparse%Nelements=0
    sparse%status=.false.
  end subroutine free_matrix

  integer function nnz_matrix(sparse) result(nnz)
    !Return the cached stored-entry count for a metadata-only parent.
    !Concrete children override this binding to count their actual storage.
    !
    class(sparse_matrix),intent(in) :: sparse !Parent object; overridden bindings dispatch on its concrete dynamic type.
    nnz=sparse%Nelements
  end function nnz_matrix

  ! Local payload estimate in decimal kB; excludes descriptors and padding.
  function bytes_matrix(self) result(kbytes)
    !Return the parent metadata payload in decimal kilobytes.
    !Excludes padding and descriptors. Children include their numerical storage.
    !This is a local estimate, available without MPI.
    !
    class(sparse_matrix),intent(in) :: self !Parent object or a concrete extension.
    real(8) :: kbytes
    kbytes=(3d0*storage_size(self%Nrow)+storage_size(self%status))/8000d0
  end function bytes_matrix

#ifdef _MPI
  subroutine bcast_matrix(self,comm,root)
    !Collectively broadcast metadata for a plain parent object.
    !Defaults are MPI_COMM_WORLD and root zero in the supplied communicator.
    !MPI must be initialized. Concrete children override this binding to transfer
    !complete storage; all ranks must already hold the same concrete dynamic type.
    !
    class(sparse_matrix),intent(inout) :: self !Parent object or a concrete extension.
    integer,intent(in),optional :: comm,root !Communicator and source rank; defaults MPI_COMM_WORLD and zero.
    integer :: comm_,root_,ncpu,ierr
    if(.not.check_MPI())error stop "sparse bcast: MPI not initialized"
    comm_=MPI_COMM_WORLD
    if(present(comm))comm_=comm
    root_=0
    if(present(root))root_=root
    call MPI_Comm_size(comm_,ncpu,ierr)
    if(root_<0.or.root_>=ncpu)call MPI_Abort(comm_,1,ierr)
    call Bcast_MPI(comm_,self%status,root=root_)
    call Bcast_MPI(comm_,self%Nrow,root=root_)
    call Bcast_MPI(comm_,self%Ncol,root=root_)
    call Bcast_MPI(comm_,self%Nelements,root=root_)
  end subroutine bcast_matrix
#endif

  
  !+------------------------------------------------------------------+
  !PURPOSE:  Return the shape of a sparse matrix
  !+------------------------------------------------------------------+
   function shape_matrix(self) result(shape)
    !Return physical dimensions :code:`[Nrow,Ncol]`.
    !Accepts the parent and every concrete extension through :code:`class(sparse_matrix)`.
    !Available as :code:`shape(A)` and :code:`A%shape()`.
    !
     class(sparse_matrix),intent(in) :: self
     integer,dimension(2)            :: shape
     shape = [self%Nrow,self%Ncol]
   end function shape_matrix


  


  !+------------------------------------------------------------------+
  !PURPOSE  : Sort an array, gives the new ordering of the label.
  !+------------------------------------------------------------------+
  subroutine sort_array(array,order)
    !Sort an integer array in place and return its permutation.
    !After sorting, each output element equals the original element at the
    !corresponding position of :f:var:`order`. This legacy quicksort uses random
    !pivots and advances the Fortran random-number generator. Sparse insertion
    !does not call this routine or automatically sort stored values.
    !
    integer,dimension(:)                    :: array !Integer array to sort in place.
    integer,dimension(size(array))          :: order !Permutation mapping sorted positions to original positions.
    integer,dimension(size(array))          :: backup
    integer                                 :: i
    forall(i=1:size(array))order(i)=i
    call qsort_sort(array, order,1, size(array))
    do i=1,size(array)
       backup(i)=array(order(i))
    enddo
    array=backup
  contains
    recursive subroutine qsort_sort( array, order, left, right )
      integer, dimension(:) :: array
      integer, dimension(:) :: order
      integer               :: left
      integer               :: right
      integer               :: i
      integer               :: last
      if ( left .ge. right ) return
      call qsort_swap( order, left, qsort_rand(left,right) )
      last = left
      do i = left+1, right
         if ( compare(array(order(i)), array(order(left)) ) .lt. 0 ) then
            last = last + 1
            call qsort_swap( order, last, i )
         endif
      enddo
      call qsort_swap( order, left, last )
      call qsort_sort( array, order, left, last-1 )
      call qsort_sort( array, order, last+1, right )
    end subroutine qsort_sort
    !---------------------------------------------!
    subroutine qsort_swap( order, first, second )
      integer, dimension(:) :: order
      integer               :: first, second
      integer               :: tmp
      tmp           = order(first)
      order(first)  = order(second)
      order(second) = tmp
    end subroutine qsort_swap
    !---------------------------------------------!
    integer function qsort_rand( lower, upper )
      integer               :: lower, upper
      real(8)               :: r
      call random_number(r)
      qsort_rand =  lower + nint(r * (upper-lower))
    end function qsort_rand
    !---------------------------------------------!
    function compare(f,g)
      implicit none
      integer               :: f,g
      integer               :: compare
      if(f<g) then
         compare=-1
      else
         compare=1
      endif
    end function compare
  end subroutine sort_array


  !##################################################################
  !              AUXILIARY COMPUTATIONAL ROUTINES
  !##################################################################
  ! Input indices must be sorted increasingly. Returns zero when absent.
  integer function binary_search(Ain,value) result(pos)
    !Search a strictly increasing integer index array without modifying it.
    !Returns its one-based position, or zero when the value is absent, including
    !for an empty array. No allocation or sorting is performed.
    !
    integer,intent(in) :: Ain(:),value
    integer :: lo,hi,mid
    lo=1
    hi=size(Ain)
    pos=0
    do while(lo<=hi)
      mid=lo+(hi-lo)/2
      if(Ain(mid)==value)then
        pos=mid
        return
      elseif(Ain(mid)<value)then
        lo=mid+1
      else
        hi=mid-1
      endif
    enddo
  end function binary_search


  subroutine append_I(vec,val)
    !Append one scalar to an allocatable vector, preserving existing entries.
    !An unallocated vector is allocated. Storage may be reallocated; intended
    !for dynamic assembly rather than bulk numerical kernels.
    !
    integer,dimension(:),allocatable,intent(inout) :: vec !Allocatable destination vector; existing entries are retained.
    integer,intent(in)                             :: val !Scalar to append.
    integer,dimension(:),allocatable               :: tmp
    integer                                        :: n
    !
    if (allocated(vec)) then
       n = size(vec)
       allocate(tmp(n+1))
       tmp(:n) = vec
       call move_alloc(tmp,vec)
       n = n + 1
    else
       n = 1
       allocate(vec(n))
    end if
    !
    !Put val as last entry:
    vec(n) = val
    !
    if(allocated(tmp))deallocate(tmp)
  end subroutine append_I

  subroutine append_D(vec,val)
    !Append one scalar to an allocatable vector, preserving existing entries.
    !An unallocated vector is allocated. Storage may be reallocated; intended
    !for dynamic assembly rather than bulk numerical kernels.
    !
    real(8),dimension(:),allocatable,intent(inout) :: vec !Allocatable destination vector; existing entries are retained.
    real(8),intent(in)                             :: val !Scalar to append.
    real(8),dimension(:),allocatable               :: tmp
    integer                                        :: n
    !
    if (allocated(vec)) then
       n = size(vec)
       allocate(tmp(n+1))
       tmp(:n) = vec
       call move_alloc(tmp,vec)
       n = n + 1
    else
       n = 1
       allocate(vec(n))
    end if
    !
    !Put val as last entry:
    vec(n) = val
    !
    if(allocated(tmp))deallocate(tmp)
  end subroutine append_D

  subroutine append_Z(vec,val)
    !Append one scalar to an allocatable vector, preserving existing entries.
    !An unallocated vector is allocated. Storage may be reallocated; intended
    !for dynamic assembly rather than bulk numerical kernels.
    !
    complex(8),dimension(:),allocatable,intent(inout) :: vec !Allocatable destination vector; existing entries are retained.
    complex(8),intent(in)                             :: val !Scalar to append.
    complex(8),dimension(:),allocatable               :: tmp
    integer                                           :: n
    !
    if (allocated(vec)) then
       n = size(vec)
       allocate(tmp(n+1))
       tmp(:n) = vec
       call move_alloc(tmp,vec)
       n = n + 1
    else
       n = 1
       allocate(vec(n))
    end if
    !
    !Put val as last entry:
    vec(n) = val
    !
    if(allocated(tmp))deallocate(tmp)
  end subroutine append_Z






  subroutine append_IA(vec,values)
    !Append an array to an allocatable vector, preserving entry order.
    !An unallocated destination is allocated. This helper does not sort indices,
    !combine duplicates, or prune zeros.
    !
    integer,allocatable,intent(inout) :: vec(:) !Allocatable destination vector; existing entries are retained.
    integer,intent(in) :: values(:) !Array to append in its existing order.
    if(allocated(vec))then
      vec=[vec,values]
    else
      vec=values
    endif
  end subroutine append_IA

  subroutine append_DA(vec,values)
    !Append an array to an allocatable vector, preserving entry order.
    !An unallocated destination is allocated. This helper does not sort indices,
    !combine duplicates, or prune zeros.
    !
    real(8),allocatable,intent(inout) :: vec(:) !Allocatable destination vector; existing entries are retained.
    real(8),intent(in) :: values(:) !Array to append in its existing order.
    if(allocated(vec))then
      vec=[vec,values]
    else
      vec=values
    endif
  end subroutine append_DA

  subroutine append_ZA(vec,values)
    !Append an array to an allocatable vector, preserving entry order.
    !An unallocated destination is allocated. This helper does not sort indices,
    !combine duplicates, or prune zeros.
    !
    complex(8),allocatable,intent(inout) :: vec(:) !Allocatable destination vector; existing entries are retained.
    complex(8),intent(in) :: values(:) !Array to append in its existing order.
    if(allocated(vec))then
      vec=[vec,values]
    else
      vec=values
    endif
  end subroutine append_ZA

END MODULE SF_SPARSE_COMMON
