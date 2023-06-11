module SF_SP_LINALG
   USE SF_MISC,   only: assert_shape
   USE SF_RANDOM, only: mt_random
   USE SF_LINALG, only: eye,eigh
#ifdef _MPI
   USE SF_MPI
#endif
   implicit none
#ifdef _MPI
   include 'mpif.h'
#endif
   private


   interface sp_eigh
      module procedure :: lanczos_arpack_d
      module procedure :: lanczos_arpack_c
#ifdef _MPI
      module procedure :: lanczos_parpack_d
      module procedure :: lanczos_parpack_c
#endif
   end interface sp_eigh


   interface sp_lanc_eigh
      module procedure :: lanczos_eigh_d
      module procedure :: lanczos_eigh_c
#ifdef _MPI
      module procedure :: mpi_lanczos_eigh_d
      module procedure :: mpi_lanczos_eigh_c
#endif
   end interface sp_lanc_eigh


   interface sp_lanc_tridiag
      module procedure :: lanczos_tridiag_d
      module procedure :: lanczos_tridiag_c
#ifdef _MPI
      module procedure :: mpi_lanczos_tridiag_d
      module procedure :: mpi_lanczos_tridiag_c
#endif
   end interface sp_lanc_tridiag


   interface sp_dvdson_eigh
      module procedure :: dvdson_eigh_d
   end interface sp_dvdson_eigh


   complex(8),parameter              :: zero=(0d0,0d0)
   complex(8),parameter              :: one=(1d0,0d0)
   complex(8),parameter              :: xi=(0d0,1d0)
   integer,allocatable               :: seed_random(:)
   integer                           :: nrandom
   logical                           :: verb=.false.
   real(8)                           :: threshold_=1.d-12
   integer                           :: ncheck_=10



   !****************************************************************************************
   !                                      PUBLIC
   !****************************************************************************************
   public :: sp_eigh
   public :: sp_lanc_eigh
   public :: sp_lanc_tridiag
   public :: sp_dvdson_eigh
   !****************************************************************************************




contains


   !##################################################################
   ! ARPACK METHOD for LOWEST part of the spectrum of a of a SPARSE
   !   MATRIX (defined via H*v)
   ! - DBLE and CMPLX versions included
   ! - SERIAL and PARALLEL-MPI versions included
   !        [COMM    MPI  Communicator for the processor grid.  (INPUT)]
   include "arpack_d.f90"
   include "arpack_c.f90"

#ifdef _MPI
   subroutine lanczos_parpack_d(MpiComm,MatVec,eval,evec,Nblock,Nitermax,v0,tol,iverbose,vrandom)
      !Arguments
      integer                   :: MpiComm
      !Interface to Matrix-Vector routine:
      interface
         subroutine MatVec(nchunk,vin,vout)
            integer                   :: nchunk
            real(8),dimension(nchunk) :: vin,vout
         end subroutine MatVec
      end interface
      !Arguments
      real(8)                   :: eval(:)![Neigen]
      real(8)                   :: evec(:,:)![Ns,Neigen]
      integer,optional          :: Nblock
      integer,optional          :: Nitermax
      ! character(len=2),optional :: which
      real(8),optional          :: v0(size(evec,1))
      real(8),optional          :: tol
      logical,optional          :: iverbose
      logical,optional          :: vrandom
      !Dimensions:
      integer                   :: Ns
      integer                   :: Neigen
      integer                   :: maxn,maxnev,maxncv,ldv
      integer                   :: n,nconv,ncv,nev
      !Arrays:
      real(8),allocatable       :: evec_tmp(:) ![Nloc] see above
      real(8),allocatable       :: ax(:),d(:,:)
      real(8),allocatable       :: resid(:),vec(:)
      real(8),allocatable       :: workl(:),workd(:)
      real(8),allocatable       :: v(:,:)
      logical,allocatable       :: select(:)
      integer                   :: iparam(11)
      integer                   :: ipntr(11)
      !Control Vars:
      integer                   :: ido,ierr,info,ishfts,j,lworkl,maxitr,mode1
      logical                   :: verb,vran
      integer                   :: i
      real(8)                   :: sigma,norm,norm_tmp
      real(8)                   :: tol_
      character                 :: bmat
      character(len=2)          :: which_
      real(8),external          :: dnrm2
      !MPI
      logical                   :: mpi_master
      !
      integer ::  logfil, ndigit, mgetv0,&
         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,&
         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,&
         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
      common /debug/logfil, ndigit, mgetv0,&
         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,&
         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,&
         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
      !
      if(MpiComm == MPI_COMM_NULL)return
      !
      !MPI setup:
      mpi_master= Get_master_MPI(MpiComm)
      !
      Ns     = size(evec,1)
      Neigen = size(eval)
      call assert_shape(Evec,[Ns,Neigen],"P_arpack_d","Evec")
      !
      maxn   = Ns
      maxnev = Neigen
      !
      maxncv = 10*Neigen ; if(present(Nblock))maxncv = Nblock
      maxitr = 512       ; if(present(Nitermax))maxitr = Nitermax
      which_ = 'SA'     ! ; if(present(which))which_=which
      tol_   = 0d0       ; if(present(tol))tol_=tol
      verb   = .false.   ; if(present(iverbose))verb=iverbose
      vran   = .true.    ; if(present(vrandom))vran=vrandom
      if(verb)then
         ndigit=-4
         logfil = 6
         mcaupd=1;mnaupd=1
         mcaup2=1;mnaup2=1
         mceupd=4;mneupd=4
      endif
      if(maxncv>Ns)then
         maxncv=Ns
         print*,"PARPACK WARNING Ncv > Ns: reset block size to ",Ns
      endif
      !BUG FIX FOR THE BLOCK RESIZE STUCK BEHAVIOR, from Xuanyu Long
      call MPI_ALLREDUCE(MPI_IN_PLACE,maxncv,1,MPI_INTEGER,MPI_MIN,MpiComm,ierr)
      !
      ldv    = maxn
      n      = maxn
      nev    = maxnev
      ncv    = maxncv
      bmat   = 'I'
      !
      !=========================================================================
      ! Setup distribution of data to nodes:
      if ( ldv > maxn ) then
         stop ' ERROR with _SDRV1: NLOC is greater than MAXNLOC '
      else if ( nev > maxnev ) then
         stop ' ERROR with _SDRV1: NEV is greater than MAXNEV '
      else if ( ncv > maxncv ) then
         stop ' ERROR with _SDRV1: NCV is greater than MAXNCV '
      end if
      !
      allocate(ax(ldv))
      allocate(resid(ldv))
      allocate(workd(3*ldv))
      allocate(v(ldv,ncv))
      allocate(d(ncv,2))
      allocate(workl(ncv*(ncv+8)))
      allocate(select(ncv))
      !
      ax     =0.d0
      d      =0.d0
      resid  =0.d0
      workl  =0.d0
      workd  =0.d0
      v      =0.d0
      lworkl = ncv*(ncv+8)+10
      info   = 1
      ido    = 0
      select = .true.
      !
      if(present(v0))then
         resid = v0
      else
         allocate(vec(ldv))
         if(vran)then
            ! call random_seed(size=nrandom)
            ! if(allocated(seed_random))deallocate(seed_random)
            ! allocate(seed_random(nrandom))
            ! seed_random=1234567
            ! call random_seed(put=seed_random)
            ! deallocate(seed_random)
            ! call random_number(vec)
            call mt_random(vec)
         else
            vec = 1d0              !start with unitary vector 1/sqrt(Ndim)
         endif
         norm_tmp = dot_product(vec,vec)
         norm = 0d0
         call AllReduce_MPI(MpiComm,norm_tmp,norm)
         resid=vec/sqrt(norm)
         deallocate(vec)
      endif
      !
      ishfts    = 1
      mode1     = 1
      iparam(1) = ishfts
      iparam(3) = maxitr
      iparam(7) = mode1
      !
      !  MAIN LOOP (Reverse communication loop)
      do
         call pdsaupd(MpiComm,ido,bmat,ldv,which_,nev,tol_,resid,&
            ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,info )
         !
         if(ido/=-1.AND.ido/=1)exit
         !  Perform matrix vector multiplication: y <--- OP*x ; workd(ipntr(1))=input, workd(ipntr(2))=output
         call MatVec(ldv,workd(ipntr(1)),workd(ipntr(2)))
      end do
      !
      !POST PROCESSING:
      if(info/=0)then
         if(mpi_master)then
            write(*,'(a,i6)')'Warning/Error in PDSAUPD, info = ', info
            include "error_msg_arpack.h90"
         endif
      else
         !
         call pdseupd (MpiComm,.true.,'All',select,d,v,ldv,sigma,bmat,&
            ldv,which_,nev,tol_,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,ierr)
         do j=1,neigen
            eval(j)=d(j,1)
         enddo
         !
         evec=0d0
         do j=1,neigen
            do i=1,ldv
               evec(i,j)=v(i,j)
            enddo
         enddo
         !
         !  Compute the residual norm ||  A*x - lambda*x ||
         !  for the NCONV accurately computed eigenvalues and eigenvectors.
         if(ierr/=0)then
            write(*,'(a,i6)')'Error with PDSEUPD, IERR = ',ierr
            write(*,'(a)')'Check the documentation of PDSEUPD.'
            stop
         else
            nconv =  iparam(5)
         end if
         !
         if(mpi_master)then
            include "info_msg_arpack.h90"
         endif
         !
         ! if(mpi_master.and.nconv==0.and.verb)stop "None of the required values was found."
      endif
      deallocate(ax,resid,workd,v,d,workl,select)
   end subroutine lanczos_parpack_d
#endif

#ifdef _MPI

   subroutine lanczos_parpack_c(MpiComm,MatVec,eval,evec,Nblock,Nitermax,v0,tol,iverbose,vrandom)
      !Arguments
      integer                    :: MpiComm
      !Interface to Matrix-Vector routine:
      interface
         subroutine MatVec(nchunk,vin,vout)
            integer                      :: nchunk
            complex(8),dimension(nchunk) :: vin,vout
         end subroutine MatVec
      end interface
      !Arguments
      real(8)                   :: eval(:)   ![Neigen]
      complex(8)                :: evec(:,:) ![Nloc,Neigen]
      integer,optional          :: Nblock
      integer,optional          :: Nitermax
      ! character(len=2),optional :: which
      complex(8),optional       :: v0(size(evec,1))
      real(8),optional          :: tol
      logical,optional          :: iverbose
      logical,optional          :: vrandom
      !Dimensions:
      integer                   :: Ns
      integer                   :: Neigen
      integer                   :: maxn,maxnev,maxncv,ldv,nloc
      integer                   :: n,nconv,ncv,nev
      !Arrays:
      complex(8),allocatable    :: evec_tmp(:) ![Nloc] see above
      complex(8),allocatable    :: ax(:),d(:)
      complex(8),allocatable    :: v(:,:)
      complex(8),allocatable    :: workl(:),workd(:),workev(:)
      complex(8),allocatable    :: resid(:),vec(:)
      real(8),allocatable       :: rwork(:),rd(:,:)
      logical,allocatable       :: select(:)
      integer                   :: iparam(11)
      integer                   :: ipntr(14)
      !Control Vars:
      integer                   :: ido,ierr,info,ishfts,j,lworkl,maxitr,mode1
      logical                   :: verb,vran
      integer                   :: i
      real(8)                   :: sigma,norm,norm_tmp
      real(8)                   :: tol_
      character                 :: bmat
      character(len=2)          :: which_
      real(8),external          :: dznrm2,dlapy2
      real(8),allocatable       :: reV(:),imV(:)
      integer,allocatable       :: Eorder(:)
      !MPI
      logical                   :: mpi_master

      integer ::  logfil, ndigit, mgetv0,&
         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,&
         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,&
         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
      common /debug/logfil, ndigit, mgetv0,&
         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,&
         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,&
         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
      !
      if(MpiComm == MPI_COMM_NULL)return
      !
      !MPI setup:
      mpi_master= Get_master_MPI(MpiComm)
      !
      Ns     = size(evec,1)
      Neigen = size(eval)
      call assert_shape(Evec,[Ns,Neigen],"P_arpack_d","Evec")
      !
      maxn   = Ns
      maxnev = Neigen
      !
      maxncv = 10*Neigen ; if(present(Nblock))maxncv = Nblock
      maxitr = 512       ; if(present(Nitermax))maxitr = Nitermax
      which_ = 'SR'      !; if(present(which))which_=which
      tol_   = 0d0       ; if(present(tol))tol_=tol
      verb   = .false.   ; if(present(iverbose))verb=iverbose
      vran   = .true.    ; if(present(vrandom))vran=vrandom
      if(verb)then
         ndigit=-4
         logfil = 6
         mcaupd=1;mnaupd=1
         mcaup2=1;mnaup2=1
         mceupd=4;mneupd=4
      endif
      !
      !
      if(maxncv>Ns)then
         maxncv=Ns
         print*,"PARPACK WARNING Ncv > Ns: reset block size to ",Ns
      endif
      !BUG FIX FOR THE BLOCK RESIZE STUCK BEHAVIOR, from Xuanyu Long
      call MPI_ALLREDUCE(MPI_IN_PLACE,maxncv,1,MPI_INTEGER,MPI_MIN,MpiComm,ierr)
      !

      !
      n      = maxn
      nev    = maxnev
      ncv    = maxncv
      bmat   = 'I'
      !
      !=========================================================================
      ! Setup distribution of data to nodes:
      ldv = size(evec,1)             !ldv is the SMALL dimension
      ! if(ldv>0.AND.ldv < ncv)stop "LANCZOS_PARPACK_C error: ldv < maxNblock. Unstable! Find a way to increase ldv... (less cpu?)"
      if ( ldv > maxn ) then
         stop ' ERROR with _SDRV1: NLOC is greater than MAXNLOC '
      else if ( nev > maxnev ) then
         stop ' ERROR with _SDRV1: NEV is greater than MAXNEV '
      else if ( ncv > maxncv ) then
         stop ' ERROR with _SDRV1: NCV is greater than MAXNCV '
      end if
      !
      allocate(ax(ldv))
      allocate(d(ncv))
      allocate(resid(ldv))
      allocate(v(ldv,ncv))
      allocate(workd(3*ldv))
      allocate(workev(3*ncv))
      allocate(workl(ncv*(3*ncv+5) + 10))
      allocate(rwork(ncv))
      allocate(rd(ncv,3))
      allocate(select(ncv))
      !
      ax     = zero
      d      = zero
      v      = zero
      workl  = zero
      workd  = zero
      resid  = zero
      workev = zero
      rwork  = zero
      rd     = zero
      lworkl = ncv*(3*ncv+5) + 10
      info   = 1
      ido    = 0
      select = .true.
      !
      !
      if(present(v0))then
         resid = v0
      else
         ! allocate(Vec(ldv))
         if(vran)then
            ! call random_seed(size=nrandom)
            ! if(allocated(seed_random))deallocate(seed_random)
            ! allocate(seed_random(nrandom))
            ! seed_random=1234567
            ! call random_seed(put=seed_random)
            ! allocate(reV(ldv),imV(ldv))
            ! call random_number(reV)
            ! call random_number(imV)
            ! call mt_random(reV)
            ! call mt_random(imV)
            ! vec=dcmplx(reV,imV)
            ! deallocate(reV,imV)
            call mt_random(resid)
         else
            ! vec = one               !start with unitary vector 1/sqrt(Ndim)
            resid = one               !start with unitary vector 1/sqrt(Ndim)
         endif
         norm_tmp = dot_product(resid,resid)
         norm = 0d0
         call AllReduce_MPI(MpiComm,norm_tmp,norm)
         resid=resid/sqrt(norm)
         ! deallocate(Vec)
      endif
      !
      ishfts    = 1
      mode1     = 1
      iparam(1) = ishfts
      iparam(3) = maxitr
      iparam(7) = mode1
      !
      !
      !  MAIN LOOP (Reverse communication loop)
      do
         call pznaupd(MpiComm,ido,bmat,ldv,which_,nev,tol_,resid,&
            ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,rwork,info)
         !
         if(ido/=-1.AND.ido/=1)exit
         !  Perform matrix vector multiplication:
         !    y <--- OP*x ; workd(ipntr(1))=input, workd(ipntr(2))=output
         call MatVec(ldv,workd(ipntr(1)),workd(ipntr(2)))
      end do
      !
      !POST PROCESSING:
      if(info/=0)then
         if(mpi_master)then
            write(*,'(a,i6)')'Warning/Error in PZNAUPD, info = ', info
            include "error_msg_arpack.h90"
         endif
      else
         !
         call pzneupd (MpiComm,.true.,'All',select,d,v,ldv,sigma,workev,bmat,&
            ldv,which_,nev,tol_,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,rwork,ierr)
         !
         do j=1,neigen
            eval(j)=dreal(d(j))
         enddo
         allocate(Eorder(Neigen))
         call sort_array(Eval,Eorder)
         !
         evec=zero
         do j=1,neigen
            do i=1,ldv
               evec(i,j)=v(i,Eorder(j))
            enddo
         enddo
         !
         !  Compute the residual norm ||  A*x - lambda*x ||
         !  for the NCONV accurately computed eigenvalues and eigenvectors.
         if(ierr/=0)then
            write(*,'(a,i6)')'Error with PZNEUPD, IERR = ',ierr
            write(*,'(a)')'Check the documentation of PZNEUPD.'
         else
            nconv =  iparam(5)
         end if
         !
         if(mpi_master)then
            include "info_msg_arpack.h90"
         end if
         !
         !
         ! if(mpi_master.and.nconv==0.and.verb)stop "None of the required values was found."
      endif
      deallocate(ax,d,resid,v,workd,workev,workl,rwork,rd,select)
      !
      !
   contains
      !
      !+------------------------------------------------------------------+
      !PURPOSE  : Sort an array, gives the new ordering of the label.
      !+------------------------------------------------------------------+
      subroutine sort_array(array,order2,no_touch_array)
         implicit none
         real(8),dimension(:)                    :: array
         real(8),dimension(size(array))          :: backup
         integer,dimension(size(array))          :: order
         integer,dimension(size(array)),optional :: order2
         integer                                 :: i
         logical,optional                        :: no_touch_array
         do i=1,size(order)
            order(i)=i
         enddo
         call qsort_sort( array, order, 1, size(array) )
         if(.not.present(no_touch_array))then
            do i=1,size(order)
               backup(i)=array(order(i))
            enddo
            array=backup
         endif
         if(present(order2)) order2=order
      end subroutine sort_array
      !---------------------------------------------!
      recursive subroutine qsort_sort( array, order, left, right )
         implicit none
         real(8), dimension(:)                 :: array
         integer, dimension(:)                 :: order
         integer                               :: left
         integer                               :: right
         integer                               :: i
         integer                               :: last
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
         implicit none
         integer, dimension(:)                 :: order
         integer                               :: first, second
         integer                               :: tmp
         tmp           = order(first)
         order(first)  = order(second)
         order(second) = tmp
      end subroutine qsort_swap
      !---------------------------------------------!
      integer function qsort_rand( lower, upper )
         implicit none
         integer                               :: lower, upper
         real(4)                               :: r
         r=drand()
         qsort_rand =  lower + nint(r * (upper-lower))
      end function qsort_rand
      !---------------------------------------------!
      function drand()
         implicit none
         real(8)                               :: drand
         real(4)                               :: r
         call random_number(r)
         drand=dble(r)
      end function drand
      !---------------------------------------------!
      function compare(f,g)
         implicit none
         real(8)                               :: f,g
         integer                               :: compare
         if(f<g) then
            compare=-1
         else
            compare=1
         endif
      end function compare
      !
      !
   end subroutine lanczos_parpack_c

#endif


   !##################################################################
   ! LANCZOS METHOD for LOWEST EigenSolution OR tri-diagonalization
   !    of a SPARSE MATRIX (defined via H*v)
   ! - DBLE and CMPLX versions included
   ! - SERIAL and PARALLEL-MPI versions included
   include "lanczos_d.f90"
   include "lanczos_c.f90"
#ifdef _MPI
   include "mpi_lanczos_d.f90"
   include "mpi_lanczos_c.f90"
#endif


   !##################################################################
   ! DAVIDSON METHOD for LOWEST EigenSolution of a SPARSE MATRIX (defined via H*v)
   include "dvdson_serial.f90"



end module SF_SP_LINALG
























! !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! !++++++++++++++++++COMPUTATIONAL ROUTINE: TQL2++++++++++++++++++++++++
! !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! !---------------------------------------------------------------------
! ! PURPOSE computes all eigenvalues/vectors, real symmetric tridiagonal matrix.
! !    This subroutine finds the eigenvalues and eigenvectors of a symmetric
! !    tridiagonal matrix by the QL method.  The eigenvectors of a full
! !    symmetric matrix can also be found if TRED2 has been used to reduce this
! !    full matrix to tridiagonal form.
! !  Parameters:
! !    Input, integer ( kind = 4 ) N, the order of the matrix.
! !
! !    Input/output, real ( kind = 8 ) D(N).  On input, the diagonal elements of
! !    the matrix.  On output, the eigenvalues in ascending order.  If an error
! !    exit is made, the eigenvalues are correct but unordered for indices
! !    1,2,...,IERR-1.
! !
! !    Input/output, real ( kind = 8 ) E(N).  On input, E(2:N) contains the
! !    subdiagonal elements of the input matrix, and E(1) is arbitrary.
! !    On output, E has been destroyed.
! !
! !    Input, real ( kind = 8 ) Z(N,N).  On input, the transformation matrix
! !    produced in the reduction by TRED2, if performed.  If the eigenvectors of
! !    the tridiagonal matrix are desired, Z must contain the identity matrix.
! !    On output, Z contains the orthonormal eigenvectors of the symmetric
! !    tridiagonal (or full) matrix.  If an error exit is made, Z contains
! !    the eigenvectors associated with the stored eigenvalues.
! !
! !    Output, integer ( kind = 4 ) IERR, error flag.
! !    0, normal return,
! !    J, if the J-th eigenvalue has not been determined after
! !    30 iterations.
! !
! !---------------------------------------------------------------------
! subroutine tql2 ( n, d, e, z, ierr )
!   integer :: n
!   real(8) :: c
!   real(8) :: c2
!   real(8) :: c3
!   real(8) :: d(n)
!   real(8) :: dl1
!   real(8) :: e(n)
!   real(8) :: el1
!   real(8) :: f
!   real(8) :: g
!   real(8) :: h
!   integer ( kind = 4 ) i
!   integer ( kind = 4 ) ierr
!   integer ( kind = 4 ) ii
!   integer ( kind = 4 ) j
!   integer ( kind = 4 ) k
!   integer ( kind = 4 ) l
!   integer ( kind = 4 ) l1
!   integer ( kind = 4 ) l2
!   integer ( kind = 4 ) m
!   integer ( kind = 4 ) mml
!   real(8) :: p
!   real(8) :: r
!   real(8) :: s
!   real(8) :: s2
!   real(8) :: tst1
!   real(8) :: tst2
!   real(8) :: z(n,n)
!   ierr = 0
!   if ( n == 1 ) then
!      return
!   end if
!   do i = 2, n
!      e(i-1) = e(i)
!   end do
!   f = 0.0D+00
!   tst1 = 0.0D+00
!   e(n) = 0.0D+00
!   do l = 1, n
!      j = 0
!      h = abs ( d(l) ) + abs ( e(l) )
!      tst1 = max ( tst1, h )
!      !
!      !  Look for a small sub-diagonal element.
!      !
!      do m = l, n
!         tst2 = tst1 + abs ( e(m) )
!         if ( tst2 == tst1 ) then
!            exit
!         end if
!      end do
!      if ( m == l ) then
!         go to 220
!      end if
! 130  continue
!      if ( 30 <= j ) then
!         ierr = l
!         return
!      end if
!      j = j + 1
!      !
!      !  Form shift.
!      !
!      l1 = l + 1
!      l2 = l1 + 1
!      g = d(l)
!      p = ( d(l1) - g ) / ( 2.0D+00 * e(l) )
!      r = pythag ( p, 1.0D+00 )
!      d(l) = e(l) / ( p + sign ( r, p ) )
!      d(l1) = e(l) * ( p + sign ( r, p ) )
!      dl1 = d(l1)
!      h = g - d(l)
!      d(l2:n) = d(l2:n) - h
!      f = f + h
!      !
!      !  QL transformation.
!      !
!      p = d(m)
!      c = 1.0D+00
!      c2 = c
!      el1 = e(l1)
!      s = 0.0D+00
!      mml = m - l
!      do ii = 1, mml
!         c3 = c2
!         c2 = c
!         s2 = s
!         i = m - ii
!         g = c * e(i)
!         h = c * p
!         r = pythag ( p, e(i) )
!         e(i+1) = s * r
!         s = e(i) / r
!         c = p / r
!         p = c * d(i) - s * g
!         d(i+1) = h + s * ( c * g + s * d(i) )
!         !
!         !  Form vector.
!         !
!         do k = 1, n
!            h = z(k,i+1)
!            z(k,i+1) = s * z(k,i) + c * h
!            z(k,i) = c * z(k,i) - s * h
!         end do
!      end do
!      p = - s * s2 * c3 * el1 * e(l) / dl1
!      e(l) = s * p
!      d(l) = c * p
!      tst2 = tst1 + abs ( e(l) )
!      if ( tst2 > tst1 ) then
!         go to 130
!      end if
! 220  continue
!      d(l) = d(l) + f
!   end do
!   !
!   !  Order eigenvalues and eigenvectors.
!   !
!   do ii = 2, n
!      i = ii - 1
!      k = i
!      p = d(i)
!      do j = ii, n
!         if ( d(j) < p ) then
!            k = j
!            p = d(j)
!         end if
!      end do
!      if ( k /= i ) then
!         d(k) = d(i)
!         d(i) = p
!         do j = 1, n
!            call r8_swap ( z(j,i), z(j,k) )
!         end do
!      end if
!   end do
!   return
! end subroutine tql2


! !---------------------------------------------------------------------
! ! PURPOSE: computes SQRT ( A * A + B * B ) carefully.
! !    The formula
! !    PYTHAG = sqrt ( A * A + B * B )
! !    is reasonably accurate, but can fail if, for example, A**2 is larger
! !    than the machine overflow.  The formula can lose most of its accuracy
! !    if the sum of the squares is very large or very small.
! !  Parameters:
! !    Input, real(8) :: A, B, the two legs of a right triangle.
! !    Output, real(8) :: PYTHAG, the length of the hypotenuse.
! !---------------------------------------------------------------------
! function pythag ( a, b )
!   implicit none
!   real(8) :: a
!   real(8) :: b
!   real(8) :: p
!   real(8) :: pythag
!   real(8) :: r
!   real(8) :: s
!   real(8) :: t
!   real(8) :: u
!   p = max ( abs ( a ), abs ( b ) )
!   if ( p /= 0.0D+00 ) then
!      r = ( min ( abs ( a ), abs ( b ) ) / p )**2
!      do
!         t = 4.0D+00 + r
!         if ( t == 4.0D+00 ) then
!            exit
!         end if
!         s = r / t
!         u = 1.0D+00 + 2.0D+00 * s
!         p = u * p
!         r = ( s / u )**2 * r
!      end do
!   end if
!   pythag = p
!   return
! end function pythag

! !---------------------------------------------------------------------
! ! PURPOSE: swaps two R8's.
! !  Parameters:
! !    Input/output, real(8) :: X, Y.  On output, the values of X and
! !    Y have been interchanged.
! !---------------------------------------------------------------------
! subroutine r8_swap ( x, y )
!   real(8) :: x
!   real(8) :: y
!   real(8) :: z
!   z = x
!   x = y
!   y = z
!   return
! end subroutine r8_swap



