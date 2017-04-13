module SF_SP_LINALG
#ifdef _MPI
  USE MPI
#endif
  implicit none
  private 


  abstract interface 
     subroutine lanc_htimesv_d(nloc,vin,vout)
       integer :: nloc
       real(8) :: vin(nloc)
       real(8) :: vout(nloc)
     end subroutine lanc_htimesv_d
     !
     subroutine lanc_htimesv_c(nloc,vin,vout)
       integer    :: nloc
       complex(8) :: vin(nloc)
       complex(8) :: vout(nloc)
     end subroutine lanc_htimesv_c
  end interface


  interface sp_lanc_eigh
     module procedure :: lanczos_plain_d
     module procedure :: lanczos_plain_c
#ifdef _MPI
     module procedure :: p_lanczos_plain_d
     module procedure :: p_lanczos_plain_c
#endif
  end interface sp_lanc_eigh


  interface sp_lanc_tridiag
     module procedure :: lanczos_plain_tridiag_d
     module procedure :: lanczos_plain_tridiag_c
  end interface sp_lanc_tridiag


  interface lanczos_iteration
     module procedure :: lanczos_plain_iteration_d
     module procedure :: lanczos_plain_iteration_c
  end interface lanczos_iteration


  interface sp_eigh
     module procedure :: lanczos_arpack_d
     module procedure :: lanczos_arpack_c
#ifdef _MPI
     module procedure :: lanczos_parpack_d
     module procedure :: lanczos_parpack_c
#endif
  end interface sp_eigh


  complex(8),parameter              :: zero=(0d0,0d0),one=(1d0,0d0),xi=(0d0,1d0)
  integer,allocatable               :: seed_random(:)
  integer                           :: nrandom

  procedure(lanc_htimesv_d),pointer :: dp_hprod
  procedure(lanc_htimesv_c),pointer :: cp_hprod
  logical                           :: verb=.false.
  real(8)                           :: threshold_=1.d-13
  integer                           :: ncheck_=10

  integer                           :: lanc_mpi_comm
  integer                           :: lanc_mpi_rank=0
  integer                           :: lanc_mpi_size=1
  logical                           :: lanc_mpi_master=.true.





  !****************************************************************************************
  !                                      PUBLIC 
  !****************************************************************************************
  public :: sp_eigh
  !
  public :: sp_lanc_eigh
  public :: sp_lanc_tridiag
  !****************************************************************************************




contains





  subroutine lanczos_arpack_d(MatVec,Ns,Neigen,Nblock,Nitermax,eval,evec,which,v0,tol,iverbose)
    !Interface to Matrix-Vector routine:
    interface
       subroutine MatVec(Nloc,vin,vout)
         integer                 :: Nloc
         real(8),dimension(Nloc) :: vin
         real(8),dimension(Nloc) :: vout
       end subroutine MatVec
    end interface
    !Arguments
    integer                   :: Ns
    integer                   :: Neigen
    integer                   :: Nblock
    integer                   :: Nitermax
    real(8)                   :: eval(Neigen)
    real(8)                   :: evec(Ns,Neigen)
    character(len=2),optional :: which
    real(8),optional          :: v0(ns)
    real(8),optional          :: tol
    logical,optional          :: iverbose
    !Dimensions:
    integer                   :: maxn,maxnev,maxncv,ldv
    integer                   :: n,nconv,ncv,nev
    !Arrays:
    real(8),allocatable       :: ax(:),d(:,:)
    real(8),allocatable       :: resid(:)
    real(8),allocatable       :: workl(:),workd(:)
    real(8),allocatable       :: v(:,:)
    logical,allocatable       :: select(:)
    integer                   :: iparam(11)
    integer                   :: ipntr(11)
    !Control Vars:
    integer                   :: ido,ierr,info,ishfts,j,lworkl,maxitr,mode1
    logical                   :: rvec,verb
    integer                   :: i
    real(8)                   :: sigma
    real(8)                   :: tol_
    character                 :: bmat  
    character(len=2)          :: which_
    real(8),external          :: dnrm2
    !
    which_='SA';if(present(which))which_=which
    tol_=1d-12;if(present(tol))tol_=tol
    verb=.false.;if(present(iverbose))verb=iverbose
    !
    maxn   = Ns
    maxnev = Neigen
    maxncv = Nblock
    ldv    = Ns
    if(maxncv>Ns)maxncv=Ns
    n      = maxn
    nev    = maxnev
    ncv    = maxncv
    bmat   = 'I'
    maxitr = Nitermax
    ! 
    allocate(ax(n))
    allocate(d(ncv,2))
    allocate(resid(n))
    allocate(workl(ncv*(ncv+8)))
    allocate(workd(3*n))
    allocate(v(ldv,ncv))
    allocate(select(ncv))
    ax     =0d0
    d      =0d0
    resid  =0d0
    workl  =0d0
    workd  =0d0
    v      =0d0
    lworkl = ncv*(ncv+8)
    info   = 1
    ido    = 0
    ishfts    = 1
    mode1     = 1
    iparam(1) = ishfts
    iparam(3) = maxitr
    iparam(7) = mode1
    if(present(v0))then
       resid=v0
    else
       call random_seed(size=nrandom)
       if(allocated(seed_random))deallocate(seed_random)
       allocate(seed_random(nrandom))
       seed_random=1234567
       call random_seed(put=seed_random)
       call random_number(resid)
    endif
    resid=resid/sqrt(dot_product(resid,resid))
    !
    !MAIN LOOP, REVERSE COMMUNICATION
    do
       call dsaupd(ido,bmat,n,which_,nev,tol_,resid,ncv,v,ldv,&
            iparam,ipntr,workd,workl,lworkl,info)
       if(ido/=-1.AND.ido/=1)then
          exit
       end if
       !  Perform matrix vector multiplication
       !    y <--- OP*x ; workd(ipntr(1))=input, workd(ipntr(2))=output
       call MatVec(N,workd(ipntr(1)),workd(ipntr(2)) )
    end do
    !
    !POST PROCESSING:
    if(info<0)then
       write(*,'(a,i6)')'Error in DSAUPD, info = ', info
       stop
    else
       rvec = .true.
       call dseupd(rvec,'All',select,d,v,ldv,sigma,bmat,n,which_,&
            nev,tol_,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,ierr)
       do j=1,neigen
          eval(j)=d(j,1)
          ! if(present(evec))then
          do i=1,ns
             evec(i,j)=v(i,j)
          enddo
          ! endif
       enddo
       !
       !  Compute the residual norm
       !    ||  A*x - lambda*x ||
       !  for the NCONV accurately computed eigenvalues and 
       !  eigenvectors.  (iparam(5) indicates how many are 
       !  accurate to the requested tolerance)
       if(ierr/=0)then
          write(*,'(a,i6)')'Error with DSEUPD (get Evec), ierr = ',ierr
          ! else
          !    nconv =  iparam(5)
          !    do j = 1, nconv
          !       call MatVec(n, v(1,j), ax )
          !       call daxpy( n, -d(j,1), v(1,j), 1, ax, 1 )
          !       d(j,2) = dnrm2(n,ax,1)
          !       d(j,2) = d(j,2) / abs ( d(j,1) )
          !    end do
          !    if(verb)call dmout(6,nconv,2,d,maxncv,-6,'Ritz values and relative residuals')
       end if
       !
       if(info==1) then
          write(*,'(a)' ) ' '
          write(*,'(a)' ) '  Maximum number of iterations reached.'
       elseif(info==3) then
          write(*,'(a)' ) ' '
          write(*,'(a)' ) '  No shifts could be applied during implicit '&
               //'Arnoldi update, try increasing NCV.'
       end if
    endif
  end subroutine lanczos_arpack_d

  subroutine lanczos_arpack_c(MatVec,Ns,Neigen,Nblock,Nitermax,eval,evec,which,v0,tol,iverbose)
    !Interface to Matrix-Vector routine:
    interface
       subroutine MatVec(Nloc,vin,vout)
         integer                    :: Nloc
         complex(8),dimension(Nloc) :: vin
         complex(8),dimension(Nloc) :: vout
       end subroutine MatVec
    end interface
    !Arguments
    integer                      :: Ns
    integer                      :: Neigen
    integer                      :: Nblock
    integer                      :: Nitermax
    real(8)                      :: eval(Neigen)
    complex(8)                   :: evec(Ns,Neigen)
    character(len=2),optional    :: which
    complex(8),optional          :: v0(ns)
    real(8),optional             :: tol
    logical,optional             :: iverbose
    !Dimensions:
    integer                      :: maxn,maxnev,maxncv,ldv
    integer                      :: n,nconv,ncv,nev
    !Arrays:
    complex(8),allocatable       :: ax(:),d(:)
    complex(8),allocatable       :: v(:,:)
    complex(8),allocatable       :: workl(:),workd(:),workev(:)
    complex(8),allocatable       :: resid(:)
    real(8),allocatable          :: rwork(:),rd(:,:)
    logical,allocatable          :: select(:)
    integer                      :: iparam(11)
    integer                      :: ipntr(14)
    !Control Vars:
    integer                      :: ido,ierr,info,ishfts,j,lworkl,maxitr,mode1
    logical                      :: rvec,verb
    integer                      :: i
    real(8)                      :: sigma
    real(8)                      :: tol_
    character                    :: bmat  
    character(len=2)             :: which_
    real(8),external             :: dznrm2,dlapy2
    real(8),allocatable          :: reV(:),imV(:)
    !
    which_='SR';if(present(which))which_=which
    tol_=1d-12;if(present(tol))tol_=tol
    verb=.false.;if(present(iverbose))verb=iverbose
    !
    maxn   = Ns
    maxnev = Neigen
    maxncv = Nblock
    ldv    = Ns
    if(maxncv>Ns)maxncv=Ns
    n      = maxn
    nev    = maxnev
    ncv    = maxncv
    bmat   = 'I'
    maxitr = Nitermax
    ! 
    allocate(ax(n))
    allocate(d(ncv))
    allocate(resid(n))
    allocate(v(ldv,ncv))
    allocate(workd(3*n))
    allocate(workev(3*ncv))
    allocate(workl(ncv*(3*ncv+5) + 10))
    allocate(rwork(ncv))
    allocate(rd(ncv,3))
    allocate(select(ncv))
    ax     = zero
    d      = zero
    v      = zero
    workl  = zero
    workd  = zero
    resid  = zero
    workev = zero
    rwork  = 0.d0
    rd     = 0.d0
    lworkl = ncv*(3*ncv+5) + 10
    info   = 1
    ido    = 0
    ishfts    = 1
    mode1     = 1
    iparam(1) = ishfts
    iparam(3) = maxitr
    iparam(7) = mode1
    if(present(v0))then
       resid=v0
    else
       call random_seed(size=nrandom)
       if(allocated(seed_random))deallocate(seed_random)
       allocate(seed_random(nrandom))
       seed_random=1234567
       call random_seed(put=seed_random)
       allocate(reV(size(resid)),imV(size(resid)))
       call random_number(reV)
       call random_number(imV)
       resid=dcmplx(reV,imV)
       deallocate(reV,imV)
       resid=resid/sqrt(dot_product(resid,resid))
    endif
    resid=resid/sqrt(dot_product(resid,resid))
    !
    !MAIN LOOP, REVERSE COMMUNICATION
    do
       call znaupd(ido,bmat,n,which_,nev,tol_,resid,ncv,v,ldv,&
            iparam,ipntr,workd,workl,lworkl,rwork,info)
       if(ido/=-1.AND.ido/=1)exit
       !  Perform matrix vector multiplication
       !    y <--- OP*x ; workd(ipntr(1))=input, workd(ipntr(2))=output
       call MatVec(N,workd(ipntr(1)),workd(ipntr(2)) )
    end do
    !
    !POST PROCESSING:
    if(info<0)then
       write(*,'(a,i6)')'Error in ZSAUPD, info = ', info
       stop
    else
       rvec = .true.
       call zneupd  (rvec,'A',select,d,v,ldv,sigma,workev,bmat,n,which_,&
            nev,tol_,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,rwork,ierr)
       !  Eigenvalues are returned in the first column of the two dimensional 
       !  array D and the corresponding eigenvectors are returned in the first 
       !  NCONV (=IPARAM(5)) columns of the two dimensional array V if requested.
       !  Otherwise, an orthogonal basis for the invariant subspace corresponding 
       !  to the eigenvalues in D is returned in V.
       do j=1,neigen
          eval(j)=dreal(d(j))
          ! if(present(evec))then
          do i=1,ns
             evec(i,j)=v(i,j)
          enddo
          ! endif
       enddo
       !
       !  Compute the residual norm
       !    ||  A*x - lambda*x ||
       !  for the NCONV accurately computed eigenvalues and 
       !  eigenvectors.  (iparam(5) indicates how many are 
       !  accurate to the requested tolerance)
       if(ierr/=0)then
          write(*,'(a,i6)')'Error with DSEUPD (get Evec), ierr = ',ierr
          ! else
          !    nconv =  iparam(5)
          !    do j = 1, nconv
          !       call hprod(1, n, v(1,j), ax )
          !       call zaxpy( n, -d(j), v(1,j), 1, ax, 1 )
          !       rd(j,1) = dble (d(j))
          !       rd(j,2) = dimag (d(j))
          !       rd(j,3) = dznrm2 (n, ax, 1)
          !       rd(j,3) = rd(j,3) / dlapy2 (rd(j,1),rd(j,2))
          !    end do
          !    if(verb)call dmout(6,nconv,3,rd,maxncv,-6,'Ritz values and relative residuals')
       end if
       !
       if(info==1) then
          write(*,'(a)' ) ' '
          write(*,'(a)' ) '  Maximum number of iterations reached.'
       elseif(info==3) then
          write(*,'(a)' ) ' '
          write(*,'(a)' ) '  No shifts could be applied during implicit '&
               //'Arnoldi update, try increasing NCV.'
       end if
    endif
  end subroutine lanczos_arpack_c





#ifdef _MPI
  subroutine lanczos_parpack_d(MpiComm,MatVec,Ns,Neigen,Nblock,Nitermax,eval,evec,which,v0,tol,iverbose)
    !Arguments
    integer                      :: MpiComm
    !Interface to Matrix-Vector routine:
    interface
       subroutine MatVec(nchunk,vin,vout)
         integer                   :: nchunk
         real(8),dimension(nchunk) :: vin,vout
       end subroutine MatVec
    end interface
    !Arguments
    integer                      :: Ns
    integer                      :: Neigen
    integer                      :: Nblock
    integer                      :: Nitermax
    real(8)                      :: eval(Neigen)
    real(8)                      :: evec(Ns,Neigen)
    character(len=2),optional    :: which
    real(8),optional             :: v0(ns)
    real(8),optional             :: tol
    logical,optional             :: iverbose
    !Dimensions:
    integer                      :: maxn,maxnev,maxncv,ldv
    integer                      :: n,nconv,ncv,nev
    !Arrays:
    real(8)                      :: evec_tmp(Ns)
    real(8),allocatable          :: ax(:),d(:,:)
    real(8),allocatable          :: resid(:),vec(:)
    real(8),allocatable          :: workl(:),workd(:)
    real(8),allocatable          :: v(:,:)
    logical,allocatable          :: select(:)
    integer                      :: iparam(11)
    integer                      :: ipntr(11)
    !Control Vars:
    integer                      :: ido,ierr,info,ishfts,j,lworkl,maxitr,mode1
    logical                      :: rvec,verb
    integer                      :: i
    real(8)                      :: sigma
    real(8)                      :: tol_
    character                    :: bmat  
    character(len=2)             :: which_
    real(8),external             :: dnrm2
    !MPI
    integer                      :: mpi_ierr
    integer                      :: mpi_rank
    integer                      :: mpi_size
    logical                      :: mpi_master
    integer                      :: mpiQ
    integer                      :: mpiR
    !
    which_='SA';if(present(which))which_=which
    tol_=1d-12;if(present(tol))tol_=tol
    verb=.false.;if(present(iverbose))verb=iverbose
    !
    !MPI setup:
    mpi_size  = MPI_Get_size(MpiComm)
    mpi_rank  = MPI_Get_rank(MpiComm)
    mpi_master= MPI_Get_master(MpiComm)
    mpiQ = Ns/mpi_size
    mpiR = 0
    if(mpi_rank == mpi_size-1)mpiR=mod(Ns,mpi_size)
    maxn   = Ns
    maxnev = Neigen
    maxncv = Nblock
    ldv    = maxn
    if(maxncv>Ns)maxncv=Ns
    !
    n      = maxn
    nev    = maxnev
    ncv    = maxncv
    bmat   = 'I'
    maxitr = Nitermax
    !=========================================================================
    ! Setup distribution of data to nodes:
    ldv = mpiQ+mpiR             !ldv is the SMALL dimension
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
    lworkl = ncv*(ncv+8)
    info   = 1
    ido    = 0
    allocate(vec(n))
    if(present(v0))then
       vec=v0
    else
       call random_seed(size=nrandom)
       if(allocated(seed_random))deallocate(seed_random)
       allocate(seed_random(nrandom))
       seed_random=1234567
       call random_seed(put=seed_random)
       call random_number(vec)
    endif
    vec=vec/sqrt(dot_product(vec,vec))
    do i=1 + mpi_rank*mpiQ, (mpi_rank+1)*mpiQ + mpiR
       resid(i-mpi_rank*mpiQ)=vec(i)
    enddo
    ishfts    = 1
    mode1     = 1
    iparam(1) = ishfts
    iparam(3) = maxitr
    iparam(7) = mode1

    !  MAIN LOOP (Reverse communication loop)
    do
       call pdsaupd(MpiComm,ido,bmat,ldv,which_,nev,tol_,resid,&
            ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,info )
       if(ido/=-1.AND.ido/=1)exit
       !  Perform matrix vector multiplication
       !    y <--- OP*x ; workd(ipntr(1))=input, workd(ipntr(2))=output
       call MatVec(ldv,workd(ipntr(1)),workd(ipntr(2)))
    end do

    if(info<0)then
       if(mpi_master)write(*,'(a,i6)')'Fatal Error in PDSAUPD, info = ', info
       stop
    else
       rvec = .true.
       call pdseupd (MpiComm,rvec,'All',select,d,v,ldv,sigma,bmat,&
            ldv,which_,nev,tol_,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,ierr)
       do j=1,neigen
          eval(j)=d(j,1)
       enddo
       ! if(present(evec))then
       evec=0d0
       do j=1,neigen
          evec_tmp=0d0
          do i=1 + mpi_rank*mpiQ , (mpi_rank+1)*mpiQ+mpiR
             evec_tmp(i)=v(i-mpi_rank*mpiQ,j)
          enddo
          call MPI_Allreduce(evec_tmp,evec(:,j),Ns,MPI_DOUBLE_PRECISION,MPI_SUM,MpiComm,mpi_ierr)
       enddo
       ! endif
       nconv =  iparam(5)
       !
       if(mpi_master)then
          if(info==1) then
             write(*,'(a)' ) ' '
             write(*,'(a)' ) '  Maximum number of iterations reached.'
          elseif(info==3) then
             write(*,'(a)' ) ' '
             write(*,'(a)' ) '  No shifts could be applied during implicit '&
                  //'Arnoldi update, try increasing NCV.'
          end if
       endif
    endif
    call mpi_barrier(MpiComm,mpi_ierr)
    deallocate(ax,resid,workd,v,d,workl,select)
  end subroutine lanczos_parpack_d


  subroutine lanczos_parpack_c(MpiComm,MatVec,Ns,Neigen,Nblock,Nitermax,eval,evec,which,v0,tol,iverbose)
    !Arguments
    integer                           :: MpiComm
    !Interface to Matrix-Vector routine:
    interface
       subroutine MatVec(nchunk,vin,vout)
         integer                      :: nchunk
         complex(8),dimension(nchunk) :: vin,vout
       end subroutine MatVec
    end interface
    !Arguments
    integer                           :: Ns
    integer                           :: Neigen
    integer                           :: Nblock
    integer                           :: Nitermax
    real(8)                           :: eval(Neigen)
    complex(8)                        :: evec(Ns,Neigen)
    character(len=2),optional         :: which
    complex(8),optional               :: v0(ns)
    real(8),optional                  :: tol
    logical,optional                  :: iverbose
    !Dimensions:
    integer                           :: maxn,maxnev,maxncv,ldv
    integer                           :: n,nconv,ncv,nev
    !Arrays:
    complex(8)                        :: evec_tmp(Ns)
    complex(8),allocatable            :: ax(:),d(:)
    complex(8),allocatable            :: v(:,:)
    complex(8),allocatable            :: workl(:),workd(:),workev(:)
    complex(8),allocatable            :: resid(:),vec(:)
    real(8),allocatable               :: rwork(:),rd(:,:)
    logical,allocatable               :: select(:)
    integer                           :: iparam(11)
    integer                           :: ipntr(14)
    !Control Vars:
    integer                           :: ido,ierr,info,ishfts,j,lworkl,maxitr,mode1
    logical                           :: rvec,verb
    integer                           :: i
    real(8)                           :: sigma
    real(8)                           :: tol_
    character                         :: bmat  
    character(len=2)                  :: which_
    real(8),external                  :: dznrm2,dlapy2
    real(8),allocatable               :: reV(:),imV(:)
    !MPI
    integer                           :: mpi_ierr
    integer                           :: mpi_rank
    integer                           :: mpi_size
    logical                           :: mpi_master
    integer                           :: mpiQ
    integer                           :: mpiR
    !
    which_='SR';if(present(which))which_=which
    tol_=1d-12;if(present(tol))tol_=tol
    verb=.false.;if(present(iverbose))verb=iverbose
    !
    !MPI setup:
    mpi_size  = MPI_Get_size(MpiComm)
    mpi_rank  = MPI_Get_rank(MpiComm)
    mpi_master= MPI_Get_master(MpiComm)
    mpiQ = Ns/mpi_size
    mpiR = 0
    if(mpi_rank == mpi_size-1)mpiR=mod(Ns,mpi_size)
    !
    maxn   = Ns
    maxnev = Neigen
    maxncv = Nblock
    ldv    = maxn
    if(maxncv>Ns)maxncv=Ns
    !
    n      = maxn
    nev    = maxnev
    ncv    = maxncv
    bmat   = 'I'
    maxitr = Nitermax
    !=========================================================================
    ! Setup distribution of data to nodes:
    ldv = mpiQ+mpiR             !ldv is the SMALL dimension
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
    rwork  = 0.d0
    rd     = 0.d0
    lworkl = ncv*(3*ncv+5) + 10
    info   = 1
    ido    = 0
    allocate(vec(n))
    if(present(v0))then
       vec=v0
    else
       call random_seed(size=nrandom)
       if(allocated(seed_random))deallocate(seed_random)
       allocate(seed_random(nrandom))
       seed_random=1234567
       call random_seed(put=seed_random)
       allocate(reV(n),imV(n))
       call random_number(reV)
       call random_number(imV)
       vec=dcmplx(reV,imV)
       deallocate(reV,imV)
    endif
    vec=vec/sqrt(dot_product(vec,vec))
    do i=1 + mpi_rank*mpiQ, (mpi_rank+1)*mpiQ + mpiR
       resid(i-mpi_rank*mpiQ)=vec(i)
    enddo
    ishfts    = 1
    mode1     = 1
    iparam(1) = ishfts
    iparam(3) = maxitr
    iparam(7) = mode1
    !  MAIN LOOP (Reverse communication loop)
    do
       call pznaupd(MpiComm,ido,bmat,ldv,which_,nev,tol_,resid,&
            ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,rwork,info)
       if(ido/=-1.AND.ido/=1)exit
       !  Perform matrix vector multiplication:
       !    y <--- OP*x ; workd(ipntr(1))=input, workd(ipntr(2))=output
       call MatVec(ldv,workd(ipntr(1)),workd(ipntr(2)))
    end do
    if(info<0)then
       if(mpi_master)write(*,'(a,i6)')'Fatal Error in PZNAUPD, info = ', info
       stop
    else
       rvec = .true.
       call pzneupd (MpiComm,rvec,'All',select,d,v,ldv,sigma,workev,bmat,&
            mpiQ+mpiR,which_,nev,tol_,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,rwork,ierr)
       do j=1,neigen
          eval(j)=dreal(d(j))
       enddo
       ! if(present(evec))then
       evec=zero
       do j=1,neigen
          evec_tmp=zero
          do i=1 + mpi_rank*mpiQ , (mpi_rank+1)*mpiQ+mpiR
             evec_tmp(i)=v(i-mpi_rank*mpiQ,j)
          enddo
          call MPI_Allreduce(evec_tmp,evec(:,j),Ns,MPI_DOUBLE_COMPLEX,MPI_SUM,MpiComm,mpi_ierr)
       enddo
       ! endif
       nconv =  iparam(5)
       !
       if(mpi_master)then
          if(info==1) then
             write(*,'(a)' ) ' '
             write(*,'(a)' ) '  Maximum number of iterations reached.'
          elseif(info==3) then
             write(*,'(a)' ) ' '
             write(*,'(a)' ) '  No shifts could be applied during implicit '&
                  //'Arnoldi update, try increasing NCV.'
          end if
       endif
    endif
    call mpi_barrier(MpiComm,mpi_ierr)
    deallocate(ax,resid,workd,v,d,workl,select)
  end subroutine lanczos_parpack_c
#endif



  !##################################################################
  !        LANCZOS METHOD for GROUND STATE of SPARSE MATRIX
  !##################################################################
  !---------------------------------------------------------------------
  !Purpose: use plain lanczos to get the groundstate energy
  !---------------------------------------------------------------------
  subroutine lanczos_plain_d(Hprod,Ns,Nitermax,Egs,Vect,Nlanc,iverbose,threshold,ncheck)
    procedure(lanc_htimesv_d)            :: Hprod
    integer                              :: ns,nitermax
    real(8)                              :: egs
    real(8),dimension(ns)                :: vect
    real(8),dimension(ns)                :: vin,vout
    integer                              :: iter,nlanc
    real(8),dimension(nitermax+1)        :: alanc,blanc
    real(8),dimension(Nitermax,Nitermax) :: Z
    real(8),dimension(Nitermax)          :: diag,subdiag,esave
    real(8)                              :: a_,b_,norm,diff
    integer                              :: i,ierr
    real(8),optional                     :: threshold
    integer,optional                     :: ncheck
    logical,optional                     :: iverbose
    if(present(iverbose))verb=iverbose
    if(present(threshold))threshold_=threshold
    if(present(ncheck))ncheck_=ncheck
    if(associated(dp_hprod))nullify(dp_hprod)
    dp_hprod=>Hprod
    if(.not.associated(dp_hprod))then
       print*,"LANCZOS_PLAIN: dp_hprod is not set. call lanczos_plain_set_htimesv"
       stop
    endif
    norm=dot_product(vect,vect)
    if(norm==0.d0)then
       call random_number(vect)
       vect=vect/sqrt(dot_product(vect,vect))
       if(verb)write(*,*)"LANCZOS_PLAIN: random initial vector generated:"
    endif
    !
    !============= LANCZOS LOOP =====================
    !
    vin = vect
    vout= 0.d0
    alanc=0.d0
    blanc=0.d0
    nlanc=0
    lanc_loop: do iter=1,Nitermax
       if(verb)then
          print*,""
          write(*,*)"Lanczos iteration:",iter
       endif
       call lanczos_plain_iteration_d(iter,vin,vout,a_,b_)
       if(abs(b_)<threshold_)exit lanc_loop
       if(verb)print*,"alanc,blanc=",a_,b_
       nlanc=nlanc+1
       alanc(iter) = a_ ; blanc(iter+1) = b_
       diag = 0.d0 ; subdiag = 0.d0 ; Z = 0.d0
       forall(i=1:Nlanc)Z(i,i)=1.d0
       diag(1:Nlanc)    = alanc(1:Nlanc)
       subdiag(2:Nlanc) = blanc(2:Nlanc)
       call tql2(Nlanc,diag,subdiag,Z,ierr)
       if(verb)then
          print *,'---> lowest eigenvalue  <---'
          write(*,*)"E_lowest    = ",diag(1)
          open(10,file="lanc_eigenvals.dat")
          do i=1,Nlanc
             write(10,*)i,diag(i)
          enddo
          close(10)
       endif
       if(nlanc >= ncheck_)then
          esave(nlanc-(Ncheck_-1))=diag(1)
          if(nlanc >= (Ncheck_+1))then
             diff=esave(Nlanc-(Ncheck_-1))-esave(Nlanc-(Ncheck_-1)-1)
             if(verb)write(*,*)'test deltaE = ',diff
             if(abs(diff).le.threshold_)exit lanc_loop
          endif
       endif
    enddo lanc_loop
    if(verb)write(*,*)'lanczos deltaE  = ',diff
    if(nlanc==nitermax)print*,"LANCZOS_SIMPLE: reach Nitermax"
    !
    !============== END LANCZOS LOOP ======================
    !
    diag=0.d0 ; subdiag=0.d0 ; Z=0.d0
    forall(i=1:Nlanc)Z(i,i)=1.d0
    diag(1:Nlanc)    = alanc(1:Nlanc)
    subdiag(2:Nlanc) = blanc(2:Nlanc)
    call tql2(Nlanc,diag,subdiag,Z,ierr)
    !
    !Get the Eigenvalues:
    egs = diag(1)
    !
    !Get the Eigenvector:
    vin =vect
    vout=0.d0
    vect=0.d0
    do iter=1,nlanc
       call lanczos_plain_iteration_d(iter,vin,vout,alanc(iter),blanc(iter))
       vect = vect + vin*Z(iter,1)
    end do
    norm=sqrt(dot_product(vect,vect))
    vect=vect/norm
  end subroutine lanczos_plain_d
  !
  subroutine lanczos_plain_c(Hprod,Ns,Nitermax,Egs,Vect,Nlanc,iverbose,threshold,ncheck)
    procedure(lanc_htimesv_c)            :: Hprod
    integer                              :: ns,nitermax
    real(8)                              :: egs
    complex(8),dimension(ns)             :: vect
    complex(8),dimension(ns)             :: vin,vout
    integer                              :: iter,nlanc
    real(8),dimension(nitermax+1)        :: alanc,blanc
    real(8),dimension(Nitermax,Nitermax) :: Z
    real(8),dimension(Nitermax)          :: diag,subdiag,esave
    real(8)                              :: a_,b_,norm,diff,ran(2)
    integer                              :: i,ierr
    real(8),optional                     :: threshold
    integer,optional                     :: ncheck
    logical,optional                     :: iverbose
    if(present(iverbose))verb=iverbose
    if(present(threshold))threshold_=threshold
    if(present(ncheck))ncheck_=ncheck
    if(associated(cp_hprod))nullify(cp_hprod)
    cp_hprod=>Hprod
    if(.not.associated(cp_hprod))then
       print*,"LANCZOS_PLAIN: cp_hprod is not set. call lanczos_plain_set_htimesv"
       stop
    endif
    norm=dot_product(vect,vect)
    if(norm==0.d0)then
       do i=1,ns
          call random_number(ran)
          vect(i)=cmplx(ran(1),ran(2),8)
       enddo
       vect=vect/sqrt(dot_product(vect,vect))
       if(verb)write(*,*)"LANCZOS_PLAIN: random initial vector generated:"
    endif
    !
    !============= LANCZOS LOOP =====================
    !
    vin = vect
    vout= zero
    alanc=0.d0
    blanc=0.d0
    nlanc=0
    lanc_loop: do iter=1,Nitermax
       if(verb)then
          print*,""
          write(*,*)"Lanczos iteration:",iter
       endif
       call lanczos_plain_iteration_c(iter,vin,vout,a_,b_)
       if(abs(b_)<threshold_)exit lanc_loop
       if(verb)print*,"alanc,blanc=",a_,b_
       nlanc=nlanc+1
       alanc(iter)  =a_
       blanc(iter+1)=b_   
       diag=0.d0 ; subdiag=0.d0 ; Z=0.d0
       forall(i=1:Nlanc)Z(i,i)=1.d0
       diag(1:Nlanc)    = alanc(1:Nlanc)
       subdiag(2:Nlanc) = blanc(2:Nlanc)
       call tql2(Nlanc,diag,subdiag,Z,ierr)
       if(verb)then
          print *,'---> lowest eigenvalue  <---'
          write(*,*)"E_lowest    = ",diag(1)
          open(10,file="lanc_eigenvals.dat")
          do i=1,Nlanc
             write(10,*)i,diag(i)
          enddo
          close(10)
       endif
       if(nlanc >= ncheck_)then
          esave(nlanc-(Ncheck_-1))=diag(1)
          if(nlanc >= (Ncheck_+1))then
             diff=esave(Nlanc-(Ncheck_-1))-esave(Nlanc-(Ncheck_-1)-1)
             if(verb)write(*,*)'test deltaE = ',diff
             if(abs(diff).le.threshold_)exit lanc_loop
          endif
       endif
    enddo lanc_loop
    if(verb)write(*,*)'lanczos deltaE  = ',diff
    if(nlanc==nitermax)print*,"LANCZOS_SIMPLE: reach Nitermax"
    !
    !============== END LANCZOS LOOP ======================
    !
    diag=0.d0 ; subdiag=0.d0 ; Z=0.d0
    forall(i=1:Nlanc)Z(i,i)=1.d0
    diag(1:Nlanc)    = alanc(1:Nlanc)
    subdiag(2:Nlanc) = blanc(2:Nlanc)
    call tql2(Nlanc,diag,subdiag,Z,ierr)
    !
    !Get the Eigenvalues:
    egs = diag(1)
    !
    !Get the Eigenvector:
    vin =vect
    vout=0.d0
    vect=0.d0
    do iter=1,nlanc
       call lanczos_plain_iteration_c(iter,vin,vout,alanc(iter),blanc(iter))
       vect = vect + vin*Z(iter,1)
    end do
    norm=sqrt(dot_product(vect,vect))
    vect=vect/norm
  end subroutine lanczos_plain_c


#ifdef _MPI
  subroutine p_lanczos_plain_d(MpiComm,Hprod,Ns,Nitermax,Egs,Vect,Nlanc,iverbose,threshold,ncheck)
    integer                              :: MpiComm
    procedure(lanc_htimesv_d)            :: Hprod
    integer                              :: ns,nitermax
    real(8)                              :: egs
    real(8),dimension(ns)                :: vect
    real(8),dimension(ns)                :: vin,vout
    integer                              :: iter,nlanc
    real(8),dimension(nitermax+1)        :: alanc,blanc
    real(8),dimension(Nitermax,Nitermax) :: Z
    real(8),dimension(Nitermax)          :: diag,subdiag,esave
    real(8)                              :: a_,b_,norm,diff
    integer                              :: i,ierr
    real(8),optional                     :: threshold
    integer,optional                     :: ncheck
    logical,optional                     :: iverbose
    !
    lanc_mpi_comm=MpiComm
    lanc_mpi_rank=Mpi_Get_rank(Lanc_Mpi_Comm)
    lanc_mpi_size=Mpi_Get_size(Lanc_Mpi_Comm)
    lanc_mpi_master=MPI_Get_master(Lanc_Mpi_Comm)
    !
    if(present(iverbose))verb=iverbose
    if(present(threshold))threshold_=threshold
    if(present(ncheck))ncheck_=ncheck
    if(associated(dp_hprod))nullify(dp_hprod)
    dp_hprod=>Hprod
    if(.not.associated(dp_hprod))&
         stop "LANCZOS_PLAIN: dp_hprod is not set. call lanczos_plain_set_htimesv"
    norm=dot_product(vect,vect)
    if(norm==0.d0)then
       call random_number(vect)
       vect=vect/sqrt(dot_product(vect,vect))
       if(verb.AND.lanc_mpi_master)write(*,*)"LANCZOS_PLAIN: random initial vector generated:"
    endif
    !
    !============= LANCZOS LOOP =====================
    !
    vin = vect
    vout= 0.d0
    alanc=0.d0
    blanc=0.d0
    nlanc=0
    lanc_loop: do iter=1,Nitermax
       if(verb.AND.lanc_mpi_master)then
          print*,""
          write(*,*)"Lanczos iteration:",iter
       endif
       call lanczos_plain_iteration_d(iter,vin,vout,a_,b_)
       if(abs(b_)<threshold_)exit lanc_loop
       if(verb.AND.lanc_mpi_master)print*,"alanc,blanc=",a_,b_
       nlanc=nlanc+1
       alanc(iter) = a_ ; blanc(iter+1) = b_
       diag = 0.d0 ; subdiag = 0.d0 ; Z = 0.d0
       forall(i=1:Nlanc)Z(i,i)=1.d0
       diag(1:Nlanc)    = alanc(1:Nlanc)
       subdiag(2:Nlanc) = blanc(2:Nlanc)
       call tql2(Nlanc,diag,subdiag,Z,ierr)
       if(verb.AND.lanc_mpi_master)then
          print *,'---> lowest eigenvalue  <---'
          write(*,*)"E_lowest    = ",diag(1)
          open(10,file="lanc_eigenvals.dat")
          do i=1,Nlanc
             write(10,*)i,diag(i)
          enddo
          close(10)
       endif
       if(nlanc >= ncheck_)then
          esave(nlanc-(Ncheck_-1))=diag(1)
          if(nlanc >= (Ncheck_+1))then
             diff=esave(Nlanc-(Ncheck_-1))-esave(Nlanc-(Ncheck_-1)-1)
             if(verb.AND.lanc_mpi_master)write(*,*)'test deltaE = ',diff
             if(abs(diff).le.threshold_)exit lanc_loop
          endif
       endif
    enddo lanc_loop
    if(verb.AND.lanc_mpi_master)write(*,*)'Lanczos deltaE = ',diff
    if(nlanc==nitermax)print*,"LANCZOS_SIMPLE: reach Nitermax"
    !
    !============== END LANCZOS LOOP ======================
    !
    diag=0.d0 ; subdiag=0.d0 ; Z=0.d0
    forall(i=1:Nlanc)Z(i,i)=1.d0
    diag(1:Nlanc)    = alanc(1:Nlanc)
    subdiag(2:Nlanc) = blanc(2:Nlanc)
    call tql2(Nlanc,diag,subdiag,Z,ierr)
    !
    !Get the Eigenvalues:
    egs = diag(1)
    !
    !Get the Eigenvector:
    vin =vect
    vout=0.d0
    vect=0.d0
    do iter=1,nlanc
       call lanczos_plain_iteration_d(iter,vin,vout,alanc(iter),blanc(iter))
       vect = vect + vin*Z(iter,1)
    end do
    norm=sqrt(dot_product(vect,vect))
    vect=vect/norm
  end subroutine p_lanczos_plain_d
  !
  subroutine p_lanczos_plain_c(MpiComm,Hprod,Ns,Nitermax,Egs,Vect,Nlanc,iverbose,threshold,ncheck)
    integer                              :: MpiComm
    procedure(lanc_htimesv_c)            :: Hprod
    integer                              :: ns,nitermax
    real(8)                              :: egs
    complex(8),dimension(ns)             :: vect
    complex(8),dimension(ns)             :: vin,vout
    integer                              :: iter,nlanc
    real(8),dimension(nitermax+1)        :: alanc,blanc
    real(8),dimension(Nitermax,Nitermax) :: Z
    real(8),dimension(Nitermax)          :: diag,subdiag,esave
    real(8)                              :: a_,b_,norm,diff,ran(2)
    integer                              :: i,ierr
    real(8),optional                     :: threshold
    integer,optional                     :: ncheck
    logical,optional                     :: iverbose
    !
    lanc_mpi_comm=MpiComm
    lanc_mpi_rank=Mpi_Get_rank(Lanc_Mpi_Comm)
    lanc_mpi_size=Mpi_Get_size(Lanc_Mpi_Comm)
    lanc_mpi_master=MPI_Get_master(Lanc_Mpi_Comm)
    !
    if(present(iverbose))verb=iverbose
    if(present(threshold))threshold_=threshold
    if(present(ncheck))ncheck_=ncheck
    if(associated(cp_hprod))nullify(cp_hprod)
    cp_hprod=>Hprod
    if(.not.associated(cp_hprod))&
         stop "LANCZOS_PLAIN: cp_hprod is not set. call lanczos_plain_set_htimesv"
    norm=dot_product(vect,vect)
    if(norm==0.d0)then
       do i=1,ns
          call random_number(ran)
          vect(i)=cmplx(ran(1),ran(2),8)
       enddo
       vect=vect/sqrt(dot_product(vect,vect))
       if(verb.AND.lanc_mpi_master)write(*,*)"LANCZOS_PLAIN: random initial vector generated:"
    endif
    !
    !============= LANCZOS LOOP =====================
    !
    vin = vect
    vout= zero
    alanc=0.d0
    blanc=0.d0
    nlanc=0
    lanc_loop: do iter=1,Nitermax
       if(verb.AND.lanc_mpi_master)then
          print*,""
          write(*,*)"Lanczos iteration:",iter
       endif
       call lanczos_plain_iteration_c(iter,vin,vout,a_,b_)
       if(abs(b_)<threshold_)exit lanc_loop
       if(verb.AND.lanc_mpi_master)print*,"alanc,blanc=",a_,b_
       nlanc=nlanc+1
       alanc(iter)  =a_
       blanc(iter+1)=b_   
       diag=0.d0 ; subdiag=0.d0 ; Z=0.d0
       forall(i=1:Nlanc)Z(i,i)=1.d0
       diag(1:Nlanc)    = alanc(1:Nlanc)
       subdiag(2:Nlanc) = blanc(2:Nlanc)
       call tql2(Nlanc,diag,subdiag,Z,ierr)
       if(verb.AND.lanc_mpi_master)then
          print *,'---> lowest eigenvalue  <---'
          write(*,*)"E_lowest    = ",diag(1)
          open(10,file="lanc_eigenvals.dat")
          do i=1,Nlanc
             write(10,*)i,diag(i)
          enddo
          close(10)
       endif
       if(nlanc >= ncheck_)then
          esave(nlanc-(Ncheck_-1))=diag(1)
          if(nlanc >= (Ncheck_+1))then
             diff=esave(Nlanc-(Ncheck_-1))-esave(Nlanc-(Ncheck_-1)-1)
             if(verb.AND.lanc_mpi_master)write(*,*)'test deltaE = ',diff
             if(abs(diff).le.threshold_)exit lanc_loop
          endif
       endif
    enddo lanc_loop
    if(verb.AND.lanc_mpi_master)write(*,*)'lanczos deltaE = ',diff
    if(nlanc==nitermax)print*,"LANCZOS_SIMPLE: reach Nitermax"
    !
    !============== END LANCZOS LOOP ======================
    !
    diag=0.d0 ; subdiag=0.d0 ; Z=0.d0
    forall(i=1:Nlanc)Z(i,i)=1.d0
    diag(1:Nlanc)    = alanc(1:Nlanc)
    subdiag(2:Nlanc) = blanc(2:Nlanc)
    call tql2(Nlanc,diag,subdiag,Z,ierr)
    !
    !Get the Eigenvalues:
    egs = diag(1)
    !
    !Get the Eigenvector:
    vin =vect
    vout=0.d0
    vect=0.d0
    do iter=1,nlanc
       call lanczos_plain_iteration_c(iter,vin,vout,alanc(iter),blanc(iter))
       vect = vect + vin*Z(iter,1)
    end do
    norm=sqrt(dot_product(vect,vect))
    vect=vect/norm
  end subroutine p_lanczos_plain_c
#endif






  !##################################################################
  !        LANCZOS METHOD for TRI-DIAGONALIZATION OF A MATRIX
  !                          (KRYLOV METHOD)
  !##################################################################
  !---------------------------------------------------------------------
  !Purpose: use simple Lanczos to tri-diagonalize a matrix H (defined 
  ! in the H*v function).
  !---------------------------------------------------------------------
  subroutine lanczos_plain_tridiag_d(hprod,vin,alanc,blanc,nitermax,iverbose,threshold)
    real(8),dimension(:),intent(inout)        :: vin
    real(8),dimension(size(vin))              :: vout
    real(8),dimension(nitermax),intent(inout) :: alanc
    real(8),dimension(nitermax),intent(inout) :: blanc
    procedure(lanc_htimesv_d)                 :: Hprod
    integer                                   :: nitermax
    integer                                   :: iter
    real(8)                                   :: a_,b_
    real(8),optional                          :: threshold
    logical,optional                          :: iverbose
    if(present(iverbose))verb=iverbose
    if(present(threshold))threshold_=threshold
    if(associated(dp_hprod))nullify(dp_hprod)
    dp_hprod=>Hprod
    a_=0.d0
    b_=0.d0
    vout=0.d0
    do iter=1,nitermax
       call lanczos_plain_iteration_d(iter,vin,vout,a_,b_)
       if(verb.AND.lanc_mpi_master)print*,iter,a_,b_
       alanc(iter)=a_
       if(abs(b_)<threshold_)exit
       if(iter<nitermax)blanc(iter+1)=b_
    enddo
    if(iter==nitermax.AND.lanc_mpi_master)print*,"LANCZOS_SIMPLE: reach Nitermax"
  end subroutine lanczos_plain_tridiag_d
  !
  subroutine lanczos_plain_tridiag_c(hprod,vin,alanc,blanc,nitermax,iverbose,threshold)
    complex(8),dimension(:),intent(inout)     :: vin
    complex(8),dimension(size(vin))           :: vout
    real(8),dimension(nitermax),intent(inout) :: alanc
    real(8),dimension(nitermax),intent(inout) :: blanc
    procedure(lanc_htimesv_c)                 :: Hprod
    integer                                   :: nitermax
    integer                                   :: iter
    real(8)                                   :: a_,b_
    real(8),optional                          :: threshold
    logical,optional                          :: iverbose
    if(present(iverbose))verb=iverbose
    if(present(threshold))threshold_=threshold
    if(associated(cp_hprod))nullify(cp_hprod)
    cp_hprod=>Hprod
    a_=0.d0
    b_=0.d0
    vout=zero
    do iter=1,nitermax
       call lanczos_plain_iteration_c(iter,vin,vout,a_,b_)
       if(verb.AND.lanc_mpi_master)print*,iter,a_,b_
       alanc(iter)=a_
       if(iter<nitermax)blanc(iter+1)=b_
       if(abs(b_)<threshold_)exit
    enddo
    if(iter==nitermax.AND.lanc_mpi_master)print*,"LANCZOS_SIMPLE: reach Nitermax"
  end subroutine lanczos_plain_tridiag_c







  !##################################################################
  !                     LANCZOS ITERATION
  !##################################################################
  !---------------------------------------------------------------------
  !Purpose: plain homebrew lanczos iteration (no orthogonalization)
  !note: the a,b variables are real, even in the complex matrix case
  !to understand why check out the Gollub-Van Loan textbook.
  !a it is easy: hermiticity->diag\in\RRR
  !b: is fixed by requiring |b|^2 = <v,v> thus you can only fix the 
  !the absolute value. A lemma shows that the phase can be chosen 
  !identically zero
  !---------------------------------------------------------------------
  subroutine lanczos_plain_iteration_d(iter,vin,vout,a,b)
    real(8),dimension(:),intent(inout)            :: vin
    real(8),dimension(size(vin)),intent(inout)    :: vout
    real(8),dimension(size(vin))                  :: tmp
    real(8),intent(inout)                         :: a,b
    integer                                       :: iter,ndim,nloc
    real(8)                                       :: norm
    ndim=size(vin)
    if(iter==1)then
       norm=sqrt(dot_product(vin,vin))
       if(norm==0.d0)stop "lanczos_plain_iteration: norm =0!!"
       vin=vin/norm
       b=0.d0
    end if
    call dp_Hprod(ndim,vin,tmp) !Ndim = size(vin)
    tmp = tmp-b*vout
    a   = dot_product(vin,tmp)
    tmp = tmp-a*vin
    b   = sqrt(dot_product(tmp,tmp))
    vout= vin
    vin = tmp/b
  end subroutine lanczos_plain_iteration_d
  !
  subroutine lanczos_plain_iteration_c(iter,vin,vout,a,b)
    complex(8),dimension(:),intent(inout)         :: vin
    complex(8),dimension(size(vin)),intent(inout) :: vout
    complex(8),dimension(size(vin))               :: tmp
    real(8),intent(inout)                         :: a,b
    integer                                       :: iter,ndim,nloc
    real(8)                                       :: norm
    ndim=size(vin)
    if(iter==1)then
       norm=sqrt(dot_product(vin,vin))
       if(norm==0.d0)stop "lanczos_plain_iteration: norm =0!!"
       vin=vin/norm
       b=0.d0
    end if
    call cp_Hprod(ndim,vin,tmp) !Ndim = size(vin)
    tmp = tmp-b*vout
    a   = dot_product(vin,tmp)
    tmp = tmp-a*vin
    b   = sqrt(dot_product(tmp,tmp))
    vout= vin
    vin = tmp/b
  end subroutine lanczos_plain_iteration_c
  !










  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !++++++++++++++++++COMPUTATIONAL ROUTINE: TQL2++++++++++++++++++++++++ 
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !---------------------------------------------------------------------
  ! PURPOSE computes all eigenvalues/vectors, real symmetric tridiagonal matrix.
  !    This subroutine finds the eigenvalues and eigenvectors of a symmetric
  !    tridiagonal matrix by the QL method.  The eigenvectors of a full
  !    symmetric matrix can also be found if TRED2 has been used to reduce this
  !    full matrix to tridiagonal form.
  !  Parameters:
  !    Input, integer ( kind = 4 ) N, the order of the matrix.
  !
  !    Input/output, real ( kind = 8 ) D(N).  On input, the diagonal elements of
  !    the matrix.  On output, the eigenvalues in ascending order.  If an error
  !    exit is made, the eigenvalues are correct but unordered for indices
  !    1,2,...,IERR-1.
  !
  !    Input/output, real ( kind = 8 ) E(N).  On input, E(2:N) contains the
  !    subdiagonal elements of the input matrix, and E(1) is arbitrary.
  !    On output, E has been destroyed.
  !
  !    Input, real ( kind = 8 ) Z(N,N).  On input, the transformation matrix
  !    produced in the reduction by TRED2, if performed.  If the eigenvectors of
  !    the tridiagonal matrix are desired, Z must contain the identity matrix.
  !    On output, Z contains the orthonormal eigenvectors of the symmetric
  !    tridiagonal (or full) matrix.  If an error exit is made, Z contains
  !    the eigenvectors associated with the stored eigenvalues.
  !
  !    Output, integer ( kind = 4 ) IERR, error flag.
  !    0, normal return,
  !    J, if the J-th eigenvalue has not been determined after
  !    30 iterations.
  !
  !---------------------------------------------------------------------
  subroutine tql2 ( n, d, e, z, ierr )
    integer :: n
    real(8) :: c
    real(8) :: c2
    real(8) :: c3
    real(8) :: d(n)
    real(8) :: dl1
    real(8) :: e(n)
    real(8) :: el1
    real(8) :: f
    real(8) :: g
    real(8) :: h
    integer ( kind = 4 ) i
    integer ( kind = 4 ) ierr
    integer ( kind = 4 ) ii
    integer ( kind = 4 ) j
    integer ( kind = 4 ) k
    integer ( kind = 4 ) l
    integer ( kind = 4 ) l1
    integer ( kind = 4 ) l2
    integer ( kind = 4 ) m
    integer ( kind = 4 ) mml
    real(8) :: p
    real(8) :: r
    real(8) :: s
    real(8) :: s2
    real(8) :: tst1
    real(8) :: tst2
    real(8) :: z(n,n)
    ierr = 0
    if ( n == 1 ) then
       return
    end if
    do i = 2, n
       e(i-1) = e(i)
    end do
    f = 0.0D+00
    tst1 = 0.0D+00
    e(n) = 0.0D+00
    do l = 1, n
       j = 0
       h = abs ( d(l) ) + abs ( e(l) )
       tst1 = max ( tst1, h )
       !
       !  Look for a small sub-diagonal element.
       !
       do m = l, n
          tst2 = tst1 + abs ( e(m) )
          if ( tst2 == tst1 ) then
             exit
          end if
       end do
       if ( m == l ) then
          go to 220
       end if
130    continue
       if ( 30 <= j ) then
          ierr = l
          return
       end if
       j = j + 1
       !
       !  Form shift.
       !
       l1 = l + 1
       l2 = l1 + 1
       g = d(l)
       p = ( d(l1) - g ) / ( 2.0D+00 * e(l) )
       r = pythag ( p, 1.0D+00 )
       d(l) = e(l) / ( p + sign ( r, p ) )
       d(l1) = e(l) * ( p + sign ( r, p ) )
       dl1 = d(l1)
       h = g - d(l)
       d(l2:n) = d(l2:n) - h
       f = f + h
       !
       !  QL transformation.
       !
       p = d(m)
       c = 1.0D+00
       c2 = c
       el1 = e(l1)
       s = 0.0D+00
       mml = m - l
       do ii = 1, mml
          c3 = c2
          c2 = c
          s2 = s
          i = m - ii
          g = c * e(i)
          h = c * p
          r = pythag ( p, e(i) )
          e(i+1) = s * r
          s = e(i) / r
          c = p / r
          p = c * d(i) - s * g
          d(i+1) = h + s * ( c * g + s * d(i) )
          !
          !  Form vector.
          !
          do k = 1, n
             h = z(k,i+1)
             z(k,i+1) = s * z(k,i) + c * h
             z(k,i) = c * z(k,i) - s * h
          end do
       end do
       p = - s * s2 * c3 * el1 * e(l) / dl1
       e(l) = s * p
       d(l) = c * p
       tst2 = tst1 + abs ( e(l) )
       if ( tst2 > tst1 ) then
          go to 130
       end if
220    continue
       d(l) = d(l) + f
    end do
    !
    !  Order eigenvalues and eigenvectors.
    !
    do ii = 2, n
       i = ii - 1
       k = i
       p = d(i)
       do j = ii, n
          if ( d(j) < p ) then
             k = j
             p = d(j)
          end if
       end do
       if ( k /= i ) then
          d(k) = d(i)
          d(i) = p
          do j = 1, n
             call r8_swap ( z(j,i), z(j,k) )
          end do
       end if
    end do
    return
  end subroutine tql2


  !---------------------------------------------------------------------
  ! PURPOSE: computes SQRT ( A * A + B * B ) carefully.
  !    The formula
  !    PYTHAG = sqrt ( A * A + B * B )
  !    is reasonably accurate, but can fail if, for example, A**2 is larger
  !    than the machine overflow.  The formula can lose most of its accuracy
  !    if the sum of the squares is very large or very small.
  !  Parameters:
  !    Input, real(8) :: A, B, the two legs of a right triangle.
  !    Output, real(8) :: PYTHAG, the length of the hypotenuse.
  !---------------------------------------------------------------------
  function pythag ( a, b )
    implicit none
    real(8) :: a
    real(8) :: b
    real(8) :: p
    real(8) :: pythag
    real(8) :: r
    real(8) :: s
    real(8) :: t
    real(8) :: u
    p = max ( abs ( a ), abs ( b ) )
    if ( p /= 0.0D+00 ) then
       r = ( min ( abs ( a ), abs ( b ) ) / p )**2
       do
          t = 4.0D+00 + r
          if ( t == 4.0D+00 ) then
             exit
          end if
          s = r / t
          u = 1.0D+00 + 2.0D+00 * s
          p = u * p
          r = ( s / u )**2 * r
       end do
    end if
    pythag = p
    return
  end function pythag

  !---------------------------------------------------------------------
  ! PURPOSE: swaps two R8's.
  !  Parameters:
  !    Input/output, real(8) :: X, Y.  On output, the values of X and
  !    Y have been interchanged.
  !---------------------------------------------------------------------
  subroutine r8_swap ( x, y )
    real(8) :: x
    real(8) :: y
    real(8) :: z
    z = x
    x = y
    y = z
    return
  end subroutine r8_swap





#ifdef _MPI
  function MPI_Get_size(comm) result(size)
    integer :: comm
    integer :: size,ierr
    call MPI_Comm_size(comm,size,ierr)
  end function MPI_Get_size



  function MPI_Get_rank(comm) result(rank)
    integer :: comm
    integer :: rank,ierr
    call MPI_Comm_rank(comm,rank,ierr)
  end function MPI_Get_rank



  function MPI_Get_master(comm) result(master)
    integer :: comm
    logical :: master
    integer :: rank,ierr
    call MPI_Comm_rank(comm,rank,ierr)
    master=.false.
    if(rank==0)master=.true.
  end function MPI_Get_master
#endif

end module SF_SP_LINALG
