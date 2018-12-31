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
  tol_=0d0;if(present(tol))tol_=tol
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
  if(info/=0)then
     write(*,'(a,i6)')'Warning/Error in DSAUPD, info = ', info
     select case(info)
     case(1)
        write(*,'(a)')'Maximum number of iterations reached.'
        write(*,'(a)')'All possible eigenvalues of OP has been found. '
        write(*,'(a)')'IPARAM(5) returns the number of wanted converged Ritz values.'
        write(*,'(a,I0)')'IPARAM(5) = ',Iparam(5)              
     case(3)
        write(*,'(a)') ' No shifts could be applied during implicit '&
             //'Arnoldi update, try increasing NCV.'
        stop
     case(-1)
        write(*,'(a)')'N must be positive.'
        stop
     case(-2)
        write(*,'(a)')'NEV must be positive.'
        stop
     case(-3)
        write(*,'(a)')'NCV must be greater than NEV and less than or equal to N.'
        stop
     case(-4)
        write(*,'(a)')'The maximum number of Arnoldi update iterations allowed must be greater than zero.'
        stop
     case(-5)
        write(*,'(a)')'WHICH must be one of LM, SM, LA, SA or BE.'
        stop
     case(-6)
        write(*,'(a)')'BMAT must be one of I or G.'
        stop
     case(-7)
        write(*,'(a)')'Length of private work array WORKL is not sufficient.'
        stop
     case(-8)
        write(*,'(a)')'Error return from trid. eigenvalue calculation; Informatinal error from LAPACK routine dsteqr .'
        stop
     case(-9)
        write(*,'(a)')'Starting vector is zero.'
        stop
     case(-10)
        write(*,'(a)')'IPARAM(7) must be 1,2,3,4,5.'
        stop
     case(-11)
        write(*,'(a)')'IPARAM(7) = 1 and BMAT = G are incompatable.'
        stop
     case(-12)
        write(*,'(a)')'IPARAM(1) must be equal to 0 or 1.'
        stop
     case(-13)
        write(*,'(a)')'NEV and WHICH = BE are incompatable.'
        stop
     case(-9999)
        write(*,'(a)')'Could not build an Arnoldi factorization.'
        write(*,'(a)')'IPARAM(5) returns the size of the current Arnoldi factorization.'
        write(*,'(a)')'The user is advised to check that enough workspace and array storage has been allocated.'
        stop
     end select
  else
     rvec = .true.
     call dseupd(rvec,'All',select,d,v,ldv,sigma,bmat,n,which_,&
          nev,tol_,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,ierr)
     do j=1,neigen
        eval(j)=d(j,1)
        do i=1,ns
           evec(i,j)=v(i,j)
        enddo
     enddo
     !
     !=========================================================================
     !  Compute the residual norm
     !    ||  A*x - lambda*x ||
     !  for the NCONV accurately computed eigenvalues and 
     !  eigenvectors.  (iparam(5) indicates how many are 
     !  accurate to the requested tolerance)
     if(ierr/=0)then
        write(*,'(a,i6)')'Error with DSEUPD, IERR = ',ierr
        write(*,'(a)')'Check the documentation of DSEUPD.'
        stop
     else
        nconv =  iparam(5)
     end if
     !
     if(verb)then
        write(*,'(a)') ''
        write(*,'(a)') 'ARPACK::'
        write(*,'(a)') ''
        write(*,'(a,i6)') '  Size of the matrix is:                      ', n
        write(*,'(a,i6)') '  Number of Ritz values requested is:         ', nev
        write(*,'(a,i6)') '  Number of Arnoldi vectors generated is:     ', ncv
        write(*,'(a)')    '  Portion of the spectrum:                        '//trim(which_)
        write(*,'(a,i6)') '  Number of converged Ritz values is:         ', nconv
        write(*,'(a,i6)') '  Number of Implicit Arnoldi iterations is:   ', iparam(3)
        write(*,'(a,i6)') '  Number of OP*x is:                          ', iparam(9)
        write(*,'(a,ES14.6)') '  The convergence criterion is:           ', tol
     end if
     !
     !if(nconv == 0) stop "ARPACK:: no converged eigenvalues have been found."
     !
  endif
  deallocate(ax,d,resid,workl,workd,v,select)
end subroutine lanczos_arpack_d
