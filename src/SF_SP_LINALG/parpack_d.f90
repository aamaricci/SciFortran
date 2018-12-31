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
  real(8)                      :: evec(:,:) ![Nloc,Neigen]
  character(len=2),optional    :: which
  real(8),optional             :: v0(size(evec,1))
  real(8),optional             :: tol
  logical,optional             :: iverbose
  !Dimensions:
  integer                      :: maxn,maxnev,maxncv,ldv
  integer                      :: n,nconv,ncv,nev
  !Arrays:
  real(8),allocatable          :: evec_tmp(:) ![Nloc] see above
  real(8),allocatable          :: ax(:),d(:,:)
  real(8),allocatable          :: resid(:),vec(:)
  real(8),allocatable          :: workl(:),workd(:)
  real(8),allocatable          :: v(:,:)
  logical,allocatable          :: select(:)
  integer                      :: iparam(11)
  integer                      :: ipntr(11)
  !Control Vars:
  integer                      :: ido,ierr,info,ishfts,j,lworkl,maxitr,mode1
  logical                      :: verb
  integer                      :: i
  real(8)                      :: sigma,norm,norm_tmp
  real(8)                      :: tol_
  character                    :: bmat  
  character(len=2)             :: which_
  real(8),external             :: dnrm2
  !MPI
  logical                      :: mpi_master
  !
  if(MpiComm /= MPI_COMM_NULL)then
     !
     which_ = 'SA'   ; if(present(which))which_=which
     tol_   = 1d-12  ; if(present(tol))tol_=tol
     verb   = .false.; if(present(iverbose))verb=iverbose
     !
     !MPI setup:
     mpi_master= Get_master_MPI(MpiComm)
     !
     maxn   = Ns
     maxnev = Neigen
     maxncv = Nblock
     ldv    = maxn
     if(maxncv>Ns)then
        maxncv=Ns
        print*,"PARPACK WARNING Ncv > Ns: reset block size to ",Ns
     endif
     !
     n      = maxn
     nev    = maxnev
     ncv    = maxncv
     bmat   = 'I'
     maxitr = Nitermax
     !=========================================================================
     ! Setup distribution of data to nodes:
     ldv = size(evec,1)            !ldv is the SMALL dimension
     !if(ldv>0.AND.ldv < ncv)stop "LANCZOS_PARPACK_D error: ldv < maxNblock. Unstable! Find a way to increase ldv... (less cpu?)"
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
        call random_seed(size=nrandom)
        if(allocated(seed_random))deallocate(seed_random)
        allocate(seed_random(nrandom))
        seed_random=1234567
        call random_seed(put=seed_random)
        !
        allocate(vec(ldv))
        call random_number(vec)
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
        !  Perform matrix vector multiplication
        !    y <--- OP*x ; workd(ipntr(1))=input, workd(ipntr(2))=output
        call MatVec(ldv,workd(ipntr(1)),workd(ipntr(2)))
     end do
     !
     !POST PROCESSING:
     if(info/=0)then
        if(mpi_master)then
           write(*,'(a,i6)')'Warning/Error in PDSAUPD, info = ', info
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
        endif
     else
        !
        call pdseupd (MpiComm,.true.,'All',select,d,v,ldv,sigma,bmat,&
             ldv,which_,nev,tol_,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,ierr)        
        do j=1,neigen
           eval(j)=d(j,1)
        enddo
        !
        if(any(shape(evec)/=[ldv,Neigen]))stop "arpack_mpi ERROR: evec has wrong dimension. Reduce=T"
        evec=0d0
        do j=1,neigen
           do i=1,ldv
              evec(i,j)=v(i,j)
           enddo
        enddo
        !
        !
        !=========================================================================
        !  Compute the residual norm
        !    ||  A*x - lambda*x ||
        !  for the NCONV accurately computed eigenvalues and 
        !  eigenvectors.  (iparam(5) indicates how many are 
        !  accurate to the requested tolerance)
        if(ierr/=0)then
           write(*,'(a,i6)')'Error with PDSEUPD, IERR = ',ierr
           write(*,'(a)')'Check the documentation of PDSEUPD.'
        else
           nconv =  iparam(5)
        end if
        !
        if(mpi_master.AND.verb)then
           write(*,'(a)') ''
           write(*,'(a)') 'ARPACK::'
           write(*,'(a)') ''
           write(*,'(a,i12)') '  Size of the matrix is:                      ', n
           write(*,'(a,i6)') '  Number of Ritz values requested is:         ', nev
           write(*,'(a,i6)') '  Number of Arnoldi vectors generated is:     ', ncv
           write(*,'(a)')    '  Portion of the spectrum:                        '//trim(which_)
           write(*,'(a,i6)') '  Number of converged Ritz values is:         ', nconv
           write(*,'(a,i6)') '  Number of Implicit Arnoldi iterations is:   ', iparam(3)
           write(*,'(a,i6)') '  Number of OP*x is:                          ', iparam(9)
           write(*,'(a,ES14.6)') '  The convergence criterion is:           ', tol
        end if
        !
        ! if(mpi_master.and.nconv==0.and.verb)stop "None of the required values was found."
     endif
     call Barrier_MPI(MpiComm)
     deallocate(ax,resid,workd,v,d,workl,select)
  endif
end subroutine lanczos_parpack_d
