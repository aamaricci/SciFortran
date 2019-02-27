subroutine lanczos_parpack_d(MpiComm,MatVec,Ns,Neigen,Nblock,Nitermax,eval,evec,which,v0,tol,iverbose,vrandom)
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
  logical,optional             :: vrandom
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
  logical                      :: verb,vran
  integer                      :: i
  real(8)                      :: sigma,norm,norm_tmp
  real(8)                      :: tol_
  character                    :: bmat  
  character(len=2)             :: which_
  real(8),external             :: dnrm2
  !MPI
  logical                      :: mpi_master
  !
  if(MpiComm == MPI_COMM_NULL)return
  !
  which_ = 'SA'   ; if(present(which))which_=which
  tol_   = 0d0    ; if(present(tol))tol_=tol
  verb   = .false.; if(present(iverbose))verb=iverbose
  vran   = .true. ;if(present(vrandom))vran=vrandom
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
        call random_seed(size=nrandom)
        if(allocated(seed_random))deallocate(seed_random)
        allocate(seed_random(nrandom))
        seed_random=1234567
        call random_seed(put=seed_random)
        deallocate(seed_random)
        call random_number(vec)
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
     !  Perform matrix vector multiplication
     !    y <--- OP*x ; workd(ipntr(1))=input, workd(ipntr(2))=output
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
     if(any(shape(evec)/=[ldv,Neigen]))stop "arpack_mpi ERROR: evec has wrong dimension. Reduce=T"
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
  call Barrier_MPI(MpiComm)
  deallocate(ax,resid,workd,v,d,workl,select)
end subroutine lanczos_parpack_d
