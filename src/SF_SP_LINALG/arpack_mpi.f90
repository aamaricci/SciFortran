subroutine lanczos_parpack_d(MpiComm,MatVec,Ns,Neigen,Nblock,Nitermax,eval,evec,which,v0,tol,iverbose,mpi_reduce)
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
  real(8)                      :: evec(:,:) ![Ns,Neigen]if Reduce=F (default) OR [Nloc,Neigen] if Reduce=T
  character(len=2),optional    :: which
  real(8),optional             :: v0(ns)
  real(8),optional             :: tol
  logical,optional             :: iverbose,mpi_reduce
  !Dimensions:
  integer                      :: maxn,maxnev,maxncv,ldv
  integer                      :: n,nconv,ncv,nev
  !Arrays:
  real(8),allocatable          :: evec_tmp(:) ![Ns] OR [Nloc] see above
  real(8),allocatable          :: ax(:),d(:,:)
  real(8),allocatable          :: resid(:),vec(:)
  real(8),allocatable          :: workl(:),workd(:)
  real(8),allocatable          :: v(:,:)
  logical,allocatable          :: select(:)
  integer                      :: iparam(11)
  integer                      :: ipntr(11)
  !Control Vars:
  integer                      :: ido,ierr,info,ishfts,j,lworkl,maxitr,mode1
  logical                      :: verb,reduce_
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
  which_ = 'SA'   ; if(present(which))which_=which
  tol_   = 0d0    ; if(present(tol))tol_=tol
  verb   = .false.; if(present(iverbose))verb=iverbose
  reduce_= .true.; if(present(mpi_reduce))reduce_=mpi_reduce
  !
  !MPI setup:
  mpi_size  = Get_size_MPI(MpiComm)
  mpi_rank  = Get_rank_MPI(MpiComm)
  mpi_master= Get_master_MPI(MpiComm)
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
  !
  if(present(v0))then
     do i=1 + mpi_rank*mpiQ, (mpi_rank+1)*mpiQ + mpiR
        resid(i-mpi_rank*mpiQ)=v0(i)
     enddo
  else
     call random_seed(size=nrandom)
     if(allocated(seed_random))deallocate(seed_random)
     allocate(seed_random(nrandom))
     seed_random=1234567
     call random_seed(put=seed_random)
     !
     allocate(vec(n))
     call random_number(vec)
     vec=vec/sqrt(dot_product(vec,vec))
     do i=1 + mpi_rank*mpiQ, (mpi_rank+1)*mpiQ + mpiR
        resid(i-mpi_rank*mpiQ)=vec(i)
     enddo
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
     if(ido/=-1.AND.ido/=1)exit
     !  Perform matrix vector multiplication
     !    y <--- OP*x ; workd(ipntr(1))=input, workd(ipntr(2))=output
     call MatVec(ldv,workd(ipntr(1)),workd(ipntr(2)))
  end do
  !
  if(info<0)then
     if(mpi_master)write(*,'(a,i6)')'Fatal Error in PDSAUPD, info = ', info
     stop
  else
     call pdseupd (MpiComm,.true.,'All',select,d,v,ldv,sigma,bmat,&
          ldv,which_,nev,tol_,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,ierr)
     do j=1,neigen
        eval(j)=d(j,1)
     enddo
     !
     ! if(present(evec))then
     select case(reduce_)
     case(.false.)
        if(any(shape(evec)/=[ldv,Neigen]))stop "arpack_mpi ERROR: evec has wrong dimension. Reduce=T"
        evec=0d0
        do j=1,neigen
           do i=1,mpiQ+mpiR
              evec(i,j)=v(i,j)
           enddo
        enddo
        !
     case(.true.)
        if(any(shape(evec)/=[Ns,Neigen]))stop "arpack_mpi ERROR: evec has wrong dimension. Reduce=F"
        evec=0d0
        allocate(evec_tmp(Ns))
        do j=1,neigen
           evec_tmp=0d0
           do i=1 + mpi_rank*mpiQ , (mpi_rank+1)*mpiQ+mpiR
              evec_tmp(i)=v(i-mpi_rank*mpiQ,j)
           enddo
           call MPI_Allreduce(evec_tmp,evec(:,j),Ns,MPI_DOUBLE_PRECISION,MPI_SUM,MpiComm,mpi_ierr)
        enddo
        deallocate(evec_tmp)
        !
     end select
     ! endif
     !
     !=========================================================================
     !  Compute the residual norm
     !    ||  A*x - lambda*x ||
     !  for the NCONV accurately computed eigenvalues and 
     !  eigenvectors.  (iparam(5) indicates how many are 
     !  accurate to the requested tolerance)
     if(ierr/=0)then
        write(*,'(a,i6)')'Error with PDSEUPD, IERR = ',ierr
        write(*,'(a)')'Check the documentation of SSEUPD.'
     else
        nconv =  iparam(5)
        !
        ! do j = 1, nconv
        !    call hprod( n, v(1,j), ax )
        !    call zaxpy( n, -d(j), v(1,j), 1, ax, 1 )
        !    rd(j,1) = dble (d(j))
        !    rd(j,2) = dimag (d(j))
        !    rd(j,3) = dznrm2 (n, ax, 1)
        !    rd(j,3) = rd(j,3) / dlapy2 (rd(j,1),rd(j,2))
        ! end do
        ! if(mpi_master.AND.verb)call dmout(6,nconv,3,rd,maxncv,-6,'Ritz values and relative residuals')
     end if

     if(mpi_master.AND.verb)then
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










subroutine lanczos_parpack_c(MpiComm,MatVec,Ns,Neigen,Nblock,Nitermax,eval,evec,which,v0,tol,iverbose,mpi_reduce)
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
  complex(8)                        :: evec(:,:) ![Ns,Neigen]if Reduce=F (default) OR [Nloc,Neigen] if Reduce=T
  character(len=2),optional         :: which
  complex(8),optional               :: v0(ns)
  real(8),optional                  :: tol
  logical,optional                  :: iverbose,mpi_reduce
  !Dimensions:
  integer                           :: maxn,maxnev,maxncv,ldv,nloc
  integer                           :: n,nconv,ncv,nev
  !Arrays:
  complex(8),allocatable            :: evec_tmp(:)
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
  logical                           :: verb,reduce_
  integer                           :: i
  real(8)                           :: sigma
  real(8)                           :: tol_
  character                         :: bmat  
  character(len=2)                  :: which_
  real(8),external                  :: dznrm2,dlapy2
  real(8),allocatable               :: reV(:),imV(:)
  integer,allocatable               :: Eorder(:)
  !MPI
  integer                           :: mpi_ierr
  integer                           :: mpi_rank
  integer                           :: mpi_size
  logical                           :: mpi_master
  integer                           :: mpiQ
  integer                           :: mpiR
  !
  which_ = 'SR'   ; if(present(which))which_=which
  tol_   = 1d-12  ; if(present(tol))tol_=tol
  verb   = .false.; if(present(iverbose))verb=iverbose
  reduce_= .true.; if(present(mpi_reduce))reduce_=mpi_reduce
  !
  !MPI setup:
  mpi_size  = Get_size_MPI(MpiComm)
  mpi_rank  = Get_rank_MPI(MpiComm)
  mpi_master= Get_master_MPI(MpiComm)
  mpiQ = Ns/mpi_size
  mpiR = 0
  if(mpi_rank == mpi_size-1)mpiR=mod(Ns,mpi_size)
  !
  maxn   = Ns
  maxnev = Neigen
  maxncv = Nblock
  if(maxncv>Ns)maxncv=Ns
  !
  n      = maxn
  nev    = maxnev
  ncv    = maxncv
  bmat   = 'I'
  maxitr = Nitermax
  !
  ! Setup distribution of data to nodes:
  ldv = mpiQ+mpiR             !ldv is the SMALL dimension
  if(ldv < ncv)stop "LANCZOS_PARPACK_C error: ldv < maxNblock. Unstable! Find a way to increase ldv... (less cpu?)"
  !
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
  ishfts    = 1
  mode1     = 1
  iparam(1) = ishfts
  iparam(3) = maxitr
  iparam(7) = mode1
  !
  !
  if(present(v0))then
     do i=1 + mpi_rank*mpiQ, (mpi_rank+1)*mpiQ + mpiR
        resid(i-mpi_rank*mpiQ)=v0(i)
     enddo
  else
     call random_seed(size=nrandom)
     if(allocated(seed_random))deallocate(seed_random)
     allocate(seed_random(nrandom))
     seed_random=1234567
     call random_seed(put=seed_random)
     allocate(vec(n),reV(n),imV(n))
     call random_number(reV)
     call random_number(imV)
     vec=dcmplx(reV,imV)
     vec=vec/sqrt(dot_product(vec,vec))
     do i=1 + mpi_rank*mpiQ, (mpi_rank+1)*mpiQ + mpiR
        resid(i-mpi_rank*mpiQ)=vec(i)
     enddo
     deallocate(reV,imV)
  endif
  !
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
  if(info<0)then
     if(mpi_master)write(*,'(a,i6)')'Fatal Error in PZNAUPD, info = ', info
     stop
     !
  else
     !
     call pzneupd (MpiComm,.true.,'All',select,d,v,ldv,sigma,workev,bmat,&
          mpiQ+mpiR,which_,nev,tol_,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,rwork,ierr)
     !
     do j=1,neigen
        eval(j)=dreal(d(j))
     enddo
     allocate(Eorder(Neigen))
     call sort_array(Eval,Eorder)
     !
     ! if(present(evec))then
     select case(reduce_)
     case(.false.)
        if(any(shape(evec)/=[ldv,Neigen]))stop "arpack_mpi ERROR: evec has wrong dimension. Reduce=T"
        evec=zero
        do j=1,neigen
           do i=1,mpiQ+mpiR
              evec(i,j)=v(i,Eorder(j))!j)
           enddo
        enddo
        !
     case(.true.)
        if(any(shape(evec)/=[Ns,Neigen]))stop "arpack_mpi ERROR: evec has wrong dimension. Reduce=F"
        evec=zero
        allocate(evec_tmp(Ns))
        do j=1,neigen
           evec_tmp=zero
           do i=1 + mpi_rank*mpiQ , (mpi_rank+1)*mpiQ+mpiR
              evec_tmp(i)=v(i-mpi_rank*mpiQ,Eorder(j))!j)
           enddo
           call MPI_Allreduce(evec_tmp,evec(:,j),Ns,MPI_DOUBLE_COMPLEX,MPI_SUM,MpiComm,mpi_ierr)
        enddo
        deallocate(evec_tmp)
        !
     end select
     ! endif
     !

     !=========================================================================
     !  Compute the residual norm
     !    ||  A*x - lambda*x ||
     !  for the NCONV accurately computed eigenvalues and 
     !  eigenvectors.  (iparam(5) indicates how many are 
     !  accurate to the requested tolerance)
     if(ierr/=0)then
        write(*,'(a,i6)')'Error with PDSEUPD, IERR = ',ierr
        write(*,'(a)')'Check the documentation of SSEUPD.'
     else
        nconv =  iparam(5)
        !
        ! do j = 1, nconv
        !    call hprod( n, v(1,j), ax )
        !    call zaxpy( n, -d(j), v(1,j), 1, ax, 1 )
        !    rd(j,1) = dble (d(j))
        !    rd(j,2) = dimag (d(j))
        !    rd(j,3) = dznrm2 (n, ax, 1)
        !    rd(j,3) = rd(j,3) / dlapy2 (rd(j,1),rd(j,2))
        ! end do
        ! if(mpi_master.AND.verb)call dmout(6,nconv,3,rd,maxncv,-6,'Ritz values and relative residuals')
     end if
     !
     if(mpi_master.AND.verb)then
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



contains

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


end subroutine lanczos_parpack_c
