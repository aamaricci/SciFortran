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
  integer,allocatable             :: Eorder(:)
  !
  which_='SR';if(present(which))which_=which
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
     end do
     allocate(Eorder(Neigen))
     call sort_array(Eval,Eorder)
     do j=1,neigen
        evec(:,j)=v(:,Eorder(j))
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

end subroutine lanczos_arpack_c
