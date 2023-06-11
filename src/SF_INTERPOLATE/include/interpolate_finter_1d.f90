subroutine init_finter_d(func,Xin,Fin,N)
  type(finter_type) :: func
  real(8)           :: xin(:)
  real(8)           :: fin(size(xin))
  integer           :: N,Lin
  if(func%status)call delete_finter(func)
  Lin=size(xin)
  allocate(func%x(Lin),func%f(Lin))
  func%X    = Xin
  func%F    = Fin
  func%Imax = Lin
  func%Imin = 1
  func%N    = N
  func%status=.true.
end subroutine init_finter_d
!
subroutine init_finter_c(func,Xin,Fin,N)
  type(finter_type) :: func
  real(8)           :: xin(:)
  complex(8)        :: fin(size(xin))
  integer           :: N,Lin
  if(func%status)call delete_finter(func)
  Lin=size(xin)
  allocate(func%x(Lin),func%f(Lin),func%g(Lin))
  func%X    = Xin
  func%F    = dreal(fin)
  func%G    = dimag(fin)
  func%Imax = Lin
  func%Imin = 1
  func%N    = N
  func%status=.true.
end subroutine init_finter_c




!*******************************************************************
!*******************************************************************
!*******************************************************************




subroutine delete_finter(func)
  type(finter_type) :: func
  if(allocated(func%x))deallocate(func%x)
  if(allocated(func%f))deallocate(func%f)
  if(allocated(func%g))deallocate(func%g)
  func%imax=0
  func%imin=0
  func%N   =0
  func%status=.false.
end subroutine delete_finter





!*******************************************************************
!*******************************************************************
!*******************************************************************



function finter(func,x)
  real(8)           :: x
  type(finter_type) :: func
  real(8)           :: finter
  real(8)           :: y,dy
  integer           :: j,k,k0,k1
  integer           :: n
  N=func%N    !order of polynomial interpolation
  finter=0.d0
  j=locate(func%X(func%Imin:func%Imax),x)
  !k = min(max(j-(N-1)/2,1),func%Imax+1-N)
  k=max(j-(N-1)/2,1)
  k0=k
  if(k0 < func%Imin)k0=func%Imin
  k1=k+N+1
  if(k1 > func%Imax)then
     k1=func%Imax
     k0=k1-N-1
  endif
  call polint(func%X(k0:k1),func%F(k0:k1),x,y,dy)
  !call polint(func%X(k:k+n),func%F(k:k+n),x,y,dy)
  finter=y
  return
end function finter

function cinter(func,x)
  real(8)           :: x
  type(finter_type) :: func
  complex(8)        :: cinter
  real(8)           :: ry,dry
  real(8)           :: iy,diy
  integer           :: j,k,k0,k1
  integer           :: n
  N=func%N    !order of polynomial interpolation
  cinter=dcmplx(0d0,0d0)
  !
  j=locate(func%X(func%Imin:func%Imax),x)
  k=max(j-(N-1)/2,1)
  k0=k
  if(k0 < func%Imin)k0=func%Imin
  k1=k+N+1
  if(k1 > func%Imax)then
     k1=func%Imax
     k0=k1-N-1
  endif
  !
  call polint(func%X(k0:k1),func%F(k0:k1),x,ry,dry)
  call polint(func%X(k0:k1),func%G(k0:k1),x,iy,diy)
  cinter=dcmplx(ry,iy)
  return
end function cinter


