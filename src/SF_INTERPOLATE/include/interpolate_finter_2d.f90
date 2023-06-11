subroutine init_finter2d(func,Xin,Yin,Fin,N)
  type(finter2d_type) :: func
  real(8)       :: xin(:),yin(:)
  real(8)       :: fin(size(xin),size(yin))
  integer       :: N,Lx,Ly
  if(func%status)deallocate(func%x,func%y,func%f)
  Lx=size(xin) ; Ly=size(yin)
  allocate(func%x(Lx),func%y(Ly),func%f(Lx,Ly))
  func%X    = Xin
  func%Y    = Yin
  func%F    = Fin
  func%Imin = 1
  func%Jmin = 1
  func%Imax = Lx
  func%Jmax = Ly
  func%N    = N
  func%status=.true.
end subroutine init_finter2d





subroutine delete_finter2d(func)
  type(finter2d_type) :: func
  if(allocated(func%x))deallocate(func%x)
  if(allocated(func%y))deallocate(func%y)
  if(allocated(func%f))deallocate(func%f)
  func%imin=0
  func%jmin=0
  func%imax=0
  func%jmax=0
  func%N   =0
  func%status=.false.
end subroutine delete_finter2d








function finter2d(func,x,y)
  real(8)         :: x,y
  type(finter2d_type) :: func
  real(8)         :: finter2d
  real(8)         :: f,df
  integer         :: itmp,jtmp,kx,ky,k0x,k0y,k1x,k1y
  integer         :: n
  N=func%N    !order of polynomial interpolation
  finter2d=0.d0
  itmp=locate(func%X(func%Imin:func%Imax),x)
  jtmp=locate(func%Y(func%Jmin:func%Jmax),y)
  kx=max(itmp-(N-1)/2,1)
  ky=max(jtmp-(N-1)/2,1)
  k0x = kx ; if(k0x < func%Imin)k0x=func%Imin
  k0y = ky ; if(k0y < func%Jmin)k0y=func%Jmin
  k1x = kx+N+1
  if(k1x > func%Imax)then         
     k1x=func%Imax
     k0x=k1x-N-1
  endif
  k1y = ky+N+1
  if(k1y > func%Jmax)then
     k1y=func%Jmax
     k0y=k1y-N-1
  endif
  call polin2(func%X(k0x:k1x),func%Y(k0y:k1y),func%F(k0x:k1x,k0y:k1y),x,y,f,df)
  finter2d=f
  return
end function finter2d

