program test_SF_ARRAY
  USE SF_ARRAYS
  USE ASSERTING
  implicit none

  integer,parameter :: L=10,P=3,Q=10,N=P*Q+1
  integer           :: irange(L)
  real(8)           :: grid(L),upmgrid(N),upmint(2*P*Q+1),fmesh

  grid = linspace(0d0,1d0,L)
  call assert([grid(1),grid(L)],[0d0,1d0],"LINSPACE")
  write(*,"(10F12.8)")grid

  grid = linspace(0d0,1d0,L,iend=.false.)
  call assert([grid(1),grid(L)],[0d0,0.9d0],"LINSPACE IEND=F")
  write(*,"(10F12.8)")grid

  grid = linspace(0d0,1d0,L,istart=.false.)
  call assert([grid(1),grid(L)],[0.1d0,1d0],"LINSPACE ISTART=F")
  write(*,"(10F12.8)")grid

  grid = linspace(0d0,1d0,L,istart=.false.,mesh=fmesh)
  call assert(fmesh,0.1d0,"LINSPACE FMESH")
  write(*,"(10F12.8)")grid


  grid = logspace(1d-1,1d0,L)
  call assert([grid(1),grid(L)],[0.1d0,1d0],"LOGSPACE")
  write(*,"(10F12.8)")grid

  grid = logspace(1d-1,1d0,L,base=2d0)
  call assert([grid(1),grid(L)],[0.1d0,1d0],"LOGSPACE BASE=2")
  write(*,"(10F12.8)")grid


  grid = powspace(0d0,1d0,L,base=2d0)
  call assert([grid(1),grid(L)],[0d0,1d0],"POWSPACE BASE=2")
  write(*,"(10F12.8)")grid


  irange = arange(1,L)
  call assert([irange(1),irange(5),irange(L)],[1,5,L],"ARANGE")
  write(*,"(10I12)")irange


  upmgrid = upmspace(0d0,10d0,P,Q,N,base=10d0)
  call assert([upmgrid(10),upmgrid(20),upmgrid(30)],[0.09d0,0.91d0,9.1d0],"UPMSPACE P=10, Q=3")
  write(*,"(10F12.8)")upmgrid

  upmint = upminterval(0d0,10d0,5d0,P,Q,type=1,base=10d0)
  call assert([upmint(30),upmint(31),upmint(32)],[4.995d0,5d0,5.005d0],"UPMINTERVAL P=10, Q=3")
  write(*,"(10F12.8)")upmint

end program test_SF_ARRAY
