!+-------------------------------------------------------------------+
!PROGRAM  : 
!TYPE     : Subroutine
!PURPOSE  :  
!+-------------------------------------------------------------------+
subroutine Interp_Gtau(FG1,FG2,L1,L2)
  integer             :: L1, L2
  real(8)             :: FG1(0:L1), FG2(0:L2)
  real(8)             :: FctL1(-L1:L1), FctL2(-L2:L2)
  real(8),allocatable :: xa(:), ya(:,:)!, y2(:)
  integer             :: L11, L12, L13
  real(8)             :: x, y
  integer             :: i, j
  FctL1(0:L1)=FG1 ; forall(i=1:L1)FctL1(-i)=-FctL1(L1-i)
  L11 = L1 + 1
  L12 = L1 + 2
  L13 = L1 + 3
  allocate(xa(L11),ya(4,L11))
  do i=1, L11
     xa(i)=dble(i-1)/dble(L1)
     ya(1,i)=FctL1(i-1)
  enddo
  call CUBSPL(xa,ya,L11,0,0)
  do i=1, L2
     x=dble(i)/dble(L2)
     FctL2(i)=PPVALU(xa,ya,L1,4,x,0)
  enddo
  FctL2(0)=FctL1(0)
  do i=1, L2
     FctL2(-i)=-FctL2(L2-i)
  enddo
  FG2=FctL2(0:L2)
end subroutine Interp_Gtau
!*******************************************************************
!*******************************************************************
!*******************************************************************


!+-------------------------------------------------------------------+
!PROGRAM  : EXTRACT
!TYPE     : Subroutine
!PURPOSE  : Sample a given function G(tau) over Nfak < N points.
!COMMENTS : Incoming function is expected to have N+1 points (tau=0,beta)
!this is modification with respect to the std extract routine used in
!HFqmc.
! g0 has N+1 points
! g00 will have Nfak+1 (tau_fak=0,beta)
!+-------------------------------------------------------------------+
subroutine extract_gtau(g0,g00)
  real(8),dimension(0:)  :: g0 !0:L
  real(8),dimension(0:) :: g00 !0:Lfak
  integer :: N,Nfak
  integer :: i,ip
  real(8) :: p,mismatch
  N=size(g0)-1
  Nfak=size(g00)-1
  g00(0)=g0(0)
  if(g0(0) > 0.d0) g00(Nfak)=1.d0-g0(0)
  if(g0(0) < 0.d0) g00(Nfak)=-(g0(0)+1.d0)
  mismatch=dble(N)/dble(Nfak)
  do i=1,Nfak-1
     p=dble(i)*mismatch
     ip=int(p)
     g00(i)=g0(ip)
  enddo
end subroutine extract_gtau
!*******************************************************************
!*******************************************************************
!*******************************************************************
