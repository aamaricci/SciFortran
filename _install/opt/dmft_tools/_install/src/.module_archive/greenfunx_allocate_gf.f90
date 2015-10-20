subroutine allocate_matsubara_gf(fgk,L)
  type(matsubara_gf),intent(inout) :: fgk
  integer,intent(in) :: L
  allocate(fgk%iw(L),fgk%tau(0:L))
  fgk%status=.true.
end subroutine allocate_matsubara_gf
subroutine allocate_matsubara_gf_(fgk,L)
  type(matsubara_gf),intent(inout) :: fgk(:)
  integer,intent(in) :: L
  integer            :: L1,i
  L1=size(fgk,1)
  do i=1,L1
     allocate(fgk(i)%iw(L),fgk(i)%tau(0:L))
     fgk(i)%status=.true.
  enddo
end subroutine allocate_matsubara_gf_
subroutine allocate_matsubara_gf__(fgk,L)
  type(matsubara_gf),intent(inout) :: fgk(:,:)
  integer,intent(in) :: L
  integer            :: L1,L2,i,j
  L1=size(fgk,1);L2=size(fgk,2)
  do i=1,L1
     do j=1,L2
        allocate(fgk(i,j)%iw(L),fgk(i,j)%tau(0:L))
        fgk(i,j)%status=.true.
     enddo
  enddo
end subroutine allocate_matsubara_gf__



subroutine allocate_real_gf(fgk,L)
  type(real_gf),intent(inout) :: fgk
  integer,intent(in)          :: L
  allocate(fgk%w(2*L),fgk%t(-L:L))
  fgk%status=.true.
end subroutine allocate_real_gf
subroutine allocate_real_gf_(fgk,L)
  type(real_gf),intent(inout) :: fgk(:)
  integer,intent(in)          :: L
  integer                     :: L1,i
  L1=size(fgk,1)
  do i=1,L1
     allocate(fgk(i)%w(2*L),fgk(i)%t(-L:L))
     fgk(i)%status=.true.
  enddo
end subroutine allocate_real_gf_
subroutine allocate_real_gf__(fgk,L)
  type(real_gf),intent(inout) :: fgk(:,:)
  integer,intent(in)          :: L
  integer                     :: L1,L2,i,j
  L1=size(fgk,1);L2=size(fgk,2)
  do i=1,L1
     do j=1,L2
        allocate(fgk(i,j)%w(2*L),fgk(i,j)%t(-L:L))
        fgk(i,j)%status=.true.
     enddo
  enddo
end subroutine allocate_real_gf__



subroutine allocate_keldysh_equilibrium_component(fgk,L)
  type(keldysh_equilibrium_component),intent(inout) :: fgk
  integer,intent(in)                                :: L
  allocate(fgk%w(2*L),fgk%t(-L:L))
end subroutine allocate_keldysh_equilibrium_component

subroutine allocate_keldysh_equilibrium_gf(fgk,L)
  type(keldysh_equilibrium_gf),intent(inout)        :: fgk
  integer,intent(in)                                :: L
  call allocate_keldysh_equilibrium_component(fgk%less,L)
  call allocate_keldysh_equilibrium_component(fgk%gtr,L)
  call allocate_keldysh_equilibrium_component(fgk%ret,L)
  fgk%status=.true.
end subroutine allocate_keldysh_equilibrium_gf
