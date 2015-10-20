subroutine matsubara_gf_identity(fgk1,fgk2)
  type(matsubara_gf),intent(inout) :: fgk1
  type(matsubara_gf),intent(in)    :: fgk2
  fgk1%iw  = fgk2%iw
  fgk1%tau = fgk2%tau
end subroutine matsubara_gf_identity
subroutine matsubara_gf_identity_(fgk1,fgk2)
  type(matsubara_gf),intent(inout) :: fgk1(:)
  type(matsubara_gf),intent(in)    :: fgk2(:)
  integer            :: L1,i
  L1=size(fgk1,1)
  do i=1,L1
     fgk1(i)%iw  = fgk2(i)%iw
     fgk1(i)%tau = fgk2(i)%tau
  enddo
end subroutine matsubara_gf_identity_
subroutine matsubara_gf_identity__(fgk1,fgk2)
  type(matsubara_gf),intent(inout) :: fgk1(:,:)
  type(matsubara_gf),intent(in)    :: fgk2(:,:)
  integer                          :: L1,L2,i,j
  L1=size(fgk1,1);L2=size(fgk1,2)
  do i=1,L1
     do j=1,L2
        fgk1(i,j)%iw  = fgk2(i,j)%iw
        fgk1(i,j)%tau = fgk2(i,j)%tau
     end do
  end do
end subroutine matsubara_gf_identity__


subroutine real_gf_identity(fgk1,fgk2)
  type(real_gf),intent(inout) :: fgk1
  type(real_gf),intent(in)    :: fgk2
  fgk1%w = fgk2%w
  fgk1%t  = fgk2%t
end subroutine real_gf_identity
subroutine real_gf_identity_(fgk1,fgk2)
  type(real_gf),intent(inout) :: fgk1(:)
  type(real_gf),intent(in)    :: fgk2(:)
  integer            :: L1,i
  L1=size(fgk1,1)
  do i=1,L1
     fgk1(i)%w  = fgk2(i)%w
     fgk1(i)%t  = fgk2(i)%t
  end do
end subroutine real_gf_identity_
subroutine real_gf_identity__(fgk1,fgk2)
  type(real_gf),intent(inout) :: fgk1(:,:)
  type(real_gf),intent(in)    :: fgk2(:,:)
  integer            :: L1,L2,i,j
  L1=size(fgk1,1);L2=size(fgk1,2)
  do i=1,L1
     do j=1,L2
        fgk1(i,j)%w  = fgk2(i,j)%w
        fgk1(i,j)%t  = fgk2(i,j)%t
     end do
  end do
end subroutine real_gf_identity__



subroutine keldysh_eq_component_gf_identity(fgk1,fgk2)
  type(keldysh_equilibrium_component),intent(inout) :: fgk1
  type(keldysh_equilibrium_component),intent(in)    :: fgk2
  fgk1%t = fgk2%t
  fgk1%w = fgk2%w
end subroutine keldysh_eq_component_gf_identity

subroutine keldysh_eq_gf_identity(fgk1,fgk2)
  type(keldysh_equilibrium_gf),intent(inout) :: fgk1
  type(keldysh_equilibrium_gf),intent(in)    :: fgk2
  fgk1%less%t = fgk2%less%t
  fgk1%gtr%t  = fgk2%gtr%t
  fgk1%ret%t  = fgk2%ret%t
  fgk1%less%w = fgk2%less%w
  fgk1%gtr%w  = fgk2%gtr%w
  fgk1%ret%w  = fgk2%ret%w
end subroutine keldysh_eq_gf_identity
