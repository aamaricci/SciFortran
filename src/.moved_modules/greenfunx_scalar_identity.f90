subroutine matsubara_scalar_identity(fgk,C)
  type(matsubara_gf),intent(inout) :: fgk
  complex(8),intent(in)            :: C
  fgk%iw = C
  fgk%tau= real(C,8)
end subroutine matsubara_scalar_identity
subroutine matsubara_scalar_identity_(fgk,C)
  type(matsubara_gf),intent(inout) :: fgk(:)
  complex(8),intent(in)            :: C
  integer                          :: L1,i
  L1=size(fgk,1)
  do i=1,L1
     fgk(i)%iw = C
     fgk(i)%tau= real(C,8)
  enddo
end subroutine matsubara_scalar_identity_
subroutine matsubara_scalar_identity__(fgk,C)
  type(matsubara_gf),intent(inout) :: fgk(:,:)
  complex(8),intent(in)            :: C
  integer            :: L1,L2,i,j
  L1=size(fgk,1);L2=size(fgk,2)
  do i=1,L1
     do j=1,L2
        fgk(i,j)%iw = C
        fgk(i,j)%tau= real(C,8)
     end do
  end do
end subroutine matsubara_scalar_identity__


subroutine real_scalar_identity(fgk,C)
  type(real_gf),intent(inout) :: fgk
  complex(8),intent(in)        :: C
  fgk%w = C
  fgk%t  = C
end subroutine real_scalar_identity
subroutine real_scalar_identity_(fgk,C)
  type(real_gf),intent(inout) :: fgk(:)
  complex(8),intent(in)       :: C
  integer                     :: L1,i
  L1=size(fgk,1)
  do i=1,L1
     fgk(i)%w = C
     fgk(i)%t  = C
  end do
end subroutine real_scalar_identity_
subroutine real_scalar_identity__(fgk,C)
  type(real_gf),intent(inout) :: fgk(:,:)
  complex(8),intent(in)        :: C
  integer            :: L1,L2,i,j
  L1=size(fgk,1);L2=size(fgk,2)
  do i=1,L1
     do j=1,L2
        fgk(i,j)%w = C
        fgk(i,j)%t  = C
     end do
  end do
end subroutine real_scalar_identity__


subroutine keldysh_eq_component_scalar_identity(fgk,C)
  type(keldysh_equilibrium_component),intent(inout) :: fgk
  complex(8),intent(in)                      :: C
  fgk%t = C
  fgk%w = C
end subroutine keldysh_eq_component_scalar_identity

subroutine keldysh_eq_scalar_identity(fgk,C)
  type(keldysh_equilibrium_gf),intent(inout) :: fgk
  complex(8),intent(in)                      :: C
  fgk%less%t = C
  fgk%gtr%t  = C
  fgk%ret%t  = C
  fgk%less%w = C
  fgk%gtr%w  = C
  fgk%ret%w  = C
end subroutine keldysh_eq_scalar_identity
