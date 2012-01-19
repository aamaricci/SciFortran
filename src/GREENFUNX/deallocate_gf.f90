subroutine deallocate_matsubara_gf(fgk)
  type(matsubara_gf),intent(inout) :: fgk
  deallocate(fgk%iw,fgk%tau)
  fgk%status=.false.
end subroutine deallocate_matsubara_gf
subroutine deallocate_matsubara_gf_(fgk)
  type(matsubara_gf),intent(inout) :: fgk(:)
  integer                          :: i,L1
  L1=size(fgk,1)
  do i=1,L1
     deallocate(fgk(i)%iw,fgk(i)%tau)
     fgk(i)%status=.false.
  enddo
end subroutine deallocate_matsubara_gf_
subroutine deallocate_matsubara_gf__(fgk)
  type(matsubara_gf),intent(inout) :: fgk(:,:)
  integer                          :: L1,L2,i,j
  L1=size(fgk,1);L2=size(fgk,2)
  do i=1,L1
     do j=1,L2
        deallocate(fgk(i,j)%iw,fgk(i,j)%tau)
        fgk(i,j)%status=.false.
     enddo
  enddo
end subroutine deallocate_matsubara_gf__


subroutine deallocate_real_gf(fgk)
  type(real_gf),intent(inout)      :: fgk
  deallocate(fgk%w,fgk%t)
  fgk%status=.false.
end subroutine deallocate_real_gf
subroutine deallocate_real_gf_(fgk)
  type(real_gf),intent(inout)      :: fgk(:)
  integer                          :: L1,i
  L1=size(fgk,1)
  do i=1,L1
     deallocate(fgk(i)%w,fgk(i)%t)
     fgk(i)%status=.false.
  enddo
end subroutine deallocate_real_gf_
subroutine deallocate_real_gf__(fgk)
  type(real_gf),intent(inout)      :: fgk(:,:)
  integer                          :: L1,L2,i,j
  L1=size(fgk,1);L2=size(fgk,2)
  do i=1,L1
     do j=1,L2
        deallocate(fgk(i,j)%w,fgk(i,j)%t)
        fgk(i,j)%status=.false.
     enddo
  enddo
end subroutine deallocate_real_gf__



subroutine deallocate_keldysh_equilibrium_component(fgk)
  type(keldysh_equilibrium_component),intent(inout) :: fgk
  deallocate(fgk%w,fgk%t)
end subroutine deallocate_keldysh_equilibrium_component

subroutine deallocate_keldysh_equilibrium_gf(fgk)
  type(keldysh_equilibrium_gf),intent(inout)        :: fgk
  call deallocate_keldysh_equilibrium_component(fgk%less)
  call deallocate_keldysh_equilibrium_component(fgk%gtr)
  call deallocate_keldysh_equilibrium_component(fgk%ret)
  fgk%status=.false.
end subroutine deallocate_keldysh_equilibrium_gf
