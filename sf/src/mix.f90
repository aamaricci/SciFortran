subroutine mix(fin,weight)
  complex(8),dimension(:),intent(inout) :: fin
  real(8),intent(in)                    :: weight
  complex(8),dimension(:),allocatable   :: old_fin
  if(.not.allocated(old_fin))then
     allocate(old_fin(size(fin)))
     old_fin=0.d0
     return
  endif
  fin=weight*fin + (1.d0-weight)*old_fin
end subroutine mix
