!-------------------------------------
elemental function modulo3d(V) result(res)
  type(vect3D),intent(in)  :: V
  real(8) :: res
  res = sqrt(V%x**2 + V%y**2 + V%z**2)
end function modulo3d
!-------------------------------------


!-------------------------------------
elemental function add3d(V,W) result(res)
  type(vect3D),intent(in)  :: V,W
  type(vect3D) :: res
  res%x = V%x + W%x
  res%y = V%y + W%y
  res%z = V%z + W%z
end function add3d
!-------------------------------------


!-------------------------------------
elemental function subtract3d(V,W) result(res)
  type(vect3D),intent(in)  :: V,W
  type(vect3D) :: res
  res%x = V%x - W%x
  res%y = V%y - W%y
  res%z = V%z - W%y
end function subtract3d
!-------------------------------------


!-------------------------------------
elemental function prodL3d(C,V) result(res)
  real(8),intent(in) :: C
  type(vect3D),intent(in)  :: V
  type(vect3D) :: res
  res%x = C*V%x 
  res%y = C*V%y
  res%z = C*V%z
end function prodL3d
!-------------------------------------


!-------------------------------------
elemental function prodR3d(V,C) result(res)
  real(8),intent(in) :: C
  type(vect3D),intent(in)  :: V
  type(vect3D) :: res
  res%x = V%x*C 
  res%y = V%y*C
  res%z = V%z*C
end function prodR3d
!-------------------------------------


!-------------------------------------
elemental function dot3d(V,W) result(res)
  type(vect3D),intent(in)  :: V,W
  real(8) :: res
  res = V%x*W%x + V%y*W%y
end function dot3d
!-------------------------------------

!-------------------------------------
elemental subroutine scalar_eq3d(V,C)
  real(8),intent(in)  :: C
  type(vect3D),intent(inout) :: V
  V%x = C
  V%y = C
  V%z = C
end subroutine  scalar_eq3d
!-------------------------------------

!-------------------------------------
elemental subroutine vector_eq3d(V,W)
  type(vect3D),intent(in)  :: W
  type(vect3D),intent(inout) :: V
  V%x = W%x
  V%y = W%y
  V%z = W%z
end subroutine  vector_eq3d
