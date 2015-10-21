!-------------------------------------
elemental function modulo2d(V) result(res)
  type(vect2D),intent(in)  :: V
  real(8) :: res
  res = sqrt(V%x**2 + V%y**2)
end function modulo2d
!-------------------------------------


!-------------------------------------
elemental function add2d(V,W) result(res)
  type(vect2D),intent(in)  :: V,W
  type(vect2D) :: res
  res%x = V%x + W%x
  res%y = V%y + W%y
end function add2d
!-------------------------------------


!-------------------------------------
elemental function subtract2d(V,W) result(res)
  type(vect2D),intent(in)  :: V,W
  type(vect2D) :: res
  res%x = V%x - W%x
  res%y = V%y - W%y
end function subtract2d
!-------------------------------------

!-------------------------------------
elemental function prodL2d(C,V) result(res)
  real(8),intent(in) :: C
  type(vect2D),intent(in)  :: V
  type(vect2D) :: res
  res%x = C*V%x 
  res%y = C*V%y
end function prodL2d
!-------------------------------------

!-------------------------------------
elemental function prodR2d(V,C) result(res)
  real(8),intent(in) :: C
  type(vect2D),intent(in)  :: V
  type(vect2D) :: res
  res%x = V%x*C 
  res%y = V%y*C
end function prodR2d
!-------------------------------------

!-------------------------------------
elemental function dot2d(V,W) result(res)
  type(vect2D),intent(in)  :: V,W
  real(8) :: res
  res = V%x*W%x + V%y*W%y
end function dot2d
!-------------------------------------

!-------------------------------------
elemental subroutine scalar_eq2d(V,C)
  real(8),intent(in)  :: C
  type(vect2D),intent(inout) :: V
  V%x = C
  V%y = C
end subroutine  scalar_eq2d
!-------------------------------------

!-------------------------------------
elemental subroutine vector_eq2d(V,W)
  type(vect2D),intent(in)  :: W
  type(vect2D),intent(inout) :: V
  V%x = W%x
  V%y = W%y
end subroutine  vector_eq2d
