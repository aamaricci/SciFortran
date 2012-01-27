!###############################################################
! PROGRAM  : FMINLN
! TYPE     : Module
! PURPOSE  : 
!###############################################################
MODULE BROYDN_FMINLN
  implicit none
  REAL(8), DIMENSION(:), POINTER :: fmin_fvecp
CONTAINS
  FUNCTION fmin(x)
    REAL(8), DIMENSION(:), INTENT(IN) :: x
    REAL(8)                           :: fmin
    INTERFACE
       FUNCTION funcv(x)
         REAL(8), DIMENSION(:), INTENT(IN) :: x
         REAL(8), DIMENSION(size(x))       :: funcv
       END FUNCTION funcv
    END INTERFACE
    if (.not. associated(fmin_fvecp)) then
       print*,'fmin: problem with pointer for returned values'
       stop
    endif
    fmin_fvecp=funcv(x)
    fmin=0.5d0*dot_product(fmin_fvecp,fmin_fvecp)
  END FUNCTION fmin
END MODULE BROYDN_FMINLN
