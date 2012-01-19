MODULE CGFIT_F1DIM_MOD
  implicit none
  INTEGER                        :: ncom
  REAL(8), DIMENSION(:), POINTER :: pcom,xicom

CONTAINS


  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  FUNCTION f1dim(x)
    REAL(8), INTENT(IN)                :: x
    real(8)                            :: f1dim
    REAL(8), DIMENSION(:), ALLOCATABLE :: xt
    INTERFACE
       FUNCTION func(x)
         REAL(8), DIMENSION(:),INTENT(IN) :: x
         REAL(8)                          :: func
       END FUNCTION func
    END INTERFACE
    allocate(xt(ncom))
    xt(:)=pcom(:)+x*xicom(:)
    f1dim=func(xt)
    deallocate(xt)
  END FUNCTION f1dim

  !********************************************************************
  !********************************************************************
  !********************************************************************

END MODULE CGFIT_F1DIM_MOD
