MODULE PADE
  USE COMMON_VARS
  USE TOOLS
  implicit none
  private

  integer                 :: Npade
  complex(16),allocatable :: p(:,:),z(:)
  public :: pade_analytic_continuation

contains


  function pade_analytic_continuation(x,n,gm,wm) result(gr)
    complex(8),dimension(:),intent(in)        :: x
    integer,intent(in)                        :: n
    complex(8),dimension(:),intent(in)        :: gm
    real(8),dimension(size(gm)),intent(in)    :: wm
    complex(8),dimension(size(x))             :: gr
    integer                                   :: i
    Npade=n ; if(size(gm) < Npade)call error("error in pade_analytic_continuation. STOP")
    call set_pade_coefficient(gm,wm)
    do i=1,size(x)
       gr(i) = gpade(x(i))
    enddo
    call end_pade_coefficient
  end function pade_analytic_continuation


  subroutine set_pade_coefficient(gm,wm)
    integer                    :: i,j
    complex(8),dimension(:)    :: gm
    real(8),dimension(size(gm)):: wm
    allocate(p(Npade,Npade)) ; p=zero
    allocate(z(Npade))       ; z=xi*wm(1:Npade)
    p(1,1:Npade) = gm(1:Npade)
    do j=2,Npade
       do i=2,j
          p(i,j)=(p(i-1,i-1)-p(i-1,j))/(z(j)-z(i-1))/p(i-1,j)
       enddo
    enddo
  end subroutine set_pade_coefficient


  function gpade(x) result(ge)
    integer    :: i
    complex(8) :: x,ge
    complex(8) :: aw(0:Npade),b(0:Npade)
    if(size(z) < Npade)call error("error in gpade")
    aw(0)=zero ; aw(1)=p(1,1)
    b(0)=one   ; b(1)=one
    do i=1,Npade-1
       aw(i+1)=aw(i)+(x-z(i))*p(i+1,i+1)*aw(i-1)
       b(i+1)=b(i)+(x-z(i))*p(i+1,i+1)*b(i-1)
    enddo
    ge=aw(Npade)/b(Npade)
  end function gpade


  subroutine end_pade_coefficient
    deallocate(p,z)
  end subroutine end_pade_coefficient


END MODULE PADE
