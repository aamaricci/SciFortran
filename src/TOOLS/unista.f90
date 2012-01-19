!   unista = (stable unique) removes duplicates from an array,
!            leaving unique entries in the order of their first
!            appearance in the initial set.
MODULE M_UNISTA
  USE M_UNIINV
  implicit none
  private

  interface unista
     module procedure d_unista, r_unista, i_unista
  end interface unista

  public  :: unista

contains

  subroutine d_unista (xdont, nuni, mask)
    real(kind=8), dimension (:), intent (inout) :: xdont
    integer, intent (out)                       :: nuni
    integer, dimension (size(xdont)) :: iwrkt
    logical, dimension (size(xdont)) :: ifmptyt
    logical, dimension (size(xdont)),optional :: mask
    integer :: icrs
    call uniinv (xdont, iwrkt)

    ifmptyt = .true.
    nuni = 0
    do icrs = 1, size(xdont)
       if(present(mask))mask(icrs)=ifmptyt(iwrkt(icrs))
       if (ifmptyt(iwrkt(icrs))) then
          ifmptyt(iwrkt(icrs)) = .false. !modify the mask to that next acces is false
          nuni = nuni + 1
          xdont (nuni) = xdont (icrs)
       end if
    end do
    return
  end subroutine d_unista

  subroutine r_unista (xdont, nuni, mask)
    real, dimension (:), intent (inout) :: xdont
    integer, intent (out)               :: nuni
    integer, dimension (size(xdont)) :: iwrkt
    logical, dimension (size(xdont)) :: ifmptyt
    logical, dimension (size(xdont)),optional :: mask
    integer :: icrs
    call uniinv (xdont, iwrkt)
    ifmptyt = .true.
    nuni = 0
    do icrs = 1, size(xdont)
       if(present(mask))mask(icrs)=ifmptyt(iwrkt(icrs))
       if (ifmptyt(iwrkt(icrs))) then
          ifmptyt(iwrkt(icrs)) = .false.
          nuni = nuni + 1
          xdont (nuni) = xdont (icrs)
       end if
    end do
    return
  end subroutine r_unista

  subroutine i_unista (xdont, nuni, mask)
    integer, dimension (:), intent (inout)  :: xdont
    integer, intent (out) :: nuni
    integer, dimension (size(xdont)) :: iwrkt
    logical, dimension (size(xdont)) :: ifmptyt
    logical, dimension (size(xdont)),optional :: mask
    integer :: icrs
    call uniinv (xdont, iwrkt)
    ifmptyt = .true.
    nuni = 0
    do icrs = 1, size(xdont)
       if(present(mask))mask(icrs)=ifmptyt(iwrkt(icrs))
       if (ifmptyt(iwrkt(icrs))) then
          ifmptyt(iwrkt(icrs)) = .false.
          nuni = nuni + 1
          xdont (nuni) = xdont (icrs)
       end if
    end do
    return
  end subroutine i_unista

END MODULE M_UNISTA
