!----------------------------------------------------------------------
!   UNIINV = Merge-sort inverse ranking of an array, with removal of
!   duplicate entries.
!   The routine is similar to pure merge-sort ranking, but on
!   the last pass, it sets indices in IGOEST to the rank
!   of the value in the ordered set with duplicates removed.
!   For performance reasons, the first 2 passes are taken
!   out of the standard loop, and use dedicated coding.
!----------------------------------------------------------------------
Module m_uniinv
  implicit none
  private 

  interface uniinv
     module procedure d_uniinv, r_uniinv, i_uniinv
  end interface uniinv

  interface nearless
     module procedure D_nearless, R_nearless, I_nearless
  end interface nearless

  public :: uniinv

contains

  subroutine d_uniinv (xdont, igoest)
    real (kind=8), dimension (:), intent (in) :: xdont
    integer, dimension (:), intent (out)      :: igoest
    real (kind=8) :: xtst, xdona, xdonb
    integer, dimension (size(igoest)) :: jwrkt, irngt
    integer :: lmtna, lmtnc, irng, irng1, irng2, nuni
    integer :: nval, iind, iwrkd, iwrk, iwrkf, jinda, iinda, iindb
    !
    nval = min (size(xdont), size(igoest))
    !
    select case (nval)
    case (:0)
       return
    case (1)
       igoest (1) = 1
       return
    case default
       continue
    end select
    !
    !  fill-in the index array, creating ordered couples
    !
    do iind = 2, nval, 2
       if (xdont(iind-1) < xdont(iind)) then
          irngt (iind-1) = iind - 1
          irngt (iind) = iind
       else
          irngt (iind-1) = iind
          irngt (iind) = iind - 1
       end if
    end do
    if (modulo (nval, 2) /= 0) then
       irngt (nval) = nval
    end if
    !
    !  we will now have ordered subsets a - b - a - b - ...
    !  and merge a and b couples into     c   -   c   - ...
    !
    lmtna = 2
    lmtnc = 4
    !
    !  first iteration. the length of the ordered subsets goes from 2 to 4
    !
    do
       if (nval <= 4) exit
       !
       !   loop on merges of a and b into c
       !
       do iwrkd = 0, nval - 1, 4
          if ((iwrkd+4) > nval) then
             if ((iwrkd+2) >= nval) exit
             !
             !   1 2 3
             !
             if (xdont(irngt(iwrkd+2)) <= xdont(irngt(iwrkd+3))) exit
             !
             !   1 3 2
             !
             if (xdont(irngt(iwrkd+1)) <= xdont(irngt(iwrkd+3))) then
                irng2 = irngt (iwrkd+2)
                irngt (iwrkd+2) = irngt (iwrkd+3)
                irngt (iwrkd+3) = irng2
                !
                !   3 1 2
                !
             else
                irng1 = irngt (iwrkd+1)
                irngt (iwrkd+1) = irngt (iwrkd+3)
                irngt (iwrkd+3) = irngt (iwrkd+2)
                irngt (iwrkd+2) = irng1
             end if
             exit
          end if
          !
          !   1 2 3 4
          !
          if (xdont(irngt(iwrkd+2)) <= xdont(irngt(iwrkd+3))) cycle
          !
          !   1 3 x x
          !
          if (xdont(irngt(iwrkd+1)) <= xdont(irngt(iwrkd+3))) then
             irng2 = irngt (iwrkd+2)
             irngt (iwrkd+2) = irngt (iwrkd+3)
             if (xdont(irng2) <= xdont(irngt(iwrkd+4))) then
                !   1 3 2 4
                irngt (iwrkd+3) = irng2
             else
                !   1 3 4 2
                irngt (iwrkd+3) = irngt (iwrkd+4)
                irngt (iwrkd+4) = irng2
             end if
             !
             !   3 x x x
             !
          else
             irng1 = irngt (iwrkd+1)
             irng2 = irngt (iwrkd+2)
             irngt (iwrkd+1) = irngt (iwrkd+3)
             if (xdont(irng1) <= xdont(irngt(iwrkd+4))) then
                irngt (iwrkd+2) = irng1
                if (xdont(irng2) <= xdont(irngt(iwrkd+4))) then
                   !   3 1 2 4
                   irngt (iwrkd+3) = irng2
                else
                   !   3 1 4 2
                   irngt (iwrkd+3) = irngt (iwrkd+4)
                   irngt (iwrkd+4) = irng2
                end if
             else
                !   3 4 1 2
                irngt (iwrkd+2) = irngt (iwrkd+4)
                irngt (iwrkd+3) = irng1
                irngt (iwrkd+4) = irng2
             end if
          end if
       end do
       !
       !  the cs become as and bs
       !
       lmtna = 4
       exit
    end do
    !
    !  iteration loop. each time, the length of the ordered subsets
    !  is doubled.
    !
    do
       if (2*lmtna >= nval) exit
       iwrkf = 0
       lmtnc = 2 * lmtnc
       !
       !   loop on merges of a and b into c
       !
       do
          iwrk = iwrkf
          iwrkd = iwrkf + 1
          jinda = iwrkf + lmtna
          iwrkf = iwrkf + lmtnc
          if (iwrkf >= nval) then
             if (jinda >= nval) exit
             iwrkf = nval
          end if
          iinda = 1
          iindb = jinda + 1
          !
          !  one steps in the c subset, that we create in the final rank array
          !
          !  make a copy of the rank array for the iteration
          !
          jwrkt (1:lmtna) = irngt (iwrkd:jinda)
          xdona = xdont (jwrkt(iinda))
          xdonb = xdont (irngt(iindb))
          !
          do
             iwrk = iwrk + 1
             !
             !  we still have unprocessed values in both a and b
             !
             if (xdona > xdonb) then
                irngt (iwrk) = irngt (iindb)
                iindb = iindb + 1
                if (iindb > iwrkf) then
                   !  only a still with unprocessed values
                   irngt (iwrk+1:iwrkf) = jwrkt (iinda:lmtna)
                   exit
                end if
                xdonb = xdont (irngt(iindb))
             else
                irngt (iwrk) = jwrkt (iinda)
                iinda = iinda + 1
                if (iinda > lmtna) exit! only b still with unprocessed values
                xdona = xdont (jwrkt(iinda))
             end if
             !
          end do
       end do
       !
       !  the cs become as and bs
       !
       lmtna = 2 * lmtna
    end do
    !
    !   last merge of a and b into c, with removal of duplicates.
    !
    iinda = 1
    iindb = lmtna + 1
    nuni = 0
    !
    !  one steps in the c subset, that we create in the final rank array
    !
    jwrkt (1:lmtna) = irngt (1:lmtna)
    if (iindb <= nval) then
       xtst = nearless (min(xdont(jwrkt(1)), xdont(irngt(iindb))))
    else
       xtst = nearless (xdont(jwrkt(1)))
    endif
    do iwrk = 1, nval
       !
       !  we still have unprocessed values in both a and b
       !
       if (iinda <= lmtna) then
          if (iindb <= nval) then
             if (xdont(jwrkt(iinda)) > xdont(irngt(iindb))) then
                irng = irngt (iindb)
                iindb = iindb + 1
             else
                irng = jwrkt (iinda)
                iinda = iinda + 1
             end if
          else
             !
             !  only a still with unprocessed values
             !
             irng = jwrkt (iinda)
             iinda = iinda + 1
          end if
       else
          !
          !  only b still with unprocessed values
          !
          irng = irngt (iwrk)
       end if
       if (xdont(irng) > xtst) then
          xtst = xdont (irng)
          nuni = nuni + 1
       end if
       igoest (irng) = nuni
       !
    end do
    !
    return
    !
  end subroutine d_uniinv

  subroutine r_uniinv (xdont, igoest)
    real, dimension (:), intent (in) :: xdont
    integer, dimension (:), intent (out) :: igoest
    real    :: xtst, xdona, xdonb
    integer, dimension (size(igoest)) :: jwrkt, irngt
    integer :: lmtna, lmtnc, irng, irng1, irng2, nuni
    integer :: nval, iind, iwrkd, iwrk, iwrkf, jinda, iinda, iindb
    !
    nval = min (size(xdont), size(igoest))
    !
    select case (nval)
    case (:0)
       return
    case (1)
       igoest (1) = 1
       return
    case default
       continue
    end select
    !
    !  fill-in the index array, creating ordered couples
    !
    do iind = 2, nval, 2
       if (xdont(iind-1) < xdont(iind)) then
          irngt (iind-1) = iind - 1
          irngt (iind) = iind
       else
          irngt (iind-1) = iind
          irngt (iind) = iind - 1
       end if
    end do
    if (modulo (nval, 2) /= 0) then
       irngt (nval) = nval
    end if
    !
    !  we will now have ordered subsets a - b - a - b - ...
    !  and merge a and b couples into     c   -   c   - ...
    !
    lmtna = 2
    lmtnc = 4
    !
    !  first iteration. the length of the ordered subsets goes from 2 to 4
    !
    do
       if (nval <= 4) exit
       !
       !   loop on merges of a and b into c
       !
       do iwrkd = 0, nval - 1, 4
          if ((iwrkd+4) > nval) then
             if ((iwrkd+2) >= nval) exit
             !
             !   1 2 3
             !
             if (xdont(irngt(iwrkd+2)) <= xdont(irngt(iwrkd+3))) exit
             !
             !   1 3 2
             !
             if (xdont(irngt(iwrkd+1)) <= xdont(irngt(iwrkd+3))) then
                irng2 = irngt (iwrkd+2)
                irngt (iwrkd+2) = irngt (iwrkd+3)
                irngt (iwrkd+3) = irng2
                !
                !   3 1 2
                !
             else
                irng1 = irngt (iwrkd+1)
                irngt (iwrkd+1) = irngt (iwrkd+3)
                irngt (iwrkd+3) = irngt (iwrkd+2)
                irngt (iwrkd+2) = irng1
             end if
             exit
          end if
          !
          !   1 2 3 4
          !
          if (xdont(irngt(iwrkd+2)) <= xdont(irngt(iwrkd+3))) cycle
          !
          !   1 3 x x
          !
          if (xdont(irngt(iwrkd+1)) <= xdont(irngt(iwrkd+3))) then
             irng2 = irngt (iwrkd+2)
             irngt (iwrkd+2) = irngt (iwrkd+3)
             if (xdont(irng2) <= xdont(irngt(iwrkd+4))) then
                !   1 3 2 4
                irngt (iwrkd+3) = irng2
             else
                !   1 3 4 2
                irngt (iwrkd+3) = irngt (iwrkd+4)
                irngt (iwrkd+4) = irng2
             end if
             !
             !   3 x x x
             !
          else
             irng1 = irngt (iwrkd+1)
             irng2 = irngt (iwrkd+2)
             irngt (iwrkd+1) = irngt (iwrkd+3)
             if (xdont(irng1) <= xdont(irngt(iwrkd+4))) then
                irngt (iwrkd+2) = irng1
                if (xdont(irng2) <= xdont(irngt(iwrkd+4))) then
                   !   3 1 2 4
                   irngt (iwrkd+3) = irng2
                else
                   !   3 1 4 2
                   irngt (iwrkd+3) = irngt (iwrkd+4)
                   irngt (iwrkd+4) = irng2
                end if
             else
                !   3 4 1 2
                irngt (iwrkd+2) = irngt (iwrkd+4)
                irngt (iwrkd+3) = irng1
                irngt (iwrkd+4) = irng2
             end if
          end if
       end do
       !
       !  the cs become as and bs
       !
       lmtna = 4
       exit
    end do
    !
    !  iteration loop. each time, the length of the ordered subsets
    !  is doubled.
    !
    do
       if (2*lmtna >= nval) exit
       iwrkf = 0
       lmtnc = 2 * lmtnc
       !
       !   loop on merges of a and b into c
       !
       do
          iwrk = iwrkf
          iwrkd = iwrkf + 1
          jinda = iwrkf + lmtna
          iwrkf = iwrkf + lmtnc
          if (iwrkf >= nval) then
             if (jinda >= nval) exit
             iwrkf = nval
          end if
          iinda = 1
          iindb = jinda + 1
          !
          !  one steps in the c subset, that we create in the final rank array
          !
          !  make a copy of the rank array for the iteration
          !
          jwrkt (1:lmtna) = irngt (iwrkd:jinda)
          xdona = xdont (jwrkt(iinda))
          xdonb = xdont (irngt(iindb))
          !
          do
             iwrk = iwrk + 1
             !
             !  we still have unprocessed values in both a and b
             !
             if (xdona > xdonb) then
                irngt (iwrk) = irngt (iindb)
                iindb = iindb + 1
                if (iindb > iwrkf) then
                   !  only a still with unprocessed values
                   irngt (iwrk+1:iwrkf) = jwrkt (iinda:lmtna)
                   exit
                end if
                xdonb = xdont (irngt(iindb))
             else
                irngt (iwrk) = jwrkt (iinda)
                iinda = iinda + 1
                if (iinda > lmtna) exit! only b still with unprocessed values
                xdona = xdont (jwrkt(iinda))
             end if
             !
          end do
       end do
       !
       !  the cs become as and bs
       !
       lmtna = 2 * lmtna
    end do
    !
    !   last merge of a and b into c, with removal of duplicates.
    !
    iinda = 1
    iindb = lmtna + 1
    nuni = 0
    !
    !  one steps in the c subset, that we create in the final rank array
    !
    jwrkt (1:lmtna) = irngt (1:lmtna)
    if (iindb <= nval) then
       xtst = nearless (min(xdont(jwrkt(1)), xdont(irngt(iindb))))
    else
       xtst = nearless (xdont(jwrkt(1)))
    endif
    do iwrk = 1, nval
       !
       !  we still have unprocessed values in both a and b
       !
       if (iinda <= lmtna) then
          if (iindb <= nval) then
             if (xdont(jwrkt(iinda)) > xdont(irngt(iindb))) then
                irng = irngt (iindb)
                iindb = iindb + 1
             else
                irng = jwrkt (iinda)
                iinda = iinda + 1
             end if
          else
             !
             !  only a still with unprocessed values
             !
             irng = jwrkt (iinda)
             iinda = iinda + 1
          end if
       else
          !
          !  only b still with unprocessed values
          !
          irng = irngt (iwrk)
       end if
       if (xdont(irng) > xtst) then
          xtst = xdont (irng)
          nuni = nuni + 1
       end if
       igoest (irng) = nuni
       !
    end do
    !
    return
    !
  end subroutine r_uniinv
  subroutine i_uniinv (xdont, igoest)
    ! __________________________________________________________
    !   uniinv = merge-sort inverse ranking of an array, with removal of
    !   duplicate entries.
    !   the routine is similar to pure merge-sort ranking, but on
    !   the last pass, it sets indices in igoest to the rank
    !   of the value in the ordered set with duplicates removed.
    !   for performance reasons, the first 2 passes are taken
    !   out of the standard loop, and use dedicated coding.
    ! __________________________________________________________
    ! __________________________________________________________
    integer, dimension (:), intent (in)  :: xdont
    integer, dimension (:), intent (out) :: igoest
    ! __________________________________________________________
    integer :: xtst, xdona, xdonb
    !
    ! __________________________________________________________
    integer, dimension (size(igoest)) :: jwrkt, irngt
    integer :: lmtna, lmtnc, irng, irng1, irng2, nuni
    integer :: nval, iind, iwrkd, iwrk, iwrkf, jinda, iinda, iindb
    !
    nval = min (size(xdont), size(igoest))
    !
    select case (nval)
    case (:0)
       return
    case (1)
       igoest (1) = 1
       return
    case default
       continue
    end select
    !
    !  fill-in the index array, creating ordered couples
    !
    do iind = 2, nval, 2
       if (xdont(iind-1) < xdont(iind)) then
          irngt (iind-1) = iind - 1
          irngt (iind) = iind
       else
          irngt (iind-1) = iind
          irngt (iind) = iind - 1
       end if
    end do
    if (modulo (nval, 2) /= 0) then
       irngt (nval) = nval
    end if
    !
    !  we will now have ordered subsets a - b - a - b - ...
    !  and merge a and b couples into     c   -   c   - ...
    !
    lmtna = 2
    lmtnc = 4
    !
    !  first iteration. the length of the ordered subsets goes from 2 to 4
    !
    do
       if (nval <= 4) exit
       !
       !   loop on merges of a and b into c
       !
       do iwrkd = 0, nval - 1, 4
          if ((iwrkd+4) > nval) then
             if ((iwrkd+2) >= nval) exit
             !
             !   1 2 3
             !
             if (xdont(irngt(iwrkd+2)) <= xdont(irngt(iwrkd+3))) exit
             !
             !   1 3 2
             !
             if (xdont(irngt(iwrkd+1)) <= xdont(irngt(iwrkd+3))) then
                irng2 = irngt (iwrkd+2)
                irngt (iwrkd+2) = irngt (iwrkd+3)
                irngt (iwrkd+3) = irng2
                !
                !   3 1 2
                !
             else
                irng1 = irngt (iwrkd+1)
                irngt (iwrkd+1) = irngt (iwrkd+3)
                irngt (iwrkd+3) = irngt (iwrkd+2)
                irngt (iwrkd+2) = irng1
             end if
             exit
          end if
          !
          !   1 2 3 4
          !
          if (xdont(irngt(iwrkd+2)) <= xdont(irngt(iwrkd+3))) cycle
          !
          !   1 3 x x
          !
          if (xdont(irngt(iwrkd+1)) <= xdont(irngt(iwrkd+3))) then
             irng2 = irngt (iwrkd+2)
             irngt (iwrkd+2) = irngt (iwrkd+3)
             if (xdont(irng2) <= xdont(irngt(iwrkd+4))) then
                !   1 3 2 4
                irngt (iwrkd+3) = irng2
             else
                !   1 3 4 2
                irngt (iwrkd+3) = irngt (iwrkd+4)
                irngt (iwrkd+4) = irng2
             end if
             !
             !   3 x x x
             !
          else
             irng1 = irngt (iwrkd+1)
             irng2 = irngt (iwrkd+2)
             irngt (iwrkd+1) = irngt (iwrkd+3)
             if (xdont(irng1) <= xdont(irngt(iwrkd+4))) then
                irngt (iwrkd+2) = irng1
                if (xdont(irng2) <= xdont(irngt(iwrkd+4))) then
                   !   3 1 2 4
                   irngt (iwrkd+3) = irng2
                else
                   !   3 1 4 2
                   irngt (iwrkd+3) = irngt (iwrkd+4)
                   irngt (iwrkd+4) = irng2
                end if
             else
                !   3 4 1 2
                irngt (iwrkd+2) = irngt (iwrkd+4)
                irngt (iwrkd+3) = irng1
                irngt (iwrkd+4) = irng2
             end if
          end if
       end do
       !
       !  the cs become as and bs
       !
       lmtna = 4
       exit
    end do
    !
    !  iteration loop. each time, the length of the ordered subsets
    !  is doubled.
    !
    do
       if (2*lmtna >= nval) exit
       iwrkf = 0
       lmtnc = 2 * lmtnc
       !
       !   loop on merges of a and b into c
       !
       do
          iwrk = iwrkf
          iwrkd = iwrkf + 1
          jinda = iwrkf + lmtna
          iwrkf = iwrkf + lmtnc
          if (iwrkf >= nval) then
             if (jinda >= nval) exit
             iwrkf = nval
          end if
          iinda = 1
          iindb = jinda + 1
          !
          !  one steps in the c subset, that we create in the final rank array
          !
          !  make a copy of the rank array for the iteration
          !
          jwrkt (1:lmtna) = irngt (iwrkd:jinda)
          xdona = xdont (jwrkt(iinda))
          xdonb = xdont (irngt(iindb))
          !
          do
             iwrk = iwrk + 1
             !
             !  we still have unprocessed values in both a and b
             !
             if (xdona > xdonb) then
                irngt (iwrk) = irngt (iindb)
                iindb = iindb + 1
                if (iindb > iwrkf) then
                   !  only a still with unprocessed values
                   irngt (iwrk+1:iwrkf) = jwrkt (iinda:lmtna)
                   exit
                end if
                xdonb = xdont (irngt(iindb))
             else
                irngt (iwrk) = jwrkt (iinda)
                iinda = iinda + 1
                if (iinda > lmtna) exit! only b still with unprocessed values
                xdona = xdont (jwrkt(iinda))
             end if
             !
          end do
       end do
       !
       !  the cs become as and bs
       !
       lmtna = 2 * lmtna
    end do
    !
    !   last merge of a and b into c, with removal of duplicates.
    !
    iinda = 1
    iindb = lmtna + 1
    nuni = 0
    !
    !  one steps in the c subset, that we create in the final rank array
    !
    jwrkt (1:lmtna) = irngt (1:lmtna)
    if (iindb <= nval) then
       xtst = nearless (min(xdont(jwrkt(1)), xdont(irngt(iindb))))
    else
       xtst = nearless (xdont(jwrkt(1)))
    endif
    do iwrk = 1, nval
       !
       !  we still have unprocessed values in both a and b
       !
       if (iinda <= lmtna) then
          if (iindb <= nval) then
             if (xdont(jwrkt(iinda)) > xdont(irngt(iindb))) then
                irng = irngt (iindb)
                iindb = iindb + 1
             else
                irng = jwrkt (iinda)
                iinda = iinda + 1
             end if
          else
             !
             !  only a still with unprocessed values
             !
             irng = jwrkt (iinda)
             iinda = iinda + 1
          end if
       else
          !
          !  only b still with unprocessed values
          !
          irng = irngt (iwrk)
       end if
       if (xdont(irng) > xtst) then
          xtst = xdont (irng)
          nuni = nuni + 1
       end if
       igoest (irng) = nuni
       !
    end do
    !
    return
    !
  end subroutine i_uniinv

  function d_nearless (xval) result (d_nl)
    !  nearest value less than given value
    real (kind=8), intent (in) :: xval
    real (kind=8) :: d_nl
    d_nl = nearest (xval, -1.0d0)
    return
  end function d_nearless
  function r_nearless (xval) result (r_nl)
    !  nearest value less than given value
    real, intent (in) :: xval
    real :: r_nl
    r_nl = nearest (xval, -1.0)
    return
  end function r_nearless
  function i_nearless (xval) result (i_nl)
    !  nearest value less than given value
    integer, intent (in) :: xval
    integer :: i_nl
    i_nl = xval - 1
    return
  end function i_nearless

end module m_uniinv
