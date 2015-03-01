module SF_MISC
  implicit none
  private
  real(8),parameter    :: pi    = 3.14159265358979323846264338327950288419716939937510d0

  !UNIINV = Merge-sort inverse ranking of an array, with removal of duplicate entries. 
  interface uniinv
     module procedure d_uniinv, r_uniinv, i_uniinv
  end interface uniinv

  !UNISTA = (stable unique) removes duplicates from an array, leaving unique entries
  ! in the order of their first appearance in the initial set.
  interface unista
     module procedure d_unista, r_unista, i_unista
  end interface unista


  interface uniq
     module procedure i_uniq,d_uniq
  end interface uniq


  interface nearless
     module procedure D_nearless, R_nearless, I_nearless
  end interface nearless


  !SORT,UNIQ,SHUFFLE:
  public :: sort
  public :: sort_array
  public :: reshuffle
  public :: uniq
  public :: uniinv
  public :: unista

contains



  ! SORTING 1D:
  !###################################################################
  !+-----------------------------------------------------------------+
  !PURPOSE  :   
  !+-----------------------------------------------------------------+
  subroutine i_uniq(AIN,AOUT,MASK)
    integer,dimension(:),intent(INOUT)                    :: AIN
    integer,dimension(:),allocatable,intent(OUT)          :: AOUT
    integer                                               :: NDIM
    logical,dimension(:),allocatable,intent(OUT),optional :: MASK
    if(present(MASK))then
       allocate(MASK(size(AIN)))
       call i_unista(AIN,ndim,MASK)
    else
       call i_unista(AIN,ndim)
    endif
    allocate(AOUT(NDIM))
    AOUT(1:NDIM)=AIN(1:NDIM)
  end subroutine i_uniq
  subroutine d_uniq(AIN,AOUT,MASK)
    real(8),dimension(:),intent(INOUT)                    :: AIN
    real(8),dimension(:),allocatable,intent(OUT)          :: AOUT
    integer                                               :: NDIM
    logical,dimension(:),allocatable,intent(OUT),optional :: MASK
    if(present(MASK))then
       allocate(MASK(size(AIN)))
       call d_unista(AIN,ndim,MASK)
    else
       call d_unista(AIN,ndim)
    endif
    allocate(AOUT(NDIM))
    AOUT(1:NDIM)=AIN(1:NDIM)
  end subroutine d_uniq



  !+-----------------------------------------------------------------+
  !PURPOSE  :   
  !+-----------------------------------------------------------------+
  subroutine reshuffle(array,order) 
    real(8),dimension(:)                    :: array
    integer,dimension(size(array))          :: order
    real(8),dimension(size(array))          :: dummy
    integer                                 :: i,Lk
    Lk=size(array)
    forall(i=1:Lk)dummy(i)=array(order(i))
    array=dummy
  end subroutine reshuffle




  !+-----------------------------------------------------------------+
  !PURPOSE  :   
  !+-----------------------------------------------------------------+
  subroutine sort(ARR,M)
    implicit none
    integer :: i,j,M
    real(8) :: a
    real(8),dimension(M) :: ARR
    do j=2, M
       a=ARR(j)
       do i=j-1,1,-1
          if (ARR(i)<=a) goto 10
          ARR(i+1)=ARR(i)
       enddo
       i=0
10     ARR(i+1)=a
    enddo
    return
  end subroutine sort




  !+------------------------------------------------------------------+
  !PURPOSE  : Sort an array, gives the new ordering of the label.
  !+------------------------------------------------------------------+
  subroutine sort_array(array,order2,no_touch_array)
    implicit none
    real(8),dimension(:)                    :: array
    real(8),dimension(size(array))          :: backup
    integer,dimension(size(array))          :: order
    integer,dimension(size(array)),optional :: order2
    integer                                 :: i
    logical,optional                        :: no_touch_array
    do i=1,size(order)
       order(i)=i
    enddo
    call qsort_sort( array, order, 1, size(array) )
    if(.not.present(no_touch_array))then
       do i=1,size(order)
          backup(i)=array(order(i))
       enddo
       array=backup
    endif
    if(present(order2)) order2=order
  contains
    recursive subroutine qsort_sort( array, order, left, right )
      implicit none
      real(8), dimension(:)                 :: array
      integer, dimension(:)                 :: order
      integer                               :: left
      integer                               :: right
      integer                               :: i
      integer                               :: last
      if ( left .ge. right ) return
      call qsort_swap( order, left, qsort_rand(left,right) )
      last = left
      do i = left+1, right
         if ( compare(array(order(i)), array(order(left)) ) .lt. 0 ) then
            last = last + 1
            call qsort_swap( order, last, i )
         endif
      enddo
      call qsort_swap( order, left, last )
      call qsort_sort( array, order, left, last-1 )
      call qsort_sort( array, order, last+1, right )
    end subroutine qsort_sort
    !---------------------------------------------!
    subroutine qsort_swap( order, first, second )
      implicit none
      integer, dimension(:)                 :: order
      integer                               :: first, second
      integer                               :: tmp
      tmp           = order(first)
      order(first)  = order(second)
      order(second) = tmp
    end subroutine qsort_swap
    !---------------------------------------------!
    integer function qsort_rand( lower, upper )
      implicit none
      integer                               :: lower, upper
      real(4)                               :: r
      r=drand()
      qsort_rand =  lower + nint(r * (upper-lower))
    end function qsort_rand
    !---------------------------------------------!
    function drand()
      implicit none
      real(8)                               :: drand
      real(4)                               :: r
      call random_number(r)
      drand=dble(r)
    end function drand
    !---------------------------------------------!
    function compare(f,g)
      implicit none
      real(8)                               :: f,g
      integer                               :: compare
      if(f<g) then
         compare=-1
      else
         compare=1
      endif
    end function compare
  end subroutine sort_array





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
  !
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




  !*******************************************************************
  !*******************************************************************
  !*******************************************************************




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




  !*******************************************************************
  !*******************************************************************
  !*******************************************************************



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



end module SF_MISC
