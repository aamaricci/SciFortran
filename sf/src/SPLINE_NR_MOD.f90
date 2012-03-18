module SPLINE_NR_MOD
  implicit none

  INTERFACE assert_eq
     MODULE PROCEDURE assert_eq2,assert_eq3,assert_eq4,assert_eqn
  END INTERFACE assert_eq


contains

  function locate(xx,x)
    REAL(8), DIMENSION(:), INTENT(IN) :: xx
    REAL(8), INTENT(IN) :: x
    INTEGER :: locate
    INTEGER :: n,jl,jm,ju
    LOGICAL :: ascnd
    n=size(xx)
    ascnd = (xx(n) >= xx(1))
    jl=0
    ju=n+1
    do
       if (ju-jl <= 1) exit
       jm=(ju+jl)/2
       if (ascnd .eqv. (x >= xx(jm))) then
          jl=jm
       else
          ju=jm
       end if
    end do
    if (x == xx(1)) then
       locate=1
    else if (x == xx(n)) then
       locate=n-1
    else
       locate=jl
    end if
  END FUNCTION locate




  SUBROUTINE polint(xa,ya,x,y,dy)
    REAL(8), DIMENSION(:), INTENT(IN) :: xa,ya
    REAL(8), INTENT(IN)          :: x
    REAL(8), INTENT(OUT)         :: y,dy
    INTEGER                      :: m,n,ns
    REAL(8), DIMENSION(size(xa)) :: c,d,den,ho
    n=assert_eq(size(xa),size(ya),'polint')
    c=ya
    d=ya
    ho=xa-x
    ns=iminloc(abs(x-xa))
    y=ya(ns)
    ns=ns-1
    do m=1,n-1
       den(1:n-m)=ho(1:n-m)-ho(1+m:n)
       if (any(den(1:n-m) == 0.0))then
          print*,'polint: calculation failure'
          stop
       endif
       den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
       d(1:n-m)=ho(1+m:n)*den(1:n-m)
       c(1:n-m)=ho(1:n-m)*den(1:n-m)
       if (2*ns < n-m) then
          dy=c(ns+1)
       else
          dy=d(ns)
          ns=ns-1
       end if
       y=y+dy
    end do
  END SUBROUTINE polint


  !   subroutine locate(xx,n,x,j) 
  !     INTEGER j,n 
  !     DOUBLE PRECISION x,xx(n) 
  !     !Given an array xx(1:n), and given a value x, returns a value j such that x is between xx(j) and xx(j+1). xx(1:n) must be monotonic, either increasing or decreasing. j=0 or j=n is returned to indicate that x is out of range.

  !     INTEGER jl,jm,ju
  !     jl=0                      !Initialize lower 
  !     ju=n+1                    !and upper limits. 
  ! 10  if(ju-jl.gt.1)then        ! If we are not yet done,
  !        jm=(ju+jl)/2           !compute a midpoint,
  !        if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm)))then
  !           jl=jm               !and replace either the lower limit
  !        else
  !           ju=jm               !or the upper limit, as appropriate. 
  !        endif
  !        goto 10                !Repeat until
  !     endif                     !the test condition 10 is satis ed.
  !     if(x.eq.xx(1))then        !Then set the output
  !        j=1
  !     else if(x.eq.xx(n))then 
  !        j=n-1
  !     else 
  !        j=jl
  !     endif
  !     return                    !and return.
  !   end subroutine locate
  !   !-------------------------------------------------------!


  ! !-------------------------------------------------------!
  ! subroutine polint(xa,ya,n,x,y,dy) !Num Recipes
  !   INTEGER n,NMAX
  !   DOUBLE PRECISION dy,x,y,xa(n),ya(n)
  !   PARAMETER (NMAX=10)
  !   !     Largest anticipated value of n. Given arrays xa and ya, each of length n,
  !   !     and given a value x, this routine returns a value y, and an error estimate
  !   !     dy. If P(x) is the polynomial of degree N   1 such that P(xai) = yai, i = 1,
  !   !     . . . , n, then the returned value y = P(x).
  !   INTEGER i,m,ns
  !   DOUBLE PRECISION den,dif,dift,ho,hp,w,c(NMAX),d(NMAX) 
  !   ns=1 
  !   dif=abs(x-xa(1))
  !   do i=1,n
  !      dift=abs(x-xa(i))
  !      if (dift.lt.dif) then
  !         ns=i
  !         dif=dift
  !      endif
  !      c(i)=ya(i)             !and initialize the tableau of c s and d s. 
  !      d(i)=ya(i)
  !   enddo
  !   y=ya(ns)
  !   ! This is the initial approximation to y.
  !   ns=ns-1 
  !   !For each column of the tableau,
  !   do m=1,n-1
  !      do i=1,n-m             !we loop over the current c s and d s and update them.
  !         ho=xa(i)-x 
  !         hp=xa(i+m)-x 
  !         w=c(i+1)-d(i)
  !         den=ho-hp 
  !         !            if(den.eq.0.)pause  'failure in polint'
  !         if(den.eq.0.) write(*,*) 'failure in polint at x',x
  !         den=w/den
  !         d(i)=hp*den
  !         !     Here the c s and d s are updated.
  !         c(i)=ho*den
  !      enddo
  !      if (2*ns.lt.n-m)then
  !         ! After each column in the tableau is completed, we decide which correction, c or d,
  !         ! we want to add to our accumulating value of y, i.e., which path to take through the tableau
  !         ! forking up or down. We do this in such a way as to take the most  straight line  route through
  !         ! the tableau to its apex, updating ns accordingly to keep track of where we are. This route keeps
  !         ! the partial approximations centered (insofar as possible) on the target x. T he last dy added is
  !         ! thus the error indication. 
  !         dy=c(ns+1) 
  !      else
  !         dy=d(ns)
  !         ns=ns-1
  !      endif
  !      y=y+dy 
  !   enddo
  !   return
  ! END SUBROUTINE polint
  ! !-------------------------------------------------------!

  FUNCTION iminloc(arr)
    REAL(8), DIMENSION(:), INTENT(IN) :: arr
    INTEGER, DIMENSION(1) :: imin
    INTEGER :: iminloc
    imin=minloc(arr(:))
    iminloc=imin(1)
  END FUNCTION iminloc

  FUNCTION assert_eq2(n1,n2,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2
    INTEGER :: assert_eq2
    if (n1 == n2) then
       assert_eq2=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq2'
    end if
  END FUNCTION assert_eq2

  FUNCTION assert_eq3(n1,n2,n3,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2,n3
    INTEGER :: assert_eq3
    if (n1 == n2 .and. n2 == n3) then
       assert_eq3=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq3'
    end if
  END FUNCTION assert_eq3

  FUNCTION assert_eq4(n1,n2,n3,n4,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2,n3,n4
    INTEGER :: assert_eq4
    if (n1 == n2 .and. n2 == n3 .and. n3 == n4) then
       assert_eq4=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq4'
    end if
  END FUNCTION assert_eq4

  FUNCTION assert_eqn(nn,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, DIMENSION(:), INTENT(IN) :: nn
    INTEGER :: assert_eqn
    if (all(nn(2:) == nn(1))) then
       assert_eqn=nn(1)
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eqn'
    end if
  END FUNCTION assert_eqn

end module SPLINE_NR_MOD
