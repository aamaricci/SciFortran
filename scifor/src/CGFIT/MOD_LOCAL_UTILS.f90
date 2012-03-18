MODULE CGFIT_LOCAL_UTILS
  !USE FUNCTION
  implicit none

contains

  !+-------------------------------------------------------------------+
  !PURPOSE : Given a function FUNC, and given distinct initial points 
  ! AX and BX, this routine searches in the downhill direction (defined 
  ! by the function as evaluated at the initial points) and returns new 
  ! points AX, BX, CX which bracket a minimum of the function.  
  ! Also returned are the function values at the three points: FA, FB, FC.
  !+-------------------------------------------------------------------+
  SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func)
    REAL(8), INTENT(INOUT) :: ax,bx
    REAL(8), INTENT(OUT)   :: cx,fa,fb,fc
    !...The first parameter is the default ratio by which successive intervals
    !   are magnified; the second is the maximum magnification allowed for a
    !   parabolic-fit step
    REAL(8), PARAMETER     :: GOLD=1.618034d0,GLIMIT=100.0d0,TINY=1.0d-20
    REAL(8)                :: fu,q,r,u,ulim
    INTERFACE
       FUNCTION func(x)
         REAL(8), INTENT(IN) :: x
         REAL(8)             :: func
       END FUNCTION func
    END INTERFACE
    fa=func(ax)
    fb=func(bx)
    if (fb > fa) then
       call swap(ax,bx)
       call swap(fa,fb)
    end if
    cx=bx+GOLD*(bx-ax)
    fc=func(cx)
    do
       if (fb < fc) RETURN
       r=(bx-ax)*(fb-fc)
       q=(bx-cx)*(fb-fa)
       u=bx-((bx-cx)*q-(bx-ax)*r)/(2.d0*sign(max(abs(q-r),TINY),q-r))
       ulim=bx+GLIMIT*(cx-bx)
       if ((bx-u)*(u-cx) > 0.0) then
          fu=func(u)
          if (fu < fc) then
             ax=bx
             fa=fb
             bx=u
             fb=fu
             RETURN
          else if (fu > fb) then
             cx=u
             fc=fu
             RETURN
          end if
          u=cx+GOLD*(cx-bx)
          fu=func(u)
       else if ((cx-u)*(u-ulim) > 0.0) then
          fu=func(u)
          if (fu < fc) then
             bx=cx
             cx=u
             u=cx+GOLD*(cx-bx)
             call shft(fb,fc,fu,func(u))
          end if
       else if ((u-ulim)*(ulim-cx) >= 0.0) then
          u=ulim
          fu=func(u)
       else
          u=cx+GOLD*(cx-bx)
          fu=func(u)
       end if
       call shft(ax,bx,cx,u)
       call shft(fa,fb,fc,fu)
    end do
  CONTAINS
    SUBROUTINE swap(a,b)
      REAL(8), INTENT(INOUT) :: a,b
      REAL(8) :: dum
      dum=a
      a=b
      b=dum
    END SUBROUTINE swap
    !-------------------
    !-------------------
    !-------------------
    SUBROUTINE shft(a,b,c,d)
      REAL(8), INTENT(OUT) :: a
      REAL(8), INTENT(INOUT) :: b,c
      REAL(8), INTENT(IN) :: d
      a=b
      b=c
      c=d
    END SUBROUTINE shft
  END SUBROUTINE mnbrak



  !********************************************************************
  !********************************************************************
  !********************************************************************




  !+-------------------------------------------------------------------+
  !PURPOSE  : Given a function F, and given a bracketing triplet of 
  ! abscissas AX, BX, CX (such that BX is between AX and CX, 
  ! and F(BX) is less than both F(AX) and F(CX)), this routine 
  ! isolates the minimum to a fractional precision of about TOL 
  ! using Brent's method. The abscissa of the minimum is returned 
  ! as XMIN, and the minimum function value is returned as BRENT, 
  ! the returned function value.
  !+-------------------------------------------------------------------+
  FUNCTION brent(ax,bx,cx,func,tol,xmin)
    REAL(8), INTENT(IN)  :: ax,bx,cx,tol
    REAL(8), INTENT(OUT) :: xmin
    REAL(8)              :: brent
    INTEGER, PARAMETER   :: ITMAX=100
    REAL(8), PARAMETER   :: CGOLD=0.3819660d0,ZEPS=1.0d-3*epsilon(ax)
    INTEGER              :: iter
    REAL(8)              :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
    INTERFACE
       FUNCTION func(x)
         REAL(8), INTENT(IN) :: x
         REAL(8)             :: func
       END FUNCTION func
    END INTERFACE
    a=min(ax,cx)
    b=max(ax,cx)
    v=bx
    w=v
    x=v
    e=0.d0
    fx=func(x)
    fv=fx
    fw=fx
    do iter=1,ITMAX
       xm=0.5d0*(a+b)
       tol1=tol*abs(x)+ZEPS
       tol2=2.d0*tol1
       if (abs(x-xm) <= (tol2-0.5d0*(b-a))) then
          xmin=x
          brent=fx
          RETURN
       end if
       if (abs(e) > tol1) then
          r=(x-w)*(fx-fv)
          q=(x-v)*(fx-fw)
          p=(x-v)*q-(x-w)*r
          q=2.0*(q-r)
          if (q > 0.d0) p=-p
          q=abs(q)
          etemp=e
          e=d
          if (abs(p) >= abs(0.5*q*etemp) .or. &
               p <= q*(a-x) .or. p >= q*(b-x)) then
             e=merge(a-x,b-x, x >= xm )
             d=CGOLD*e
          else
             d=p/q
             u=x+d
             if (u-a < tol2 .or. b-u < tol2) d=sign(tol1,xm-x)
          end if
       else
          e=merge(a-x,b-x, x >= xm )
          d=CGOLD*e
       end if
       u=merge(x+d,x+sign(tol1,d), abs(d) >= tol1 )
       fu=func(u)
       if (fu <= fx) then
          if (u >= x) then
             a=x
          else
             b=x
          end if
          call shft(v,w,x,u)
          call shft(fv,fw,fx,fu)
       else
          if (u < x) then
             a=u
          else
             b=u
          end if
          if (fu <= fw .or. w == x) then
             v=w
             fv=fw
             w=u
             fw=fu
          else if (fu <= fv .or. v == x .or. v == w) then
             v=u
             fv=fu
          end if
       end if
    end do
    print*,'brent: exceed maximum iterations'
  CONTAINS
    SUBROUTINE shft(a,b,c,d)
      REAL(8), INTENT(OUT) :: a
      REAL(8), INTENT(INOUT) :: b,c
      REAL(8), INTENT(IN) :: d
      a=b
      b=c
      c=d
    END SUBROUTINE shft
  END FUNCTION brent
  !********************************************************************
  !********************************************************************
  !********************************************************************
end module CGFIT_LOCAL_UTILS
