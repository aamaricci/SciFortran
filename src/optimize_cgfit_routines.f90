module CGFIT_FUNC_INTERFACE
  abstract interface
     function cgfit_func(a)
       real(8),dimension(:)  ::  a
       real(8)               ::  cgfit_func
     end function cgfit_func
     function cgfit_fjac(a)
       real(8),dimension(:)       :: a
       real(8),dimension(size(a)) :: cgfit_fjac
     end function cgfit_fjac
  end interface
end module CGFIT_FUNC_INTERFACE




!********************************************************************
!********************************************************************
!********************************************************************
!********************************************************************
!********************************************************************
!********************************************************************




MODULE CGFIT_ROUTINES
  USE CGFIT_FUNC_INTERFACE
  implicit none
  integer                        :: ncom
  real(8), dimension(:), pointer :: pcom,xicom
  procedure(cgfit_func),pointer :: func
  procedure(cgfit_fjac),pointer :: dfunc

contains

  !+-------------------------------------------------------------------+
  ! PURPOSE  : Minimization routine
  ! Given an N dimensional point P and an N dimensional direction 
  ! XI, LINMIN moves and resets P to where the function FUNC(P) takes 
  ! on a minimum
  ! along the direction XI from P, and replaces XI by the actual vector
  ! displacement that P was moved.  Also returns FRET the value of 
  ! FUNC at the returned location P.  
  ! This is actually all accomplished by calling the routines 
  ! MNBRAK and BRENT.
  !+-------------------------------------------------------------------+
  SUBROUTINE linmin(p,xi,fret)
    real(8), intent(out) :: fret
    real(8), dimension(:), target, intent(inout) :: p,xi
    real(8), parameter :: tol=1.0e-4
    real(8) :: ax,bx,fa,fb,fx,xmin,xx
    ncom=size(p) ; if(ncom /= size(xi))stop "Error in LinMin"
    pcom=>p
    xicom=>xi
    ax=0.0
    xx=1.0
    call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
    fret=brent_(ax,xx,bx,f1dim,TOL,xmin)
    !...construct the vector results to return
    xi=xmin*xi
    p=p+xi
    return
  end subroutine linmin



  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  function f1dim(x)
    real(8), intent(in)                :: x
    real(8)                            :: f1dim
    real(8), dimension(:), allocatable :: xt
    allocate(xt(ncom))
    xt(:)=pcom(:)+x*xicom(:)
    f1dim=func(xt)
    deallocate(xt)
  end function f1dim


  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !     Given a function FUNC, and given distinct initial points AX and BX,
  !     this routine searches in the downhill direction (defined by the 
  !     function as evaluated at the initial points) and returns new points
  !     AX, BX, CX which bracket a minimum of the function.  
  !     Also returned are the function values at the three points, 
  !     FA, FB, and FC.
  !+-------------------------------------------------------------------+
  subroutine mnbrak(ax,bx,cx,fa,fb,fc,func)
    real(8), intent(inout) :: ax,bx
    real(8), intent(out) :: cx,fa,fb,fc
    !...the first parameter is the default ratio by which successive intervals
    !   are magnified; the second is the maximum magnification allowed for a
    !   parabolic-fit step
    real(8), parameter :: gold=1.618034,glimit=100.0,tiny=1.0e-20
    real(8) :: fu,q,r,u,ulim
    interface
       function func(x)
         real(8), intent(in) :: x
         real(8)             :: func
       end function func
    end interface
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
       u=bx-((bx-cx)*q-(bx-ax)*r)/(2.0*sign(max(abs(q-r),TINY),q-r))
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
  contains
    subroutine swap(a,b)
      real(8), intent(inout) :: a,b
      real(8) :: dum
      dum=a
      a=b
      b=dum
    end subroutine swap
    !-------------------
    subroutine shft(a,b,c,d)
      real(8), intent(out) :: a
      real(8), intent(inout) :: b,c
      real(8), intent(in) :: d
      a=b
      b=c
      c=d
    end subroutine shft
  end subroutine mnbrak


  !+-------------------------------------------------------------------+
  !purpose  : 
  !     given a function f, and given a bracketing triplet of 
  !     abscissas ax, bx, cx (such that bx is between ax and cx, 
  !     and f(bx) is less than both f(ax) and f(cx)), this routine 
  !     isolates the minimum to a fractional precision of about tol 
  !     using brents method. the abscissa of the minimum is returned 
  !     as xmin, and the minimum function value is returned as brent, 
  !     the returned function value.
  !+-------------------------------------------------------------------+
  function brent_(ax,bx,cx,func,tol,xmin)
    real(8), intent(in)  :: ax,bx,cx,tol
    real(8), intent(out) :: xmin
    real(8)              :: brent_
    integer, parameter   :: itmax=100
    real(8), parameter   :: cgold=0.3819660d0,zeps=1.0d-3*epsilon(ax)
    integer              :: iter
    real(8)              :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
    interface
       function func(x)
         real(8), intent(in) :: x
         real(8)             :: func
       end function func
    end interface
    a=min(ax,cx)
    b=max(ax,cx)
    v=bx
    w=v
    x=v
    e=0.0
    fx=func(x)
    fv=fx
    fw=fx
    do iter=1,itmax
       xm=0.5*(a+b)
       tol1=tol*abs(x)+zeps
       tol2=2.0*tol1
       if (abs(x-xm) <= (tol2-0.5*(b-a))) then
          xmin=x
          brent_=fx
          return
       end if
       if (abs(e) > tol1) then
          r=(x-w)*(fx-fv)
          q=(x-v)*(fx-fw)
          p=(x-v)*q-(x-w)*r
          q=2.0*(q-r)
          if (q > 0.0) p=-p
          q=abs(q)
          etemp=e
          e=d
          if (abs(p) >= abs(0.5*q*etemp) .or. &
               p <= q*(a-x) .or. p >= q*(b-x)) then
             e=merge(a-x,b-x, x >= xm )
             d=cgold*e
          else
             d=p/q
             u=x+d
             if (u-a < tol2 .or. b-u < tol2) d=sign(tol1,xm-x)
          end if
       else
          e=merge(a-x,b-x, x >= xm )
          d=cgold*e
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
    !pause 'brent: exceed maximum iterations'
  contains
    subroutine shft(a,b,c,d)
      real(8), intent(out) :: a
      real(8), intent(inout) :: b,c
      real(8), intent(in) :: d
      a=b
      b=c
      c=d
    end subroutine shft
  end function brent_

end module CGFIT_ROUTINES

