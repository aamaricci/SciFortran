MODULE CGFIT_ROUTINES
  implicit none

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



  integer                        :: ncom
  real(8), dimension(:), pointer :: pcom,xicom
  procedure(cgfit_func),pointer  :: func
  procedure(cgfit_fjac),pointer  :: fjac

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
  SUBROUTINE linmin(p,xi,fret,ftol)
    real(8), intent(out)                         :: fret
    real(8), dimension(:), target, intent(inout) :: p,xi
    real(8),optional                             :: ftol
    real(8)                                      :: tol
    real(8)                                      :: ax,bx,fa,fb,fx,xmin,xx
    tol=1.d-6;if(present(ftol))tol=ftol
    ncom=size(p) ; if(ncom /= size(xi))stop "Error in LinMin"
    pcom=>p
    xicom=>xi
    ax=0.d0
    xx=1.d0
    call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
    fret=brent_(ax,xx,bx,f1dim,tol,xmin)
    xi=xmin*xi
    p=p+xi
    return
  end subroutine linmin

  SUBROUTINE dlinmin(p,xi,fret,ftol,itmax)
    real(8), intent(out)                         :: fret
    real(8), dimension(:), target, intent(inout) :: p,xi
    real(8),optional                             :: ftol
    integer,optional                             :: itmax
    real(8)                                      :: tol
    integer                                      :: itmax_
    real(8)                                      :: ax,bx,fa,fb,fx,xmin,xx
    tol=1.d-6;if(present(ftol))tol=ftol
    itmax_=100;if(present(itmax))itmax_=itmax
    ncom=size(p) ; if(ncom /= size(xi))stop "Error in DLinMin"
    pcom=>p
    xicom=>xi
    ax=0.d0
    xx=1.d0
    call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
    fret=dbrent_(ax,xx,bx,f1dim,df1dim,tol,xmin,itmax_)
    xi=xmin*xi
    p=p+xi
    return
  end subroutine dlinmin


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

  function df1dim(x)
    real(8), intent(in)                :: x
    real(8)                            :: df1dim
    real(8), dimension(:), allocatable :: xt,df
    allocate(xt(ncom),df(ncom))
    xt(:)=pcom(:)+x*xicom(:)
    df(:)=fjac(xt)
    df1dim=dot_product(df,xicom)
    deallocate(xt,df)
  end function df1dim

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
    real(8), intent(out)   :: cx,fa,fb,fc
    !...the first parameter is the default ratio by which successive intervals
    !   are magnified; the second is the maximum magnification allowed for a
    !   parabolic-fit step
    real(8), parameter     :: gold=(1d0+sqrt(5d0))/2d0,glimit=10d0,tiny=1.d-20
    real(8)                :: fu,q,r,u,ulim
    interface
       function func(x)
         real(8), intent(in) :: x
         real(8)             :: func
       end function func
    end interface
    fa=func(ax)
    fb=func(bx)
    if (fb > fa) then           !switch the role of a and b to that we can go downhill
       call swap(ax,bx)
       call swap(fa,fb)
    end if
    cx=bx+GOLD*(bx-ax)
    fc=func(cx)
    ! do while (fb.ge.fc)
    do                          !do while”: keep returning here until we bracket.
       ! print*,"mnbrak count",ax,cx,bx
       if(isnan(ax).OR.isnan(cx).OR.isnan(bx))stop "MNBRAK error: ax/bx/cx are Nan!"
       if (fb < fc) RETURN
       r=(bx-ax)*(fb-fc)        !Compute u by parabolic extrapolation from a, b, c. TINY
       q=(bx-cx)*(fb-fa)        !is used to prevent any possible division by zero
       u=bx-((bx-cx)*q-(bx-ax)*r)/(2d0*sign(max(abs(q-r),TINY),q-r))
       ulim=bx+GLIMIT*(cx-bx)
       if ((bx-u)*(u-cx) > 0.d0) then !Parabolic u is between b and c: try it.
          fu=func(u)
          if (fu < fc) then     !Got a minimum between b and c.
             ax=bx
             fa=fb
             bx=u
             fb=fu
             RETURN
          else if (fu > fb) then !Got a minimum between between a and u.
             cx=u
             fc=fu
             RETURN
          end if
          u=cx+GOLD*(cx-bx)     !Parabolic fit was no use. Use default magnification.
          fu=func(u)
       else if ((cx-u)*(u-ulim) > 0.d0) then !Parabolic fit is between c and its allowed limit
          fu=func(u)
          if (fu < fc) then
             bx=cx
             cx=u
             u=cx+GOLD*(cx-bx)
             call shft(fb,fc,fu,func(u))
          end if
       else if ((u-ulim)*(ulim-cx) >= 0.d0) then !Limit parabolic u to maximum allowed value.
          u=ulim
          fu=func(u)
       else                     !Reject parabolic u, use default magnification.
          u=cx+GOLD*(cx-bx)
          fu=func(u)
       end if                   !Eliminate oldest point and continue.
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
    real(8), parameter   :: zeps=1d0-3d0*epsilon(ax)
    real(8), parameter   :: cgold=(3d0-sqrt(5d0))/2d0 !0.3819660d0
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
    e=0.d0
    fx=func(x)
    fv=fx
    fw=fx
    do iter=1,itmax
       xm=0.5d0*(a+b)
       tol1=tol*abs(x)+zeps
       tol2=2.0*tol1
       if (abs(x-xm) <= (tol2-0.5d0*(b-a))) then
          xmin=x
          brent_=fx
          return
       end if
       if (abs(e) > tol1) then
          r=(x-w)*(fx-fv)
          q=(x-v)*(fx-fw)
          p=(x-v)*q-(x-w)*r
          q=2.d0*(q-r)
          if (q > 0.d0) p=-p
          q=abs(q)
          etemp=e
          e=d
          if (abs(p) >= abs(0.5d0*q*etemp) .or. &
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
    ! pause 'brent: exceed maximum iterations'
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




  function dbrent_(ax,bx,cx,func,fjac,tol,xmin,itmax) result(dbrent)
    real(8),intent(in)  :: ax,bx,cx,tol
    real(8),intent(out) :: xmin
    real(8)             :: dbrent
    integer,optional    :: itmax
    integer             :: itmax_=100
    real(8), parameter  :: zeps=1.d-3*epsilon(ax)
    integer             :: iter
    real(8)             :: a,b,d,d1,d2,du,dv,dw,dx,e,fu,fv,fw,fx,olde,tol1,tol2
    real(8)             :: u,u1,u2,v,w,x,xm
    logical             :: ok1,ok2
    interface
       function func(x)
         real(8), intent(in) :: x
         real(8)             :: func
       end function func
       function fjac(x)
         real(8), intent(in) :: x
         real(8)             :: fjac
       end function fjac
    end interface
    itmax_=100;if(present(itmax))itmax_=itmax    
    a=min(ax,cx)
    b=max(ax,cx)
    v=bx
    w=v
    x=v
    e=0.d0
    fx=func(x)
    fv=fx
    fw=fx
    dx=fjac(x)
    dv=dx
    dw=dx
    do iter=1,ITMAX
       xm=0.5d0*(a+b)
       tol1=tol*abs(x)+ZEPS
       tol2=2.0d0*tol1
       if (abs(x-xm) <= (tol2-0.5d0*(b-a))) exit
       if (abs(e) > tol1) then
          d1=2.0d0*(b-a)
          d2=d1
          if (dw /= dx) d1=(w-x)*dx/(dx-dw)
          if (dv /= dx) d2=(v-x)*dx/(dx-dv)
          u1=x+d1
          u2=x+d2
          ok1=((a-u1)*(u1-b) > 0.d0) .and. (dx*d1 <= 0.d0)
          ok2=((a-u2)*(u2-b) > 0.d0) .and. (dx*d2 <= 0.d0)
          olde=e
          e=d
          if (ok1 .or. ok2) then
             if (ok1 .and. ok2) then
                d=merge(d1,d2, abs(d1) < abs(d2))
             else
                d=merge(d1,d2,ok1)
             end if
             if (abs(d) <= abs(0.5d0*olde)) then
                u=x+d
                if (u-a < tol2 .or. b-u < tol2) &
                     d=sign(tol1,xm-x)
             else
                e=merge(a,b, dx >= 0.d0)-x
                d=0.5d0*e
             end if
          else
             e=merge(a,b, dx >= 0.d0)-x
             d=0.5d0*e
          end if
       else
          e=merge(a,b, dx >= 0.d0)-x
          d=0.5d0*e
       end if
       if (abs(d) >= tol1) then
          u=x+d
          fu=func(u)
       else
          u=x+sign(tol1,d)
          fu=func(u)
          if (fu > fx) exit
       end if
       du=fjac(u)
       if (fu <= fx) then
          if (u >= x) then
             a=x
          else
             b=x
          end if
          call mov3(v,fv,dv,w,fw,dw)
          call mov3(w,fw,dw,x,fx,dx)
          call mov3(x,fx,dx,u,fu,du)
       else
          if (u < x) then
             a=u
          else
             b=u
          end if
          if (fu <= fw .or. w == x) then
             call mov3(v,fv,dv,w,fw,dw)
             call mov3(w,fw,dw,u,fu,du)
          else if (fu <= fv .or. v == x .or. v == w) then
             call mov3(v,fv,dv,u,fu,du)
          end if
       end if
    end do
    if (iter > ITMAX) stop 'dbrent: exceeded maximum iterations'
    xmin=x
    dbrent=fx
  contains
    !bl
    subroutine mov3(a,b,c,d,e,f)
      real(8), intent(in) :: d,e,f
      real(8), intent(out) :: a,b,c
      a=d
      b=e
      c=f
    end subroutine mov3
  end function dbrent_




  function isinfty(a) result(bool)
    real(8) :: a
    logical :: bool
    bool = (a-1 == a)
  end function isinfty

  function isnan(a) result(bool)
    real(8) :: a
    logical :: bool
    bool = (a /= a)
  end function isnan


end module CGFIT_ROUTINES

