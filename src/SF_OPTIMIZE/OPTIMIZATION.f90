MODULE OPTIMIZE_MINIMIZE
  USE CGFIT_ROUTINES
  implicit none
  private

  interface fmin_cg
     module procedure fmin_cg_df,fmin_cg_f
  end interface fmin_cg

  interface fmin_cgplus
     module procedure fmin_cgplus_df,fmin_cgplus_f
  end interface fmin_cgplus

  interface fmin_cgminimize
     module procedure fmin_cgminimize_func,fmin_cgminimize_sub
  end interface fmin_cgminimize

  interface leastsq
     module procedure :: leastsq_lmdif_func
     module procedure :: leastsq_lmdif_sub
     module procedure :: leastsq_lmder_func
     module procedure :: leastsq_lmder_sub
  end interface leastsq

  interface curvefit
     module procedure :: curvefit_lmdif_func
     module procedure :: curvefit_lmdif_sub
     module procedure :: curvefit_lmder_func
     module procedure :: curvefit_lmder_sub
  end interface curvefit

  interface dbrent
     module procedure :: dbrent_wgrad
     module procedure :: dbrent_nograd
  end interface dbrent

  interface fmin_bfgs
     module procedure :: bfgs_with_grad
     module procedure :: bfgs_no_grad
  end interface fmin_bfgs


  !General-purpose
  public   :: fmin         !Minimize a function using the downhill simplex algorithm.
  ! public :: fmin_powell  !Minimize a function using modified Powell’s method. This method
  public   :: fmin_cg
  public   :: fmin_cgplus
  public   :: fmin_cgminimize
  public :: fmin_bfgs    !Minimize a function using the BFGS algorithm.
  ! public :: fmin_ncg     !Unconstrained minimization of a function using the Newton-CG method.
  public   :: leastsq      !Minimize the sum of squares of a set of equations. a wrapper around MINPACKs lmdif and lmder algorithms.
  public   :: curvefit     !Use non-linear least squares to fit a function, f, to data.

  !Global
  ! public :: anneal       !Minimize a function using simulated annealing.
  ! public :: basinhopping ! Find the global minimum of a function using the basin-hopping algorithm ..


  !Scalar function minimizers
  public   :: brent         !minimize a given a function of one-variable with a possible bracketing interval without using derivative information
  public   :: dbrent        !minimize a given a function of one-variable with a possible bracketing interval  using derivative information
  public   :: bracket       !Bracket the minimum of the function.


  real(8)  :: df_eps=0d0


contains




  !##################################################################
  !##################################################################
  !                       SCALAR FUNCTIONS
  !##################################################################
  !##################################################################

  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  ! Given a function f, and given a bracketing triplet of abscissas
  ! ax, bx, cx (such that bx is between ax and cx, and f(bx) is less
  ! than both f(ax) and f(cx)), this routine isolates the minimum to a
  ! fractional precision of about tol using Brent’s method. The abscissa
  ! of the minimum is returned as xmin, and the minimum function value
  ! is returned as brent, the returned function value.
  ! Parameters: Maximum allowed number of iterations; golden ratio;
  ! and a small number that protects against trying to achieve
  ! fractional accuracy for a minimum that happens to be exactly zero.
  !+-------------------------------------------------------------------+
  subroutine brent(func,xmin,brack,tol,niter)
    interface
       function func(x)
         real(8) :: x
         real(8) :: func
       end function func
    end interface
    real(8),intent(inout)         :: xmin
    real(8),dimension(:),optional :: brack
    real(8),optional              :: tol
    integer,optional              :: niter
    real(8)                       :: tol_
    integer                       :: niter_
    integer                       :: iter
    real(8)                       :: ax,xx,bx,fa,fx,fb,fret
    !
    tol_=1d-9;if(present(tol))tol_=tol
    Niter_=200;if(present(Niter))Niter_=Niter
    !
    if(present(brack))then
       select case(size(brack))
       case(1)
          stop "Brent error: calling brent with size(brack)==1. None or two points are necessary."
       case(2)
          ax = brack(1)
          xx = brack(2)
          call bracket(ax,xx,bx,fa,fx,fb,func)
       case (3)
          ax = brack(1)
          xx = brack(2)
          bx = brack(3)
       end select
    else
       ax=0d0
       xx=1d0
       call bracket(ax,xx,bx,fa,fx,fb,func)
    endif
    fret=brent_optimize(ax,xx,bx,func,tol_,niter_,xmin)
  end subroutine brent
  !

  !
  function brent_optimize(ax,bx,cx,func,tol,itmax,xmin)
    real(8), intent(in)  :: ax,bx,cx,tol
    real(8), intent(out) :: xmin
    real(8)              :: brent_optimize
    integer              :: itmax
    real(8), parameter   :: cgold=0.3819660d0,zeps=1.0d-3*epsilon(ax)
    integer              :: iter
    real(8)              :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
    interface
       function func(x)
         real(8) :: x
         real(8) :: func
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
          brent_optimize=fx
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
  end function brent_optimize






  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !  Given a function f and its derivative function df, and given a
  !  bracketing triplet of abscissas ax, bx, cx [such that bx is between
  !  ax and cx, and f(bx) is less than both f(ax) and f(cx)], this
  !  routine isolates the minimum to a fractional precision of about
  !  tol using a modification of Brent’s method that uses derivatives.
  !  The abscissa of the minimum is returned as xmin, and the minimum
  !  function value is returned as dbrent, the returned function value.
  !+-------------------------------------------------------------------+ 
  subroutine dbrent_wgrad(func,dfunc,xmin,brack,tol,niter)
    interface
       function func(x)
         real(8) :: x
         real(8) :: func
       end function func
       function dfunc(x)
         real(8) :: x
         real(8) :: dfunc
       end function dfunc
    end interface
    real(8),intent(inout)         :: xmin
    real(8),dimension(:),optional :: brack
    real(8),optional              :: tol
    integer,optional              :: niter
    real(8)                       :: tol_
    integer                       :: niter_
    integer                       :: iter
    real(8)                       :: ax,xx,bx,fa,fx,fb,fret
    !
    tol_=1d-9;if(present(tol))tol_=tol
    Niter_=200;if(present(Niter))Niter_=Niter
    !
    if(present(brack))then
       select case(size(brack))
       case(1)
          stop "Brent error: calling brent with size(brack)==1. None or two points are necessary."
       case(2)
          ax = brack(1)
          xx = brack(2)
          call bracket(ax,xx,bx,fa,fx,fb,func)
       case (3)
          ax = brack(1)
          xx = brack(2)
          bx = brack(3)
       end select
    else
       ax=0d0
       xx=1d0
       call bracket(ax,xx,bx,fa,fx,fb,func)
    endif
    fret=dbrent_optimize(ax,xx,bx,func,dfunc,tol_,niter_,xmin)
  end subroutine dbrent_wgrad
  !
  subroutine dbrent_nograd(func,xmin,brack,tol,niter)
    interface
       function func(x)
         real(8) :: x
         real(8) :: func
       end function func
    end interface
    real(8),intent(inout)         :: xmin
    real(8),dimension(:),optional :: brack
    real(8),optional              :: tol
    integer,optional              :: niter
    real(8)                       :: tol_
    integer                       :: niter_
    integer                       :: iter
    real(8)                       :: ax,xx,bx,fa,fx,fb,fret
    !
    tol_=1d-9;if(present(tol))tol_=tol
    Niter_=200;if(present(Niter))Niter_=Niter
    !
    if(present(brack))then
       select case(size(brack))
       case(1)
          stop "Brent error: calling brent with size(brack)==1. None or two points are necessary."
       case(2)
          ax = brack(1)
          xx = brack(2)
          call bracket(ax,xx,bx,fa,fx,fb,func)
       case (3)
          ax = brack(1)
          xx = brack(2)
          bx = brack(3)
       end select
    else
       ax=0d0
       xx=1d0
       call bracket(ax,xx,bx,fa,fx,fb,func)
    endif
    fret=dbrent_optimize(ax,xx,bx,func,dfunc,tol_,niter_,xmin)
  contains
    function dfunc(x)
      real(8) :: x
      real(8) :: dfunc
      call fgradient_func(func,x,dfunc)
    end function dfunc
    !
    subroutine fgradient_func(funcv,x,fjac,epsfcn)
      implicit none
      interface 
         function funcv(x)
           real(8) :: x
           real(8) :: funcv
         end function funcv
      end interface
      integer          ::  n
      real(8),intent(in) ::  x
      real(8)            ::  x_
      real(8)          ::  fvec
      real(8)          ::  fjac
      real(8),optional ::  epsfcn
      real(8)          ::  eps,eps_
      real(8)          ::  epsmch
      real(8)          ::  h,temp
      real(8)          ::  wa1
      real(8)          ::  wa2
      integer          :: i,j,k
      x_ = x
      eps_= 0.d0; if(present(epsfcn))eps_=epsfcn
      epsmch = epsilon(epsmch)
      eps  = sqrt(max(eps_,epsmch))
      !  Evaluate the function
      fvec = funcv(x_)
      temp = x_
      h    = eps*abs(temp)
      if(h==0d0) h = eps
      x_   = temp + h
      wa1  = funcv(x_)
      x_   = temp
      fjac = (wa1 - fvec)/h
    end subroutine fgradient_func
  end subroutine dbrent_nograd



  function dbrent_optimize(ax,bx,cx,func,fjac,tol,itmax,xmin) result(dbrent)
    real(8),intent(in)  :: ax,bx,cx,tol
    real(8),intent(out) :: xmin
    real(8)             :: dbrent
    integer             :: itmax
    real(8), parameter  :: zeps=1.d-3*epsilon(ax)
    integer             :: iter
    real(8)             :: a,b,d,d1,d2,du,dv,dw,dx,e,fu,fv,fw,fx,olde,tol1,tol2
    real(8)             :: u,u1,u2,v,w,x,xm
    logical             :: ok1,ok2
    interface
       function func(x)
         real(8) :: x
         real(8) :: func
       end function func
       function fjac(x)
         real(8) :: x
         real(8) :: fjac
       end function fjac
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
  end function dbrent_optimize




  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !     Given a function FUNC, and given distinct initial points AX and BX,
  !     this routine searches in the downhill direction (defined by the 
  !     function as evaluated at the initial points) and returns new points
  !     AX, BX, CX which bracket a minimum of the function.  
  !     Also returned are the function values at the three points, 
  !     FA, FB, and FC.
  !+-------------------------------------------------------------------+
  subroutine bracket(ax,bx,cx,fa,fb,fc,func)
    real(8), intent(inout) :: ax,bx
    real(8), intent(out) :: cx,fa,fb,fc
    !...the first parameter is the default ratio by which successive intervals
    !   are magnified; the second is the maximum magnification allowed for a
    !   parabolic-fit step
    real(8), parameter :: gold=1.618034d0,glimit=100.d0,tiny=1.d-20
    real(8) :: fu,q,r,u,ulim
    interface
       function func(x)
         real(8) :: x
         real(8) :: func
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
       if ((bx-u)*(u-cx) > 0.d0) then
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
       else if ((cx-u)*(u-ulim) > 0.d0) then
          fu=func(u)
          if (fu < fc) then
             bx=cx
             cx=u
             u=cx+GOLD*(cx-bx)
             call shft(fb,fc,fu,func(u))
          end if
       else if ((u-ulim)*(ulim-cx) >= 0.d0) then
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
  end subroutine bracket





  !##################################################################
  !##################################################################
  !                 MULTI-DIMENSIONAL FUNCTIONS
  !##################################################################
  !##################################################################



  !+-------------------------------------------------------------------+
  !LMDIF INTERFACE:
  !solve M nonlinear equations in N unknowns with M>N
  !so f(x)=0 can NOT be solved.
  !This look for a solution x so that the norm
  ! transpose(f(x))*f(x) is minimized.
  !+-------------------------------------------------------------------+
  subroutine leastsq_lmdif_func(func,a,xdata,ydata,tol,info)
    interface
       function func(a,xdata,ydata)
         real(8),dimension(:)           :: a
         real(8),dimension(:)           :: xdata
         real(8),dimension(size(xdata)) :: ydata
         real(8),dimension(size(xdata)) :: func
       end function func
    end interface
    real(8),dimension(:)           :: a
    real(8),dimension(:)           :: xdata
    real(8),dimension(size(xdata)) :: ydata
    integer                        :: m
    real(8),optional               :: tol
    integer,optional               :: info
    real(8)                        :: tol_
    integer                        :: info_
    integer                        :: n
    real(8),dimension(size(xdata)) :: fvec
    tol_ = 1.d-15;if(present(tol))tol_=tol
    n=size(a)
    m=size(xdata)
    call lmdif1(leastsq_lmdif1_func2sub,m,n,a,fvec,tol_,info_)
    if(present(info))info=info_
  contains
    subroutine leastsq_lmdif1_func2sub(m,n,a,fvec,iflag)
      integer ::  m
      integer ::  n
      real(8) ::  a(n)
      real(8) ::  fvec(m)
      integer ::  iflag
      fvec = func(a,xdata,ydata)
      if(iflag<0)stop "LEASTSQ_lmdif1_func2sub ERROR: iflag < 0 "
    end subroutine leastsq_lmdif1_func2sub
  end subroutine leastsq_lmdif_func

  subroutine leastsq_lmdif_sub(func,a,xdata,ydata,tol,info)
    interface
       subroutine func(a,xdata,ydata,f)
         real(8),dimension(:)           :: a
         real(8),dimension(:)           :: xdata
         real(8),dimension(size(xdata)) :: ydata
         real(8),dimension(size(xdata)) :: f
       end subroutine func
    end interface
    !
    real(8),dimension(:)           :: a
    real(8),dimension(:)           :: xdata
    real(8),dimension(size(xdata)) :: ydata
    integer                        :: m
    real(8),optional               :: tol
    integer,optional               :: info
    real(8)                        :: tol_
    integer                        :: info_
    integer                        :: n
    real(8),dimension(size(xdata)) :: fvec
    tol_ = 1.d-15;if(present(tol))tol_=tol
    n=size(a)
    m=size(xdata)
    call lmdif1(leastsq_lmdif1_sub2sub,m,n,a,fvec,tol_,info_)
    if(present(info))info=info_
  contains
    subroutine leastsq_lmdif1_sub2sub(m,n,a,fvec,iflag)
      integer ::  m
      integer ::  n
      real(8) ::  a(n)
      real(8) ::  fvec(m)
      integer ::  iflag
      call func(a,xdata,ydata,fvec)
      if(iflag<0)stop "LEASTSQ_LMDIF1_sub2sub ERROR: iflag < 0 "
    end subroutine leastsq_lmdif1_sub2sub
  end subroutine leastsq_lmdif_sub



  !LMDER INTERFACE:
  subroutine leastsq_lmder_func(func,dfunc,a,xdata,ydata,tol,info)
    interface
       function func(a,xdata,ydata)
         real(8),dimension(:)           :: a
         real(8),dimension(:)           :: xdata
         real(8),dimension(size(xdata)) :: ydata
         real(8),dimension(size(xdata)) :: func
       end function func
       !
       function dfunc(a,xdata,ydata)
         real(8),dimension(:)                   :: a
         real(8),dimension(:)                   :: xdata
         real(8),dimension(size(xdata))         :: ydata
         real(8),dimension(size(xdata),size(a)) :: dfunc
       end function dfunc
    end interface
    real(8),dimension(:)                   :: a
    real(8),dimension(:)                   :: xdata
    real(8),dimension(size(xdata))         :: ydata
    integer                                :: m
    real(8),optional                       :: tol
    integer,optional                       :: info
    real(8)                                :: tol_
    integer                                :: info_
    integer                                :: n
    real(8),dimension(size(xdata))         :: fvec
    real(8),dimension(size(xdata),size(a)) :: fjac
    tol_ = 1.d-15;if(present(tol))tol_=tol
    n=size(a)
    m=size(xdata)
    call lmder1(leastsq_lmder1_func2sub,m,n,a,fvec,fjac,m,tol_,info_)
    if(present(info))info=info_
  contains
    subroutine leastsq_lmder1_func2sub(m,n,a,fvec,fjac,ldfjac,iflag)
      integer ::  m
      integer ::  n
      integer ::  ldfjac
      real(8) ::  a(n)
      real(8) ::  fvec(m)
      real(8) ::  fjac(ldfjac,n)
      integer ::  iflag
      if(iflag==1)then
         fvec = func(a,xdata,ydata)
      elseif(iflag==2)then
         fjac = dfunc(a,xdata,ydata)
      endif
      if(iflag<0)stop "LEASTSQ_LMDER1_func2sub ERROR: iflag < 0 "
    end subroutine leastsq_lmder1_func2sub
  end subroutine leastsq_lmder_func

  subroutine leastsq_lmder_sub(func,dfunc,a,xdata,ydata,tol,info)
    interface
       subroutine func(a,xdata,ydata,f)
         real(8),dimension(:)           :: a
         real(8),dimension(:)           :: xdata
         real(8),dimension(size(xdata)) :: ydata
         real(8),dimension(size(xdata)) :: f
       end subroutine func
       !
       subroutine dfunc(a,xdata,ydata,df)
         real(8),dimension(:)                   :: a
         real(8),dimension(:)                   :: xdata
         real(8),dimension(size(xdata))         :: ydata
         real(8),dimension(size(xdata),size(a)) :: df
       end subroutine dfunc
    end interface
    real(8),dimension(:)                   :: a
    real(8),dimension(:)                   :: xdata
    real(8),dimension(size(xdata))         :: ydata
    integer                                :: m
    real(8),optional                       :: tol
    integer,optional                       :: info
    real(8)                                :: tol_
    integer                                :: info_
    integer                                :: n
    real(8),dimension(size(xdata))         :: fvec
    real(8),dimension(size(xdata),size(a)) :: fjac
    tol_ = 1.d-15;if(present(tol))tol_=tol
    n=size(a)
    m=size(xdata)
    call lmder1(leastsq_lmder1_sub2sub,m,n,a,fvec,fjac,m,tol_,info_)
    if(present(info))info=info_
  contains
    subroutine leastsq_lmder1_sub2sub(m,n,a,fvec,fjac,ldfjac,iflag)
      integer ::  m
      integer ::  n
      integer ::  ldfjac
      real(8) ::  a(n)
      real(8) ::  fvec(m)
      real(8) ::  fjac(ldfjac,n)
      integer ::  iflag
      if(iflag==1)then
         call func(a,xdata,ydata,fvec)
      elseif(iflag==2)then
         call dfunc(a,xdata,ydata,fjac)
      endif
      if(iflag<0)stop "LEASTSQ_LMDER1_sub2sub ERROR: iflag < 0 "
    end subroutine leastsq_lmder1_sub2sub
  end subroutine leastsq_lmder_sub










  !+-------------------------------------------------------------------+
  !PURPOSE  : Use non-linear least squares to fit a function, f, to data.
  !+-------------------------------------------------------------------+
  !LMDIF INTERFACE
  subroutine curvefit_lmdif_func(model_func,a,xdata,ydata,tol,info)
    interface
       function model_func(x,a)
         real(8),dimension(:)       :: x
         real(8),dimension(:)       :: a
         real(8),dimension(size(x)) :: model_func
       end function model_func
    end interface
    real(8),dimension(:)           :: a
    real(8),dimension(:)           :: xdata
    real(8),dimension(size(xdata)) :: ydata
    integer                        :: m
    real(8),optional               :: tol
    integer,optional               :: info
    real(8)                        :: tol_
    integer                        :: info_
    integer                        :: n
    real(8),dimension(size(xdata)) :: fvec
    !
    tol_ = 1.d-15;if(present(tol))tol_=tol
    !
    n=size(a)
    m=size(xdata)
    !
    call lmdif1(curvefit_lmdif_func2sub,m,n,a,fvec,tol_,info_)
    !
    if(present(info))info=info_
  contains
    subroutine curvefit_lmdif_func2sub(m,n,a,fvec,iflag)
      integer ::  m
      integer ::  n
      real(8) ::  a(n)
      real(8) ::  fvec(m)
      integer ::  iflag
      fvec = model_func(xdata,a) - ydata
      if(iflag<0)stop "CURVEFIT_LMDIF_func2sub ERROR: iflag < 0 "
    end subroutine curvefit_lmdif_func2sub
  end subroutine curvefit_lmdif_func

  subroutine curvefit_lmdif_sub(model_func,a,xdata,ydata,tol,info)
    interface
       subroutine model_func(x,a,f)
         real(8),dimension(:)       :: x
         real(8),dimension(:)       :: a
         real(8),dimension(size(x)) :: f
       end subroutine model_func
    end interface
    real(8),dimension(:)           :: a
    real(8),dimension(:)           :: xdata
    real(8),dimension(size(xdata)) :: ydata
    integer                        :: m
    real(8),optional               :: tol
    integer,optional               :: info
    real(8)                        :: tol_
    integer                        :: info_
    integer                        :: n
    real(8),dimension(size(xdata)) :: fvec
    !
    tol_ = 1.d-15;if(present(tol))tol_=tol
    !
    n=size(a)
    m=size(xdata)
    !
    call lmdif1(curvefit_lmdif_sub2sub,m,n,a,fvec,tol_,info_)
    !
    if(present(info))info=info_
  contains
    subroutine curvefit_lmdif_sub2sub(m,n,a,fvec,iflag)
      integer ::  m
      integer ::  n
      real(8) ::  a(n)
      real(8) ::  fvec(m),fvec_(m)
      integer ::  iflag
      call model_func(xdata,a,fvec_)
      fvec = fvec_ - ydata
      if(iflag<0)stop "CURVEFIT_LMDIF_sub2sub ERROR: iflag < 0 "
    end subroutine curvefit_lmdif_sub2sub
  end subroutine curvefit_lmdif_sub


  !LMDER INTERFACE:
  subroutine curvefit_lmder_func(model_func,model_dfunc,a,xdata,ydata,tol,info)
    interface
       function model_func(x,a)
         real(8),dimension(:)       :: x
         real(8),dimension(:)       :: a
         real(8),dimension(size(x)) :: model_func
       end function model_func
       !
       function model_dfunc(x,a)
         real(8),dimension(:)               :: x
         real(8),dimension(:)               :: a
         real(8),dimension(size(x),size(a)) :: model_dfunc
       end function model_dfunc
    end interface
    real(8),dimension(:)                   :: a
    real(8),dimension(:)                   :: xdata
    real(8),dimension(size(xdata))         :: ydata
    integer                                :: m
    real(8),optional                       :: tol
    integer,optional                       :: info
    real(8)                                :: tol_
    integer                                :: info_
    integer                                :: n
    real(8),dimension(size(xdata))         :: fvec
    real(8),dimension(size(xdata),size(a)) :: fjac
    !
    tol_ = 1.d-15;if(present(tol))tol_=tol
    n=size(a)
    m=size(xdata)
    !
    call lmder1(curvefit_lmder1_func2sub,m,n,a,fvec,fjac,m,tol_,info_)
    !
    if(present(info))info=info_
  contains
    subroutine curvefit_lmder1_func2sub(m,n,a,fvec,fjac,ldfjac,iflag)
      integer ::  m
      integer ::  n
      integer ::  ldfjac
      real(8) ::  a(n)
      real(8) ::  fvec(m)
      real(8) ::  fjac(ldfjac,n)
      integer ::  iflag
      if(iflag==1)then
         fvec = model_func(xdata,a) - ydata
      elseif(iflag==2)then
         fjac = model_dfunc(xdata,a)
      endif
      if(iflag<0)stop "CURVEFIT_LMDER1_func2sub ERROR: iflag < 0 "
    end subroutine curvefit_lmder1_func2sub
  end subroutine curvefit_lmder_func

  subroutine curvefit_lmder_sub(model_func,model_dfunc,a,xdata,ydata,tol,info)
    interface
       subroutine model_func(x,a,f)
         real(8),dimension(:)       :: x
         real(8),dimension(:)       :: a
         real(8),dimension(size(x)) :: f
       end subroutine model_func
       !
       subroutine model_dfunc(x,a,df)
         real(8),dimension(:)               :: x
         real(8),dimension(:)               :: a
         real(8),dimension(size(x),size(a)) :: df
       end subroutine model_dfunc
    end interface
    real(8),dimension(:)                   :: a
    real(8),dimension(:)                   :: xdata
    real(8),dimension(size(xdata))         :: ydata
    integer                                :: m
    real(8),optional                       :: tol
    integer,optional                       :: info
    real(8)                                :: tol_
    integer                                :: info_
    integer                                :: n
    real(8),dimension(size(xdata))         :: fvec
    real(8),dimension(size(xdata),size(a)) :: fjac
    tol_ = 1.d-15;if(present(tol))tol_=tol
    n=size(a)
    m=size(xdata)
    call lmder1(curvefit_lmder1_sub2sub,m,n,a,fvec,fjac,m,tol_,info_)
    if(present(info))info=info_
  contains
    subroutine curvefit_lmder1_sub2sub(m,n,a,fvec,fjac,ldfjac,iflag)
      integer ::  m
      integer ::  n
      integer ::  ldfjac
      real(8) ::  a(n)
      real(8) ::  fvec(m),fvec_(m)
      real(8) ::  fjac(ldfjac,n)
      integer ::  iflag
      if(iflag==1)then
         call model_func(xdata,a,fvec_)
         fvec = fvec_ - ydata
      elseif(iflag==2)then
         call model_dfunc(xdata,a,fjac)
      endif
      if(iflag<0)stop "CURVEFIT_LMDER1_sub2sub ERROR: iflag < 0 "
    end subroutine curvefit_lmder1_sub2sub
  end subroutine curvefit_lmder_sub




  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  ! NELMIN minimizes a function using the Nelder-Mead algorithm.
  !    This routine seeks the minimum value of a user-specified function.
  !    Simplex function minimisation procedure due to Nelder and Mead (1965),
  !    as implemented by O'Neill(1971, Appl.Statist. 20, 338-45), with
  !    subsequent comments by Chambers+Ertel(1974, 23, 250-1), Benyon(1976,
  !    25, 97) and Hill(1978, 27, 380-2)
  !
  !    Input, external FN, the name of the function which evaluates
  !    the function to be minimized.
  !
  !    Input/output, real ( kind = 8 ) START(N).  On input, a starting point
  !    for the iteration.  On output, this data may have been overwritten.
  !
  ! !    Output, real ( kind = 8 ) XMIN(N), the coordinates of the point which
  ! !    is estimated to minimize the function.
  !
  !  OPTIONAL:
  !
  !    Input, real ( kind = 8 ) REQMIN, the terminating limit for the variance
  !    of the function values.  0 < REQMIN is required.
  !
  !    Input, real ( kind = 8 ) STEP(N), determines the size and shape of the
  !    initial simplex.  The relative magnitudes of its elements should reflect
  !    the units of the variables.
  !
  !    Input, integer ( kind = 4 ) KONVGE, the convergence check is carried out
  !    every KONVGE iterations. 0 < KONVGE is required.
  !
  !    Input, integer ( kind = 4 ) KCOUNT, the maximum number of function
  !    evaluations.
  !
  !    Output, integer ( kind = 4 ) ICOUNT, the number of function evaluations
  !    used.
  !
  !    Output, integer ( kind = 4 ) NUMRES, the number of restarts.
  !
  !    Output, integer ( kind = 4 ) IFAULT, error indicator.
  !    0, no errors detected.
  !    1, REQMIN, N, or KONVGE has an illegal value.
  !    2, iteration terminated because KCOUNT was exceeded without convergence.
  !+-------------------------------------------------------------------+
  subroutine fmin(fn,start,&
       lambda,tol,conv_check,max_fun_calls,fun_calls,num_restart,ierr)
    interface
       function fn(x)
         real(8),dimension(:) :: x
         real(8)              :: fn
       end function fn
    end interface
    real(8)            :: start(:)
    real(8),optional   :: lambda(size(start)) !--> step
    real(8),optional   :: tol                 !--> reqmin 
    integer,optional   :: conv_check          !--> konvge
    integer,optional   :: max_fun_calls       !--> kcount
    integer,optional   :: fun_calls           !--> icount
    integer,optional   :: num_restart         !--> numres
    integer,optional   :: ierr                !--> ifault
    !
    real(8)            :: step(size(start))
    real(8)            :: reqmin
    integer            :: konvge   
    integer            :: kcount
    integer            :: icount
    integer            :: numres    
    integer            :: ifault
    !
    real(8)            :: xmin(size(start))
    integer            :: n
    real(8), parameter :: ccoeff = 0.5D+00
    real(8)            :: del
    real(8), parameter :: ecoeff = 2.0D+00
    real(8), parameter :: eps = 0.001D+00
    integer            :: i
    integer            :: ihi
    integer            :: ilo
    integer            :: j
    integer            :: jcount
    integer            :: l
    real(8)            :: p(size(start),size(start)+1)
    real(8)            :: p2star(size(start))
    real(8)            :: pbar(size(start))
    real(8)            :: pstar(size(start))
    real(8),parameter  :: rcoeff = 1.0D+00
    real(8)            :: rq
    real(8)            :: x
    real(8)            :: y(size(start)+1)
    real(8)            :: y2star
    real(8)            :: ylo
    real(8)            :: ynewlo
    real(8)            :: ystar
    real(8)            :: z
    !
    n = size(start)
    !
    step=1d0;if(present(lambda))step=lambda
    reqmin=1d-8;if(present(tol))reqmin=tol
    konvge=10;if(present(conv_check))konvge=conv_check
    kcount=500;if(present(max_fun_calls))kcount=max_fun_calls    
    !
    !  Check the input parameters.
    !
    if ( reqmin <= 0.0D+00 ) then
       ifault = 1
       return
    end if
    if ( n < 1 ) then
       ifault = 1
       return
    end if
    if ( konvge < 1 ) then
       ifault = 1
       return
    end if
    !
    !  Initialization.
    !
    icount = 0
    numres = 0
    jcount = konvge
    del = 1.0D+00
    rq = reqmin * real ( n, kind = 8 )
    !
    !  Initial or restarted loop.
    !
    do
       p(1:n,n+1) = start(1:n)
       y(n+1) = fn ( start )
       icount = icount + 1
       !
       !  Define the initial simplex.
       !
       do j = 1, n
          x = start(j)
          start(j) = start(j) + step(j) * del
          p(1:n,j) = start(1:n)
          y(j) = fn ( start )
          icount = icount + 1
          start(j) = x
       end do
       !
       !  Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
       !  the vertex of the simplex to be replaced.
       !
       ilo = minloc ( y(1:n+1), 1 )
       ylo = y(ilo)
       !
       !  Inner loop.
       !
       do while ( icount < kcount )
          !
          !  YNEWLO is, of course, the HIGHEST value???
          !
          ihi = maxloc ( y(1:n+1), 1 )
          ynewlo = y(ihi)
          !
          !  Calculate PBAR, the centroid of the simplex vertices
          !  excepting the vertex with Y value YNEWLO.
          !
          do i = 1, n
             pbar(i) = ( sum ( p(i,1:n+1) ) - p(i,ihi) ) / real ( n, kind = 8 )
          end do
          !
          !  Reflection through the centroid.
          !
          pstar(1:n) = pbar(1:n) + rcoeff * ( pbar(1:n) - p(1:n,ihi) )
          ystar = fn ( pstar )
          icount = icount + 1
          !
          !  Successful reflection, so extension.
          !
          if ( ystar < ylo ) then
             p2star(1:n) = pbar(1:n) + ecoeff * ( pstar(1:n) - pbar(1:n) )
             y2star = fn ( p2star )
             icount = icount + 1
             !
             !  Retain extension or contraction.
             !
             if ( ystar < y2star ) then
                p(1:n,ihi) = pstar(1:n)
                y(ihi) = ystar
             else
                p(1:n,ihi) = p2star(1:n)
                y(ihi) = y2star
             end if
             !
             !  No extension.
             !
          else
             l = 0
             do i = 1, n + 1
                if ( ystar < y(i) ) then
                   l = l + 1
                end if
             end do
             if ( 1 < l ) then
                p(1:n,ihi) = pstar(1:n)
                y(ihi) = ystar
                !
                !  Contraction on the Y(IHI) side of the centroid.
                !
             else if ( l == 0 ) then
                p2star(1:n) = pbar(1:n) + ccoeff * ( p(1:n,ihi) - pbar(1:n) )
                y2star = fn ( p2star )
                icount = icount + 1
                !
                !  Contract the whole simplex.
                !
                if ( y(ihi) < y2star ) then
                   do j = 1, n + 1
                      p(1:n,j) = ( p(1:n,j) + p(1:n,ilo) ) * 0.5D+00
                      xmin(1:n) = p(1:n,j)
                      y(j) = fn ( xmin )
                      icount = icount + 1
                   end do
                   ilo = minloc ( y(1:n+1), 1 )
                   ylo = y(ilo)
                   cycle
                   !
                   !  Retain contraction.
                   !
                else
                   p(1:n,ihi) = p2star(1:n)
                   y(ihi) = y2star
                end if
                !
                !  Contraction on the reflection side of the centroid.
                !
             else if ( l == 1 ) then
                p2star(1:n) = pbar(1:n) + ccoeff * ( pstar(1:n) - pbar(1:n) )
                y2star = fn ( p2star )
                icount = icount + 1
                !
                !  Retain reflection?
                !
                if ( y2star <= ystar ) then
                   p(1:n,ihi) = p2star(1:n)
                   y(ihi) = y2star
                else
                   p(1:n,ihi) = pstar(1:n)
                   y(ihi) = ystar
                end if
             end if
          end if
          !
          !  Check if YLO improved.
          !
          if ( y(ihi) < ylo ) then
             ylo = y(ihi)
             ilo = ihi
          end if
          jcount = jcount - 1
          if ( 0 < jcount ) then
             cycle
          end if
          !
          !  Check to see if minimum reached.
          !
          if ( icount <= kcount ) then
             jcount = konvge
             x = sum ( y(1:n+1) ) / real ( n + 1, kind = 8 )
             z = sum ( ( y(1:n+1) - x )**2 )
             if ( z <= rq ) then
                exit
             end if
          end if
       end do
       !
       !  Factorial tests to check that YNEWLO is a local minimum.
       !
       xmin(1:n) = p(1:n,ilo)
       ynewlo = y(ilo)
       if ( kcount < icount ) then
          ifault = 2
          exit
       end if
       ifault = 0
       do i = 1, n
          del = step(i) * eps
          xmin(i) = xmin(i) + del
          z = fn ( xmin )
          icount = icount + 1
          if ( z < ynewlo ) then
             ifault = 2
             exit
          end if
          xmin(i) = xmin(i) - del - del
          z = fn ( xmin )
          icount = icount + 1
          if ( z < ynewlo ) then
             ifault = 2
             exit
          end if
          xmin(i) = xmin(i) + del
       end do
       if ( ifault == 0 ) then
          exit
       end if
       !
       !  Restart the procedure.
       !
       start(1:n) = xmin(1:n)
       del = eps
       numres = numres + 1
    end do
    if(present(fun_calls))fun_calls=icount
    if(present(num_restart))num_restart=numres
    if(present(ierr))ierr=ifault
    start = xmin
    return
  end subroutine fmin





  !+-------------------------------------------------------------------+
  !  PURPOSE  : Minimize the Chi^2 distance using conjugate gradient
  !     Adapted by FRPRM subroutine from NumRec (10.6)
  !     Given a starting point P that is a vector of length N, 
  !     the Fletcher-Reeves-Polak-Ribiere minimisation is performed 
  !     n a functin FUNC,using its gradient as calculated by a 
  !     routine DFUNC. The convergence tolerance on the function 
  !     value is input as FTOL.  
  !     Returned quantities are: 
  !     - P (the location of the minimum), 
  !     - ITER (the number of iterations that were performed), 
  !     - FRET (the minimum value of the function). 
  !     The routine LINMIN is called to perform line minimisations.
  !     Minimisation routines: DFPMIN, D/LINMIN, MNBRAK, D/BRENT and D/F1DIM
  !     come from Numerical Recipes.
  !  NOTE: this routine makes use of abstract interface to communicate 
  !     with routines contained elsewhere. an easier way would be to include
  !     the routines inside each of the two following fmin_cg routines. 
  !+-------------------------------------------------------------------+
  subroutine fmin_cg_df(p,f,df,iter,fret,ftol,itmax,eps,istop,type,iverbose)
    procedure(cgfit_func)                :: f
    procedure(cgfit_fjac)                :: df
    real(8), dimension(:), intent(inout) :: p
    integer, intent(out)                 :: iter
    real(8), intent(out)                 :: fret
    real(8),optional                     :: ftol,eps
    real(8)                              :: ftol_,eps_
    integer, optional                    :: itmax,type,istop
    integer                              :: itmax_,type_,istop_
    integer                              :: its
    real(8)                              :: dgg,fp,gam,gg,err_
    real(8), dimension(size(p))          :: g,h,xi
    logical,optional                     :: iverbose
    logical                              :: iverbose_
    !
    if(associated(func))nullify(func) ; func=>f
    if(associated(fjac))nullify(fjac) ; fjac=>df
    !
    iverbose_=.false.;if(present(iverbose))iverbose_=iverbose
    ftol_=1.d-5
    if(present(ftol))then
       ftol_=ftol
       if(iverbose_)write(*,"(A,ES9.2)")"CG: ftol updated to:",ftol
    endif
    eps_=1.d-4
    if(present(eps))then
       eps_=eps
       if(iverbose_)write(*,"(A,ES9.2)")"CG: eps updated to:",eps
    endif
    itmax_=500
    if(present(itmax))then
       itmax_=itmax
       if(iverbose_)write(*,"(A,I5)")"CG: itmax updated to:",itmax
    endif
    istop_=0
    if(present(istop))then
       istop_=istop
       if(iverbose_)write(*,"(A,I3)")"CG: istop update to:",istop
    endif
    type_=0
    if(present(type))then
       type_=type
       if(iverbose_)write(*,"(A,I3)")"CG: type update to:",type
    endif
    !
    fp=func(p)
    xi=fjac(p)
    g=-xi
    h=g
    xi=h
    do its=1,itmax_
       iter=its
       call dlinmin(p,xi,fret,ftol_)
       select case(istop_)
       case default
          err_ = abs(fret-fp)/(abs(fret)+abs(fp)+eps_)
       case(1)
          err_ = abs(fret-fp)/(abs(fp)+eps_)
       case(2)
          err_ = abs(fret-fp)
       end select
       if( err_ <= ftol_ )return
       fp = fret
       xi = fjac(p)
       gg=dot_product(g,g)
       select case(type_)
       case default             
          dgg=dot_product(xi+g,xi)  !polak-ribiere
       case (1)
          dgg=dot_product(xi,xi)   !fletcher-reeves.
       end select
       if (gg == 0.d0) return
       gam=dgg/gg
       g=-xi
       h=g+gam*h
       xi=h
    end do
    if(iverbose_)write(*,*)"CG: MaxIter",itmax_," exceeded."
    nullify(func)
    nullify(fjac)
    return
  end subroutine fmin_cg_df
  !
  !
  !NUMERICAL EVALUATION OF THE GRADIENT:
  subroutine fmin_cg_f(p,f,iter,fret,ftol,itmax,eps,istop,type,deps,iverbose)
    procedure(cgfit_func)                :: f
    real(8), dimension(:), intent(inout) :: p
    integer, intent(out)                 :: iter
    real(8), intent(out)                 :: fret
    real(8),optional                     :: ftol,eps,deps
    real(8)                              :: ftol_,eps_,deps_
    integer, optional                    :: itmax,type,istop
    integer                              :: itmax_,type_,istop_
    integer                              :: its
    real(8)                              :: dgg,fp,gam,gg,err_
    real(8), dimension(size(p))          :: g,h,xi
    logical,optional                     :: iverbose
    logical                              :: iverbose_
    !
    !this is just to ensure that routine needing dfunc allocated
    !and properly definted will continue to work.
    if(associated(func))nullify(func) ; func=>f
    if(associated(fjac))nullify(fjac) ; fjac=>df
    !
    iverbose_=.false.;if(present(iverbose))iverbose_=iverbose
    ftol_=1.d-5
    if(present(ftol))then
       ftol_=ftol
       if(iverbose_)write(*,"(A,ES9.2)")"CG: ftol updated to:",ftol
    endif
    eps_=1.d-4
    if(present(eps))then
       eps_=eps
       if(iverbose_)write(*,"(A,ES9.2)")"CG: eps updated to:",eps
    endif
    itmax_=500
    if(present(itmax))then
       itmax_=itmax
       if(iverbose_)write(*,"(A,I5)")"CG: itmax updated to:",itmax
    endif
    istop_=0
    if(present(istop))then
       istop_=istop
       if(iverbose_)write(*,"(A,I3)")"CG: istop update to:",istop
    endif
    type_=0
    if(present(type))then
       type_=type
       if(iverbose_)write(*,"(A,I3)")"CG: type update to:",type
    endif
    deps_=0d0
    if(present(deps))then
       deps_=deps
       if(iverbose_)write(*,"(A,ES9.2)")"CG: deps update to:",deps
    endif
    df_eps = deps_
    !
    fp=func(p)
    xi=fjac(p)
    g=-xi;
    h=g
    xi=h
    do its=1,itmax_
       iter=its
       ! print*,"iter:",its,"--",p
       ! print*,"f(p): >>",fp
       ! print*,"grad: >>",xi
       call dlinmin(p,xi,fret,ftol_)
       select case(istop_)
       case default
          err_ = abs(fret-fp)/(abs(fret)+abs(fp)+eps_)
       case(1)
          err_ = abs(fret-fp)/(abs(fp)+eps_)
       case(2)
          err_ = abs(fret-fp)
       end select
       if( err_ <= ftol_)return
       fp=fret
       xi = fjac(p)
       ! if(ifix_)then
       !    if(dot_product(xi,xi)<1d-15)xi=xi+1d-15
       ! endif
       !
       gg=dot_product(g,g)
       !
       select case(type_)
       case default             
          dgg=dot_product(xi+g,xi)  !polak-ribiere
       case (1)
          dgg=dot_product(xi,xi)   !fletcher-reeves.
       end select
       if (gg == 0d0) return
       gam=dgg/gg
       g=-xi
       h=g+gam*h
       xi=h
    end do
    if(iverbose_)write(*,*)"CG: MaxIter",itmax_," exceeded."
    nullify(func)
    nullify(fjac)
    return
  end subroutine fmin_cg_f
  !
  function df(p) 
    real(8),dimension(:)       :: p
    real(8),dimension(size(p)) :: df
    df=f_jac_1n_func(func,size(p),p)
  end function df






  !+-------------------------------------------------------------------+
  !     PURPOSE  : Minimize the Chi^2 distance using conjugate gradient
  !     Adapted from unkown minimize.f routine.
  !     don't worry it works...
  !+-------------------------------------------------------------------+
  subroutine fmin_cgminimize_func(p,fcn,iter,fret,ftol,itmax,iverbose,mode)
    real(8),dimension(:),intent(inout) :: p
    procedure(cgfit_func)              :: fcn
    integer                            :: iter
    real(8)                            :: fret
    real(8),optional                   :: ftol
    real(8)                            :: ftol_
    integer, optional                  :: itmax,mode
    integer                            :: itmax_,mode_,iprint_
    integer                            :: n
    real(8)                            :: f
    real(8),allocatable,dimension(:)   :: x,g,h,w,xprmt
    real(8)                            :: dfn,deps,hh
    integer                            :: iexit,itn
    logical,optional                   :: iverbose
    logical                            :: iverbose_
    !
    if(associated(func))nullify(func) ; func=>fcn
    !
    iverbose_=.false.;if(present(iverbose))iverbose_=iverbose
    iprint_=0;if(iverbose_)iprint_=1
    !
    ftol_=1.d-5
    if(present(ftol))then
       ftol_=ftol
       if(iverbose_)write(*,"(A,ES9.2)")"CG-mininize: ftol updated to:",ftol
    endif
    itmax_=1000
    if(present(itmax))then
       itmax_=itmax
       if(iverbose_)write(*,"(A,I5)")"CG-minimize: itmax updated to:",itmax
    endif
    mode_ =1
    if(present(mode))then
       mode_=mode_
       if(iverbose_)write(*,"(A,I5)")"CG-minimize: mode updated to:",mode       
    endif
    !
    N=size(p)
    allocate(x(N),g(N),h(N*N),w(100*N),xprmt(N))
    dfn=-0.5d0
    hh = 1.d-5
    iexit=0
    !set initial point
    x=p
    xprmt=abs(p)+1.d-15
    call minimize_(fcn_,n,x,f,g,h,w,&
         dfn,xprmt,hh,ftol_,mode_,itmax_,iprint_,iexit,itn)
    !set output variables
    iter=itn
    fret=f
    p=x
    deallocate(x,g,h,w,xprmt)
  end subroutine fmin_cgminimize_func
  subroutine fcn_(n,x,f)
    integer :: n
    real(8) :: x(n)
    real(8) :: f
    f=func(x)
  end subroutine fcn_
  !
  subroutine fmin_cgminimize_sub(p,fcn,iter,fret,ftol,itmax,iverbose,mode)
    real(8),dimension(:),intent(inout) :: p
    interface 
       subroutine fcn(n,x_,f_)
         integer                       :: n
         real(8),dimension(n)          :: x_
         real(8)                       :: f_
       end subroutine fcn
    end interface
    integer                            :: iter
    real(8)                            :: fret
    real(8),optional                   :: ftol
    real(8)                            :: ftol_
    integer, optional                  :: itmax,mode
    integer                            :: itmax_,mode_,iprint_
    integer                            :: n
    real(8)                            :: f
    real(8),allocatable,dimension(:)   :: x,g,h,w,xprmt
    real(8)                            :: dfn,deps,hh
    integer                            :: iexit,itn
    logical,optional                   :: iverbose
    logical                            :: iverbose_
    !
    iverbose_=.false.;if(present(iverbose))iverbose_=iverbose
    iprint_=0;if(iverbose_)iprint_=1
    !
    ftol_=1.d-5
    if(present(ftol))then
       ftol_=ftol
       if(iverbose_)write(*,"(A,ES9.2)")"CG-mininize: ftol updated to:",ftol
    endif
    itmax_=1000
    if(present(itmax))then
       itmax_=itmax
       if(iverbose_)write(*,"(A,I5)")"CG-minimize: itmax updated to:",itmax
    endif
    mode_ =1
    if(present(mode))then
       mode_=mode_
       if(iverbose_)write(*,"(A,I5)")"CG-minimize: mode updated to:",mode       
    endif
    n=size(p)
    allocate(x(n),g(n),h(n*n),w(100*n),xprmt(n))
    dfn=-0.5d0
    hh = 1.d-5
    iexit=0
    !set initial point
    x=p
    xprmt=abs(p)+1.d-15
    call minimize_(fcn,n,x,f,g,h,w,&
         dfn,xprmt,hh,ftol_,mode_,itmax_,iprint_,iexit,itn)
    !set output variables
    iter=itn
    fret=f
    p=x
    deallocate(x,g,h,w,xprmt)
  end subroutine fmin_cgminimize_sub







  !--------------------------------------------------------------------
  ! Conjugate Gradient methods for solving unconstrained nonlinear
  !  optimization problems, as described in the paper:
  !
  ! Gilbert, J.C. and Nocedal, J. (1992). "Global Convergence Properties 
  ! of Conjugate Gradient Methods", SIAM Journal on Optimization, Vol. 2,
  ! pp. 21-42. 
  !--------------------------------------------------------------------
  subroutine fmin_cgplus_df(p,func,fjac,iter,fret,ftol,itmax,imethod,iverb1,iverb2)
    real(8),dimension(:),intent(inout) :: p
    integer                            :: N,i
    interface 
       function func(a)
         real(8),dimension(:)          ::  a
         real(8)                       ::  func
       end function func
       function fjac(a)
         real(8),dimension(:)          :: a
         real(8),dimension(size(a))    :: fjac
       end function fjac
    end interface
    integer,intent(out)                :: iter
    real(8)                            :: fret
    real(8),optional                   :: ftol
    real(8)                            :: ftol_
    integer, optional                  :: itmax,iverb1,iverb2,imethod
    integer                            :: itmax_
    real(8),allocatable,dimension(:)   :: x,g,d,gold,w
    real(8)                            :: f,eps,tlev
    logical                            :: finish
    integer                            :: iprint(2),iflag,method
    integer                            :: nfun,irest
    iprint(1)= -1;if(present(iverb1))iprint(1)=iverb1
    iprint(2)= 0;if(present(iverb2))iprint(2)=iverb2
    method   = 2;if(present(imethod))method=imethod
    ftol_=1.d-5
    if(present(ftol))then
       ftol_=ftol
       if(iprint(1)>=0)write(*,"(A,ES9.2)")"CG+: ftol updated to:",ftol
    endif
    itmax_=1000
    if(present(itmax))then
       itmax_=itmax
       if(iprint(1)>=0)write(*,"(A,I5)")"CG+: itmax updated to:",itmax
    endif
    n     = size(p)
    finish= .false. 
    irest = 1
    allocate(x(n),g(n),d(n),gold(n),w(n))
    x     = p
    iflag = 0
    fgloop: do
       !calculate the function and gradient values here
       f = func(x)
       g = fjac(x)
       cgloop: do
          !call the CG code
          !iflag= 0 : successful termination
          !       1 : return to evaluate f and g
          !       2 : return with a new iterate, try termination test
          !      -i : error
          call cgfam(n,x,f,g,d,gold,iprint,ftol_,w,iflag,irest,method,finish,iter,nfun)
          if(iflag <= 0 .OR. iter > itmax_) exit fgloop
          if(iflag == 1) cycle fgloop
          if(iflag == 2) then
             ! termination test.  the user may replace it by some other test. however,
             ! the parameter 'finish' must be set to 'true' when the test is satisfied.
             tlev= ftol_*(1.d0 + dabs(f))
             i=0
             iloop: do 
                i=i+1
                if(i > n) then
                   finish = .true. 
                   cycle cgloop
                endif
                if(dabs(g(i)) > tlev) then
                   cycle cgloop
                else
                   cycle iloop
                endif
             enddo iloop
          endif
       enddo cgloop
    enddo fgloop
    p=x
    fret=f
    if(iflag<0)stop "CG+ error: iflag < 0"
    if(iprint(1)>=0.AND.iter>=itmax_)write(0,*)"CG+ exit with iter >= itmax"
  end subroutine fmin_cgplus_df

  subroutine fmin_cgplus_f(p,fcn,iter,fret,ftol,itmax,imethod,deps,iverb1,iverb2)
    real(8),dimension(:),intent(inout) :: p
    integer                            :: N,i
    ! interface 
    !    function fcn(a)
    !      real(8),dimension(:)          ::  a
    !      real(8)                       ::  fcn
    !    end function fcn
    ! end interface
    procedure(cgfit_func)              :: fcn
    integer,intent(out)                :: iter
    real(8)                            :: fret
    real(8),optional                   :: ftol,deps
    real(8)                            :: ftol_,deps_
    integer, optional                  :: itmax,iverb1,iverb2,imethod
    integer                            :: itmax_
    real(8),allocatable,dimension(:)   :: x,g,d,gold,w
    real(8)                            :: f,eps,tlev
    logical                            :: finish
    integer                            :: iprint(2),iflag,method
    integer                            :: nfun,irest
    if(associated(func))nullify(func) ; func=>fcn
    if(associated(fjac))nullify(fjac) ; fjac=>dfcn
    iprint(1)= -1;if(present(iverb1))iprint(1)=iverb1
    iprint(2)= 0;if(present(iverb2))iprint(2)=iverb2
    method   = 2;if(present(imethod))method=imethod
    ftol_=1.d-5
    if(present(ftol))then
       ftol_=ftol
       if(iprint(1)>=0)write(*,"(A,ES9.2)")"CG+: ftol updated to:",ftol
    endif
    itmax_=1000
    if(present(itmax))then
       itmax_=itmax
       if(iprint(1)>=0)write(*,"(A,I5)")"CG+: itmax updated to:",itmax
    endif
    deps_=0d0
    if(present(deps))then
       deps_=deps
       if(iprint(1)>=0)write(*,"(A,ES9.2)")"CG: deps update to:",deps
    endif
    df_eps = deps_
    !
    n     = size(p)
    finish= .false. 
    irest = 1
    allocate(x(n),g(n),d(n),gold(n),w(n))
    x     = p
    iflag = 0
    fgloop: do
       !calculate the function and gradient values here
       f = func(x)
       g = fjac(x)
       cgloop: do
          !call the CG code
          !iflag= 0 : successful termination
          !       1 : return to evaluate f and g
          !       2 : return with a new iterate, try termination test
          !      -i : error
          call cgfam(n,x,f,g,d,gold,iprint,ftol_,w,iflag,irest,method,finish,iter,nfun)
          if(iflag <= 0 .OR. iter > itmax_) exit fgloop
          if(iflag == 1) cycle fgloop
          if(iflag == 2) then
             ! termination test.  the user may replace it by some other test. however,
             ! the parameter 'finish' must be set to 'true' when the test is satisfied.
             tlev= ftol_*(1.d0 + dabs(f))
             i=0
             iloop: do 
                i=i+1
                if(i > n) then
                   finish = .true. 
                   cycle cgloop
                endif
                if(dabs(g(i)) > tlev) then
                   cycle cgloop
                else
                   cycle iloop
                endif
             enddo iloop
          endif
       enddo cgloop
    enddo fgloop
    p=x
    fret=f
    if(iflag<0)stop "CG+ error: iflag < 0"
    if(iprint(1)>=0.AND.iter>=itmax_)write(0,*)"CG+ exit with iter >= itmax"
    nullify(func)
    nullify(fjac)
    return
  end subroutine fmin_cgplus_f
  !
  function dfcn(p) 
    real(8),dimension(:)       :: p
    real(8),dimension(size(p)) :: dfcn
    dfcn=f_jac_1n_func(func,size(p),p)
  end function dfcn


  !--------------------------------------------------------------------
  ! Constrained BFGS (L-BFGS_B)
  !  optimization problems, as described in the paper:
  !
  ! Ciyou Zhu , Richard H. Byrd , Peihuang Lu and Jorge Nocedal: "L-BFGS-B: 
  ! FORTRAN SUBROUTINES FOR LARGE-SCALE BOUND CONSTRAINED OPTIMIZATION"
  !--------------------------------------------------------------------



  subroutine bfgs_with_grad(func,grad,x,l_,u_,nbd_)
    interface
       function func(x)
         real(8),dimension(:) :: x
         real(8) :: func
       end function func
    end interface
    interface
       function grad(x)
         real(8),dimension(:) :: x
         real(8),dimension(size(x)) :: grad
       end function grad
    end interface
    integer                                   :: n,m = 5, iprint = -1
    real(8), parameter                        :: factr  = 1.0d+7, pgtol  = 1.0d-5
    character(len=60)                         :: task, csave
    logical                                   :: lsave(4)
    integer                                   :: isave(44)
    real(8)                                   :: dsave(29),f
    integer,dimension(:),allocatable          :: iwa,nbd
    integer,dimension(:),allocatable,optional :: nbd_
    real(8),dimension(:),allocatable          :: x,l,u,g,wa
    real(8),dimension(:),allocatable,optional :: l_,u_
!
    n=size(x)
    allocate( g(n))
    allocate ( iwa(3*n) )
    allocate ( wa(2*m*n + 5*n + 11*m*m + 8*m) )
    task = 'START'
    if(present(nbd_))then
      nbd=nbd_
      l=l_
      u=u_
    else
      allocate ( nbd(n), l(n), u(n))
      nbd=0
      l=0.d0
      u=0.d0
    endif
!     The beginning of the loop
 
    do while(task(1:2).eq.'FG'.or.task.eq.'NEW_X'.or.task.eq.'START') 
         
!     This is the call to the L-BFGS-B code.
         
      call setulb ( n, m, x, l, u, nbd, f, g, factr, pgtol, &
                       wa, iwa, task, iprint,&
                       csave, lsave, isave, dsave )
      if (task(1:2) .eq. 'FG') then  
        f=func(x)
        g=grad(x)
      endif
    end do
    deallocate (nbd,l,u,iwa,wa,g)
  end subroutine bfgs_with_grad



  subroutine bfgs_no_grad(func,x,l_,u_,nbd_)
    interface
       function func(x)
         real(8),dimension(:) :: x
         real(8) :: func
       end function func
    end interface
    integer                                   :: n,m = 5, iprint = -1
    real(8), parameter                        :: factr  = 1.0d+7, pgtol  = 1.0d-5
    character(len=60)                         :: task, csave
    logical                                   :: lsave(4)
    integer                                   :: isave(44)
    real(8)                                   :: dsave(29),f
    integer,dimension(:),allocatable          :: iwa,nbd
    integer,dimension(:),allocatable,optional :: nbd_
    real(8),dimension(:),allocatable          :: x,l,u,g,wa
    real(8),dimension(:),allocatable,optional :: l_,u_
!
    n=size(x)
    allocate( g(n))
    allocate ( iwa(3*n) )
    allocate ( wa(2*m*n + 5*n + 11*m*m + 8*m) )
    task = 'START'
    if(present(nbd_))then
      nbd=nbd_
      l=l_
      u=u_
    else
      allocate ( nbd(n), l(n), u(n))
      nbd=0
      l=0.d0
      u=0.d0
    endif
!     The beginning of the loop
 
    do while(task(1:2).eq.'FG'.or.task.eq.'NEW_X'.or.task.eq.'START') 
         
!     This is the call to the L-BFGS-B code.
         
      call setulb ( n, m, x, l, u, nbd, f, g, factr, pgtol, &
                       wa, iwa, task, iprint,&
                       csave, lsave, isave, dsave )
      if (task(1:2) .eq. 'FG') then  
        f=func(x)
        g=f_jac_1n_func(func,size(x),x)
      endif
    end do
    deallocate (nbd,l,u,iwa,wa,g)
  end subroutine bfgs_no_grad




  !           AUXILIARY JACOBIAN/GRADIENT CALCULATIONS
  !
  !          1 x N Jacobian (df_i/dx_j for i=1;j=1,...,N)
  !-----------------------------------------------------------------------
  subroutine fdjac_1n_func(funcv,x,fjac,epsfcn)
    implicit none
    interface 
       function funcv(x)
         implicit none
         real(8),dimension(:) :: x
         real(8)              :: funcv
       end function funcv
    end interface
    integer          ::  n
    real(8)          ::  x(:)
    real(8)          ::  fvec
    real(8)          ::  fjac(size(x))
    real(8),optional ::  epsfcn
    real(8)          ::  eps,eps_
    real(8)          ::  epsmch
    real(8)          ::  h,temp
    real(8)          ::  wa1
    real(8)          ::  wa2
    integer          :: i,j,k
    n=size(x)
    eps_= df_eps; if(present(epsfcn))eps_=epsfcn
    epsmch = epsilon(epsmch)
    eps  = sqrt(max(eps_,epsmch))
    !  Evaluate the function
    fvec = funcv(x)
    do j=1,n
       temp = x(j)
       h    = eps*abs(temp)
       if(h==0.d0) h = eps
       x(j) = temp + h
       wa1  = funcv(x)
       x(j) = temp
       fjac(j) = (wa1 - fvec)/h
    enddo
  end subroutine fdjac_1n_func

  subroutine fdjac_1n_sub(funcv,x,fjac,epsfcn)
    implicit none
    interface 
       subroutine funcv(n,x,y)
         implicit none
         integer              :: n
         real(8),dimension(n) :: x
         real(8)              :: y
       end subroutine funcv
    end interface
    integer          ::  n
    real(8)          ::  x(:)
    real(8)          ::  fvec
    real(8)          ::  fjac(size(x))
    real(8),optional ::  epsfcn
    real(8)          ::  eps,eps_
    real(8)          ::  epsmch
    real(8)          ::  h,temp
    real(8)          ::  wa1
    real(8)          ::  wa2
    integer          :: i,j,k
    n=size(x)
    eps_= df_eps; if(present(epsfcn))eps_=epsfcn
    epsmch = epsilon(epsmch)
    eps  = sqrt(max(eps_,epsmch))
    !  Evaluate the function
    call funcv(n,x,fvec)
    !  Computation of dense approximate jacobian.
    do j=1,n
       temp = x(j)
       h    = eps*abs(temp)
       if(h==0.d0) h = eps
       x(j) = temp + h
       call funcv(n,x,wa1)
       x(j) = temp
       fjac(j) = (wa1-fvec)/h
    enddo
    return
  end subroutine fdjac_1n_sub

  function f_jac_1n_func(funcv,n,x) result(df)
    interface
       function funcv(x)
         implicit none
         real(8),dimension(:) :: x
         real(8)              :: funcv
       end function funcv
    end interface
    integer               :: n
    real(8), dimension(n) :: x
    real(8), dimension(n) :: df
    call fdjac_1n_func(funcv,x,df)
  end function f_jac_1n_func

  function f_jac_1n_sub(funcv,n,x) result(df)
    interface
       subroutine funcv(n,x,y)
         implicit none
         integer               :: n
         real(8), dimension(n) :: x
         real(8)               :: y
       end subroutine funcv
    end interface
    integer               :: n
    real(8), dimension(n) :: x
    real(8), dimension(n) :: df
    call fdjac_1n_sub(funcv,x,df)
  end function f_jac_1n_sub







END MODULE OPTIMIZE_MINIMIZE
