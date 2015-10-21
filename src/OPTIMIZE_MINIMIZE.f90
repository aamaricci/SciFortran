MODULE OPTIMIZE_MINIMIZE
  USE CGFIT_FUNC_INTERFACE
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

  public :: fmin_cg
  public :: fmin_cgplus
  public :: fmin_cgminimize


contains

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
    logical,optional :: iverbose
    logical           :: iverbose_
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
  subroutine fmin_cg_f(p,f,iter,fret,ftol,itmax,eps,istop,type,iverbose)
    procedure(cgfit_func)                :: f
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
    logical,optional :: iverbose
    logical           :: iverbose_
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
    !
    fp=func(p)
    xi=fjac(p)!f_dgradient(func,size(p),p)
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
       if( err_ <= ftol_)return
       fp=fret
       xi = fjac(p)!f_dgradient(func,size(p),p)        
       gg=dot_product(g,g)
       select case(type_)
       case default             
          dgg=dot_product(xi+g,xi)  !polak-ribiere
       case (1)
          dgg=dot_product(xi,xi)   !fletcher-reeves.
       end select
       if (gg == 0.0) return
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
  subroutine fmin_cgminimize_func(p,fcn,iter,fret,ftol,itmax,iprint,mode)
    real(8),dimension(:),intent(inout) :: p
    procedure(cgfit_func)              :: fcn
    integer                            :: iter
    real(8)                            :: fret
    real(8),optional                   :: ftol
    real(8)                            :: ftol_
    integer, optional                  :: itmax,mode,iprint
    integer                            :: itmax_,mode_,iprint_
    integer                            :: n
    real(8)                            :: f
    real(8),allocatable,dimension(:)   :: x,g,h,w,xprmt
    real(8)                            :: dfn,deps,hh
    integer                            :: iexit,itn
    if(associated(func))nullify(func) ; func=>fcn
    iprint_=0;if(present(iprint))iprint_=iprint
    ftol_=1.d-5
    if(present(ftol))then
       ftol_=ftol
       if(iprint_>0)write(*,"(A,ES9.2)")"CG-mininize: ftol updated to:",ftol
    endif
    itmax_=1000
    if(present(itmax))then
       itmax_=itmax
       if(iprint_>0)write(*,"(A,I5)")"CG-minimize: itmax updated to:",itmax
    endif
    mode_ =1
    if(present(mode))then
       mode_=mode_
       if(iprint_>0)write(*,"(A,I5)")"CG-minimize: mode updated to:",mode       
    endif
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
  subroutine fmin_cgminimize_sub(p,fcn,iter,fret,ftol,itmax,iprint,mode)
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
    integer, optional                  :: itmax,mode,iprint
    integer                            :: itmax_,mode_,iprint_
    integer                            :: n
    real(8)                            :: f
    real(8),allocatable,dimension(:)   :: x,g,h,w,xprmt
    real(8)                            :: dfn,deps,hh
    integer                            :: iexit,itn
    iprint_=0;if(present(iprint))iprint_=iprint
    ftol_=1.d-5
    if(present(ftol))then
       ftol_=ftol
       if(iprint_>0)write(*,"(A,ES9.2)")"CG-mininize: ftol updated to:",ftol
    endif
    itmax_=1000
    if(present(itmax))then
       itmax_=itmax
       if(iprint_>0)write(*,"(A,I5)")"CG-minimize: itmax updated to:",itmax
    endif
    mode_ =1
    if(present(mode))then
       mode_=mode_
       if(iprint_>0)write(*,"(A,I5)")"CG-minimize: mode updated to:",mode       
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

  subroutine fmin_cgplus_f(p,fcn,iter,fret,ftol,itmax,imethod,iverb1,iverb2)
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
    real(8),optional                   :: ftol
    real(8)                            :: ftol_
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
    eps_= 0.d0; if(present(epsfcn))eps_=epsfcn
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
    eps_= 0.d0; if(present(epsfcn))eps_=epsfcn
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
