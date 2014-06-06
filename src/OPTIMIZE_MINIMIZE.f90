  include "optimize_cgfit_routines.f90"
  MODULE MINIMIZE
    USE CGFIT_FUNC_INTERFACE
    USE CGFIT_ROUTINES
    USE DERIVATE, only:f_dgradient
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


    subroutine fmin_cgminimize_func(p,func,iter,fret,ftol,itmax,iprint,mode)
      real(8),dimension(:),intent(inout) :: p
      interface 
         function func(x_)
           real(8),dimension(:)          :: x_
           real(8)                       :: func
         end function func
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
    contains
      subroutine fcn_(n,x,f)
        integer :: n
        real(8) :: x(n)
        real(8) :: f
        f=func(x)
      end subroutine fcn_
    end subroutine fmin_cgminimize_func



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

    subroutine fmin_cgplus_f(p,func,iter,fret,ftol,itmax,imethod,iverb1,iverb2)
      real(8),dimension(:),intent(inout) :: p
      integer                            :: N,i
      interface 
         function func(a)
           real(8),dimension(:)          ::  a
           real(8)                       ::  func
         end function func
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
         g = f_dgradient(func,size(x),x)
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
    end subroutine fmin_cgplus_f









    !+-------------------------------------------------------------------+
    !     PURPOSE  : Minimize the Chi^2 distance using conjugate gradient
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
      if(associated(dfunc))nullify(dfunc) ; dfunc=>df
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
      xi=dfunc(p)
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
         xi = dfunc(p)
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
      nullify(dfunc)
      return
    end subroutine fmin_cg_df

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
      if(associated(func))nullify(func) ; func=>f
      !this is just to ensure that routine needing dfunc allocated
      !and properly definted will continue to work.
      if(associated(dfunc))nullify(dfunc) ; dfunc=>dfunc_
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
      xi=f_dgradient(func,size(p),p)
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
         xi = f_dgradient(func,size(p),p)        
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
      nullify(dfunc)
      return
    end subroutine fmin_cg_f

    function dfunc_(p) 
      real(8),dimension(:)       :: p
      real(8),dimension(size(p)) :: dfunc_
      dfunc_=f_dgradient(func,size(p),p)
    end function dfunc_

  END MODULE MINIMIZE
