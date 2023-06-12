MODULE SF_OPTIMIZE
   USE CGFIT_ROUTINES
   USE BROYDEN_ROUTINES
   !
   USE SF_CONSTANTS
   USE SF_LINALG, only: inv_sym
   private

   !!! PUBLIC API: namespace

   !OPTIMIZATION:
   public   :: brent         !minimize a given a function of one-variable with a possible bracketing interval without using derivative information
   public   :: dbrent        !minimize a given a function of one-variable with a possible bracketing interval  using derivative information
   public   :: bracket       !Bracket the minimum of the function.

   !General purpose
   public   :: fmin                !Minimize a function using the Nelder-Mead downhill simplex algorithm.
   public   :: fmin_cg             !Conjugate-Gradient 1
   public   :: fmin_cgplus         !Conjugate-Gradient 2
   public   :: fmin_cgminimize     !Conjugate-Gradient 3 (very old f77)

   !Constrained (multivariate)
   public   :: fmin_bfgs    !Minimize a function using the BFGS algorithm.
   public   :: leastsq      !Minimize the sum of squares of a set of equations. Wrap MINPACK: lmdif/lmder
   public   :: curvefit     !Use non-linear least squares to fit a function, f, to data.

   !> TODO:
   ! public :: fmin_powell  !Minimize a function using modified Powell’s method. This method
   ! public :: fmin_ncg     !Unconstrained minimization of a function using the Newton-CG method.
   ! public :: anneal       !Minimize a function using simulated annealing.
   ! public :: basinhopping ! Find the global minimum of a function using the basin-hopping algorithm ..

   !ROOT FINDING:
   public :: brentq
   public :: bisect
   public :: newton
   public :: fzero
   !Multidimensional
   !General nonlinear solvers:
   public :: fsolve              !
   public :: broyden1
   !Large-scale nonlinear solvers:

   !Fixed points accelerators:
   public :: linear_mix
   public :: adaptive_mix
   public :: broyden_mix
   ! public :: broyden2 !Find a root of a function, using Broyden’s second Jacobian approximation.
   ! public :: newton_krylov !Find a root of a function, using Krylov approximation for inverse Jacobian.
   ! public :: anderson !Find a root of a function, using (extended) Anderson mixing.

   !!! PUBLIC API: interfaces

   interface fmin_cg
      module subroutine fmin_cg_df(p,f,df,iter,fret,ftol,itmax,istop,iverbose)
         procedure(cgfit_func)                :: f
         procedure(cgfit_fjac)                :: df
         real(8), dimension(:), intent(inout) :: p
         integer, intent(out)                 :: iter
         real(8), intent(out)                 :: fret
         real(8),optional                     :: ftol
         integer, optional                    :: itmax,istop
         real(8), dimension(size(p))          :: g,h,xi,p_prev
         logical,optional                     :: iverbose
         ! IMPLEMENTED IN SUBMODULE FMIN_CG.f90
      end subroutine fmin_cg_df
      !
      module subroutine fmin_cg_f(p,f,iter,fret,ftol,itmax,istop,deps,iverbose)
         procedure(cgfit_func)                :: f
         real(8), dimension(:), intent(inout) :: p
         integer, intent(out)                 :: iter
         real(8), intent(out)                 :: fret
         real(8),optional                     :: ftol,deps
         integer, optional                    :: itmax,istop
         logical,optional                     :: iverbose
         ! IMPLEMENTED IN SUBMODULE FMIN_CG.f90
      end subroutine fmin_cg_f
   end interface fmin_cg

   interface fmin_cgplus
      !
      module subroutine fmin_cgplus_df(p,func,fjac,iter,fret,ftol,itmax,imethod,iverb1,iverb2)
         real(8),dimension(:),intent(inout) :: p
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
         integer, optional                  :: itmax,imethod,iverb1,iverb2
         ! IMPLEMENTED IN SUBMODULE FMIN_CG_CGPLUS.f90
      end subroutine fmin_cgplus_df
      !
      module subroutine fmin_cgplus_f(p,fcn,iter,fret,ftol,itmax,imethod,deps,iverb1,iverb2)
         real(8),dimension(:),intent(inout) :: p
         procedure(cgfit_func)              :: fcn
         integer,intent(out)                :: iter
         real(8)                            :: fret
         real(8),optional                   :: ftol,deps
         integer, optional                  :: itmax,imethod,iverb1,iverb2
         ! IMPLEMENTED IN SUBMODULE FMIN_CG_CGPLUS.f90
      end subroutine fmin_cgplus_f
      !
   end interface fmin_cgplus

   interface fmin_cgminimize
      !
      module subroutine fmin_cgminimize_func(p,fcn,iter,fret,ftol,itmax,iverbose,mode,new_version,hh_par)
         real(8),dimension(:),intent(inout) :: p
         procedure(cgfit_func)              :: fcn
         integer                            :: iter
         real(8)                            :: fret
         real(8),optional                   :: ftol
         integer, optional                  :: itmax
         logical,optional                   :: iverbose
         integer,optional                   :: mode
         logical,optional                   :: new_version
         real(8),optional                   :: hh_par
         ! IMPLEMENTATION IN SUBMODULE FMIN_CG_MINIMIZE.f90
      end subroutine fmin_cgminimize_func
      !
      module subroutine fmin_cgminimize_sub(p,fcn,iter,fret,ftol,itmax,iverbose,mode,new_version,hh_par)
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
         integer, optional                  :: itmax
         logical,optional                   :: iverbose
         integer,optional                   :: mode
         logical,optional                   :: new_version
         real(8),optional                   :: hh_par
         ! IMPLEMENTATION IN SUBMODULE FMIN_CG_MINIMIZE.f90
      end subroutine fmin_cgminimize_sub
      !
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

   interface brent
      module subroutine brent(func,xmin,brack,tol,niter)
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
         ! IMPLEMENTED IN SUBMODULE BRENT.f90
      end subroutine
   end interface brent

   interface bracket
      module subroutine bracket(ax,bx,cx,fa,fb,fc,func)
         real(8), intent(inout) :: ax,bx
         real(8), intent(out) :: cx,fa,fb,fc
         interface
            function func(x)
               real(8) :: x
               real(8) :: func
            end function func
         end interface
         ! IMPLEMENTED IN SUBMODULE BRENT.f90
      end subroutine bracket
   end interface bracket

   interface dbrent
      module subroutine dbrent_wgrad(func,dfunc,xmin,brack,tol,niter)
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
         ! IMPLEMENTED IN SUBMODULE BRENT.f90
      end subroutine
      module subroutine dbrent_nograd(func,xmin,brack,tol,niter)
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
         ! IMPLEMENTED IN SUBMODULE BRENT.f90
      end subroutine
   end interface dbrent

   interface brentq
      module function brentq(func,a,b,tol) result(fzero)
         interface
            function func(x)
               real(8),intent(in) :: x
               real(8)            :: func
            end function func
         end interface
         real(8),intent(in) :: a,b
         real(8),optional   :: tol
         real(8)            :: fzero
         ! IMPLEMENTATION IN SUBMODULE FROOT_SCALAR.f90
      end function brentq
   end interface brentq

   interface bisect
      module subroutine bisect(f,x1,x2,eps,niter,flag)
         interface
            function f(x)
               real(8) :: x
               real(8) :: f
            end function f
         end interface
         real(8)          :: x1, x2
         real(8),optional :: eps
         integer,optional :: Niter
         integer,optional :: flag
         ! IMPLEMENTATION IN SUBMODULE FROOT_SCALAR.f90
      end subroutine bisect
   end interface bisect

   interface newton
      module subroutine newton(f,xinit,eps,niter)
         interface
            function f(x)
               real(8) :: x
               real(8) :: f
            end function f
         end interface
         real(8), intent(inout) :: xinit
         real(8),optional       :: eps
         integer,optional       :: Niter
         ! IMPLEMENTATION IN SUBMODULE FROOT_SCALAR.f90
      end subroutine newton
   end interface newton

   interface fzero
      module subroutine fzero(f,b,c,iflag,rguess,tol_rel,tol_abs)
         interface
            function f(x)
               real(8),intent(in) :: x
               real(8)            :: f
            end function f
         end interface
         real(8) :: b
         real(8) :: c
         integer :: iflag
         real(8),optional :: rguess
         real(8),optional :: tol_rel
         real(8),optional :: tol_abs
         ! IMPLEMENTATION IN SUBMODULE FROOT_SCALAR.f90
      end subroutine fzero
   end interface fzero

   interface broyden1
      module subroutine broyden1(ff,x,check,maxits,tol,tol1,tolmin,stpmx,noexit)
         procedure(broydn_func)               :: ff
         real(8), dimension(:), intent(inout) :: x
         logical, optional                    :: noexit
         logical, optional                    :: check
         integer, optional                    :: maxits
         real(8), optional                    :: tol,tol1,tolmin,stpmx
         ! IMPLEMENTATION IN SUBMODULE BROYDEN1.f90
      end subroutine broyden1
   end interface broyden1

   interface fmin_bfgs
      module procedure :: bfgs_with_grad
      module procedure :: bfgs_no_grad
   end interface fmin_bfgs

   interface linear_mix
      module procedure :: d_linear_mix_1
      module procedure :: d_linear_mix_2
      module procedure :: d_linear_mix_3
      module procedure :: d_linear_mix_4
      module procedure :: d_linear_mix_5
      module procedure :: d_linear_mix_6
      module procedure :: d_linear_mix_7
      module procedure :: c_linear_mix_1
      module procedure :: c_linear_mix_2
      module procedure :: c_linear_mix_3
      module procedure :: c_linear_mix_4
      module procedure :: c_linear_mix_5
      module procedure :: c_linear_mix_6
      module procedure :: c_linear_mix_7
   end interface linear_mix


   interface adaptive_mix
      module procedure :: d_adaptive_mix
      module procedure :: c_adaptive_mix
   end interface adaptive_mix

   interface broyden_mix
      module procedure :: d_broyden_mix
      module procedure :: c_broyden_mix
   end interface broyden_mix


   interface fsolve
      module procedure :: fsolve_hybrd_func
      module procedure :: fsolve_hybrd_sub
      !
      module procedure :: fsolve_hybrj_func
      module procedure :: fsolve_hybrj_sub
   end interface fsolve


   !!! MODULE VARIABLES

   real(8)                         :: df_eps=tiny(1d0)
   ! procedure(hybrd_func),pointer :: hybrd_funcv
   real(8), dimension(:),pointer   :: fmin_fvecp

contains

   !!! LOCAL/INCLUDED IMPLEMENTATIONS

   ! INTERFACES TO MINPACK lmder/lmdif
   include "leastsq.f90"
   include "curvefit.f90"


   ! Minimizes a function using the Nelder-Mead algorithm.
   !    This routine seeks the minimum value of a user-specified function.
   !    Simplex function minimisation procedure due to Nelder and Mead (1965),
   !    as implemented by O'Neill(1971, Appl.Statist. 20, 338-45), with
   !    subsequent comments by Chambers+Ertel(1974, 23, 250-1), Benyon(1976,
   !    25, 97) and Hill(1978, 27, 380-2)
   include "fmin_Nelder_Mead.f90"



   ! Constrained BFGS (L-BFGS_B) optimization problems:
   ! Ciyou Zhu , Richard H. Byrd , Peihuang Lu and Jorge Nocedal: "L-BFGS-B:
   ! FORTRAN SUBROUTINES FOR LARGE-SCALE BOUND CONSTRAINED OPTIMIZATION"
   include "fmin_bfgs.f90"




   ! Mixing and Acceleration:
   include "linear_mix.f90"
   include "adaptive_mix.f90"
   include "broyden_mix.f90"


   ! Interface to MINPACK hybrd/hybrj: FSOLVE
   include "fsolve.f90"








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









END MODULE SF_OPTIMIZE


!!! IMPLEMENTATIONS MOVED INTO SUBMODULES

! Conjugate Gradient methods for solving unconstrained nonlinear
!  optimization problems:
! Gilbert, J.C. and Nocedal, J. (1992). "Global Convergence Properties
! of Conjugate Gradient Methods", SIAM Journal on Optimization, Vol. 2,
! pp. 21-42.
include "fmin_cg_cgplus.f90"

!  Minimize the Chi^2 distance using conjugate gradient
!     Wraps legacy fixed-form minimize_*.f routines.
include "fmin_cg_minimize.f90"

! Minimize the Chi^2 distance using conjugate gradient
!     Adapted by FRPRM subroutine from NumRec (10.6),,
!     the Fletcher-Reeves-Polak-Ribiere minimisation is performed
include "fmin_cg.f90"

! Brent methods, including bracket
include "brent.f90"

! Find root of scalar functions
include "froot_scalar.f90"

! Broyden root finding method:
include "broyden1.f90"
