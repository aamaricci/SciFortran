  include "MOD_F1DIM.f90"
  include "MOD_LOCAL_UTILS.f90"
  module CGFIT
    USE CGFIT_F1DIM_MOD
    USE CGFIT_LOCAL_UTILS
    implicit none
    private

    public  :: conjugate

  contains

    !+-------------------------------------------------------------------+
    !     PURPOSE  : Minimize the a function FUNC using conjugate gradient
    !     Adapted by FRPRM subroutine from NumRec (10.6).
    !     Given a starting point P(N) the Fletcher-Reeves/Polak-Ribiere 
    !     minimisation is performed on  FUNC, using its gradient 
    !     stored in DFUNC. The convergence tolerance on the function 
    !     value is input as FTOL.  
    !     Returned quantities are: 
    !     - P (the location of the minimum), 
    !     - ITER (the number of iterations that were performed), 
    !     - FRET (the minimum value of the function). 
    !     The routine LINMIN is called to perform line minimisations.
    !     Minimisation routines: DFPMIN, LINMIN, MNBRAK, BRENT and F1DIM
    !     come from Numerical Recipes.
    !+-------------------------------------------------------------------+
    subroutine conjugate(p,ftol,iter,fret)
      IMPLICIT NONE
      REAL(8), DIMENSION(:), INTENT(INOUT) :: p
      REAL(8), INTENT(IN)                  :: ftol
      INTEGER, INTENT(OUT)                 :: iter
      REAL(8), INTENT(OUT)                 :: fret
      INTEGER, PARAMETER                   :: ITMAX=10000
      REAL(8), PARAMETER                   :: EPS=1.0e-6
      INTEGER                              :: its
      REAL(8)                              :: dgg,fp,gam,gg
      REAL(8), DIMENSION(size(p))          :: g,h,xi
      INTERFACE
         FUNCTION FUNC(X)
           REAL(8),DIMENSION(:),INTENT(IN) :: X
           REAL(8)                         :: FUNC
         END FUNCTION FUNC
         FUNCTION DFUNC(X)
           REAL(8),DIMENSION(:),INTENT(IN) :: X
           REAL(8),DIMENSION(SIZE(X))      :: DFUNC
         END FUNCTION DFUNC
      END INTERFACE
      fp=func(p)
      xi=dfunc(p)
      g=-xi
      h=g
      xi=h
      do its=1,ITMAX
         iter=its
         call linmin(p,xi,fret)
         if (2.0*abs(fret-fp) <= ftol*(abs(fret)+abs(fp)+EPS)) RETURN
         !fp=fret
         fp = func(p) !========MODIFICATION=======
         xi = dfunc(p)        
         gg=dot_product(g,g)
         !dgg=dot_product(xi,xi)   !Fletcher-Reeves.
         dgg=dot_product(xi+g,xi)  !Polak-Ribiere
         if (gg == 0.d0) RETURN
         gam=dgg/gg
         g=-xi
         h=g+gam*h
         xi=h
      end do
      print*, 'Too many iteration in CG'
      return
    END SUBROUTINE conjugate



    !********************************************************************
    !********************************************************************
    !********************************************************************





    !+-------------------------------------------------------------------+
    ! PURPOSE: Given an N dimensional point P and an N dimensional 
    ! direction XI, LINMIN moves and resets P to where the function 
    ! FUNC(P) takes on a minimum along the direction XI from P, and 
    ! replaces XI by the actual vector displacement that P was moved.  
    ! Also returns FRET the value of FUNC at the returned location P.  
    ! This is actually all accomplished by calling the routines MNBRAK 
    ! and BRENT.
    !+-------------------------------------------------------------------+
    SUBROUTINE linmin(p,xi,fret)
      REAL(8), INTENT(OUT)                         :: fret
      REAL(8), DIMENSION(:), TARGET, INTENT(INOUT) :: p,xi
      REAL(8), PARAMETER                           :: TOL=1.0e-4
      REAL(8)                                      :: ax,bx,fa,fb,fx,xmin,xx
      ncom=size(p) ; if(ncom /= size(xi))stop "Error in LinMin"
      pcom=>p
      xicom=>xi
      ax=0.0
      xx=1.0
      call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
      fret=brent(ax,xx,bx,f1dim,TOL,xmin)
      !...construct the vector results to return
      xi=xmin*xi
      p=p+xi
      return
    end subroutine linmin


    !********************************************************************
    !********************************************************************
    !********************************************************************

  end module CGFIT
