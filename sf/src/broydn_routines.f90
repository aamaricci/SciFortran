!###############################################################
! PROGRAM  : ROUTINES
! TYPE     : Module
! PURPOSE  : 
!###############################################################
MODULE BROYDN_ROUTINES
  USE BROYDN_UTILS, ONLY : assert_eq
  implicit none

contains

  SUBROUTINE fdjac(x,fvec,df)
    REAL(8), DIMENSION(:), INTENT(IN)    :: fvec
    REAL(8), DIMENSION(:), INTENT(INOUT) :: x
    REAL(8), DIMENSION(:,:), INTENT(OUT) :: df
    INTERFACE
       FUNCTION funcv(x)
         REAL(8), DIMENSION(:), INTENT(IN) :: x
         REAL(8), DIMENSION(size(x))       :: funcv
       END FUNCTION funcv
    END INTERFACE
    REAL(8), PARAMETER :: EPS=1.0d-4
    INTEGER :: j,n
    REAL(8), DIMENSION(size(x)) :: xsav,xph,h
    n=assert_eq(size(x),size(fvec),size(df,1),size(df,2),'fdjac')
    xsav=x
    h=EPS*abs(xsav)
    where (h == 0.0) h=EPS
    xph=xsav+h
    h=xph-xsav
    do j=1,n
       x(j)=xph(j)
       df(:,j)=(funcv(x)-fvec(:))/h(j)
       x(j)=xsav(j)
    end do
  END SUBROUTINE fdjac


  SUBROUTINE lnsrch(xold,fold,g,p,x,f,stpmax,check,func)
    USE BROYDN_UTILS, ONLY : assert_eq,vabs
    REAL(8), DIMENSION(:), INTENT(IN) :: xold,g
    REAL(8), DIMENSION(:), INTENT(INOUT) :: p
    REAL(8), INTENT(IN) :: fold,stpmax
    REAL(8), DIMENSION(:), INTENT(OUT) :: x
    REAL(8), INTENT(OUT) :: f
    LOGICAL, INTENT(OUT) :: check
    INTERFACE
       FUNCTION func(x)
         REAL(8)                          :: func
         REAL(8), DIMENSION(:), INTENT(IN):: x
       END FUNCTION func
    END INTERFACE
    REAL(8), PARAMETER :: ALF=1.0e-4,TOLX=epsilon(x)
    INTEGER :: ndum
    REAL(8) :: a,alam,alam2,alamin,b,disc,f2,pabs,rhs1,rhs2,slope,tmplam
    ndum=assert_eq(size(g),size(p),size(x),size(xold),'lnsrch')
    check=.false.
    pabs=vabs(p(:))
    if (pabs > stpmax) p(:)=p(:)*stpmax/pabs
    slope=dot_product(g,p)
    if (slope >= 0.0) then
       print*,'roundoff problem in lnsrch'
       stop
    endif
    alamin=TOLX/maxval(abs(p(:))/max(abs(xold(:)),1.0d0))
    alam=1.0
    do
       x(:)=xold(:)+alam*p(:)
       f=func(x)
       if (alam < alamin) then
          x(:)=xold(:)
          check=.true.
          RETURN
       else if (f <= fold+ALF*alam*slope) then
          RETURN
       else
          if (alam == 1.0) then
             tmplam=-slope/(2.0d0*(f-fold-slope))
          else
             rhs1=f-fold-alam*slope
             rhs2=f2-fold-alam2*slope
             a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
             b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/&
                  (alam-alam2)
             if (a == 0.0) then
                tmplam=-slope/(2.0d0*b)
             else
                disc=b*b-3.0d0*a*slope
                if (disc < 0.0) then
                   tmplam=0.5d0*alam
                else if (b <= 0.0) then
                   tmplam=(-b+sqrt(disc))/(3.0d0*a)
                else
                   tmplam=-slope/(b+sqrt(disc))
                end if
             end if
             if (tmplam > 0.5d0*alam) tmplam=0.5d0*alam
          end if
       end if
       alam2=alam
       f2=f
       alam=max(tmplam,0.1d0*alam)
    end do
  END SUBROUTINE lnsrch

  SUBROUTINE qrdcmp(a,c,d,sing)
    USE BROYDN_UTILS, ONLY : assert_eq,outerprod,vabs
    REAL(8), DIMENSION(:,:), INTENT(INOUT) :: a
    REAL(8), DIMENSION(:), INTENT(OUT) :: c,d
    LOGICAL, INTENT(OUT) :: sing
    INTEGER :: k,n
    REAL(8) :: scale,sigma
    n=assert_eq(size(a,1),size(a,2),size(c),size(d),'qrdcmp')
    sing=.false.
    do k=1,n-1
       scale=maxval(abs(a(k:n,k)))
       if (scale == 0.0) then
          sing=.true.
          c(k)=0.0
          d(k)=0.0
       else
          a(k:n,k)=a(k:n,k)/scale
          sigma=sign(vabs(a(k:n,k)),a(k,k))
          a(k,k)=a(k,k)+sigma
          c(k)=sigma*a(k,k)
          d(k)=-scale*sigma
          a(k:n,k+1:n)=a(k:n,k+1:n)-outerprod(a(k:n,k),&
               matmul(a(k:n,k),a(k:n,k+1:n)))/c(k)
       end if
    end do
    d(n)=a(n,n)
    if (d(n) == 0.0) sing=.true.
  END SUBROUTINE qrdcmp


  SUBROUTINE qrupdt(r,qt,u,v)
    USE BROYDN_UTILS, ONLY : assert_eq,ifirstloc
    REAL(8), DIMENSION(:,:), INTENT(INOUT) :: r,qt
    REAL(8), DIMENSION(:), INTENT(INOUT) :: u
    REAL(8), DIMENSION(:), INTENT(IN) :: v
    INTEGER :: i,k,n
    n=assert_eq((/size(r,1),size(r,2),size(qt,1),size(qt,2),size(u),&
         size(v)/),'qrupdt')
    k=n+1-ifirstloc(u(n:1:-1) /= 0.0)
    if (k < 1) k=1
    do i=k-1,1,-1
       call rotate(r,qt,i,u(i),-u(i+1))
       u(i)=pythag(u(i),u(i+1))
    end do
    r(1,:)=r(1,:)+u(1)*v
    do i=1,k-1
       call rotate(r,qt,i,r(i,i),-r(i+1,i))
    end do
  contains
    SUBROUTINE rotate(r,qt,i,a,b)
      USE BROYDN_UTILS, ONLY : assert_eq
      REAL(8), DIMENSION(:,:), TARGET, INTENT(INOUT) :: r,qt
      INTEGER, INTENT(IN) :: i
      REAL(8), INTENT(IN) :: a,b
      REAL(8), DIMENSION(size(r,1)) :: temp
      INTEGER :: n
      REAL(8) :: c,fact,s
      n=assert_eq(size(r,1),size(r,2),size(qt,1),size(qt,2),'rotate')
      if (a == 0.0) then
         c=0.0
         s=sign(1.0d0,b)
      else if (abs(a) > abs(b)) then
         fact=b/a
         c=sign(1.0d0/sqrt(1.0d0+fact**2),a)
         s=fact*c
      else
         fact=a/b
         s=sign(1.0d0/sqrt(1.0d0+fact**2),b)
         c=fact*s
      end if
      temp(i:n)=r(i,i:n)
      r(i,i:n)=c*temp(i:n)-s*r(i+1,i:n)
      r(i+1,i:n)=s*temp(i:n)+c*r(i+1,i:n)
      temp=qt(i,:)
      qt(i,:)=c*temp-s*qt(i+1,:)
      qt(i+1,:)=s*temp+c*qt(i+1,:)
    END SUBROUTINE rotate

    FUNCTION pythag(a,b)
      REAL(8), INTENT(IN) :: a,b
      REAL(8) :: pythag
      REAL(8) :: absa,absb
      absa=abs(a)
      absb=abs(b)
      if (absa > absb) then
         pythag=absa*sqrt(1.0d0+(absb/absa)**2)
      else
         if (absb == 0.0) then
            pythag=0.0
         else
            pythag=absb*sqrt(1.0d0+(absa/absb)**2)
         end if
      end if
    END FUNCTION pythag
  END SUBROUTINE qrupdt


  SUBROUTINE rsolv(a,d,b)
    USE BROYDN_UTILS, ONLY : assert_eq
    REAL(8), DIMENSION(:,:), INTENT(IN) :: a
    REAL(8), DIMENSION(:), INTENT(IN) :: d
    REAL(8), DIMENSION(:), INTENT(INOUT) :: b
    INTEGER :: i,n
    n=assert_eq(size(a,1),size(a,2),size(b),size(d),'rsolv')
    b(n)=b(n)/d(n)
    do i=n-1,1,-1
       b(i)=(b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/d(i)
    end do
  END SUBROUTINE rsolv
END MODULE BROYDN_ROUTINES



