MODULE BROYDN_UTILS
  implicit none
  INTEGER, PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8
  INTERFACE get_diag
     MODULE PROCEDURE get_diag_rv, get_diag_dv
  END INTERFACE get_diag
  INTERFACE outerdiff
     MODULE PROCEDURE outerdiff_r,outerdiff_d,outerdiff_i
  END INTERFACE outerdiff
  INTERFACE outerprod
     MODULE PROCEDURE outerprod_r,outerprod_d
  END INTERFACE outerprod
  INTERFACE put_diag
     MODULE PROCEDURE put_diag_rv, put_diag_r
  END INTERFACE put_diag
  INTERFACE assert_eq
     MODULE PROCEDURE assert_eq2,assert_eq3,assert_eq4,assert_eqn
  END INTERFACE assert_eq
contains

  FUNCTION get_diag_rv(mat)
    REAL(4),DIMENSION(:,:),INTENT(IN) :: mat
    REAL(4), DIMENSION(size(mat,1)) :: get_diag_rv
    INTEGER :: j
    j=assert_eq2(size(mat,1),size(mat,2),'get_diag_rv')
    do j=1,size(mat,1)
       get_diag_rv(j)=mat(j,j)
    end do
  END FUNCTION get_diag_rv
  !BL
  FUNCTION get_diag_dv(mat)
    REAL(8),DIMENSION(:,:),INTENT(IN) :: mat
    REAL(8), DIMENSION(size(mat,1)) :: get_diag_dv
    INTEGER :: j
    j=assert_eq2(size(mat,1),size(mat,2),'get_diag_dv')
    do j=1,size(mat,1)
       get_diag_dv(j)=mat(j,j)
    end do
  END FUNCTION get_diag_dv
  !
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
  !BL
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
  !BL
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
  !BL
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


  FUNCTION lower_triangle(j,k,extra)
    INTEGER, INTENT(IN) :: j,k
    INTEGER, OPTIONAL, INTENT(IN) :: extra
    LOGICAL, DIMENSION(j,k) :: lower_triangle
    INTEGER :: n
    n=0;if (present(extra)) n=extra
    lower_triangle=(outerdiff(arth_i(1,1,j),arth_i(1,1,k)) > -n)
  contains
    function arth_i(first,increment,n)
      INTEGER, INTENT(IN) :: first,increment,n
      INTEGER, DIMENSION(n) :: arth_i
      INTEGER :: k,k2,temp
      if (n > 0) arth_i(1)=first
      if (n <= NPAR_ARTH) then
         do k=2,n
            arth_i(k)=arth_i(k-1)+increment
         end do
      else
         do k=2,NPAR2_ARTH
            arth_i(k)=arth_i(k-1)+increment
         end do
         temp=increment*NPAR2_ARTH
         k=NPAR2_ARTH
         do
            if (k >= n) exit
            k2=k+k
            arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
            temp=temp+temp
            k=k2
         end do
      end if
    END FUNCTION arth_i
  END FUNCTION lower_triangle


  FUNCTION outerdiff_r(a,b)
    REAL(4), DIMENSION(:), INTENT(IN) :: a,b
    REAL(4), DIMENSION(size(a),size(b)) :: outerdiff_r
    outerdiff_r = spread(a,dim=2,ncopies=size(b)) - spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerdiff_r
  !BL
  FUNCTION outerdiff_d(a,b)
    REAL(8), DIMENSION(:), INTENT(IN) :: a,b
    REAL(8), DIMENSION(size(a),size(b)) :: outerdiff_d
    outerdiff_d = spread(a,dim=2,ncopies=size(b)) - spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerdiff_d
  !BL
  FUNCTION outerdiff_i(a,b)
    INTEGER, DIMENSION(:), INTENT(IN) :: a,b
    INTEGER, DIMENSION(size(a),size(b)) :: outerdiff_i
    outerdiff_i = spread(a,dim=2,ncopies=size(b)) - spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerdiff_i

  FUNCTION outerprod_r(a,b)
    REAL(4), DIMENSION(:), INTENT(IN) :: a,b
    REAL(4), DIMENSION(size(a),size(b)) :: outerprod_r
    outerprod_r = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerprod_r
  !BL
  FUNCTION outerprod_d(a,b)
    REAL(8), DIMENSION(:), INTENT(IN) :: a,b
    REAL(8), DIMENSION(size(a),size(b)) :: outerprod_d
    outerprod_d = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerprod_d

  SUBROUTINE put_diag_rv(diagv,mat)
    REAL(8), DIMENSION(:), INTENT(IN) :: diagv
    REAL(8), DIMENSION(:,:), INTENT(INOUT) :: mat
    INTEGER :: j,n
    n=assert_eq2(size(diagv),min(size(mat,1),size(mat,2)),'put_diag_rv')
    do j=1,n
       mat(j,j)=diagv(j)
    end do
  END SUBROUTINE put_diag_rv
  !BL
  SUBROUTINE put_diag_r(scal,mat)
    REAL(8), INTENT(IN) :: scal
    REAL(8), DIMENSION(:,:), INTENT(INOUT) :: mat
    INTEGER :: j,n
    n = min(size(mat,1),size(mat,2))
    do j=1,n
       mat(j,j)=scal
    end do
  END SUBROUTINE put_diag_r

  SUBROUTINE unit_matrix(mat)
    REAL(8), DIMENSION(:,:), INTENT(OUT) :: mat
    INTEGER :: i,n
    n=min(size(mat,1),size(mat,2))
    mat(:,:)=0.0
    do i=1,n
       mat(i,i)=1.0
    end do
  END SUBROUTINE unit_matrix

  FUNCTION vabs(v)
    REAL(8), DIMENSION(:), INTENT(IN) :: v
    REAL(8) :: vabs
    vabs=sqrt(dot_product(v,v))
  END FUNCTION vabs

  FUNCTION ifirstloc(mask)
    LOGICAL, DIMENSION(:), INTENT(IN) :: mask
    INTEGER               :: ifirstloc
    INTEGER, DIMENSION(1) :: loc
    loc=maxloc(merge(1,0,mask))
    ifirstloc=loc(1)
    if (.not. mask(ifirstloc)) ifirstloc=size(mask)+1
  END FUNCTION ifirstloc
END MODULE BROYDN_UTILS









MODULE BROYDN_ROUTINES
  USE BROYDN_UTILS
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



