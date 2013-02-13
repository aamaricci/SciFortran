  include "broydn_utils.f90"
  include "broydn_routines.f90"
  include "broydn_fminln.f90"
  MODULE BROYDEN
    USE BROYDN_UTILS
    USE BROYDN_ROUTINES
    USE BROYDN_FMINLN
    implicit none
    private 

    public :: broydn

  contains

    SUBROUTINE broydn(x,check,noexit)
      REAL(8), DIMENSION(:), INTENT(INOUT) :: x
      logical,optional     :: noexit
      LOGICAL, INTENT(OUT) :: check
      INTEGER, PARAMETER :: MAXITS=50
      REAL(8), PARAMETER :: EPS=epsilon(x),TOLF=1.0e-5,TOLMIN=1.0e-7,TOLX=EPS,STPMX=100.0
      INTEGER :: i,its,k,n
      REAL(8) :: f,fold,stpmax
      REAL(8), DIMENSION(size(x)), TARGET :: fvec
      REAL(8), DIMENSION(size(x)) :: c,d,fvcold,g,p,s,t,w,xold
      REAL(8), DIMENSION(size(x),size(x)) :: qt,r
      LOGICAL :: restrt,sing,noexit_

      noexit_=.false.                   !AA: exit on error, set to true makes broyden not exiting on error!
      if(present(noexit))noexit_=noexit 

      fmin_fvecp=>fvec
      n=size(x)
      f=fmin(x)
      if (maxval(abs(fvec(:))) < 0.01d0*TOLF) then
         check=.false.
         RETURN
      end if
      stpmax=STPMX*max(vabs(x(:)),real(n,8))
      restrt=.true.
      do its=1,MAXITS
         if (restrt) then
            call fdjac(x,fvec,r)
            call qrdcmp(r,c,d,sing)
            if (sing) then
               print*,'singular Jacobian in broydn'
               if(noexit_)then
                  noexit=.false.
                  return
               else
                  stop
               endif
            endif
            call unit_matrix(qt)
            do k=1,n-1
               if (c(k) /= 0.0) then
                  qt(k:n,:)=qt(k:n,:)-outerprod(r(k:n,k),&
                       matmul(r(k:n,k),qt(k:n,:)))/c(k)
               end if
            end do
            where (lower_triangle(n,n)) r(:,:)=0.0
            call put_diag(d(:),r(:,:))
         else
            s(:)=x(:)-xold(:)
            do i=1,n
               t(i)=dot_product(r(i,i:n),s(i:n))
            end do
            w(:)=fvec(:)-fvcold(:)-matmul(t(:),qt(:,:))
            where (abs(w(:)) < EPS*(abs(fvec(:))+abs(fvcold(:)))) &
                 w(:)=0.0
            if (any(w(:) /= 0.0)) then
               t(:)=matmul(qt(:,:),w(:))
               s(:)=s(:)/dot_product(s,s)
               call qrupdt(r,qt,t,s)
               d(:)=get_diag(r(:,:))
               if (any(d(:) == 0.0))then
                  print*,'r singular in broydn'
                  if(noexit_)then
                     noexit=.false.
                     return
                  else
                     stop
                  endif
               endif
            end if
         end if
         p(:)=-matmul(qt(:,:),fvec(:))
         do i=1,n
            g(i)=-dot_product(r(1:i,i),p(1:i))
         end do
         xold(:)=x(:)
         fvcold(:)=fvec(:)
         fold=f
         call rsolv(r,d,p)
         call lnsrch(xold,fold,g,p,x,f,stpmax,check,fmin)
         if (maxval(abs(fvec(:))) < TOLF) then
            check=.false.
            RETURN
         end if
         if (check) then
            if (restrt .or. maxval(abs(g(:))*max(abs(x(:)), &
                 1.0d0)/max(f,0.5d0*n)) < TOLMIN) RETURN
            restrt=.true.
         else
            restrt=.false.
            if (maxval((abs(x(:)-xold(:)))/max(abs(x(:)), &
                 1.0d0)) < TOLX) RETURN
         end if
      end do
      print*,'MAXITS exceeded in broydn'
      if(noexit_)then
         noexit=.false.
         return
      else
         stop
      endif
    END SUBROUTINE broydn
  END MODULE BROYDEN
