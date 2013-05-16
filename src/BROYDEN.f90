  include "broydn_routines.f90"
  MODULE BROYDEN    
    USE BROYDN_ROUTINES
    implicit none
    private 
    public :: broydn
    real(8), dimension(:), pointer :: fmin_fvecp

  contains

    function fmin(x)
      real(8), dimension(:), intent(in) :: x
      real(8)                           :: fmin
      interface
         function funcv(x)
           real(8), dimension(:), intent(in) :: x
           real(8), dimension(size(x))       :: funcv
         end function funcv
      end interface
      if (.not. associated(fmin_fvecp)) then
         print*,'fmin: problem with pointer for returned values'
         stop
      endif
      fmin_fvecp=funcv(x)
      fmin=0.5d0*dot_product(fmin_fvecp,fmin_fvecp)
    end function fmin


    subroutine broydn(x,check,maxits,tolf,tolmin,stpmx,noexit)
      real(8), dimension(:), intent(inout) :: x
      logical, optional                    :: noexit
      logical, intent(out)                 :: check
      integer, optional                    :: maxits
      integer                              :: maxits_=200
      real(8), optional                    :: tolf,tolmin,stpmx
      real(8)                              :: tolf_=1.0e-5,tolmin_=1.0e-7,stpmx_=100.0
      real(8),parameter                    :: eps=epsilon(x),tolx=eps
      integer                              :: i,its,k,n
      real(8)                              :: f,fold,stpmax
      real(8), dimension(size(x)), target  :: fvec
      real(8), dimension(size(x))          :: c,d,fvcold,g,p,s,t,w,xold
      real(8), dimension(size(x),size(x))  :: qt,r
      logical                              :: restrt,sing,noexit_=.false. 
      if(present(noexit))   noexit_  = noexit
      !
      if(present(MAXITS))   MAXITS_  = MAXITS  
      !
      if(present(TOLF))     TOLF_    = TOLF
      !
      if(present(TOLMIN))   TOLMIN_  = TOLMIN  
      !
      if(present(STPMX))    STPMX_   = STPMX  
      !
      fmin_fvecp=>fvec
      !
      n=size(x)
      f=fmin(x)
      if (maxval(abs(fvec(:))) < 0.01d0*TOLF_) then
         check=.false.
         RETURN
      end if
      stpmax=STPMX_*max(vabs(x(:)),real(n,8))
      restrt=.true.
      do its=1,MAXITS_
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
         if (maxval(abs(fvec(:))) < TOLF_) then
            check=.false.
            RETURN
         end if
         if (check) then
            if (restrt .or. maxval(abs(g(:))*max(abs(x(:)), &
                 1.0d0)/max(f,0.5d0*n)) < TOLMIN_) RETURN
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
