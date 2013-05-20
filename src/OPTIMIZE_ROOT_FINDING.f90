  include "optimize_broyden_routines.f90"
  module HYBRD_INTERFACE
    implicit none
    abstract interface 
       function hybrd_func(x)
         real(8),dimension(:)       :: x
         real(8),dimension(size(x)) :: hybrd_func
       end function hybrd_func
    end interface
  end module HYBRD_INTERFACE

  MODULE ROOT_FINDING
    USE BROYDEN_ROUTINES
    USE HYBRD_INTERFACE
    implicit none
    private 
    interface broydn
       module procedure fzero_broyden
    end interface broydn

    public :: fzero_brentq
    public :: zbrent

    public :: fzero_hybrd
    public :: fzero_broyden
    public :: broydn            !backward compatibility

    public :: fsolve

    procedure(hybrd_func),pointer :: hybrd_funcv
    real(8), dimension(:),pointer :: fmin_fvecp

  contains

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !SCALAR FUNCTIONS:
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    function fzero_brentq(func,a,b) result(fzero)
      interface
         function func(x)
           real(8),intent(in) :: x
           real(8)            :: func
         end function func
      end interface
      real(8),intent(in) :: a,b
      real(8)            :: fzero
      real(8)            :: tol
      tol=epsilon(a)
      fzero = zbrent(func,a,b,tol)
    end function fzero_brentq

    function zbrent(func,x1,x2,tol)
      implicit none
      interface
         function func(x)
           implicit none
           real(8), intent(in) :: x
           real(8)             :: func
         end function func
      end interface
      real(8), intent(in) :: x1,x2,tol
      real(8)             :: zbrent
      integer, parameter  :: itmax=100
      real(8), parameter  :: eps=epsilon(x1)
      integer             :: iter
      real(8)             :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
      a=x1       ; b=x2
      fa=func(a) ; fb=func(b)
      if ((fa > 0.d0 .AND. fb > 0.d0) .OR. (fa < 0.d0 .AND. fb < 0.d0))then
         write(*,"(A)")'ROOT MUST BE BRACKETED FOR ZBRENT'
         stop
      endif
      c=b ; fc=fb
      do iter=1,itmax
         if ((fb > 0.d0 .AND. fc > 0.d0) .OR. (fb < 0.d0 .AND. fc < 0.d0)) then
            c=a
            fc=fa
            d=b-a
            e=d
         end if
         if (abs(fc) < abs(fb)) then
            a=b
            b=c
            c=a
            fa=fb
            fb=fc
            fc=fa
         end if
         tol1=2.d0*eps*abs(b)+0.5d0*tol
         xm=0.5d0*(c-b)
         if (abs(xm) <= tol1 .or. fb == 0.d0) then
            zbrent=b
            return
         end if
         !
         if (abs(e) >= tol1 .AND. abs(fa) > abs(fb)) then
            s=fb/fa
            if (a == c) then
               p=2.d0*xm*s
               q=1.d0-s
            else
               q=fa/fc
               r=fb/fc
               p=s*(2.d0*xm*q*(q-r)-(b-a)*(r-1.d0))
               q=(q-1.d0)*(r-1.d0)*(s-1.d0)
            end if
            if (p > 0.d0) q=-q
            p=abs(p)
            if (2.d0*p  <  min(3.d0*xm*q-abs(tol1*q),abs(e*q))) then
               e=d
               d=p/q
            else
               d=xm
               e=d
            end if
         else
            d=xm
            e=d
         end if
         a=b
         fa=fb
         b=b+merge(d,sign(tol1,xm), abs(d) > tol1 )
         fb=func(b)
      end do
      write(*,"(A)")'zbrent: exceeded maximum iterations'
      zbrent=b
    end function zbrent






    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !MULTI-DIMENSIONAL:
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine fsolve(ff,x,tol,info)
      procedure(hybrd_func)      :: ff
      real(8),dimension(:)       :: x      
      real(8),optional           :: tol
      integer,optional           :: info
      real(8)                    :: tol_
      integer                    :: info_
      integer                    :: n
      real(8),dimension(size(x)) :: fvec
      if(associated(hybrd_funcv))nullify(hybrd_funcv)
      hybrd_funcv=>ff
      tol_ = 1.d-15;if(present(tol))tol_=tol
      n=size(x)
      call hybrd1(func,n,x,fvec,tol_,info_)
      if(present(info))info=info_
    end subroutine fsolve
    !
    subroutine func(n,x,fvec,iflag)
      integer ::  n
      real(8) ::  x(n)
      real(8) ::  fvec(n)
      integer ::  iflag
      fvec(:) = hybrd_funcv(x)
    end subroutine func


    subroutine fzero_hybrd(func,x,tol,info)
      real(8),dimension(:)       :: x
      real(8),dimension(size(x)) :: fvec
      integer                    :: n
      real(8),optional           :: tol
      integer,optional           :: info
      real(8)                    :: tol_
      integer                    :: info_
      external func
      tol_ = 1.d-15;if(present(tol))tol_=tol
      n=size(x)
      call hybrd1(func,n,x,fvec,tol_,info_)
      if(present(info))info=info_
    end subroutine fzero_hybrd



    subroutine fzero_broyden(ff,x,check,maxits,tolf,tolmin,stpmx,noexit)
      procedure(broydn_func)               :: ff
      real(8), dimension(:), intent(inout) :: x
      logical, optional                    :: noexit
      logical, optional                    :: check
      integer, optional                    :: maxits
      real(8), optional                    :: tolf,tolmin,stpmx
      logical                              :: check_
      integer                              :: maxits_=200
      real(8)                              :: tolf_=1.0e-10,tolmin_=1.0e-7,stpmx_=100.0
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
      if(associated(funcv))nullify(funcv)
      funcv=>ff
      fmin_fvecp=>fvec
      !
      n=size(x)
      f=fmin(x)
      if (maxval(abs(fvec(:))) < 0.01d0*TOLF_) then
         check_=.false.
         if(present(check))check=check_
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
         call lnsrch(xold,fold,g,p,x,f,stpmax,check_,fmin)
         if (maxval(abs(fvec(:))) < TOLF_) then
            check_=.false.
            if(present(check))check=check_
            RETURN
         end if
         if (check_) then
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
         if(present(check))check=check_
         return
      else
         if(present(check))check=check_
         stop
      endif
    end subroutine fzero_broyden


    function fmin(x)
      real(8), dimension(:), intent(in) :: x
      real(8)                           :: fmin
      if (.not. associated(fmin_fvecp)) then
         print*,'fmin: problem with pointer for returned values'
         stop
      endif
      fmin_fvecp=funcv(x)
      fmin=0.5d0*dot_product(fmin_fvecp,fmin_fvecp)
    end function fmin


  END MODULE ROOT_FINDING
