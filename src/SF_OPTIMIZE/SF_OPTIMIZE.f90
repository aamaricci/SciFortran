MODULE SF_OPTIMIZE
  USE OPTIMIZE_ROOT_FINDING
  USE OPTIMIZE_MINIMIZE
  !
  USE SF_CONSTANTS
  USE SF_LINALG, only: inv_sym
  private

  !OPTIMIZATION:
  !
  !General purpose
  public :: fmin                !Nelder-Mead
  public :: fmin_cg             !Conjugate-Gradient 1
  public :: fmin_cgplus         !Conjugate-Gradient 2
  public :: fmin_cgminimize     !Conjugate-Gradient 3 (very old f77)
  !
  public :: leastsq             
  public :: curvefit            
  !
  !Constrained (multivariate)
  public :: fmin_bfgs           !BFGS (constrained and not)
  !Global
  ! 
  !Scalar function minimizers
  public :: brent
  public :: dbrent
  public :: bracket


  !ROOT FINDING:
  !Scalar functions
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


contains


  subroutine d_linear_mix_1(x,Fx,alpha)
    real(8),intent(inout),dimension(:)      :: x
    real(8),intent(in),dimension(size(x))   :: Fx
    real(8),intent(in)                      :: alpha
    x = x + alpha*Fx
  end subroutine d_linear_mix_1

  subroutine d_linear_mix_2(x,Fx,alpha)
    real(8),intent(inout),dimension(:,:)              :: x
    real(8),intent(in),dimension(size(x,1),size(x,2)) :: Fx
    real(8),intent(in)                                :: alpha
    x = x + alpha*Fx
  end subroutine d_linear_mix_2

  subroutine d_linear_mix_3(x,Fx,alpha)
    real(8),intent(inout),dimension(:,:,:)                      :: x
    real(8),intent(in),dimension(size(x,1),size(x,2),size(x,3)) :: Fx
    real(8),intent(in)                                          :: alpha
    x = x + alpha*Fx
  end subroutine d_linear_mix_3

  subroutine d_linear_mix_4(x,Fx,alpha)
    real(8),intent(inout),dimension(:,:,:,:)                              :: x
    real(8),intent(in),dimension(size(x,1),size(x,2),size(x,3),size(x,4)) :: Fx
    real(8),intent(in)                                                    :: alpha
    x = x + alpha*Fx
  end subroutine d_linear_mix_4

  subroutine d_linear_mix_5(x,Fx,alpha)
    real(8),intent(inout),dimension(:,:,:,:,:)                                      :: x
    real(8),intent(in),dimension(size(x,1),size(x,2),size(x,3),size(x,4),size(x,5)) :: Fx
    real(8),intent(in)                                                              :: alpha
    x = x + alpha*Fx
  end subroutine d_linear_mix_5

  subroutine d_linear_mix_6(x,Fx,alpha)
    real(8),intent(inout),dimension(:,:,:,:,:,:)                                               :: x
    real(8),intent(in),dimension(size(x,1),size(x,2),size(x,3),size(x,4),size(x,5),size(x,6)) :: Fx
    real(8),intent(in)                                                                         :: alpha
    x = x + alpha*Fx
  end subroutine d_linear_mix_6

  subroutine d_linear_mix_7(x,Fx,alpha)
    real(8),intent(inout),dimension(:,:,:,:,:,:,:)                                                      :: x
    real(8),intent(in),dimension(size(x,1),size(x,2),size(x,3),size(x,4),size(x,5),size(x,6),size(x,7)) :: Fx
    real(8),intent(in)                                                                                  :: alpha
    x = x + alpha*Fx
  end subroutine d_linear_mix_7


  subroutine c_linear_mix_1(x,Fx,alpha)
    complex(8),intent(inout),dimension(:)    :: x
    complex(8),intent(in),dimension(size(x)) :: Fx
    real(8),intent(in)                       :: alpha
    x = x + alpha*Fx
  end subroutine c_linear_mix_1

  subroutine c_linear_mix_2(x,Fx,alpha)
    complex(8),intent(inout),dimension(:,:)              :: x
    complex(8),intent(in),dimension(size(x,1),size(x,2)) :: Fx
    real(8),intent(in)                                   :: alpha
    x = x + alpha*Fx
  end subroutine c_linear_mix_2

  subroutine c_linear_mix_3(x,Fx,alpha)
    complex(8),intent(inout),dimension(:,:,:)                      :: x
    complex(8),intent(in),dimension(size(x,1),size(x,2),size(x,3)) :: Fx
    real(8),intent(in)                                             :: alpha
    x = x + alpha*Fx
  end subroutine c_linear_mix_3

  subroutine c_linear_mix_4(x,Fx,alpha)
    complex(8),intent(inout),dimension(:,:,:,:)                              :: x
    complex(8),intent(in),dimension(size(x,1),size(x,2),size(x,3),size(x,4)) :: Fx
    real(8),intent(in)                                                       :: alpha
    x = x + alpha*Fx
  end subroutine c_linear_mix_4

  subroutine c_linear_mix_5(x,Fx,alpha)
    complex(8),intent(inout),dimension(:,:,:,:,:)                                      :: x
    complex(8),intent(in),dimension(size(x,1),size(x,2),size(x,3),size(x,4),size(x,5)) :: Fx
    real(8),intent(in)                                                                 :: alpha
    x = x + alpha*Fx
  end subroutine c_linear_mix_5

  subroutine c_linear_mix_6(x,Fx,alpha)
    complex(8),intent(inout),dimension(:,:,:,:,:,:)                                               :: x
    complex(8),intent(in),dimension(size(x,1),size(x,2),size(x,3),size(x,4),size(x,5),size(x,6)) :: Fx
    real(8),intent(in)                                                                            :: alpha
    x = x + alpha*Fx
  end subroutine c_linear_mix_6

  subroutine c_linear_mix_7(x,Fx,alpha)
    complex(8),intent(inout),dimension(:,:,:,:,:,:,:)                                                      :: x
    complex(8),intent(in),dimension(size(x,1),size(x,2),size(x,3),size(x,4),size(x,5),size(x,6),size(x,7)) :: Fx
    real(8),intent(in)                                                                                     :: alpha
    x = x + alpha*Fx
  end subroutine c_linear_mix_7








  subroutine d_adaptive_mix(x,Fx,alpha,iter)
    real(8),intent(inout),dimension(:)      :: x
    real(8),intent(in),dimension(size(x))   :: Fx
    real(8),intent(in)                      :: alpha
    integer,intent(in)                      :: iter
    integer                                 :: N
    real(8),parameter                       :: alpha_max = 1d0 ! mixing parameter
    real(8),allocatable,dimension(:),save   :: Fx_prev
    real(8),allocatable,dimension(:),save   :: Beta
    integer                                 :: j
    !
    N = size(x)
    !
    if(iter==1)then
       if(allocated(Fx_prev))deallocate(Fx_prev)
       allocate(Fx_prev(N))
       Fx_prev = 0d0
       !
       if(allocated(Beta))deallocate(Beta)
       allocate(Beta(N))
       Beta = alpha
    endif
    !
    X = X + beta*Fx             !(Fx-x)
    !
    if (iter > 1) then
       do j=1,N
          if(Fx_prev(j) * Fx(j) > 0) then
             beta(j) = beta(j) + alpha
             if (beta(j) > alpha_max) beta(j) = alpha_max
          else
             beta(j) = alpha
          end if
       end do
    end if
    !
    Fx_prev = Fx
    !
  end subroutine d_adaptive_mix

  subroutine c_adaptive_mix(x,Fx,alpha,iter)
    complex(8),intent(inout),dimension(:)    :: x
    complex(8),intent(in),dimension(size(x)) :: Fx
    real(8),intent(in)                       :: alpha
    integer,intent(in)                       :: iter
    integer                                  :: N
    real(8),parameter                        :: alpha_max = 1d0 ! mixing parameter
    complex(8),allocatable,dimension(:),save :: Fx_prev
    real(8),allocatable,dimension(:),save    :: Beta
    integer                                  :: j
    !
    N = size(x)
    !
    if(iter==1)then
       if(allocated(Fx_prev))deallocate(Fx_prev)
       allocate(Fx_prev(N))
       Fx_prev = 0d0
       !
       if(allocated(Beta))deallocate(Beta)
       allocate(Beta(N))
       Beta = alpha
    endif
    !
    X = X + beta*Fx             !(Fx-x)
    !
    if (iter > 1) then
       do j=1,N
          if(abs(Fx_prev(j) * Fx(j)) > 0) then
             beta(j) = beta(j) + alpha
             if (beta(j) > alpha_max) beta(j) = alpha_max
          else
             beta(j) = alpha
          end if
       end do
    end if
    !
    Fx_prev = Fx
    !
  end subroutine c_adaptive_mix


  !X =v[n]
  !Fx=v[n+1]
  subroutine d_broyden_mix(X,Fx,alpha,M,iter,w0)
    real(8),intent(inout),dimension(:)      :: x
    real(8),intent(in),dimension(size(x))   :: Fx
    real(8),intent(in)                      :: alpha
    integer,intent(in)                      :: M
    integer,intent(in)                      :: iter
    real(8),optional                        :: w0
    real(8)                                 :: wg0
    real(8),dimension(M)                    :: wg
    integer                                 :: N
    real(8)                                 :: Xsave(size(X))
    integer                                 :: iter_used,ipos,inext
    real(8),allocatable,dimension(:,:),save :: Df,Dv
    real(8),dimension(M,M)                  :: Beta
    real(8),dimension(M)                    :: Work
    real(8)                                 :: norm,GammaMix
    integer                                 :: i,j
    N   = size(X)
    wg0 = 0.01d0 ; if(present(w0))wg0=w0
    wg  = min(M,8)*1d0
    if(iter==1)then
       if(allocated(Df))deallocate(Df)
       if(allocated(Dv))deallocate(Dv)
       allocate(Df(M,N))
       allocate(Dv(M,N))
       Df   = 0d0
       Dv   = 0d0
    endif
    !
    Xsave = X
    !
    !linear mixing if M=0
    if(M==0)then
       X=X+alpha*Fx
       return
    endif
    !
    !Get broyden mixing pointer
    iter_used = min(iter-1,M)
    ipos      = iter-1-((iter-2)/M)*M
    !
    !DeltaF^(n) = F^(n+1)-F^(n)/|F^(n+1)-F^(n)| 
    !DeltaV^(n) = V^(n+1)-V^(n)/|F^(n+1)-F^(n)| 
    if(iter==1)then
       X=X+alpha*Fx               !(Fx-X)
       return
    else
       Df(ipos,:) = Fx - Df(ipos,:)
       Dv(ipos,:) = X  - Dv(ipos,:)
       norm = dot_product(Df(ipos,:),Df(ipos,:))
       norm = sqrt(norm)
       Df(ipos,:)=Df(ipos,:)/norm
       Dv(ipos,:)=Dv(ipos,:)/norm
    endif
    !
    !Build Beta
    !beta  = [w(0)*w(0)*delta(i,j) + w(i)*w(j)*a(i,j)]^-1
    beta=0d0
    do i=1,iter_used
       do j=i+1,iter_used
          beta(i,j) = wg(i)*wg(j)*dot_product(Df(j,:),Df(i,:))
       enddo
       beta(i,i) = wg0**2 + wg(i)**2
    enddo
    !
    call inv_sym(beta(1:iter_used,1:iter_used))
    !
    !Work(i) = Df(i)*f
    do i=1,iter_used
       work(i) = dot_product(Df(i,:),Fx)
    enddo
    !
    X=X+alpha*Fx
    !
    !v^{i+1} = v^i + alpha*f - sum_i g(i)*w(i)*u(i)
    !where:
    ! - g(i) = [sum_j work(j)*b(j,i)*w(j)]
    ! - u(i) = (alpha*DeltaF(i) + DeltaV(i))
    do i=1,iter_used
       GammaMix=0d0
       do j=1,iter_used
          GammaMix=GammaMix + beta(j,i)*wg(j)*Work(j)
       enddo
       !
       X = X - wg(i)*GammaMix*(alpha*Df(i,:) + Dv(i,:))
    enddo
    !
    !Store the next entries for DeltaF and DeltaV
    inext      = iter-((iter-1)/M)*M
    Df(inext,:)= Fx
    Dv(inext,:)= Xsave
  end subroutine d_broyden_mix

  subroutine c_broyden_mix(X,Fx,alpha,M,iter,w0)
    complex(8),intent(inout),dimension(:)      :: x
    complex(8),intent(in),dimension(size(x))   :: Fx
    real(8),intent(in)                         :: alpha
    integer,intent(in)                         :: M
    integer,intent(in)                         :: iter
    real(8),optional                           :: w0
    real(8)                                    :: wg0
    real(8),dimension(M)                       :: wg
    integer                                    :: N
    complex(8)                                 :: Xsave(size(X))
    integer                                    :: iter_used,ipos,inext
    complex(8),allocatable,dimension(:,:),save :: Df,Dv
    complex(8),dimension(M,M)                  :: Beta
    complex(8),dimension(M)                    :: Work
    complex(8)                                 :: GammaMix
    real(8)                                    :: norm
    integer                                    :: i,j
    N   = size(X)
    wg0 = 0.01d0 ; if(present(w0))wg0=w0
    wg  = min(M,8)*1d0
    if(iter==1)then
       if(allocated(Df))deallocate(Df)
       if(allocated(Dv))deallocate(Dv)
       allocate(Df(M,N))
       allocate(Dv(M,N))
       Df   = 0d0
       Dv   = 0d0
    endif
    !
    Xsave = X
    !
    !linear mixing if M=0
    if(M==0)then
       X=X+alpha*Fx
       return
    endif
    !
    !Get broyden mixing pointer
    iter_used = min(iter-1,M)
    ipos      = iter-1-((iter-2)/M)*M
    !
    !DeltaF^(n) = F^(n+1)-F^(n)/|F^(n+1)-F^(n)| 
    !DeltaV^(n) = V^(n+1)-V^(n)/|F^(n+1)-F^(n)| 
    if(iter==1)then
       X=X+alpha*Fx               !(Fx-X)
       return
    else
       Df(ipos,:) = Fx - Df(ipos,:)
       Dv(ipos,:) = X  - Dv(ipos,:)
       norm = dot_product(Df(ipos,:),Df(ipos,:))
       norm = sqrt(norm)
       Df(ipos,:)=Df(ipos,:)/norm
       Dv(ipos,:)=Dv(ipos,:)/norm
    endif
    !
    !Build Beta
    !beta  = [w(0)*w(0)*delta(i,j) + w(i)*w(j)*a(i,j)]^-1
    beta=0d0
    do i=1,iter_used
       do j=i+1,iter_used
          beta(i,j) = wg(i)*wg(j)*dot_product(Df(j,:),Df(i,:))
       enddo
       beta(i,i) = wg0**2 + wg(i)**2
    enddo
    !
    call inv_sym(beta(1:iter_used,1:iter_used))
    !
    !Work(i) = Df(i)*f
    do i=1,iter_used
       work(i) = dot_product(Df(i,:),Fx)
    enddo
    !
    X=X+alpha*Fx
    !
    !v^{i+1} = v^i + alpha*f - sum_i g(i)*w(i)*u(i)
    !where:
    ! - g(i) = [sum_j work(j)*b(j,i)*w(j)]
    ! - u(i) = (alpha*DeltaF(i) + DeltaV(i))
    do i=1,iter_used
       GammaMix=0d0
       do j=1,iter_used
          GammaMix=GammaMix + beta(j,i)*wg(j)*Work(j)
       enddo
       !
       X = X - wg(i)*GammaMix*(alpha*Df(i,:) + Dv(i,:))
    enddo
    !
    !Store the next entries for DeltaF and DeltaV
    inext      = iter-((iter-1)/M)*M
    Df(inext,:)= Fx
    Dv(inext,:)= Xsave
  end subroutine c_broyden_mix

END MODULE SF_OPTIMIZE
