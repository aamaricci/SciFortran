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
