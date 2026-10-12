module SF_CHEBYCHEV
  USE SF_CONSTANTS, only: pi, zero
  USE SF_IOTOOLS, only: to_lower, str
  implicit none
  private

  interface cheb_coefficients
    module procedure :: coefficients_d
    module procedure :: coefficients_c
  end interface cheb_coefficients

  interface cheb_set
    module procedure :: set_d
    module procedure :: set_c
  end interface cheb_set

  interface cheb_projection
    module procedure :: projection_rows_d
    module procedure :: projection_rows_c
    module procedure :: projection_vectors_c
    module procedure :: projection_vectors_d
  end interface cheb_projection

  interface cheb_moments
    module procedure :: moments_rows_d
    module procedure :: moments_rows_c
    module procedure :: moments_vectors_d
    module procedure :: moments_vectors_c
  end interface cheb_moments

  public :: cheb_coefficients
  public :: cheb_set
  public :: cheb_projection
  public :: cheb_moments

  !Variables.
  !shift and rescale global Factors:
  real(8), save                 :: cheb_a = 1d0, cheb_b = 0d0
  logical, save                 :: scale_ready = .false.
  logical, save                 :: complex_coeff = .false.
  real(8), allocatable, save    :: cheb_coeff(:)
  complex(8), allocatable, save :: cheb_zcoeff(:)

contains





  !Get Chebyshev coefficients for a given function. A=scale, B=shift.
  subroutine coefficients_d(fun, a, b, coeff, kind, tol, activate)
    interface
      function fun(x) result(y)
        real(8), intent(in) :: x
        real(8)             :: y
      end function fun
    end interface
    logical, optional, intent(in) :: activate
    logical                       :: activate_
    integer, parameter            :: MaxIter = 128
    real(8), intent(in)           :: a, b
    real(8), intent(out)          :: coeff(0:)
    character(len=*), optional    :: kind
    real(8), optional, intent(in) :: tol
    integer                       :: m, order, j, nq, iter
    real(8)                       :: theta, eps, g, value
    real(8)                       :: previous(0:size(coeff) - 1)
    character(len=16)             :: kernel
    !
    activate_ = .true.; if(present(activate)) activate_ = activate
    kernel    = 'none'; if(present(kind)) kernel = kind
    eps       = 1d-12 ; if(present(tol)) eps = tol
    !
    !Checks:
    select case (to_lower(str(kernel)))
     case ('jackson', 'j'); kernel = 'jackson'
     case ('none', 'n'); kernel = 'none'
     case default; stop 'cheb_coefficients: unknown kernel'
    end select
    !
    if(eps <= 0d0) stop 'cheb_coefficients: invalid tolerance'
    if(a <= 0d0 .OR. size(coeff) == 0) stop 'cheb_coefficients: invalid inputs a,b,coeff'
    !
    order = size(coeff) - 1
    theta = pi/(order + 2)
    nq    = max(64, 2*size(coeff))
    previous = 0d0
    !
    do iter=1,MaxIter
      coeff = 0d0
      do j=1,nq
        theta = pi*(j - 0.5d0)/nq
        value = fun(a*cos(theta) + b) !Any check on value?
        do m=0,order
          coeff(m) = coeff(m) + value*cos(m*theta)
        enddo
      enddo
      !
      coeff    = 2d0*coeff/nq
      coeff(0) = coeff(0)/2d0
      if(iter > 1)then
        if( maxval(abs(coeff - previous)) <= eps*max(1d0, maxval(abs(coeff))) ) exit
      end if
      previous = coeff
      if( nq > ishft(huge(nq), -1) ) stop 'cheb_coefficients: quadrature size overflow'
      nq = 2*nq
    end do
    if(iter > MaxIter) stop 'cheb_coefficients: quadratures not converged'
    theta = pi/(order + 2)
    !
    if(kernel == 'jackson')then
      do m=0,order
        g = ((order + 2 - m)*cos(m*theta) + sin(m*theta)/tan(theta))/(order + 2)
        coeff(m) = coeff(m)*g
      enddo
    endif
    !
    if(activate_) call cheb_set(a, b, coeff)
    !
  end subroutine coefficients_d

  subroutine coefficients_c(fun, a, b, coeff, kind, tol, activate)
    interface
      function fun(x) result(y)
        real(8), intent(in)         :: x
        complex(8)                  :: y
      end function fun
    end interface
    logical, optional, intent(in) :: activate
    logical                       :: activate_
    integer, parameter            :: MaxIter = 128
    real(8), intent(in)           :: a, b
    complex(8), intent(out)       :: coeff(0:)
    character(len=*), optional    :: kind
    real(8), optional, intent(in) :: tol
    integer                       :: m, order, j, nq, iter
    real(8)                       :: theta, eps, g
    complex(8)                    :: value
    complex(8)                    :: previous(0:size(coeff) - 1)
    character(len=16)             :: kernel
    !
    activate_ = .true.; if(present(activate)) activate_ = activate
    kernel    = 'none'; if(present(kind)) kernel = kind
    eps       = 1d-12 ; if(present(tol)) eps = tol
    !
    !Checks:
    select case (to_lower(str(kernel)))
     case ('jackson', 'j'); kernel = 'jackson'
     case ('none', 'n'); kernel = 'none'
     case default; stop 'cheb_coefficients: unknown kernel'
    end select
    !
    if(eps <= 0d0) stop 'cheb_coefficients: invalid tolerance'
    !
    if(a <= 0d0 .OR. size(coeff) == 0) stop 'cheb_coefficients: invalid inputs a,b,coeff'
    !
    order = size(coeff) - 1
    theta = pi/(order + 2)
    nq    = max(64, 2*size(coeff))
    previous = 0d0
    !
    do iter=1,MaxIter
      coeff = 0d0
      do j=1,nq
        theta = pi*(j - 0.5d0)/nq
        value = fun(a*cos(theta) + b)
        do m=0,order
          coeff(m) = coeff(m) + value*cos(m*theta)
        end do
      end do
      !
      coeff    = 2d0*coeff/nq
      coeff(0) = coeff(0)/2d0
      if(iter > 1)then
        if( maxval(abs(coeff - previous)) <= eps*max(1d0, maxval(abs(coeff))) ) exit
      end if
      previous = coeff
      if(nq > ishft(huge(nq), -1)) stop 'cheb_coefficients: quadrature size overflow'
      nq = 2*nq
    end do
    if( iter > MaxIter ) stop 'cheb_coefficients: quadratures not converged'
    theta = pi/(order + 2)
    !
    if( kernel == 'jackson' ) then
      do m=0,order
        g = ((order + 2 - m)*cos(m*theta) + sin(m*theta)/tan(theta))/(order + 2)
        coeff(m) = coeff(m)*g
      enddo
    endif
    !
    if(activate_) call cheb_set(a, b, coeff)
    !
  end subroutine coefficients_c









  !Set the internal memory for the Chebyshev expansion.
  !With no coefficients, prepare only the scale and invalidate the old expansion.
  subroutine set_d(a, b, coeff, info)
    real(8), intent(in)            :: a, b
    real(8), optional, intent(in)  :: coeff(0:)
    integer, optional, intent(out) :: info
    integer                        :: status !error code
    !
    status = 1
    !
    !Use block (yeah I learned some new Fortran feature... ) to implement teh ckecks in a compact way
    setup: block
      if(a <= 0d0) exit setup
      if(present(coeff)) then
        if(size(coeff) == 0) exit setup
      end if
      if(allocated(cheb_coeff)) deallocate (cheb_coeff)
      if(allocated(cheb_zcoeff)) deallocate (cheb_zcoeff)
      if(present(coeff)) then
        allocate (cheb_coeff(0:size(coeff) - 1), cheb_zcoeff(0:size(coeff) - 1))
        cheb_coeff = coeff
        cheb_zcoeff = cmplx(coeff, 0d0, 8)
      end if
      cheb_a = a
      cheb_b = b
      complex_coeff = .false.
      scale_ready = .true.
      status = 0
    end block setup
    if(present(info)) info = status
    if(status /= 0 .and. .not. present(info)) stop 'cheb_set: invalid inputs'
  end subroutine set_d

  subroutine set_c(a, b, coeff, info)
    real(8), intent(in)            :: a, b
    complex(8), intent(in)         :: coeff(0:)
    integer, optional, intent(out) :: info
    integer                        :: status
    status = 1
    setup: block
      if(a <= 0d0 .or. size(coeff) == 0) exit setup
      call set_d(a, b)
      allocate (cheb_coeff(0:size(coeff) - 1), cheb_zcoeff(0:size(coeff) - 1))
      cheb_coeff = real(coeff, 8)
      cheb_zcoeff = coeff
      complex_coeff = maxval(abs(aimag(coeff))) > 0d0
      status = 0
    end block setup
    if(present(info)) info = status
    if(status /= 0 .and. .not. present(info)) stop 'cheb_set: invalid inputs'
  end subroutine set_c











  !
  !
  !The WORKHORSE Chebyshev recurrence method
  !
  !This is the core of the algorithm: it approximates fun(A)*Vseed, knowing 
  !only Vseed\in\RRR^N and the action of A on \RRR^N vectors. 
  !The coefficiens of the Chebyshev polynomial are stored in coeff.
  !

  !Optionally it takes: the rows indices of the A matrix to be evaluated, 
  !or a matrix U containing the vector to contract with fun(A)*Vseed.
  !
  !If moments are requested, then a,b,coeff might not be set.
  !
  !Rows or U should be present, but not both. 
  !
  !DBLE:
  subroutine recurrence_d(MatVec, Vseed, values, status, nop, rows, U, moments)
    interface
      subroutine MatVec(Nloc, vin, vout)
        integer :: Nloc
        real(8) :: vin(Nloc), vout(Nloc)
      end subroutine MatVec
    end interface
    real(8), intent(in)            :: Vseed(:)
    real(8), intent(out)           :: values(:)
    integer, intent(out)           :: status, nop
    integer, optional, intent(in)  :: rows(:)
    real(8), optional, intent(in)  :: U(:, :)
    real(8), optional, intent(out) :: moments(:, 0:)
    real(8)                        :: t0(size(Vseed)), t1(size(Vseed)), t2(size(Vseed)), projected(size(values))
    integer                        :: n, m, k, order
    !
    N      = size(Vseed)
    status = 1
    nop    = 0
    values = 0d0
    if(present(moments)) moments = 0d0
    !
    !Use a block to implement checks and early exit
    evaluation: block
      if(.not. scale_ready) then
        status = 3
        exit evaluation
      end if
      if(n == 0) exit evaluation
      !
      if(present(rows)) then
        if(size(rows) /= size(values)) exit evaluation
        if(any(rows < 1) .or. any(rows > n)) exit evaluation
      else if(present(U)) then
        if(size(U, 1) /= n .or. size(U, 2) /= size(values)) exit evaluation
      else
        exit evaluation
      end if
      !
      if(present(moments)) then
        if(size(moments, 1) /= size(values) .or. size(moments, 2) == 0) exit evaluation
        order = size(moments,2)-1
      else
        if(.not. allocated(cheb_coeff)) then
          status = 3
          exit evaluation
        end if
        if(complex_coeff) exit evaluation
        order = size(cheb_coeff) - 1
      end if
      !
      !
      !Now, the actual recurrence
      status = 2
      t0     = Vseed
      do m=0,order
        !
        if(m > 0) then
          call MatVec(n, t0, t2)
          nop = nop + 1
          t2  = (t2 - cheb_b*t0)/cheb_a
          if( m > 1 )t2 = 2d0*t2 - t1
          t1 = t0
          t0 = t2
        end if
        !
        if(present(rows)) then
          projected = t0(rows)
        else
          do k = 1, size(values)
            projected(k) = dot_product(U(:, k), t0)
          end do
        end if
        !
        if(present(moments)) then
          moments(:, m) = projected
        else
          values = values + cheb_coeff(m)*projected
        end if
      end do
      !
      status = 0
    end block evaluation
    !
    !
    if(status /= 0) then
      values = 0d0
      if(present(moments)) moments = 0d0
    end if
  end subroutine recurrence_d

  !CMPLX  
  subroutine recurrence_c(MatVec, Vseed, values, status, nop, rows, U, moments)
    interface
      subroutine MatVec(Nloc, vin, vout)
        integer                     :: Nloc
        complex(8)                  :: vin(Nloc), vout(Nloc)
      end subroutine MatVec
    end interface
    complex(8), intent(in)            :: Vseed(:)
    complex(8), intent(out)           :: values(:)
    integer, intent(out)              :: status, nop
    integer, optional, intent(in)     :: rows(:)
    complex(8), optional, intent(in)  :: U(:, :)
    complex(8), optional, intent(out) :: moments(:, 0:)
    complex(8)                        :: t0(size(Vseed)), t1(size(Vseed)), t2(size(Vseed)), projected(size(values))
    integer                           :: n, m, k, order

    n      = size(Vseed)
    status = 1
    nop    = 0
    values = 0d0
    if(present(moments)) moments = 0d0
    !
    evaluation: block
      if(.not. scale_ready) then
        status = 3
        exit evaluation
      end if
      if(n == 0) exit evaluation
      !
      if(present(rows)) then
        if(size(rows) /= size(values)) exit evaluation
        if(any(rows < 1) .or. any(rows > n)) exit evaluation
      else if(present(U)) then
        if(size(U, 1) /= n .or. size(U, 2) /= size(values)) exit evaluation
      else
        exit evaluation
      end if
      !
      if(present(moments)) then
        if(size(moments, 1) /= size(values) .or. size(moments, 2) == 0) exit evaluation
        order = size(moments, 2) - 1
      else
        if(.not. allocated(cheb_coeff)) then
          status = 3
          exit evaluation
        end if
        order = size(cheb_coeff) - 1
      end if
      !
      status = 2
      t0     = Vseed
      do m=0,order
        !
        if(m > 0) then
          call MatVec(n, t0, t2)
          nop = nop + 1
          t2  = (t2 - cheb_b*t0)/cheb_a
          if(m > 1) t2 = 2d0*t2 - t1
          t1 = t0
          t0 = t2
        end if
        !
        if(present(rows)) then
          projected = t0(rows)
        else
          do k = 1, size(values)
            projected(k) = dot_product(U(:, k), t0)
          end do
        end if
        !
        if(present(moments)) then
          moments(:, m) = projected
        else
          values = values + cheb_zcoeff(m)*projected
        end if
        !
      end do
      status = 0
    end block evaluation
    !
    if(status /= 0) then
      values = 0d0
      if(present(moments)) moments = 0d0
    end if
    !
  end subroutine recurrence_c









  !New use simpler interfaces for the most common cases:
  !
  !Passing Rows:
  subroutine projection_rows_d(MatVec, Vseed, rows, values, info, NumOp)
    interface
      subroutine MatVec(Nloc, vin, vout)
        integer :: Nloc
        real(8) :: vin(Nloc), vout(Nloc)
      end subroutine MatVec
    end interface
    real(8), intent(in)            :: Vseed(:)
    integer, intent(in)            :: rows(:)
    real(8), intent(out)           :: values(:)
    integer, optional, intent(out) :: info, NumOp
    integer                        :: status, nop
    call recurrence_d(MatVec, Vseed, values, status, nop, rows=rows)
    if(present(info)) info = status
    if(present(NumOp)) NumOp = nop
    if(status /= 0 .and. .not. present(info)) stop 'cheb_projection: invalid input, state or recurrence'
  end subroutine projection_rows_d

  subroutine projection_rows_c(MatVec, Vseed, rows, values, info, NumOp)
    interface
      subroutine MatVec(Nloc, vin, vout)
        integer                     :: Nloc
        complex(8)                  :: vin(Nloc), vout(Nloc)
      end subroutine MatVec
    end interface
    complex(8), intent(in)           :: Vseed(:)
    integer, intent(in)              :: rows(:)
    complex(8), intent(out)          :: values(:)
    integer, optional, intent(out)    :: info, NumOp
    integer                         :: status, nop
    call recurrence_c(MatVec, Vseed, values, status, nop, rows=rows)
    if(present(info)) info = status
    if(present(NumOp)) NumOp = nop
    if(status /= 0 .and. .not. present(info)) stop 'cheb_projection: invalid input, state or recurrence'
  end subroutine projection_rows_c





  !Passing vectors in U:
  subroutine projection_vectors_d(MatVec, Vseed, U, values, info, NumOp)
    interface
      subroutine MatVec(Nloc, vin, vout)
        integer :: Nloc
        real(8) :: vin(Nloc), vout(Nloc)
      end subroutine MatVec
    end interface
    real(8), intent(in)            :: Vseed(:)
    real(8), intent(in)            :: U(:, :)
    real(8), intent(out)           :: values(:)
    integer, optional, intent(out) :: info, NumOp
    integer                        :: status, nop
    call recurrence_d(MatVec, Vseed, values, status, nop, U=U)
    if(present(info)) info = status
    if(present(NumOp)) NumOp = nop
    if(status /= 0 .and. .not. present(info)) stop 'cheb_projection: invalid input, state or recurrence'
  end subroutine projection_vectors_d

  subroutine projection_vectors_c(MatVec, Vseed, U, values, info, NumOp)
    interface
      subroutine MatVec(Nloc, vin, vout)
        integer                     :: Nloc
        complex(8)                  :: vin(Nloc), vout(Nloc)
      end subroutine MatVec
    end interface
    complex(8), intent(in)           :: Vseed(:)
    complex(8), intent(in)           :: U(:, :)
    complex(8), intent(out)          :: values(:)
    integer, optional, intent(out)    :: info, NumOp
    integer                         :: status, nop
    call recurrence_c(MatVec, Vseed, values, status, nop, U=U)
    if(present(info)) info = status
    if(present(NumOp)) NumOp = nop
    if(status /= 0 .and. .not. present(info)) stop 'cheb_projection: invalid input, state or recurrence'
  end subroutine projection_vectors_c





  !Get the moments for selected rows (independent of a,b,coeff)
  subroutine moments_rows_d(MatVec, Vseed, rows, moments, info, NumOp)
    interface
      subroutine MatVec(Nloc, vin, vout)
        integer :: Nloc
        real(8) :: vin(Nloc), vout(Nloc)
      end subroutine MatVec
    end interface
    real(8), intent(in)            :: Vseed(:)
    integer, intent(in)            :: rows(:)
    real(8), intent(out)           :: moments(:, 0:)
    integer, optional, intent(out) :: info, NumOp
    real(8)                        :: values(size(moments, 1))
    integer                        :: status, nop
    call recurrence_d(MatVec, Vseed, values, status, nop, rows=rows, moments=moments)
    if(present(info)) info = status
    if(present(NumOp)) NumOp = nop
    if(status /= 0 .and. .not. present(info)) stop 'cheb_moments: invalid input, state or recurrence'
  end subroutine moments_rows_d

  subroutine moments_rows_c(MatVec, Vseed, rows, moments, info, NumOp)
    interface
      subroutine MatVec(Nloc, vin, vout)
        integer                     :: Nloc
        complex(8)                  :: vin(Nloc), vout(Nloc)
      end subroutine MatVec
    end interface
    complex(8), intent(in)           :: Vseed(:)
    integer, intent(in)              :: rows(:)
    complex(8), intent(out)          :: moments(:, 0:)
    integer, optional, intent(out)    :: info, NumOp
    complex(8)                      :: values(size(moments, 1))
    integer                         :: status, nop
    call recurrence_c(MatVec, Vseed, values, status, nop, rows=rows, moments=moments)
    if(present(info)) info = status
    if(present(NumOp)) NumOp = nop
    if(status /= 0 .and. .not. present(info)) stop 'cheb_moments: invalid input, state or recurrence'
  end subroutine moments_rows_c




  !Get the moments for given matrix U
  subroutine moments_vectors_d(MatVec, Vseed, U, moments, info, NumOp)
    interface
      subroutine MatVec(Nloc, vin, vout)
        integer :: Nloc
        real(8) :: vin(Nloc), vout(Nloc)
      end subroutine MatVec
    end interface
    real(8), intent(in)            :: Vseed(:)
    real(8), intent(in)            :: U(:, :)
    real(8), intent(out)           :: moments(:, 0:)
    integer, optional, intent(out) :: info, NumOp
    real(8)                        :: values(size(moments, 1))
    integer                        :: status, nop
    call recurrence_d(MatVec, Vseed, values, status, nop, U=U, moments=moments)
    if(present(info)) info = status
    if(present(NumOp)) NumOp = nop
    if(status /= 0 .and. .not. present(info)) stop 'cheb_moments: invalid input, state or recurrence'
  end subroutine moments_vectors_d

  subroutine moments_vectors_c(MatVec, Vseed, U, moments, info, NumOp)
    interface
      subroutine MatVec(Nloc, vin, vout)
        integer                     :: Nloc
        complex(8)                  :: vin(Nloc), vout(Nloc)
      end subroutine MatVec
    end interface
    complex(8), intent(in)           :: Vseed(:)
    complex(8), intent(in)           :: U(:, :)
    complex(8), intent(out)          :: moments(:, 0:)
    integer, optional, intent(out)    :: info, NumOp
    complex(8)                      :: values(size(moments, 1))
    integer                         :: status, nop
    call recurrence_c(MatVec, Vseed, values, status, nop, U=U, moments=moments)
    if(present(info)) info = status
    if(present(NumOp)) NumOp = nop
    if(status /= 0 .and. .not. present(info)) stop 'cheb_moments: invalid input, state or recurrence'
  end subroutine moments_vectors_c





end module SF_CHEBYCHEV
