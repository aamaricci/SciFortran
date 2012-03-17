module GREENFUNX
  !###############################################################
  !     PROGRAM  : GREENFUNX
  !     TYPE     : Module
  !     PURPOSE  : a library for Green's functions type and related 
  !operations. This comes really handy when FFT must be used as it
  !automatically set the correct dimensions.
  !###############################################################
  implicit none
  private 

  type,public :: matsubara_gf
     complex(8),dimension(:),pointer :: iw
     real(8),dimension(:),pointer    :: tau
     logical                         :: status=.false.
  end type matsubara_gf

  type,public ::  real_gf
     complex(8),dimension(:),pointer :: w
     complex(8),dimension(:),pointer :: t
     logical                         :: status=.false.
  end type real_gf

  type,public ::  keldysh_equilibrium_component
     complex(8),dimension(:),pointer :: w
     complex(8),dimension(:),pointer :: t
  end type keldysh_equilibrium_component
  type,public ::  keldysh_equilibrium_gf
     type(keldysh_equilibrium_component) :: less
     type(keldysh_equilibrium_component) :: gtr
     type(keldysh_equilibrium_component) :: ret
     logical                             :: status=.false.
  end type keldysh_equilibrium_gf


  interface allocate_gf
     module procedure allocate_matsubara_gf,allocate_matsubara_gf_,allocate_matsubara_gf__,&
          allocate_real_gf,allocate_real_gf_,allocate_real_gf__,&
          allocate_keldysh_equilibrium_component,allocate_keldysh_equilibrium_gf
  end interface allocate_gf

  interface deallocate_gf
     module procedure deallocate_matsubara_gf,deallocate_matsubara_gf_,deallocate_matsubara_gf__,&
          deallocate_real_gf, deallocate_real_gf_, deallocate_real_gf__,&          
          deallocate_keldysh_equilibrium_component,deallocate_keldysh_equilibrium_gf
  end interface deallocate_gf

  interface assignment(=)
     module procedure matsubara_gf_identity,matsubara_gf_identity_,matsubara_gf_identity__,&
          real_gf_identity,real_gf_identity_,real_gf_identity__,&
          keldysh_eq_component_gf_identity,keldysh_eq_gf_identity,&          
          matsubara_scalar_identity,matsubara_scalar_identity_,matsubara_scalar_identity__,&
          real_scalar_identity,real_scalar_identity_,real_scalar_identity__,&
          keldysh_eq_component_scalar_identity,keldysh_eq_scalar_identity
  end interface assignment(=)

  interface operator(+)
     module procedure add_matsubara_gf,add_real_gf,&
          add_keldysh_equilibrium_component,add_keldysh_equilibrium_gf
  end interface operator(+)

  public :: allocate_gf
  public :: deallocate_gf
  public :: assignment(=)
  public :: operator(+)
  public :: ret_component_t
  public :: less_component_w
  public :: gtr_component_w

contains

  !+----------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : subroutine
  !PURPOSE  : 
  !+----------------------------------------------------------------+
  function ret_component_t(fgkgtr,fgkless,t,L) result(fgkret)
    complex(8),dimension(-L:L),intent(in)  :: fgkless,fgkgtr
    complex(8),dimension(-L:L)             :: fgkret
    real(8),dimension(-L:L),intent(in)     :: t
    integer :: i,L
    forall(i=-L:L)fgkret(i)=step(t(i))*(fgkgtr(i)-fgkless(i))
  contains
    pure function step(x)
      real(8),intent(in) :: x
      real(8) :: step
      if(x < 0.d0) then
         step = 0.0d0
      elseif(x==0.d0)then
         step = 0.50d0
      else
         step = 1.0d0
      endif
    end function step
  end function ret_component_t
  !--------------------------------------------------------!


  !--------------------------------------------------------!
  function less_component_w(fret,wr,beta) result(fless)
    complex(8),dimension(:),intent(in)  :: fret
    complex(8),dimension(size(fret))    :: fless
    real(8),dimension(size(fret))       :: wr
    real(8) :: pi,A,beta,w
    integer :: i,L
    pi=acos(-1.d0)
    L=size(fret)
    do i=1,L
       w       = wr(i)
       A       = -dimag(fret(i))/pi
       fless(i)= fermi(w,beta)*A
    enddo
  contains
    function fermi(x,beta)
      real(8) :: fermi, x, beta
      if(x*beta > 100.d0)then
         fermi=0.d0
         return
      endif
      fermi = 1.d0/(1.d0+exp(beta*x))
    end function fermi
  end function less_component_w
  !--------------------------------------------------------!


  !--------------------------------------------------------!
  function gtr_component_w(fret,wr,beta) result(fgtr)
    complex(8),dimension(:),intent(in)  :: fret
    complex(8),dimension(size(fret))    :: fgtr
    real(8),dimension(size(fret))       :: wr
    real(8) :: pi,A,beta,w
    integer :: i,L
    pi=acos(-1.d0)
    L=size(fret)
    do i=1,L
       w      = wr(i)
       A      = -dimag(fret(i))/pi
       fgtr(i)= (fermi(w,beta)-1.d0)*A
    enddo
  contains
    function fermi(x,beta)
      real(8) :: fermi, x, beta
      if(x*beta > 100.d0)then
         fermi=0.d0
         return
      endif
      fermi = 1.d0/(1.d0+exp(beta*x))
    end function fermi
  end function gtr_component_w
  !--------------------------------------------------------!



  !+----------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : subroutine
  !PURPOSE  : 
  !+----------------------------------------------------------------+
  include "allocate_gf.f90"
  include "deallocate_gf.f90"
  !******************************************************************
  !******************************************************************
  !******************************************************************


  !+----------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : subroutine
  !PURPOSE  : 
  !+----------------------------------------------------------------+
  include "gf_identity.f90"
  include "scalar_identity.f90"
  !******************************************************************
  !******************************************************************
  !******************************************************************


  !+----------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : subroutine
  !PURPOSE  : 
  !+----------------------------------------------------------------+
  elemental function add_matsubara_gf(fgk1,fgk2) result(fgk3)
    type(matsubara_gf),intent(in) :: fgk1,fgk2
    type(matsubara_gf)            :: fgk3
    fgk3%iw  = fgk1%iw  + fgk2%iw
    fgk3%tau = fgk1%tau + fgk2%tau
  end function add_matsubara_gf

  elemental function add_real_gf(fgk1,fgk2) result(fgk3)
    type(real_gf),intent(in) :: fgk1,fgk2
    type(real_gf)            :: fgk3
    fgk3%w = fgk1%w + fgk2%w
    fgk3%t  = fgk1%t  + fgk2%t
  end function add_real_gf

  elemental function add_keldysh_equilibrium_component(fgk1,fgk2) result(fgk3)
    type(keldysh_equilibrium_component),intent(in) :: fgk1,fgk2
    type(keldysh_equilibrium_component)            :: fgk3
    fgk3%t = fgk1%t + fgk2%t
    fgk3%w = fgk1%w  + fgk2%w
  end function add_keldysh_equilibrium_component

  elemental function add_keldysh_equilibrium_gf(fgk1,fgk2) result(fgk3)
    type(keldysh_equilibrium_gf),intent(in) :: fgk1,fgk2
    type(keldysh_equilibrium_gf)            :: fgk3
    fgk3%less%t = fgk1%less%t + fgk2%less%t
    fgk3%gtr%t  = fgk1%gtr%t  + fgk2%gtr%t
    fgk3%ret%t  = fgk1%ret%t  + fgk2%ret%t
    fgk3%less%w = fgk1%less%w + fgk2%less%w
    fgk3%gtr%w  = fgk1%gtr%w  + fgk2%gtr%w
    fgk3%ret%w  = fgk1%ret%w  + fgk2%ret%w
  end function add_keldysh_equilibrium_gf
  !******************************************************************
  !******************************************************************
  !******************************************************************

end module GREENFUNX
