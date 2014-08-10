module MPI_VARS
  implicit none
  integer          :: MPIID=0
  integer          :: MPISIZE=1
  integer          :: MPIERR
  character(len=3) :: MPICHAR
  ! public :: init_mpi
  ! public :: finalize_mpi
  ! contains
  !   subroutine init_mpi
  !     call MPI_INIT(mpiERR)
  !     call MPI_COMM_RANK(MPI_COMM_WORLD,mpiID,mpiERR)
  !     call MPI_COMM_SIZE(MPI_COMM_WORLD,mpiSIZE,mpiERR)
  !     write(*,"(A,I4,A,I4,A)")'Processor ',mpiID,' of ',mpiSIZE,' is alive'
  !     call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  !   end subroutine init_mpi
  !   subroutine finalize_mpi
  !     call MPI_FINALIZE(mpiERR)
  !   end subroutine finalize_mpi
end module MPI_VARS

module OMP_VARS
  implicit none
  integer :: omp_num_threads
  integer :: omp_id
  integer :: omp_size
end module OMP_VARS


module CONSTANTS
  implicit none
  private

  !PARAMETERS
  !===============================================================
  complex(8),parameter,public :: zero=(0.d0,0.d0)
  complex(8),parameter,public :: xi=(0.d0,1.d0)
  complex(8),parameter,public :: one=(1.d0,0.d0)
  real(8),parameter,public    :: sqrt2 = 1.41421356237309504880169d0
  real(8),parameter,public    :: sqrt3 = 1.73205080756887729352745d0
  real(8),parameter,public    :: sqrt6 = 2.44948974278317809819728d0
  real(8),parameter,public    :: pi    = 3.14159265358979323846264338327950288419716939937510d0
  real(8),parameter,public    :: pi2   = 6.28318530717959d0
  real(8),parameter,public    :: gamma_euler = 0.57721566490153286060d0  !euler s constant
  real(8),parameter,public    :: euler= 2.7182818284590452353602874713526624977572470936999595749669676277240766303535d0
  integer,parameter,public    :: max_int  = huge(1) 
  real(8),parameter,public    :: max_real = huge(1.d0)
  real(8),parameter,public    :: epsilonr=epsilon(1.d0),epsilonq=1.d-30
  integer,parameter,public    :: dbl=8,dp=8        ! "double" precision
  integer,parameter,public    :: ddp=16            ! "quad"   precision
  integer,parameter,public    :: sp = kind(1.0)    ! "single" precision


  real(8),parameter,public ::                                 Avogadro_constant=  0.602214129000D+24
  real(8),parameter,public ::                                     Bohr_magneton=  0.927400968000D-23
  real(8),parameter,public ::                             Bohr_magneton_in_eVoT=  0.578838180660D-04
  real(8),parameter,public ::                             Bohr_magneton_in_HzoT=  0.139962455500D+11
  real(8),parameter,public ::         Bohr_magneton_in_inverse_meters_per_tesla=     46.686449800000
  real(8),parameter,public ::                              Bohr_magneton_in_KoT=      0.671713880000
  real(8),parameter,public ::                                       Bohr_radius=  0.529177210920D-10
  real(8),parameter,public ::                                Boltzmann_constant=  0.138064880000D-22
  real(8),parameter,public ::                        Boltzmann_constant_in_eVoK=  0.861733240000D-04
  real(8),parameter,public ::                        Boltzmann_constant_in_HzoK=  0.208366180000D+11
  real(8),parameter,public ::   Boltzmann_constant_in_inverse_meters_per_kelvin=     69.503476000000
  real(8),parameter,public ::                                Compton_wavelength=  0.242631023890D-11
  real(8),parameter,public ::                      Compton_wavelength_over_2_pi=  0.386159268000D-12
  real(8),parameter,public ::                                 electric_constant=  0.885418781700D-11
  real(8),parameter,public ::                  electron_charge_to_mass_quotient= -0.175882008800D+12
  real(8),parameter,public ::                                 electron_g_factor= -0.200231930436D+01
  real(8),parameter,public ::                           electron_gyromag__ratio=  0.176085970800D+12
  real(8),parameter,public ::                 electron_gyromag__ratio_over_2_pi=  0.280249526600D+05
  real(8),parameter,public ::                                electron_mag__mom_= -0.928476430000D-23
  real(8),parameter,public ::         electron_mag__mom__to_Bohr_magneton_ratio= -0.100115965218D+01
  real(8),parameter,public ::                                     electron_mass=  0.910938291000D-30
  real(8),parameter,public ::                   electron_mass_energy_equivalent=  0.818710506000D-13
  real(8),parameter,public ::            electron_mass_energy_equivalent_in_MeV=      0.510998928000
  real(8),parameter,public ::                                     electron_volt=  0.160217656500D-18
  real(8),parameter,public ::       electron_volt_atomic_mass_unit_relationship=  0.107354415000D-08
  real(8),parameter,public ::                electron_volt_hartree_relationship=      0.036749323790
  real(8),parameter,public ::                  electron_volt_hertz_relationship=  0.241798934800D+15
  real(8),parameter,public ::          electron_volt_inverse_meter_relationship=  0.806554429000D+06
  real(8),parameter,public ::                  electron_volt_joule_relationship=  0.160217656500D-18
  real(8),parameter,public ::                 electron_volt_kelvin_relationship=  0.116045190000D+05
  real(8),parameter,public ::               electron_volt_kilogram_relationship=  0.178266184500D-35
  real(8),parameter,public ::                                 elementary_charge=  0.160217656500D-18
  real(8),parameter,public ::                          elementary_charge_over_h=  0.241798934800D+15
  real(8),parameter,public ::                                  Faraday_constant=  0.964853365000D+05
  real(8),parameter,public ::Faraday_constant_for_conventional_electric_current=  0.964853321000D+05
  real(8),parameter,public ::                           fine_structure_constant=  0.729735256980D-02
  real(8),parameter,public ::                                Josephson_constant=  0.483597870000D+15
  real(8),parameter,public ::                  joule_electron_volt_relationship=  0.624150934000D+19
  real(8),parameter,public ::                          joule_hertz_relationship=  0.150919031100D+34
  real(8),parameter,public ::                  joule_inverse_meter_relationship=  0.503411701000D+25
  real(8),parameter,public ::                         joule_kelvin_relationship=  0.724297160000D+23
  real(8),parameter,public ::                       joule_kilogram_relationship=  0.111265005600D-16
  real(8),parameter,public ::              kelvin_atomic_mass_unit_relationship=  0.925108680000D-13
  real(8),parameter,public ::                 kelvin_electron_volt_relationship=  0.861733240000D-04
  real(8),parameter,public ::                       kelvin_hartree_relationship=  0.316681140000D-05
  real(8),parameter,public ::                         kelvin_hertz_relationship=  0.208366180000D+11
  real(8),parameter,public ::                 kelvin_inverse_meter_relationship=     69.503476000000
  real(8),parameter,public ::                         kelvin_joule_relationship=  0.138064880000D-22
  real(8),parameter,public ::                      kelvin_kilogram_relationship=  0.153617900000D-39
  real(8),parameter,public ::            kilogram_atomic_mass_unit_relationship=  0.602214129000D+27
  real(8),parameter,public ::               kilogram_electron_volt_relationship=  0.560958885000D+36
  real(8),parameter,public ::                     kilogram_hartree_relationship=  0.206148596800D+35
  real(8),parameter,public ::                       kilogram_hertz_relationship=  0.135639260800D+50
  real(8),parameter,public ::               kilogram_inverse_meter_relationship=  0.452443873000D+42
  real(8),parameter,public ::                       kilogram_joule_relationship=  0.898755178700D+17
  real(8),parameter,public ::                      kilogram_kelvin_relationship=  0.650965820000D+40
  real(8),parameter,public ::                      lattice_parameter_of_silicon=  0.543102050400D-09
  real(8),parameter,public ::                            natural_unit_of_action=  0.105457172600D-33
  real(8),parameter,public ::                    natural_unit_of_action_in_eV_s=  0.658211928000D-15
  real(8),parameter,public ::                            natural_unit_of_energy=  0.818710506000D-13
  real(8),parameter,public ::                     natural_unit_of_energy_in_MeV=      0.510998928000
  real(8),parameter,public ::                            natural_unit_of_length=  0.386159268000D-12
  real(8),parameter,public ::                              natural_unit_of_mass=  0.910938291000D-30
  real(8),parameter,public ::                            natural_unit_of_mom_um=  0.273092429000D-21
  real(8),parameter,public ::                   natural_unit_of_mom_um_in_MeVoc=      0.510998928000
  real(8),parameter,public ::                              natural_unit_of_time=  0.128808866833D-20
  real(8),parameter,public ::                          natural_unit_of_velocity=  0.299792458000D+09
  real(8),parameter,public ::                 Newtonian_constant_of_gravitation=  0.667384000000D-10
  real(8),parameter,public ::                                   Planck_constant=  0.662606957000D-33
  real(8),parameter,public ::                           Planck_constant_in_eV_s=  0.413566751600D-14
  real(8),parameter,public ::                         Planck_constant_over_2_pi=  0.105457172600D-33
  real(8),parameter,public ::                                  Rydberg_constant=  0.109737315685D+08
  real(8),parameter,public ::                    Rydberg_constant_times_c_in_Hz=  0.328984196036D+16
  real(8),parameter,public ::                   Rydberg_constant_times_hc_in_eV=     13.605692530000
  real(8),parameter,public ::                    Rydberg_constant_times_hc_in_J=  0.217987217100D-17
  real(8),parameter,public ::                          speed_of_light_in_vacuum=  0.299792458000D+09
  real(8),parameter,public ::                  standard_acceleration_of_gravity=      9.806650000000
  real(8),parameter,public ::                         Stefan_Boltzmann_constant=  0.567037300000D-07





  public :: timestamp



contains



  !+-------------------------------------------------------------------+
  !PURPOSE  : prints the current YMDHMS date as a time stamp.
  ! Example: 31 May 2001   9:45:54.872 AM
  !+-------------------------------------------------------------------+
  subroutine timestamp(unit)
    integer,optional        :: unit
    integer                 :: unit_
    integer(4),dimension(8) :: data
    unit_=6;if(present(unit))unit_=unit
    call date_and_time(values=data)
    call print_date(data,unit_)
  end subroutine timestamp



  !+-------------------------------------------------------------------+
  !PURPOSE  : print actual date
  !+-------------------------------------------------------------------+
  subroutine print_date(dummy,unit)
    integer(4),dimension(8) :: dummy
    integer                 :: unit
    integer(4)                          :: year
    integer(4)                          :: mese
    integer(4)                          :: day
    integer(4)                          :: h
    integer(4)                          :: m
    integer(4)                          :: s
    integer(4)                          :: ms
    character(len=9),parameter,dimension(12) :: month = (/ &
         'January  ', 'February ', 'March    ', 'April    ', &
         'May      ', 'June     ', 'July     ', 'August   ', &
         'September', 'October  ', 'November ', 'December ' /)
    year = dummy(1)
    mese = dummy(2)
    day  = dummy(3)
    h    = dummy(5)
    m    = dummy(6)
    s    = dummy(7)
    ms   = dummy(8)
    write(unit,"(A,i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3)")&
         "Timestamp: +",day,trim(month(mese)),year, h,':',m,':',s,'.',ms
    write(unit,*)""
  end subroutine print_date



  !******************************************************************
  !******************************************************************
  !******************************************************************



END MODULE CONSTANTS


MODULE FONTS
  implicit none
  private


  public :: bold,underline,highlight,erased
  public :: red,green,yellow,blue
  public :: bold_red,bold_green,bold_yellow,bold_blue
  public :: bg_red,bg_green,bg_yellow,bg_blue

contains

  !******************************************************************
  !******************************************************************
  !******************************************************************

  function bold(text) result(textout)
    character(len=*) :: text
    character(len=8+len(text)) :: textout
    textout=achar(27)//"[1m"//text//achar(27)//"[0m"
  end function bold

  function underline(text) result(textout)
    character(len=*) :: text
    character(len=8+len(text)) :: textout
    textout=achar(27)//"[4m"//text//achar(27)//"[0m"
  end function underline

  function highlight(text) result(textout)
    character(len=*) :: text
    character(len=8+len(text)) :: textout
    textout=achar(27)//"[7m"//text//achar(27)//"[0m"
  end function highlight

  function erased(text) result(textout)
    character(len=*) :: text
    character(len=8+len(text)) :: textout
    textout=achar(27)//"[9m"//text//achar(27)//"[0m"
  end function erased

  function red(text) result(textout)
    character(len=*) :: text
    character(len=9+len(text)) :: textout
    textout=achar(27)//"[91m"//text//achar(27)//"[0m"
  end function red

  function green(text) result(textout)
    character(len=*) :: text
    character(len=9+len(text)) :: textout
    textout=achar(27)//"[92m"//text//achar(27)//"[0m"
  end function green

  function yellow(text) result(textout)
    character(len=*) :: text
    character(len=9+len(text)) :: textout
    textout=achar(27)//"[93m"//text//achar(27)//"[0m"
  end function yellow

  function blue(text) result(textout)
    character(len=*) :: text
    character(len=9+len(text)) :: textout
    textout=achar(27)//"[94m"//text//achar(27)//"[0m"
  end function blue

  function bold_red(text) result(textout)
    character(len=*) :: text
    character(len=11+len(text)) :: textout
    textout=achar(27)//"[1;91m"//text//achar(27)//"[0m"
  end function bold_red

  function bold_green(text) result(textout)
    character(len=*) :: text
    character(len=11+len(text)) :: textout
    textout=achar(27)//"[1;92m"//text//achar(27)//"[0m"
  end function bold_green

  function bold_yellow(text) result(textout)
    character(len=*) :: text
    character(len=11+len(text)) :: textout
    textout=achar(27)//"[1;93m"//text//achar(27)//"[0m"
  end function bold_yellow

  function bold_blue(text) result(textout)
    character(len=*) :: text
    character(len=11+len(text)) :: textout
    textout=achar(27)//"[1;94m"//text//achar(27)//"[0m"
  end function bold_blue

  function bg_red(text) result(textout)
    character(len=*) :: text
    character(len=9+len(text)) :: textout
    textout=achar(27)//"[41m"//text//achar(27)//"[0m"
  end function bg_red

  function bg_green(text) result(textout)
    character(len=*) :: text
    character(len=9+len(text)) :: textout
    textout=achar(27)//"[42m"//text//achar(27)//"[0m"
  end function bg_green

  function bg_yellow(text) result(textout)
    character(len=*) :: text
    character(len=9+len(text)) :: textout
    textout=achar(27)//"[43m"//text//achar(27)//"[0m"
  end function bg_yellow

  function bg_blue(text) result(textout)
    character(len=*) :: text
    character(len=9+len(text)) :: textout
    textout=achar(27)//"[44m"//text//achar(27)//"[0m"
  end function bg_blue

END MODULE FONTS
