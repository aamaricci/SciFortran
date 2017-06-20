module IOREAD
  USE IOFILE
  implicit none
  private

  integer            :: unit
  character(len=128) :: fmt
  logical            :: control

  interface sread
     module procedure :: sreadA1_RR
     module procedure :: sreadA1_RC
     module procedure :: sreadA2_RR
     module procedure :: sreadA2_RC
     module procedure :: sreadA3_RR
     module procedure :: sreadA3_RC
     module procedure :: sreadA4_RR
     module procedure :: sreadA4_RC
     module procedure :: sreadA5_RR
     module procedure :: sreadA5_RC
     module procedure :: sreadA6_RR
     module procedure :: sreadA6_RC
     module procedure :: sreadA7_RR
     module procedure :: sreadA7_RC
  end interface sread


  interface read_data
     module procedure :: data_readA1_R
     module procedure :: data_readA1_C
     module procedure :: data_readA2_R
     module procedure :: data_readA2_C
     module procedure :: data_readA3_R
     module procedure :: data_readA3_C
     module procedure :: data_readA4_R
     module procedure :: data_readA4_C
     module procedure :: data_readA5_R
     module procedure :: data_readA5_C
     module procedure :: data_readA6_R
     module procedure :: data_readA6_C
     module procedure :: data_readA7_R
     module procedure :: data_readA7_C
  end interface read_data


  public :: sread
  public :: read_data

contains



  ! SPLOT ararys (1--7)
  include "ioread_sread.f90"


  ! STORE_DATA arrays (1--7)
  include "ioread_read_data.f90"


end module IOREAD
