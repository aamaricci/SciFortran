module precision
  integer, parameter :: kr4 = selected_real_kind(6,37)       ! single precision real
  integer, parameter :: kr8 = selected_real_kind(15,307)     ! double precision real
  integer, parameter :: ki4 = selected_int_kind(9)           ! single precision integer
  integer, parameter :: ki8 = selected_int_kind(18)          ! double precision integer
  integer, parameter :: kc4 = kr4                            ! single precision complex
  integer, parameter :: kc8 = kr8                            ! double precision complex
end module precision
