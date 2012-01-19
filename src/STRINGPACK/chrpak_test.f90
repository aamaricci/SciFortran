program main

!*****************************************************************************80
!
!! MAIN is the main program for CHRPAK_PRB.
!
!  Discussion:
!
!    CHRPAK_PRB tests routines from the CHRPAK library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CHRPAK_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test the routines in the CHRPAK library.'

  call test001 ( )
  call test002 ( )
  call test003 ( )
  call test005 ( )
  call test006 ( )
  call test007 ( )
  call test008 ( )

  call test011 ( )
  call test012 ( )
  call test013 ( )
  call test014 ( )
  call test015 ( )
  call test016 ( )
  call test017 ( )
  call test018 ( )
  call test019 ( )

  call test020 ( )
  call test021 ( )
  call test022 ( )
  call test023 ( )
  call test024 ( )
  call test025 ( )
  call test026 ( )
  call test027 ( )
  call test028 ( )
  call test029 ( )

  call test030 ( )
  call test031 ( )
  call test032 ( )
  call test033 ( )
  call test034 ( )
  call test035 ( )
  call test036 ( )
  call test037 ( )
  call test038 ( )
  call test039 ( )

  call test040 ( )
  call test041 ( )
  call test042 ( )
  call test043 ( )
  call test045 ( )
  call test046 ( )
  call test047 ( )
  call test048 ( )
  call test049 ( )

  call test050 ( )
  call test051 ( )
  call test052 ( )
  call test054 ( )
  call test055 ( )
  call test056 ( )
  call test057 ( )
  call test058 ( )
  call test059 ( )

  call test060 ( )
  call test061 ( )
  call test062 ( )
  call test063 ( )
  call test064 ( )
  call test065 ( )
  call test066 ( )
  call test067 ( )
  call test068 ( )

  call test070 ( )
  call test071 ( )
  call test072 ( )
  call test073 ( )
  call test074 ( )
  call test075 ( )
  call test076 ( )
  call test077 ( )
  call test078 ( )
  call test079 ( )

  call test080 ( )
  call test081 ( )
  call test082 ( )
  call test083 ( )
  call test085 ( )
  call test086 ( )
  call test087 ( )
  call test088 ( )
  call test089 ( )

  call test090 ( )
  call test091 ( )
  call test092 ( )
  call test093 ( )
  call test094 ( )
  call test095 ( )
  call test096 ( )
  call test097 ( )
  call test0975 ( )
  call test098 ( )
  call test099 ( )

  call test100 ( )
  call test101 ( )
  call test1013 ( )
  call test1015 ( )
  call test102 ( )
  call test103 ( )
  call test104 ( )
  call test105 ( )
  call test1055 ( )
  call test106 ( )
  call test107 ( )
  call test108 ( )
  call test109 ( )

  call test110 ( )
  call test111 ( )
  call test112 ( )
  call test113 ( )
  call test114 ( )
  call test115 ( )
  call test116 ( )
  call test117 ( )
  call test118 ( )
  call test119 ( )

  call test120 ( )
  call test121 ( )
  call test122 ( )
  call test1225 ( )
  call test123 ( )
  call test124 ( )
  call test125 ( )
  call test1255 ( )
  call test126 ( )
  call test127 ( )
  call test128 ( )
  call test129 ( )

  call test130 ( )
  call test131 ( )
  call test132 ( )
  call test133 ( )
  call test134 ( )
  call test135 ( )
  call test136 ( )
  call test137 ( )
  call test138 ( )
  call test139 ( )

  call test140 ( )
  call test141 ( )
  call test142 ( )
  call test143 ( )
  call test144 ( )
  call test145 ( )
  call test146 ( )
  call test147 ( )
  call test148 ( )
  call test149 ( )

  call test150 ( )
  call test152 ( )
  call test153 ( )
  call test154 ( )
  call test155 ( )
  call test156 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CHRPAK_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test001 ( )

!*****************************************************************************80
!
!! TEST001 tests A_TO_I4 and I4_TO_A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character a
  integer ( kind = 4 ) a_to_i4
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  character i4_to_a

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST001'
  write ( *, '(a)' ) '  A_TO_I4: Alphabetic character => I'
  write ( *, '(a)' ) '  I4_TO_A: I => Alphabetic character'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  1:26 = A:Z'
  write ( *, '(a)' ) '  27:52 = a:z'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   I  ==>  A  ==>  I'
  write ( *, '(a)' ) ' '

  do i = 0, 55, 3

    a = i4_to_a ( i )
    i2 = a_to_i4 ( a )

    write ( *, '(i8,5x,a1,5x,i8)' ) i, a, i2

  end do

  return
end
subroutine test002 ( )

!*****************************************************************************80
!
!! TEST002 tests B4_IEEE_TO_R4 and R4_TO_B4_IEEE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 16

  character ( len = 32 ) bits
  integer ( kind = 4 ) e
  integer ( kind = 4 ) f
  character ( len = 32 ) i4_to_s32
  real    ( kind = 4 ) r1
  real    ( kind = 4 ) r2
  real    ( kind = 4 ), dimension ( test_num ) :: r4_test = (/ &
      0.25E+00,  0.5E+00,   1.0E+00,   2.0E+00,  4.0E+00, &
      1.5E+00,   1.75E+00,  1.875E+00, 6.5E+00, -6.5E+00, &
     99.0E+00, 100.0E+00, 101.0E+00,   0.0E+00, -1.0E+00, &
    huge ( 1.0E+00 ) /)
  integer ( kind = 4 ) s
  integer ( kind = 4 ) test
  integer ( kind = 4 ) word

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST002'
  write ( *, '(a)' ) '  B4_IEEE_TO_R4:  32 bit string => R4'
  write ( *, '(a)' ) '  R4_TO_B4_IEEE:  R4 => 32 bit string'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    R1               --------------Word--------------      R2'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    r1 = r4_test(test)
    call r4_to_b4_ieee ( r1, word )
    bits = i4_to_s32 ( word )
    call b4_ieee_to_r4 ( word, r2 )
    write ( *, '(g20.12,2x,a32,2x,g20.12)' ) r1, bits, r2
  end do
!
!  Extra test values, some of which are unnormalized real quantities.
!
   s = 0
   e = -125
   f = 3

   call sef_to_r4 ( s, e, f, r1 )
   call r4_to_b4_ieee ( r1, word )
   bits = i4_to_s32 ( word )
   call b4_ieee_to_r4 ( word, r2 )
   write ( *, '(g20.12,2x,a32,2x,g20.12)' ) r1, bits, r2

   s = 0
   e = -127
   f = 3

   call sef_to_r4 ( s, e, f, r1 )
   call r4_to_b4_ieee ( r1, word )
   bits = i4_to_s32 ( word )
   call b4_ieee_to_r4 ( word, r2 )
   write ( *, '(g20.12,2x,a32,2x,g20.12)' ) r1, bits, r2

   s = 0
   e = -129
   f = 3

   call sef_to_r4 ( s, e, f, r1 )
   call r4_to_b4_ieee ( r1, word )
   bits = i4_to_s32 ( word )
   call b4_ieee_to_r4 ( word, r2 )
   write ( *, '(g20.12,2x,a32,2x,g20.12)' ) r1, bits, r2

   s = 0
   e = -132
   f = 7

   call sef_to_r4 ( s, e, f, r1 )
   call r4_to_b4_ieee ( r1, word )
   bits = i4_to_s32 ( word )
   call b4_ieee_to_r4 ( word, r2 )
   write ( *, '(g20.12,2x,a32,2x,g20.12)' ) r1, bits, r2

   s = 0
   e = -135
   f = 15

   call sef_to_r4 ( s, e, f, r1 )
   call r4_to_b4_ieee ( r1, word )
   bits = i4_to_s32 ( word )
   call b4_ieee_to_r4 ( word, r2 )
   write ( *, '(g20.12,2x,a32,2x,g20.12)' ) r1, bits, r2

  return
end
subroutine test003 ( )

!*****************************************************************************80
!
!! TEST003 tests B4_IEEE_TO_SEF and SEF_TO_B4_IEEE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 26

  character ( len = 32 ) bits
  integer ( kind = 4 ) e
  integer ( kind = 4 ) e2
  integer ( kind = 4 ), parameter, dimension ( test_num) :: etest = (/ &
      -2,   -1,    0,    1,    2, &
      -1,   -2,   -3,   -1,   -1, &
       0,    2,    0,    0,    0, &
     104, -125, -127, -129, -132, &
    -135,    0,    0,  128,  128, &
     128 /)
  integer ( kind = 4 ) f
  integer ( kind = 4 ) f2
  integer ( kind = 4 ), parameter, dimension ( test_num) :: ftest = (/ &
           1,  1,   1,  1,  1, &
           3,  7,  15, 13, 13, &
          99, 25, 101,  0,  1, &
    16777215,  3,   3,  3,  7, &
          15,  0,   0,  1,  1, &
           0  /)
  character ( len = 32 ) i4_to_s32
  integer ( kind = 4 ) s
  integer ( kind = 4 ) s2
  integer ( kind = 4 ), parameter, dimension ( test_num) :: s_test = (/ &
    0, 0, 0, 0, 0, &
    0, 0, 0, 0, 1, &
    0, 0, 0, 0, 1, &
    0, 0, 0, 0, 0, &
    0, 0, 1, 0, 1, &
    0 /)
  integer ( kind = 4 ) test
  integer ( kind = 4 ) word

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST003'
  write ( *, '(a)' ) '  B4_IEEE_TO_SEF converts a real IEEE word to SEF form.'
  write ( *, '(a)' ) '  SEF_TO_B4_IEEE converts SEF form to a real IEEE word.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  S is the sign bit (0 = positive, 1 = negative)'
  write ( *, '(a)' ) '  E is the exponent base 2'
  write ( *, '(a)' ) '  F is the mantissa'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '   S    E         F  SEEEEEEEEFFFFFFFFFFFFFFFFFFFFFFF  S2   E2        F2'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    s = s_test(test)
    e = etest(test)
    f = ftest(test)

    call sef_to_b4_ieee ( s, e, f, word )
    bits = i4_to_s32 ( word )
    call b4_ieee_to_sef ( word, s2, e2, f2 )
    write ( *, '(2x,i2,i5,i10,2x,a32,2x,i2,i5,i10)' ) s, e, f, bits, s2, e2, f2

  end do

  return
end
subroutine test005 ( )

!*****************************************************************************80
!
!! TEST005 tests BASE_TO_I4 and I4_TO_BASE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 6

  integer ( kind = 4 ) base
  integer ( kind = 4 ), dimension ( test_num ) :: base_test = (/ &
    -1, 1, 2, 3, 4, 8 /)
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ), dimension ( test_num ) :: i4_test = (/ &
    5, 5, 21, -243, 16, 15 /)
  character ( len = 20 ) s
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST005'
  write ( *, '(a)' ) '  BASE_TO_I4 converts an integer in some other'
  write ( *, '(a)' ) '    base into base 10.'
  write ( *, '(a)' ) '  I4_TO_BASE converts an integer base 10 to '
  write ( *, '(a)' ) '    its representation in another base;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  BASE, I, I4_TO_BASE(I), BASE_TO_I4(I4_TO_BASE(I))'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    i1 = i4_test(test)
    base = base_test(test)

    call i4_to_base ( i1, base, s )
    call base_to_i4 ( s, base, i2 )

    write ( *, '(i8,2x,i8,2x,a,i8)' ) base, i1, s, i2

  end do

  return
end
subroutine test006 ( )

!*****************************************************************************80
!
!! TEST006 tests BINARY_TO_I4 and I4_TO_BINARY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 4

  integer ( kind = 4 ) i4
  integer ( kind = 4 ), dimension ( test_num ) :: i4_test = (/ 21, -32, 2, 128 /)
  integer ( kind = 4 ) j4
  character ( len = 10 ) s
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST006'
  write ( *, '(a)' ) '  BINARY_TO_I4 converts a binary to an integer.'
  write ( *, '(a)' ) '  I4_TO_BINARY converts an integer to binary,'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        I4  ==> BINARY    ==> I4'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    i4 = i4_test(test)
    call i4_to_binary ( i4, s )
    call binary_to_i4 ( s, j4 )

    write ( *, '(2x,i8,2x,a,2x,i8)' ) i4, s, j4

  end do

  return
end
subroutine test007 ( )

!*****************************************************************************80
!
!! TEST007 tests BINARY_TO_R4 and R4_TO_BINARY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 3

  real    ( kind = 4 ) r1
  real    ( kind = 4 ) r2
  real    ( kind = 4 ), dimension ( test_num ) :: r4_test = (/ &
    -10.75E+00, 0.4078125E+00, 0.666666E+00 /)
  character ( len = 20 ) s
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST007'
  write ( *, '(a)' ) '  BINARY_TO_R4: binary string => R4.'
  write ( *, '(a)' ) '  R4_TO_BINARY: R4 => binary string;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    R4   =>    S   =>   R4'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    r1 = r4_test(test)
    call r4_to_binary ( r1, s )
    call binary_to_r4 ( s, r2 )
    write ( *, '(f12.6, 2x, a, 2x, f12.6)' ) r1, s, r2
  end do

  return
end
subroutine test008 ( )

!*****************************************************************************80
!
!! TEST008 tests BINARY_TO_R8 and R8_TO_BINARY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 3

  real    ( kind = 8 ) r1
  real    ( kind = 8 ) r2
  real    ( kind = 8 ), dimension ( test_num ) :: r8_test = (/ &
    -10.75D+00, 0.4078125D+00, 0.666666D+00 /)
  character ( len = 20 ) s
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST008'
  write ( *, '(a)' ) '  BINARY_TO_R8: binary string => R8.'
  write ( *, '(a)' ) '  R8_TO_BINARY: R8 => binary string;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    R8   =>    S   =>   R8'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    r1 = r8_test(test)
    call r8_to_binary ( r1, s )
    call binary_to_r8 ( s, r2 )
    write ( *, '(f12.6, 2x, a, 2x, f12.6)' ) r1, s, r2
  end do

  return
end
subroutine test011 ( )

!*****************************************************************************80
!
!! TEST011 tests CH_CAP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: test_num = 5

  character c
  character, dimension ( test_num ) :: c_test = (/ &
    'F', 'f', '1', 'b', 'B' /)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST011'
  write ( *, '(a)' ) '  CH_CAP uppercases a character.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  C  CH_CAP(C)'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    c = c_test(test)
    call ch_cap ( c )
    write ( *, '(2x,a,2x,a)' ) c_test(test), c
  end do

  return
end
subroutine test012 ( )

!*****************************************************************************80
!
!! TEST012 tests CH_COUNT_FILE_ADD.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) count(0:255)
  character ( len = 80 ) :: file_name = 'chrpak_prb.f90'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST012'
  write ( *, '(a)' ) '  CH_COUNT_FILE_ADD adds the characters in a file'
  write ( *, '(a)' ) '  to a character count.'

  call ch_count_init ( count )

  call ch_count_file_add ( file_name, count )

  call ch_count_print ( count, 'Raw character count data:' )

  call ch_count_histogram_print ( count, file_name )

  return
end
subroutine test013 ( )

!*****************************************************************************80
!
!! TEST013 tests CH_EXTRACT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character c
  character ( len = 80 ) s

  s = '  A  bc $ '

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST013'
  write ( *, '(a)' ) '  CH_EXTRACT extracts characters from a string.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The string: "' // trim ( s ) // '"'
  write ( *, '(a)' ) ' '

  do

    call ch_extract ( s, c )

    if ( c == ' ' ) then
      exit
    end if

    write ( *, '(4x,a)' ) c

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Reached the last character.'

  return
end
subroutine test014 ( )

!*****************************************************************************80
!
!! TEST014 tests CH_INDEX_FIRST.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character              c
  integer   ( kind = 4 ) ch_index_first
  integer   ( kind = 4 ) iloc
  character ( len = 40 ) s

  c = 'g'
  s = 'Joel prefers graphics to graphs.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST014'
  write ( *, '(a)' ) '  CH_INDEX_FIRST searches a string for a character.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  String =    "' // trim ( s ) // '"'
  write ( *, '(a)' ) '  Character = ' // c

  iloc = ch_index_first ( s, c )

  write ( *, '(a,i8)' ) '  Character occurs at location ', iloc

  return
end
subroutine test015 ( )

!*****************************************************************************80
!
!! TEST015 tests CH_INDEX_LAST.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character c
  integer ( kind = 4 ) ch_index_last
  integer ( kind = 4 ) j
  character ( len = 40 ) s

  c = 'o'
  s = 'HELLO World, how ARE you?'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST015'
  write ( *, '(a)' ) '  CH_INDEX_LAST finds the LAST occurrence of a character.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  String =    "' // trim ( s ) // '"'
  write ( *, '(a)' ) '  Character = ' // c

  j = ch_index_last ( s, c )

  write ( *, '(a,i8)' ) '  Character occurs last at location ', j

  return
end
subroutine test016 ( )

!*****************************************************************************80
!
!! TEST016 tests CH_LOW.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: test_num = 5

  character c
  character, dimension ( test_num ) :: c_test = (/ &
    'F', 'f', '1', 'b', 'B' /)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST016'
  write ( *, '(a)' ) '  CH_LOW lowercases a character.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  C  CH_LOW(C)'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    c = c_test(test)
    call ch_low ( c )
    write ( *, '(2x,a,2x,a)' ) c_test(test), c
  end do

  return
end
subroutine test017 ( )

!*****************************************************************************80
!
!! TEST017 tests CH_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character c
  logical done
  character ( len = 20 ) s

  s = 'A B, C  DE  F'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST017'
  write ( *, '(a)' ) '  CH_NEXT returns characters from a string.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Input: "' // trim ( s ) // '"'

  done = .true.

  do

    call ch_next ( s, c, done )

    if ( done ) then
      write ( *, '(a)' ) '  No more characters.'
      exit
    end if

    write ( *, '(2x,a)' ) c

  end do

  return
end
subroutine test018 ( )

!*****************************************************************************80
!
!! TEST018 tests CH_ROMAN_TO_I4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) ch_roman_to_i4
  character c
  logical done
  integer ( kind = 4 ) ival
  character ( len = 20 ) s

  s = 'IJVXLCDMijvxlcdm0 W%'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST018'
  write ( *, '(a)' ) '  CH_ROMAN_TO_I4 converts a Roman numeral character'
  write ( *, '(a)' ) '  to its corresponding integer value.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Input: "' // trim ( s ) // '"'

  done = .true.

  do

    call ch_next ( s, c, done )

    if ( done ) then
      exit
    end if

    ival = ch_roman_to_i4 ( c )

    write ( *, '(2x,a,2x,i8)' ) c, ival

  end do

  return
end
subroutine test019 ( )

!*****************************************************************************80
!
!! TEST019 tests CH_TO_BRAILLE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) ncol2
  character ( len = 6 ) braille(3)
  character ( len = 12 ) :: s = 'SOS Titanic!'
  character ( len = 100 ) string2(3)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST019'
  write ( *, '(a)' ) '  CH_TO_BRAILLE converts a character to Braille.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here is the string to be converted:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  "' // trim ( s ) // '"'

  string2(1:3) = ' '

  ncol2 = 0

  do i = 1, len_trim ( s )

    call ch_to_braille ( s(i:i), ncol, braille )

    if ( 0 < ncol ) then

      do j = 1, 3
        string2(j)(ncol2+1:ncol2+ncol) = braille(j)(1:ncol)
      end do

      ncol2 = ncol2 + ncol

    end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Braille translation:'
  write ( *, '(a)' ) ' '

  do i = 1, 3
    write ( *, '(4x,a)' ) string2(i)(1:ncol2)
  end do

  return
end
subroutine test020 ( )

!*****************************************************************************80
!
!! TEST020 tests CH_TO_AMINO_NAME, CH_TO_CH3_AMINO, CH3_TO_CH_AMINO, I4_TO_AMINO_CODE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 27 ) amino_name
  character c
  character ch_back
  character ( len = 3 ) c3
  integer ( kind = 4 ) i
  character i4_to_a

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST020'
  write ( *, '(a)' ) '  CH_TO_CH3_AMINO converts a 1 character amino'
  write ( *, '(a)' ) '    acid code to 3 characters,'
  write ( *, '(a)' ) '  CH3_TO_CH_AMINO converts a 3 character amino'
  write ( *, '(a)' ) '    acid code to 1 character.'
  write ( *, '(a)' ) '  CH_TO_AMINO_NAME converts a 1 character amino'
  write ( *, '(a)' ) '    acid code to an amino acid name.'
  write ( *, '(a)' ) '  I4_TO_AMINO_CODE converts an integer to an'
  write ( *, '(a)' ) '    amino code.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   I -> A -> CCC -> C'
  write ( *, '(a)' ) ' '

  do i = 1, 26
    c = i4_to_a ( i )
    call ch_to_ch3_amino ( c, c3 )
    call ch3_to_ch_amino ( c3, ch_back )
    write ( *, '(2x,i2,4x,a1,4x,a3,4x,a1)' ) i, c, c3, ch_back
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   I -> Alpha -> AMINO_NAME'
  write ( *, '(a)' ) ' '

  do i = 1, 26
    c = i4_to_a ( i )
    call ch_to_amino_name ( c, amino_name )
    write ( *, '(2x,i2,4x,a1,4x,a27)' ) i, c, amino_name
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   I -> AMINO_CODE -> AMINO_NAME'
  write ( *, '(a)' ) ' '

  do i = 1, 23
    call i4_to_amino_code ( i, c )
    call ch_to_amino_name ( c, amino_name )
    write ( *, '(2x,i2,4x,a1,4x,a27)' ) i, c, amino_name
  end do

  return
end
subroutine test021 ( )

!*****************************************************************************80
!
!! TEST021 tests CH_TO_DIGIT and DIGIT_TO_CH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST021'
  write ( *, '(a)' ) '  CH_TO_DIGIT: character -> decimal digit'
  write ( *, '(a)' ) '  DIGIT_TO_C: decimal digit -> character.'
  write ( *, '(a)' ) ' '

  do i = -2, 11
    call digit_to_ch ( i, c )
    call ch_to_digit ( c, i2 )
    write ( *, '(2x,i8,a6,i8)' ) i, c, i2
  end do

  return
end
subroutine test022 ( )

!*****************************************************************************80
!
!! TEST022 tests HEX_DIGIT_TO_I4 and I4_TO_HEX_DIGIT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character              c
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) i2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST022'
  write ( *, '(a)' ) '  HEX_DIGIT_TO_I4: hexadecimal digit -> I4,'
  write ( *, '(a)' ) '  I4_TO_HEX_DIGIT: I4 -> hexadecimal digit.'
  write ( *, '(a)' ) ' '

  do i = -2, 17
    call i4_to_hex_digit ( i, c )
    call hex_digit_to_i4 ( c, i2 )
    write ( *, '(2x,i8,a6,i8)' ) i, c, i2
  end do

  return
end
subroutine test023 ( )

!*****************************************************************************80
!
!! TEST023 tests CH_TO_DIGIT_OCT and DIGIT_OCT_TO_C.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character              c
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) i2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST023'
  write ( *, '(a)' ) '  CH_TO_DIGIT_OCT: character -> hexadecimal'
  write ( *, '(a)' ) '  DIGIT_OCT_TO_C: hexadecimal -> character.'
  write ( *, '(a)' ) ' '

  do i = -2, 9
    call digit_oct_to_ch ( i, c )
    call ch_to_digit_oct ( c, i2 )
    write ( *, '(2x,i8,a6,i8)' ) i, c, i2
  end do

  return
end
subroutine test024 ( )

!*****************************************************************************80
!
!! TEST024 tests CH_TO_MILITARY and MILITARY_TO_CH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character c
  character ch_back
  character ( len = 8 ) c8
  integer ( kind = 4 ) i
  character i4_to_a

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST024'
  write ( *, '(a)' ) '  CH_TO_MILITARY converts a character to military code.'
  write ( *, '(a)' ) '  MILITARY_TO_CH converts a military code to a character.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   I -> C -> Mil -> C'
  write ( *, '(a)' ) ' '

  do i = 1, 52, 4
    c = i4_to_a ( i )
    call ch_to_military ( c, c8 )
    call military_to_ch ( c8, ch_back )
    write ( *, '(4x,i2,4x,a1,4x,a8,4x,a1)' ) i, c, c8, ch_back
  end do

  return
end
subroutine test025 ( )

!*****************************************************************************80
!
!! TEST025 tests CH_TO_MORSE and S_CAT1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  character ( len = 6 ) morse
  character ( len = 20 ) s
  character ( len = 80 ) s2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST025'
  write ( *, '(a)' ) '  CH_TO_MORSE converts ASCII to Morse.'
  write ( *, '(a)' ) '  S_CAT1 concatenates strings with a blank separator.'

  s = 'SOS Titanic!'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The string to be converted:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  "' // trim ( s ) // '"'

  s2 = ' '
  do i = 1, len_trim ( s )
    call ch_to_morse ( s(i:i), morse )
    call s_cat1 ( s2, morse, s2 )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Morse translation:'
  write ( *, '(a)' ) ' '
  write ( *, '(4x,a)' ) trim ( s2 )

  return
end
subroutine test026 ( )

!*****************************************************************************80
!
!! TEST026 tests CH_TO_ROT13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ch_to_rot13
  integer ( kind = 4 ) i
  character ( len = 80 ) s1
  integer ( kind = 4 ) s1_length
  character ( len = 80 ) s2
  character ( len = 80 ) s3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST026'
  write ( *, '(a)' ) '  CH_TO_ROT13 "encodes" a character using ROT13.'

  s1 = 'ABCDEFGHIJKLMNOPQRSTUVQXYZ'
  s1_length = len_trim ( s1 )

  s2 = ' '
  s3 = ' '

  do i = 1, s1_length
    s2(i:i) = ch_to_rot13 ( s1(i:i) )
    s3(i:i) = ch_to_rot13 ( s2(i:i) )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             CH  :' // trim ( s1 )
  write ( *, '(a)' ) '       ROT13(CH) :' // trim ( s2 )
  write ( *, '(a)' ) ' ROT13(ROT13(CH)):' // trim ( s3 )

  s1 = '  CH_TO_ROT13 "encodes" a character using ROT13.'
  s1_length = len_trim ( s1 )

  s2 = ' '
  s3 = ' '

  do i = 1, s1_length
    s2(i:i) = ch_to_rot13 ( s1(i:i) )
    s3(i:i) = ch_to_rot13 ( s2(i:i) )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             CH  :' // trim ( s1 )
  write ( *, '(a)' ) '       ROT13(CH) :' // trim ( s2 )
  write ( *, '(a)' ) ' ROT13(ROT13(CH)):' // trim ( s3 )

  return
end
subroutine test027 ( )

!*****************************************************************************80
!
!! TEST027 tests CH_TO_SOUNDEX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  character ( len = 30 ) s1
  character ( len = 30 ) s2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST027'
  write ( *, '(a)' ) '  CH_TO_SOUNDEX converts ASCII characters'
  write ( *, '(a)' ) '  to Soundex characters (digits).'

  s1 = 'SOS - Titanic & Mayflower!'
  s2 = ' '

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here is the string to be converted:'
  write ( *, '(a)' ) ' '
  write ( *, '(4x,a)' ) trim ( s1 )

  do i = 1, len_trim ( s1 )
    call ch_to_soundex ( s1(i:i), s2(i:i) )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Soundex translation:'
  write ( *, '(a)' ) ' '
  write ( *, '(4x,a)' ) trim ( s2 )

  return
end
subroutine test028 ( )

!*****************************************************************************80
!
!! TEST028 tests CH_TO_SYM and SYM_TO_CH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ch
  character ch2
  character ( len = 4 ) failok
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  logical ch_is_printable
  character ( len = 4 ) sym

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST028'
  write ( *, '(a)' ) '  CH_TO_SYM converts ANY charcter to a printable symbol.'
  write ( *, '(a)' ) '  SYM_TO_CH converts a printable symbol to a character.'
  write ( *, '(a)' ) ' '

  do i = 0, 255

    ch = char ( i )
    call ch_to_sym ( ch, sym )
    call sym_to_ch ( sym, ch2, ihi )

    if ( ch == ch2 ) then
      failok = 'OK'
    else
      failok = 'FAIL'
    end if

    if ( ch_is_printable ( ch ) ) then
      write ( *, '(2x,a4,2x,i3,2x,a1,4x,a4,4x,a1)' ) failok, i, ch, sym, ch2
    else
      write ( *, '(2x,a4,2x,i3,2x,1x,4x,a4,4x,1x)' ) failok, i,     sym
    end if

  end do

  return
end
subroutine test029 ( )

!*****************************************************************************80
!
!! TEST029 tests CH_UNIFORM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) a_to_i4
  character ch
  character ch_uniform
  character chi
  character clo
  integer ( kind = 4 ) count(26)
  integer ( kind = 4 ) i
  character i4_to_a
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST029'
  write ( *, '(a)' ) '  CH_UNIFORM returns a random character.'

  count(1:26) = 0

  clo = 'D'
  chi = 'W'
  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   I  A  Count'
  write ( *, '(a)' ) ' '

  do i = 1, 100000

    ch = ch_uniform ( clo, chi, seed )

    j = a_to_i4 ( ch )
    count(j) = count(j) + 1

  end do

  do i = 1, 26
    write ( *, '(2x,i2,2x,a1,2x,i5)' ) i, i4_to_a(i), count(i)
  end do

  return
end
subroutine test030 ( )

!*****************************************************************************80
!
!! TEST030 tests CH4_TO_I4 and I4_TO_CH4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 4

  integer ( kind = 4 ) intval
  integer ( kind = 4 ) test
  character ( len = 4 ), dimension ( test_num ) :: word = (/ &
    'Adam', &
    'Bill', &
    'Crow', &
    'Dave' /)
  character ( len = 4 ) word2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST030'
  write ( *, '(a)' ) '  I4_TO_CH4: Integer -> 4 characters;'
  write ( *, '(a)' ) '  CH4_TO_I4: 4 characters -> Integer.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  CH4 --> CH4_TO_I4(CH4) --> I4_TO_CH4(CH4_TO_I4(CH4))'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    call ch4_to_i4 ( word(test), intval )

    call i4_to_ch4 ( intval, word2 )

    write ( *, '(2x,a4,2x,i12,2x,a4)' ) word(test), intval, word2

  end do

  do test = 1, test_num
    call s_reverse ( word(test) )
  end do

  do test = 1, test_num

    call ch4_to_i4 ( word(test), intval )

    call i4_to_ch4 ( intval, word2 )

    write ( *, '(2x,a4,2x,i12,2x,a4)' ) word(test), intval, word2

  end do

  return
end
subroutine test031 ( )

!*****************************************************************************80
!
!! TEST031 tests CH4_TO_R4 and R4_TO_CH4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 4

  integer ( kind = 4 ) i
  real    ( kind = 4 ) rval
  integer ( kind = 4 ) test
  character ( len = 4 ), dimension ( test_num ) :: word = (/ &
    'Adam', &
    'Bill', &
    'Crow', &
    'Dave' /)
  character ( len = 4 ) word2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST031'
  write ( *, '(a)' ) '  CH4_TO_R4: 4 character => R4.'
  write ( *, '(a)' ) '  R4_TO_CH4: R4 => 4 character.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  word --> CH4_TO_R4(word) --> R4_TO_CH4(CH4_TO_R4(word))'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    call ch4_to_r4 ( word(test), rval )

    call r4_to_ch4 ( rval, word2 )

    write ( *, '(2x,a4,2x,g14.6,2x,a4)' ) word(test), rval, word2

  end do

  return
end
subroutine test032 ( )

!*****************************************************************************80
!
!! TEST032 tests CH4VEC_TO_I4VEC and I4VEC_TO_CH4VEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer   ( kind = 4 ), parameter :: n = 11

  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) i4vec(n)
  integer   ( kind = 4 ) i4vec2(n)
  character ( len = 4 * n ) s

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST032'
  write ( *, '(a)') '   For vectors of integers and character*4 strings:'
  write ( *, '(a)' ) '  CH4VEC_TO_I4VEC: CH4 => I.'
  write ( *, '(a)' ) '  I4VEC_TO_CH4VEC: I => CH4.'

  do i = 1, n
    i4vec(i) = i - ( n / 2 )
  end do

  call i4vec_to_ch4vec ( n, i4vec, s )

  call ch4vec_to_i4vec ( n, s, i4vec2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I  Input  Output'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i3,2x,i8,2x,i8)' ) i, i4vec(i), i4vec2(i)
  end do

  return
end
subroutine test033

!*****************************************************************************80
!
!! TEST033 tests CH4VEC_TO_I4VEC and I4VEC_TO_CH4VEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  integer ( kind = 4 ) i4vec(n)
  character ( len = 4*n ) s
  character ( len = 4*n ) s2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST033'
  write ( *, '(a)' ) '   For vectors of integers and character*4 strings:'
  write ( *, '(a)' ) '  CH4VEC_TO_I4VEC: CH4 => I4.'
  write ( *, '(a)' ) '  I4VEC_TO_CH4VEC: I4 => CH4.'

  s = 'Bartleby !'

  call ch4vec_to_i4vec ( n, s, i4vec )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Input string:  ' // trim ( s(1:4*n) )

  call i4vec_print ( n, i4vec, '  Integer vector:' )

  call i4vec_to_ch4vec ( n, i4vec, s2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Output string: ' // trim ( s2(1:4*n) )

  return
end
subroutine test034 ( )

!*****************************************************************************80
!
!! TEST034 tests CHR4_TO_8 and CHR8_TO_4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character chrtmp
  character chrtmp2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ichr
  integer ( kind = 4 ) j
  character ( len = 256 ) s1
  character ( len = 512 ) s2
  character ( len = 256 ) s3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST034'
  write ( *, '(a)' ) '  CHR8_TO_4 convert characters to pairs of hexadecimals.'
  write ( *, '(a)' ) '  CHR4_TO_8 converts pairs of hexadecimals to characters.'
  write ( *, '(a)' ) ' '

  do i = 1, 256
    s1(i:i) = char(i-1)
  end do

  call chr8_to_4 ( s1, s2 )

  call chr4_to_8 ( s2, s3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Coded characters that can''t be printed are shown as blanks.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   ASCII  Coded  Decoded'
  write ( *, '(a)' ) ' '

  do i = 1, 256

    ichr = i - 1
    j = 2 * i - 1

    if ( 33 <= ichr .and. ichr <= 127 ) then
      chrtmp = s1(i:i)
      chrtmp2 = s3(i:i)
    else
      chrtmp = ' '
      chrtmp2 = ' '
    end if

    write ( *, '(2x,i3,1x,a1,6x,a2,7x,a1)' ) ichr, chrtmp, s2(j:j+1), chrtmp2

  end do

  return
end
subroutine test035 ( )

!*****************************************************************************80
!
!! TEST035 tests CHRASS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 8

  character ( len = 20 ) lhs
  character ( len = 20 ) rhs
  character ( len = 20 ) s
  character ( len = 20 ), dimension ( test_num ) :: s_test = (/ &
    'a = 1.0             ', &
    'n = -17             ', &
    'scale = +5.3E-2     ', &
    'filename = myprog.f ', &
    ' =  A pot of gold   ', &
    'Fred                ', &
    ' =  Bob             ', &
    '1 = 2, 2 = 3, 3 = 4 ' /)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST035'
  write ( *, '(a)' ) '  CHRASS parses an assignment statement.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  S                  LHS(S)                   RHS(S)'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    s = s_test(test)
    call chrass ( s, lhs, rhs )
    write ( *, '(2x,a20,2x,a20,2x,a20)' ) s, lhs, rhs
  end do

  return
end
subroutine test036 ( )

!*****************************************************************************80
!
!! TEST036 tests CHRCTP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 10

  complex cval
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) length
  character ( len = 20 ) string(test_num)
  integer ( kind = 4 ) test

  string ( 1) = '(1,1)'
  string ( 2) = '(,)'
  string ( 3) = '( 20 , 99 )'
  string ( 4) = '(-1.2E+2, +30E-2)'
  string ( 5) = '(1)'
  string ( 6) = '(1,2,3)'
  string ( 7) = '(4,5('
  string ( 8) = '(6,)'
  string ( 9) = '(7;8)'
  string (10) = '9'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST036'
  write ( *, '(a)' ) '  CHRCTP accepts a string of characters'
  write ( *, '(a)' ) '  and extracts a complex value from them,'
  write ( *, '(a)' ) '  assuming the format (A,B) for complex numbers.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  STRING                        CVAL    IERROR   LENGTH'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    call chrctp ( string(test), cval, ierror, length )

    write ( *, '(2x,a20,2x,2f8.1,2x,i2,6x,i2)' ) &
      string(test), cval, ierror, length

  end do

  return
end
subroutine test037 ( )

!*****************************************************************************80
!
!! TEST037 tests CHVEC_PERMUTE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  character chvec(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST037'
  write ( *, '(a)' ) '  CHVEC_PERMUTE permutes a character vector.'

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Using random number seed SEED = ', seed

  call perm_uniform ( n, seed, p )

  call i4vec_print ( n, p, '  The random permutation:' )

  do i = 1, n
    chvec(i) = char ( ichar ( 'A' ) + i - 1 )
  end do

  call chvec_print ( n, chvec, '  CHVEC before permutation:' )

  call chvec_permute ( n, chvec, p )

  call chvec_print ( n, chvec, '  CHVEC after permutation:' )

  return
end
subroutine test038 ( )

!*****************************************************************************80
!
!! TEST038 tests CHVEC_TO_S and S_TO_CHVEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 20

  character chvec(20)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  character ( len = 20 ) s

  s = 'Yabba Blabba'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST038'
  write ( *, '(a)' ) '  CHVEC_TO_S: character vector => string;'
  write ( *, '(a)' ) '  S_TO_CHVEC: string to character vector.'
  write ( *, '(a)' ) ' '

  n = 0
  call s_to_chvec ( s, n, chvec )

  write ( *, '(a)' ) '  String: ' // trim ( s )
  write ( *, '(a)' ) ' '
  write ( *, '(a,20(1x,a1))' ) '  CHVEC: ', ( chvec(i), i = 1, n )

  s = ' '

  call chvec_to_s ( n, chvec, s )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Recovered string: "' // trim ( s ) // '"'

  return
end
subroutine test039 ( )

!*****************************************************************************80
!
!! TEST039 tests COMMA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 30 ) input
  character ( len = 30 ) output

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST039'
  write ( *, '(a)' ) '  COMMA shifts commas left through blanks.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  --------Input-------  -------Output-------'
  write ( *, '(a)' ) ' '

  input = ' To Henry , our dog ,'
  output = input
  call comma ( output )
  write ( *, '(2x,a,2x,a)' ) input, output

  input = ' 14 , 15  , 16   ,'
  output = input
  call comma ( output )
  write ( *, '(2x,a,2x,a)' ) input, output

  input = ' ,  ,  , '
  output = input
  call comma ( output )
  write ( *, '(2x,a,2x,a)' ) input, output

  return
end
subroutine test040 ( )

!*****************************************************************************80
!
!! TEST040 tests DEC_TO_S_LEFT and S_TO_DEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 20

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ), dimension ( test_num ) :: itest = (/ &
       0,   21,   -3,  -31,  147,   16,   34,  123,  123,  123, &
     123,  123, -123, -123, -123, -123, -123,   34,   99,   99 /)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ), dimension ( test_num ) :: jtest = (/ &
      0,   3,   0, -2,   -2,  -5,  30, -19, -20, -21, &
    -22, -23, -19, -20, -21, -22, -23, -30, -99,  99 /)
  integer ( kind = 4 ) length
  character ( len = 22 ) s
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST040'
  write ( *, '(a)' ) '  For decimals I * 10**J,'
  write ( *, '(a)' ) '  DEC_TO_S_LEFT: -> decimal to left string;'
  write ( *, '(a)' ) '  S_TO_DEC: string to decimal.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   I         J                 S_LEFT  ' &
       // '      LENGTH     I2     J2'
  write ( *, '(a)' ) '--------- ---------  ' // &
       '----------------------  ------  --------------'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    i = itest(test)
    j = jtest(test)

    call dec_to_s_left ( i, j, s )
    call s_to_dec ( s, i2, j2, length )

    write ( *, '(2x,i10,i10,2x,a22,2x,i3,2x,i10,i10)' ) i, j, s, length, i2, j2

  end do

  return
end
subroutine test041 ( )

!*****************************************************************************80
!
!! TEST041 tests DEC_TO_S_RIGHT and S_TO_DEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 20

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ), dimension ( test_num ) :: itest = (/ &
       0,   21,   -3,  -31,  147,   16,   34,  123,  123,  123, &
     123,  123, -123, -123, -123, -123, -123,   34,   99,   99 /)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ), dimension ( test_num ) :: jtest = (/ &
      0,   3,   0, -2,   -2,  -5,  30, -19, -20, -21, &
    -22, -23, -19, -20, -21, -22, -23, -30, -99,  99 /)
  integer ( kind = 4 ) length
  character ( len = 22 ) s
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST041'
  write ( *, '(a)' ) '  For decimals I * 10**J,'
  write ( *, '(a)' ) '  DEC_TO_S_RIGHT: -> decimal to right string.'
  write ( *, '(a)' ) '  S_TO_DEC: string to decimal.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   I         J                S_RIGHT  ' &
       // '      LENGTH     I2     J2'
  write ( *, '(a)' ) '--------- ---------  ' // &
       '----------------------  ------  --------------'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    i = itest(test)
    j = jtest(test)

    call dec_to_s_right ( i, j, s )
    call s_to_dec ( s, i2, j2, length )

    write ( *, '(2x,i10,i10,2x,a22,2x,i3,2x,i10,i10)' ) i, j, s, length, i2, j2

  end do

  return
end
subroutine test042 ( )

!*****************************************************************************80
!
!! TEST042 tests DEC_TO_S_LEFT and S_TO_DEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 11

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) length
  character ( len = 10 ) s
  character ( len = 22 ) s2
  character ( len = 10 ), dimension ( test_num ) :: s_test = (/ &
    '1         ', '1A        ', '+12,34,56 ', '  34 7    ', &
    '-1 E2ABCD ', '-1 X2ABCD ', ' 2E-1     ', '23.45     ', &
    'Inf       ', 'NaN       ', '  &#99    ' /)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST042'
  write ( *, '(a)' ) '  For decimals I * 10**J,'
  write ( *, '(a)' ) '  DEC_TO_S_LEFT: -> decimal to left string;'
  write ( *, '(a)' ) '  S_TO_DEC: string to decimal.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          S in               I       J' &
    // '  LENGTH           S out       '
  write ( *, '(a)' ) '----------------------  ------  ------' &
    // '  ------  ---------------------'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    s = s_test(test)

    call s_to_dec ( s, i, j, length )
    call dec_to_s_left ( i, j, s2 )

    write ( *, '(2x,a,2x,i8,2x,i8,2x,i8,2x,a)' ) s, i, j, length, s2

  end do

  return
end
subroutine test043 ( )

!*****************************************************************************80
!
!! TEST043 tests EBCDIC_TO_S.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 13 ) s

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST043'
  write ( *, '(a)' ) '  EBCDIC_TO_S converts a EBCDIC string to ASCII.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We will not print out the "before" string!'
  write ( *, '(a)' ) ' '

  s = char(200) // char(133) // char(147) // char(147) // char(150) // &
      char(107) // char( 64) // char(166) // char(150) // char(153) // &
      char(147) // char(132) // char( 90)

  call ebcdic_to_s ( s )

  write ( *, '(a)' ) '  After conversion:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  "' // trim ( s ) // '"'

  return
end
subroutine test045 ( )

!*****************************************************************************80
!
!! TEST045 tests FLT_TO_S and R4_TO_FLT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 September 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 10

  integer ( kind = 4 ) iexp
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) mant
  integer ( kind = 4 ) ndig
  real    ( kind = 4 ), dimension ( test_num ) :: r4_test = (/ &
    1.0E+00, 10.0E+00, 100.0E+00, 101.0E+00, 99.0E+00, &
    0.0E+00, -1.0E+00, -123.456E+00, -0.123456E+00, 0.000000123456E+00 /)
  real    ( kind = 4 ) rval
  character ( len = 40 ) s
  integer ( kind = 4 ) test

  ndig = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST045'
  write ( *, '(a)' ) '  R4_TO_FLT: real -> scientific representation;'
  write ( *, '(a)' ) '  FLT_TO_S: scientific representation -> string:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of digits used is ', ndig
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      RVAL      ISGN  MANT  IEXP  S'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    rval = r4_test(test)
    call r4_to_flt ( rval, isgn, mant, iexp, ndig )
    mant = isgn * mant

    call flt_to_s ( mant, iexp, ndig, s )

    write ( *, '(2x,g14.6,2x,i2,2x,i8,2x,i8,2x,a40)' ) rval, isgn, mant, iexp, s

  end do

  return
end
subroutine test046 ( )

!*****************************************************************************80
!
!! TEST046 tests HEX_TO_I4 and I4_TO_HEX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 3

  character ( len = 8 ) hex
  integer ( kind = 4 ) i4
  integer ( kind = 4 ), dimension (test_num) :: i4_test = (/ 21, -32, 1776 /)
  integer ( kind = 4 ) j4
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST046'
  write ( *, '(a)' ) '  HEX_TO_I4, hexadecimal->integer.'
  write ( *, '(a)' ) '  I4_TO_HEX, integer->hexadecimal'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I  I4_TO_HEX(I)  HEX_TO_I4(I4_TO_HEX(I)) '
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    i4 = i4_test(test)
    call i4_to_hex ( i4, hex )
    call hex_to_i4 ( hex, j4 )

    write ( *, '(2x,i8,2x,a8,2x,i8)' ) i4, hex, j4

  end do

  return
end
subroutine test047 ( )

!*****************************************************************************80
!
!! TEST047 tests HEX_TO_S and S_TO_HEX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 3

  character ( len = 5 ) chrval(test_num)
  character ( len = 5 ) chrval2
  character ( len = 10 ) hexstr
  integer ( kind = 4 ) test

  chrval(1) = 'ABC'
  chrval(2) = 'Wow!!'
  chrval(3) = '1234'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST047'
  write ( *, '(a)' ) '  S_TO_HEX: string -> hexadecimal;'
  write ( *, '(a)' ) '  HEX_TO_S: hexadecimal -> string.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  String  Hexadecimal  Recovered string'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    call s_to_hex ( chrval(test), hexstr )
    call hex_to_s ( hexstr, chrval2 )

    write ( *, '(2x,a5,2x,a10,2x,a5)' ) chrval(test), hexstr, chrval2

  end do

  return
end
subroutine test048 ( )

!*****************************************************************************80
!
!! TEST048 tests I4_BYTE_SWAP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  real    ( kind = 4 ) pi
  real    ( kind = 4 ) temp
  real    ( kind = 4 ) x(n)

  pi = 4.0E+00 * atan2 ( 1.0E+00, 1.0E+00 )

  temp = 1.0E+00

  do i = 1, n
    temp = - pi * temp
    x(i) = temp
  end do
!
!  Tell the user our data.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST048'
  write ( *, '(a)' ) '  I4_BYTE_SWAP swaps bytes in a 4 byte word.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data from a different computer can be'
  write ( *, '(a)' ) '  read this way, if necessary.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here is the data written to the file:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(g14.6)' ) x(i)
  end do
!
!  Write the data to a fixed length record file.
!
  open ( unit = 1, file = 'chrprb.dat', form = 'unformatted', &
    access = 'direct', recl = 4, iostat = ios, status = 'replace' )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  Error while opening unit 1.'
    stop
  end if

  do i = 1, n
    write ( 1, rec = i ) x(i)
  end do

  close ( unit = 1 )

  return
end
subroutine test049 ( )

!*****************************************************************************80
!
!! TEST049 tests I4_BYTE_SWAP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) bytes(4)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  real    ( kind = 4 ) temp
  real    ( kind = 4 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST049'
  write ( *, '(a)' ) '  I4_BYTE_SWAP swaps bytes.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read the data in CHRPRB.DAT.'
!
!  Read the data from a fixed length record file.
!
  open ( unit = 1, file = 'chrprb.dat', form = 'unformatted', &
    access = 'direct', recl = 4, iostat = ios, status = 'old' )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  Error while opening unit 1.'
    stop
  end if

  do i = 1, n
    read ( 1, rec = i ) x(i)
  end do

  close ( unit = 1, status = 'delete' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Here is the plain data from the file:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(g14.6)' ) x(i)
  end do

  bytes = (/ 4, 3, 2, 1 /)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Using byte order:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,4i1)' ) bytes
  write ( *, '(a)' ) '  our data becomes:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    temp = x(i)
    call i4_byte_swap ( temp, bytes )
    write ( *, '(g14.6)' ) temp
  end do

  bytes = (/ 2, 1, 4, 3 /)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Using byte order:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,4i1)' ) bytes
  write ( *, '(a)' ) '  our data becomes:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    temp = x(i)
    call i4_byte_swap ( temp, bytes )
    write ( *, '(g14.6)' ) temp
  end do

  bytes = (/ 3, 4, 1, 2 /)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Using byte order:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,4i1)' ) bytes
  write ( *, '(a)' ) '  our data becomes:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    temp = x(i)
    call i4_byte_swap ( temp, bytes )
    write ( *, '(g14.6)' ) temp
  end do

  bytes = (/ 2, 2, 2, 4 /)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Using byte order:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,4i1)' ) bytes
  write ( *, '(a)' ) '  our data becomes:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    temp = x(i)
    call i4_byte_swap ( temp, bytes )
    write ( *, '(g14.6)' ) temp
  end do

  return
end
subroutine test050 ( )

!*****************************************************************************80
!
!! TEST050 tests I4_EXTRACT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i4
  integer ( kind = 4 ) ierror
  character ( len = 80 ) s

  s = '  123    45   789'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST050'
  write ( *, '(a)' ) '  I4_EXTRACT extracts integers from a string.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  "' // trim ( s ) // '"'
  write ( *, '(a)' ) ' '

  do

    call i4_extract ( s, i4, ierror )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Reached the last integer.'
      exit
    end if

    write ( *, '(2x,i8)' ) i4

  end do

  return
end
subroutine test051 ( )

!*****************************************************************************80
!
!! TEST051 tests I4_LENGTH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 6

  integer ( kind = 4 ) i4
  integer ( kind = 4 ) i4_length
  integer ( kind = 4 ), dimension ( test_num ) :: i4_test = (/ &
    0, 1, -1, 140, -1952, 123456 /)
  integer ( kind = 4 ) j4
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST051'
  write ( *, '(a)' ) '  I4_LENGTH computes an integer''s "length".'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        I4    Length'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    i4 = i4_test(test)

    j4 = i4_length ( i4_test(test) )

    write ( *, '(2x,i8,2x,i8)' ) i4, j4

  end do

  return
end
subroutine test052 ( )

!*****************************************************************************80
!
!! TEST052 tests I4_NEXT_READ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) intval
  character ( len = 80 ) s

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST052'
  write ( *, '(a)' ) '  I4_NEXT_READ extracts integers from a string.'

  s = 'Data set #12 extends from (5,-43) and is worth $4.56'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  String to be analyzed:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  "' // trim ( s ) // '"'

  ierror = -1
  i = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  #  Integer'
  write ( *, '(a)' ) ' '

  do

    call i4_next_read ( s, intval, ierror )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of integers found was ', i
      exit
    end if

    i = i + 1
    write ( *, '(2x,i3,2x,i10)' ) i, intval

    if ( 99 <= i ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Fatal error!'
      write ( *, '(a)' ) '  Reading phantom data from string.'
      stop
    end if

  end do

  return
end
subroutine test054 ( )

!*****************************************************************************80
!
!! TEST054 tests I4_TO_BINHEX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  character i4_to_binhex
  character ( len = 64 ) s

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST054'
  write ( *, '(a)' ) '  I4_TO_BINHEX: I => BINHEX character'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The BINHEX alphabet'
  write ( *, '(a)' ) ' '

  do i = 1, 64
    s(i:i) = i4_to_binhex ( i )
  end do

  write ( *, '(2x,a)' ) s

  return
end
subroutine test055 ( )

!*****************************************************************************80
!
!! TEST055 tests I4_TO_MONTH_NAME and MONTH_NAME_TO_I4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 9

  integer ( kind = 4 ) month
  character ( len = 9 ) month_name
  character ( len = 9 ) name_test(test_num)
  integer ( kind = 4 ) test

  name_test(1) = 'J'
  name_test(2) = 'Febooary'
  name_test(3) = 'Dec.'
  name_test(4) = 'April'
  name_test(5) = 'Aug'
  name_test(6) = 'Mar'
  name_test(7) = 'May'
  name_test(8) = 'o'
  name_test(9) = 'nO'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST055'
  write ( *, '(a)' ) '  I4_TO_MONTH_NAME: I => Month_Name'
  write ( *, '(a)' ) '  MONTH_NAME_TO_I4: Month_Name => I.'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    call month_name_to_i4 ( name_test(test), month )
    call i4_to_month_name ( month, month_name )
    write ( *, '(2x,a3,2x,i2,2x,a9)' ) name_test(test), month, month_name
  end do

  return
end
subroutine test056 ( )

!*****************************************************************************80
!
!! TEST056 tests I4_TO_NUNARY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 3

  integer i4
  integer ( kind = 4 ), dimension ( test_num ) :: i4_test = (/ -5, 0, 7 /)
  character ( len = 20 ) s
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST056'
  write ( *, '(a)' ) '  I4_TO_NUNARY converts an integer to negative unary.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      I4    NUNARY'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    i4 = i4_test(test)

    call i4_to_nunary ( i4, s )

    write ( *, '(2x,i8,2x,a)' ) i4, s

  end do

  return
end
subroutine test057 ( )

!*****************************************************************************80
!
!! TEST057 tests I4_TO_OCT and OCT_TO_I4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 3

  integer ( kind = 4 ) i4
  integer ( kind = 4 ), dimension ( test_num ) :: i4_test = (/ 21, -32, 1776 /)
  integer ( kind = 4 ) j4
  character ( len = 10 ) s
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST057'
  write ( *, '(a)' ) '  I4_TO_OCT, integer->octal'
  write ( *, '(a)' ) '  OCT_TO_I4, octal->integer.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      I4  ==>  OCT   ==>  I4'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    i4 = i4_test(test)
    call i4_to_oct ( i4, s )
    call oct_to_i4 ( s, j4 )

    write ( *, '(2x,i8,2x,a10,2x,i8)' ) i4, s, j4

  end do

  return
end
subroutine test058 ( )

!*****************************************************************************80
!
!! TEST058 tests I4_TO_S_LEFT and S_TO_I4;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 4

  integer ( kind = 4 ) i2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) length
  character ( len = 20 ) s1
  character ( len = 20 ) s2
  character ( len = 20 ), dimension ( test_num ) :: s_test = (/ &
    ' -124 56 AbC        ', &
    '25,50,5             ', &
    '+15.9               ', &
    '  123abc            ' /)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST058'
  write ( *, '(a)' ) '  I4_TO_S_LEFT: I4 -> left-justified string;'
  write ( *, '(a)' ) '  S_TO_I4:      string->I4.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    STRING ==> S_TO_I4 ==> I4_TO_S_LEFT'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    s1 = s_test(test)

    call s_to_i4 ( s1, i2, ierror, length )

    call i4_to_s_left ( i2, s2 )

    write ( *, '(2x,a,2x,i8,2x,a)' ) s1, i2, s2

  end do

  return
end
subroutine test059 ( )

!*****************************************************************************80
!
!! TEST059 tests I4_TO_S_LEFT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 7

  integer ( kind = 4 ) i4
  integer ( kind = 4 ), dimension ( test_num) :: i4_test = (/ &
    0, 1, -1, 140, -1952, 123456, 1234567 /)
  character ( len = 6 ) s
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST059'
  write ( *, '(a)' ) '  I4_TO_S_LEFT:  I4 -> Left-justified string;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I4       S_LEFT'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    i4 = i4_test(test)

    call i4_to_s_left ( i4, s )

    write ( *, '(2x,i8,2x,a)' ) i4, '"' // s // '"'

  end do

  return
end
subroutine test060 ( )

!*****************************************************************************80
!
!! TEST060 tests I4_TO_S_RIGHT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 7

  integer ( kind = 4 ) i4
  integer ( kind = 4 ), dimension ( test_num) :: i4_test = (/ &
    0, 1, -1, 140, -1952, 123456, 1234567 /)
  character ( len = 6 ) s
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST060'
  write ( *, '(a)' ) '  I4_TO_S_RIGHT:  I4 -> Right-justified string;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I4       S_RIGHT'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    i4 = i4_test(test)

    call i4_to_s_right ( i4, s )

    write ( *, '(2x,i8,2x,a)' ) i4, '"' // s // '"'

  end do

  return
end
subroutine test061 ( )

!*****************************************************************************80
!
!! TEST061 tests I4_TO_S_RIGHT_COMMA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 9

  integer ( kind = 4 ) i4
  integer ( kind = 4 ), dimension ( test_num) :: i4_test = (/ &
    0, 1, -1, 140, -1952, 123456, 1234567, 123456789, 1234567890 /)
  character ( len = 15 ) s
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST061'
  write ( *, '(a)' ) '  I4_TO_S_RIGHT_COMMA:'
  write ( *, '(a)' ) '  I4 -> Right-justified string with commas;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I4       S_RIGHT_COMMA'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    i4 = i4_test(test)

    call i4_to_s_right_comma ( i4, s )

    write ( *, '(2x,i10,2x,a)' ) i4, '"' // s // '"'

  end do

  return
end
subroutine test062 ( )

!*****************************************************************************80
!
!! TEST062 tests I4_TO_S_ROMAN and S_ROMAN_TO_I4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 5

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ), dimension ( test_num ) :: i_test = (/ 99, 157, 486, 1999, 4999 /)
  character ( len = 20 ) s
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST062'
  write ( *, '(a)' ) '  I4_TO_S_ROMAN: Integer -> Roman Numerals'
  write ( *, '(a)' ) '  S_ROMAN_TO_I4: Roman Numerals -> Integer.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    I4 ==> S'
  write ( *, '(a)' ) ' '

  do test = -5, 10
    i = test
    call i4_to_s_roman ( i, s )
    call s_roman_to_i4 ( s, i2 )
    write ( *, '(2x,i8,2x,a,2x,i8)' ) i, s, i2
  end do

  do test = 1, test_num
    i = i_test(test)
    call i4_to_s_roman ( i, s )
    call s_roman_to_i4 ( s, i2 )
    write ( *, '(2x,i8,2x,a,2x,i8)' ) i, s, i2
  end do

  return
end
subroutine test063 ( )

!*****************************************************************************80
!
!! TEST063 tests I4_TO_S_ZERO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 7

  integer ( kind = 4 ) i4
  integer ( kind = 4 ), dimension ( test_num) :: i4_test = (/ &
    0, 1, -1, 140, -1952, 123456, 1234567 /)
  character ( len = 6 ) s
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST063'
  write ( *, '(a)' ) '  I4_TO_S_ZERO:  I4 -> Zero-padded string;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I4       S_ZERO'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    i4 = i4_test(test)

    call i4_to_s_zero ( i4, s )

    write ( *, '(2x,i8,2x,a)' ) i4, '"' // s // '"'

  end do

  return
end
subroutine test064 ( )

!*****************************************************************************80
!
!! TEST064 tests I4_TO_S32 and S32_TO_I4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 4

  integer ( kind = 4 ) i4
  integer ( kind = 4 ), dimension ( test_num ) :: i4_test = (/ &
    0, 1, -1, 15 /)
  character ( len = 32 ) i4_to_s32
  integer ( kind = 4 ) j4
  character ( len = 32 ) s32
  integer ( kind = 4 ) s32_to_i4
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST064'
  write ( *, '(a)' ) '  I4_TO_S32: integer => character ( len = 32 );'
  write ( *, '(a)' ) '  S32_TO_I4: character ( len = 32 ) => integer.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '            I4  ---------------S32--------------            I4'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    i4 = i4_test(test)
    s32 = i4_to_s32 ( i4 )
    j4 = s32_to_i4 ( s32 )

    write ( *, '( 2x, i12, 2x, a32, 2x, i12 )' ) i4, s32, j4

  end do

  return
end
subroutine test065 ( )

!*****************************************************************************80
!
!! TEST065 tests I4_TO_UNARY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 3

  integer ( kind = 4 ) i4
  integer ( kind = 4 ), dimension ( test_num ) :: i4_test = (/ -5, 0, 7 /)
  character ( len = 10 ) s
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST065'
  write ( *, '(a)' ) '  I4_TO_UNARY converts an integer to unary.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        I4     UNARY'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    i4 = i4_test(test)

    call i4_to_unary ( i4, s )

    write ( *, '(2x,i8,2x,a)' ) i4, s

  end do

  return
end
subroutine test066 ( )

!*****************************************************************************80
!
!! TEST066 tests I4_TO_UUDECODE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  character i4_to_uudecode
  character ( len = 64 ) s

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST066'
  write ( *, '(a)' ) '  I4_TO_UUDECODE: I => UUDECODE character'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The UUDECODE alphabet'
  write ( *, '(a)' ) ' '

  do i = 1, 64
    s(i:i) = i4_to_uudecode ( i )
  end do

  write ( *, '(2x,a)' ) s

  return
end
subroutine test067 ( )

!*****************************************************************************80
!
!! TEST067 tests I4_TO_XXDECODE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
  integer ( kind = 4 ) i
  character i4_to_xxdecode
  character ( len = 64 ) s

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST067'
  write ( *, '(a)' ) '  I4_TO_XXDECODE: I => XXDECODE character'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The XXDECODE alphabet'
  write ( *, '(a)' ) ' '

  do i = 1, 64
    s(i:i) = i4_to_xxdecode ( i )
  end do

  write ( *, '(2x,a)' ) s

  return
end
subroutine test068 ( )

!*****************************************************************************80
!
!! TEST068 tests ISTRCMP and ISTRNCMP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 5

  integer ( kind = 4 ) istrcmp
  integer ( kind = 4 ) istrncmp
  integer ( kind = 4 ) itemp1
  integer ( kind = 4 ) itemp2
  integer ( kind = 4 ) nchar
  character ( len = 15 ) s1(test_num)
  character ( len = 15 ) s2(test_num)
  integer ( kind = 4 ) test

  nchar = 5

  s1(1) = 'Alex'
  s1(2) = 'Barney'
  s1(3) = 'Cray YMP'
  s1(4) = 'ZULU'
  s1(5) = 'BeHanna'

  s2(1) = 'Alexander'
  s2(2) = 'Babushka'
  s2(3) = 'Zulu'
  s2(4) = 'Zulu'
  s2(5) = 'BeHanna'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST068'
  write ( *, '(a)' ) '  ISTRCMP, C-like string comparison.'
  write ( *, '(a)' ) '  ISTRNCMP, C-like string comparisons.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'String1        String2      ISTRNCMP          ISTRCMP'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    itemp1 = istrncmp ( s1(test), s2(test), nchar )
    itemp2 = istrcmp ( s1(test), s2(test) )
    write ( *, '(2x,a,2x,a,2x,i2,2x,i2)' ) s1(test), s2(test), itemp1, itemp2
  end do

  return
end
subroutine test070 ( )

!*****************************************************************************80
!
!! TEST070 tests NAMEFL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 4

  character ( len = 30 ), dimension ( test_num ) :: name_test = (/ &
    'Brown, Charlie                ', &
    'Cher                          ', &
    'Howell, James Thurston        ', &
    'Shakespeare Joe Bob           ' /)
  character ( len = 30 ) s
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST070'
  write ( *, '(a)' ) '  NAMEFL takes a name in the '
  write ( *, '(a)' ) '  last name, first name order and restores the'
  write ( *, '(a)' ) '  first name, last name order.'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    s = name_test(test)
    call namefl ( s )
    write ( *, '(2x,a30,2x,a30)' ) name_test(test), s
  end do

  return
end
subroutine test071 ( )

!*****************************************************************************80
!
!! TEST071 tests NAMELF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 3

  character ( len = 30 ), dimension ( test_num ) :: s_test = (/ &
    'Charlie Brown                 ', &
    'Cher                          ', &
    'James Thurston Howell         ' /)
  character ( len = 30 ) s
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST071'
  write ( *, '(a)' ) '  NAMELF moves a last name first.'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    s = s_test(test)
    call namelf ( s )
    write ( *, '(2x,a30,2x,a30)' ) s_test(test), s
  end do

  return
end
subroutine test072 ( )

!*****************************************************************************80
!
!! TEST072 tests NEXCHR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 6

  character chr
  integer ( kind = 4 ) ichr
  integer ( kind = 4 ) jchr
  character ( len = 16 ), dimension ( test_num ) :: s_test = (/ &
    'Here I am!      ', &
    '   O  !         ', &
    'D o u b l e     ', &
    'T  r  i  p  l  e', &
    'F           a  r', &
    '         1      ' /)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST072'
  write ( *, '(a)' ) '  NEXCHR finds the next nonblank in a string.'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    jchr = 0

    do

      if ( jchr == 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  String is "' // trim ( s_test(test) ) // '"'
        write ( *, '(a)' ) ' '
      end if

      call nexchr ( s_test(test)(jchr+1:), ichr, chr )

      if ( ichr <= 0 ) then
        write ( *, '(a)' ) '  No more nonblanks!'
        exit
      end if

      jchr = jchr + ichr
      write ( *, '(a)' ) '  Next character is "' // chr // '".'

    end do

  end do

  return
end
subroutine test073 ( )

!*****************************************************************************80
!
!! TEST073 tests NEXSTR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 6

  integer ( kind = 4 ) isub
  integer ( kind = 4 ) jsub
  integer ( kind = 4 ) nsub
  character ( len = 16 ), dimension ( test_num ) :: s_test = (/ &
    'Here I am!      ', &
    '   O  !         ', &
    'D o u b l e     ', &
    'T  r  i  p  l  e', &
    'F           a  r', &
    '         1      ' /)
  character ( len = 2 ) sub
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST073'
  write ( *, '(a)' ) '  NEXSTR finds the next several characters in a string.'
  write ( *, '(a)' ) ' '

  nsub = 2

  do test = 1, test_num

    jsub = 0

    do

      if ( jsub == 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  String is "' // trim ( s_test(test) ) // '"'
        write ( *, '(a)' ) ' '
      end if

      call nexstr ( s_test(test)(jsub+1:), nsub, isub, sub )

      if ( isub <= 0 ) then
        write ( *, '(a)' ) '  No more nonblanks!'
        exit
      end if

      write ( *, '(a)' ) '  Next substring: ' // trim ( sub )
      jsub = jsub + isub

    end do

  end do

  return
end
subroutine test074 ( )

!*****************************************************************************80
!
!! TEST074 tests R4_TO_FLT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 10

  integer ( kind = 4 ) iexp
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) mant
  integer ( kind = 4 ) ndig
  real    ( kind = 4 ), dimension ( test_num ) :: r4_test = (/ &
    1.0E+00,   10.0E+00,  100.0E+00,  101.0E+00,   99.0E+00, &
    0.0E+00,   -1.0E+00, -123.456E+00,   -0.123456E+00,   0.000000123456E+00 /)
  real    ( kind = 4 ) rval
  real    ( kind = 4 ) sval
  integer ( kind = 4 ) test

  ndig = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST074'
  write ( *, '(a)' ) '  R4_TO_FLT computes the scientific representation'
  write ( *, '(a)' ) '  (floating point, base 10) of a real number.'
  write ( *, '(a)' ) ' '

  do ndig = 1, 6

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  The number of digits used is ', ndig
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '     RVAL        ISGN    MANT    IEXP      SVAL'
    write ( *, '(a)' ) ' '

    do test = 1, test_num

      rval = r4_test(test)
      call r4_to_flt ( rval, isgn, mant, iexp, ndig )
      sval = isgn * mant * 10.0E+00**iexp

      write ( *, '(g14.6,3i8,g14.6)' ) rval, isgn, mant, iexp, sval

    end do

  end do

  return
end
subroutine test075 ( )

!*****************************************************************************80
!
!! TEST075 tests R4_TO_S_LEFT, R4_TO_S_RIGHT, S_BLANKS_INSERT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real    ( kind = 4 ) rval
  character ( len = 40 ) s
  character ( len = 14 ) s2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST075'
  write ( *, '(a)' ) '  R4_TO_S_LEFT: real -> left justified string;'
  write ( *, '(a)' ) '  R4_TO_S_RIGHT: real -> right justified string.'
  write ( *, '(a)' ) '  S_BLANKS_INSERT inserts blanks in a string;'
  write ( *, '(a)' ) ' '

  s = 'There were guests.'
  write ( *, '(a)' ) '  Before call, STRING1 = "' // trim ( s ) // '"'

  call s_blanks_insert ( s, 11, 25 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  After inserting blanks into positions 11 through 25,'
  write ( *, '(a)' ) '  STRING1 = "' // trim ( s ) // '"'

  rval = 78.25

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now use R4_TO_S_RIGHT to turn the real value'
  write ( *, '(a,g14.6,a)' ) '  R = ', rval, ' into a right-justified string:'

  call r4_to_s_right ( rval, s2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  STRING2 = "' // trim ( s2 ) // '"'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now insert STRING2 into STRING1.'
  s(12:25) = s2
  write ( *, '(a)' ) '  The resulting string is "' // trim ( s ) // '"'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Repeating for R4_TO_S_LEFT:'

  s = 'There were guests.'
  call s_blanks_insert ( s, 11, 25 )

  rval = 78.25
  call r4_to_s_left ( rval, s2 )
  s(12:25) = s2

  write ( *, '(a)' ) '  The resulting string is "' // trim ( s ) // '"'

  return
end
subroutine test076 ( )

!*****************************************************************************80
!
!! TEST076 tests R4_TO_S_LEFT and S_TO_R4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 3

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) length
  real    ( kind = 4 ) r
  character ( len = 14 ) s
  character ( len = 14 ) s_test(test_num)
  character ( len = 14 ) s2
  integer ( kind = 4 ) test

  s_test(1) = ' 52.134ABCDE'
  s_test(2) = ' 8.0/2.0'
  s_test(3) = '12E1, 34, 56'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST076'
  write ( *, '(a)' ) '  S_TO_R4, string -> real number;'
  write ( *, '(a)' ) '  R4_TO_S_LEFT, real number -> string.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  S --> S_TO_R4 --> R4_TO_S_LEFT'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    s = s_test(test)
    write ( *, '(2x,a14,g14.6,a14)' ) s
    call s_to_r4 ( s, r, ierror, length )
    write ( *, '(2x,a14,g14.6,a14)' ) s, r
    call r4_to_s_left ( r, s2 )

    write ( *, '(2x,a14,g14.6,a14)' ) s, r, s2

  end do

  return
end
subroutine test077 ( )

!*****************************************************************************80
!
!! TEST077 tests R4_TO_S32 and S32_TO_R4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 4

  integer ( kind = 4 ) i
  real    ( kind = 4 ) r4
  character ( len = 32 ) r4_to_s32
  real    ( kind = 4 ), dimension ( test_num ) :: r4_test = (/ &
    0.0E+00, 1.0E+00, 7.0E+00, 15.0E+00 /)
  real    ( kind = 4 ) rval2
  character ( len = 32 ) s
  real    ( kind = 4 ) s32_to_r4
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST077'
  write ( *, '(a)' ) '  R4_TO_S32 converts a real to a string'
  write ( *, '(a)' ) '  S32_TO_R4 converts a string to a real.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  R   R4_TO_S32(R)   S32_TO_R4(R4_TO_S32(R))'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    r4 = r4_test(test)
    s = r4_to_s32 ( r4 )
    rval2 = s32_to_r4 ( s )

    write ( *, '( 2x, g14.6, 2x, a32, 2x, g14.6 )' ) r4, s, rval2

  end do

  return
end
subroutine test078 ( )

!*****************************************************************************80
!
!! TEST078 tests R4_TO_SEF and SEF_TO_R4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 16

  integer ( kind = 4 ) e
  integer ( kind = 4 ) e2
  integer ( kind = 4 ) f
  integer ( kind = 4 ) f2
  real    ( kind = 4 ) r
  real    ( kind = 4 ) r2
  real    ( kind = 4 ), dimension ( test_num ) :: r4_test = (/ &
      0.25E+00,  0.5E+00,   1.0E+00,   2.0E+00,  4.0E+00, &
      1.5E+00,   1.75E+00,  1.875E+00, 6.5E+00, -6.5E+00, &
     99.0E+00, 100.0E+00, 101.0E+00,   0.0E+00, -1.0E+00, &
    huge ( 1.0E+00 ) /)
  integer ( kind = 4 ) s
  integer ( kind = 4 ) s2
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST078'
  write ( *, '(a)' ) '  R4_TO_SEF converts an R4 to SEF form.'
  write ( *, '(a)' ) '  SEF_TO_R4 converts SEF form to an R4.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  S is the sign bit (0 = positive, 1 = negative)'
  write ( *, '(a)' ) '  E is the exponent base 2'
  write ( *, '(a)' ) '  F is the mantissa'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    R           S       E           F     R2'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    r = r4_test(test)

    call r4_to_sef ( r, s, e, f )

    call sef_to_r4 ( s, e, f, r2 )

    write ( *, '(2x,g16.8,i2,i8,i12,g16.8)' ) r, s, e, f, r2

  end do
!
!  Extra test values, some of which are unnormalized real quantities.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  S    E         F     R        S2   E2        F2'
  write ( *, '(a)' ) ' '

  s = 0
  e = -125
  f = 3

  call sef_to_r4 ( s, e, f, r )
  call r4_to_sef ( r, s2, e2, f2 )
  write ( *, '(2x,i2,i5,i10,g14.6,i2,i5,i10)' ) s, e, f, r, s2, e2, f2

  s = 0
  e = -127
  f = 3

  call sef_to_r4 ( s, e, f, r )
  call r4_to_sef ( r, s2, e2, f2 )
  write ( *, '(2x,i2,i5,i10,g14.6,i2,i5,i10)' ) s, e, f, r, s2, e2, f2

  s = 0
  e = -129
  f = 3

  call sef_to_r4 ( s, e, f, r )
  call r4_to_sef ( r, s2, e2, f2 )
  write ( *, '(2x,i2,i5,i10,g14.6,i2,i5,i10)' ) s, e, f, r, s2, e2, f2

  s = 0
  e = -132
  f = 7

  call sef_to_r4 ( s, e, f, r )
  call r4_to_sef ( r, s2, e2, f2 )
  write ( *, '(2x,i2,i5,i10,g14.6,i2,i5,i10)' ) s, e, f, r, s2, e2, f2

  s = 0
  e = -135
  f = 15

  call sef_to_r4 ( s, e, f, r )
  call r4_to_sef ( r, s2, e2, f2 )
  write ( *, '(2x,i2,i5,i10,g14.6,i2,i5,i10)' ) s, e, f, r, s2, e2, f2

  return
end
subroutine test079 ( )

!*****************************************************************************80
!
!! TEST079 tests R8_EXTRACT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) ierror
  character ( len = 80 ) s
  real    ( kind = 8 ) r

  s = '  12.3    45   -0.789'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST079'
  write ( *, '(a)' ) '  R8_EXTRACT extracts reals from a string.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Our string: "' // trim ( s ) // '"'
  write ( *, '(a)' ) ' '

  do

    call r8_extract ( s, r, ierror )

    if ( ierror /= 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6)' ) r

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Reached the last real value.'

  return
end
subroutine test080 ( )

!*****************************************************************************80
!
!! TEST080 tests R8_TO_S_LEFT, R8_TO_S_RIGHT and S_TO_R8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 4

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) length
  real    ( kind = 8 ) r
  character ( len = 20 ) s
  character ( len = 20 ) s_test(test_num)
  character ( len = 14 ) s2
  integer ( kind = 4 ) test

  s_test(1) = '52.134ABCDE'
  s_test(2) = ' 2.0/6.0'
  s_test(3) = '   12D-1, 34, 56'
  s_test(4) = '0.0001234'

  write ( *, '(a)') ' '
  write ( *, '(a)' ) 'TEST080'
  write ( *, '(a)' ) '  S_TO_R8,       string -> R8;'
  write ( *, '(a)' ) '  R8_TO_S_LEFT,  R8 -> left string.'
  write ( *, '(a)' ) '  R8_TO_S_RIGHT, R8 -> right string.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  S --> S_TO_R8 --> R8_TO_S_LEFT'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    s = s_test(test)
    call s_to_r8 ( s, r, ierror, length )
    call r8_to_s_left ( r, s2 )

    write ( *, '(2x,a20,2x,g14.6,2x,a14)' ) s, r, s2

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  S --> S_TO_R8 --> R8_TO_S_RIGHT'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    s = s_test(test)
    call s_to_r8 ( s, r, ierror, length )
    call r8_to_s_right ( r, s2 )

    write ( *, '(2x,a20,2x,g14.6,2x,a14)' ) s, r, s2

  end do
  return
end
subroutine test081

!*****************************************************************************80
!
!! TEST081 tests R8VEC_TO_S.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 6

  integer ( kind = 4 ) i
  character ( len = 100 ) s
  real ( kind = 8 ), dimension ( n ) :: x = (/ &
    1234.56D+00, &
    -0.00125D+00, &
    0.0D+00, &
    10203040506.0D+00, &
    77.0D+00, &
    1.5D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST081'
  write ( *, '(a)' ) '  R8VEC_TO_S writes an R8VEC to a string.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The real vector data:'
  write ( *, '(a)' ) ' '
  do i = 1, n
     write ( *, '(2x,i8,2x,g14.6)' ) i, x(i)
  end do

  call r8vec_to_s ( n, x, s )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The string: "' // trim ( s ) // '"'

  return
end
subroutine test082

!*****************************************************************************80
!
!! TEST082 tests RANGER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxval = 30

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ival(maxval)
  integer ( kind = 4 ) nval
  character ( len = 40 ) s

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST082'
  write ( *, '(a)' ) '  RANGER interprets a range description.'
  write ( *, '(a)' ) ' '

  s = ' 4:8 2 14:20 2:-1 81:81 10'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The input string is "' // trim ( s ) // '"'

  call ranger ( s, maxval, nval, ival )

  if ( nval <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  RANGER found no integers.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8,a)' ) '  RANGER found ', nval, ' integers:'
    write ( *, '(a)' ) ' '
    do i = 1, nval
      write ( *, '(2x,i8)' ) ival(i)
    end do

  end if

  return
end
subroutine test083

!*****************************************************************************80
!
!! TEST083 tests RAT_TO_S_LEFT and RAT_TO_S_RIGHT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 8

  integer ( kind = 4 ) i
  integer ( kind = 4 ) itest(test_num)
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) jtest(test_num)
  integer ( kind = 4 ) jval
  character ( len = 22 ) s1
  character ( len = 22 ) s2
  integer ( kind = 4 ) test

  itest(1) = 12
  jtest(1) = 10

  itest(2) = 48
  jtest(2) = -96

  itest(3) = -44
  jtest(3) = -44

  itest(4) = 23
  jtest(4) = 0

  itest(5) = -99
  jtest(5) = 0

  itest(6) = 0
  jtest(6) = 0

  itest(7) = 123456789
  jtest(7) = 987654321

  itest(8) = 0
  jtest(8) = 909

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST083'
  write ( *, '(a)' ) '  RAT_TO_S_LEFT prints a ratio left justified,'
  write ( *, '(a)' ) '  RAT_TO_S_RIGHT prints it right justified.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   IVAL      JVAL            Right          ' // &
    '         Left         '
  write ( *, '(a)' ) '  --------- ---------  ----------------------  ' // &
    '----------------------'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    ival = itest(test)
    jval = jtest(test)
    call rat_to_s_right ( ival, jval, s1 )
    call rat_to_s_left ( ival, jval, s2 )
    write ( *, '(2x,i10,i10,2x,a22,2x,a22)' ) ival, jval, s1, s2
  end do

  return
end
subroutine test085

!*****************************************************************************80
!
!! TEST085 tests S_ADJUSTL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 3

  character ( len = 10 ) s_test(test_num)
  character ( len = 10 ) s2
  integer ( kind = 4 ) test

  s_test(1) = '  Hello!  '
  s_test(2) = 'Ouch!'
  s_test(3) = '  A B C   '

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST085'
  write ( *, '(a)' ) '  S_ADJUSTL justifies a string to the left;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' )'   Original   S_ADJUSTL'
  write ( *, '(a)' )'  ----------  ---------- '
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    s2 = s_test(test)
    call s_adjustl ( s2 )

    write ( *, '(2x,a10,2x,a10)' ) s_test(test), s2

  end do

  return
end
subroutine test086

!*****************************************************************************80
!
!! TEST086 tests S_ADJUSTR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 3

  character ( len = 10 ) s_test(test_num)
  character ( len = 10 ) s2
  integer ( kind = 4 ) test

  s_test(1) = '  Hello!  '
  s_test(2) = 'Ouch!'
  s_test(3) = '  A B C   '

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST086'
  write ( *, '(a)' ) '  S_ADJUSTR justifies a string to the right.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' )'   Original   S_ADJUSTR'
  write ( *, '(a)' )'  ----------  ----------'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    s2 = s_test(test)
    call s_adjustr ( s2 )

    write ( *, '(2x,a10,2x,a10)' ) s_test(test), s2

  end do

  return
end
subroutine test087

!*****************************************************************************80
!
!! TEST087 tests S_AFTER_SS_COPY and S_BEFORE_SS_COPY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 6

  character ( len = 3 ) her
  integer ( kind = 4 ) ii
  character paren
  character ( len = 30 ) s_test(test_num)
  character ( len = 30 ) s2
  integer ( kind = 4 ) test

  paren = '('
  her = 'her'
  s_test(1) = 'John (or Jack)'
  s_test(2) = 'Jill St John (her real name)'
  s_test(3) = 'Jeff is OK (Rather!)'
  s_test(4) = 'FUNCTION SDOT(N,X,INCX,Y,INCY)'
  s_test(5) = 'Another remarkable string.'
  s_test(6) = 'On the (other (hand!!)'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST087'
  write ( *, '(a)' ) '  S_BEFORE_SS_COPY copies a string'
  write ( *, '(a)' ) '    before the first occurrence of a substring.'
  write ( *, '(a)' ) '  S_AFTER_SS_COPY copies a string'
  write ( *, '(a)' ) '    after the first occurrence of a substring.'
  write ( *, '(a)' ) ' '

  do ii = 1, 2

    write ( *, '(a)' ) ' '

    if ( ii == 1 ) then
      write ( *, '(a)' ) '  Our flag string is ' // paren
    else
      write ( *, '(a)' ) '  Our flag string is ' // her
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  String                          Copy'
    write ( *, '(a)' ) ' '

    do test = 1, test_num

      if ( ii == 1 ) then
        call s_before_ss_copy ( s_test(test), paren, s2 )
      else
        call s_before_ss_copy ( s_test(test), her, s2 )
      end if

      write ( *, '(2x,a30,2x,a30)' ) s_test(test), s2

    end do

  end do

  do ii = 1, 2

    write ( *, '(a)' ) ' '

    if ( ii == 1 ) then
      write ( *, '(a)' ) '  Our flag string is ' // paren
    else
      write ( *, '(a)' ) '  Our flag string is ' // her
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  String                          Copy'
    write ( *, '(a)' ) ' '

    do test = 1, test_num

      if ( ii == 1 ) then
        call s_after_ss_copy ( s_test(test), paren, s2 )
      else
        call s_after_ss_copy ( s_test(test), her, s2 )
      end if

      write ( *, '(2x,a30,2x,a30)' ) s_test(test), s2

    end do

  end do

  return
end
subroutine test088

!*****************************************************************************80
!
!! TEST088 tests S_ALPHA_LAST
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 4

  integer ( kind = 4 ) iloc
  character ( len = 20 ) s_test(test_num)
  integer ( kind = 4 ) test

  s_test(1) = 'HELLO World   !! !  '
  s_test(2) = '12345678901234567890'
  s_test(3) = '0.314159E+01'
  s_test(4) = '!@#$%a^&A(){}[]\\|<>?'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST088'
  write ( *, '(a)' ) '  S_ALPHA_LAST returns the location of the '
  write ( *, '(a)' ) '  last alphabetic character;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  ------String------  S_ALPHA_LAST'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    call s_alpha_last ( s_test(test), iloc )
    write ( *, '(2x,a20,2x,i8)' ) s_test(test), iloc
  end do

  return
end
subroutine test089

!*****************************************************************************80
!
!! TEST089 tests S_ANY_ALPHA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 4

  logical s_any_alpha
  character ( len = 20 ) s_test(test_num)
  integer ( kind = 4 ) test

  s_test(1) = 'HELLO World   !! !  '
  s_test(2) = '12345678901234567890'
  s_test(3) = '0.314159E+01'
  s_test(4) = '!@#$%a^&A(){}[]\\|<>?'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST089'
  write ( *, '(a)' ) '  S_ANY_ALPHA reports if a string'
  write ( *, '(a)' ) '  contains any alphabetic characters'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  ------String------  --S_ANY_ALPHA--'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    write ( *, '(2x,a20,2x,l1)' ) s_test(test), s_any_alpha ( s_test(test) )
  end do

  return
end
subroutine test090

!*****************************************************************************80
!
!! TEST090 tests S_BEGIN.
!
!  Discussion:
!
!    'Bob'          'BOB'     TRUE
!    '  B  o b '    ' bo b'   TRUE
!    'Bob'          'Bobby'   TRUE
!    'Bobo'         'Bobb'    FALSE
!    ' '            'Bob'     TRUE  (because blank matches anything)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 6

  logical s_begin
  character ( len = 12 ) s1
  character ( len = 12 ) s2
  character ( len = 12 ) s_test1(test_num)
  character ( len = 12 ) s_test2(test_num)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST090'
  write ( *, '(a)' ) '  S_BEGIN checks the beginning of a string for a'
  write ( *, '(a)' ) '  substring, ignoring case and spaces.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     S1          S2     S_BEGIN(S1,S2)'
  write ( *, '(a)' ) ' '

  s_test1(1) = 'Bob'
  s_test1(2) = ' B  o b'
  s_test1(3) = 'Bob'
  s_test1(4) = 'Bobo'
  s_test1(5) = ' '
  s_test1(6) = 'cubic meter'

  s_test2(1) = 'BOB'
  s_test2(2) = ' bo b'
  s_test2(3) = 'BOBBY'
  s_test2(4) = 'Bobb'
  s_test2(5) = 'Bob'
  s_test2(6) = 'cubic meter'

  do test = 1, test_num

    s1 = s_test1(test)
    s2 = s_test2(test)

    write ( *, '(2x,a,2x,a,2x,l1)' ) s1, s2, s_begin ( s1, s2 )

  end do

  return
end
subroutine test091

!*****************************************************************************80
!
!! TEST091 tests S_BEHEAD_SUBSTRING
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 4

  character ( len = 20 ) s_test(test_num)
  character ( len = 20 ) s_old
  character ( len = 20 ) sub(test_num)
  integer ( kind = 4 ) test

  s_test(1) = '    HELLO World!'
  sub(1) = 'HELLO'
  s_test(2) = '12345678901234567890'
  sub(2) = '12345'
  s_test(3) = '0.314159E+01'
  sub(3) = '314'
  s_test(4) = '!@#$%a^&A(){}[]\\|<>?'
  sub(4) = '!@#$%a^&A(){}[]\\|<>?'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST091'
  write ( *, '(a)' ) '  S_BEHEAD_SUBSTRING removes an initial substring from a '
  write ( *, '(a)' ) '  string, if it occurs'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  ------String--------  -----SUB------------  ---Beheaded----'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    s_old = s_test(test)
    call s_behead_substring ( s_test(test), sub(test) )
    write ( *, '(2x,a20,2x,a20,2x,a20)' ) s_old, sub(test), s_test(test)
  end do

  return
end
subroutine test092

!*****************************************************************************80
!
!! TEST092 tests S_BLANK_DELETE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 20 ) s

  s = 'HELLO World   !! !  '

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST092'
  write ( *, '(a)' ) '  S_BLANK_DELETE removes all blanks.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Input:  "' // trim ( s ) // '"'

  call s_blank_delete ( s )

  write ( *, '(a)' ) '  Output: "' // trim ( s ) // '"'

  return
end
subroutine test093

!*****************************************************************************80
!
!! TEST093 tests S_BLANKS_DELETE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 20 ) s

  s = 'HELLO World   !! !  '

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST093'
  write ( *, '(a)' ) '  S_BLANKS_DELETE removes double blanks.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Input:  ' // trim ( s )

  call s_blanks_delete ( s )

  write ( *, '(a)' ) '  Output: ' // trim ( s )

  return
end
subroutine test094

!*****************************************************************************80
!
!! TEST094 tests S_CAP, S_LOW and S_WORD_CAP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 5

  character ( len = 20 ) s_test(test_num)
  character ( len = 20 ) s1
  character ( len = 20 ) s2
  character ( len = 20 ) s3
  integer ( kind = 4 ) test

  s_test(1) = 'HELLO World   !! !  '
  s_test(2) = '12345678901234567890'
  s_test(3) = 'Abc Def Ghi Jkl Mno '
  s_test(4) = '!@#$%a^&A(){}[]\\|<>?'
  s_test(5) = 'a waste is a terrible thing to mind.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST094'
  write ( *, '(a)' ) '  S_CAP capitalizes all characters in a string;'
  write ( *, '(a)' ) '  S_LOW lowercases all characters;'
  write ( *, '(a)' ) '  S_WORD_CAP initial-capitalizes words in a string;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  ------Original------  -----Capitalized-----' // &
    '-----Lower Cased-----  -----Word_Caps-----'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    s1 = s_test(test)
    call s_cap ( s1 )
    s2 = s_test(test)
    call s_low ( s2 )
    s3 = s_test(test)
    call s_word_cap ( s3 )

    write ( *, '(2x,a20,2x,a20,2x,a20,2x,a20)' ) s_test(test), s1, s2, s3

  end do

  return
end
subroutine test095

!*****************************************************************************80
!
!! TEST095 tests S_CAT and S_CAT1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 5 ) s1
  character ( len = 5 ) s2
  character ( len = 10 ) s3
  character ( len = 10 ) s4
  character ( len = 10 ) s5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST095'
  write ( *, '(a)' ) '  // concatenates two strings;'
  write ( *, '(a)' ) '  S_CAT concatenates two strings, trimming blanks;'
  write ( *, '(a)' ) '  S_CAT1 concatenates two strings with a'
  write ( *, '(a)' ) '    single blank separator.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  --S1-  --S2-   --S1//S2--   --S_CAT--  --S_CAT1--'
  write ( *, '(a)' ) ' '

  s1 = 'Cat'
  s2 = 'fish'

  s3 = s1 // s2
  call s_cat ( s1, s2, s4 )
  call s_cat1 ( s1, s2, s5 )

  write ( *, '(2x,a,5x,a,5x,a,5x,a,5x,a)' ) s1, s2, s3, s4, s5

  return
end
subroutine test096

!*****************************************************************************80
!
!! TEST096 tests S_CENTER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 14 ) string1
  character ( len = 14 ) string2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST096'
  write ( *, '(a)' ) '  S_CENTER centers a string.'

  string1 = 'A'
  string2 = string1

  call s_center ( string2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   12345677654321'
  write ( *, '(a)' ) '  "' // string1 // '"'
  write ( *, '(a)' ) '  "' // string2 // '"'
  write ( *, '(a)' ) '   12345677654321'

  string1 = '  B  C  '
  string2 = string1

  call s_center ( string2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   12345677654321'
  write ( *, '(a)' ) '  "' // string1 // '"'
  write ( *, '(a)' ) '  "' // string2 // '"'
  write ( *, '(a)' ) '   12345677654321'

  string1 = '     67   4   '
  string2 = string1

  call s_center ( string2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   12345677654321'
  write ( *, '(a)' ) '  "' // string1 // '"'
  write ( *, '(a)' ) '  "' // string2 // '"'
  write ( *, '(a)' ) '   12345677654321'

  return
end
subroutine test097

!*****************************************************************************80
!
!! TEST097 tests S_CENTER_INSERT, S_LEFT_INSERT and S_RIGHT_INSERT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 10 ) string1
  character ( len = 30 ) string2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST097'
  write ( *, '(a)' ) '  S_LEFT_INSERT inserts a string left of another;'
  write ( *, '(a)' ) '  S_CENTER_INSERT inserts it in the center;'
  write ( *, '(a)' ) '  S_RIGHT_INSERT inserts it to the right.'
  write ( *, '(a)' ) ' '

  string1 = 'ZOWIE'
  string2 = '123456789012345678901234567890'

  write ( *, '(a)' ) '  The string to be inserted is: ' // trim ( string1 )
  write ( *, '(a)' ) '  The string in which we insert is: ' // trim ( string2 )

  call s_left_insert ( string1, string2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  After calling S_LEFT_INSERT:'
  write ( *, '(a)' ) '  "' // trim ( string2 ) // '"'

  string1 = 'ZOWIE'
  string2 = '123456789012345678901234567890'

  call s_center_insert ( string1, string2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  After calling S_CENTER_INSERT: '
  write ( *, '(a)' ) '  "' // trim ( string2 ) // '"'

  string1 = 'ZOWIE'
  string2 = '123456789012345678901234567890'

  call s_right_insert ( string1, string2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  After calling S_RIGHT_INSERT:  '
  write ( *, '(a)' ) '  "' // trim ( string2 ) // '"'

  return
end
subroutine test0975 ( )

!*****************************************************************************80
!
!! TEST0975 tests S_CH_COUNT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character              ch
  integer   ( kind = 4 ) ch_count
  character ( len = 30 ) s

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0975'
  write ( *, '(a)' ) '  S_CH_COUNT counts occurrences of a character.'

  s = 'Bob is debobbing the bobber!'
  ch = 'b'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  String =     "' // trim ( s ) // '".'
  write ( *, '(a)' ) '  Character is "' // ch // '".'

  call s_ch_count ( s, ch, ch_count )

  write ( *, '(a,i8)' ) '  Number of occurrences = ', ch_count

  return
end
subroutine test098

!*****************************************************************************80
!
!! TEST098 tests S_CH_DELETE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 4

  character c(test_num)
  character ( len = 35 ) s_test(test_num)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST098'
  write ( *, '(a)' ) '  S_CH_DELETE removes a character from a string.'
  write ( *, '(a)' ) ' '

  s_test(1) = 'A man, a plan, a canal, Panama!'
  c(1) = ' '

  s_test(2) = 'A man, a plan, a canal, Panama!'
  c(2) = 'a'

  s_test(3) = 'A man, a plan, a canal, Panama!'
  c(3) = 'n'

  s_test(4) = 'aaaaannnnnQ!'
  c(4) = 'n'

  do test = 1, test_num

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Remove "' // c(test) // '" from "' &
      // trim ( s_test(test) ) // '"'
    write ( *, '(a)' )
    call s_ch_delete ( s_test(test), c(test) )

    write ( *, '(a)' ) '  Result: ' // trim ( s_test(test) )

  end do

  return
end
subroutine test099

!*****************************************************************************80
!
!! TEST099 tests S_CH_LAST.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 5

  character ( len = 46 ), dimension ( test_num ) :: s_test = (/ &
    'HELLO World   !! !                  ', &
    '12345678901234567890                ', &
    'Abc Def Ghi Jkl Mno                 ', &
    '!@#$%a^&A(){}[]\\|<>?               ', &
    'a taste is a wearable thing to mind.' /)
  character s_ch_last
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST099'
  write ( *, '(a)' ) '  S_CH_LAST returns the last nonblank in a string.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  ------String------       Last'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    write ( *, '(2x,a20,10x,a1)' ) s_test(test), s_ch_last ( s_test(test) )
  end do

  return
end
subroutine test100

!*****************************************************************************80
!
!! TEST100 tests S_CHOP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  character ( len = 30 ) s

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST100'
  write ( *, '(a)' ) '  S_CHOP chops out part of a string.'
  write ( *, '(a)' ) ' '

  s = 'CHRPAK is not working today!'
  write ( *, '(a)' ) '  Original string = "' // trim ( s ) // '"'
  ilo = 11
  ihi = 14
  write ( *, '(2x,a,i8,a,i8)' ) '  We delete entries ', ilo, ' to ', ihi
  call s_chop ( s, ilo, ihi )
  write ( *, '(a)' ) '  Chopped string = "' // trim ( s ) // '"'

  return
end
subroutine test101

!*****************************************************************************80
!
!! TEST101 tests S_DETAG.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 60 ) s

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST101'
  write ( *, '(a)' ) '  S_DETAG removes HTML tags from a string.'

  s = 'This is <I>italic</I> whereas this <B>boldly</B> goes on!'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Original string:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  "' // trim ( s ) // '"'

  call s_detag ( s )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Detagged string:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  "' // trim ( s ) // '"'

  s = 'This is an example <A HREF = "stuff.html">of a link </A>.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Original string:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  "' // trim ( s ) // '"'

  call s_detag ( s )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Detagged string:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  "' // trim ( s ) // '"'

  return
end
subroutine test1013 ( )

!*****************************************************************************80
!
!! TEST1013 tests S_DETROFF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character, parameter :: bs = char ( 8 )
  integer   ( kind = 4 ) i
  character ( len = 80 ) s
  integer   ( kind = 4 ) s_length

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1013'
  write ( *, '(a)' ) '  S_DETROFF deletes CHARACTER+Backspace pairs.'

  s = '_#T_#e_#x_#t ##is## B#B#B#Bo#o#o#ol#l#l#ld#d#d#d'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  String with symbolic backspaces:'
  write ( *, '(a)' ) '    "' // trim ( s ) // '".'
  write ( *, '(a,i8)' ) '  Length is ', len_trim ( s )

  s_length = len ( s )

  do i = 1, s_length
    if ( s(i:i) == '#' ) then
      s(i:i) = bs
    end if
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  String with actual backspaces'
  write ( *, '(a)' ) '    "' // trim ( s ) // '".'
  write ( *, '(a,i8)' ) '  Length is ', len_trim ( s )


  call s_detroff ( s )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  De-TROFF''ed string:'
  write ( *, '(a)' ) '    "' // trim ( s ) // '".'
  write ( *, '(a,i8)' ) '  Length is ', len_trim ( s )

  return
end
subroutine test1015

!*****************************************************************************80
!
!! TEST1015 tests S_EQI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  logical s_eqi
  character ( len = 10 ) s1
  character ( len = 10 ) s2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1015'
  write ( *, '(a)' ) '  S_EQI compares two strings for equality,'
  write ( *, '(a)' ) '  ignoring case and trailing blanks.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   A      B    S_EQI(A,B)'
  write ( *, '(a)' ) ' '

  s1 = 'NixoN'
  s2 = 'niXon'
  write ( *, '(2x,a10,2x,a10,2x,l1)' ) s1, s2, s_eqi ( s1, s2 )

  s1 = 'animal'
  s2 = 'CRACKER'
  write ( *, '(2x,a10,2x,a10,2x,l1)' ) s1, s2, s_eqi ( s1, s2 )

  s1 = 'Yes'
  s2 = 'y'
  write ( *, '(2x,a10,2x,a10,2x,l1)' ) s1, s2, s_eqi ( s1, s2 )

  s1 = 'ALPHA'
  s2 = 'zeta'
  write ( *, '(2x,a10,2x,a10,2x,l1)' ) s1, s2, s_eqi ( s1, s2 )

  s1 = 'NIX on'
  s2 = 'Nixon'
  write ( *, '(2x,a10,2x,a10,2x,l1)' ) s1, s2, s_eqi ( s1, s2 )

  s1 = 'blank'
  s2 = 'blank     '
  write ( *, '(2x,a10,2x,a10,2x,l1)' ) s1, s2, s_eqi ( s1, s2 )

  return
end
subroutine test102

!*****************************************************************************80
!
!! TEST102 tests S_ESCAPE_TEX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 80 ) s1
  character ( len = 80 ) s2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST102'
  write ( *, '(a)' ) '  S_ESCAPE_TEX "protects" characters in a string'
  write ( *, '(a)' ) '  that might otherwise be interpreted as TeX'
  write ( *, '(a)' ) '  escape characters.'

  s1 = 'The file A_B.TXT is {I think__so} of size 2^8 or C\B.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Original string:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  "' // trim ( s1 ) // '"'

  call s_escape_tex ( s1, s2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  De-escaped string:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  "' // trim ( s2 ) // '"'

  return
end
subroutine test103

!*****************************************************************************80
!
!! TEST103 tests S_FILL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character c
  character ( len = 10 ) s

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST103'
  write ( *, '(a)' ) '  S_FILL fills a string with a character.'
  write ( *, '(a)' ) ' '

  s = 'My word!'
  write ( *, '(2x,a,a)' ) '  Before: ', trim ( s )
  c = '$'
  call s_fill ( s, c )
  write ( *, '(2x,a,a)' ) '  After:  ', trim ( s )

  return
end
subroutine test104

!*****************************************************************************80
!
!! TEST104 tests S_GEI, S_GTI, S_LEI, S_LTI, S_NEQI, S_EQI, S_EQIDB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 5

  character ( len = 7 ) a(test_num)
  character ( len = 7 ) b(test_num)
  logical comp(test_num,14)
  integer ( kind = 4 ) i
  logical s_eqi
  logical s_eqidb
  logical s_gei
  logical s_gti
  logical s_lei
  logical s_lti
  logical s_neqi
  integer ( kind = 4 ) test

  a(1) = 'NixoN'
  b(1) = 'niXon'
  a(2) = 'animal'
  b(2) = 'CRACKER'
  a(3) = 'Yes'
  b(3) = 'y'
  a(4) = 'ALPHA'
  b(4) = 'zeta'
  a(5) = 'NIX on'
  b(5) = 'Nixon'

  do test = 1, test_num

    comp(test,1) = a(test) == b(test)
    comp(test,2) = a(test) == b(test)
    comp(test,3) = lge ( a(test), b(test) )
    comp(test,4) = lgt ( a(test), b(test) )
    comp(test,5) = lle ( a(test), b(test) )
    comp(test,6) = llt ( a(test), b(test) )
    comp(test,7) = a(test) /= b(test)

    comp(test,8) = s_eqi ( a(test), b(test) )
    comp(test,9) = s_eqidb ( a(test), b(test) )
    comp(test,10) = s_gei ( a(test), b(test) )
    comp(test,11) = s_gti ( a(test), b(test) )
    comp(test,12) = s_lei ( a(test), b(test) )
    comp(test,13) = s_lti ( a(test), b(test) )
    comp(test,14) = s_neqi ( a(test), b(test) )

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST104'
  write ( *, '(a)' ) '  For implicitly capitalized strings S1 and S2'
  write ( *, '(a)' ) '  S_EQI,   S1 =  S2'
  write ( *, '(a)' ) '  S_EQIDB, S1 =  S2, blank insensitive'
  write ( *, '(a)' ) '  S_GEI    S1 >= S2'
  write ( *, '(a)' ) '  S_GTI    S1 >  S2'
  write ( *, '(a)' ) '  S_LEI    S1 <= S2'
  write ( *, '(a)' ) '  S_LTI    S1 <  S2'
  write ( *, '(a)' ) '  S_NEQI   S1 != S2'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Results of    "A compare B"'
  write ( *, '(a)' ) '  First line is FORTRAN (case sensitive)'
  write ( *, '(a)' ) '  Second line is CHRPAK (case insensitive)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   A      B    =   =  = >  >   < =   <    = / = '
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    write ( *, '(2x,a7,2x,a7,7(3x,l1))' ) a(test), b(test), comp(test,1:7)
    write ( *, '(2x,7x,2x,7x,7(3x,l1))' ) comp(test,8:14)
    write ( *, '(a)' ) ' '
  end do

  return
end
subroutine test105

!*****************************************************************************80
!
!! TEST105 tests S_INC_C.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  character ( len = 30 ) s

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST105'
  write ( *, '(a)' ) '  S_INC_C can "increment" the characters in a string.'

  s = 'Tax'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Starting string: "' // trim ( s ) // '"'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Incremented forms:'
  write ( *, '(a)' ) ' '

  do i = 1, 5
    call s_inc_c ( s )
    write ( *, '(2x,a)' ) '  "' // trim ( s ) // '"'
  end do

  s = 'aB34c* 8zY'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Starting string: "' // trim ( s ) // '"'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Incremented forms:'
  write ( *, '(a)' ) ' '

  do i = 1, 5
    call s_inc_c ( s )
    write ( *, '(2x,a)' ) '  "' // trim ( s ) // '"'
  end do

  return
end
subroutine test1055

!*****************************************************************************80
!
!! TEST1055 tests S_INC_N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  character ( len = 30 ) s

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1055'
  write ( *, '(a)' ) '  S_INC_N can "increment" the numeric part'
  write ( *, '(a)' ) '  of a file name.'

  s = 'data01.txt'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Starting string: "' // trim ( s ) // '"'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Incremented forms:'
  write ( *, '(a)' ) ' '

  do i = 1, 5
    call s_inc_n ( s )
    write ( *, '(2x,a)' ) '  "' // trim ( s ) // '"'
  end do

  s = 'mat09lab98.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Starting string: "' // trim ( s ) // '"'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Incremented forms:'
  write ( *, '(a)' ) ' '

  do i = 1, 5
    call s_inc_n ( s )
    write ( *, '(2x,a)' ) '  "' // trim ( s ) // '"'
  end do

  return
end
subroutine test106

!*****************************************************************************80
!
!! TEST106 tests S_INDEX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) i4
  integer ( kind = 4 ) s_index
  character ( len = 30 ) s
  character ( len = 10 ) substring

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST106'
  write ( *, '(a)' ) '  S_INDEX reports the first occurrence of a substring.'
  write ( *, '(a)' ) '    The comparison ignores trailing blanks.'

  s = 'Bob is debobbing the bobber!'
  substring = 'bob'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  String =     ' // trim ( s )
  write ( *, '(a)' ) '  Substring is ' // trim ( substring )

  i1 =   index ( s, trim ( substring ) )
  i2 = s_index ( s, trim ( substring ) )
  i3 =   index ( s,        substring )
  i4 = s_index ( s,        substring )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '    INDEX ( S, TRIM ( SUBSTRING ) ) = ', i1
  write ( *, '(a,i8)' ) '  S_INDEX ( S, TRIM ( SUBSTRING ) ) = ', i2
  write ( *, '(a,i8)' ) '    INDEX ( S,        SUBSTRING   ) = ', i3
  write ( *, '(a,i8)' ) '  S_INDEX ( S,        SUBSTRING   ) = ', i4

  return
end
subroutine test107

!*****************************************************************************80
!
!! TEST107 tests S_INDEX_SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character blank
  character hat
  integer ( kind = 4 ) i
  integer ( kind = 4 ) loc_new
  integer ( kind = 4 ) loc_old
  character ( len = 40 ) s
  character ( len = 10 ) s2
  integer ( kind = 4 ) s_index_set

  blank = ' '
  hat = '^'
  s2 = '0123456789'
  s = '1 way 4 U 2 deb8 of10 is 2 Rgu!'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST107'
  write ( *, '(a)' ) '  S_INDEX_SET searches a string for any character'
  write ( *, '(a)' ) '    in a given set.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  String: ' // trim ( s )
  write ( *, '(a)' ) '  Set:	' // trim ( s2 )
  write ( *, '(a)' ) ' '

  loc_new = 0

  do

    loc_old = loc_new

    loc_new = s_index_set ( s(loc_old+1:), s2 ) + loc_old

    if ( loc_new == loc_old ) then
      exit
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  "' // trim ( s ) // '"'
    write ( *, '(40a)' ) ( blank, i = 1, loc_new-1 ), hat

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  No more matches found.'

  return
end
subroutine test108

!*****************************************************************************80
!
!! TEST108 tests S_INDEX_LAST and S_INDEXI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) i4
  integer ( kind = 4 ) s_indexi
  integer ( kind = 4 ) s_index_last
  character ( len = 30 ) s
  character ( len = 10 ) substring

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST108'
  write ( *, '(a)' ) '  S_INDEXI reports the first occurrence of a'
  write ( *, '(a)' ) '    substring, case and trailing space'
  write ( *, '(a)' ) '    insensitive.'
  write ( *, '(a)' ) '  S_INDEX_LAST reports the LAST occurrence'
  write ( *, '(a)' ) '    of a substring.'
  write ( *, '(a)' ) '  INDEX is a case and trailing space sensitive'
  write ( *, '(a)' ) '    routine which reports the first occurrence'
  write ( *, '(a)' ) '    of a substring.'

  s = 'Bob is debobbing the bobber!'
  substring = 'bob'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  String =     ' // trim ( s )
  write ( *, '(a)' ) '  Substring is ' // trim ( substring )

  i1 = index ( s, substring )
  i2 = index ( s, trim ( substring ) )
  i3 = s_indexi ( s, substring )
  i4 = s_index_last ( s, substring )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  INDEX =              ', i1
  write ( *, '(a,i8)' ) '  INDEX (restricted) = ', i2
  write ( *, '(a,i8)' ) '  INDEXI =             ', i3
  write ( *, '(a,i8)' ) '  S_INDEX_LAST =       ', i4

  return
end
subroutine test109

!*****************************************************************************80
!
!! TEST109 tests S_INDEX_LAST_C.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character c
  integer ( kind = 4 ) i
  character i4_to_a
  integer ( kind = 4 ) j
  character ( len = 60 ) s
  integer ( kind = 4 ) s_index_last_c

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST109'
  write ( *, '(a)' ) '  S_INDEX_LAST_C reports the LAST occurrence'
  write ( *, '(a)' ) '    of a character.'

  s = 'The quick brown fox jumps right over the big lazy dog!'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  String =     ' // trim ( s )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I     C     J'
  write ( *, '(a)' ) ' '

  do i = 27, 52
    c = i4_to_a ( i )
    j = s_index_last_c ( s, c )
    write ( *, '(2x,i8,5x,a1,i8)' ) i, c, j
  end do

  return
end
subroutine test110

!*****************************************************************************80
!
!! TEST110 tests S_IS_DIGIT and S_IS_I.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 6

  integer ( kind = 4 ) ival
  logical lval1
  logical lval2
  logical s_is_digit
  logical s_is_i
  character ( len = 10 ) s_test(test_num)
  integer ( kind = 4 ) test

  s_test(1) = '123 '
  s_test(2) = ' 1.2 - 3'
  s_test(3) = ' A4'
  s_test(4) = '-3.14E+2'
  s_test(5) = ' 2 3 4 '
  s_test(6) = ' +2, '

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST110'
  write ( *, '(a)' ) '  S_IS_DIGIT reports whether a string'
  write ( *, '(a)' ) '    contains only digits.'
  write ( *, '(a)' ) '  S_IS_I reports whether a string'
  write ( *, '(a)' ) '    represents a single integer.'
  write ( *, '(a)' ) ' '
  ival = 0

  do test = 1, test_num
    lval1 = s_is_digit ( s_test(test) )
    lval2 = s_is_i ( s_test(test), ival )
    write ( *, '(2x,a10,2x,l1,2x,l1,2x,i8)' ) s_test(test), lval1, lval2, ival
  end do

  return
end
subroutine test111

!*****************************************************************************80
!
!! TEST111 tests S_IS_F77_NAME and S_IS_F90_NAME.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 9

  logical s_is_f77_name
  logical s_is_f90_name
  character ( len = 10 ), dimension ( test_num ) :: s_test = (/ &
    'arthur    ', &
    'art hur   ', &
    '    Mario ', &
    '3.14159   ', &
    'zo#wy     ', &
    '          ', &
    'R2D2      ', &
    'A_1       ', &
    '_A1       ' /)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST111'
  write ( *, '(a)' ) '  S_IS_F77_NAME reports if a string is a'
  write ( *, '(a)' ) '    legal FORTRAN-77 identifier.'
  write ( *, '(a)' ) '  S_IS_F90_NAME reports if a string is a'
  write ( *, '(a)' ) '    legal FORTRAN-90 identifier.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  -------String-------  F77?     F90?'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    write ( *, '(2x,a,5x,l1,9x,l1)' ) &
      s_test(test), s_is_f77_name ( s_test(test) ), &
      s_is_f90_name ( s_test(test) )

  end do

  return
end
subroutine test112

!*****************************************************************************80
!
!! TEST112 tests S_IS_R.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 6

  logical lval
  real    ( kind = 4 ) rval
  character ( len = 10 ), dimension ( test_num ) :: s_test = (/ &
    '123       ', &
    ' 1.2 - 3  ', &
    ' A4.5     ', &
    '-3.14E+2  ', &
    ' 2 3 4    ', &
    ' +2.3,    ' /)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST112'
  write ( *, '(a)' ) '  S_IS_R reports whether a string'
  write ( *, '(a)' ) '    represents a single real value.'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    call s_is_r ( s_test(test), rval, lval )
    write ( *, '(2x,a10,2x,l1,2x,g14.6)' ) s_test(test), lval, rval
  end do

  return
end
subroutine test113

!*****************************************************************************80
!
!! TEST113 tests S_ONLY_ALPHAB and S_ONLY_DIGITB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 9

  character ( len = 4 ), dimension ( test_num ) :: s_test = (/ &
    '1984', &
    'Fred', &
    'C3PO', &
    '/#4D', &
    ' Bc ', &
    '2 34', &
    '-198', &
    '8 +4', &
    '10*8' /)
  logical s_only_alphab
  logical s_only_digitb
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST113'
  write ( *, '(a)' ) '  S_ONLY_ALPHAB reports if a string is only'
  write ( *, '(a)' ) '    alphabetic and blanks.'
  write ( *, '(a)' ) '  S_ONLY_DIGITB reports if a string is only digits and blanks.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  S       S_ONLY_DIGITB  S_ONLY_ALPHAB'
  write ( *, '(a)') ' '

  do test = 1, test_num
    write ( *, '(2x,3x,a4,5x,l1,5x,l1)' ) &
      s_test(test), s_only_digitb( s_test(test) ), &
      s_only_alphab( s_test(test) )
  end do

  return
end
subroutine test114

!*****************************************************************************80
!
!! TEST114 tests S_OVERLAP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 5

  integer ( kind = 4 ) overlap
  character ( len = 10 ) s1
  character ( len = 10 ), save, dimension ( test_num ) :: s1_test = (/ &
    'timber    ', &
    'timber    ', &
    'beret     ', &
    'beret     ', &
    'beret     ' /)
  character ( len = 10 ) s2
  character ( len = 10 ), save, dimension ( test_num ) :: s2_test = (/ &
    'beret     ', &
    'timber    ', &
    'timber    ', &
    'berets    ', &
    'berth     ' /)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST114'
  write ( *, '(a)' ) '  S_OVERLAP measures the overlap between two strings.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  S1          S2          Overlap'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    s1 = s1_test(test)
    s2 = s2_test(test)
    call s_overlap ( s1, s2, overlap )
    write ( *, '(2x,a,3x,a,3x,i2)' ) s1, s2, overlap
  end do

  return
end
subroutine test115

!*****************************************************************************80
!
!! TEST115 tests S_REPLACE_CH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character c1
  character c2
  character ( len = 15 ) s
  character ( len = 15 ) :: s_old = 'No pennies now.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST115'
  write ( *, '(a)' ) '  S_REPLACE_CH replaces one character by another;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    C1  C2  Original String  Modified String'
  write ( *, '(a)' ) ' '

  c1 = 'n'
  c2 = 't'

  s = s_old
  call s_replace_ch ( s, c1, c2 )

  write ( *, '(5x,a1,3x,a1,2x,a,2x,a)' ) c1, c2, s_old, s

  return
end
subroutine test116

!*****************************************************************************80
!
!! TEST116 tests S_REPLACE_ONE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 35 ) s1
  character ( len = 35 ) s2
  character ( len = 2 ) :: sub1 = 'an'
  character ( len = 4 ) :: sub2 = 'ORK '

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST116'
  write ( *, '(a)' ) '  S_REPLACE_ONE replaces one occurrence of a string.'

  s1 = 'A man, a plan, a canal - Panama!'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Replace the first occurrence of '
  write ( *, '(4x,a)' ) '"' // trim ( sub1 ) // ' by "' // trim ( sub2 ) &
    // '" in "' // trim ( s1 ) // '"'

  call s_replace_one ( s1, sub1, sub2, s2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Result:'
  write ( *, '(4x,a)' ) '"' // trim ( s2 ) // '"'

  return
end
subroutine test117

!*****************************************************************************80
!
!! TEST117 tests S_REPLACE_REC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) irep
  character ( len = 35 ) s
  character ( len = 2 ) sub1a
  character sub2a

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST117'
  write ( *, '(a)' ) '  S_REPLACE_REC recursively replaces a string.'
  write ( *, '(a)' ) ' '

  s = 'aaaaannnnnBC'
  sub1a = 'an'
  sub2a = 'a'
  write ( *, '(a)' ) '  Replace all occurrences of '
  write ( *, '(4x,a)' ) trim ( sub1a ) // ' by ' // trim ( sub2a ) &
    // ' in ' // trim ( s )
  write ( *, '(a)' ) ' '
  call s_replace_rec ( s, sub1a, sub2a, irep )
  write ( *, '(a)' ) '  Result "' // trim ( s ) // '"'
  write ( *, '(2x,i8,a)' ) irep, ' replacements were made.'

  return
end
subroutine test118

!*****************************************************************************80
!
!! TEST118 tests S_REPLACE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) irep
  character ( len = 35 ) string
  character ( len = 3 ) sub1
  character ( len = 3 ) sub2
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST118'
  write ( *, '(a)' ) '  S_REPLACE replaces a pattern in a string.'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    string = 'A man, a plan, a canal, Panama!'

    if ( test == 1 ) then
      sub1 = 'an'
      sub2 = '&@'
    else if ( test == 2 ) then
      sub1 = 'an,'
      sub2 = '8'
    else if ( test == 3 ) then
      sub1 = 'a'
      sub2 = 'oro'
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Replace all occurrences of '
    write ( *, '(4x,a)' ) trim ( sub1 ) // ' by ' // trim ( sub2 ) &
      // ' in ' // trim ( string )
    write ( *, '(a)' )
    call s_replace ( string, sub1, sub2, irep )
    write ( *, '(a)' ) '  Result: ' // trim ( string )
    write ( *, '(2x,i8,a)' ) irep, ' replacements were made.'

  end do

  return
end
subroutine test119

!*****************************************************************************80
!
!! TEST119 tests S_REVERSE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 35 ) s

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST119'
  write ( *, '(a)' ) '  S_REVERSE reverses a string.'
  write ( *, '(a)' ) ' '

  s = 'A man, a plan, a canal, Panama!'
  write ( *, '(2x,a,a)' ) '  Before: "' // trim ( s ) // '"'
  call s_reverse ( s )
  write ( *, '(2x,a,a)' ) '  After: "' // trim ( s ) // '"'

  return
end
subroutine test120

!*****************************************************************************80
!
!! TEST120 tests S_S_DELETE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 4

  integer ( kind = 4 ) irep
  character ( len = 31 ) s
  character ( len = 31 ), dimension ( test_num ) :: s_test = (/ &
    'A man, a plan, a canal, Panama!', &
    'A man, a plan, a canal, Panama!', &
    'A man, a plan, a canal, Panama!', &
    'aaaaannnnnQ!                   ' /)
  character ( len = 5 ) sub
  character ( len = 5 ), dimension ( test_num ) :: sub_test = (/ &
    ',    ', &
    'an   ', &
    'canal', &
    'an   ' /)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST120'
  write ( *, '(a)' ) '  S_S_DELETE removes a substring;'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    s = s_test(test)
    sub = sub_test(test)

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Remove "' // &
      trim ( sub ) // '" from "' // trim ( s ) // '"'

    call s_s_delete ( s, trim ( sub ), irep )

    write ( *, '(a)' )
    write ( *, '(a)' ) '  Result: ' // trim ( s_test(test) )
    write ( *, '(2x,i8,a)' ) irep, ' removals'

  end do

  return
end
subroutine test121

!*****************************************************************************80
!
!! TEST121 tests S_S_DELETE2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 4

  integer ( kind = 4 ) irep
  character ( len = 31 ) s
  character ( len = 31 ), dimension ( test_num ) :: s_test = (/ &
    'A man, a plan, a canal, Panama!', &
    'A man, a plan, a canal, Panama!', &
    'A man, a plan, a canal, Panama!', &
    'aaaaannnnnQ!                   ' /)
  character ( len = 5 ) sub
  character ( len = 5 ), dimension ( test_num ) :: sub_test = (/ &
    ',    ', &
    'an   ', &
    'canal', &
    'an   ' /)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST121'
  write ( *, '(a)' ) '  S_S_DELETE2 recursively removes a substring;'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    s = s_test(test)
    sub = sub_test(test)

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) &
      '  Remove "' // trim ( sub ) // '" from "' // trim ( s ) // '"'
    write ( *, '(a)' ) ' '

    call s_s_delete2 ( s, trim ( sub ), irep )

    write ( *, '(a)' ) '  Result: ' // trim ( s )
    write ( *, '(2x,i8,a)' ) irep, ' removals'

  end do

  return
end
subroutine test122

!*****************************************************************************80
!
!! TEST122 tests S_S_INSERT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) ipos
  character ( len = 40 ) s
  character ( len = 4 ) s2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST122'
  write ( *, '(a)' ) '  S_S_INSERT inserts one string into another.'
  write ( *, '(a)' ) ' '

  s = 'Where are the snows of yesteryear?'
  s2 = 'plow'
  ipos = 19

  write ( *, '(a,i8,a)' ) '  Insert ''' // trim ( s2 ) // '" into position ', &
    ipos, ' of '
  write ( *, '(a)' ) '  "' // trim ( s ) // '"'

  call s_s_insert ( s, ipos, s2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Result:'
  write ( *, '(a)' ) '  "' // trim ( s ) // '"'

  return
end
subroutine test1225

!*****************************************************************************80
!
!! TEST1225 tests S_S_SUBANAGRAM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 4

  logical s_s_subanagram
  character ( len = 13 ) s1
  character ( len = 13 ), dimension ( test_num ) :: s1_test = (/ &
    'Get a priest!', &
    'Get a priest!', &
    'Get a priest!', &
    'Get a priest!' /)
  character ( len = 6 ) s2
  character ( len = 6 ), dimension ( test_num ) :: s2_test = (/ &
    'stripe', &
    'pastor', &
    'a sip ', &
    'tag!  ' /)
  integer ( kind = 4 ) test
  logical value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1225'
  write ( *, '(a)' ) '  S_S_SUBANAGRAM is TRUE if S2 is a "subanagram"'
  write ( *, '(a)' ) '  of S1.'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    s1 = s1_test(test)
    s2 = s2_test(test)

    value = s_s_subanagram ( trim ( s1 ), trim ( s2 ) )

    write ( *, '(2x,a,2x,a,2x,l1)' ) &
      '"' // s1_test(test) // '"', '"' // trim ( s2_test(test) ) // '"', value

  end do

  return
end
subroutine test123

!*****************************************************************************80
!
!! TEST123 tests S_SET_DELETE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 20 ) s
  character ( len = 10 ) s2

  s2 = '0123456789'
  s = '1 way 4 U 2 deb8'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST123'
  write ( *, '(a)' ) '  S_SET_DELETE removes all occurrences of a set'
  write ( *, '(a)' ) '    of characters.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  String: "' // trim ( s ) // '"'
  write ( *, '(a)' ) '  Set:    "' // trim ( s2 ) // '"'

  call s_set_delete ( s, s2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Result:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  "' // trim ( s ) // '"'

  return
end
subroutine test124

!*****************************************************************************80
!
!! TEST124 tests S_SHIFT_CIRCULAR, S_SHIFT_LEFT and S_SHIFT_RIGHT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) ishft
  character ( len = 6 ) string
  character ( len = 6 ) string1
  character ( len = 6 ) string2
  character ( len = 6 ) string3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST124'
  write ( *, '(a)' ) '  S_SHIFT_CIRCULAR, right circular shift.'
  write ( *, '(a)' ) '  S_SHIFT_LEFT, left shift, blank pad.'
  write ( *, '(a)' ) '  S_SHIFT_RIGHT, right shift, blank pad.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  String   Shift  Shift_Circular   Shift_Right  Shift_Left'
  write ( *, '(a)' ) ' '

  ishft = 2
  string = 'Abcde '

  string1 = string
  call s_shift_circular ( string1, ishft )

  string2 = string
  call s_shift_right ( string2, ishft )

  string3 = string
  call s_shift_left ( string3, ishft )

  write ( *, '(2x,a6,2x,i8,2x,a6,2x,a6,2x,a6)' ) &
    string, ishft, string1, string2, string3

  ishft = 3
  string = '123456'

  string1 = string
  call s_shift_circular ( string1, ishft )

  string2 = string
  call s_shift_right ( string2, ishft )

  string3 = string
  call s_shift_left ( string3, ishft )

  write ( *, '(2x,a6,2x,i8,2x,a6,2x,a6,2x,a6)' ) &
    string, ishft, string1, string2, string3

  ishft = -2
  string = 'Shazam'

  string1 = string
  call s_shift_circular ( string1, ishft )

  string2 = string
  call s_shift_right ( string2, ishft )

  string3 = string
  call s_shift_left ( string3, ishft )

  write ( *, '(2x,a6,2x,i8,2x,a6,2x,a6,2x,a6)' ) &
    string, ishft, string1, string2, string3

  return
end
subroutine test125

!*****************************************************************************80
!
!! TEST125 tests S_SKIP_SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character blank
  character hat
  integer ( kind = 4 ) i
  integer ( kind = 4 ) locnew
  integer ( kind = 4 ) locold
  character ( len = 20 ) s
  character ( len = 10 ) s2
  integer ( kind = 4 ) s_skip_set

  blank = ' '
  hat = '^'
  s2 = '0123456789'
  s = '1 way 4 U 2 deb8!'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST125'
  write ( *, '(a)' ) '  S_SKIP_SET finds the next character that '
  write ( *, '(a)' ) '    IS NOT part of a given set of characters;'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Our string is'
  write ( *, '(a)' ) '  "' // trim ( s ) // '"'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Our character set is'
  write ( *, '(a)' ) '  "' // trim ( s2 ) // '"'

  locnew = 0

  do

    locold = locnew

    locnew = s_skip_set ( s(locold+1:), s2 ) + locold

    if ( locnew == locold ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  No more matches.'
      exit
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  "' // trim ( s ) // '"'
    write ( *, '(40a)' ) ( blank, i = 1, locnew-1 ), hat

  end do

  return
end
subroutine test1255

!*****************************************************************************80
!
!! TEST1255 tests S_SORT_A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 5

  character ( len = 20 ) s
  character ( len = 20 ), dimension ( test_num ) :: s_test = (/ &
    'HELLO World   !! !  ', &
    '12345678901234567890', &
    'Abc Def Ghi Jkl Mno ', &
    'AbleBakerCharlieDelt', &
    'What? You have seen?' /)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1255'
  write ( *, '(a)' ) '  S_SORT_A ascending sorts a string.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  -------String-------  -------Sorted-------'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    s = s_test(test)
    call s_sort_a ( s )
    write ( *, '(2x,a22,2x,a22)' ) &
      '"' // s_test(test) // '"', '"' // s // '"'
  end do

  return
end
subroutine test126

!*****************************************************************************80
!
!! TEST126 tests S_SPLIT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 4

  character ( len = 250 ) output
  character ( len = 50 ) s
  character ( len = 50 ) s_test(test_num)
  character ( len = 50 ) s1
  character ( len = 50 ) s2
  character ( len = 50 ) s3
  character ( len = 50 ) sub
  character ( len = 50 ) sub_test(test_num)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST126'
  write ( *, '(a)' ) '  S_SPLIT splits a string at a substring.'
  write ( *, '(a)' ) ' '

  s_test(1) = '      REAL FUNCTION GRAMMA ( X, Y, Z )'
  sub_test(1) = 'real function'

  s_test(2) = '      real  function   gramma ( x, y, z )'
  sub_test(2) = 'real function'

  s_test(3) = '      REAL FUNCTION GRAMMA ( X, Y, Z )'
  sub_test(3) = 'unc'

  s_test(4) = '      real  function   gramma ( x, y, z )'
  sub_test(4) = 'lemon'

  do test = 1, test_num

    s = s_test(test)
    sub = sub_test(test)

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  String = "' // trim ( s ) // '"'
    write ( *, '(a)' ) '  Substring = '// trim ( sub )

    call s_split ( s, sub, s1, s2, s3 )

    if ( s2 == ' ' ) then
      write ( *, '(a)' ) '  No match'
    else
      output = s1 // ' // ' // s2 // ' // ' // s3
      call s_blanks_delete ( output )
      write ( *, '(2x,a)' ) trim ( output )
    end if

  end do

  return
end
subroutine test127

!*****************************************************************************80
!
!! TEST127 tests S_TAB_BLANKS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 4
  character, parameter :: TAB = char ( 9 )

  character ( len = 80 ) s
  character ( len = 80 ) s_test(test_num)
  integer ( kind = 4 ) test

  s_test(1) = 'No tabs in me.'
  s_test(2) = 'I''ve got one' // TAB // 'tab here!'
  s_test(3) = 'I' // TAB // 'have' // TAB // 'three' // TAB // '!'
  s_test(4) = TAB // 'I begin and end with them!' // TAB

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST127'
  write ( *, '(a)' ) '  S_TAB_BLANKS replaces TAB''s by 6 spaces.'

  do test = 1, test_num
    write ( *, '(a)' ) ' '
    s = s_test(test)
    write ( *, '(2x,a)' ) '  "' // trim ( s ) // '"'
    call s_tab_blanks ( s )
    write ( *, '(2x,a)' ) '  "' // trim ( s ) // '"'
  end do

  return
end
subroutine test128

!*****************************************************************************80
!
!! TEST128 tests S_TO_C4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 8

  complex ( kind = 4 ) cval
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) length
  character ( len = 10 ), dimension ( test_num ) :: s_test = (/ &
    '1         ', &
    '2+I       ', &
    '3 + 4 I   ', &
    '5 + 6*I   ', &
    'I         ', &
    '7 I       ', &
    '-8 * I    ', &
    '44 * 99   ' /)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST128'
  write ( *, '(a)' ) '  S_TO_C4 accepts a string of characters'
  write ( *, '(a)' ) '  and extracts a complex value from them,'
  write ( *, '(a)' ) '  assuming a format of A+BI for complex values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  String              CVAL    IERROR   LENGTH'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    call s_to_c4 ( s_test(test), cval, ierror, length )

    write ( *, '(2x,a10,2x,2f8.1,2x,i2,6x,i2)' ) &
      s_test(test), cval, ierror, length

  end do

  return
end
subroutine test129

!*****************************************************************************80
!
!! TEST129 tests S_TO_FORMAT
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 7

  character c
  integer ( kind = 4 ) m
  integer ( kind = 4 ) r
  character ( len = 20 ) s_test(test_num)
  integer ( kind = 4 ) test
  integer ( kind = 4 ) w

  s_test(1) = 'a80'
  s_test(2) = 'f8.4'
  s_test(3) = '3g14.6'
  s_test(4) = 'i12'
  s_test(5) = '12l1'
  s_test(6) = '(10o11)'
  s_test(7) = ' ( 5 z 11.7  )'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST129'
  write ( *, '(a)' ) '  S_TO_FORMAT, string -> FORTRAN format RcW.M;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  --------String------     R  c     W      M'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    call s_to_format ( s_test(test), r, c, w, m )

    write ( *, '(2x,a20,i8,2x,a1,i8,i8)' ) s_test(test), r, c, w, m

  end do

  return
end
subroutine test130

!*****************************************************************************80
!
!! TEST130 tests S_TO_L.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 December 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 10

  logical l
  character ( len = 10 ) s
  logical s_to_l
  character ( len = 10 ) string(test_num)
  integer ( kind = 4 ) test

  string(1) = '0'
  string(2) = 'F'
  string(3) = 'f'
  string(4) = '1'
  string(5) = 'T'
  string(6) = 't'
  string(7) = '  0'
  string(8) = '  1  0'
  string(9) = '  01'
  string(10) = '  Talse'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST130'
  write ( *, '(a)' ) '  S_TO_L reads logical data from a string.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  S   L'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    s = string(test)
    l = s_to_l ( s )
    write ( *, '(2x,a10,2x,l1,4x,i2,4x,i2)' ) s, l
  end do

  return
end
subroutine test131

!*****************************************************************************80
!
!! TEST131 tests S_TO_R4VEC
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3
  integer ( kind = 4 ), parameter :: test_num = 3

  integer ( kind = 4 ) ierror
  real    ( kind = 4 ) rvec(n)
  character ( len = 20 ) s_test(test_num)
  character ( len = 20 ) s
  integer ( kind = 4 ) test

  s_test(1) = ' 1 2 3'
  s_test(2) = '1.5 2.25 3.75'
  s_test(3) = '10, 21.0, 32.0, 43.0'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST131'
  write ( *, '(a)' ) '  S_TO_R4VEC, string -> R4VEC;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  --------String------   R(1)      R(2)      R(3)'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    s = s_test(test)
    call s_to_r4vec ( s, n, rvec, ierror )

    write ( *, '(2x,a20,3f10.4)' ) s, rvec(1:n)

  end do

  return
end
subroutine test132

!*****************************************************************************80
!
!! TEST132 tests S_TO_ROT13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 3

  character ( len = 30 ), dimension ( test_num ) :: s_test = (/ &
    'abcdefghijklmnopqrstuvwxyz    ', &
    'Cher                          ', &
    'James Thurston Howell         ' /)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST132'
  write ( *, '(a)' ) '  S_TO_ROT13 encrypts a string.'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    write ( *, '(a)' ) ' '
    write ( *, '(a,a30)' )'  Original:      ', trim ( s_test(test) )
    call s_to_rot13 ( s_test(test) )
    write ( *, '(a,a30)' )'  Rotated once:  ', trim ( s_test(test) )
    call s_to_rot13 ( s_test(test) )
    write ( *, '(a,a30)' )'  Rotated twice: ', trim ( s_test(test) )
  end do

  return
end
subroutine test133

!*****************************************************************************80
!
!! TEST133 tests S_TO_SOUNDEX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 14

  character ( len = 4 ) code
  character ( len = 15 ) s_test(test_num)
  integer ( kind = 4 ) test

  s_test(1) = 'Ellery'
  s_test(2) = 'Euler'
  s_test(3) = 'Gauss'
  s_test(4) = 'Ghosh'
  s_test(5) = 'Heilbronn'
  s_test(6) = 'hi-lo-ball'
  s_test(7) = 'Hilbert'
  s_test(8) = 'Kant'
  s_test(9) = 'Knuth'
  s_test(10) = 'Ladd'
  s_test(11) = 'Lloyd'
  s_test(12) = 'Lissajous'
  s_test(13) = 'Lukasiewicz'
  s_test(14) = 'Bob'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST133'
  write ( *, '(a)' ) '  S_TO_SOUNDEX converts a string to a Soundex code.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' Test      String      Code'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    call s_to_soundex ( s_test(test), code )
    write ( *, '(2x,i3,2x,a15,2x,a4)' ) test, s_test(test), code
  end do

  return
end
subroutine test134

!*****************************************************************************80
!
!! TEST134 tests S_TO_W.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 3

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) length
  integer ( kind = 4 ) j
  character ( len = 30 ) s_test(test_num)
  integer ( kind = 4 ) test
  character ( len = 30 ) w

  s_test(1) = 'This is simple.'
  s_test(2) = '  HERE''s a  har*der_one!'
  s_test(3) = '  what       now?'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST134'
  write ( *, '(a)' ) '  S_TO_W accepts a string of characters'
  write ( *, '(a)' ) '  and extracts blank-delimited words.'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Test string S = "' // trim ( s_test(test) ) // '"'
    write ( *, '(a)' ) ' '

    j = 1

    do

      call s_to_w ( s_test(test), w, ierror, length )

      if ( ierror /= 0 ) then
        exit
      end if

      write ( *, '(2x,i3,2x,a)' ) j, trim ( w )

      s_test(test) = s_test(test)(length+1:)

    end do

  end do

  return
end
subroutine test135

!*****************************************************************************80
!
!! TEST135 tests S_TOKEN_EQUAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nset = 5

  integer ( kind = 4 ) i
  integer ( kind = 4 ) iset
  character ( len = 10 ) set(nset)
  character ( len = 10 ) s

  set(1) = 'Bob'
  set(2) = 'Archie'
  set(3) = 'Veronica'
  set(4) = 'Jughead'
  set(5) = 'Betty'

  s = 'verONICa'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST135'
  write ( *, '(a)' ) '  S_TOKEN_EQUAL searches for whether'
  write ( *, '(a)' ) '  a string is in a set.  Here, the string is'
  write ( *, '(a)' ) '  "' // trim ( s ) // '"'
  write ( *, '(a)' ) '  and the set is'

  do i = 1, nset
    write ( *, '(2x,a)' ) trim ( set(i) )
  end do

  call s_token_equal ( s, set, nset, iset )

  write ( *, '(a)' ) ' '

  if ( iset /= 0 ) then
    write ( *, '(a)' ) '  The matching entry is ' // trim ( set(iset) )
  else
    write ( *, '(a)' ) '  No match.'
  end if

  return
end
subroutine test136

!*****************************************************************************80
!
!! TEST136 tests S_TOKEN_MATCH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: token_num = 4

  integer ( kind = 4 ) match
  character ( len = 40 ) s
  character ( len = 20 ) token(token_num)
  integer ( kind = 4 ) token_i

  s = 'TommyGun'

  token(1) = 'Tom'
  token(2) = 'Zebra'
  token(3) = 'TommY'
  token(4) = 'TommyKnocker'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST136'
  write ( *, '(a)' ) '  S_TOKEN_MATCH finds longest token match.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' Our string is'
  write ( *, '(a)' ) '  "' // trim ( s ) // '"'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Our tokens are:'
  write ( *, '(a)' ) ' '
  do token_i = 1, token_num
    write ( *, '(2x,a)' ) token(token_i)
  end do

  call s_token_match ( s, token_num, token, match )

  write ( *, '(a)' ) ' '
  if ( match == 0 ) then
    write ( *, '(a)' ) '  No matching token was found.'
  else
    write ( *, '(a,i8)' ) '  Maximum match occurs with token ', match
  end if

  return
end
subroutine test137

!*****************************************************************************80
!
!! TEST137 tests S_WORD_COUNT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 4

  integer ( kind = 4 ) nword
  character ( len = 32 ) s_test(test_num)
  integer ( kind = 4 ) test

  s_test(1) = '?'
  s_test(2) = 'A man, a plan, a canal - Panama!'
  s_test(3) = ' justone!word,-@#$ '
  s_test(4) = 'How about a day in the park?'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST137'
  write ( *, '(a)' ) '  S_WORD_COUNT counts the words in a string'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  STRING                                    Words'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    call s_word_count ( s_test(test), nword )

    write ( *, '( 2x, a32, 2x, i12 )' ) s_test(test), nword

  end do

  return
end
subroutine test138

!*****************************************************************************80
!
!! TEST138 tests S_WORD_EXTRACT_FIRST.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 80 ) s
  character ( len = 80 ) word

  s = 'Just an incontrovertible sample of text!'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST138'
  write ( *, '(a)' ) '  S_WORD_EXTRACT_FIRST extracts the first word'
  write ( *, '(a)' ) '  from a string.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Our input string is:'
  write ( *, '(a)' ) '  "' // trim ( s ) // '"'
  write ( *, '(a)' ) ' '

  do

    call s_word_extract_first ( s, word )

    if ( len_trim ( word ) <= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Reached the last word.'
      exit
    end if

    write ( *, '(a)' ) '  "' // trim ( word ) // '"'

  end do

  return
end
subroutine test139

!*****************************************************************************80
!
!! TEST139 tests S_WORD_FIND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) iword
  integer ( kind = 4 ) nword
  character ( len = 30 ) s
  character ( len = 10 ) word

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST139'
  write ( *, '(a)' ) '  S_WORD_FIND looks for a particular word in a string.'
  write ( *, '(a)' ) ' '

  s = 'Fred is following me around!'
  write ( *, '(a)' ) '  string = ' // s
  iword = 4
  write ( *, '(a,i8)' ) '  We want to find word number ', iword

  call s_word_find ( s, iword, word, nword )

  if ( nword == 0 ) then
    write ( *, '(a)' ) '  S_WORD_FIND could not find the requested word.'
  else
    write ( *, '(a,i8)' ) '  Word has length ', nword
    write ( *, '(a)' ) '  The requested word is ' // trim ( word )
  end if

  return
end
subroutine test140

!*****************************************************************************80
!
!! TEST140 tests S_WORD_INDEX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) iword
  character ( len = 30 ) s

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST140'
  write ( *, '(a)' ) '  S_WORD_INDEX finds the Nth word in a string.'
  write ( *, '(a)' ) ' '

  s = 'Fred is following me around!'
  write ( *, '(a)' ) '  String = ' // trim ( s )
  iword = 4
  write ( *, '(a,i8)' ) '  We want to find word number ', iword

  call s_word_index ( s, iword, ilo, ihi )

  if ( ilo == 0 .and. ihi == 0 ) then
    write ( *, '(a)' ) '  S_WORD_INDEX could not find the requested word.'
  else
    write ( *, '(a,i8,a,i8)' ) '  Word lies between locations ', ilo, &
      ' and ', ihi
    write ( *, '(a)' ) '  The requested word is ' // s(ilo:ihi)
  end if

  return
end
subroutine test141

!*****************************************************************************80
!
!! TEST141 tests S_WORD_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  logical done
  integer ( kind = 4 ) i
  character ( len = 80 ) s
  character ( len = 80 ) word

  s = 'Just an incontrovertible sample of text!'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST141'
  write ( *, '(a)' ) '  S_WORD_NEXT returns each word '
  write ( *, '(a)' ) '  in order, from a string.'

  do i = 1, 2

    if ( i == 1 ) then
      s = 'Just, an incontrovertible (sample of) text!'
    else
      s = 'A "second" string.'
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Input string:'
    write ( *, '(a)' ) '  "' // trim ( s ) // '"'

    done = .true.

    do

      call s_word_next ( s, word, done )

      if ( done ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  No more words in the string.'
        exit
      end if

      write ( *, '(4x,a)' ) trim ( word )

    end do

  end do

  return
end
subroutine test142

!*****************************************************************************80
!
!! TEST142 tests S_WORD_PERMUTE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 6

  integer ( kind = 4 ), dimension ( n ) :: perm = (/ 3, 4, 1, 2, 6, 5 /)
  character ( len = 80 ) s1
  character ( len = 80 ) s2

  s1 = 'Just an incontrovertible sample of text!'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST142'
  write ( *, '(a)' ) '  S_WORD_PERMUTE permutes the words in a string.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The permutation is:'
  write ( *, '(a)' ) ' '
  write ( *, '(6(2x,i2))' ) (/ 1, 2, 3, 4, 5, 6 /)
  write ( *, '(6(2x,i2))' ) perm(1:n)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Our input string is:'
  write ( *, '(a)' ) '  "' // trim ( s1 ) // '"'
  write ( *, '(a)' ) ' '

  call s_word_permute ( s1, n, perm, s2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Our output string is:'
  write ( *, '(a)' ) '  "' // trim ( s2 ) // '"'

  return
end
subroutine test143

!*****************************************************************************80
!
!! TEST143 tests SVEC_LAB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 15

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ident(n)
  integer ( kind = 4 ) nuniq
  character ( len = 20 ) string(n)

  string(1) = 'ALPHA'
  string(2) = 'BETA'
  string(3) = ' '
  string(4) = 'ALPHA'
  string(5) = 'Alpha'
  string(6) = 'GAMMA'
  string(7) = 'BETA'
  string(8) = 'BETA'
  string(9) = 'ALPHA'
  string(10) = 'GAMMA'
  string(11) = ' '
  string(12) = ' '
  string(13) = 'RHO'
  string(14) = 'EPSILON'
  string(15) = 'Alpha'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST143'
  write ( *, '(a)' ) '  SVEC_LAB marks unique strings in a list.'

  call svec_lab ( n, nuniq, string, ident )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of unique entries = ', nuniq
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  String, ID'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,a20,2x,i8)' ) string(i), ident(i)
  end do

  return
end
subroutine test144

!*****************************************************************************80
!
!! TEST144 tests SVEC_MERGE_A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: na = 10
  integer ( kind = 4 ), parameter :: nb = 10

  character ( len = 4 ) a(na)
  character ( len = 4 ) b(nb)
  character ( len = 4 ) c(na+nb)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nc

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST144'
  write ( *, '(a)' ) '  SVEC_MERGE_A merges two sorted character arrays.'
  write ( *, '(a)' ) ' '

  a(1) = 'Adam'
  a(2) = 'Bill'
  a(3) = 'Bob'
  a(4) = 'Carl'
  a(5) = 'Carl'
  a(6) = 'Earl'
  a(7) = 'Fred'
  a(8) = 'Jean'
  a(9) = 'Lynn'
  a(10) = 'Zeke'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Input vector A:'
  write ( *, '(a)' ) ' '

  do i = 1, na
    write ( *, '(2x,a4)' ) a(i)
  end do

  b(1) = 'Ada'
  b(2) = 'Barb'
  b(3) = 'Cath'
  b(4) = 'Deb'
  b(5) = 'Eve'
  b(6) = 'Inez'
  b(7) = 'Jane'
  b(8) = 'Jean'
  b(9) = 'Jill'
  b(10) = 'Lynn'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Input vector B:'
  write ( *, '(a)' ) ' '

  do i = 1, nb
    write ( *, '(2x,a4)' ) b(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call SVEC_MERGE_A to merge the two lists.'

  call svec_merge_a ( na, a, nb, b, nc, c )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Merged output vector C:'
  write ( *, '(a)' ) ' '

  do i = 1, nc
    write ( *, '(2x,a4)' ) c(i)
  end do

  return
end
subroutine test145

!*****************************************************************************80
!
!! TEST145 tests SVEC_REVERSE and SVEC_SORT_HEAP_A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  character ( len = 10 ) carray(n)
  integer ( kind = 4 ) i

  carray(1) = 'FRED'
  carray(2) = 'fred'
  carray(3) = 'Abacus'
  carray(4) = 'beetles'
  carray(5) = 'XYLOPHONE'
  carray(6) = 'banana'
  carray(7) = 'goofball'
  carray(8) = 'abbot'
  carray(9) = 'BARBECUE'
  carray(10) = 'abbots'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST145'
  write ( *, '(a)' ) '  SVEC_SORT_HEAP_A sorts a string vector.'
  write ( *, '(a)' ) '  SVEC_REVERSE reverses a string vector.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Unsorted list:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(4x,a)' ) carray(i)
  end do

  call svec_sort_heap_a ( n, carray )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sorted list:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(4x,a)' ) carray(i)
  end do

  call svec_reverse ( n, carray )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Reversed sorted list:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(4x,a)' ) carray(i)
  end do

  return
end
subroutine test146

!*****************************************************************************80
!
!! TEST146 tests SVEC_SORT_HEAP_A, SVEC_MERGE_A and SVEC_SEARCH_BINARY_A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: na = 10
  integer ( kind = 4 ), parameter :: nb = 10

  character ( len = 4 ) a(na)
  character ( len = 4 ) b(nb)
  character ( len = 4 ) c(na+nb)
  character ch_uniform
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nc
  integer ( kind = 4 ) seed
  character ( len = 4 ) string

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST146'
  write ( *, '(a)' ) '  For ascending order:'
  write ( *, '(a)' ) '  SVEC_SORT_HEAP_A sorts a character array;'
  write ( *, '(a)' ) '  SVEC_MERGE_A merges two sorted character '
  write ( *, '(a)' ) '    arrays into a single sorted array.'
  write ( *, '(a)' ) '  SVEC_SEARCH_BINARY_A searches a string array for'
  write ( *, '(a)' ) '    a particular value.'
  write ( *, '(a)' ) ' '

  seed = 123456789

  do i = 1, na
    do j = 1, 4
      a(i)(j:j) = ch_uniform ( 'A', 'E', seed )
    end do
  end do

  call svec_sort_heap_a ( na, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sorted vector A:'
  write ( *, '(a)' ) ' '

  do i = 1, na
    write ( *, '(4x,a)' ) a(i)
  end do

  do i = 1, nb
    do j = 1, 4
      b(i)(j:j) = ch_uniform ( 'B', 'F', seed )
    end do
  end do

  call svec_sort_heap_a ( nb, b )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sorted vector B:'
  write ( *, '(a)' ) ' '

  do i = 1, nb
    write ( *, '(4x,a)' ) b(i)
  end do

  call svec_merge_a ( na, a, nb, b, nc, c )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Merged output vector C = A + B:'
  write ( *, '(a)' ) ' '

  do i = 1, nc
    write ( *, '(4x,a)' ) c(i)
  end do

  string = a(2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Search C for value ' // string
  write ( *, '(a)' ) ' '

  call svec_search_binary_a ( nc, c, string, indx )

  if ( indx == 0 ) then
    write ( *, '(a)' ) '  The value does not occur'
  else
    write ( *, '(a,i8)' ) '  The value occurs at index ', indx
  end if

  return
end
subroutine test147

!*****************************************************************************80
!
!! TEST147 tests SVEC_SORT_HEAP_A and SVEC_SORTED_UNIQUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  character ( len = 3 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nuniq

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST147'
  write ( *, '(a)' ) '  SVEC_SORT_HEAP_A sorts a character array;'
  write ( *, '(a)' ) '  SVEC_SORTED_UNIQUE finds the unique entries in a'
  write ( *, '(a)' ) '    sorted character array.'
  write ( *, '(a)' ) ' '

  a(1) = 'Cat'
  a(2) = 'Bat'
  a(3) = 'Mat'
  a(4) = 'Tab'
  a(5) = 'Ax'
  a(6) = 'Ax'
  a(7) = 'Tab'
  a(8) = 'Pyx'
  a(9) = 'Ax'
  a(10) = 'Bat'

  call svec_sort_heap_a ( n, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Input vector A:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(4x,a)' ) a(i)
  end do

  call svec_sorted_unique ( n, a, nuniq )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Unique entries:'
  write ( *, '(a)' ) ' '

  do i = 1, nuniq
    write ( *, '(4x,a)' ) a(i)
  end do

  return
end
subroutine test148

!*****************************************************************************80
!
!! TEST148 tests SVEC_SORT_HEAP_A and SVECI_SORT_HEAP_A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 14

  character ( len = 10 ) svec(n)
  integer ( kind = 4 ) i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST148'
  write ( *, '(a)' ) '  Sort an array of character strings:'
  write ( *, '(a)' ) '  SVEC_SORT_HEAP_A, case-sensitive;'
  write ( *, '(a)' ) '  SVECI_SORT_HEAP_A, case-insensitive.'
  write ( *, '(a)' ) ' '

  svec(1) = 'FRED'
  svec(2) = 'fred'
  svec(3) = 'Abacus'
  svec(4) = 'beetles'
  svec(5) = 'XYLOPHONE'
  svec(6) = 'banana'
  svec(7) = 'goofball'
  svec(8) = 'abbot'
  svec(9) = 'BARBECUE'
  svec(10) = 'abbots'
  svec(11) = ' indented'
  svec(12) = '123456'
  svec(13) = 'beetles'
  svec(14) = 'Abacus'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Unsorted list:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,a)' ) svec(i)
  end do

  call svec_sort_heap_a ( n, svec )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sorted list:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,a)' ) svec(i)
  end do

  svec(1) = 'FRED'
  svec(2) = 'fred'
  svec(3) = 'Abacus'
  svec(4) = 'beetles'
  svec(5) = 'XYLOPHONE'
  svec(6) = 'banana'
  svec(7) = 'goofball'
  svec(8) = 'abbot'
  svec(9) = 'BARBECUE'
  svec(10) = 'abbots'
  svec(11) = ' indented'
  svec(12) = '123456'
  svec(13) = 'beetles'
  svec(14) = 'Abacus'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now do a case-insensitive sort:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Unsorted list:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,a)' ) svec(i)
  end do

  call sveci_sort_heap_a ( n, svec )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sorted list:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,a)' ) svec(i)
  end do

  return
end
subroutine test149

!*****************************************************************************80
!
!! TEST149 tests SVEC_SORT_HEAP_A_INDEX and SVECI_SORT_HEAP_A_INDEX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 14

  character ( len = 10 ) carray(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST149'
  write ( *, '(a)' ) '  Indexed heap sort of strings:'
  write ( *, '(a)' ) '  SVEC_SORT_HEAP_A_INDEX, case-sensitive;'
  write ( *, '(a)' ) '  SVECI_SORT_HEAP_A_INDEX, case-insensitive.'
  write ( *, '(a)' ) ' '

  carray(1) = 'FRED'
  carray(2) = 'fred'
  carray(3) = 'Abacus'
  carray(4) = 'beetles'
  carray(5) = 'XYLOPHONE'
  carray(6) = 'banana'
  carray(7) = 'goofball'
  carray(8) = 'abbot'
  carray(9) = 'BARBECUE'
  carray(10) = 'abbots'
  carray(11) = ' indented'
  carray(12) = '123456'
  carray(13) = 'beetles'
  carray(14) = 'Abacus'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Unsorted list:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,a)' ) carray(i)
  end do

  call svec_sort_heap_a_index ( n, carray, indx )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sorted list:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,a)' ) carray(indx(i))
  end do

  carray(1) = 'FRED'
  carray(2) = 'fred'
  carray(3) = 'Abacus'
  carray(4) = 'beetles'
  carray(5) = 'XYLOPHONE'
  carray(6) = 'banana'
  carray(7) = 'goofball'
  carray(8) = 'abbot'
  carray(9) = 'BARBECUE'
  carray(10) = 'abbots'
  carray(11) = ' indented'
  carray(12) = '123456'
  carray(13) = 'beetles'
  carray(14) = 'Abacus'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now do a case-insensitive sort:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Unsorted list:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,a)' ) carray(i)
  end do

  call sveci_sort_heap_a_index ( n, carray, indx )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sorted list:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,a)' ) carray(indx(i))
  end do

  return
end
subroutine test150

!*****************************************************************************80
!
!! TEST150 tests SVECI_SEARCH_BINARY_A and SVECI_SORT_HEAP_A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  character ( len = 10 ) carray(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  character ( len = 10 ) string

  carray(1) = 'FRED'
  carray(2) = 'fred'
  carray(3) = 'Abacus'
  carray(4) = 'beetles'
  carray(5) = 'XYLOPHONE'
  carray(6) = 'banana'
  carray(7) = 'goofball'
  carray(8) = 'abbot'
  carray(9) = 'BARBECUE'
  carray(10) = 'abbots'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST150'
  write ( *, '(a)') '  For implicitly capitalized strings,'
  write ( *, '(a)' ) '  SVECI_SORT_HEAP_A sorts;'
  write ( *, '(a)' ) '  SVECI_SEARCH_BINARY_A searches.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Unsorted list:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(4x,a)' ) carray(i)
  end do

  call sveci_sort_heap_a ( n, carray )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sorted list:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(4x,a)' ) trim ( carray(i) )
  end do

  string = 'ABBoT'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now search for the string ' // trim ( string )

  call sveci_search_binary_a ( n, carray, string, indx )

  write ( *, '(a)' ) ' '
  if ( indx == 0 ) then
    write ( *, '(a)' ) '  The search string does not occur.'
  else
    write ( *, '(a,i8)' ) '  The search string occurs in index ', indx
  end if

  return
end
subroutine test152

!*****************************************************************************80
!
!! TEST152 tests WORD_LAST_READ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 80 ) s
  integer ( kind = 4 ) test
  character ( len = 80 ) word

  s = 'Just an incontrovertible sample of text!'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST152'
  write ( *, '(a)' ) '  WORD_LAST_READ returns the last word from a string.'

  do test = 1, 2

    if ( test == 1 ) then
      s = 'Just, an incontrovertible (sample of) text!'
    else
      s = 'A "second" string.'
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Input string:'
    write ( *, '(a)' ) '  "' // trim ( s ) // '"'

    call word_last_read ( s, word )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Last word:  "' // trim ( word ) // '"'

  end do

  return
end
subroutine test153

!*****************************************************************************80
!
!! TEST153 tests WORD_NEXT.
!
!  Discussion:
!
!    Thanks to Bill Richmond for pointing out that ILO and IHI must
!    be initialized to 0, in order for this routine to work properly!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  character ( len = 80 ) s

  s = '  Just an incontrovertible ,   sample of text!'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST153'
  write ( *, '(a)' ) '  WORD_NEXT returns each "word" from a string.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We use the following string:'
  write ( *, '(a)' ) '    "' // trim ( s ) // '"'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here are the individual words:'
  write ( *, '(a)' ) ' '

  ilo = 0
  ihi = 0

  do

    call word_next ( s, ilo, ihi )

    if ( ilo <= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  There are no more words in the string.'
      exit
    end if

    write ( *, '(2x,a)' ) '"' // s(ilo:ihi) // '"'

  end do

  return
end
subroutine test154

!*****************************************************************************80
!
!! TEST154 tests WORD_NEXT_READ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  logical done
  character ( len = 80 ) s
  character ( len = 80 ) w
  integer ( kind = 4 ) word_num

  s = '  Here is a string, (you see) with x[1] = {gamma}!'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST154'
  write ( *, '(a)' ) '  WORD_NEXT_READ returns each "word" from a string.'
  write ( *, '(a)' ) '  It pays attention to various parentheses and brackets.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We use the following string:'
  write ( *, '(a)' ) '    "' // trim ( s ) // '"'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here are the individual words:'
  write ( *, '(a)' ) ' '

  done = .true.
  word_num = 0

  do

    call word_next_read ( s, w, done )

    if ( done ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of words was ', word_num
      exit
    end if

    word_num = word_num + 1

    write ( *, '(2x,i8,2x,a)' ) word_num, '"' // trim ( w ) // '"'

  end do

  return
end
subroutine test155

!*****************************************************************************80
!
!! TEST155 tests WORD_NEXT2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 80 ) first
  character ( len = 80 ) s
  character ( len = 80 ) last

  s = 'Just an incontrovertible sample of text!'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST155'
  write ( *, '(a)' ) '  WORD_NEXT2 returns each word from a string.'
  write ( *, '(a)' ) ' '
  write ( *, '(4x,a)' ) s
  write ( *, '(a)' ) ' '

  do

    call word_next2 ( s, first, last )

    if ( len_trim ( first ) <= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(2x,a)' ) 'Reached the last word.'
      exit
    end if

    write ( *, '(4x,a)' ) trim ( first )
    s = last

  end do

  return
end
subroutine test156

!*****************************************************************************80
!
!! TEST156 tests WORD_SWAP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) iword1
  integer ( kind = 4 ) iword2
  character ( len = 80 ) line

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST156'
  write ( *, '(a)' ) '  WORD_SWAP swaps two words in a string'

  line = 'This is the true story of six roommates who '
  write ( *, '(4x,a)' ) trim ( line )

  iword1 = 4
  iword2 = 8

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a,i8)' ) '  Now swap words ', iword1, ' and ', iword2

  call word_swap ( line, iword1, iword2 )

  write ( *, '(a)' ) ' '
  write ( *, '(4x,a)' ) trim ( line )

  return
end

