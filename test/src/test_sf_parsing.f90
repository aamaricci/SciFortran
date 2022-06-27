program test_parsing
  use asserting
  use scifor
  implicit none
  
  integer                  :: Le, Le_r=1000
  real(8)                  :: wmixing,wmixing_r=0.5d0
  character(len=16)        :: finput, finput_r="inputAHM.conf"
  logical                  :: phsym, phsym_r=.false.

  logical                  :: info(4)


  call parse_cmd_variable(finput,"FINPUT",default='inputAHM.conf')
  call parse_input_variable(wmixing,"wmixing",finput,default=1.d0,comment="Mixing bath parameter")
  call parse_input_variable(Le,"LE",finput,default=500)
  call parse_input_variable(phsym,"phsym",finput,default=.false.,comment="Flag to enforce p-h symmetry of the bath.")

  call assert_c(finput,finput_r,info(1))
  call assert_d(wmixing,wmixing_r,info(2))
  call assert_i(Le,Le_r,info(3))
  info(4) = (phsym.eqv.phsym_r)

  if(all(info)) then
     write(*,"(A)") "test_parsing exit status is "//achar(27)//"[32m POSITIVE "//achar(27)//"[0m."
  else
     write(*,"(A)") "test_parsing exit status is "//achar(27)//"[31m NEGATIVE "//achar(27)//"[0m."
     error stop 2
  endif
  
end program test_parsing
