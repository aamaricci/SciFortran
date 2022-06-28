program test_parsing
  USE SCIFOR, only: parse_cmd_variable,parse_input_variable
  USE ASSERTING
  implicit none

  integer                  :: L, L_r=1000
  real(8)                  :: dble,dble_r=0.5d0
  real(8),dimension(3)     :: dble_vec,dble_vec_r=[1d0,1d0,0d0]
  character(len=16)        :: string, string_r="input.conf"
  logical                  :: bool, bool_r=.false.

  !This parse from cmd line only:
  call parse_cmd_variable(string,"string",default='input.conf')
  !These parse from file and cmd line, case insesitive (as in Fortran), order irrelevant
  call parse_input_variable(dble,"Dble",string,default=1d0,comment="A dble")
  call parse_input_variable(L,"L",string,default=500) !comment are optional
  call parse_input_variable(bool,"bool",string,default=.false.,comment="A bool")
  call parse_input_variable(dble_vec,"Dble_vec",string,default=[0d0,0d0,0d0],comment="A dble array rank 1")

  call assert(string,string_r,"PARSE_STRING")
  call assert(dble,dble_r,"PARSE_DBLE")
  call assert(dble_vec,dble_vec_r,"PARSE_DBLE_ARRAY")
  call assert(L,L_r,"PARSE_INT")
  call assert(bool,bool_r,"PARSE_BOOL")

end program test_parsing
