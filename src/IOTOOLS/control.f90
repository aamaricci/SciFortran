  inquire(file=reg_filename(pname),exist=control)
  if(.not.control)then
     write(*,"(A,A,A)")"I can not read : ",reg_filename(pname),"SKIP"
     return
  else
     write(*,"(A,A)")"read:     ",reg_filename(pname)
  endif
