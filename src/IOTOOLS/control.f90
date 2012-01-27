  inquire(file=reg_filename(pname),exist=control)
  if(.not.control)then
     call msg(bold_red("I can not read : +"//reg_filename(pname)//"SKIP"))
     return
  else
     write(*,"(A,A)")"read:     ",reg_filename(pname)
  endif
