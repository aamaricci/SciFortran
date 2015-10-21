  !if(isnan(err))call abort("check_convergence: error is NaN. EXIT...")
  open(10,file=reg(file_),position="append")
  write(10,"(I5,ES15.7)")check,err
  close(10)
