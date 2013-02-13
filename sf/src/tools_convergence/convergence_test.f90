  if(success > N1)convergence=.true.
  if(check>=N2)then
     open(10,file="ERROR.README")
     write(10,*)""
     close(10)
     write(*,"(A)")bg_red("Convergence not achieved after N2 iterations!! exiting..")
     convergence=.true.
  endif
  if(convergence)then
     check=0
     deallocate(Xold)
  endif
