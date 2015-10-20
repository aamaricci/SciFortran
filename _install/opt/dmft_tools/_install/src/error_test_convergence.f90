  if(success > N1)convergence=.true.
  if(check>=N2)then
     open(10,file="ERROR.README")
     write(10,*)""
     close(10)
     write(*,"(A,I4,A)")"Not converged after",N2," iterations."
     if(extend_)then
        write(*,"(A,I4,A3,I4)")        "Increasing    N2    :",N2,"-->",N2+3*N1
        write(*,"(A,ES15.7,A3,ES15.7)")"Increasing Threshold:",eps,"-->",10.d0*eps
        !increase error_threshold by 1 order of magnitude:
        eps=10.d0*eps
        !increase the max number of iterations by twice the required number of successes.
        N2=N2+3*N1
        !this allowes the system to perform another 2*Nsuccess iterations and to converge and 
        !exit properly.
        !This might be a problem when performing a fixed density calculation which has to be 
        !handled differently.
     endif
     if(.not.extend_)convergence=.true.
  endif
  if(convergence.and.reset_)then
     check=0
     deallocate(Xold)
  endif
