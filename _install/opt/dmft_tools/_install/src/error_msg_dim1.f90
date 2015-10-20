  if(convergence)then
     write(*,"(A,ES15.7)")bold_green("max error="),error(1)
     write(*,"(A,ES15.7)")bold_green("    error="),err
     write(*,"(A,ES15.7)")bold_green("min error="),error(2)
  else
     if(err < eps)then
        write(*,"(A,ES15.7)")bold_yellow("max error="),error(1)
        write(*,"(A,ES15.7)")bold_yellow("    error="),err
        write(*,"(A,ES15.7)")bold_yellow("min error="),error(2)
     else
        write(*,"(A,ES15.7)")bold_red("max error="),error(1)
        write(*,"(A,ES15.7)")bold_red("    error="),err
        write(*,"(A,ES15.7)")bold_red("min error="),error(2)
     endif
  endif
