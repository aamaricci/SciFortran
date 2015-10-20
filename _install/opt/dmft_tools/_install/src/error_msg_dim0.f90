  if(convergence)then
     write(*,"(A,ES15.7)")bold_green("error="),err
  else
     if(err < eps)then
        write(*,"(A,ES15.7)")bold_yellow("error="),err
     else
        write(*,"(A,ES15.7)")bold_red("error="),err
     endif
  endif
