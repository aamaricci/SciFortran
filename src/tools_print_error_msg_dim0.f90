  if(convergence)then
     write(*,"(A,ES18.9)")bold_green("error="),err
  else
     if(err < eps)then
        write(*,"(A,ES18.9)")bold_yellow("error="),err
     else
        write(*,"(A,ES18.9)")bold_red("error="),err
     endif
  endif
