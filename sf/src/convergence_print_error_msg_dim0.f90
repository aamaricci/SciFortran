  if(convergence)then
     write(*,"(A,F18.12)")bold_green("error="),err
  else
     if(err < eps)then
        write(*,"(A,F18.12)")bold_yellow("error="),err
     else
        write(*,"(A,F18.12)")bold_red("error="),err
     endif
  endif
