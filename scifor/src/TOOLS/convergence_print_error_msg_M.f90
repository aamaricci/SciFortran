  if(convergence)then
     write(*,"(A,F18.12)")bold_green("max error="),error(1)
     write(*,"(A,F18.12)")bold_green("    error="),err
     write(*,"(A,F18.12)")bold_green("min error="),error(2)
  else
     if(err < eps)then
        write(*,"(A,F18.12)")bold_yellow("max error="),error(1)
        write(*,"(A,F18.12)")bold_yellow("    error="),err
        write(*,"(A,F18.12)")bold_yellow("min error="),error(2)
     else
        write(*,"(A,F18.12)")bold_red("max error="),error(1)
        write(*,"(A,F18.12)")bold_red("    error="),err
        write(*,"(A,F18.12)")bold_red("min error="),error(2)
     endif
  endif
