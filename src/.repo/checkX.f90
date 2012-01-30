  if(present(X))then
     Nx=size(X)
     if(Nx < Ny1) then
        write(*,*)"big problem while printing "//trim(pname)//": skipping"
        return
     endif
     if(Nx/=Ny1 .OR. Nx/=Ny2) write(*,"(a,1x,I6,I6,I6)")"problem while printing "//trim(pname)//" Nx,Ny1,Ny2",Nx,Ny1,Ny2
  endif
