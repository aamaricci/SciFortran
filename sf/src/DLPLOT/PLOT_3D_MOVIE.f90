!+-----------------------------------------------------------------+
!PROGRAM  : PLOT3DSURFACE_MOVIE
!TYPE     : Subroutine
!PURPOSE  : realizes an animated .gif out of the .png color plots
!+-----------------------------------------------------------------+
subroutine plot_dislin_3D_movie(pname,xname,yname,zname,X,Y,GF)
  character(len=*) :: xname,yname,zname,pname
  character(len=4) :: char_temp
  character(len=32):: xsuff,ysuff
  character(len=128):: cmd
  integer :: Nt,Nx,Ny,ixstep,iystep,i,time
  real(4) :: xstep,ystep
  real(4) :: xmin,xmax,dex
  real(4) :: ymin,ymax,dey
  real(4) :: zmin,zmax,dez
  real(8),dimension(:) :: X,Y
  real(8),dimension(:,:,:) :: GF

  Nx=size(GF,1) ;Ny=size(GF,2); Nt=size(GF,3)
  xmin=minval(X);xmax=maxval(X)
  if(xmin==xmax)then
     xmin=xmin-0.25
     xmax=xmax+0.25
  endif
  dex=abs(xmax-xmin)/5.

  ymin=minval(Y);ymax=maxval(Y)
  if(ymin==ymax)then
     ymin=ymin-0.25
     ymax=ymax+0.25
  endif
  dey=abs(ymax-ymin)/5.

  !Set the steps for the functions to plot:
  zmin=minval((GF));zmax=maxval((GF))
  GF=GF+abs(zmin) !shift to 0
  if(maxval(abs(GF))>1.d-12) GF=GF/maxval(abs(GF)) !normalization
  zmin=minval((GF));zmax=maxval((GF))
  if(zmin==zmax)zmax=zmin+0.01d0
  dez=max(abs(zmin),abs(zmax))/2.
  write(*,"(A,A)")"print movie: ",adjustl(trim(pname))
  call system("if [ -d "//trim(pname)//"_movie ]; then rm -rf "//trim(pname)//"_movie ; fi")
  call system("mkdir "//trim(pname)//"_movie ")

  do time=1,Nt
     write(char_temp,'(i4)')time
     CALL METAFL('PNG')         !DEFINES THE METAFILE FORMAT (\ie THE OUTPUT TYPE)
     CALL WINSIZ(600,400)       !16:9, SET THE DIMENSION/RESOLUTION OF THE IMAGE
     CALL SCRMOD('REVERSE')     !SET BACKGROUND TO WHITE
     CALL IMGFMT('RGB')         !SET PALETTE TO R.G.B. def:8palette
     CALL SETFIL(trim(adjustl( trim(char_temp)//trim(pname)//".png")))
     CALL FILMOD('DELETE')      !ACTION OF A PRE-EXISENT FILE def:COUNT
     !CALL UNITS('POINTS')      !CHANGE UNITS TO POINTS. def:CM
     CALL SETPAG('DA4L')        !SET PAGE SIZE: DA4L/P DIN A4, USAL/P US LETTER 
     CALL ERRDEV('APPEND')      !SET BEHAVIOR ON ERROR MESSAGE
     CALL ERRFIL('dislin.out')  !SET THE FILE NAME FOR ERROR AND OUTPUT MESSAGES
     !INIT DISLIN:
     CALL DISINI
     CALL ERRMOD('PROTOCOL','FILE')
     CALL COMPLX
     !SET SOME AXIS PROPERTIES:
     CALL TICKS(4,"XY")        !SET THE NUMBER OF TICKS per AXIS
     CALL TICLEN(45,15)         !LENGTH OF MAJ/MIN TICKS
     CALL LABDIG(3,'XY')        !SET NUMBER OF DIGITS: -2=AUTO, -1=INT, n=#DIGITS
     CALL FRAME(3)              !SET THE THICKNESS OF THE FRAME
     !SET SOME PAGE PROPERTIES: 
     CALL PAGORG('BOTTOM')     !SETS THE ORIGIN AT "BOTTOM"-LEFT
     !COLOR:
     CALL SETVLT("RAIN")        !SET COLOR PALETTE: (R)RAIN,(R)TEMP,(R)GREY,VGA
     !CHARACTER:
     CALL HEIGHT(50)            !SET CHARACHTER HEIGHT
     CALL HWFONT                !SET TO HARDWARE/SYSTEM FONT
     CALL TEXMOD('ON')          !SET LaTeX MODE
     !TITLE:
     CALL TITLIN(trim(pname)//" time:"//trim(adjustl(char_temp)),2)
     !AXIS:
     CALL NAME(trim(xname),'X')
     CALL NAME(trim(yname),'Y')
     CALL NAME(trim(zname),'Z')
     !PLOT:
     CALL INTAX
     CALL AUTRES(Nx,Ny)
     CALL GRAF3(xmin,xmax,xmin,dex,ymin,ymax,ymin,dey,0.,1.05,0.,1.)
     CALL CRVMAT(real(GF(:,:,time),4),Nx,Ny,20,20)
     !DISLIN FINALIZE:
     CALL TITLE
     CALL DISFIN
     call system("mv *"//trim(pname)//".png "//trim(pname)//"_movie/ 2>/dev/null")
  enddo
  call system("echo '#!/bin/bash' > make_gif.sh")
  call system("echo 'echo 'Converting to .gif'' >> make_gif.sh")
  call system("echo 'convert -delay 50 `ls *.png|sort -n` -loop 0 "//trim(pname)//".gif' >> make_gif.sh")
  call system("echo 'if [ ! -d "//trim(pname)//".gif ]; then rm -f *.png ; fi' >> make_gif.sh") 
  call system("chmod u+x make_gif.sh")
  call system("mv make_gif.sh "//trim(pname)//"_movie/ 2>/dev/null")
  return
end subroutine plot_dislin_3D_movie
!********************************************************************
!********************************************************************
!********************************************************************




!+-----------------------------------------------------------------+
!PROGRAM  : PLOT3DSURFACE_MOVIE
!TYPE     : Subroutine
!PURPOSE  : realizes an animated .gif out of the .png surface plots
!+-----------------------------------------------------------------+
subroutine plot_3D_surface_movie(pname,xname,yname,zname,X,Y,GF)
  character(len=*) :: xname,yname,zname,pname
  character(len=4) :: char_temp
  character(len=32):: xsuff,ysuff
  character(len=128):: cmd
  integer :: Nt,Nx,Ny,ixstep,iystep,i,time
  real(4) :: xstep,ystep
  real(4) :: xmin,xmax,dex
  real(4) :: ymin,ymax,dey
  real(4) :: zmin,zmax,dez
  real(8),dimension(:) :: X,Y
  real(8),dimension(:,:,:) :: GF

  Nx=size(GF,1) ;Ny=size(GF,2); Nt=size(GF,3)
  xmin=minval(X);xmax=maxval(X)
  if(xmin==xmax)then
     xmin=xmin-0.25
     xmax=xmax+0.25
  endif
  dex=abs(xmax-xmin)/5.

  ymin=minval(Y);ymax=maxval(Y)
  if(ymin==ymax)then
     ymin=ymin-0.25
     ymax=ymax+0.25
  endif
  dey=abs(ymax-ymin)/5.


  zmin=minval((GF));zmax=maxval((GF))
  if(zmin==zmax)zmax=zmin+0.01d0
  dez=max(abs(zmin),abs(zmax))    
  !dez=abs(zmax-zmin)
  dez=dez/5.

  write(*,"(A,A)")"print movie: ",adjustl(trim(pname))
  call system("if [ -d "//trim(pname)//"_surface_movie ]; then rm -rf "//trim(pname)//"_surface_movie ; fi")
  call system("mkdir "//trim(pname)//"_surface_movie ")
  do time=1,Nt
     write(char_temp,'(i4)')time
     CALL METAFL('PNG')         !DEFINES THE METAFILE FORMAT (\ie THE OUTPUT TYPE)
     CALL WINSIZ(800,600)       !16:9, SET THE DIMENSION/RESOLUTION OF THE IMAGE
     CALL SCRMOD('REVERSE')     !SET BACKGROUND TO WHITE
     CALL IMGFMT('RGB')         !SET PALETTE TO R.G.B. def:8palette
     CALL SETFIL(trim(adjustl( trim(char_temp)//trim(pname)//".png")))
     CALL FILMOD('DELETE')      !ACTION OF A PRE-EXISENT FILE def:COUNT

     !CALL UNITS('POINTS')      !CHANGE UNITS TO POINTS. def:CM
     !CALL SETPAG('DA4L')        !SET PAGE SIZE: DA4L/P DIN A4, USAL/P US LETTER 
     CALL ERRDEV('APPEND')      !SET BEHAVIOR ON ERROR MESSAGE
     CALL ERRFIL('dislin.out')  !SET THE FILE NAME FOR ERROR AND OUTPUT MESSAGES
     !INIT DISLIN:
     CALL DISINI
     CALL ERRMOD('PROTOCOL','FILE')
     CALL COMPLX
     !SET SOME AXIS PROPERTIES:
     CALL TICKS(4,"XY")        !SET THE NUMBER OF TICKS per AXIS
     CALL TICLEN(45,15)         !LENGTH OF MAJ/MIN TICKS
     !SET SOME PAGE PROPERTIES: 
     CALL PAGORG('BOTTOM')     !SETS THE ORIGIN AT "BOTTOM"-LEFT
     !COLOR:
     CALL SETVLT("RAIN")        !SET COLOR PALETTE: (R)RAIN,(R)TEMP,(R)GREY,VGA
     !CHARACTER:
     CALL HEIGHT(50)            !SET CHARACHTER HEIGHT
     CALL PAGERA
     CALL HWFONT                !SET TO HARDWARE/SYSTEM FONT
     CALL TEXMOD('ON')          !SET LaTeX MODE
     !TITLE:
     CALL TITLIN(trim(pname)//" time:"//trim(adjustl(char_temp)),2)
     !AXIS:
     CALL NAME(trim(xname),'X')
     CALL NAME(trim(yname),'Y')
     CALL NAME(trim(zname),'Z')
     CALL LABDIG(6,'Z')
     CALL LABDIG(0,'XY')

     !PLOT:
     CALL VIEW3D(210.0,30.,5.5,'ANGLE')
     CALL SETVLT('RAIN') !TEMP
     CALL GRAF3D(xmin,xmax,xmin,dex,ymin,ymax,ymin,dey,zmin,zmax,zmin,dez)
     CALL SHDMOD("SMOOTH","SURFACE")
     CALL SURSHD(real(X,4),Nx,real(Y,4),Ny,real(GF(:,:,time),4))
     CALL SURMAT(real(GF(:,:,time),4),Nx,Ny,1,1)

     !DISLIN FINALIZE:
     CALL TITLE
     CALL DISFIN
     call system("mv *"//trim(pname)//".png "//trim(pname)//"_surface_movie/ 2>/dev/null")    
  enddo
  call system("echo '#!/bin/bash' > make_gif.sh")
  call system("echo 'echo 'Converting to .gif'' >> make_gif.sh")
  call system("echo 'convert -delay 30 `ls *.png|sort -n` -loop 0 "//trim(pname)//".gif' >> make_gif.sh")
  call system("echo 'if [ ! -d "//trim(pname)//".gif ]; then rm -f *.png ; fi' >> make_gif.sh") 
  call system("chmod u+x make_gif.sh")
  call system("mv make_gif.sh "//trim(pname)//"_surface_movie/ 2>/dev/null")
  return
end subroutine plot_3D_surface_movie
!********************************************************************
!********************************************************************
!********************************************************************
