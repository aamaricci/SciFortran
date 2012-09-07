
!+-----------------------------------------------------------------+
!PROGRAM  : PLOT_DISLIN
!TYPE     : Subroutine
!PURPOSE  : 
!+-----------------------------------------------------------------+
subroutine plot_dislin2Dmovie_(pname,X,Y,&
     Xlabel,Ylabel, &
     wlp,           &
     Xmin,Xmax,     &
     Ymin,Ymax)

  character(len=*)              :: pname
  real(8),dimension(:)          :: X
  real(8),dimension(:,:)        :: Y
  character(len=*),optional     :: Xlabel,Ylabel
  character(len=*),optional     :: wlp    
  real(8),optional              :: Xmin,Xmax,Ymin,Ymax
  character(32)                 :: labelX,labelY
  real(4)                       :: minX,maxX,minY,maxY    
  real(4)                       :: dex,dey
  character(len=50)             :: cbuf
  character(len=3)              :: withlp
  integer                       :: narray
  integer                       :: I,K,Nx,Ny,Nt
  real(4)                       :: dt       
  character(len=7)              :: char_temp

  Nx=size(X);Nt=size(Y,1);Ny=size(Y,2)

  minX=minval(X); if(present(Xmin))minX=real(Xmin)
  maxX=maxval(X); if(present(Xmax))maxX=real(Xmax)
  if(maxX==minX)then
     minX = minX - 0.25
     maxX = maxX + 0.25
  endif

  minY=minval(Y);if(present(Ymin))minY=real(Ymin)
  maxY=maxval(Y);if(present(Ymax))maxY=real(Ymax)
  if(maxY==minY)then
     minY = minY - 0.25
     maxY = maxY + 0.25
  endif

  dex=abs(maxX-minX)/4.d0
  dey=abs(maxY-minY)/4.d0

  labelX = "X"  ; if(present(Xlabel))labelX = Xlabel
  labelY = "Y"  ; if(present(Ylabel))labelY = Ylabel
  withlp = "wl" ; if(present(wlp))   withlp = wlp

  write(*,"(a,a,a)")"print movie: ",adjustl(trim(pname))
  do i=1,Nt
     dt=dble(i)
     write(char_temp,'(f4.0)')dt
     CALL METAFL('PNG')         !DEFINES THE METAFILE FORMAT (\ie THE OUTPUT TYPE)
     CALL WINSIZ(1024,768)     !16:9, SET THE DIMENSION/RESOLUTION OF THE IMAGE
     CALL SCRMOD('REVERSE')     !SET BACKGROUND TO WHITE
     CALL IMGFMT('RGB')         !SET PALETTE TO R.G.B. def:8palette
     CALL SETFIL(trim(adjustl(trim(char_temp)//trim(pname)//".png")))!SET THE FILE NAME
     CALL FILMOD('DELETE')      !ACTION OF A PRE-EXISENT FILE def:COUNT
     CALL SETPAG('DA4L')        !SET PAGE SIZE: DA4L/P DIN A4, USAL/P US LETTER 
     CALL ERRDEV('APPEND')     !SET BEHAVIOR ON ERROR MESSAGE
     CALL ERRFIL('dislin.out') !SET THE FILE NAME FOR ERROR AND OUTPUT MESSAGES

     !INIT DISLIN:
     CALL DISINI
     CALL ERRMOD('PROTOCOL','FILE')!SET MOD ON ERROR MESSAGE
     CALL COMPLX
     CALL TICKS(5,"XY")        !SET THE NUMBER OF TICKS per AXIS
     CALL TICLEN(45,15)         !LENGTH OF MAJ/MIN TICKS
     CALL LABDIG(2,'XY')        !SET NUMBER OF DIGITS: -2=AUTO, -1=INT, n=#DIGITS
     CALL FRAME(3)              !SET THE THICKNESS OF THE FRAME
     !COLOR:
     CALL SETVLT("RAIN")        !SET COLOR PALETTE: (R)RAIN,(R)TEMP,(R)GREY,VGA
     !CHARACTER:
     CALL HEIGHT(50)            !SET CHARACHTER HEIGHT
     CALL HWFONT                !SET TO HARDWARE/SYSTEM FONT
     CALL TEXMOD('ON')          !SET LaTeX MODE
     !TITLE:
     CALL TITLIN(trim(pname)//" time:"//trim(adjustl(char_temp)),2)!SET THE TITLE AND ITS LINE
     !AXIS:
     CALL NAME(trim(labelX),'X')
     CALL NAME(trim(labelY),'Y')

     !PLOT:
     CALL GRAF(minX,maxX,minX,dex,minY,maxY,minY,dey)
     CALL GRID(2,2)

     CALL TITLE

     if(trim(withlp)=="wp")  CALL INCMRK(-1)
     if(trim(withlp)=="wl")  CALL INCMRK(0)
     if(trim(withlp)=="wlp") CALL INCMRK(1)

     CALL THKCRV(5)

     !CALL CHNCRV("COLOR")
     call PLOT2D(X,Y(I,1:Nx))
     CALL DISFIN

  enddo

  call system("if [ ! -d "//trim(pname)//"_movie ]; then mkdir "//trim(pname)//"_movie ; fi")
  call system("mv *"//trim(pname)//".png "//trim(pname)//"_movie/ 2>/dev/null")    
  call system("echo '#!/bin/bash' > make_gif.sh")
  call system("echo 'echo 'Converting to .gif'' >> make_gif.sh")
  call system("echo 'convert -delay 50 `ls *.*.png|sort -n` -loop 0 "//trim(pname)//".gif' >> make_gif.sh")
  call system("echo 'if [ ! -d "//trim(pname)//".gif ]; then rm -f *.*.png ; fi' >> make_gif.sh")
  call system("chmod u+x make_gif.sh")
  call system("mv make_gif.sh "//trim(pname)//"_movie/ 2>/dev/null")

contains
  subroutine plot2D(X,Y,color_)
    character(len=*),optional :: color_
    real(8),dimension(:) :: X,Y
    if(present(color_))CALL COLOR(color_)
    CALL CURVE(real(X),real(Y),size(X))
    CALL COLOR("FORE")
  end subroutine plot2D

end subroutine plot_dislin2Dmovie_
!********************************************************************
!********************************************************************
!********************************************************************





!+-----------------------------------------------------------------+
!PROGRAM  : PLOT_DISLIN
!TYPE     : Subroutine
!PURPOSE  : 
!+-----------------------------------------------------------------+
subroutine plot_dislin2Dmovie__(pname,X,Y,&
     Xlabel,Ylabel, &
     wlp,           &
     Xmin,Xmax,     &
     Ymin,Ymax)

  character(len=*)              :: pname
  real(8),dimension(:,:)        :: X
  real(8),dimension(:,:)        :: Y
  character(len=*),optional     :: Xlabel,Ylabel
  character(len=*),optional     :: wlp    
  real(8),optional              :: Xmin,Xmax,Ymin,Ymax
  character(32)                 :: labelX,labelY
  real(4)                       :: minX,maxX,minY,maxY    
  real(4)                       :: dex,dey
  character(len=50)             :: cbuf
  character(len=3)              :: withlp
  integer                       :: narray
  integer                       :: I,K,Nx,Ny,Nt
  real(4)                       :: dt       
  character(len=7)              :: char_temp

  Nx=size(X,2);Nt=size(Y,1);Ny=size(Y,2)
  if(Nx/=Ny)stop "Warning in plot_2Dmovie: Nx .ne. Ny"
  if(size(X,1)/=Nt)stop "Error in plot_2Dmovie: Xtime .NE. Ytime"

  minX=minval(X); if(present(Xmin))minX=real(Xmin)
  maxX=maxval(X); if(present(Xmax))maxX=real(Xmax)
  if(maxX==minX)then
     minX = minX - 0.25
     maxX = maxX + 0.25
  endif

  minY=minval(Y);if(present(Ymin))minY=real(Ymin)
  maxY=maxval(Y);if(present(Ymax))maxY=real(Ymax)
  if(maxY==minY)then
     minY = minY - 0.25
     maxY = maxY + 0.25
  endif

  dex=abs(maxX-minX)/4.d0
  dey=abs(maxY-minY)/4.d0

  labelX = "X"  ; if(present(Xlabel))labelX = Xlabel
  labelY = "Y"  ; if(present(Ylabel))labelY = Ylabel
  withlp = "wl" ; if(present(wlp))   withlp = wlp

  write(*,"(a,a,a)")"print movie: ",adjustl(trim(pname))
  do i=1,Nt
     dt=float(i)
     write(char_temp,'(f4.0)')dt
     CALL METAFL('PNG')         !DEFINES THE METAFILE FORMAT (\ie THE OUTPUT TYPE)
     CALL WINSIZ(1024,768)     !16:9, SET THE DIMENSION/RESOLUTION OF THE IMAGE
     CALL SCRMOD('REVERSE')     !SET BACKGROUND TO WHITE
     CALL IMGFMT('RGB')         !SET PALETTE TO R.G.B. def:8palette
     CALL SETFIL(trim(adjustl(trim(char_temp)//trim(pname)//".png")))!SET THE FILE NAME
     CALL FILMOD('DELETE')      !ACTION OF A PRE-EXISENT FILE def:COUNT
     CALL SETPAG('DA4L')        !SET PAGE SIZE: DA4L/P DIN A4, USAL/P US LETTER 
     CALL ERRDEV('APPEND')     !SET BEHAVIOR ON ERROR MESSAGE
     CALL ERRFIL('dislin.out') !SET THE FILE NAME FOR ERROR AND OUTPUT MESSAGES

     !INIT DISLIN:
     CALL DISINI
     CALL ERRMOD('PROTOCOL','FILE')!SET MOD ON ERROR MESSAGE
     CALL COMPLX
     CALL TICKS(5,"XY")        !SET THE NUMBER OF TICKS per AXIS
     CALL TICLEN(45,15)         !LENGTH OF MAJ/MIN TICKS
     CALL LABDIG(2,'XY')        !SET NUMBER OF DIGITS: -2=AUTO, -1=INT, n=#DIGITS
     CALL FRAME(3)              !SET THE THICKNESS OF THE FRAME
     !COLOR:
     CALL SETVLT("RAIN")        !SET COLOR PALETTE: (R)RAIN,(R)TEMP,(R)GREY,VGA
     !CHARACTER:
     CALL HEIGHT(50)            !SET CHARACHTER HEIGHT
     CALL HWFONT                !SET TO HARDWARE/SYSTEM FONT
     CALL TEXMOD('ON')          !SET LaTeX MODE
     !TITLE:
     CALL TITLIN(trim(pname)//" time:"//trim(adjustl(char_temp)),2)!SET THE TITLE AND ITS LINE
     !AXIS:
     CALL NAME(trim(labelX),'X')
     CALL NAME(trim(labelY),'Y')

     !PLOT:
     CALL GRAF(minX,maxX,minX,dex,minY,maxY,minY,dey)
     CALL GRID(2,2)

     CALL TITLE

     if(trim(withlp)=="wp")  CALL INCMRK(-1)
     if(trim(withlp)=="wl")  CALL INCMRK(0)
     if(trim(withlp)=="wlp") CALL INCMRK(1)

     CALL THKCRV(5)

     !CALL CHNCRV("COLOR")
     call PLOT2D(X(I,:),Y(I,1:Nx))
     CALL DISFIN

  enddo


  call system("if [ ! -d "//trim(pname)//"_movie ]; then mkdir "//trim(pname)//"_movie ; fi")
  call system("mv *"//trim(pname)//".png "//trim(pname)//"_movie/ 2>/dev/null")    
  call system("echo '#!/bin/bash' > make_gif.sh")
  call system("echo 'echo 'Converting to .gif'' >> make_gif.sh")
  call system("echo 'convert -delay 50 `ls *.*.png|sort -n` -loop 0 "//trim(pname)//".gif'  >> make_gif.sh")
  call system("echo 'if [ ! -d "//trim(pname)//".gif ]; then rm -f *.*.png ; fi' >> make_gif.sh")
  call system("chmod u+x make_gif.sh")
  call system("mv make_gif.sh "//trim(pname)//"_movie/ 2>/dev/null")

  !=========================================
contains
  subroutine plot2D(X,Y,color_)
    character(len=*),optional :: color_
    real(8),dimension(:) :: X,Y
    if(present(color_))CALL COLOR(color_)
    CALL CURVE(real(X),real(Y),size(X))
    CALL COLOR("FORE")
  end subroutine plot2D
end subroutine plot_dislin2Dmovie__
!********************************************************************
!********************************************************************
!********************************************************************
