!+-----------------------------------------------------------------+
!PROGRAM  : 
!TYPE     : Subroutine
!PURPOSE  : 
!+-----------------------------------------------------------------+
subroutine plot_dislinVF_(pname,X,Y,Vx,Vy,Xlabel,Ylabel)
  character(len=*)                   :: pname
  real(8),dimension(:)               :: X,Y
  real(8),dimension(size(X),size(Y)) :: Vx,Vy
  real(4),dimension(size(X),size(Y)) :: WMAT
  integer,dimension(size(X),size(Y)) :: ITMAT,IWMAT
  character(len=*),optional     :: Xlabel,Ylabel
  character(len=50)             :: labelX,labelY
  real(8)                       :: Xmin,Xmax,dex
  real(8)                       :: Ymin,Ymax,dey
  real(8)                       :: ZXmin,ZXmax,ZYmin,ZYmax,Zmin,Zmax,dez
  integer                       :: Nx,Ny
  integer                       :: i,j
  Nx=size(X);Ny=size(Y)

  Xmin=minval(X); Xmax=maxval(X)
  if(Xmin==Xmax)then
     Xmin = Xmin - 0.25
     Xmax = Xmax + 0.25
  endif

  Ymin=minval(Y); Ymax=maxval(Y)
  if(Ymax==Ymin)then
     Ymin = Ymin - 0.25
     Ymax = Ymax + 0.25
  endif

  dex=abs(Xmax-Xmin)/4.d0
  dey=abs(Ymax-Ymin)/4.d0

  labelX = "X"  ; if(present(Xlabel))labelX = Xlabel
  labelY = "Y"  ; if(present(Ylabel))labelY = Ylabel

  write(*,"(A,A)")"print: ",adjustl(trim(pname))

  CALL METAFL('PNG')         !DEFINES THE METAFILE FORMAT (\ie THE OUTPUT TYPE)
  CALL WINSIZ(1920,1280)     !16:9, SET THE DIMENSION/RESOLUTION OF THE IMAGE
  CALL SCRMOD('REVERSE')     !SET BACKGROUND TO WHITE
  CALL IMGFMT('RGB')         !SET PALETTE TO R.G.B. def:8palette
  CALL SETFIL(adjustl(trim(pname))//".png")   !SET THE FILE NAME
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
  CALL LABDIG(-2,'XY')        !SET NUMBER OF DIGITS: -2=AUTO, -1=INT, n=#DIGITS
  CALL FRAME(3)              !SET THE THICKNESS OF THE FRAME
  !COLOR:
  CALL SETVLT("RAIN")        !SET COLOR PALETTE: (R)RAIN,(R)TEMP,(R)GREY,VGA
  !CHARACTER:
  CALL HEIGHT(50)            !SET CHARACHTER HEIGHT
  CALL HWFONT                !SET TO HARDWARE/SYSTEM FONT
  CALL TEXMOD('ON')          !SET LaTeX MODE
  !TITLE:
  CALL TITLIN(adjustl(trim(pname)),1) !SET THE TITLE AND ITS LINE
  !AXIS:
  CALL NAME(trim(labelX),'X')
  CALL NAME(trim(labelY),'Y')

  !Vector plot OLD type: small arrows at each point:
  !===========================================================
  CALL GRAF(real(Xmin,4),real(Xmax,4),real(Xmin,4),real(dex,4),&
       real(Ymin,4),real(Ymax,4),real(Ymin,4),real(dey,4))
  CALL VECCLR(-1)
  !#CALL VECOPT(0.,"SCALE")
  CALL VECMAT(real(Vx,4),real(Vy,4),Nx,Ny,real(X,4),real(Y,4),-1)
  !#CALL STREAM(real(V%x,4),real(V%y,4),N,N,X,Y,NULLX,NULLY,0)


  !Vector plot NEW type: vector intensity plot:
  !===========================================================
  !CALL TXTURE(ITMAT, Nx, Ny) !Generatet a random texture:
  !We could do something better and faster using MKL random library
  !This is not implemented yet to keep this module DISLIN dependent only
  ! do i=1,Nx
  !    do j=1,Ny
  !       ITMAT(i,j)=int(drand()*10)
  !    enddo
  ! enddo
  ! CALL LICMOD('ON','SCALE')
  ! CALL LICPTS(real(Vx,4),real(Vy,4),Nx,Ny,ITMAT,IWMAT,WMAT)
  ! CALL AUTRES(Nx,Ny)
  ! CALL GRAF3(real(Xmin,4),real(Xmax,4),real(Xmin,4),real(dex,4),&
  !      real(Ymin,4),real(Ymax,4),real(Ymin,4),real(dey,4),&
  !      0.,10.,0.,2.5)
  ! CALL CRVMAT(WMAT,NX,NY,20,20)

  CALL TITLE
  CALL DISFIN

  call system("if [ ! -d PNG ]; then mkdir PNG; fi")
  call system("mv *.png PNG/")
end subroutine plot_dislinVF_
!-----------------------------------
!-----------------------------------
!-----------------------------------
subroutine plot_dislinVF__(pname,X,Y,Vx,Vy,Xlabel,Ylabel)
  character(len=*)                     :: pname
  real(8),dimension(:)                 :: X,Y
  real(8),dimension(:,:,:)             :: Vx,Vy
  real(4),dimension(size(X),size(Y))   :: WMAT
  integer,dimension(size(X),size(Y))   :: ITMAT,IWMAT
  character(len=*),optional     :: Xlabel,Ylabel
  character(len=50)             :: labelX,labelY
  real(8)                       :: Xmin,Xmax,dex
  real(8)                       :: Ymin,Ymax,dey
  integer                       :: I,Nx,Ny,Nt,ki,kj
  real(4)          :: dt       
  character(len=7) :: char_temp

  Nx=size(X);Ny=size(Y)
  Nt=size(Vx,1)
  if(size(Vy,1) < Nt)Nt=size(Vy,1)

  Xmin=minval(X); Xmax=maxval(X)
  if(Xmin==Xmax)then
     Xmin = Xmin - 0.25
     Xmax = Xmax + 0.25
  endif

  Ymin=minval(Y); Ymax=maxval(Y)
  if(Ymax==Ymin)then
     Ymin = Ymin - 0.25
     Ymax = Ymax + 0.25
  endif

  dex=abs(Xmax-Xmin)/4.d0
  dey=abs(Ymax-Ymin)/4.d0

  labelX = "X"  ; if(present(Xlabel))labelX = Xlabel
  labelY = "Y"  ; if(present(Ylabel))labelY = Ylabel

  write(*,"(A,A)")"print movie: ",adjustl(trim(pname))
  do I=1,Nt
     dt=float(I)
     write(char_temp,'(f4.0)')dt       
     !print*,char_temp
     CALL METAFL('PNG')         !DEFINES THE METAFILE FORMAT (\ie THE OUTPUT TYPE)
     CALL WINSIZ(800,600)       !16:9, SET THE DIMENSION/RESOLUTION OF THE IMAGE
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
     CALL LABDIG(-2,'XY')        !SET NUMBER OF DIGITS: -2=AUTO, -1=INT, n=#DIGITS
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

     !Vector plot OLD type: small arrows at each point:
     !===========================================================
     CALL GRAF(real(Xmin,4),real(Xmax,4),real(Xmin,4),real(dex,4),&
          real(Ymin,4),real(Ymax,4),real(Ymin,4),real(dey,4))
     CALL VECCLR(-1)
     !#CALL VECOPT(0.,"SCALE")
     CALL VECMAT(real(Vx(I,:,:),4),real(Vy(I,:,:),4),&
          Nx,Ny,real(X,4),real(Y,4),-1)
     !#CALL STREAM(real(V%x,4),real(V%y,4),N,N,X,Y,NULLX,NULLY,0)

     !Vector plot NEW type: vector intensity plot:
     !===========================================================
     !CALL TXTURE(ITMAT, Nx, Ny) !Generatet a random texture:
     ! ITMAT=0
     ! do ki=1,Nx
     !    do kj=1,Ny
     !       ITMAT(ki,kj)=int(drand()*10)
     !    enddo
     ! enddo
     ! CALL LICMOD('ON','SCALE')
     ! CALL LICPTS(real(Vx(I,:,:),4),real(Vy(I,:,:),4),Nx,Ny,ITMAT,IWMAT,WMAT)
     ! CALL AUTRES(Nx,Ny)
     ! CALL GRAF3(real(Xmin,4),real(Xmax,4),real(Xmin,4),real(dex,4),&
     !      real(Ymin,4),real(Ymax,4),real(Ymin,4),real(dey,4),&
     !      0.,10.,0.,2.5)
     ! CALL CRVMAT(WMAT,NX,NY,20,20)

     CALL TITLE
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
  return
end subroutine plot_dislinVF__
!********************************************************************
!********************************************************************
!********************************************************************
