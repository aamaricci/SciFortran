!+-----------------------------------------------------------------+
!PROGRAM  : PLOT
!TYPE     : Subroutine
!PURPOSE  : 
!+-----------------------------------------------------------------+
subroutine D_plot(pname,X,Y1,&
     Y1label,       &
     Y2,Y2label,    &
     Y3,Y3label,    &
     Y4,Y4label,    &
     Y5,Y5label,    &
     Y6,Y6label,    &
     Xlabel,Ylabel, &
     wlp,           &
     Xmin,Xmax,     &
     Ymin,Ymax)
  character(len=*)              :: pname
  real(8),dimension(:)          :: X
  real(8),dimension(:)          :: Y1
  real(8),dimension(:),optional :: Y2,Y3,Y4,Y5,Y6
  character(len=*),optional     :: Y1label,Y2label,Y3label,Y4label,Y5label,Y6label
  character(len=*),optional     :: Xlabel,Ylabel
  character(len=*),optional     :: wlp
  real(8),optional              :: Xmin,Xmax,Ymin,Ymax
  character(len=50)             :: labelX,labelY
  character(len=3)              :: withlp
  real(8)                       :: minX,maxX,dex
  real(8)                       :: minY,maxY,dey
  character(len=50)             :: cbuf
  integer                       :: narray
  character(len=64)             :: labelY1,labelY2,labelY3,labelY4,labelY5,labelY6
  !Get all the variables to pass, with names changed w/respect to input vars of plot_dislin,plot_xmgrace
  minX=minval(X); if(present(Xmin))minX=real(Xmin)
  maxX=maxval(X); if(present(Xmax))maxX=real(Xmax)
  if(maxX==minX)then
     minX = minX - 0.25
     maxX = maxX + 0.25
  endif

  minY=minval(Y1)
  if(present(Y2))minY=min(minval(Y1),minval(Y2))
  if(present(Y3))minY=min(minval(Y1),minval(Y2),minval(Y3))
  if(present(Y4))minY=min(minval(Y1),minval(Y2),minval(Y3),minval(Y4))
  if(present(Y5))minY=min(minval(Y1),minval(Y2),minval(Y3),minval(Y4),minval(Y5))
  if(present(Ymin))minY=real(Ymin)
  maxY=maxval(Y1)
  if(present(Y2))maxY=max(maxval(Y1),maxval(Y2))
  if(present(Y3))maxY=max(maxval(Y1),maxval(Y2),maxval(Y3))
  if(present(Y4))maxY=max(maxval(Y1),maxval(Y2),maxval(Y3),maxval(Y4))
  if(present(Y5))maxY=max(maxval(Y1),maxval(Y2),maxval(Y3),maxval(Y4),maxval(Y5))
  if(present(Ymax))maxY=real(Ymax)
  if(maxY==minY)then
     minY = minY - 0.25
     maxY = maxY + 0.25
  endif

  labelX = "X"  ; if(present(Xlabel))labelX = Xlabel
  labelY = "Y"  ; if(present(Ylabel))labelY = Ylabel
  withlp = "wl" ; if(present(wlp))   withlp = wlp

  labelY1 = "Y1" ; if(present(Y1label))labelY1 = Y1label
  labelY2 = "Y2" ; if(present(Y2label))labelY2 = Y2label
  labelY3 = "Y3" ; if(present(Y3label))labelY3 = Y3label
  labelY4 = "Y4" ; if(present(Y4label))labelY4 = Y4label
  labelY5 = "Y5" ; if(present(Y5label))labelY5 = Y5label
  labelY6 = "Y6" ; if(present(Y6label))labelY6 = Y6label

  narray=1
  if(present(Y2))  narray=narray+1
  if(present(Y3))  narray=narray+1
  if(present(Y4))  narray=narray+1
  if(present(Y5))  narray=narray+1
  if(present(Y6))  narray=narray+1

  write(*,"(A,A)")"print: ",adjustl(trim(pname))
  if(narray==1)then
     call plot_dislin(pname, X,Y1,labelY1, &
          Xlabel=labelX, Ylabel=labelY, wlp=withlp, Xmin=minX, Xmax=maxX, Ymin=minY, Ymax=maxY)
     call plot_xmgrace(pname,X,Y1,labelY1, &
          Xlabel=labelX, Ylabel=labelY)

  elseif(narray==2)then
     call plot_dislin(pname, X,Y1,labelY1,Y2,labelY2, &
          Xlabel=labelX, Ylabel=labelY, wlp=withlp, Xmin=minX, Xmax=maxX, Ymin=minY, Ymax=maxY)
     call plot_xmgrace(pname,X,Y1,labelY1,Y2,labelY2, &
          Xlabel=labelX, Ylabel=labelY)

  elseif(narray==3)then
     call plot_dislin(pname, X,Y1,labelY1,Y2,labelY2,Y3,labelY3, &
          Xlabel=labelX, Ylabel=labelY, wlp=withlp, Xmin=minX, Xmax=maxX, Ymin=minY, Ymax=maxY)
     call plot_xmgrace(pname,X,Y1,labelY1,Y2,labelY2,Y3,labelY3, &
          Xlabel=labelX, Ylabel=labelY)

  elseif(narray==4)then
     call plot_dislin(pname, X,Y1,labelY1,Y2,labelY2,Y3,labelY3,Y4,labelY4, &
          Xlabel=labelX, Ylabel=labelY, wlp=withlp, Xmin=minX, Xmax=maxX, Ymin=minY, Ymax=maxY)
     call plot_xmgrace(pname,X,Y1,labelY1,Y2,labelY2,Y3,labelY3,Y4,labelY4, &
          Xlabel=labelX, Ylabel=labelY)

  elseif(narray==5)then
     call plot_dislin(pname, X,Y1,labelY1,Y2,labelY2,Y3,labelY3,Y4,labelY4,Y5,labelY5, &
          Xlabel=labelX,Ylabel=labelY, wlp=withlp, Xmin=minX, Xmax=maxX, Ymin=minY, Ymax=maxY)
     call plot_xmgrace(pname,X,Y1,labelY1,Y2,labelY2,Y3,labelY3,Y4,labelY4,Y5,labelY5, &
          Xlabel=labelX, Ylabel=labelY)

  elseif(narray==6)then
     call plot_dislin(pname, X,Y1,labelY1,Y2,labelY2,Y3,labelY3,Y4,labelY4,Y5,labelY5,Y6,labelY6, &
          Xlabel=labelX, Ylabel=labelY, wlp=withlp, Xmin=minX, Xmax=maxX, Ymin=minY, Ymax=maxY)
     call plot_xmgrace(pname,X,Y1,labelY1,Y2,labelY2,Y3,labelY3,Y4,labelY4,Y5,labelY5,Y6,labelY6, &
          Xlabel=labelX, Ylabel=labelY)
  end if

  call system("if [ ! -d AGR ]; then mkdir AGR; fi")
  call system("mv *.agr AGR/")
  call system("if [ ! -d PNG ]; then mkdir PNG; fi")
  call system("mv *.png PNG/")
end subroutine D_plot

!-----------------------------------
!-----------------------------------
!-----------------------------------
subroutine Z_plot(pname,X,Y1,Y1label,       &
     Y2,Y2label,    &
     Y3,Y3label,    &
     Xlabel,Ylabel, &
     wlp,           &
     Xmin,Xmax,     &
     Ymin,Ymax)
  character(len=*)              :: pname
  real(8),dimension(:)          :: X
  complex(8),dimension(:)          :: Y1
  complex(8),dimension(:),optional :: Y2,Y3
  character(len=*),optional     :: Y1label,Y2label,Y3label
  character(len=*),optional     :: Xlabel,Ylabel
  character(len=*),optional     :: wlp
  real(8),optional              :: Xmin,Xmax,Ymin,Ymax
  character(len=50)             :: labelX,labelY
  character(len=3)              :: withlp
  real(8)                       :: minX,maxX,dex
  real(8)                       :: minY,R_minY,R_maxY,R_dey
  real(8)                       :: maxY,I_minY,I_maxY,I_dey
  character(len=50)             :: cbuf
  integer                       :: narray
  character(len=64)             :: labelY1,labelY2,labelY3
  !Get all the variables to pass, with names changed w/respect to input vars of plot_dislin,plot_xmgrace
  minX=minval(X); if(present(Xmin))minX=real(Xmin)
  maxX=maxval(X); if(present(Xmax))maxX=real(Xmax)
  if(maxX==minX)then
     minX = minX - 0.25
     maxX = maxX + 0.25
  endif

  R_minY=minval(real(Y1))
  if(present(Y2))R_minY=min(minval(real(Y1)),minval(real(Y2)))
  if(present(Y3))R_minY=min(minval(real(Y1)),minval(real(Y2)),minval(real(Y3)))
  if(present(Ymin))R_minY=real(Ymin)
  R_maxY=maxval(real(Y1))
  if(present(Y2))R_maxY=max(maxval(real(Y1)),maxval(real(Y2)))
  if(present(Y3))R_maxY=max(maxval(real(Y1)),maxval(real(Y2)),maxval(real(Y3)))
  if(present(Ymax))R_maxY=real(Ymax)
  if(R_maxY==R_minY)then
     R_minY = R_minY - 0.25
     R_maxY = R_maxY + 0.25
  endif

  I_minY=minval(aimag(Y1))
  if(present(Y2))I_minY=min(minval(aimag(Y1)),minval(aimag(Y2)))
  if(present(Y3))I_minY=min(minval(aimag(Y1)),minval(aimag(Y2)),minval(aimag(Y3)))
  if(present(Ymin))I_minY=real(Ymin)
  I_maxY=maxval(aimag(Y1))
  if(present(Y2))I_maxY=max(maxval(aimag(Y1)),maxval(aimag(Y2)))
  if(present(Y3))I_maxY=max(maxval(aimag(Y1)),maxval(aimag(Y2)),maxval(aimag(Y3)))
  if(present(Ymax))I_maxY=real(Ymax)
  if(I_maxY==I_minY)then
     I_minY = I_minY - 0.25
     I_maxY = I_maxY + 0.25
  endif

  labelX = "X"  ; if(present(Xlabel))labelX = Xlabel
  labelY = "Y"  ; if(present(Ylabel))labelY = Ylabel
  withlp = "wl" ; if(present(wlp))   withlp = wlp

  labelY1 = "Y1" ; if(present(Y1label))labelY1 = Y1label
  labelY2 = "Y2" ; if(present(Y2label))labelY2 = Y2label
  labelY3 = "Y3" ; if(present(Y3label))labelY3 = Y3label

  narray=1
  if(present(Y2))  narray=narray+1
  if(present(Y3))  narray=narray+1

  write(*,"(A,A)")"print: ",adjustl(trim(pname))
  if(narray==1)then
     call plot_dislin("Re_"//pname, X,real(Y1),"Re_"//labelY1, &
          Xlabel=labelX, Ylabel=labelY, wlp=withlp, Xmin=minX, Xmax=maxX, Ymin=R_minY, Ymax=R_maxY)
     call plot_dislin("Im_"//pname, X,aimag(Y1),"Im_"//labelY1, &
          Xlabel=labelX, Ylabel=labelY, wlp=withlp, Xmin=minX, Xmax=maxX, Ymin=I_minY, Ymax=I_maxY)
     call plot_xmgrace("Re_"//pname,X,real(Y1),"Re_"//labelY1, &
          Xlabel=labelX, Ylabel=labelY)
     call plot_xmgrace("Im_"//pname,X,aimag(Y1),"Im_"//labelY1, &
          Xlabel=labelX, Ylabel=labelY)

  elseif(narray==2)then
     call plot_dislin("Re_"//pname, X,real(Y1),"Re_"//labelY1,real(Y2),"Re_"//labelY2, &
          Xlabel=labelX, Ylabel=labelY, wlp=withlp, Xmin=minX, Xmax=maxX, Ymin=R_minY, Ymax=R_maxY)
     call plot_dislin("Im_"//pname, X,aimag(Y1),"Im_"//labelY1,aimag(Y2),"Im_"//labelY2, &
          Xlabel=labelX, Ylabel=labelY, wlp=withlp, Xmin=minX, Xmax=maxX, Ymin=I_minY, Ymax=I_maxY)
     call plot_xmgrace("Re_"//pname,X,real(Y1),"Re_"//labelY1,real(Y2),"Re_"//labelY2, &
          Xlabel=labelX, Ylabel=labelY)
     call plot_xmgrace("Im_"//pname,X,aimag(Y1),"Im_"//labelY1,aimag(Y2),"Im_"//labelY2, &
          Xlabel=labelX, Ylabel=labelY)

  elseif(narray==3)then
     call plot_dislin("Re_"//pname, X,real(Y1),"Re_"//labelY1,real(Y2),"Re_"//labelY2,real(Y3),"Re_"//labelY3, &
          Xlabel=labelX, Ylabel=labelY, wlp=withlp, Xmin=minX, Xmax=maxX, Ymin=R_minY, Ymax=R_maxY)
     call plot_dislin("Im_"//pname, X,aimag(Y1),"Im_"//labelY1,aimag(Y2),"Im_"//labelY2,aimag(Y3),"Re_"//labelY3, &
          Xlabel=labelX, Ylabel=labelY, wlp=withlp, Xmin=minX, Xmax=maxX, Ymin=I_minY, Ymax=I_maxY)
     call plot_xmgrace("Re_"//pname,X,real(Y1),"Re_"//labelY1,real(Y2),"Re_"//labelY2,real(Y3),"Re_"//labelY3, &
          Xlabel=labelX, Ylabel=labelY)
     call plot_xmgrace("Im_"//pname,X,aimag(Y1),"Im_"//labelY1,aimag(Y2),"Re_"//labelY2,aimag(Y3),"Re_"//labelY3, &
          Xlabel=labelX, Ylabel=labelY)

  end if

  call system("if [ ! -d AGR ]; then mkdir AGR; fi")
  call system("mv *.agr AGR/")
  call system("if [ ! -d PNG ]; then mkdir PNG; fi")
  call system("mv *.png PNG/")
end subroutine Z_plot
!------------------------
!------------------------
!------------------------

subroutine X_plot(pname,X,Y1,Y1label,       &
     Y2,Y2label,    &
     Y3,Y3label,    &
     Y4,Y4label,    &
     Y5,Y5label,    &
     Y6,Y6label,    &
     Xlabel,Ylabel, &
     wlp,           &
     Xmin,Xmax,     &
     Ymin,Ymax)
  character(len=*)              :: pname
  real(8)             :: X
  real(8)             :: Y1
  real(8),optional    :: Y2,Y3,Y4,Y5,Y6
  real(8),dimension(1) :: xa,ya,yb,yc,yd,ye,yf
  character(len=*),optional     :: Y1label,Y2label,Y3label,Y4label,Y5label,Y6label
  character(len=*),optional     :: Xlabel,Ylabel
  character(len=*),optional     :: wlp
  real(8),optional              :: Xmin,Xmax,Ymin,Ymax
  character(len=50)             :: labelX,labelY
  character(len=3)              :: withlp
  real(8)                       :: minX,maxX,dex
  real(8)                       :: minY,maxY,dey
  character(len=50)             :: cbuf
  integer                       :: narray
  character(len=64)             :: labelY1,labelY2,labelY3,labelY4,labelY5,labelY6
  !Get all the variables to pass, with names changed w/respect to input vars of plot_dislin,plot_xmgrace
  minX=X-0.25; if(present(Xmin))minX=real(Xmin)
  maxX=X+0.25; if(present(Xmax))maxX=real(Xmax)

  minY=Y1
  if(present(Y2))minY=min(Y1,Y2)
  if(present(Y3))minY=min(Y1,Y2,Y3)
  if(present(Y4))minY=min(Y1,Y2,Y3,Y4)
  if(present(Y5))minY=min(Y1,Y2,Y3,Y4,Y5)
  if(present(Y6))minY=min(Y1,Y2,Y3,Y4,Y5,Y6)
  if(present(Ymin))minY=real(Ymin)
  maxY=Y1
  if(present(Y2))maxY=min(Y1,Y2)
  if(present(Y3))maxY=min(Y1,Y2,Y3)
  if(present(Y4))maxY=min(Y1,Y2,Y3,Y4)
  if(present(Y5))maxY=min(Y1,Y2,Y3,Y4,Y5)
  if(present(Y6))maxY=min(Y1,Y2,Y3,Y4,Y5,Y6)
  if(present(Ymax))maxY=real(Ymax)
  minY=minY - 0.25
  maxY=maxY + 0.25

  labelX = "X"  ; if(present(Xlabel))labelX = Xlabel
  labelY = "Y"  ; if(present(Ylabel))labelY = Ylabel
  withlp = "wl" ; if(present(wlp))   withlp = wlp

  labelY1 = "Y1" ; if(present(Y1label))labelY1 = Y1label
  labelY2 = "Y2" ; if(present(Y2label))labelY2 = Y2label
  labelY3 = "Y3" ; if(present(Y3label))labelY3 = Y3label
  labelY4 = "Y4" ; if(present(Y4label))labelY4 = Y4label
  labelY5 = "Y5" ; if(present(Y5label))labelY5 = Y5label
  labelY6 = "Y6" ; if(present(Y6label))labelY6 = Y6label

  narray=1
  if(present(Y2))  narray=narray+1
  if(present(Y3))  narray=narray+1
  if(present(Y4))  narray=narray+1
  if(present(Y5))  narray=narray+1
  if(present(Y6))  narray=narray+1

  xa(1) = X
  ya(1) = Y1
  if(present(Y2))  yb(1)=Y2
  if(present(Y3))  yc(1)=Y3
  if(present(Y4))  yd(1)=Y4
  if(present(Y5))  ye(1)=Y5
  if(present(Y6))  yf(1)=Y6

  write(*,"(A,A)")"print: ",adjustl(trim(pname))
  if(narray==1)then
     call plot_dislin(pname, xa,ya,labelY1, &
          Xlabel=labelX, Ylabel=labelY, wlp=withlp, Xmin=minX, Xmax=maxX, Ymin=minY, Ymax=maxY)
     call plot_xmgrace(pname,xa,ya,labelY1, &
          Xlabel=labelX, Ylabel=labelY)

  elseif(narray==2)then
     call plot_dislin(pname, Xa,Ya,labelY1,Yb,labelY2, &
          Xlabel=labelX, Ylabel=labelY, wlp=withlp, Xmin=minX, Xmax=maxX, Ymin=minY, Ymax=maxY)
     call plot_xmgrace(pname,Xa,Ya,labelY1,Yb,labelY2, &
          Xlabel=labelX, Ylabel=labelY)

  elseif(narray==3)then
     call plot_dislin(pname, Xa,Ya,labelY1,Yb,labelY2,Yc,labelY3, &
          Xlabel=labelX, Ylabel=labelY, wlp=withlp, Xmin=minX, Xmax=maxX, Ymin=minY, Ymax=maxY)
     call plot_xmgrace(pname,Xa,Ya,labelY1,Yb,labelY2,Yc,labelY3, &
          Xlabel=labelX, Ylabel=labelY)

  elseif(narray==4)then
     call plot_dislin(pname, Xa,Ya,labelY1,Yb,labelY2,Yc,labelY3,Yd,labelY4, &
          Xlabel=labelX, Ylabel=labelY, wlp=withlp, Xmin=minX, Xmax=maxX, Ymin=minY, Ymax=maxY)
     call plot_xmgrace(pname,Xa,Ya,labelY1,Yb,labelY2,Yc,labelY3,Yd,labelY4, &
          Xlabel=labelX, Ylabel=labelY)

  elseif(narray==5)then
     call plot_dislin(pname, Xa,Ya,labelY1,Yb,labelY2,Yc,labelY3,Yd,labelY4,Ye,labelY5, &
          Xlabel=labelX,Ylabel=labelY, wlp=withlp, Xmin=minX, Xmax=maxX, Ymin=minY, Ymax=maxY)
     call plot_xmgrace(pname,Xa,Ya,labelY1,Yb,labelY2,Yc,labelY3,Yd,labelY4,Ye,labelY5, &
          Xlabel=labelX, Ylabel=labelY)

  elseif(narray==6)then
     call plot_dislin(pname, Xa,Ya,labelY1,Yb,labelY2,Yc,labelY3,Yd,labelY4,Ye,labelY5,Yf,labelY6, &
          Xlabel=labelX, Ylabel=labelY, wlp=withlp, Xmin=minX, Xmax=maxX, Ymin=minY, Ymax=maxY)
     call plot_xmgrace(pname,Xa,Ya,labelY1,Yb,labelY2,Yc,labelY3,Yd,labelY4,Ye,labelY5,Yf,labelY6, &
          Xlabel=labelX, Ylabel=labelY)
  end if

  call system("if [ ! -d AGR ]; then mkdir AGR; fi")
  call system("mv *.agr AGR/")
  call system("if [ ! -d PNG ]; then mkdir PNG; fi")
  call system("mv *.png PNG/")
end subroutine X_plot
!********************************************************************
!********************************************************************
!********************************************************************



!+-----------------------------------------------------------------+
!PROGRAM  : PLOT_DISLIN
!TYPE     : Subroutine
!PURPOSE  : 
!+-----------------------------------------------------------------+
subroutine plot_dislin(pname,X,Y1,Y1label,&
     Y2,Y2label,    &
     Y3,Y3label,    &
     Y4,Y4label,    &
     Y5,Y5label,    &
     Y6,Y6label,    &
     Xlabel,Ylabel, &
     wlp,           &
     Xmin,Xmax,     &
     Ymin,Ymax)

  character(len=*)              :: pname
  real(8),dimension(:)          :: X
  real(8),dimension(:)          :: Y1
  character(len=*)              :: Y1label
  real(8),dimension(:),optional :: Y2,Y3,Y4,Y5,Y6
  character(len=*),optional     :: Y2label,Y3label,Y4label,Y5label,Y6label
  character(len=*)              :: Xlabel,Ylabel
  character(len=*)              :: wlp    
  real(8)                       :: Xmin,Xmax,dex
  real(8)                       :: Ymin,Ymax,dey
  character(len=50)             :: cbuf
  integer                       :: narray

  dex=abs(Xmax-Xmin)/4.d0
  dey=abs(Ymax-Ymin)/4.d0


  narray=1
  if(present(Y2))  narray=narray+1
  if(present(Y3))  narray=narray+1
  if(present(Y4))  narray=narray+1
  if(present(Y5))  narray=narray+1
  if(present(Y6))  narray=narray+1

  CALL METAFL('PNG')         !DEFINES THE METAFILE FORMAT (\ie THE OUTPUT TYPE)
  CALL WINSIZ(1920,1280)     !16:9, SET THE DIMENSION/RESOLUTION OF THE IMAGE
  CALL SCRMOD('REVERSE')     !SET BACKGROUND TO WHITE
  CALL IMGFMT('RGB')         !SET PALETTE TO R.G.B. def:8palette
  CALL SETFIL(adjustl(trim(pname))//".png")   !SET THE FILE NAME
  CALL FILMOD('DELETE')      !ACTION OF A PRE-EXISENT FILE def:COUNT
  !CALL UNITS('POINTS')       !CHANGE UNITS TO POINTS. def:CM
  CALL SETPAG('DA4L')        !SET PAGE SIZE: DA4L/P DIN A4, USAL/P US LETTER 
  CALL ERRDEV('APPEND')     !SET BEHAVIOR ON ERROR MESSAGE
  CALL ERRFIL('dislin.out') !SET THE FILE NAME FOR ERROR AND OUTPUT MESSAGES

  !INIT DISLIN:
  CALL DISINI
  CALL ERRMOD('PROTOCOL','FILE')!SET MOD ON ERROR MESSAGE
  CALL COMPLX

  !SET SOME PAGE PROPERTIES: 
  !CALL PAGORG('BOTTOM')     !SETS THE ORIGIN AT "BOTTOM"-LEFT
  !CALL PAGHDR("","",3,0)    !PLOTS A HEADER WITH INFO
  !CALL PAGERA               !PLOTS THE BORDER

  !SET SOME AXIS PROPERTIES:
  !CALL AXSPOS(400,1800)     !SET THE AXIS-POSITION. don't need
  !CALL AXSLEN(1900,1400)     !SET THE LENGTH OF THE AXIS. weird units
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
  CALL TITLIN(adjustl(trim(pname)),1) !SET THE TITLE AND ITS LINE
  !AXIS:
  CALL NAME(trim(Xlabel),'X')
  CALL NAME(trim(Ylabel),'Y')

  !LEGEND:
  CALL LEGINI(CBUF,narray,10)
  CALL LEGLIN(CBUF,trim(Y1label),1)
  if(present(Y2)) CALL LEGLIN(CBUF,trim(Y2label),2)
  if(present(Y3)) CALL LEGLIN(CBUF,trim(Y3label),3)
  if(present(Y4)) CALL LEGLIN(CBUF,trim(Y4label),4)
  if(present(Y5)) CALL LEGLIN(CBUF,trim(Y5label),5)
  CALL LEGTIT('')

  !PLOT:
  CALL GRAF(real(xmin,4),real(xmax,4),real(xmin,4),real(dex,4),real(ymin,4),real(ymax,4),real(ymin,4),real(dey,4))
  CALL GRID(2,2)

  CALL TITLE

  if(trim(wlp)=="wl")  CALL INCMRK(0)
  if(trim(wlp)=="wlp") CALL INCMRK(1)
  if(trim(wlp)=="wp")  CALL INCMRK(-1)

  CALL THKCRV(5)

  !CALL CHNCRV("COLOR")
  call PLOT2D(X,Y1)
  if(present(Y2)) call PLOT2D(X,Y2,"RED")
  if(present(Y3)) call PLOT2D(X,Y3,"GREEN")
  if(present(Y4)) call PLOT2D(X,Y4,"BLUE")
  if(present(Y5)) call PLOT2D(X,Y5,"ORANGE")
  if(present(Y6)) call PLOT2D(X,Y6,"MAGENTA")

  CALL LEGEND(cbuf,3)
  CALL DISFIN
  !=========================================
contains
  subroutine plot2D(X,Y,color_)
    character(len=*),optional :: color_
    real(8),dimension(:) :: X,Y
    if(present(color_))CALL COLOR(color_)
    CALL CURVE(real(X),real(Y),size(X))
    CALL COLOR("FORE")
  end subroutine plot2D
end subroutine plot_dislin
!********************************************************************
!********************************************************************
!********************************************************************



!+-----------------------------------------------------------------+
!PROGRAM  : PLOT_DISLIN2D
!TYPE     : Subroutine
!PURPOSE  : 
!+-----------------------------------------------------------------+
subroutine plot_xmgrace(pname,X,Y1,Y1label, &
     Y2,Y2label,   &
     Y3,Y3label,   &
     Y4,Y4label,   &
     Y5,Y5label,   &
     Y6,Y6label,   &
     Xlabel,Ylabel)
  integer                       :: i,nsize,narray
  character(len=*)              :: pname,Xlabel,Ylabel
  character(len=64)             :: xname,yname
  real(8),dimension(:)          :: X
  real(8),dimension(:)          :: Y1
  character(len=*)              :: Y1label
  real(8),dimension(:),optional :: Y2,Y3,Y4,Y5,Y6
  character(len=*),optional     :: Y2label,Y3label,Y4label,Y5label,Y6label

  integer,save                  :: ITER=0
  character(len=3)              :: ITER_ 
  character(len=22),allocatable,dimension(:)   :: legend
  real(8),allocatable,dimension(:,:,:)         :: xtab


  xname=trim(Xlabel)
  yname=trim(Ylabel)

  narray=1
  if(present(Y2))  narray=narray+1
  if(present(Y3))  narray=narray+1
  if(present(Y4))  narray=narray+1
  if(present(Y5))  narray=narray+1
  if(present(Y6))  narray=narray+1

  allocate(legend(narray))
  legend(1)=Y1label
  if(present(Y2label))legend(2)=adjustl(trim(Y2label))
  if(present(Y3label))legend(3)=adjustl(trim(Y3label))
  if(present(Y4label))legend(4)=adjustl(trim(Y4label))
  if(present(Y5label))legend(5)=adjustl(trim(Y5label))
  if(present(Y6label))legend(6)=adjustl(trim(Y6label))

  nsize=size(Y1)
  if(present(Y2))nsize=min(size(Y1),size(Y2))
  if(present(Y3))nsize=min(size(Y1),size(Y2),size(Y3))
  if(present(Y4))nsize=min(size(Y1),size(Y2),size(Y3),size(Y4))
  if(present(Y5))nsize=min(size(Y1),size(Y2),size(Y3),size(Y4),size(Y5))
  if(present(Y6))nsize=min(size(Y1),size(Y2),size(Y3),size(Y4),size(Y5),size(Y6))
  if(size(X) < nsize) nsize=size(X)

  allocate(xtab(narray,nsize,3))
  xtab=0.d0
  xtab(:,:,3)=0.d0                        !Error bar 2b implemented
  forall(i=1:nsize) xtab(:,i,1)=X(i)      !All have same X-axis

  ITER=ITER+1
  do i=1,nsize
     xtab(1,i,2)=Y1(i)
     if(present(Y2))xtab(2,i,2)=Y2(i)
     if(present(Y3))xtab(3,i,2)=Y3(i)
     if(present(Y4))xtab(4,i,2)=Y4(i)
     if(present(Y5))xtab(5,i,2)=Y5(i)
     if(present(Y6))xtab(6,i,2)=Y6(i)
  enddo

  write(ITER_,'(i3)') ITER
  call dumpxmgrace(narray,nsize,xtab(1:narray,:,:),pname,xname,yname,legend, &
       adjustl(trim(pname))//"_iter"//trim(adjustl(trim(ITER_)))//'.agr',              &
       .false.,.false.,.false.,0.01d0)
  return
end subroutine plot_xmgrace
!********************************************************************
!********************************************************************
!********************************************************************
