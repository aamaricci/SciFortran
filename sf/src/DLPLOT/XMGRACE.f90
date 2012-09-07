!+-----------------------------------------------------------------+
!PROGRAM  : 
!TYPE     : Subroutine
!PURPOSE  : 
!+-----------------------------------------------------------------+
subroutine dumpxmgrace_(nset,n1,x,titlein,titlexin,titleyin,label,filename,logx,logy,&
     inset,fac2,inset_min,inset_max)
  implicit none
  integer                       :: nset
  character(len=*)              :: titlein,titlexin,titleyin
  character(len=22)             :: label(nset)
  character(len=len(titlein))   :: title
  character(len=len(titlexin))  :: titlex
  character(len=len(titleyin))  :: titley
  integer                       :: jj
  character(len=*)              :: filename
  character(len=len(filename))  :: filename2

  integer                  :: n1,i,j,k,jjj,i1,i2,kkk,kkk_
  real(8)                  :: x(nset,n1,3),fac2,fac
  real(8)                  :: r1,r2,rr1,rr2,rrr1,rrr2,rrrr1,rrrr2,tick1,tick2
  character(len=22)        :: char
  logical                  :: logx,logy,inset
  real(4),optional         :: inset_min,inset_max
  character(len(filename)) :: filename_dislin

  r1   =  (minval(x(:,:,1))) ; r2   =  (maxval(x(:,:,1)))  
  rr1  =  (minval(x(:,:,2))) ; rr2  =  (maxval(x(:,:,2)))  
  if(r1==r2)then
     r1=r1-0.25
     r2=r2+0.25
  endif
  if(rr1==rr2)then
     rr1=rr1-0.25
     rr2=rr2+0.25
  endif

  fac=r1+(r2-r1)*fac2
  if(logx) r1=max(1.d-9,r1)
  if(logy) rr1=max(1.d-9,rr1)

  i1=INT(min(1.d7,(r2-r1)/10.d0))   ;  i2=INT(min(1.d7,(rr2-rr1)/10.d0))
  tick1 = (r2-r1)/10.d0             ;  tick2 = (rr2-rr1)/10.d0

  rrr1 = r1 
  rrr2 = fac
  if(present(inset_min)) rrr1 = inset_min
  if(present(inset_max)) rrr2 = inset_max

  kkk=size(x(:,1,1))
  do j=1,n1
     if(maxval(x(:,j,1))>rrr2)then
        kkk=j; exit
     endif
  enddo
  kkk_=1
  do j=n1,1,-1
     if(minval(x(:,j,1))<rrr1)then
        kkk_=j; exit
     endif
  enddo

  rrrr1=minval(x(:,kkk_:kkk,2)); rrrr2=maxval(x(:,kkk_:kkk,2));
  title=adjustl(trim(titlein))
  titlex=adjustl(trim(titlexin))
  titley=adjustl(trim(titleyin))
  open(unit=14,file=adjustr(trim(filename)),status='new',err=102)
102 continue

  !introduction:

  write(14,*) '@version 50122'
  write(14,*) '@page size 720,523'
  write(14,*) '@page scroll 5%'
  write(14,*) '@page inout 5%'
  write(14,*) '@link page off'
  write(14,*) '@map font 8 to "Courier", "Courier"'
  write(14,*) '@map font 10 to "Courier-Bold", "Courier-Bold"'
  write(14,*) '@map font 11 to "Courier-BoldOblique", "Courier-BoldOblique"'
  write(14,*) '@map font 9 to "Courier-Oblique", "Courier-Oblique"'
  write(14,*) '@map font 4 to "Helvetica", "Helvetica"'
  write(14,*) '@map font 6 to "Helvetica-Bold", "Helvetica-Bold"'
  write(14,*) '@map font 7 to "Helvetica-BoldOblique", "Helvetica-BoldOblique"'
  write(14,*) '@map font 15 to "Helvetica-Narrow", "Helvetica-Narrow"'
  write(14,*) '@map font 16 to "Helvetica-Narrow-Bold", "Helvetica-Narrow-Bold"'
  write(14,*) '@map font 17 to "Helvetica-Narrow-BoldOblique", "Helvetica-Narrow-BoldOblique"'
  write(14,*) '@map font 18 to "Helvetica-Narrow-Oblique", "Helvetica-Narrow-Oblique"'
  write(14,*) '@map font 5 to "Helvetica-Oblique", "Helvetica-Oblique"'
  write(14,*) '@map font 20 to "NewCenturySchlbk-Bold", "NewCenturySchlbk-Bold"'
  write(14,*) '@map font 21 to "NewCenturySchlbk-BoldItalic", "NewCenturySchlbk-BoldItalic"'
  write(14,*) '@map font 22 to "NewCenturySchlbk-Italic", "NewCenturySchlbk-Italic"'
  write(14,*) '@map font 23 to "NewCenturySchlbk-Roman", "NewCenturySchlbk-Roman"'
  write(14,*) '@map font 24 to "Palatino-Bold", "Palatino-Bold"'
  write(14,*) '@map font 25 to "Palatino-BoldItalic", "Palatino-BoldItalic"'
  write(14,*) '@map font 26 to "Palatino-Italic", "Palatino-Italic"'
  write(14,*) '@map font 27 to "Palatino-Roman", "Palatino-Roman"'
  write(14,*) '@map font 12 to "Symbol", "Symbol"'
  write(14,*) '@map font 2 to "Times-Bold", "Times-Bold"'
  write(14,*) '@map font 3 to "Times-BoldItalic", "Times-BoldItalic"'
  write(14,*) '@map font 1 to "Times-Italic", "Times-Italic"'
  write(14,*) '@map font 0 to "Times-Roman", "Times-Roman"'
  write(14,*) '@map font 33 to "ZapfChancery-MediumItalic", "ZapfChancery-MediumItalic"'
  write(14,*) '@map font 13 to "ZapfDingbats", "ZapfDingbats"'
  write(14,*) '@map color 0 to (255, 255, 255), "white"'
  write(14,*) '@map color 1 to (0, 0, 0), "black"'
  write(14,*) '@map color 2 to (150, 0, 0), "red++"'
  write(14,*) '@map color 3 to (255, 0, 0), "red"'
  write(14,*) '@map color 4 to (255, 109, 0), "red--"'
  write(14,*) '@map color 5 to (0, 100, 0), "green++"'
  write(14,*) '@map color 6 to (30, 180, 0), "green"'
  write(14,*) '@map color 7 to (0, 255, 0), "green--"'
  write(14,*) '@map color 8 to (0, 0, 128), "blue++"'
  write(14,*) '@map color 9 to (0, 0, 255), "blue"'
  write(14,*) '@map color 10 to (0, 220, 255), "blue--"'
  write(14,*) '@map color 11 to (89, 0, 89), "magenta++"'
  write(14,*) '@map color 12 to (159, 0, 159), "magenta"'
  write(14,*) '@map color 13 to (255, 0, 255), "magenta--"'
  write(14,*) '@map color 14 to (77, 77, 77), "grey++"'
  write(14,*) '@map color 15 to (153, 153, 153), "grey"'
  write(14,*) '@map color 16 to (204, 204, 204), "grey--"'
  ! write(14,*) '@map color 0 to (255, 255, 255), "white"'
  ! write(14,*) '@map color 1 to (0, 0, 0), "black"'
  ! write(14,*) '@map color 2 to (255, 0, 0), "red"'
  ! write(14,*) '@map color 3 to (0, 255, 0), "green"'
  ! write(14,*) '@map color 4 to (0, 0, 255), "blue"'
  ! write(14,*) '@map color 5 to (255, 255, 0), "yellow"'
  ! write(14,*) '@map color 6 to (188, 143, 143), "brown"'
  ! write(14,*) '@map color 7 to (220, 220, 220), "grey"'
  ! write(14,*) '@map color 8 to (148, 0, 211), "violet"'
  ! write(14,*) '@map color 9 to (0, 255, 255), "cyan"'
  ! write(14,*) '@map color 10 to (255, 0, 255), "magenta"'
  ! write(14,*) '@map color 11 to (255, 165, 0), "orange"'
  ! write(14,*) '@map color 12 to (114, 33, 188), "indigo"'
  ! write(14,*) '@map color 13 to (103, 7, 72), "maroon"'
  ! write(14,*) '@map color 14 to (64, 224, 208), "turquoise"'
  ! write(14,*) '@map color 15 to (0, 139, 0), "green4"'
  write(14,*) '@reference date 0'
  write(14,*) '@date wrap off'
  write(14,*) '@date wrap year 1950'
  write(14,*) '@default linewidth 1.0'
  write(14,*) '@default linestyle 1'
  write(14,*) '@default color 1'
  write(14,*) '@default pattern 1'
  write(14,*) '@default font 0'
  write(14,*) '@default char size 1.000000'
  write(14,*) '@default symbol size 1.000000'
  write(14,*) '@default sformat "%.8g"'
  write(14,*) '@background color 0'
  write(14,*) '@page background fill on'
  write(14,*) '@timestamp off'
  write(14,*) '@timestamp 0.03, 0.03'
  write(14,*) '@timestamp color 1'
  write(14,*) '@timestamp rot 0'
  write(14,*) '@timestamp font 0'
  write(14,*) '@timestamp char size 1.000000'    
  write(14,*) '@timestamp def ""'

  write(14,*) '@g0 on'
  write(14,*) '@g0 hidden false'
  write(14,*) '@g0 type XY'
  write(14,*) '@g0 stacked false'
  write(14,*) '@g0 bar hgap 0.000000'
  write(14,*) '@g0 fixedpoint off'
  write(14,*) '@g0 fixedpoint type 0'
  write(14,*) '@g0 fixedpoint xy 0.000000, 0.000000'
  write(14,*) '@g0 fixedpoint format general general'
  write(14,*) '@g0 fixedpoint prec 6, 6'
  write(14,*) '@with g0'
  write(14,'(a11,4(f25.4,a1))') '@    world ',r1,',',rr1,',',r2,',',rr2
  write(14,*) '@    stack world 0, 0, 0, 0'
  write(14,*) '@    znorm 1'
  write(14,*) '@    view 0.150000, 0.150000, 1.150000, 0.850000'
  call plot_serie(.true.)
  call plot_set

  if(inset)then
     write(14,*) '@g1 on'
     write(14,*) '@g1 hidden false'
     write(14,*) '@g1 type XY'
     write(14,*) '@g1 stacked false'
     write(14,*) '@g1 bar hgap 0.000000'
     write(14,*) '@g1 fixedpoint off'
     write(14,*) '@g1 fixedpoint type 0'
     write(14,*) '@g1 fixedpoint xy 0.000000, 0.000000'
     write(14,*) '@g1 fixedpoint format general general'
     write(14,*) '@g1 fixedpoint prec 6, 6'
     write(14,*) '@with g1'
     write(14,'(a11,4(f22.4,a1))') '@    world ',rrr1,',',rrrr1,',',rrr2,',',rrrr2
     write(14,*) '@    stack world 0, 0, 0, 0'
     write(14,*) '@    znorm 1'
     write(14,*) '@    view 0.650000, 0.150000, 1.05, 0.450000'
     call plot_serie(.false.)
     call plot_set
  endif


  do jj=1,nset
     write(char,*) jj-1
     write(14,*) '@target G0.S'//ADJUSTL(TRIM(char))
     write(14,*) '@type xydy'
     do jjj=1,n1
        write(14,*) (x(jj,jjj,k),k=1,3)
     enddo
     write(14,*) ' & '
  enddo


  if(inset)then
     do jj=1,nset
        write(char,*) jj-1
        write(14,*) '@target G1.S'//ADJUSTL(TRIM(char))
        write(14,*) '@type xydy'
        do jjj=1,n1
           if(x(jj,jjj,1)<fac)then
              write(14,*) (x(jj,jjj,k),k=1,3)
           endif
        enddo
        write(14,*) ' & '
     enddo
  endif

  close(14,err=10)
10 continue

  return

  !****************************************************************************************
contains
  !****************************************************************************************

  !+-----------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : Subroutine
  !PURPOSE  : 
  !+-----------------------------------------------------------------+
  subroutine plot_serie(titles)
    logical :: titles

    if(titles)then
       write(14,*) '@    title "', ADJUSTL(TRIM(title)) ,'"'
    else
       write(14,*) '@    title ""'
    endif
    write(14,*) '@    title font 4'
    write(14,*) '@    title size 1.000000'
    write(14,*) '@    title color 1'

    write(14,*) '@    subtitle ""'
    write(14,*) '@    subtitle font 0'
    write(14,*) '@    subtitle size 1.000000'
    write(14,*) '@    subtitle color 1'

    if(logx)then
       write(14,*) '@    xaxes scale Logarithmic'
    else
       write(14,*) '@    xaxes scale Normal'
    endif
    if(logy)then 
       write(14,*) '@    yaxes scale Logarithmic'
    else
       write(14,*) '@    yaxes scale Normal'
    endif

    write(14,*) '@    xaxes invert off'
    write(14,*) '@    yaxes invert off'

    write(14,*) '@    xaxis  on'
    write(14,*) '@    xaxis  type zero false'
    write(14,*) '@    xaxis  offset 0.000000 , 0.000000'

    write(14,*) '@    xaxis  bar on'
    write(14,*) '@    xaxis  bar color 1'
    write(14,*) '@    xaxis  bar linestyle 1'
    write(14,*) '@    xaxis  bar linewidth 1.0'

    if(titles)then
       write(14,*) '@    xaxis  label "', ADJUSTL(TRIM(titlex)) ,'"'
    else
       write(14,*) '@    xaxis  label ""'
    endif

    write(14,*) '@    xaxis  label layout para'
    write(14,*) '@    xaxis  label place spec'
    write(14,*) '@    xaxis  label place 0.000000, 0.110000'
    write(14,*) '@    xaxis  label char size 1.00000'
    write(14,*) '@    xaxis  label font 4'
    write(14,*) '@    xaxis  label color 1'
    write(14,*) '@    xaxis  label place normal'
    write(14,*) '@    xaxis  tick on'
    write(14,*) '@    xaxis  tick major', tick1 !i1
    write(14,*) '@    xaxis  tick minor ticks 4'
    write(14,*) '@    xaxis  tick default 6'
    write(14,*) '@    xaxis  tick place rounded true'
    write(14,*) '@    xaxis  tick in'
    write(14,*) '@    xaxis  tick major size 1.640000'
    write(14,*) '@    xaxis  tick major color 1'
    write(14,*) '@    xaxis  tick major linewidth 1.0'
    write(14,*) '@    xaxis  tick major linestyle 1'
    write(14,*) '@    xaxis  tick major grid off'
    write(14,*) '@    xaxis  tick minor color 1'
    write(14,*) '@    xaxis  tick minor linewidth 1.0'
    write(14,*) '@    xaxis  tick minor linestyle 1'
    write(14,*) '@    xaxis  tick minor grid off'
    write(14,*) '@    xaxis  tick minor size 0.500000'
    write(14,*) '@    xaxis  ticklabel on'
    write(14,*) '@    xaxis  ticklabel format general'
    write(14,*) '@    xaxis  ticklabel prec 5'
    write(14,*) '@    xaxis  ticklabel formula ""'
    write(14,*) '@    xaxis  ticklabel append ""'
    write(14,*) '@    xaxis  ticklabel prepend ""'
    write(14,*) '@    xaxis  ticklabel angle 0'
    write(14,*) '@    xaxis  ticklabel skip 0'
    write(14,*) '@    xaxis  ticklabel stagger 0'
    write(14,*) '@    xaxis  ticklabel place normal'
    write(14,*) '@    xaxis  ticklabel offset auto'
    write(14,*) '@    xaxis  ticklabel offset 0.000000 , 0.010000'
    write(14,*) '@    xaxis  ticklabel start type auto'
    write(14,*) '@    xaxis  ticklabel start 0.000000'
    write(14,*) '@    xaxis  ticklabel stop type auto'
    write(14,*) '@    xaxis  ticklabel stop 0.000000'
    write(14,*) '@    xaxis  ticklabel char size 1.0000'
    write(14,*) '@    xaxis  ticklabel font 0'
    write(14,*) '@    xaxis  ticklabel color 1'
    write(14,*) '@    xaxis  tick place both'
    write(14,*) '@    xaxis  tick spec type none'
    write(14,*) '@    yaxis  on'
    write(14,*) '@    yaxis  type zero false'
    write(14,*) '@    yaxis  offset 0.000000 , 0.000000'
    write(14,*) '@    yaxis  bar on'
    write(14,*) '@    yaxis  bar color 1'
    write(14,*) '@    yaxis  bar linestyle 1'
    write(14,*) '@    yaxis  bar linewidth 1.0'

    if(titles)then
       write(14,*) '@    yaxis  label "', ADJUSTL(TRIM(titley)) ,'"'
    else
       write(14,*) '@    yaxis  label " "'
    endif

    write(14,*) '@    yaxis  label layout para'
    write(14,*) '@    yaxis  label place spec'
    write(14,*) '@    yaxis  label place 0.000000, 0.110000'
    write(14,*) '@    yaxis  label char size 1.000000'
    write(14,*) '@    yaxis  label font 4'
    write(14,*) '@    yaxis  label color 1'
    write(14,*) '@    yaxis  label place normal'
    write(14,*) '@    yaxis  tick on'
    write(14,*) '@    yaxis  tick major ', tick2 !i2
    write(14,*) '@    yaxis  tick minor ticks 4'
    write(14,*) '@    yaxis  tick default 6'
    write(14,*) '@    yaxis  tick place rounded true'
    write(14,*) '@    yaxis  tick in'
    write(14,*) '@    yaxis  tick major size 1.870000'
    write(14,*) '@    yaxis  tick major color 1'
    write(14,*) '@    yaxis  tick major linewidth 1.0'
    write(14,*) '@    yaxis  tick major linestyle 1'
    write(14,*) '@    yaxis  tick major grid off'
    write(14,*) '@    yaxis  tick minor color 1'
    write(14,*) '@    yaxis  tick minor linewidth 1.0'
    write(14,*) '@    yaxis  tick minor linestyle 1'
    write(14,*) '@    yaxis  tick minor grid off'
    write(14,*) '@    yaxis  tick minor size 0.500000'
    write(14,*) '@    yaxis  ticklabel on'
    write(14,*) '@    yaxis  ticklabel format general'
    write(14,*) '@    yaxis  ticklabel prec 5'
    write(14,*) '@    yaxis  ticklabel formula ""'
    write(14,*) '@    yaxis  ticklabel append ""'
    write(14,*) '@    yaxis  ticklabel prepend ""'
    write(14,*) '@    yaxis  ticklabel angle 0'
    write(14,*) '@    yaxis  ticklabel skip 0'
    write(14,*) '@    yaxis  ticklabel stagger 0'
    write(14,*) '@    yaxis  ticklabel place normal'
    write(14,*) '@    yaxis  ticklabel offset auto'
    write(14,*) '@    yaxis  ticklabel offset 0.000000 , 0.010000'
    write(14,*) '@    yaxis  ticklabel start type auto'
    write(14,*) '@    yaxis  ticklabel start 0.000000'
    write(14,*) '@    yaxis  ticklabel stop type auto'
    write(14,*) '@    yaxis  ticklabel stop 0.000000'
    write(14,*) '@    yaxis  ticklabel char size 1.0000'
    write(14,*) '@    yaxis  ticklabel font 0'
    write(14,*) '@    yaxis  ticklabel color 1'
    write(14,*) '@    yaxis  tick place both'
    write(14,*) '@    yaxis  tick spec type none'
    write(14,*) '@    altxaxis  off'
    write(14,*) '@    altyaxis  off'
    write(14,*) '@    legend on'
    write(14,*) '@    legend loctype view'
    write(14,*) '@    legend 0.85, 0.8'
    write(14,*) '@    legend box color 1'
    write(14,*) '@    legend box pattern 1'
    write(14,*) '@    legend box linewidth 1.0'
    write(14,*) '@    legend box linestyle 1'
    write(14,*) '@    legend box fill color 0'
    write(14,*) '@    legend box fill pattern 1'
    write(14,*) '@    legend font 4'
    write(14,*) '@    legend char size 1.000000'
    write(14,*) '@    legend color 1'
    write(14,*) '@    legend length 4'
    write(14,*) '@    legend vgap 1'
    write(14,*) '@    legend hgap 1'
    write(14,*) '@    legend invert false'
    write(14,*) '@    frame type 0'
    write(14,*) '@    frame linestyle 0'
    write(14,*) '@    frame linewidth 1.0'
    write(14,*) '@    frame color 1'
    write(14,*) '@    frame pattern 1'
    write(14,*) '@    frame background color 0'
    write(14,*) '@    frame background pattern 0'
  end subroutine plot_serie
  !********************************************************************
  !********************************************************************
  !********************************************************************


  !+-----------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : Subroutine
  !PURPOSE  : 
  !+-----------------------------------------------------------------+
  subroutine plot_set
    do jj=1,nset
       write(char,*) jj-1
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' hidden false'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' type xydy'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' symbol', jj          !2+jj 
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' symbol size 1.5000'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' symbol color ', jj   !jj+1 
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' symbol pattern 1'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' symbol fill color ', 1+jj
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' symbol fill pattern 0'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' symbol linewidth 1.0'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' symbol linestyle 1'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' symbol char 65'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' symbol char font 0'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' symbol skip 0'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' line type 1'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' line linestyle 1'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' line linewidth 3.0'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' line color ',jj    ! 1'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' line pattern 1'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' baseline type 0'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' baseline off'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' dropline off'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' fill type 0'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' fill rule 0'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' fill color 1'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' fill pattern 1'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' avalue off'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' avalue type 2'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' avalue char size 1.000000'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' avalue font 0'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' avalue color 1'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' avalue rot 0'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' avalue format general'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' avalue prec 3'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' avalue prepend ""'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' avalue append ""'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' avalue offset 0.000000 , 0.000000'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' errorbar on'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' errorbar place both'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' errorbar color 1'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' errorbar pattern 1'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' errorbar size 1.000000'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' errorbar linewidth 1.0'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' errorbar linestyle 1'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' errorbar riser linewidth 1.0'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' errorbar riser linestyle 1'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' errorbar riser clip off'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' errorbar riser clip length 0.100000'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' comment " "'
       write(14,*) '@    s'//ADJUSTL(TRIM(char))//' legend  "', label(jj),'"'

    enddo
  end subroutine plot_set
  !********************************************************************
  !********************************************************************
  !********************************************************************
end subroutine dumpxmgrace_
!********************************************************************
!********************************************************************
!********************************************************************



!+-----------------------------------------------------------------+
!PROGRAM  : 
!TYPE     : Subroutine
!PURPOSE  : 
!+-----------------------------------------------------------------+
subroutine dumpxmgrace__(nset,n2,n1,x,title,titlex,titley,label,filename,logx,logy,inset,&
     fac2,inset_min,inset_max)
  implicit none
  integer             :: nset
  character(len=*)    :: title,titlex,titley,label(nset,n2)
  character(len=22)   :: labell(nset*n2) 
  integer             :: jj,n2
  character(len=*)    :: filename
  integer             :: n1,i,j,k,jjj,i1,i2,kkk
  real(8)             :: x(nset,n2,n1,3),fac2,fac,xx(nset*n2,n1,3)
  logical             :: logx,logy,inset
  real(4),optional    :: inset_min,inset_max
  jj=0
  do i=1,nset
     do j=1,n2
        jj=jj+1
        xx(jj,:,:)=x(i,j,:,:)
        labell(jj)=label(i,j)
     enddo
  enddo
  call dumpxmgrace_(nset*n2,n1,xx,title,titlex,titley,labell,filename,logx,logy,inset,&
       fac2=fac2,inset_min=inset_min,inset_max=inset_max)
end subroutine dumpxmgrace__
!********************************************************************
!********************************************************************
!********************************************************************
