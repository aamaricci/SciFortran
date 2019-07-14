
subroutine p_deigh_simple(A,W,Nblock,blacs_end,method,jobz,uplo,vl,vu,il,iu,tol)
  real(8),dimension(:,:),intent(inout)       :: A ! M v = E v/v(i,j) = ith component of jth vec.
  real(8),dimension(size(A,2)),intent(inout) :: W ! eigenvalues
  integer                                    :: Nblock
  integer,optional                           :: blacs_end
  integer                                    :: Nb
  character(len=*),optional                  :: method
  character(len=1),optional                  :: jobz,uplo
  character(len=1)                           :: jobz_,uplo_,range
  integer                                    :: blacs_end_
  character(len=20)                          :: method_
  real(8),optional                           :: vl,vu,tol
  integer,optional                           :: il,iu
  real(8)                                    :: vL_,vU_,tol_
  integer                                    :: iL_,iU_
  integer                                    :: Ns
  integer                                    :: Qrows,Qcols
  integer                                    :: i,j,lda,info,ldz
  integer                                    :: lwork,liwork,mW,mZ
  real(8),dimension(:),allocatable           :: work,Rwork,Gap
  integer                                    :: guess_liwork(1)
  real(8)                                    :: guess_lwork(1)
  logical                                    :: boolV,boolI
  real(8),dimension(:,:),allocatable         :: Z
  integer,dimension(:),allocatable           :: Isuppz,Ifail,Iclustr
  integer,dimension(:),allocatable           :: Iwork
  !
  integer,external                           :: numroc,indxG2L,indxG2P,indxL2G
  real(8),external                           :: dlamch      
  !
  real(8),dimension(:,:),allocatable         :: A_loc,Z_loc
  integer                                    :: p_size
  integer                                    :: p_Nx,p_Ny
  integer                                    :: p_context
  integer                                    :: rank,rankX,rankY
  integer                                    :: Nrow,Ncol
  integer                                    :: sendR,sendC
  integer                                    :: myi,myj,unit,irank
  integer,dimension(9)                       :: descA,descAloc,descZloc
  real(8)                                    :: t_stop,t_start
  logical                                    :: master
  !
  method_='dsyevr';if(present(method))method_=trim(method)
  jobz_='V'      ;if(present(jobz))jobz_=jobz
  uplo_='L'      ;if(present(uplo))uplo_=uplo
  vl_  = 1d0     ;if(present(vL))vL_=vL
  vu_  = 1d0     ;if(present(vU))vU_=vU
  iL_  = 1       ;if(present(iL))iL_=iL
  iU_  = 1       ;if(present(iU))iU_=iU
  tol_ = dlamch('s')     ;if(present(tol))tol_=tol
  blacs_end_=0   ;if(present(blacs_end))blacs_end_=blacs_end
  !
  W=0d0
  !
  range='A'
  boolV=present(vL).AND.present(vU)
  boolI=present(iL).OR.present(iU)
  if(boolV.and.boolI)stop "vX and iX arguments both present. Can not set RANGE"
  if(boolV)range='V'
  if(boolI)range='I'
  !
  Ns    = max(1,size(A,1))
  if(any(shape(A)/=[Ns,Ns]))stop "my_eighD error: A has illegal shape"
  !
  !INIT SCALAPACK TREATMENT:
  !
  !< Initialize BLACS processor grid (like MPI)
  call blacs_setup(rank,p_size)  ![id, size]
  master = (rank==0)
  do i=1,int( sqrt( dble(p_size) ) + 1 )
     if(mod(p_size,i)==0) p_Nx = i
  end do
  p_Ny = p_size/p_Nx
  !
  !< Init context with p_Nx,p_Ny procs
  call sl_init(p_context,p_Nx,p_Ny)
  !
  !< Get coordinate of the processes
  call blacs_gridinfo( p_context, p_Nx, p_Ny, rankX, rankY)
  !
  if(rankX<0.AND.rankY<0)goto 100
  !
  Nb = Nblock
  !
  Qrows = numroc(Ns, Nb, rankX, 0, p_Nx)
  Qcols = numroc(Ns, Nb, rankY, 0, p_Ny)
  !
  if(master)then
     unit = 513! + rank
     open(unit,file="p_eigh.info")
     write(unit,"(A20,I8,A5,I8)")"Grid=",p_Nx,"x",p_Ny
     write(unit,"(A20,I2,I8,A5,I8)")"Qrows x Qcols=",rank,Qrows,"x",Qcols
  endif
  !
  !< allocate local distributed A
  allocate(A_loc(Qrows,Qcols))
  call descinit( descA, Ns, Ns, Nb, Nb, 0, 0, p_context, Qrows, info )
  call descinit( descAloc, Ns, Ns, Nb, Nb, 0, 0, p_context, Qrows, info )
  !
  !< Distribute A
  if(master)call cpu_time(t_start)
  do myi=1,Qrows
     i  = indxL2G(myi,Nblock,rankX,0,p_Nx)
     do myj=1,Qcols
        j  = indxL2G(myj,Nblock,rankY,0,p_Ny)
        A_loc(myi,myj) = A(i,j)
     enddo
  enddo
  if(master)call cpu_time(t_stop)
  if(master)write(unit,"(A20,F21.12)")"Time Distribute A:",t_stop-t_start
  !
  !< Allocate distributed eigenvector matrix
  allocate(Z_loc(Qrows,Qcols));Z_loc=0d0
  call descinit( descZloc, Ns, Ns, Nb, Nb, 0, 0, p_context, Qrows, info )
  !
  if(master)write(unit,"(A20,A21)")"Using Method:",method_
  if(master)call cpu_time(t_start)
  select case(method_)
  case default
     call  PDSYEVR(jobz_,range,uplo_,&
          Ns,A_loc,1,1,descAloc,vl_,vu_,il_,iu_,mW,mZ,W,Z_loc,1,1,descZloc,&
          guess_lwork,-1,guess_liwork,-1,info)
     lwork = guess_lwork(1)
     liwork= guess_liwork(1)
     allocate(work(lwork))
     allocate(iwork(liwork))
     call  PDSYEVR(jobz_,range,uplo_,&
          Ns,A_loc,1,1,descAloc,vl_,vu_,il_,iu_,mW,mZ,W,Z_loc,1,1,descZloc,&
          work, lwork, iwork, liwork,info)
     !
  case ("dsyev")
     call PDSYEV(jobz_,uplo_,&
          Ns,A_loc,1,1,descAloc,W,Z_loc,1,1,descZloc,&
          guess_lwork,-1,info)
     lwork = guess_lwork(1)
     allocate(work(lwork))
     call PDSYEV( jobz_,uplo_,&
          Ns,A_loc,1,1,descAloc,W,Z_loc,1,1,descZloc,&
          work,lwork,info)
     !
     !
  case ("dsyevd")
     call PDSYEVD(jobz_,uplo_,&
          Ns,A_loc,1,1,descAloc,W,Z_loc,1,1,descZloc,&
          guess_lwork,-1,guess_liwork,-1,INFO)
     lwork = guess_lwork(1)
     liwork= guess_liwork(1)
     allocate(work(lwork))
     allocate(iwork(liwork))
     call PDSYEVD(jobz_,uplo_,&
          Ns,A_loc,1,1,descAloc,W,Z_loc,1,1,descZloc,&
          work,lwork,iwork,liwork,info)
     !
  case ("dsyevx")
     allocate(Ifail(Ns))
     allocate(Iclustr(2*p_Nx*p_Ny))
     allocate(Gap(p_Nx*p_Ny))
     call PDSYEVX(jobz_,range,uplo_,&
          Ns,A_loc,1,1,descAloc,vl_,vu_,il_,iu_,0d0,mW, mZ,W,-1d0,Z_loc,1,1,descAloc,&
          guess_lwork,-1,guess_liwork,-1,ifail,iclustr,gap,info)
     lwork = guess_lwork(1)
     liwork= guess_liwork(1)
     allocate(work(lwork))
     allocate(iwork(liwork))
     call PDSYEVX(jobz_,range,uplo_,&
          Ns,A_loc,1,1,descAloc,vl_,vu_,il_,iu_,0d0,mW,mZ,W,-1d0,Z_loc,1,1,descZloc,&
          work,lwork,iwork,liwork,ifail,iclustr,gap,info)
  end select
  if(master)call cpu_time(t_stop)
  if(master)write(unit,"(A20,F21.12)")"Time diag A:",t_stop-t_start
  !
  if(jobz_=='V')then
     A=0d0
     if(master)call cpu_time(t_start)
     do i=1,Ns,Nb
        Nrow = Nb ; if(Ns-i<Nb-1)Nrow=Ns-i+1!;if(Nrow==0)Nrow=1
        do j=1,Ns,Nb
           Ncol = Nb ; if(Ns-j<Nb-1)Ncol=Ns-j+1!;if(Ncol==0)Ncol=1
           call infog2l(i,j,descA, p_Nx, p_Ny, rankX, rankY, myi, myj, SendR, SendC)
           if(rankX==SendR .AND. rankY==SendC)then
              call dgesd2d(p_context,Nrow,Ncol,Z_loc(myi,myj),Qrows,0,0)
           endif
           if(rank==0)then
              call dgerv2d(p_context,Nrow,Ncol,A(i,j),Ns,SendR,SendC)
           endif
        enddo
     enddo
     if(master)call cpu_time(t_stop)
     if(master)write(unit,"(A20,F21.12)")"Time gather Z:",t_stop-t_start
  endif
  !
  if(master)close(unit)
  call blacs_gridexit(p_context)
100 continue
  call blacs_exit(blacs_end_)
  return
  !
end subroutine p_deigh_simple




subroutine p_zeigh_simple(A,W,Nblock,blacs_end,method,jobz,uplo,vl,vu,il,iu,tol)
  complex(8),dimension(:,:),intent(inout)    :: A ! M v = E v/v(i,j) = ith component of jth vec.
  real(8),dimension(size(A,2)),intent(inout) :: W ! eigenvalues
  integer                                    :: Nblock
  integer,optional                           :: blacs_end
  integer                                    :: Nb
  character(len=*),optional                  :: method
  character(len=1),optional                  :: jobz,uplo
  character(len=1)                           :: jobz_,uplo_,range
  character(len=20)                          :: method_
  integer                                    :: blacs_end_
  real(8),optional                           :: vl,vu,tol
  integer,optional                           :: il,iu
  real(8)                                    :: vL_,vU_,tol_
  integer                                    :: iL_,iU_
  integer                                    :: Ns
  integer                                    :: Qrows,Qcols
  integer                                    :: i,j,lda,info,ldz
  integer                                    :: lwork,liwork,lrwork,mW,mZ
  complex(8),dimension(:),allocatable        :: Work,RRwork
  real(8),dimension(:),allocatable           :: Rwork
  integer,dimension(:),allocatable           :: Iwork
  complex(8)                                 :: guess_lwork(1)
  complex(8)                                 :: guess_lrrwork(1)
  real(8)                                    :: guess_lrwork(1)
  integer                                    :: guess_liwork(1)
  logical                                    :: boolV,boolI
  real(8),dimension(:,:),allocatable         :: Z
  integer,dimension(:),allocatable           :: Isuppz
  integer,dimension(:),allocatable           :: Ifail
  integer,dimension(:),allocatable           :: Iclustr
  real(8),dimension(:),allocatable           :: Gap
  !
  integer,external                           :: numroc,indxG2L,indxG2P,indxL2G
  real(8),external                           :: dlamch      
  !
  complex(8),dimension(:,:),allocatable      :: A_loc,Z_loc
  integer                                    :: p_size
  integer                                    :: p_Nx,p_Ny
  integer                                    :: p_context
  integer                                    :: rank,rankX,rankY
  integer                                    :: Nrow,Ncol
  integer                                    :: sendR,sendC
  integer                                    :: myi,myj,unit,irank
  integer,dimension(9)                       :: descA,descAloc,descZloc
  real(8)                                    :: t_stop,t_start
  logical                                    :: master
  !
  method_='zheevr';if(present(method))method_=trim(method)
  !
  jobz_='V'      ;if(present(jobz))jobz_=jobz
  uplo_='L'      ;if(present(uplo))uplo_=uplo
  vl_  = 1d0     ;if(present(vL))vL_=vL
  vu_  = 1d0     ;if(present(vU))vU_=vU
  iL_  = 1       ;if(present(iL))iL_=iL
  iU_  = 1       ;if(present(iU))iU_=iU
  tol_ = dlamch('s')     ;if(present(tol))tol_=tol
  blacs_end_=0   ;if(present(blacs_end))blacs_end_=blacs_end
  !
  !>DEBUG
  if(method_=="zheev" .OR. method_=="zheevd") then
     write(0,"(A)")"WARNING p_eigh: pzheev & pzheevd are bugged. Fall on pzheevr."
     method_='zheevr'
  endif
  if(uplo_=='U') then
     write(0,"(A)")"WARNING p_eigh: pzheevr OR pzheevx are bugged. Only work with with uplo=L"
     uplo_  ='L'
  endif
  !<DEBUG
  !
  W=0d0
  !
  range='A'
  boolV=present(vL).AND.present(vU)
  boolI=present(iL).OR.present(iU)
  if(boolV.and.boolI)stop "vX and iX arguments both present. Can not set RANGE"
  if(boolV)range='V'
  if(boolI)range='I'
  Ns    = max(1,size(A,1))
  if(any(shape(A)/=[Ns,Ns]))stop "my_eighD error: A has illegal shape"
  !
  !INIT SCALAPACK TREATMENT:
  !< Initialize BLACS processor grid (like MPI)
  call blacs_setup(rank,p_size)  ![id, size]
  master = (rank==0)
  do i=1,int( sqrt( dble(p_size) ) + 1 )
     if(mod(p_size,i)==0) p_Nx = i
  end do
  p_Ny = p_size/p_Nx
  !
  !< Init context with p_Nx,p_Ny procs
  call sl_init(p_context,p_Nx,p_Ny)
  !
  !< Get coordinate of the processes
  call blacs_gridinfo( p_context, p_Nx, p_Ny, rankX, rankY)
  !
  if(rankX<0.AND.rankY<0)goto 200
  !
  Nb = Nblock
  !
  Qrows = numroc(Ns, Nb, rankX, 0, p_Nx)
  Qcols = numroc(Ns, Nb, rankY, 0, p_Ny)
  !
  call descinit( descA, Ns, Ns, Nb, Nb, 0, 0, p_context, Qrows, info )
  call descinit( descAloc, Ns, Ns, Nb, Nb, 0, 0, p_context, Qrows, info )
  !
  if(master)then
     unit = 513
     open(unit,file="p_eigh.info")
     write(unit,"(A20,I8,A5,I8)")"Grid=",p_Nx,"x",p_Ny
     write(unit,"(A20,I2,I8,A5,I8)")"Qrows x Qcols=",rank,Qrows,"x",Qcols
  endif
  !
  !< allocate local distributed A
  allocate(A_loc(Qrows,Qcols));A_loc=zero
  !
  !< Distribute A
  if(master)call cpu_time(t_start)
  do myi=1,Qrows
     i  = indxL2G(myi,Nblock,rankX,0,p_Nx)
     do myj=1,Qcols
        j  = indxL2G(myj,Nblock,rankY,0,p_Ny)
        A_loc(myi,myj) = A(i,j)
     enddo
  enddo
  if(master)call cpu_time(t_stop)
  if(master)write(unit,"(A20,F21.12)")"Time Distribute A:",t_stop-t_start
  !
  !< Allocate distributed eigenvector matrix
  allocate(Z_loc(Qrows,Qcols));Z_loc=zero
  call descinit( descZloc, Ns, Ns, Nb, Nb, 0, 0, p_context, Qrows, info )
  !
  if(master)write(unit,"(A20,A21)")"Using Method:",method_
  if(master)call cpu_time(t_start)
  select case(method_)
  case default
     call  PZHEEVR(jobz_,range,uplo_,&
          Ns,A_loc,1,1,descAloc,vl_,vu_,il_,iu_,mW,mZ,W,Z_loc,1,1,descZloc,&
          guess_lwork,-1,guess_lrwork,-1,guess_liwork,-1,info)
     lwork = guess_lwork(1) ; lrwork= guess_lrwork(1)  ; liwork= guess_liwork(1)
     allocate(work(lwork))
     allocate(rwork(lrwork))
     allocate(iwork(liwork))
     call  PZHEEVR(jobz_,range,uplo_,&
          Ns,A_loc,1,1,descAloc,vl_,vu_,il_,iu_,mW,mZ,W,Z_loc,1,1,descZloc,&
          work, lwork, rwork, lrwork, iwork, liwork,info)
     !
     !>>>ACTHUNG<< BUGGED DO NOT USE: --> zheevr, w UPLO='L'
     ! case ("zheev")              !Bugged
     !    call PZHEEV(jobz_,uplo_,&
     !         Ns,A_loc,1,1,descAloc,W,Z_loc,1,1,descZloc,&
     !         guess_lwork,-1,guess_lrwork,-1,info)
     !    lwork = guess_lwork(1) ; lrwork = guess_lrwork(1)
     !    allocate(work(lwork))
     !    allocate(rwork(lrwork))
     !    call PZHEEV(jobz_,uplo_,&
     !         Ns,A_loc,1,1,descAloc,W,Z_loc,1,1,descZloc,&
     !         work,lwork,rwork,lrwork,info)
     !
     !>>>ACTHUNG<< BUGGED DO NOT USE: --> zheevr, w UPLO='L'
     ! case ("zheevd")
     !    call PZHEEVD('V',uplo_,&
     !         Ns,A_loc,1,1,descAloc,W,Z_loc,1,1,descZloc,&
     !         guess_lwork,-1,guess_lrwork,-1,guess_liwork,-1,INFO)
     !    lwork = guess_lwork(1) ; lrwork= guess_lrwork(1) ; liwork= guess_liwork(1)
     !    allocate(work(lwork))
     !    allocate(rwork(lrwork))
     !    allocate(iwork(liwork))
     !    call PZHEEVD('V',uplo_,&
     !         Ns,A_loc,1,1,descAloc,W,Z_loc,1,1,descZloc,&
     !         work,lwork,rwork,lrwork,iwork,liwork,info)
     !
  case ("zheevx")
     allocate(Ifail(Ns))
     allocate(Iclustr(2*p_Nx*p_Ny))
     allocate(Gap(p_Nx*p_Ny))
     call PZHEEVX(jobz_,range,uplo_,&
          Ns,A_loc,1,1,descAloc,vl_,vu_,il_,iu_,0d0,mW, mZ,W,-1d0,Z_loc,1,1,descAloc,&
          guess_lwork,-1,guess_lrwork,-1,guess_liwork,-1,ifail,iclustr,gap,info)
     lwork = guess_lwork(1)
     lrwork= guess_lrwork(1)
     liwork= guess_liwork(1)
     allocate(work(lwork))
     allocate(rwork(lrwork))
     allocate(iwork(liwork))
     call PZHEEVX(jobz_,range,uplo_,&
          Ns,A_loc,1,1,descAloc,vl_,vu_,il_,iu_,0d0,mW,mZ,W,-1d0,Z_loc,1,1,descZloc,&
          work,lwork,rwork,lrwork,iwork,liwork,ifail,iclustr,gap,info)
  end select
  if(master)call cpu_time(t_stop)
  if(master)write(unit,"(A20,F21.12)")"Time diag A:",t_stop-t_start
  !
  if(jobz_=='V')then
     A=zero
     if(master)call cpu_time(t_start)
     do i=1,Ns,Nb
        Nrow = Nb ; if(Ns-i<Nb-1)Nrow=Ns-i+1!;if(Nrow==0)Nrow=1
        do j=1,Ns,Nb
           Ncol = Nb ; if(Ns-j<Nb-1)Ncol=Ns-j+1!;if(Ncol==0)Ncol=1
           call infog2l(i,j,descA, p_Nx, p_Ny, rankX, rankY, myi, myj, SendR, SendC)
           if(rankX==SendR .AND. rankY==SendC)then
              call zgesd2d(p_context,Nrow,Ncol,Z_loc(myi,myj),Qrows,0,0)
           endif
           if(rank==0)then
              call zgerv2d(p_context,Nrow,Ncol,A(i,j),Ns,SendR,SendC)
           endif
        enddo
     enddo
     if(master)call cpu_time(t_stop)
     if(master)write(unit,"(A20,F21.12)")"Time gather Z:",t_stop-t_start
  endif
  !
  if(master)close(unit)
  call blacs_gridexit(p_context)
200 continue
  call blacs_exit(blacs_end_)
  return
end subroutine p_zeigh_simple

