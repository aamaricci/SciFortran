subroutine d_distribute_BLACS(M,Mloc,descMloc,unit)
  USE SF_MPI
  include 'mpif.h'
  real(8),dimension(:,:),intent(in)    :: M
  real(8),dimension(:,:),intent(inout) :: Mloc
  integer,dimension(9),intent(in)      :: descMloc
  integer,optional                     :: unit
  integer                              :: i,j
  integer                              :: myi,myj
  integer                              :: Nbi,Nbj
  integer                              :: Qrows,Qcols
  integer                              :: rankX,rankY
  integer                              :: pNx,pNy
  integer                              :: context
  real(8)                              :: t_start,t_stop
  integer,external                     :: indxG2L,indxL2G
  !
  context = descMloc(2)
  call blacs_gridinfo(context, pNx, pNy, rankX, rankY)
  !
  Qrows  = size(Mloc,1);if(Qrows/=descMloc(9))stop "Distribute_BLACS error: Qrows!=descM(9)"
  Qcols  = size(Mloc,2)
  Nbi = descMloc(5)
  Nbj = descMloc(6)
  !
  if(rankX==0.AND.rankY==0)call cpu_time(t_start)
  do myj=1,Qcols
     j  = indxL2G(myj,Nbj,rankY,0,pNy)
     do myi=1,Qrows
        i  = indxL2G(myi,Nbi,rankX,0,pNx)
        Mloc(myi,myj) = M(i,j)
     enddo
  enddo
  if(rankX==0.AND.rankY==0)call cpu_time(t_stop)
  if(present(unit))then
     if(rankX==0.AND.rankY==0)write(unit,"(A20,F21.12)")"Time Distribute :",t_stop-t_start
  endif
  return
end subroutine D_Distribute_BLACS

subroutine Z_distribute_BLACS(M,Mloc,descMloc,unit)
  USE SF_MPI
  include 'mpif.h'
  complex(8),dimension(:,:),intent(in)    :: M
  complex(8),dimension(:,:),intent(out)   :: Mloc
  integer,dimension(9),intent(in)         :: descMloc
  integer,optional                        :: unit
  integer                                 :: i,j
  integer                                 :: myi,myj
  integer                                 :: Nbi,Nbj
  integer                                 :: Qrows,Qcols
  integer                                 :: rankX,rankY
  integer                                 :: pNx,pNy
  integer                                 :: context
  real(8)                                 :: t_start,t_stop
  integer,external                        :: indxG2L,indxL2G
  !
  context = descMloc(2)
  call blacs_gridinfo(context, pNx, pNy, rankX, rankY)
  !
  Qrows  = size(Mloc,1);if(Qrows/=descMloc(9))stop "Distribute_BLACS error: Qrows!=descM(9)"
  Qcols  = size(Mloc,2)
  Nbi = descMloc(5)
  Nbj = descMloc(6)
  !
  if(rankX==0.AND.rankY==0)call cpu_time(t_start)
  do myj=1,Qcols
     j  = indxL2G(myj,Nbj,rankY,0,pNy)
     do myi=1,Qrows
        i  = indxL2G(myi,Nbi,rankX,0,pNx)
        Mloc(myi,myj) = M(i,j)
     enddo
  enddo
  if(rankX==0.AND.rankY==0)call cpu_time(t_stop)
  if(present(unit))then
     if(rankX==0.AND.rankY==0)write(unit,"(A20,F21.12)")"Time Distribute :",t_stop-t_start
  endif
  return
end subroutine Z_Distribute_BLACS











subroutine D_Gather_BLACS(Mloc,M,descM,unit)
  USE SF_MPI
  include 'mpif.h'
  real(8),dimension(:,:),intent(in)    :: Mloc
  real(8),dimension(:,:),intent(inout) :: M
  integer,dimension(9),intent(in)      :: descM
  integer,optional                     :: unit
  integer                              :: i,j
  integer                              :: myi,myj
  integer                              :: Ni,Nj
  integer                              :: Nbi,Nbj
  integer                              :: Nrow,Ncol,Qrows
  integer                              :: rankX,rankY
  integer                              :: pNx,pNy
  integer                              :: sendR,sendC
  real(8)                              :: t_start,t_stop
  integer                              :: context
  !
  context = descM(2)
  call blacs_gridinfo( context, pNx, pNy, rankX, rankY)
  !
  Ni  = size(M,1)
  Nj  = size(M,2)
  Nbi = descM(5)
  Nbj = descM(6)
  Qrows=descM(9) ; if(Qrows/=size(Mloc,1))stop "Gather_BLACS error: Qrows != size(Mloc,1)"
  !
  if(rankX==0.AND.rankY==0)call cpu_time(t_start)
  do j=1,Nj,Nbj
     Ncol = Nbj ; if(Nj-j<Nbj-1)Ncol=Nj-j+1
     do i=1,Ni,Nbi
        Nrow = Nbi ; if(Ni-i<Nbi-1)Nrow=Ni-i+1
        call infog2l(i,j,descM, pNx, pNy, rankX, rankY, myi, myj, SendR, SendC)
        if(rankX==SendR .AND. rankY==SendC)then
           call dgesd2d(context,Nrow,Ncol,Mloc(myi,myj),Qrows,0,0)
        endif
        if(rankX==0 .AND. rankY==0)then
           call dgerv2d(context,Nrow,Ncol,M(i,j),Ni,SendR,SendC)
        endif
     enddo
  enddo
  if(rankX==0.AND.rankY==0)call cpu_time(t_stop)
  if(present(unit))then
     if(rankX==0.AND.rankY==0)write(unit,"(A20,F21.12)")"Time Gather :",t_stop-t_start
  endif
  !
  if(rankX==0.AND.rankY==0)call cpu_time(t_start)
  call Bcast_MPI(MPI_COMM_WORLD,M)
  if(rankX==0.AND.rankY==0)call cpu_time(t_stop)
  if(present(unit))then
     if(rankX==0.AND.rankY==0)write(unit,"(A20,F21.12)")"Time Bcast :",t_stop-t_start
  endif
  !
  return
end subroutine D_Gather_BLACS

subroutine Z_Gather_BLACS(Mloc,M,descM,unit)
  USE SF_MPI
  include 'mpif.h'
  complex(8),dimension(:,:),intent(in)    :: Mloc
  complex(8),dimension(:,:),intent(inout) :: M
  integer,dimension(9),intent(in)         :: descM
  integer,optional                        :: unit
  integer                                 :: i,j
  integer                                 :: myi,myj
  integer                                 :: Ni,Nj
  integer                                 :: Nbi,Nbj
  integer                                 :: Nrow,Ncol,Qrows
  integer                                 :: rankX,rankY
  integer                                 :: pNx,pNy
  integer                                 :: sendR,sendC
  real(8)                                 :: t_start,t_stop
  integer                                 :: context
  !
  context = descM(2)
  call blacs_gridinfo( context, pNx, pNy, rankX, rankY)
  !
  Ni  = size(M,1)
  Nj  = size(M,2)
  Nbi = descM(5)
  Nbj = descM(6)
  Qrows=descM(9) ; if(Qrows/=size(Mloc,1))stop "Gather_BLACS error: Qrows != size(Mloc,1)"
  !
  if(rankX==0.AND.rankY==0)call cpu_time(t_start)
  do j=1,Nj,Nbj
     Ncol = Nbj ; if(Nj-j<Nbj-1)Ncol=Nj-j+1
     do i=1,Ni,Nbi
        Nrow = Nbi ; if(Ni-i<Nbi-1)Nrow=Ni-i+1
        call infog2l(i,j,descM, pNx, pNy, rankX, rankY, myi, myj, SendR, SendC)
        if(rankX==SendR .AND. rankY==SendC)then
           call zgesd2d(context,Nrow,Ncol,Mloc(myi,myj),Qrows,0,0)
        endif
        if(rankX==0 .AND. rankY==0)then
           call zgerv2d(context,Nrow,Ncol,M(i,j),Ni,SendR,SendC)
        endif
     enddo
  enddo
  if(rankX==0.AND.rankY==0)call cpu_time(t_stop)
  if(present(unit))then
     if(rankX==0.AND.rankY==0)write(unit,"(A20,F21.12)")"Time Gather :",t_stop-t_start
  endif
  !
  if(rankX==0.AND.rankY==0)call cpu_time(t_start)
  call Bcast_MPI(MPI_COMM_WORLD,M)
  if(rankX==0.AND.rankY==0)call cpu_time(t_stop)
  if(present(unit))then
     if(rankX==0.AND.rankY==0)write(unit,"(A20,F21.12)")"Time Bcast :",t_stop-t_start
  endif
  !
  return
end subroutine Z_Gather_BLACS








