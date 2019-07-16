subroutine p_d_matmul(A,B,C,Nblock,alfa,beta)
  real(8),dimension(:,:),intent(inout) :: A ![N,K]
  real(8),dimension(:,:),intent(inout) :: B ![K,M]
  real(8),dimension(:,:),intent(inout) :: C ![N,M]
  integer                              :: Nblock
  real(8),optional                     :: alfa,beta
  real(8)                              :: alfa_,beta_
  integer                              :: Nb
  integer                              :: N,K,M
  integer                              :: Nrows,Kcols,Krows,Mcols
  !
  integer                              :: i,j,lda,info
  integer,external                     :: numroc,indxG2L,indxG2P,indxL2G
  !
  integer,dimension(9)                 :: descAloc,descBloc,descCloc,descC
  real(8),dimension(:,:),allocatable   :: Aloc,Bloc,Cloc
  integer                              :: rankX,rankY
  integer                              :: sendR,sendC,Nipiv
  integer                              :: Nrow,Ncol
  integer                              :: myi,myj,unit
  real(8)                              :: t_stop,t_start
  logical                              :: master
  !
  ! C = alfa*A*B + beta*C
  !
  alfa_     = 1d0 ;if(present(alfa))alfa_=alfa
  beta_     = 0d0 ;if(present(beta))beta_=beta
  !
  !
  N = size(A,1)
  K = size(A,2) !==size(B,1)
  M = size(B,2)
  if(any(shape(B)/=[K,M]))stop "p_d_matmul error: B has illegal shape"
  if(any(shape(C)/=[N,M]))stop "p_d_matmul error: C has illegal shape"
  !
  !INIT SCALAPACK TREATMENT:
  !
  !< Get coordinate of the processes
  call blacs_gridinfo(p_context, p_Nx, p_Ny, rankX, rankY)
  master = (rankX==0).AND.(rankY==0)
  if(rankX<0.AND.rankY<0)goto 100
  !
  Nb = Nblock
  !
  Nrows = numroc(N, Nb, rankX, 0, p_Nx)
  Kcols = numroc(K, Nb, rankY, 0, p_Ny)
  Krows = numroc(K, Nb, rankX, 0, p_Nx)
  Mcols = numroc(M, Nb, rankY, 0, p_Ny)
  !
  if(master)then
     unit = free_unit()
     open(unit,file="p_dgemm.info")
     write(unit,"(A20,I8,A5,I8)")"Grid=",p_Nx,"x",p_Ny
     write(unit,"(A20,I8,A5,I8)")"Nrows x Kcols=",Nrows,"x",Kcols
     write(unit,"(A20,I8,A5,I8)")"Krows x Mcols=",Krows,"x",Mcols
  endif
  !
  !< allocate local distributed A
  !
  call descinit( descAloc, N, K, Nb, Nb, 0, 0, p_context, Nrows, info )
  call descinit( descBloc, K, M, Nb, Nb, 0, 0, p_context, Krows, info )
  call descinit( descCloc, N, M, Nb, Nb, 0, 0, p_context, Nrows, info )
  !
  !< Distribute A,B; allocate C
  allocate(Aloc(Nrows,Kcols))
  allocate(Bloc(Krows,Mcols))
  allocate(Cloc(Nrows,Mcols))
  if(master)call cpu_time(t_start)
  do myi=1,size(Aloc,1)
     i  = indxL2G(myi,Nblock,rankX,0,p_Nx)
     do myj=1,size(Aloc,2)
        j  = indxL2G(myj,Nblock,rankY,0,p_Ny)
        Aloc(myi,myj) = A(i,j)
     enddo
  enddo
  do myi=1,size(Bloc,1)
     i  = indxL2G(myi,Nblock,rankX,0,p_Nx)
     do myj=1,size(Bloc,2)
        j  = indxL2G(myj,Nblock,rankY,0,p_Ny)
        Bloc(myi,myj) = B(i,j)
     enddo
  enddo
  !
  if(master)call cpu_time(t_stop)
  if(master)write(unit,"(A20,F21.12)")"Time Distribute A,B:",t_stop-t_start
  !
  if(master)call cpu_time(t_start)
  call PDGEMM( 'No transpose', 'No transpose', N, M, K, &
       alfa_, &
       Aloc, 1, 1, descAloc, &
       Bloc, 1, 1, descBloc, &
       beta_, &
       Cloc, 1, 1, descCloc )
  !
  if(master)call cpu_time(t_stop)
  if(master)write(unit,"(A20,F21.12)")"Time matmul C_loc:",t_stop-t_start
  !
  C=0d0
  call descinit( descC, N, M, Nb, Nb, 0, 0, p_context, Nrows, info )
  if(master)call cpu_time(t_start)
  do i=1,N,Nb
     Nrow = Nb ; if(N-i<Nb-1)Nrow=N-i+1!;if(Nrow==0)Nrow=1
     do j=1,M,Nb
        Ncol = Nb ; if(M-j<Nb-1)Ncol=M-j+1!;if(Ncol==0)Ncol=1
        call infog2l(i,j,descC, p_Nx, p_Ny, rankX, rankY, myi, myj, SendR, SendC)
        if(rankX==SendR .AND. rankY==SendC)then
           call dgesd2d(p_context,Nrow,Ncol,Cloc(myi,myj),Nrows,0,0)
        endif
        if(master)then
           call dgerv2d(p_context,Nrow,Ncol,C(i,j),N,SendR,SendC)
        endif
     enddo
  enddo
  if(master)call cpu_time(t_stop)
  if(master)write(unit,"(A20,F21.12)")"Time gather C:",t_stop-t_start
  !
  if(master)close(unit)
100 continue
  return
  !
end subroutine p_d_matmul

subroutine p_z_matmul(A,B,C,Nblock,alfa,beta)
  complex(8),dimension(:,:),intent(in)                    :: A ![N,K]
  complex(8),dimension(:,:),intent(in)                    :: B ![K,M]
  complex(8),dimension(size(A,1),size(B,2)),intent(inout) :: C ![N,M]
  integer                                                 :: Nblock
  real(8),optional                                        :: alfa,beta
  real(8)                                                 :: alfa_,beta_
  integer                                                 :: Nb
  integer                                                 :: N,K,M
  integer                                                 :: Nrows,Kcols,Krows,Mcols
  !
  integer                                                 :: i,j,lda,info
  integer,external                                        :: numroc,indxG2L,indxG2P,indxL2G
  !
  integer,dimension(9)                                    :: descAloc,descBloc,descCloc,descC
  complex(8),dimension(:,:),allocatable                   :: Aloc,Bloc,Cloc
  integer                                                 :: rankX,rankY
  integer                                                 :: sendR,sendC,Nipiv
  integer                                                 :: Nrow,Ncol
  integer                                                 :: myi,myj,unit
  real(8)                                                 :: t_stop,t_start
  logical                                                 :: master
  !
  ! C = alfa*A*B + beta*C
  !
  alfa_     = dcmplx(1d0,0d0) ;if(present(alfa))alfa_=alfa
  beta_     = dcmplx(0d0,0d0) ;if(present(beta))beta_=beta
  !
  !
  N = size(A,1)
  K = size(A,2) !==size(B,1)
  M = size(B,2)
  if(any(shape(B)/=[K,M]))stop "p_d_matmul error: B has illegal shape"
  if(any(shape(C)/=[N,M]))stop "p_d_matmul error: C has illegal shape"
  !
  !INIT SCALAPACK TREATMENT:
  !
  !< Get coordinate of the processes
  call blacs_gridinfo( p_context, p_Nx, p_Ny, rankX, rankY)
  !
  master = (rankX==0).AND.(rankY==0)
  if(rankX<0.AND.rankY<0)goto 1010
  !
  Nb = Nblock
  !
  Nrows = numroc(N, Nb, rankX, 0, p_Nx)
  Kcols = numroc(K, Nb, rankY, 0, p_Ny)
  Krows = numroc(K, Nb, rankX, 0, p_Nx)
  Mcols = numroc(M, Nb, rankY, 0, p_Ny)
  !
  if(master)then
     unit = free_unit()
     open(unit,file="p_dgemm.info")
     write(unit,"(A20,I8,A5,I8)")"Grid=",p_Nx,"x",p_Ny
     write(unit,"(A20,I8,A5,I8)")"Nrows x Kcols=",Nrows,"x",Kcols
     write(unit,"(A20,I8,A5,I8)")"Krows x Mcols=",Krows,"x",Mcols
  endif
  !
  !< allocate local distributed A
  !
  call descinit( descAloc, N, K, Nb, Nb, 0, 0, p_context, Nrows, info )
  call descinit( descBloc, K, M, Nb, Nb, 0, 0, p_context, Krows, info )
  call descinit( descCloc, N, M, Nb, Nb, 0, 0, p_context, Nrows, info )
  !
  !< Distribute A,B; allocate C
  allocate(Aloc(Nrows,Kcols));Aloc=zero
  allocate(Bloc(Krows,Mcols));Bloc=zero
  allocate(Cloc(Nrows,Mcols));Cloc=zero
  if(master)call cpu_time(t_start)
  do myi=1,size(Aloc,1)
     i  = indxL2G(myi,Nblock,rankX,0,p_Nx)
     do myj=1,size(Aloc,2)
        j  = indxL2G(myj,Nblock,rankY,0,p_Ny)
        Aloc(myi,myj) = A(i,j)
     enddo
  enddo
  do myi=1,size(Bloc,1)
     i  = indxL2G(myi,Nblock,rankX,0,p_Nx)
     do myj=1,size(Bloc,2)
        j  = indxL2G(myj,Nblock,rankY,0,p_Ny)
        Bloc(myi,myj) = B(i,j)
     enddo
  enddo
  !
  if(master)call cpu_time(t_stop)
  if(master)write(unit,"(A20,F21.12)")"Time Distribute A,B:",t_stop-t_start
  !
  if(master)call cpu_time(t_start)
  call PZGEMM( 'No transpose', 'No transpose', N, M, K, &
       alfa_, &
       Aloc, 1, 1, descAloc, &
       Bloc, 1, 1, descBloc, &
       beta_, &
       Cloc, 1, 1, descCloc )
  !
  if(master)call cpu_time(t_stop)
  if(master)write(unit,"(A20,F21.12)")"Time matmul C_loc:",t_stop-t_start
  !
  C=dcmplx(0d0,0d0)
  call descinit( descC, N, M, Nb, Nb, 0, 0, p_context, N, info )
  if(master)call cpu_time(t_start)
  do i=1,N,Nb
     Nrow = Nb ; if(N-i<Nb-1)Nrow=N-i+1!;if(Nrow==0)Nrow=1
     do j=1,M,Nb
        Ncol = Nb ; if(M-j<Nb-1)Ncol=M-j+1!;if(Ncol==0)Ncol=1
        call infog2l(i,j,descC, p_Nx, p_Ny, rankX, rankY, myi, myj, SendR, SendC)
        if(rankX==SendR .AND. rankY==SendC)then
           call zgesd2d(p_context,Nrow,Ncol,Cloc(myi,myj),Nrows,0,0)
        endif
        if(rankX==0 .AND. rankY==0)then
           call zgerv2d(p_context,Nrow,Ncol,C(i,j),N,SendR,SendC)
        endif
     enddo
  enddo
  if(master)call cpu_time(t_stop)
  if(master)write(unit,"(A20,F21.12)")"Time gather C:",t_stop-t_start
  !
  if(master)close(unit)
1010 continue
  return
  !
end subroutine p_z_matmul






!###### OVERLOAD PARALLEL MATMUL OPERATOR --> .px. #########


function p_d_matmul_f(A,B) result(C)
  real(8),dimension(:,:),intent(in)      :: A ![N,K]
  real(8),dimension(:,:),intent(in)      :: B ![K,M]
  real(8),dimension(size(A,1),size(B,2)) :: C ![N,M]
  real(8)                                :: alfa_,beta_
  integer                                :: Nb
  integer                                :: N,K,M
  integer                                :: Nrows,Kcols,Krows,Mcols
  !
  integer                                :: i,j,lda,info
  integer,external                       :: numroc,indxG2L,indxG2P,indxL2G
  !
  integer,dimension(9)                   :: descAloc,descBloc,descCloc,descC
  real(8),dimension(:,:),allocatable     :: Aloc,Bloc,Cloc
  integer                                :: rankX,rankY
  integer                                :: sendR,sendC,Nipiv
  integer                                :: Nrow,Ncol
  integer                                :: myi,myj,unit
  real(8)                                :: t_stop,t_start
  logical                                :: master
  !
  ! C = alfa*A*B + beta*C
  !
  alfa_     = 1d0 
  beta_     = 0d0 
  !
  N = size(A,1)
  K = size(A,2) !==size(B,1)
  M = size(B,2)
  if(any(shape(B)/=[K,M]))stop "p_d_matmul error: B has illegal shape"
  if(any(shape(C)/=[N,M]))stop "p_d_matmul error: C has illegal shape"
  !
  !INIT SCALAPACK TREATMENT:
  !< Get coordinate of the processes
  call blacs_gridinfo( p_context, p_Nx, p_Ny, rankX, rankY)
  master = (rankX==0).AND.(rankY==0)
  if(rankX<0.AND.rankY<0)goto 201
  !
  I  = min(N,K,M)
  Nb = min(shift_dw(I)/2,64)
  !
  Nrows = numroc(N, Nb, rankX, 0, p_Nx)
  Kcols = numroc(K, Nb, rankY, 0, p_Ny)
  Krows = numroc(K, Nb, rankX, 0, p_Nx)
  Mcols = numroc(M, Nb, rankY, 0, p_Ny)
  !
  !< allocate local distributed A
  !
  call descinit( descAloc, N, K, Nb, Nb, 0, 0, p_context, Nrows, info )
  call descinit( descBloc, K, M, Nb, Nb, 0, 0, p_context, Krows, info )
  call descinit( descCloc, N, M, Nb, Nb, 0, 0, p_context, Nrows, info )
  !
  !< Distribute A,B; allocate C
  allocate(Aloc(Nrows,Kcols));Aloc=0d0
  allocate(Bloc(Krows,Mcols));Bloc=0d0
  allocate(Cloc(Nrows,Mcols));Cloc=0d0
  do myi=1,size(Aloc,1)
     i  = indxL2G(myi,Nb,rankX,0,p_Nx)
     do myj=1,size(Aloc,2)
        j  = indxL2G(myj,Nb,rankY,0,p_Ny)
        Aloc(myi,myj) = A(i,j)
     enddo
  enddo
  do myi=1,size(Bloc,1)
     i  = indxL2G(myi,Nb,rankX,0,p_Nx)
     do myj=1,size(Bloc,2)
        j  = indxL2G(myj,Nb,rankY,0,p_Ny)
        Bloc(myi,myj) = B(i,j)
     enddo
  enddo
  !
  call PDGEMM( 'No transpose', 'No transpose', N, M, K, &
       alfa_, &
       Aloc, 1, 1, descAloc, &
       Bloc, 1, 1, descBloc, &
       beta_, &
       Cloc, 1, 1, descCloc )
  !
  C=0d0
  call descinit( descC, N, M, Nb, Nb, 0, 0, p_context, N, info )
  do i=1,N,Nb
     Nrow = Nb ; if(N-i<Nb-1)Nrow=N-i+1!;if(Nrow==0)Nrow=1
     do j=1,M,Nb
        Ncol = Nb ; if(M-j<Nb-1)Ncol=M-j+1!;if(Ncol==0)Ncol=1
        call infog2l(i,j,descC, p_Nx, p_Ny, rankX, rankY, myi, myj, SendR, SendC)
        if(rankX==SendR .AND. rankY==SendC)then
           call dgesd2d(p_context,Nrow,Ncol,Cloc(myi,myj),Nrows,0,0)
        endif
        if(rankX==0 .AND. rankY==0)then
           call dgerv2d(p_context,Nrow,Ncol,C(i,j),N,SendR,SendC)
        endif
     enddo
  enddo
201 continue
  return
  !
end function p_d_matmul_f

function p_z_matmul_f(A,B) result(C)
  complex(8),dimension(:,:),intent(in)      :: A ![N,K]
  complex(8),dimension(:,:),intent(in)      :: B ![K,M]
  complex(8),dimension(size(A,1),size(B,2)) :: C ![N,M]
  real(8)                                   :: alfa_,beta_
  integer                                   :: Nb
  integer                                   :: N,K,M
  integer                                   :: Nrows,Kcols,Krows,Mcols
  !
  integer                                   :: i,j,lda,info
  integer,external                          :: numroc,indxG2L,indxG2P,indxL2G
  !
  integer,dimension(9)                      :: descAloc,descBloc,descCloc,descC
  complex(8),dimension(:,:),allocatable     :: Aloc,Bloc,Cloc
  integer                                   :: rankX,rankY
  integer                                   :: sendR,sendC,Nipiv
  integer                                   :: Nrow,Ncol
  integer                                   :: myi,myj,unit
  real(8)                                   :: t_stop,t_start
  logical                                   :: master
  !
  ! C = alfa*A*B + beta*C
  !
  alfa_     = dcmplx(1d0,0d0)
  beta_     = dcmplx(0d0,0d0)
  !
  !
  N = size(A,1)
  K = size(A,2) !==size(B,1)
  M = size(B,2)
  if(any(shape(B)/=[K,M]))stop "p_d_matmul error: B has illegal shape"
  if(any(shape(C)/=[N,M]))stop "p_d_matmul error: C has illegal shape"
  !
  !INIT SCALAPACK TREATMENT:
  !
  !< Get coordinate of the processes
  call blacs_gridinfo( p_context, p_Nx, p_Ny, rankX, rankY)
  master = (rankX==0).AND.(rankY==0)
  if(rankX<0.AND.rankY<0)goto 101
  !
  I  = min(N,K,M)
  Nb = min(shift_dw(I)/2,64)
  !
  Nrows = numroc(N, Nb, rankX, 0, p_Nx)
  Kcols = numroc(K, Nb, rankY, 0, p_Ny)
  Krows = numroc(K, Nb, rankX, 0, p_Nx)
  Mcols = numroc(M, Nb, rankY, 0, p_Ny)
  !
  !< allocate local distributed A
  !
  call descinit( descAloc, N, K, Nb, Nb, 0, 0, p_context, Nrows, info )
  call descinit( descBloc, K, M, Nb, Nb, 0, 0, p_context, Krows, info )
  call descinit( descCloc, N, M, Nb, Nb, 0, 0, p_context, Nrows, info )
  !
  !< Distribute A,B; allocate C
  allocate(Aloc(Nrows,Kcols));Aloc=zero
  allocate(Bloc(Krows,Mcols));Bloc=zero
  allocate(Cloc(Nrows,Mcols));Cloc=zero
  do myi=1,size(Aloc,1)
     i  = indxL2G(myi,Nb,rankX,0,p_Nx)
     do myj=1,size(Aloc,2)
        j  = indxL2G(myj,Nb,rankY,0,p_Ny)
        Aloc(myi,myj) = A(i,j)
     enddo
  enddo
  do myi=1,size(Bloc,1)
     i  = indxL2G(myi,Nb,rankX,0,p_Nx)
     do myj=1,size(Bloc,2)
        j  = indxL2G(myj,Nb,rankY,0,p_Ny)
        Bloc(myi,myj) = B(i,j)
     enddo
  enddo
  !
  !
  call PZGEMM( 'No transpose', 'No transpose', N, M, K, &
       alfa_, &
       Aloc, 1, 1, descAloc, &
       Bloc, 1, 1, descBloc, &
       beta_, &
       Cloc, 1, 1, descCloc )
  !
  C=dcmplx(0d0,0d0)
  call descinit( descC, N, M, Nb, Nb, 0, 0, p_context, N, info)
  do i=1,N,Nb
     Nrow = Nb ; if(N-i<Nb-1)Nrow=N-i+1!;if(Nrow==0)Nrow=1
     do j=1,M,Nb
        Ncol = Nb ; if(M-j<Nb-1)Ncol=M-j+1!;if(Ncol==0)Ncol=1
        call infog2l(i,j,descC, p_Nx, p_Ny, rankX, rankY, myi, myj, SendR, SendC)
        if(rankX==SendR .AND. rankY==SendC)then
           call zgesd2d(p_context,Nrow,Ncol,Cloc(myi,myj),Nrows,0,0)
        endif
        if(rankX==0 .AND. rankY==0)then
           call zgerv2d(p_context,Nrow,Ncol,C(i,j),N,SendR,SendC)
        endif
     enddo
  enddo
101 continue
  return
  !
end function p_z_matmul_f
