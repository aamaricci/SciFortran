!---------------------------------------------------------------------
!Purpose: use plain lanczos to get the groundstate energy
!---------------------------------------------------------------------
subroutine mpi_lanczos_eigh_d(MpiComm,MatVec,Ndim,Nitermax,Egs,Vect,iverbose,threshold,ncheck)
  integer                              :: MpiComm
  interface 
     subroutine MatVec(Nchunk,vin,vout)
       integer :: Nchunk
       real(8) :: vin(Nchunk)
       real(8) :: vout(Nchunk)
     end subroutine MatVec
  end interface
  integer                              :: Ndim
  integer                              :: Nitermax
  integer                              :: Nloc
  real(8)                              :: egs
  real(8),dimension(Ndim)              :: vect
  real(8),optional                     :: threshold
  integer,optional                     :: ncheck
  logical,optional                     :: iverbose
  !
  real(8),dimension(:),allocatable     :: vin,vout,vtmp
  integer                              :: iter,nlanc
  real(8),dimension(Nitermax+1)        :: alanc,blanc
  real(8),dimension(Nitermax,Nitermax) :: Z
  real(8),dimension(Nitermax)          :: diag,subdiag,esave
  real(8)                              :: a_,b_,norm,diff
  integer                              :: i,ierr
  !
  integer                              :: mpi_rank
  integer                              :: mpi_size
  logical                              :: mpi_master
  integer                              :: mpiQ
  integer                              :: mpiR
  !
  mpi_rank=get_rank_MPI(MpiComm)
  mpi_size=get_size_MPI(MpiComm)
  mpi_master=get_master_MPI(MpiComm)
  !
  mpiQ = Ndim/mpi_size
  mpiR = 0
  if(mpi_rank == mpi_size-1)mpiR=mod(Ndim,mpi_size)
  !
  Nloc = mpiQ+mpiR
  !
  if(present(iverbose))verb=iverbose
  if(present(threshold))threshold_=threshold
  if(present(ncheck))ncheck_=ncheck
  !
  norm=dot_product(vect,vect)
  if(norm==0d0)then
     call random_seed(size=nrandom)
     if(allocated(seed_random))deallocate(seed_random)
     allocate(seed_random(nrandom))
     seed_random=1234567
     call random_seed(put=seed_random)
     call random_number(vect)
     vect=vect/sqrt(dot_product(vect,vect))
     if(verb.AND.mpi_master)write(*,*)"MPI_LANCZOS_EIGH: random initial vector generated:"
  endif
  !
  !============= LANCZOS LOOP =====================
  !
  allocate(vin(Nloc),vout(Nloc))
  do i=1 + mpi_rank*mpiQ, (mpi_rank+1)*mpiQ + mpiR
     vin(i-mpi_rank*mpiQ)=vect(i)
  enddo
  vout= 0d0
  alanc=0d0
  blanc=0d0
  nlanc=0
  !
  lanc_loop: do iter=1,Nitermax
     !
     if(verb.AND.mpi_master)then
        print*,""
        write(*,*)"Lanczos iteration:",iter
     endif
     !
     call mpi_lanczos_iteration_d(MpiComm,MatVec,iter,vin,vout,a_,b_)
     if(abs(b_)<threshold_)exit lanc_loop
     !
     nlanc=nlanc+1
     !
     alanc(iter) = a_ ; blanc(iter+1) = b_
     !
     diag    = 0d0
     subdiag = 0.d0
     Z       = eye(Nlanc)
     diag(1:Nlanc)    = alanc(1:Nlanc)
     subdiag(2:Nlanc) = blanc(2:Nlanc)
     call tql2(Nlanc,diag,subdiag,Z,ierr)
     !
     if(verb.AND.mpi_master)then
        print *,'---> lowest eigenvalue  <---'
        write(*,*)"E_lowest    = ",diag(1)
        open(10,file="lanc_eigenvals.dat")
        do i=1,Nlanc
           write(10,*)i,diag(i)
        enddo
        close(10)
     endif
     !
     if(nlanc >= Ncheck_)then
        esave(nlanc-(Ncheck_-1))=diag(1)
        if(nlanc >= (Ncheck_+1))then
           diff=esave(Nlanc-(Ncheck_-1))-esave(Nlanc-(Ncheck_-1)-1)
           if(verb.AND.mpi_master)write(*,*)'test deltaE = ',diff
           if(abs(diff).le.threshold_)exit lanc_loop
        endif
     endif
     !
  enddo lanc_loop
  if(verb.AND.mpi_master)write(*,*)'Lanczos deltaE = ',diff
  if(nlanc==nitermax)print*,"LANCZOS_SIMPLE: reach Nitermax"
  !
  !============== END LANCZOS LOOP ======================
  !
  diag    = 0d0
  subdiag = 0.d0
  Z       = eye(Nlanc)
  diag(1:Nlanc)    = alanc(1:Nlanc)
  subdiag(2:Nlanc) = blanc(2:Nlanc)
  call tql2(Nlanc,diag,subdiag,Z,ierr)
  !
  !Get the Eigenvalues:
  egs = diag(1)
  !
  !Get the Eigenvector:
  do i=1 + mpi_rank*mpiQ, (mpi_rank+1)*mpiQ + mpiR
     vin(i-mpi_rank*mpiQ)=vect(i)
  enddo
  allocate(vtmp(Nloc))
  vtmp= 0d0
  vout= 0d0
  do iter=1,nlanc
     call mpi_lanczos_iteration_d(MpiComm,MatVec,iter,vin,vout,alanc(iter),blanc(iter))
     vtmp = vtmp + vin*Z(iter,1)
  end do
  vect= 0d0
  call AllReduce_MPI(MpiComm,vtmp,vect)
  norm=sqrt(dot_product(vect,vect))
  vect=vect/norm
end subroutine mpi_lanczos_eigh_d

subroutine mpi_lanczos_eigh_c(MpiComm,MatVec,Ndim,Nitermax,Egs,Vect,iverbose,threshold,ncheck)
  integer                              :: MpiComm
  interface 
     subroutine MatVec(Nchunk,vin,vout)
       integer    :: Nchunk
       complex(8) :: vin(Nchunk)
       complex(8) :: vout(Nchunk)
     end subroutine MatVec
  end interface
  integer                              :: Ndim
  integer                              :: Nitermax
  integer                              :: Nloc
  real(8)                              :: egs
  complex(8),dimension(Ndim)           :: vect
  real(8),optional                     :: threshold
  integer,optional                     :: ncheck
  logical,optional                     :: iverbose
  !
  complex(8),dimension(:),allocatable  :: vin,vout,vtmp
  integer                              :: iter,nlanc
  real(8),dimension(Nitermax+1)        :: alanc,blanc
  real(8),dimension(Nitermax,Nitermax) :: Z
  real(8),dimension(Nitermax)          :: diag,subdiag,esave
  real(8)                              :: a_,b_,norm,diff,ran(2)
  integer                              :: i,ierr
  !
  integer                              :: mpi_rank
  integer                              :: mpi_size
  logical                              :: mpi_master
  integer                              :: mpiQ
  integer                              :: mpiR
  !
  mpi_rank=get_rank_MPI(MpiComm)
  mpi_size=get_size_MPI(MpiComm)
  mpi_master=get_master_MPI(MpiComm)
  !
  mpiQ = Ndim/mpi_size
  mpiR = 0
  if(mpi_rank == mpi_size-1)mpiR=mod(Ndim,mpi_size)
  !
  Nloc = mpiQ+mpiR
  !
  if(present(iverbose))verb=iverbose
  if(present(threshold))threshold_=threshold
  if(present(ncheck))ncheck_=ncheck
  !
  norm=dot_product(vect,vect)
  if(norm==0d0)then
     call random_seed(size=nrandom)
     if(allocated(seed_random))deallocate(seed_random)
     allocate(seed_random(nrandom))
     seed_random=1234567
     call random_seed(put=seed_random)
     do i=1,Ndim
        call random_number(ran)
        vect(i)=dcmplx(ran(1),ran(2))
     enddo
     vect=vect/sqrt(dot_product(vect,vect))
     if(verb.AND.mpi_master)write(*,*)"MPI_LANCZOS_EIGH: random initial vector generated:"
  endif
  !
  !============= LANCZOS LOOP =====================
  !
  allocate(vin(Nloc),vout(Nloc))
  do i=1 + mpi_rank*mpiQ, (mpi_rank+1)*mpiQ + mpiR
     vin(i-mpi_rank*mpiQ)=vect(i)
  enddo
  vout= zero
  alanc=0d0
  blanc=0d0
  nlanc=0
  !
  lanc_loop: do iter=1,Nitermax
     !
     if(verb.AND.mpi_master)then
        print*,""
        write(*,*)"Lanczos iteration:",iter
     endif
     !
     call mpi_lanczos_iteration_c(MpiComm,MatVec,iter,vin,vout,a_,b_)
     if(abs(b_)<threshold_)exit lanc_loop
     !
     nlanc=nlanc+1
     !
     alanc(iter) = a_ ; blanc(iter+1) = b_
     !
     diag    = 0d0
     subdiag = 0d0
     Z       = eye(Nlanc)
     diag(1:Nlanc)    = alanc(1:Nlanc)
     subdiag(2:Nlanc) = blanc(2:Nlanc)
     call tql2(Nlanc,diag,subdiag,Z,ierr)
     !
     if(verb.AND.mpi_master)then
        print *,'---> lowest eigenvalue  <---'
        write(*,*)"E_lowest    = ",diag(1)
        open(10,file="lanc_eigenvals.dat")
        do i=1,Nlanc
           write(10,*)i,diag(i)
        enddo
        close(10)
     endif
     !
     if(nlanc >= Ncheck_)then
        esave(nlanc-(Ncheck_-1))=diag(1)
        if(nlanc >= (Ncheck_+1))then
           diff=esave(Nlanc-(Ncheck_-1))-esave(Nlanc-(Ncheck_-1)-1)
           if(verb.AND.mpi_master)write(*,*)'test deltaE = ',diff
           if(abs(diff).le.threshold_)exit lanc_loop
        endif
     endif
     !
  enddo lanc_loop
  if(verb.AND.mpi_master)write(*,*)'Lanczos deltaE = ',diff
  if(nlanc==nitermax)print*,"LANCZOS_SIMPLE: reach Nitermax"
  !
  !============== END LANCZOS LOOP ======================
  !
  diag    = 0d0
  subdiag = 0.d0
  Z       = eye(Nlanc)
  diag(1:Nlanc)    = alanc(1:Nlanc)
  subdiag(2:Nlanc) = blanc(2:Nlanc)
  call tql2(Nlanc,diag,subdiag,Z,ierr)
  !
  !Get the Eigenvalues:
  egs = diag(1)
  !
  !Get the Eigenvector:
  do i=1 + mpi_rank*mpiQ, (mpi_rank+1)*mpiQ + mpiR
     vin(i-mpi_rank*mpiQ)=vect(i)
  enddo
  allocate(vtmp(Nloc))
  vtmp= zero
  vout= zero
  do iter=1,nlanc
     call mpi_lanczos_iteration_c(MpiComm,MatVec,iter,vin,vout,alanc(iter),blanc(iter))
     vtmp = vtmp + vin*Z(iter,1)
  end do
  vect= 0d0
  call AllReduce_MPI(MpiComm,vtmp,vect)
  norm=sqrt(dot_product(vect,vect))
  vect=vect/norm
end subroutine mpi_lanczos_eigh_c






!---------------------------------------------------------------------
!Purpose: use simple Lanczos to tri-diagonalize a matrix H (defined 
! in the H*v function).
!---------------------------------------------------------------------
subroutine mpi_lanczos_tridiag_d(MpiComm,MatVec,vin,alanc,blanc,threshold)
  integer                                      :: MpiComm
  interface
     subroutine MatVec(Nchunk,vin,vout)
       integer                   :: Nchunk
       real(8),dimension(Nchunk) :: vin
       real(8),dimension(Nchunk) :: vout
     end subroutine MatVec
  end interface
  real(8),dimension(:),intent(inout)           :: vin
  real(8),dimension(:),allocatable             :: vout,vtmp
  real(8),dimension(:),intent(inout)           :: alanc
  real(8),dimension(size(alanc)),intent(inout) :: blanc
  integer                                      :: Nitermax,Nloc,i
  integer                                      :: iter
  real(8)                                      :: a_,b_
  real(8),optional                             :: threshold
  !
  integer                                      :: mpi_rank
  integer                                      :: mpi_size
  logical                                      :: mpi_master
  integer                                      :: mpiQ
  integer                                      :: mpiR
  !
  mpi_rank=get_rank_MPI(MpiComm)
  mpi_size=get_size_MPI(MpiComm)
  mpi_master=get_master_MPI(MpiComm)
  !
  mpiQ = size(vin)/mpi_size
  mpiR = 0
  if(mpi_rank == mpi_size-1)mpiR=mod(size(vin),mpi_size)
  !
  Nloc = mpiQ+mpiR
  !
  if(present(threshold))threshold_=threshold
  !
  allocate(vtmp(Nloc),vout(Nloc))
  do i=1 + mpi_rank*mpiQ, (mpi_rank+1)*mpiQ + mpiR
     vtmp(i-mpi_rank*mpiQ)=vin(i)
  enddo
  Nitermax = size(alanc)
  a_=0d0
  b_=0d0
  vout=0d0
  do iter=1,Nitermax
     call mpi_lanczos_iteration_d(MpiComm,MatVec,iter,vtmp,vout,a_,b_)
     alanc(iter)=a_
     if(abs(b_)<threshold_)exit
     if(iter<nitermax)blanc(iter+1)=b_
  enddo
  if(iter==nitermax.AND.mpi_master)write(*,"(A)")"MPI_LANCZOS_TRIDIAG_D: reach Nitermax"
  deallocate(vtmp,vout)
end subroutine mpi_lanczos_tridiag_d

subroutine mpi_lanczos_tridiag_c(MpiComm,MatVec,vin,alanc,blanc,threshold)
  integer                                      :: MpiComm
  interface
     subroutine MatVec(Nchunk,vin,vout)
       integer                      :: Nchunk
       complex(8),dimension(Nchunk) :: vin
       complex(8),dimension(Nchunk) :: vout
     end subroutine MatVec
  end interface
  complex(8),dimension(:),intent(inout)        :: vin
  complex(8),dimension(:),allocatable          :: vout,vtmp
  real(8),dimension(:),intent(inout)           :: alanc
  real(8),dimension(size(alanc)),intent(inout) :: blanc
  integer                                      :: Nitermax,Nloc,i
  integer                                      :: iter
  real(8)                                      :: a_,b_
  real(8),optional                             :: threshold
  !
  integer                                      :: mpi_rank
  integer                                      :: mpi_size
  logical                                      :: mpi_master
  integer                                      :: mpiQ
  integer                                      :: mpiR
  !
  mpi_rank=get_rank_MPI(MpiComm)
  mpi_size=get_size_MPI(MpiComm)
  mpi_master=get_master_MPI(MpiComm)
  !
  mpiQ = size(vin)/mpi_size
  mpiR = 0
  if(mpi_rank == mpi_size-1)mpiR=mod(size(vin),mpi_size)
  !
  Nloc = mpiQ+mpiR
  !
  if(present(threshold))threshold_=threshold
  !
  allocate(vtmp(Nloc),vout(Nloc))
  do i=1 + mpi_rank*mpiQ, (mpi_rank+1)*mpiQ + mpiR
     vtmp(i-mpi_rank*mpiQ)=vin(i)
  enddo
  Nitermax = size(alanc)
  a_=0d0
  b_=0d0
  vout=zero
  do iter=1,Nitermax
     call mpi_lanczos_iteration_c(MpiComm,MatVec,iter,vtmp,vout,a_,b_)
     alanc(iter)=a_
     if(abs(b_)<threshold_)exit
     if(iter<nitermax)blanc(iter+1)=b_
  enddo
  if(iter==nitermax.AND.mpi_master)write(*,"(A)")"MPI_LANCZOS_TRIDIAG_D: reach Nitermax"
  deallocate(vtmp,vout)
end subroutine mpi_lanczos_tridiag_c




!---------------------------------------------------------------------
!Purpose: plain homebrew lanczos iteration (no orthogonalization)
!note: the a,b variables are real, even in the complex matrix case
!to understand why check out the Gollub-Van Loan textbook.
!a it is easy: hermiticity->diag\in\RRR
!b: is fixed by requiring |b|^2 = <v,v> thus you can only fix the 
!the absolute value. A lemma shows that the phase can be chosen 
!identically zero
!MPI VERSION
!---------------------------------------------------------------------
subroutine mpi_lanczos_iteration_d(MpiComm,MatVec,iter,vin,vout,a,b)
  integer                                    :: MpiComm
  interface
     subroutine MatVec(Nchunk,vin,vout)
       integer                   :: Nchunk
       real(8),dimension(Nchunk) :: vin
       real(8),dimension(Nchunk) :: vout
     end subroutine MatVec
  end interface
  real(8),dimension(:),intent(inout)         :: vin
  real(8),dimension(size(vin)),intent(inout) :: vout
  real(8),dimension(size(vin))               :: tmp
  real(8),intent(inout)                      :: a,b
  real(8)                                    :: atmp,btmp
  integer                                    :: iter,ndim,nloc
  real(8)                                    :: norm,norm_tmp
  !
  logical                                    :: mpi_master
  !
  nloc=size(vin)
  !
  mpi_master=get_master_MPI(MpiComm)
  !
  if(iter==1)then
     norm = 0d0
     norm_tmp=dot_product(vin,vin)
     call AllReduce_MPI(MpiComm,norm_tmp,norm)
     if(mpi_master.AND.norm==0d0)stop "MPI_LANCZOS_ITERATION_D: norm = 0!!"
     vin=vin/sqrt(norm)
     b=0d0
  end if
  !
  call MatVec(nloc,vin,tmp)
  tmp   = tmp-b*vout
  atmp  = dot_product(vin,tmp)
  a     = 0d0
  call AllReduce_MPI(MpiComm,atmp,a)
  tmp   = tmp-a*vin
  btmp  = dot_product(tmp,tmp) !sqrt(dot_product(tmp,tmp))
  b     = 0d0
  call AllReduce_MPI(MpiComm,btmp,b)
  b     = sqrt(b)
  vout  = vin
  vin   = tmp/b
end subroutine mpi_lanczos_iteration_d

subroutine mpi_lanczos_iteration_c(MpiComm,MatVec,iter,vin,vout,a,b)
  integer                                       :: MpiComm
  interface
     subroutine MatVec(Nchunk,vin,vout)
       integer                      :: Nchunk
       complex(8),dimension(Nchunk) :: vin
       complex(8),dimension(Nchunk) :: vout
     end subroutine MatVec
  end interface
  complex(8),dimension(:),intent(inout)         :: vin
  complex(8),dimension(size(vin)),intent(inout) :: vout
  complex(8),dimension(size(vin))               :: tmp
  real(8),intent(inout)                         :: a,b
  real(8)                                       :: atmp,btmp
  integer                                       :: iter,ndim,nloc
  real(8)                                       :: norm,norm_tmp
  !
  logical                                       :: mpi_master
  !
  nloc=size(vin)
  !
  mpi_master=get_master_MPI(MpiComm)
  !
  if(iter==1)then
     norm = 0d0
     norm_tmp=dot_product(vin,vin)
     call AllReduce_MPI(MpiComm,norm_tmp,norm)
     if(mpi_master.AND.norm==0d0)stop "MPI_LANCZOS_ITERATION_D: norm = 0!!"
     vin=vin/sqrt(norm)
     b=0d0
  end if
  !
  call MatVec(nloc,vin,tmp)
  tmp   = tmp-b*vout
  atmp  = dot_product(vin,tmp)
  a     = 0d0
  call AllReduce_MPI(MpiComm,atmp,a)
  tmp   = tmp-a*vin
  btmp  = dot_product(tmp,tmp) !sqrt(dot_product(tmp,tmp))
  b     = 0d0
  call AllReduce_MPI(MpiComm,btmp,b)
  b     = sqrt(b)
  vout  = vin
  vin   = tmp/b
end subroutine mpi_lanczos_iteration_c
