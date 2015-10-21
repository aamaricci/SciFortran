subroutine splotM_II(pname,X,Y1,Y2,Y3,Y4,append)
  integer                                  :: i,j,Np,Ny1,Ny2
  character(len=*)                         :: pname
  integer,dimension(:,:)                   :: Y1
  integer,dimension(size(Y1,1),size(Y1,2)),optional :: Y2,Y3,Y4
  integer,dimension(size(Y1,2))            :: X
  logical,optional                    :: append
  logical                             :: check,rw
  rw=.false.;if(present(append))rw=append
  Ny1=size(Y1,1) ; Ny2=size(Y1,2)
  open(719,file=adjustl(trim(pname)))
  if(present(Y4))then
     do i=1,Ny1
        do j=1,Ny2
           write(719,'(5(I15))')X(j),Y1(i,j),Y2(i,j),Y3(i,j),Y4(i,j)
        enddo
        write(719,*)""
     enddo
  elseif(present(Y3))then
     do i=1,Ny1
        do j=1,Ny2
           write(719,'(4(I15))')X(j),Y1(i,j),Y2(i,j),Y3(i,j)
        enddo
        write(719,*)""
     enddo
  elseif(present(Y2))then
     do i=1,Ny1
        do j=1,Ny2
           write(719,'(3(I15))')X(j),Y1(i,j),Y2(i,j)
        enddo
        write(719,*)""
     enddo
  else
     do i=1,Ny1
        do j=1,Ny2
           write(719,*)X(j),Y1(i,j)
        enddo
        write(719,*)""
     enddo
  endif
  close(719)
end subroutine splotM_II

subroutine splotM_IR(pname,X,Y1,Y2,Y3,Y4,append)
  integer                                :: i,j,Ny1,Ny2
  character(len=*)                       :: pname
  real(8),dimension(:,:)                   :: Y1
  real(8),dimension(size(Y1,1),size(Y1,2)),optional :: Y2,Y3,Y4
  integer,dimension(size(Y1,2))            :: X
  logical,optional                    :: append
  logical                             :: check,rw
  rw=.false.;if(present(append))rw=append
  Ny1=size(Y1,1) ; Ny2=size(Y1,2)
  open(719,file=adjustl(trim(pname)))
  if(present(Y4))then
     do i=1,Ny1
        do j=1,Ny2
           write(719,"(I15,4(F21.12))")X(j),Y1(i,j),Y2(i,j),Y3(i,j),Y4(i,j)
        enddo
        write(719,*)""
     enddo
  elseif(present(Y3))then
     do i=1,Ny1
        do j=1,Ny2
           write(719,"(I15,3(F21.12))")X(j),Y1(i,j),Y2(i,j),Y3(i,j)
        enddo
        write(719,*)""
     enddo
  elseif(present(Y2))then
     do i=1,Ny1
        do j=1,Ny2
           write(719,"(I15,2(F21.12))")X(j),Y1(i,j),Y2(i,j)
        enddo
        write(719,*)""
     enddo
  else
     do i=1,Ny1
        do j=1,Ny2
           write(719,*)X(j),Y1(i,j)
        enddo
        write(719,*)""
     enddo
  endif
  close(719)
end subroutine splotM_IR

subroutine splotM_IC(pname,X,Y1,Y2,append)
  integer                                     :: i,j,Ny1,Ny2
  character(len=*)                            :: pname
  complex(8),dimension(:,:)                   :: Y1
  complex(8),dimension(size(Y1,1),size(Y1,2)),optional :: Y2
  integer,dimension(size(Y1,2))               :: X
  logical,optional                    :: append
  logical                             :: check,rw
  rw=.false.;if(present(append))rw=append
  Ny1=size(Y1,1) ; Ny2=size(Y1,2)
  open(719,file=adjustl(trim(pname)))
  if(present(Y2))then
     do i=1,Ny1
        do j=1,Ny2
           write(719,"(I15,4(F21.12))")X(j),dimag(Y1(i,j)),dreal(Y1(i,j)),&
                dimag(Y2(i,j)),dreal(Y2(i,j))
        enddo
        write(719,*)""
     enddo
  else
     do i=1,Ny1
        do j=1,Ny2
           write(719,*)X(j),dimag(Y1(i,j)),dreal(Y1(i,j))
        enddo
        write(719,*)""
     enddo
  endif
  close(719)
end subroutine splotM_IC


!----------------------------

subroutine splotM_RI(pname,X,Y1,Y2,Y3,Y4,append)
  integer                                :: i,j,Ny1,Ny2
  character(len=*)                       :: pname
  integer,dimension(:,:)                   :: Y1
  integer,dimension(size(Y1,1),size(Y1,2)),optional :: Y2,Y3,Y4
  real(8),dimension(size(Y1,2))            :: X
  logical,optional                    :: append
  logical                             :: check,rw
  rw=.false.;if(present(append))rw=append
  Ny1=size(Y1,1) ; Ny2=size(Y1,2)
  open(719,file=adjustl(trim(pname)))
  if(present(Y4))then
     do i=1,Ny1
        do j=1,Ny2
           write(719,"(F21.12,4(I15))")X(j),Y1(i,j),Y2(i,j),Y3(i,j),Y4(i,j)
        enddo
        write(719,*)""
     enddo
  elseif(present(Y3))then
     do i=1,Ny1
        do j=1,Ny2
           write(719,"(F21.12,3(I15))")X(j),Y1(i,j),Y2(i,j),Y3(i,j)
        enddo
        write(719,*)""
     enddo
  elseif(present(Y2))then
     do i=1,Ny1
        do j=1,Ny2
           write(719,*)X(j),Y1(i,j),Y2(i,j)
        enddo
        write(719,*)""
     enddo
  else
     do i=1,Ny1
        do j=1,Ny2
           write(719,*)X(j),Y1(i,j)
        enddo
        write(719,*)""
     enddo
  endif
  close(719)
end subroutine splotM_RI

subroutine splotM_RR(pname,X,Y1,Y2,Y3,Y4,append)
  integer                                :: i,j,Ny1,Ny2
  character(len=*)                       :: pname
  real(8),dimension(:,:)                   :: Y1
  real(8),dimension(size(Y1,1),size(Y1,2)),optional :: Y2,Y3,Y4
  real(8),dimension(size(Y1,2))            :: X
  logical,optional                    :: append
  logical                             :: check,rw
  rw=.false.;if(present(append))rw=append
  Ny1=size(Y1,1) ; Ny2=size(Y1,2)
  open(719,file=adjustl(trim(pname)))
  if(present(Y4))then
     do i=1,Ny1
        do j=1,Ny2
           write(719,"(5(F21.12))")X(j),Y1(i,j),Y2(i,j),Y3(i,j),Y4(i,j)
        enddo
        write(719,*)""
     enddo
  elseif(present(Y3))then
     do i=1,Ny1
        do j=1,Ny2
           write(719,"(4(F21.12))")X(j),Y1(i,j),Y2(i,j),Y3(i,j)
        enddo
        write(719,*)""
     enddo
  elseif(present(Y2))then
     do i=1,Ny1
        do j=1,Ny2
           write(719,*)X(j),Y1(i,j),Y2(i,j)
        enddo
        write(719,*)""
     enddo
  else
     do i=1,Ny1
        do j=1,Ny2
           write(719,*)X(j),Y1(i,j)
        enddo
        write(719,*)""
     enddo
  endif
  close(719)
end subroutine splotM_RR

subroutine splotM_RC(pname,X,Y1,Y2,append)
  integer                                     :: i,j,Ny1,Ny2
  character(len=*)                            :: pname
  complex(8),dimension(:,:)                   :: Y1
  complex(8),dimension(size(Y1,1),size(Y1,2)),optional :: Y2
  real(8),dimension(size(Y1,2))               :: X
  logical,optional                    :: append
  logical                             :: check,rw
  rw=.false.;if(present(append))rw=append
  Ny1=size(Y1,1) ; Ny2=size(Y1,2)
  open(719,file=adjustl(trim(pname)))
  if(present(Y2))then
     do i=1,Ny1
        do j=1,Ny2
           write(719,"(F21.12,4(F21.12))")X(j),dimag(Y1(i,j)),dreal(Y1(i,j)),&
                dimag(Y2(i,j)),dreal(Y2(i,j))
        enddo
        write(719,*)""
     enddo
  else
     do i=1,Ny1
        do j=1,Ny2
           write(719,*)X(j),dimag(Y1(i,j)),dreal(Y1(i,j))
        enddo
        write(719,*)""
     enddo
  endif
  close(719)
end subroutine splotM_RC


!----------------------------
!----------------------------
!----------------------------



subroutine splotA3_II(pname,X,Y1,Y2,append)
  integer                                                      :: i,j,k,Ny1,Ny2,Ny3
  character(len=*)                                             :: pname
  integer,dimension(:,:,:)                                     :: Y1
  integer,dimension(size(Y1,1),size(Y1,2),size(Y1,3)),optional :: Y2
  integer,dimension(size(Y1,3))                                :: X
  logical,optional                    :: append
  logical                             :: check,rw
  rw=.false.;if(present(append))rw=append
  Ny1=size(Y1,1) ; Ny2=size(Y1,2); Ny3=size(Y1,3)
  open(719,file=adjustl(trim(pname)))
  if(present(Y2))then
     do i=1,Ny1
        do j=1,Ny2
           do k=1,Ny3
              write(719,*)X(k),Y1(i,j,k),Y2(i,j,k)
           enddo
           write(719,*)""
        enddo
     enddo
  else
     do i=1,Ny1
        do j=1,Ny2
           do k=1,Ny3
              write(719,*)X(k),Y1(i,j,k)
           enddo
           write(719,*)""
        enddo
     enddo
  endif
  close(719)
end subroutine splotA3_II

subroutine splotA3_IR(pname,X,Y1,Y2,append)
  integer                                                      :: i,j,k,Ny1,Ny2,Ny3
  character(len=*)                                             :: pname
  real(8),dimension(:,:,:)                                     :: Y1
  real(8),dimension(size(Y1,1),size(Y1,2),size(Y1,3)),optional :: Y2
  integer,dimension(size(Y1,3))                                :: X
  logical,optional                    :: append
  logical                             :: check,rw
  rw=.false.;if(present(append))rw=append
  Ny1=size(Y1,1) ; Ny2=size(Y1,2); Ny3=size(Y1,3)
  open(719,file=adjustl(trim(pname)))
  if(present(Y2))then
     do i=1,Ny1
        do j=1,Ny2
           do k=1,Ny3
              write(719,*)X(k),Y1(i,j,k),Y2(i,j,k)
           enddo
           write(719,*)""
        enddo
     enddo
  else
     do i=1,Ny1
        do j=1,Ny2
           do k=1,Ny3
              write(719,*)X(k),Y1(i,j,k)
           enddo
           write(719,*)""
        enddo
     enddo
  endif
  close(719)
end subroutine splotA3_IR

subroutine splotA3_IC(pname,X,Y1,Y2,append)
  integer                                                         :: i,j,k,Ny1,Ny2,Ny3
  character(len=*)                                                :: pname
  complex(8),dimension(:,:,:)                                     :: Y1
  complex(8),dimension(size(Y1,1),size(Y1,2),size(Y1,3)),optional :: Y2
  integer,dimension(size(Y1,3))                                   :: X
  logical,optional                    :: append
  logical                             :: check,rw
  rw=.false.;if(present(append))rw=append
  Ny1=size(Y1,1) ; Ny2=size(Y1,2); Ny3=size(Y1,3)
  open(719,file=adjustl(trim(pname)))
  if(present(Y2))then
     do i=1,Ny1
        do j=1,Ny2
           do k=1,Ny3
              write(719,"(I15,4(F21.12))")X(k),dimag(Y1(i,j,k)),dreal(Y1(i,j,k)),dimag(Y2(i,j,k)),dreal(Y2(i,j,k))
           enddo
           write(719,*)""
        enddo
     enddo
  else
     do i=1,Ny1
        do j=1,Ny2
           do k=1,Ny3
              write(719,*)X(k),dimag(Y1(i,j,k)),dreal(Y1(i,j,k))
           enddo
           write(719,*)""
        enddo
     enddo
  endif
  close(719)
end subroutine splotA3_IC

!----------------------------


subroutine splotA3_RI(pname,X,Y1,Y2,append)
  integer                                                      :: i,j,k,Ny1,Ny2,Ny3
  character(len=*)                                             :: pname
  integer,dimension(:,:,:)                                     :: Y1
  integer,dimension(size(Y1,1),size(Y1,2),size(Y1,3)),optional :: Y2
  real(8),dimension(size(Y1,3))                                :: X
  logical,optional                    :: append
  logical                             :: check,rw
  rw=.false.;if(present(append))rw=append
  Ny1=size(Y1,1) ; Ny2=size(Y1,2); Ny3=size(Y1,3)
  open(719,file=adjustl(trim(pname)))
  if(present(Y2))then
     do i=1,Ny1
        do j=1,Ny2
           do k=1,Ny3
              write(719,*)X(k),Y1(i,j,k),Y2(i,j,k)
           enddo
           write(719,*)""
        enddo
     enddo
  else
     do i=1,Ny1
        do j=1,Ny2
           do k=1,Ny3
              write(719,*)X(k),Y1(i,j,k)
           enddo
           write(719,*)""
        enddo
     enddo
  endif
  close(719)
end subroutine splotA3_RI

subroutine splotA3_RR(pname,X,Y1,Y2,append)
  integer                                                      :: i,j,k,Ny1,Ny2,Ny3
  character(len=*)                                             :: pname
  real(8),dimension(:,:,:)                                     :: Y1
  real(8),dimension(size(Y1,1),size(Y1,2),size(Y1,3)),optional :: Y2
  real(8),dimension(size(Y1,3))                                :: X
  logical,optional                    :: append
  logical                             :: check,rw
  rw=.false.;if(present(append))rw=append
  Ny1=size(Y1,1) ; Ny2=size(Y1,2); Ny3=size(Y1,3)
  open(719,file=adjustl(trim(pname)))
  if(present(Y2))then
     do i=1,Ny1
        do j=1,Ny2
           do k=1,Ny3
              write(719,*)X(k),Y1(i,j,k),Y2(i,j,k)
           enddo
           write(719,*)""
        enddo
     enddo
  else
     do i=1,Ny1
        do j=1,Ny2
           do k=1,Ny3
              write(719,*)X(k),Y1(i,j,k)
           enddo
           write(719,*)""
        enddo
     enddo
  endif
  close(719)
end subroutine splotA3_RR

subroutine splotA3_RC(pname,X,Y1,Y2,append)
  integer                                                         :: i,j,k,Ny1,Ny2,Ny3
  character(len=*)                                                :: pname
  complex(8),dimension(:,:,:)                                     :: Y1
  complex(8),dimension(size(Y1,1),size(Y1,2),size(Y1,3)),optional :: Y2
  real(8),dimension(size(Y1,3))                                   :: X
  logical,optional                    :: append
  logical                             :: check,rw
  rw=.false.;if(present(append))rw=append
  Ny1=size(Y1,1) ; Ny2=size(Y1,2); Ny3=size(Y1,3)
  open(719,file=adjustl(trim(pname)))
  if(present(Y2))then
     do i=1,Ny1
        do j=1,Ny2
           do k=1,Ny3
              write(719,"(F21.12,4(F21.12))")X(k),dimag(Y1(i,j,k)),dreal(Y1(i,j,k)),dimag(Y2(i,j,k)),dreal(Y2(i,j,k))
           enddo
           write(719,*)""
        enddo
     enddo
  else
     do i=1,Ny1
        do j=1,Ny2
           do k=1,Ny3
              write(719,*)X(k),dimag(Y1(i,j,k)),dreal(Y1(i,j,k))
           enddo
           write(719,*)""
        enddo
     enddo
  endif
  close(719)
end subroutine splotA3_RC
!----------------------------
!----------------------------
!----------------------------
