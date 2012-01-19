
!+-----------------------------------------------------------------+
!PROGRAM  : 
!TYPE     : Subroutine
!PURPOSE  : 
!+-----------------------------------------------------------------+
subroutine shiftM_fw_C(Gin,Gout,step_)
  complex(8),dimension(0:,0:) :: Gin
  complex(8),dimension(:,:)   :: Gout
  integer :: i,j,Ndim1,Ndim2,step
  integer,optional :: step_
  step=1;if(present(step_))step=step_
  Ndim1=size(Gout,1) ; Ndim2=size(Gout,2)
  forall(i=1:Ndim1,j=1:Ndim2)Gout(i,j)=Gin(i-step,j-step)
end subroutine shiftM_fw_C
!-----------------------------
!-----------------------------
!-----------------------------
subroutine shiftM_fw_R(Gin,Gout,step_)
  real(8),dimension(0:,0:) :: Gin
  real(8),dimension(:,:)   :: Gout
  integer :: i,j,Ndim1,Ndim2,step
  integer,optional :: step_
  step=1;if(present(step_))step=step_
  Ndim1=size(Gout,1) ; Ndim2=size(Gout,2)
  forall(i=1:Ndim1,j=1:Ndim2)Gout(i,j)=Gin(i-step,j-step)
end subroutine shiftM_fw_R
!-----------------------------
!-----------------------------
!-----------------------------
subroutine shiftA_fw_C(Ain,Aout,step_)
  complex(8),dimension(0:) :: Ain
  complex(8),dimension(:)  :: Aout
  integer :: i,Ndim1,step
  integer,optional :: step_
  step=1;if(present(step_))step=step_
  Ndim1=size(Aout,1) 
  forall(i=1:Ndim1)Aout(i)=Ain(i-step)
end subroutine shiftA_fw_C
!-----------------------------
!-----------------------------
!-----------------------------
subroutine shiftA_fw_R(Ain,Aout,step_)
  real(8),dimension(0:) :: Ain
  real(8),dimension(:)  :: Aout
  integer :: i,Ndim1,step
  integer,optional :: step_
  step=1;if(present(step_))step=step_
  Ndim1=size(Aout,1) 
  forall(i=1:Ndim1)Aout(i)=Ain(i-step)
end subroutine shiftA_fw_R
!********************************************************************
!********************************************************************
!********************************************************************












!+-----------------------------------------------------------------+
!PROGRAM  : 
!TYPE     : Subroutine
!PURPOSE  : 
!+-----------------------------------------------------------------+
subroutine shiftM_bw_C(Gin,Gout,step_)
  complex(8),dimension(0:,0:) :: Gout
  complex(8),dimension(:,:)   :: Gin
  integer :: i,j,Ndim1,Ndim2,step
  integer,optional :: step_
  step=1;if(present(step_))step=step_
  Ndim1=size(Gin,1) ; Ndim2=size(Gin,2)
  forall(i=1:Ndim1,j=1:Ndim2)Gout(i-step,j-step)=Gin(i,j)
end subroutine shiftM_bw_C
!-----------------------------
!-----------------------------
!-----------------------------
subroutine shiftM_bw_R(Gin,Gout,step_)
  real(8),dimension(0:,0:) :: Gout
  real(8),dimension(:,:)   :: Gin
  integer :: i,j,Ndim1,Ndim2,step
  integer,optional :: step_
  step=1;if(present(step_))step=step_
  Ndim1=size(Gin,1) ; Ndim2=size(Gin,2)
  forall(i=1:Ndim1,j=1:Ndim2)Gout(i-step,j-step)=Gin(i,j)
end subroutine shiftM_bw_R
!-----------------------------
!-----------------------------
!-----------------------------
subroutine shiftA_bw_C(Ain,Aout,step_)
  complex(8),dimension(0:) :: Aout
  complex(8),dimension(:)  :: Ain
  integer :: i,Ndim1,step
  integer,optional :: step_
  step=1;if(present(step_))step=step_
  Ndim1=size(Ain,1)
  forall(i=1:Ndim1)Aout(i-step)=Ain(i)
end subroutine shiftA_bw_C
!-----------------------------
!-----------------------------
!-----------------------------
subroutine shiftA_bw_R(Ain,Aout,step_)
  real(8),dimension(0:) :: Aout
  real(8),dimension(:)  :: Ain
  integer :: i,Ndim1,step
  integer,optional :: step_
  step=1;if(present(step_))step=step_
  Ndim1=size(Ain,1)
  forall(i=1:Ndim1)Aout(i-step)=Ain(i)
end subroutine shiftA_bw_R
!********************************************************************
!********************************************************************
!********************************************************************
