!!Bool
subroutine MPI_Bcast_Bool_0(comm,data,root)
  integer,intent(in)          :: comm
  logical,intent(in)          :: data
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,1,MPI_LOGICAL,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Bool_0')
end subroutine MPI_Bcast_Bool_0
!
subroutine MPI_Bcast_Bool_1(comm,data,root)
  integer,intent(in)          :: comm
  logical,intent(in)          :: data(:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_LOGICAL,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Bool_1')
end subroutine MPI_Bcast_Bool_1
!
subroutine MPI_Bcast_Bool_2(comm,data,root)
  integer,intent(in)          :: comm
  logical,intent(in)          :: data(:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_LOGICAL,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Bool_2')
end subroutine MPI_Bcast_Bool_2
!
subroutine MPI_Bcast_Bool_3(comm,data,root)
  integer,intent(in)          :: comm
  logical,intent(in)          :: data(:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_LOGICAL,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Bool_3')
end subroutine MPI_Bcast_Bool_3
!
subroutine MPI_Bcast_Bool_4(comm,data,root)
  integer,intent(in)          :: comm
  logical,intent(in)          :: data(:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_LOGICAL,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Bool_4')
end subroutine MPI_Bcast_Bool_4
!
subroutine MPI_Bcast_Bool_5(comm,data,root)
  integer,intent(in)          :: comm
  logical,intent(in)          :: data(:,:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_LOGICAL,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Bool_5')
end subroutine MPI_Bcast_Bool_5
!
subroutine MPI_Bcast_Bool_6(comm,data,root)
  integer,intent(in)          :: comm
  logical,intent(in)          :: data(:,:,:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_LOGICAL,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Bool_6')
end subroutine MPI_Bcast_Bool_6
!
subroutine MPI_Bcast_Bool_7(comm,data,root)
  integer,intent(in)          :: comm
  logical,intent(in)          :: data(:,:,:,:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_LOGICAL,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Bool_7')
end subroutine MPI_Bcast_Bool_7





!! INTEGER
subroutine MPI_Bcast_Int_0(comm,data,root)
  integer,intent(in)          :: comm
  integer,intent(in)          :: data
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,1,MPI_INTEGER,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Int_0')
end subroutine MPI_Bcast_Int_0
!
subroutine MPI_Bcast_Int_1(comm,data,root)
  integer,intent(in)          :: comm
  integer,intent(in)          :: data(:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_INTEGER,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Int_1')
end subroutine MPI_Bcast_Int_1
!
subroutine MPI_Bcast_Int_2(comm,data,root)
  integer,intent(in)          :: comm
  integer,intent(in)          :: data(:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_INTEGER,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Int_2')
end subroutine MPI_Bcast_Int_2
!
subroutine MPI_Bcast_Int_3(comm,data,root)
  integer,intent(in)          :: comm
  integer,intent(in)          :: data(:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_INTEGER,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Int_3')
end subroutine MPI_Bcast_Int_3
!
subroutine MPI_Bcast_Int_4(comm,data,root)
  integer,intent(in)          :: comm
  integer,intent(in)          :: data(:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_INTEGER,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Int_4')
end subroutine MPI_Bcast_Int_4
!
subroutine MPI_Bcast_Int_5(comm,data,root)
  integer,intent(in)          :: comm
  integer,intent(in)          :: data(:,:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_INTEGER,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Int_5')
end subroutine MPI_Bcast_Int_5
!
subroutine MPI_Bcast_Int_6(comm,data,root)
  integer,intent(in)          :: comm
  integer,intent(in)          :: data(:,:,:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_INTEGER,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Int_6')
end subroutine MPI_Bcast_Int_6
!
subroutine MPI_Bcast_Int_7(comm,data,root)
  integer,intent(in)          :: comm
  integer,intent(in)          :: data(:,:,:,:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_INTEGER,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Int_7')
end subroutine MPI_Bcast_Int_7





!! REAL8
subroutine MPI_Bcast_Dble_0(comm,data,root)
  integer,intent(in)          :: comm
  real(8),intent(in)          :: data
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,1,MPI_DOUBLE_PRECISION,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Dble_0')
end subroutine MPI_Bcast_Dble_0
!
subroutine MPI_Bcast_Dble_1(comm,data,root)
  integer,intent(in)          :: comm
  real(8),intent(in)          :: data(:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_DOUBLE_PRECISION,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Dble_1')
end subroutine MPI_Bcast_Dble_1
!
subroutine MPI_Bcast_Dble_2(comm,data,root)
  integer,intent(in)          :: comm
  real(8),intent(in)          :: data(:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_DOUBLE_PRECISION,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Dble_2')
end subroutine MPI_Bcast_Dble_2
!
subroutine MPI_Bcast_Dble_3(comm,data,root)
  integer,intent(in)          :: comm
  real(8),intent(in)          :: data(:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_DOUBLE_PRECISION,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Dble_3')
end subroutine MPI_Bcast_Dble_3
!
subroutine MPI_Bcast_Dble_4(comm,data,root)
  integer,intent(in)          :: comm
  real(8),intent(in)          :: data(:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_DOUBLE_PRECISION,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Dble_4')
end subroutine MPI_Bcast_Dble_4
!
subroutine MPI_Bcast_Dble_5(comm,data,root)
  integer,intent(in)          :: comm
  real(8),intent(in)          :: data(:,:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_DOUBLE_PRECISION,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Dble_5')
end subroutine MPI_Bcast_Dble_5
!
subroutine MPI_Bcast_Dble_6(comm,data,root)
  integer,intent(in)          :: comm
  real(8),intent(in)          :: data(:,:,:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_DOUBLE_PRECISION,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Dble_6')
end subroutine MPI_Bcast_Dble_6
!
subroutine MPI_Bcast_Dble_7(comm,data,root)
  integer,intent(in)          :: comm
  real(8),intent(in)          :: data(:,:,:,:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_DOUBLE_PRECISION,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Dble_7')
end subroutine MPI_Bcast_Dble_7





!!CMPLX8
subroutine MPI_Bcast_Cmplx_0(comm,data,root)
  integer,intent(in)          :: comm
  complex(8),intent(in)       :: data
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,1,MPI_DOUBLE_COMPLEX,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Cmplx_0')
end subroutine MPI_Bcast_Cmplx_0
!
subroutine MPI_Bcast_Cmplx_1(comm,data,root)
  integer,intent(in)          :: comm
  complex(8),intent(in)       :: data(:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_DOUBLE_COMPLEX,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Cmplx_1')
end subroutine MPI_Bcast_Cmplx_1
!
subroutine MPI_Bcast_Cmplx_2(comm,data,root)
  integer,intent(in)          :: comm
  complex(8),intent(in)       :: data(:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_DOUBLE_COMPLEX,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Cmplx_2')
end subroutine MPI_Bcast_Cmplx_2
!
subroutine MPI_Bcast_Cmplx_3(comm,data,root)
  integer,intent(in)          :: comm
  complex(8),intent(in)       :: data(:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_DOUBLE_COMPLEX,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Cmplx_3')
end subroutine MPI_Bcast_Cmplx_3
!
subroutine MPI_Bcast_Cmplx_4(comm,data,root)
  integer,intent(in)          :: comm
  complex(8),intent(in)       :: data(:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_DOUBLE_COMPLEX,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Cmplx_4')
end subroutine MPI_Bcast_Cmplx_4
!
subroutine MPI_Bcast_Cmplx_5(comm,data,root)
  integer,intent(in)          :: comm
  complex(8),intent(in)       :: data(:,:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_DOUBLE_COMPLEX,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Cmplx_5')
end subroutine MPI_Bcast_Cmplx_5
!
subroutine MPI_Bcast_Cmplx_6(comm,data,root)
  integer,intent(in)          :: comm
  complex(8),intent(in)       :: data(:,:,:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_DOUBLE_COMPLEX,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Cmplx_6')
end subroutine MPI_Bcast_Cmplx_6
!
subroutine MPI_Bcast_Cmplx_7(comm,data,root)
  integer,intent(in)          :: comm
  complex(8),intent(in)       :: data(:,:,:,:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_DOUBLE_COMPLEX,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Cmplx_7')
end subroutine MPI_Bcast_Cmplx_7
