!!BOOL
subroutine MPI_Allgather_Bool_0(comm,send,data,root)
  integer,intent(in)          :: comm
  logical,intent(inout)       :: data
  logical,intent(in)          :: send
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLGATHER(send,1,MPI_LOGICAL,data,1,MPI_LOGICAL,comm,ierr)
  call Error_MPI(sub='MPI_Allgather_Bool_0')
end subroutine MPI_Allgather_Bool_0
!
subroutine MPI_Allgather_Bool_1(comm,send,data,root)
  integer,intent(in)          :: comm
  logical,intent(inout)       :: data(:)
  logical,intent(in)          :: send(:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLGATHER(send,size(send),MPI_LOGICAL,data,size(data),MPI_LOGICAL,comm,ierr)
  call Error_MPI(sub='MPI_Allgather_Bool_1')
end subroutine MPI_Allgather_Bool_1
!
subroutine MPI_Allgather_Bool_2(comm,send,data,root)
  integer,intent(in)          :: comm
  logical,intent(inout)       :: data(:,:)
  logical,intent(in)          :: send(:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLGATHER(send,size(send),MPI_LOGICAL,data,size(data),MPI_LOGICAL,comm,ierr)
  call Error_MPI(sub='MPI_Allgather_Bool_2')
end subroutine MPI_Allgather_Bool_2
!
subroutine MPI_Allgather_Bool_3(comm,send,data,root)
  integer,intent(in)          :: comm
  logical,intent(inout)       :: data(:,:,:)
  logical,intent(in)          :: send(:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLGATHER(send,size(send),MPI_LOGICAL,data,size(data),MPI_LOGICAL,comm,ierr)
  call Error_MPI(sub='MPI_Allgather_Bool_3')
end subroutine MPI_Allgather_Bool_3
!
subroutine MPI_Allgather_Bool_4(comm,send,data,root)
  integer,intent(in)          :: comm
  logical,intent(inout)       :: data(:,:,:,:)
  logical,intent(in)          :: send(:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLGATHER(send,size(send),MPI_LOGICAL,data,size(data),MPI_LOGICAL,comm,ierr)
  call Error_MPI(sub='MPI_Allgather_Bool_4')
end subroutine MPI_Allgather_Bool_4
!
subroutine MPI_Allgather_Bool_5(comm,send,data,root)
  integer,intent(in)          :: comm
  logical,intent(inout)       :: data(:,:,:,:,:)
  logical,intent(in)          :: send(:,:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLGATHER(send,size(send),MPI_LOGICAL,data,size(data),MPI_LOGICAL,comm,ierr)
  call Error_MPI(sub='MPI_Allgather_Bool_5')
end subroutine MPI_Allgather_Bool_5
!
subroutine MPI_Allgather_Bool_6(comm,send,data,root)
  integer,intent(in)          :: comm
  logical,intent(inout)       :: data(:,:,:,:,:,:)
  logical,intent(in)          :: send(:,:,:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLGATHER(send,size(send),MPI_LOGICAL,data,size(data),MPI_LOGICAL,comm,ierr)
  call Error_MPI(sub='MPI_Allgather_Bool_6')
end subroutine MPI_Allgather_Bool_6
!
subroutine MPI_Allgather_Bool_7(comm,send,data,root)
  integer,intent(in)          :: comm
  logical,intent(inout)       :: data(:,:,:,:,:,:,:)
  logical,intent(in)          :: send(:,:,:,:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLGATHER(send,size(send),MPI_LOGICAL,data,size(data),MPI_LOGICAL,comm,ierr)
  call Error_MPI(sub='MPI_Allgather_Bool_7')
end subroutine MPI_Allgather_Bool_7





!!INTEGER
subroutine MPI_Allgather_Int_0(comm,send,data,root)
  integer,intent(in)          :: comm
  integer,intent(inout)       :: data
  integer,intent(in)          :: send
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLGATHER(send,1,MPI_INTEGER,data,1,MPI_INTEGER,comm,ierr)
  call Error_MPI(sub='MPI_Allgather_Int_0')
end subroutine MPI_Allgather_Int_0
!
subroutine MPI_Allgather_Int_1(comm,send,data,root)
  integer,intent(in)          :: comm
  integer,intent(inout)       :: data(:)
  integer,intent(in)          :: send(:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLGATHER(send,size(send),MPI_INTEGER,data,size(data),MPI_INTEGER,comm,ierr)
  call Error_MPI(sub='MPI_Allgather_Int_1')
end subroutine MPI_Allgather_Int_1
!
subroutine MPI_Allgather_Int_2(comm,send,data,root)
  integer,intent(in)          :: comm
  integer,intent(inout)       :: data(:,:)
  integer,intent(in)          :: send(:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLGATHER(send,size(send),MPI_INTEGER,data,size(data),MPI_INTEGER,comm,ierr)
  call Error_MPI(sub='MPI_Allgather_Int_2')
end subroutine MPI_Allgather_Int_2
!
subroutine MPI_Allgather_Int_3(comm,send,data,root)
  integer,intent(in)          :: comm
  integer,intent(inout)       :: data(:,:,:)
  integer,intent(in)          :: send(:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLGATHER(send,size(send),MPI_INTEGER,data,size(data),MPI_INTEGER,comm,ierr)
  call Error_MPI(sub='MPI_Allgather_Int_3')
end subroutine MPI_Allgather_Int_3
!
subroutine MPI_Allgather_Int_4(comm,send,data,root)
  integer,intent(in)          :: comm
  integer,intent(inout)       :: data(:,:,:,:)
  integer,intent(in)          :: send(:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLGATHER(send,size(send),MPI_INTEGER,data,size(data),MPI_INTEGER,comm,ierr)
  call Error_MPI(sub='MPI_Allgather_Int_4')
end subroutine MPI_Allgather_Int_4
!
subroutine MPI_Allgather_Int_5(comm,send,data,root)
  integer,intent(in)          :: comm
  integer,intent(inout)       :: data(:,:,:,:,:)
  integer,intent(in)          :: send(:,:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLGATHER(send,size(send),MPI_INTEGER,data,size(data),MPI_INTEGER,comm,ierr)
  call Error_MPI(sub='MPI_Allgather_Int_5')
end subroutine MPI_Allgather_Int_5
!
subroutine MPI_Allgather_Int_6(comm,send,data,root)
  integer,intent(in)          :: comm
  integer,intent(inout)       :: data(:,:,:,:,:,:)
  integer,intent(in)          :: send(:,:,:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLGATHER(send,size(send),MPI_INTEGER,data,size(data),MPI_INTEGER,comm,ierr)
  call Error_MPI(sub='MPI_Allgather_Int_6')
end subroutine MPI_Allgather_Int_6
!
subroutine MPI_Allgather_Int_7(comm,send,data,root)
  integer,intent(in)          :: comm
  integer,intent(inout)       :: data(:,:,:,:,:,:,:)
  integer,intent(in)          :: send(:,:,:,:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLGATHER(send,size(send),MPI_INTEGER,data,size(data),MPI_INTEGER,comm,ierr)
  call Error_MPI(sub='MPI_Allgather_Int_7')
end subroutine MPI_Allgather_Int_7







!!REAL8
subroutine MPI_Allgather_Dble_0(comm,send,data,root)
  integer,intent(in)          :: comm
  real(8),intent(inout)       :: data
  real(8),intent(in)          :: send
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLGATHER(send,1,MPI_DOUBLE_PRECISION,data,1,MPI_DOUBLE_PRECISION,comm,ierr)
  call Error_MPI(sub='MPI_Allgather_Dble_0')
end subroutine MPI_Allgather_Dble_0
!
subroutine MPI_Allgather_Dble_1(comm,send,data,root)
  integer,intent(in)          :: comm
  real(8),intent(inout)       :: data(:)
  real(8),intent(in)          :: send(:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_PRECISION,data,size(data),MPI_DOUBLE_PRECISION,comm,ierr)
  call Error_MPI(sub='MPI_Allgather_Dble_1')
end subroutine MPI_Allgather_Dble_1
!
subroutine MPI_Allgather_Dble_2(comm,send,data,root)
  integer,intent(in)          :: comm
  real(8),intent(inout)       :: data(:,:)
  real(8),intent(in)          :: send(:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_PRECISION,data,size(data),MPI_DOUBLE_PRECISION,comm,ierr)
  call Error_MPI(sub='MPI_Allgather_Dble_2')
end subroutine MPI_Allgather_Dble_2
!
subroutine MPI_Allgather_Dble_3(comm,send,data,root)
  integer,intent(in)          :: comm
  real(8),intent(inout)       :: data(:,:,:)
  real(8),intent(in)          :: send(:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_PRECISION,data,size(data),MPI_DOUBLE_PRECISION,comm,ierr)
  call Error_MPI(sub='MPI_Allgather_Dble_3')
end subroutine MPI_Allgather_Dble_3
!
subroutine MPI_Allgather_Dble_4(comm,send,data,root)
  integer,intent(in)          :: comm
  real(8),intent(inout)       :: data(:,:,:,:)
  real(8),intent(in)          :: send(:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_PRECISION,data,size(data),MPI_DOUBLE_PRECISION,comm,ierr)
  call Error_MPI(sub='MPI_Allgather_Dble_4')
end subroutine MPI_Allgather_Dble_4
!
subroutine MPI_Allgather_Dble_5(comm,send,data,root)
  integer,intent(in)          :: comm
  real(8),intent(inout)       :: data(:,:,:,:,:)
  real(8),intent(in)          :: send(:,:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_PRECISION,data,size(data),MPI_DOUBLE_PRECISION,comm,ierr)
  call Error_MPI(sub='MPI_Allgather_Dble_5')
end subroutine MPI_Allgather_Dble_5
!
subroutine MPI_Allgather_Dble_6(comm,send,data,root)
  integer,intent(in)          :: comm
  real(8),intent(inout)       :: data(:,:,:,:,:,:)
  real(8),intent(in)          :: send(:,:,:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_PRECISION,data,size(data),MPI_DOUBLE_PRECISION,comm,ierr)
  call Error_MPI(sub='MPI_Allgather_Dble_6')
end subroutine MPI_Allgather_Dble_6
!
subroutine MPI_Allgather_Dble_7(comm,send,data,root)
  integer,intent(in)          :: comm
  real(8),intent(inout)       :: data(:,:,:,:,:,:,:)
  real(8),intent(in)          :: send(:,:,:,:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_PRECISION,data,size(data),MPI_DOUBLE_PRECISION,comm,ierr)
  call Error_MPI(sub='MPI_Allgather_Dble_7')
end subroutine MPI_Allgather_Dble_7





!!CMPLX8
subroutine MPI_Allgather_Cmplx_0(comm,send,data,root)
  integer,intent(in)          :: comm
  complex(8),intent(inout)    :: data
  complex(8),intent(in)       :: send
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLGATHER(send,1,MPI_DOUBLE_COMPLEX,data,1,MPI_DOUBLE_COMPLEX,comm,ierr)
  call Error_MPI(sub='MPI_Allgather_Cmplx_0')
end subroutine MPI_Allgather_Cmplx_0
!
subroutine MPI_Allgather_Cmplx_1(comm,send,data,root)
  integer,intent(in)          :: comm
  complex(8),intent(inout)    :: data(:)
  complex(8),intent(in)       :: send(:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_COMPLEX,data,size(data),MPI_DOUBLE_COMPLEX,comm,ierr)
  call Error_MPI(sub='MPI_Allgather_Cmplx_1')
end subroutine MPI_Allgather_Cmplx_1
!
subroutine MPI_Allgather_Cmplx_2(comm,send,data,root)
  integer,intent(in)          :: comm
  complex(8),intent(inout)    :: data(:,:)
  complex(8),intent(in)       :: send(:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_COMPLEX,data,size(data),MPI_DOUBLE_COMPLEX,comm,ierr)
  call Error_MPI(sub='MPI_Allgather_Cmplx_2')
end subroutine MPI_Allgather_Cmplx_2
!
subroutine MPI_Allgather_Cmplx_3(comm,send,data,root)
  integer,intent(in)          :: comm
  complex(8),intent(inout)    :: data(:,:,:)
  complex(8),intent(in)       :: send(:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_COMPLEX,data,size(data),MPI_DOUBLE_COMPLEX,comm,ierr)
  call Error_MPI(sub='MPI_Allgather_Cmplx_3')
end subroutine MPI_Allgather_Cmplx_3
!
subroutine MPI_Allgather_Cmplx_4(comm,send,data,root)
  integer,intent(in)          :: comm
  complex(8),intent(inout)    :: data(:,:,:,:)
  complex(8),intent(in)       :: send(:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_COMPLEX,data,size(data),MPI_DOUBLE_COMPLEX,comm,ierr)
  call Error_MPI(sub='MPI_Allgather_Cmplx_4')
end subroutine MPI_Allgather_Cmplx_4
!
subroutine MPI_Allgather_Cmplx_5(comm,send,data,root)
  integer,intent(in)          :: comm
  complex(8),intent(inout)    :: data(:,:,:,:,:)
  complex(8),intent(in)       :: send(:,:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_COMPLEX,data,size(data),MPI_DOUBLE_COMPLEX,comm,ierr)
  call Error_MPI(sub='MPI_Allgather_Cmplx_5')
end subroutine MPI_Allgather_Cmplx_5
!
subroutine MPI_Allgather_Cmplx_6(comm,send,data,root)
  integer,intent(in)          :: comm
  complex(8),intent(inout)    :: data(:,:,:,:,:,:)
  complex(8),intent(in)       :: send(:,:,:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_COMPLEX,data,size(data),MPI_DOUBLE_COMPLEX,comm,ierr)
  call Error_MPI(sub='MPI_Allgather_Cmplx_6')
end subroutine MPI_Allgather_Cmplx_6
!
subroutine MPI_Allgather_Cmplx_7(comm,send,data,root)
  integer,intent(in)          :: comm
  complex(8),intent(inout)    :: data(:,:,:,:,:,:,:)
  complex(8),intent(in)       :: send(:,:,:,:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLGATHER(send,size(send),MPI_DOUBLE_COMPLEX,data,size(data),MPI_DOUBLE_COMPLEX,comm,ierr)
  call Error_MPI(sub='MPI_Allgather_Cmplx_7')
end subroutine MPI_Allgather_Cmplx_7
