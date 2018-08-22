!!BOOL
subroutine MPI_Allreduce_Bool_0(comm,send,data,root)
  integer,intent(in)          :: comm
  logical,intent(inout)       :: data
  logical,intent(in)          :: send
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  call MPI_ALLREDUCE(send,data,1,MPI_LOGICAL,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Bool_0')
end subroutine MPI_Allreduce_Bool_0
!
subroutine MPI_Allreduce_Bool_1(comm,send,data,root)
  integer,intent(in)          :: comm
  logical,intent(inout)       :: data(:)
  logical,intent(in)          :: send(:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  call MPI_ALLREDUCE(send,data,size(data),MPI_LOGICAL,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Bool_1')
end subroutine MPI_Allreduce_Bool_1
!
subroutine MPI_Allreduce_Bool_2(comm,send,data,root)
  integer,intent(in)          :: comm
  logical,intent(inout)       :: data(:,:)
  logical,intent(in)          :: send(:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  call MPI_ALLREDUCE(send,data,size(data),MPI_LOGICAL,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Bool_2')
end subroutine MPI_Allreduce_Bool_2
!
subroutine MPI_Allreduce_Bool_3(comm,send,data,root)
  integer,intent(in)          :: comm
  logical,intent(inout)       :: data(:,:,:)
  logical,intent(in)          :: send(:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  call MPI_ALLREDUCE(send,data,size(data),MPI_LOGICAL,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Bool_3')
end subroutine MPI_Allreduce_Bool_3
!
subroutine MPI_Allreduce_Bool_4(comm,send,data,root)
  integer,intent(in)          :: comm
  logical,intent(inout)       :: data(:,:,:,:)
  logical,intent(in)          :: send(:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  call MPI_ALLREDUCE(send,data,size(data),MPI_LOGICAL,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Bool_4')
end subroutine MPI_Allreduce_Bool_4
!
subroutine MPI_Allreduce_Bool_5(comm,send,data,root)
  integer,intent(in)          :: comm
  logical,intent(inout)       :: data(:,:,:,:,:)
  logical,intent(in)          :: send(:,:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  call MPI_ALLREDUCE(send,data,size(data),MPI_LOGICAL,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Bool_5')
end subroutine MPI_Allreduce_Bool_5
!
subroutine MPI_Allreduce_Bool_6(comm,send,data,root)
  integer,intent(in)          :: comm
  logical,intent(inout)       :: data(:,:,:,:,:,:)
  logical,intent(in)          :: send(:,:,:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  call MPI_ALLREDUCE(send,data,size(data),MPI_LOGICAL,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Bool_6')
end subroutine MPI_Allreduce_Bool_6
!
subroutine MPI_Allreduce_Bool_7(comm,send,data,root)
  integer,intent(in)          :: comm
  logical,intent(inout)       :: data(:,:,:,:,:,:,:)
  logical,intent(in)          :: send(:,:,:,:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  call MPI_ALLREDUCE(send,data,size(data),MPI_LOGICAL,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Bool_7')
end subroutine MPI_Allreduce_Bool_7




!!INTEGER
subroutine MPI_Allreduce_Int_0(comm,send,data,root)
  integer,intent(in)          :: comm
  integer,intent(inout)       :: data
  integer,intent(in)          :: send
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  call MPI_ALLREDUCE(send,data,1,MPI_INTEGER,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Int_0')
end subroutine MPI_Allreduce_Int_0
!
subroutine MPI_Allreduce_Int_1(comm,send,data,root)
  integer,intent(in)          :: comm
  integer,intent(inout)       :: data(:)
  integer,intent(in)          :: send(:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  call MPI_ALLREDUCE(send,data,size(data),MPI_INTEGER,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Int_1')
end subroutine MPI_Allreduce_Int_1
!
subroutine MPI_Allreduce_Int_2(comm,send,data,root)
  integer,intent(in)          :: comm
  integer,intent(inout)       :: data(:,:)
  integer,intent(in)          :: send(:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  call MPI_ALLREDUCE(send,data,size(data),MPI_INTEGER,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Int_2')
end subroutine MPI_Allreduce_Int_2
!
subroutine MPI_Allreduce_Int_3(comm,send,data,root)
  integer,intent(in)          :: comm
  integer,intent(inout)       :: data(:,:,:)
  integer,intent(in)          :: send(:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  call MPI_ALLREDUCE(send,data,size(data),MPI_INTEGER,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Int_3')
end subroutine MPI_Allreduce_Int_3
!
subroutine MPI_Allreduce_Int_4(comm,send,data,root)
  integer,intent(in)          :: comm
  integer,intent(inout)       :: data(:,:,:,:)
  integer,intent(in)          :: send(:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  call MPI_ALLREDUCE(send,data,size(data),MPI_INTEGER,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Int_4')
end subroutine MPI_Allreduce_Int_4
!
subroutine MPI_Allreduce_Int_5(comm,send,data,root)
  integer,intent(in)          :: comm
  integer,intent(inout)       :: data(:,:,:,:,:)
  integer,intent(in)          :: send(:,:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  call MPI_ALLREDUCE(send,data,size(data),MPI_INTEGER,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Int_5')
end subroutine MPI_Allreduce_Int_5
!
subroutine MPI_Allreduce_Int_6(comm,send,data,root)
  integer,intent(in)          :: comm
  integer,intent(inout)       :: data(:,:,:,:,:,:)
  integer,intent(in)          :: send(:,:,:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  call MPI_ALLREDUCE(send,data,size(data),MPI_INTEGER,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Int_6')
end subroutine MPI_Allreduce_Int_6
!
subroutine MPI_Allreduce_Int_7(comm,send,data,root)
  integer,intent(in)          :: comm
  integer,intent(inout)       :: data(:,:,:,:,:,:,:)
  integer,intent(in)          :: send(:,:,:,:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  call MPI_ALLREDUCE(send,data,size(data),MPI_INTEGER,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Int_7')
end subroutine MPI_Allreduce_Int_7




!!REAL8
subroutine MPI_Allreduce_Dble_0(comm,send,data,root)
  integer,intent(in)          :: comm
  real(8),intent(inout)       :: data
  real(8),intent(in)          :: send
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  call MPI_ALLREDUCE(send,data,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Dble_0')
end subroutine MPI_Allreduce_Dble_0
!
subroutine MPI_Allreduce_Dble_1(comm,send,data,root)
  integer,intent(in)          :: comm
  real(8),intent(inout)       :: data(:)
  real(8),intent(in)          :: send(:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Dble_1')
end subroutine MPI_Allreduce_Dble_1
!
subroutine MPI_Allreduce_Dble_2(comm,send,data,root)
  integer,intent(in)          :: comm
  real(8),intent(inout)       :: data(:,:)
  real(8),intent(in)          :: send(:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Dble_2')
end subroutine MPI_Allreduce_Dble_2
!
subroutine MPI_Allreduce_Dble_3(comm,send,data,root)
  integer,intent(in)          :: comm
  real(8),intent(inout)       :: data(:,:,:)
  real(8),intent(in)          :: send(:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Dble_3')
end subroutine MPI_Allreduce_Dble_3
!
subroutine MPI_Allreduce_Dble_4(comm,send,data,root)
  integer,intent(in)          :: comm
  real(8),intent(inout)       :: data(:,:,:,:)
  real(8),intent(in)          :: send(:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Dble_4')
end subroutine MPI_Allreduce_Dble_4
!
subroutine MPI_Allreduce_Dble_5(comm,send,data,root)
  integer,intent(in)          :: comm
  real(8),intent(inout)       :: data(:,:,:,:,:)
  real(8),intent(in)          :: send(:,:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Dble_5')
end subroutine MPI_Allreduce_Dble_5
!
subroutine MPI_Allreduce_Dble_6(comm,send,data,root)
  integer,intent(in)          :: comm
  real(8),intent(inout)       :: data(:,:,:,:,:,:)
  real(8),intent(in)          :: send(:,:,:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Dble_6')
end subroutine MPI_Allreduce_Dble_6
!
subroutine MPI_Allreduce_Dble_7(comm,send,data,root)
  integer,intent(in)          :: comm
  real(8),intent(inout)       :: data(:,:,:,:,:,:,:)
  real(8),intent(in)          :: send(:,:,:,:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Dble_7')
end subroutine MPI_Allreduce_Dble_7





!!CMPLX8
subroutine MPI_Allreduce_Cmplx_0(comm,send,data,root)
  integer,intent(in)          :: comm
  complex(8),intent(inout)    :: data
  complex(8),intent(in)       :: send
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  call MPI_ALLREDUCE(send,data,1,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Cmplx_0')
end subroutine MPI_Allreduce_Cmplx_0
!
subroutine MPI_Allreduce_Cmplx_1(comm,send,data,root)
  integer,intent(in)          :: comm
  complex(8),intent(inout)    :: data(:)
  complex(8),intent(in)       :: send(:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Cmplx_1')
end subroutine MPI_Allreduce_Cmplx_1
!
subroutine MPI_Allreduce_Cmplx_2(comm,send,data,root)
  integer,intent(in)          :: comm
  complex(8),intent(inout)    :: data(:,:)
  complex(8),intent(in)       :: send(:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Cmplx_2')
end subroutine MPI_Allreduce_Cmplx_2
!
subroutine MPI_Allreduce_Cmplx_3(comm,send,data,root)
  integer,intent(in)          :: comm
  complex(8),intent(inout)    :: data(:,:,:)
  complex(8),intent(in)       :: send(:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Cmplx_3')
end subroutine MPI_Allreduce_Cmplx_3
!
subroutine MPI_Allreduce_Cmplx_4(comm,send,data,root)
  integer,intent(in)          :: comm
  complex(8),intent(inout)    :: data(:,:,:,:)
  complex(8),intent(in)       :: send(:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Cmplx_4')
end subroutine MPI_Allreduce_Cmplx_4
!
subroutine MPI_Allreduce_Cmplx_5(comm,send,data,root)
  integer,intent(in)          :: comm
  complex(8),intent(inout)    :: data(:,:,:,:,:)
  complex(8),intent(in)       :: send(:,:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Cmplx_5')
end subroutine MPI_Allreduce_Cmplx_5
!
subroutine MPI_Allreduce_Cmplx_6(comm,send,data,root)
  integer,intent(in)          :: comm
  complex(8),intent(inout)    :: data(:,:,:,:,:,:)
  complex(8),intent(in)       :: send(:,:,:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Cmplx_6')
end subroutine MPI_Allreduce_Cmplx_6
!
subroutine MPI_Allreduce_Cmplx_7(comm,send,data,root)
  integer,intent(in)          :: comm
  complex(8),intent(inout)    :: data(:,:,:,:,:,:,:)
  complex(8),intent(in)       :: send(:,:,:,:,:,:,:)
  integer,intent(in),optional :: root
  rank=0;if(present(root))rank=root
  call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Cmplx_7')
end subroutine MPI_Allreduce_Cmplx_7
