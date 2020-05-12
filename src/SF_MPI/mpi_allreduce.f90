!!BOOL
subroutine MPI_Allreduce_Bool_0(send,data,root,MpiComm)  
  logical,intent(inout)       :: data
  logical,intent(in)          :: send
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
  #ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLREDUCE(send,data,1,MPI_LOGICAL,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Bool_0')
#else
  return
#endif
end subroutine  MPI_Allreduce_Bool_0
!
subroutine MPI_Allreduce_Bool_1(send,data,root,MpiComm)  
  logical,intent(inout)       :: data(:)
  logical,intent(in)          :: send(:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLREDUCE(send,data,size(data),MPI_LOGICAL,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Bool_1')
#else
  return
#endif
end subroutine  MPI_Allreduce_Bool_1
!
subroutine MPI_Allreduce_Bool_2(send,data,root,MpiComm)  
  logical,intent(inout)       :: data(:,:)
  logical,intent(in)          :: send(:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLREDUCE(send,data,size(data),MPI_LOGICAL,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Bool_2')
#else
  return
#endif
end subroutine  MPI_Allreduce_Bool_2
!
subroutine MPI_Allreduce_Bool_3(send,data,root,MpiComm)  
  logical,intent(inout)       :: data(:,:,:)
  logical,intent(in)          :: send(:,:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLREDUCE(send,data,size(data),MPI_LOGICAL,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Bool_3')
#else
  return
#endif
end subroutine  MPI_Allreduce_Bool_3
!
subroutine MPI_Allreduce_Bool_4(send,data,root,MpiComm)  
  logical,intent(inout)       :: data(:,:,:,:)
  logical,intent(in)          :: send(:,:,:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLREDUCE(send,data,size(data),MPI_LOGICAL,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Bool_4')
#else
  return
#endif
end subroutine  MPI_Allreduce_Bool_4
!
subroutine MPI_Allreduce_Bool_5(send,data,root,MpiComm)  
  logical,intent(inout)       :: data(:,:,:,:,:)
  logical,intent(in)          :: send(:,:,:,:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLREDUCE(send,data,size(data),MPI_LOGICAL,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Bool_5')
#else
  return
#endif
end subroutine  MPI_Allreduce_Bool_5
!
subroutine MPI_Allreduce_Bool_6(send,data,root,MpiComm)  
  logical,intent(inout)       :: data(:,:,:,:,:,:)
  logical,intent(in)          :: send(:,:,:,:,:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLREDUCE(send,data,size(data),MPI_LOGICAL,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Bool_6')
#else
  return
#endif
end subroutine  MPI_Allreduce_Bool_6
!
subroutine MPI_Allreduce_Bool_7(send,data,root,MpiComm)  
  logical,intent(inout)       :: data(:,:,:,:,:,:,:)
  logical,intent(in)          :: send(:,:,:,:,:,:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLREDUCE(send,data,size(data),MPI_LOGICAL,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Bool_7')
#else
  return
#endif
end subroutine  MPI_Allreduce_Bool_7



!!INTEGER
subroutine MPI_Allreduce_Int_0(send,data,root,MpiComm)  
  integer,intent(inout)       :: data
  integer,intent(in)          :: send
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLREDUCE(send,data,1,MPI_INTEGER,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Int_0')
#else
  return
#endif
end subroutine  MPI_Allreduce_Int_0
!
subroutine MPI_Allreduce_Int_1(send,data,root,MpiComm)  
  integer,intent(inout)       :: data(:)
  integer,intent(in)          :: send(:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLREDUCE(send,data,size(data),MPI_INTEGER,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Int_1')
#else
  return
#endif
end subroutine  MPI_Allreduce_Int_1
!
subroutine MPI_Allreduce_Int_2(send,data,root,MpiComm)  
  integer,intent(inout)       :: data(:,:)
  integer,intent(in)          :: send(:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLREDUCE(send,data,size(data),MPI_INTEGER,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Int_2')
#else
  return
#endif
end subroutine  MPI_Allreduce_Int_2
!
subroutine MPI_Allreduce_Int_3(send,data,root,MpiComm)  
  integer,intent(inout)       :: data(:,:,:)
  integer,intent(in)          :: send(:,:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLREDUCE(send,data,size(data),MPI_INTEGER,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Int_3')
#else
  return
#endif
end subroutine  MPI_Allreduce_Int_3
!
subroutine MPI_Allreduce_Int_4(send,data,root,MpiComm)  
  integer,intent(inout)       :: data(:,:,:,:)
  integer,intent(in)          :: send(:,:,:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLREDUCE(send,data,size(data),MPI_INTEGER,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Int_4')
#else
  return
#endif
end subroutine  MPI_Allreduce_Int_4
!
subroutine MPI_Allreduce_Int_5(send,data,root,MpiComm)  
  integer,intent(inout)       :: data(:,:,:,:,:)
  integer,intent(in)          :: send(:,:,:,:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLREDUCE(send,data,size(data),MPI_INTEGER,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Int_5')
#else
  return
#endif
end subroutine  MPI_Allreduce_Int_5
!
subroutine MPI_Allreduce_Int_6(send,data,root,MpiComm)  
  integer,intent(inout)       :: data(:,:,:,:,:,:)
  integer,intent(in)          :: send(:,:,:,:,:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLREDUCE(send,data,size(data),MPI_INTEGER,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Int_6')
#else
  return
#endif
end subroutine  MPI_Allreduce_Int_6
!
subroutine MPI_Allreduce_Int_7(send,data,root,MpiComm)  
  integer,intent(inout)       :: data(:,:,:,:,:,:,:)
  integer,intent(in)          :: send(:,:,:,:,:,:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLREDUCE(send,data,size(data),MPI_INTEGER,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Int_7')
#else
  return
#endif
end subroutine  MPI_Allreduce_Int_7



!!REAL8
subroutine MPI_Allreduce_Dble_0(send,data,root,MpiComm)  
  real(8),intent(inout)       :: data
  real(8),intent(in)          :: send
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLREDUCE(send,data,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Dble_0')
#else
  return
#endif
end subroutine  MPI_Allreduce_Dble_0
!
subroutine MPI_Allreduce_Dble_1(send,data,root,MpiComm)  
  real(8),intent(inout)       :: data(:)
  real(8),intent(in)          :: send(:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Dble_1')
#else
  return
#endif
end subroutine  MPI_Allreduce_Dble_1
!
subroutine MPI_Allreduce_Dble_2(send,data,root,MpiComm)  
  real(8),intent(inout)       :: data(:,:)
  real(8),intent(in)          :: send(:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Dble_2')
#else
  return
#endif
end subroutine  MPI_Allreduce_Dble_2
!
subroutine MPI_Allreduce_Dble_3(send,data,root,MpiComm)  
  real(8),intent(inout)       :: data(:,:,:)
  real(8),intent(in)          :: send(:,:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Dble_3')
#else
  return
#endif
end subroutine  MPI_Allreduce_Dble_3
!
subroutine MPI_Allreduce_Dble_4(send,data,root,MpiComm)  
  real(8),intent(inout)       :: data(:,:,:,:)
  real(8),intent(in)          :: send(:,:,:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Dble_4')
#else
  return
#endif
end subroutine  MPI_Allreduce_Dble_4
!
subroutine MPI_Allreduce_Dble_5(send,data,root,MpiComm)  
  real(8),intent(inout)       :: data(:,:,:,:,:)
  real(8),intent(in)          :: send(:,:,:,:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Dble_5')
#else
  return
#endif
end subroutine  MPI_Allreduce_Dble_5
!
subroutine MPI_Allreduce_Dble_6(send,data,root,MpiComm)  
  real(8),intent(inout)       :: data(:,:,:,:,:,:)
  real(8),intent(in)          :: send(:,:,:,:,:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Dble_6')
#else
  return
#endif
end subroutine  MPI_Allreduce_Dble_6
!
subroutine MPI_Allreduce_Dble_7(send,data,root,MpiComm)  
  real(8),intent(inout)       :: data(:,:,:,:,:,:,:)
  real(8),intent(in)          :: send(:,:,:,:,:,:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Dble_7')
#else
  return
#endif
end subroutine  MPI_Allreduce_Dble_7





!!CMPLX8
subroutine MPI_Allreduce_Cmplx_0(send,data,root,MpiComm)  
  complex(8),intent(inout)    :: data
  complex(8),intent(in)       :: send
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLREDUCE(send,data,1,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Cmplx_0')
#else
  return
#endif
end subroutine  MPI_Allreduce_Cmplx_0
!
subroutine MPI_Allreduce_Cmplx_1(send,data,root,MpiComm)  
  complex(8),intent(inout)    :: data(:)
  complex(8),intent(in)       :: send(:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Cmplx_1')
#else
  return
#endif
end subroutine  MPI_Allreduce_Cmplx_1
!
subroutine MPI_Allreduce_Cmplx_2(send,data,root,MpiComm)  
  complex(8),intent(inout)    :: data(:,:)
  complex(8),intent(in)       :: send(:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Cmplx_2')
#else
  return
#endif
end subroutine  MPI_Allreduce_Cmplx_2
!
subroutine MPI_Allreduce_Cmplx_3(send,data,root,MpiComm)  
  complex(8),intent(inout)    :: data(:,:,:)
  complex(8),intent(in)       :: send(:,:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Cmplx_3')
#else
  return
#endif
end subroutine  MPI_Allreduce_Cmplx_3
!
subroutine MPI_Allreduce_Cmplx_4(send,data,root,MpiComm)  
  complex(8),intent(inout)    :: data(:,:,:,:)
  complex(8),intent(in)       :: send(:,:,:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Cmplx_4')
#else
  return
#endif
end subroutine  MPI_Allreduce_Cmplx_4
!
subroutine MPI_Allreduce_Cmplx_5(send,data,root,MpiComm)  
  complex(8),intent(inout)    :: data(:,:,:,:,:)
  complex(8),intent(in)       :: send(:,:,:,:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Cmplx_5')
#else
  return
#endif
end subroutine  MPI_Allreduce_Cmplx_5
!
subroutine MPI_Allreduce_Cmplx_6(send,data,root,MpiComm)  
  complex(8),intent(inout)    :: data(:,:,:,:,:,:)
  complex(8),intent(in)       :: send(:,:,:,:,:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Cmplx_6')
#else
  return
#endif
end subroutine  MPI_Allreduce_Cmplx_6
!
subroutine MPI_Allreduce_Cmplx_7(send,data,root,MpiComm)  
  complex(8),intent(inout)    :: data(:,:,:,:,:,:,:)
  complex(8),intent(in)       :: send(:,:,:,:,:,:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_ALLREDUCE(send,data,size(data),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
  call Error_MPI(sub='MPI_Allreduce_Cmplx_7')
#else
  return
#endif
end subroutine  MPI_Allreduce_Cmplx_7
