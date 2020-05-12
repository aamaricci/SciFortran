!!Bool
subroutine MPI_Bcast_Bool_0(data,root,MpiComm)
  logical,intent(in)          :: data
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,1,MPI_LOGICAL,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Bool_0')
#else
  return
#endif
end subroutine MPI_Bcast_Bool_0
!
subroutine MPI_Bcast_Bool_1(data,root,MpiComm)
  logical,intent(in)          :: data(:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_LOGICAL,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Bool_1')
#else
  return
#endif
end subroutine  MPI_Bcast_Bool_1
!
subroutine MPI_Bcast_Bool_2(data,root,MpiComm)  
  logical,intent(in)          :: data(:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_LOGICAL,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Bool_2')
#else
  return
#endif
end subroutine  MPI_Bcast_Bool_2
!
subroutine MPI_Bcast_Bool_3(data,root,MpiComm)  
  logical,intent(in)          :: data(:,:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_LOGICAL,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Bool_3')
#else
  return
#endif
end subroutine  MPI_Bcast_Bool_3
!
subroutine MPI_Bcast_Bool_4(data,root,MpiComm)  
  logical,intent(in)          :: data(:,:,:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_LOGICAL,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Bool_4')
#else
  return
#endif
end subroutine  MPI_Bcast_Bool_4
!
subroutine MPI_Bcast_Bool_5(data,root,MpiComm)  
  logical,intent(in)          :: data(:,:,:,:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_LOGICAL,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Bool_5')
#else
  return
#endif
end subroutine  MPI_Bcast_Bool_5
!
subroutine MPI_Bcast_Bool_6(data,root,MpiComm)  
  logical,intent(in)          :: data(:,:,:,:,:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_LOGICAL,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Bool_6')
#else
  return
#endif
end subroutine  MPI_Bcast_Bool_6
!
subroutine MPI_Bcast_Bool_7(data,root,MpiComm)  
  logical,intent(in)          :: data(:,:,:,:,:,:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_LOGICAL,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Bool_7')
#else
  return
#endif
end subroutine  MPI_Bcast_Bool_7




!! INTEGER
subroutine MPI_Bcast_Int_0(data,root,MpiComm)  
  integer,intent(in)          :: data
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,1,MPI_INTEGER,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Int_0')
#else
  return
#endif
end subroutine  MPI_Bcast_Int_0
!
subroutine MPI_Bcast_Int_1(data,root,MpiComm)  
  integer,intent(in)          :: data(:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_INTEGER,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Int_1')
#else
  return
#endif
end subroutine  MPI_Bcast_Int_1
!
subroutine MPI_Bcast_Int_2(data,root,MpiComm)  
  integer,intent(in)          :: data(:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_INTEGER,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Int_2')
#else
  return
#endif
end subroutine  MPI_Bcast_Int_2
!
subroutine MPI_Bcast_Int_3(data,root,MpiComm)
  integer,intent(in)          :: data(:,:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_INTEGER,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Int_3')
#else
  return
#endif
end subroutine  MPI_Bcast_Int_3
!
subroutine MPI_Bcast_Int_4(data,root,MpiComm)
  integer,intent(in)          :: data(:,:,:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_INTEGER,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Int_4')
#else
  return
#endif
end subroutine  MPI_Bcast_Int_4
!
subroutine MPI_Bcast_Int_5(data,root,MpiComm)
  integer,intent(in)          :: data(:,:,:,:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_INTEGER,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Int_5')
#else
  return
#endif
end subroutine  MPI_Bcast_Int_5
!
subroutine MPI_Bcast_Int_6(data,root,MpiComm)
  integer,intent(in)          :: data(:,:,:,:,:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_INTEGER,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Int_6')
#else
  return
#endif
end subroutine  MPI_Bcast_Int_6
!
subroutine MPI_Bcast_Int_7(data,root,MpiComm)
  integer,intent(in)          :: data(:,:,:,:,:,:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_INTEGER,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Int_7')
#else
  return
#endif
end subroutine  MPI_Bcast_Int_7



!! REAL8
subroutine MPI_Bcast_Dble_0(data,root,MpiComm)
  real(8),intent(in)          :: data
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,1,MPI_DOUBLE_PRECISION,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Dble_0')
#else
  return
#endif
end subroutine  MPI_Bcast_Dble_0
!
subroutine MPI_Bcast_Dble_1(data,root,MpiComm)
  real(8),intent(in)          :: data(:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_DOUBLE_PRECISION,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Dble_1')
#else
  return
#endif
end subroutine  MPI_Bcast_Dble_1
!
subroutine MPI_Bcast_Dble_2(data,root,MpiComm)
  real(8),intent(in)          :: data(:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_DOUBLE_PRECISION,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Dble_2')
#else
  return
#endif
end subroutine  MPI_Bcast_Dble_2
!
subroutine MPI_Bcast_Dble_3(data,root,MpiComm)
  real(8),intent(in)          :: data(:,:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_DOUBLE_PRECISION,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Dble_3')
#else
  return
#endif
end subroutine  MPI_Bcast_Dble_3
!
subroutine MPI_Bcast_Dble_4(data,root,MpiComm)
  real(8),intent(in)          :: data(:,:,:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_DOUBLE_PRECISION,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Dble_4')
#else
  return
#endif
end subroutine  MPI_Bcast_Dble_4
!
subroutine MPI_Bcast_Dble_5(data,root,MpiComm)
  real(8),intent(in)          :: data(:,:,:,:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_DOUBLE_PRECISION,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Dble_5')
#else
  return
#endif
end subroutine  MPI_Bcast_Dble_5
!
subroutine MPI_Bcast_Dble_6(data,root,MpiComm)
  real(8),intent(in)          :: data(:,:,:,:,:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_DOUBLE_PRECISION,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Dble_6')
#else
  return
#endif
end subroutine  MPI_Bcast_Dble_6
!
subroutine MPI_Bcast_Dble_7(data,root,MpiComm)
  real(8),intent(in)          :: data(:,:,:,:,:,:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_DOUBLE_PRECISION,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Dble_7')
#else
  return
#endif
end subroutine  MPI_Bcast_Dble_7




!!CMPLX8
subroutine MPI_Bcast_Cmplx_0(data,root,MpiComm)
  complex(8),intent(in)       :: data
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,1,MPI_DOUBLE_COMPLEX,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Cmplx_0')
#else
  return
#endif
end subroutine  MPI_Bcast_Cmplx_0
!
subroutine MPI_Bcast_Cmplx_1(data,root,MpiComm)
  complex(8),intent(in)       :: data(:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_DOUBLE_COMPLEX,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Cmplx_1')
#else
  return
#endif
end subroutine  MPI_Bcast_Cmplx_1
!
subroutine MPI_Bcast_Cmplx_2(data,root,MpiComm)
  complex(8),intent(in)       :: data(:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_DOUBLE_COMPLEX,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Cmplx_2')
#else
  return
#endif
end subroutine  MPI_Bcast_Cmplx_2
!
subroutine MPI_Bcast_Cmplx_3(data,root,MpiComm)
  complex(8),intent(in)       :: data(:,:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_DOUBLE_COMPLEX,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Cmplx_3')
#else
  return
#endif
end subroutine  MPI_Bcast_Cmplx_3
!
subroutine MPI_Bcast_Cmplx_4(data,root,MpiComm)
  complex(8),intent(in)       :: data(:,:,:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_DOUBLE_COMPLEX,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Cmplx_4')
#else
  return
#endif
end subroutine  MPI_Bcast_Cmplx_4
!
subroutine MPI_Bcast_Cmplx_5(data,root,MpiComm)
  complex(8),intent(in)       :: data(:,:,:,:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_DOUBLE_COMPLEX,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Cmplx_5')
#else
  return
#endif
end subroutine  MPI_Bcast_Cmplx_5
!
subroutine MPI_Bcast_Cmplx_6(data,root,MpiComm)
  complex(8),intent(in)       :: data(:,:,:,:,:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_DOUBLE_COMPLEX,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Cmplx_6')
#else
  return
#endif
end subroutine  MPI_Bcast_Cmplx_6
!
subroutine MPI_Bcast_Cmplx_7(data,root,MpiComm)
  complex(8),intent(in)       :: data(:,:,:,:,:,:,:)
  integer,intent(in),optional :: root
  integer,intent(in),optional :: MpiComm
  integer                     :: Comm
#ifdef _MPI
  comm=MPI_COMM_WORLD;if(present(MpiComm))comm=MpiComm
  rank=0;if(present(root))rank=root
  if(comm==MPI_COMM_NULL)return
  call MPI_BCAST(data,size(data),MPI_DOUBLE_COMPLEX,rank,comm,ierr)
  call Error_MPI(sub='MPI_Bcast_Cmplx_7')
#else
  return
#endif
end subroutine  MPI_Bcast_Cmplx_7
