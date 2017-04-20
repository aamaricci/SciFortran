!-------------------------------------------------------------------------------------------
!PURPOSE  : Invert a hermitian n*n matrix using MKL_LAPACK library  
! M is destroyed and replaces by its inverse M^-1
!-------------------------------------------------------------------------------------------
subroutine zinv_her(A,uplo)
  complex(8),dimension(:,:)                 :: A
  complex(8)                                :: c
  character(len=*),optional                 :: uplo
  character(len=1)                          :: uplo_
  integer                                   :: n,lda,info,lwork,i,j
  integer,dimension(:),allocatable          :: ipvt
  complex(8),dimension(:),allocatable       :: work
  complex(8),dimension(1)                   :: lwork_guess
  logical                                   :: bool
  uplo_="U";if(present(uplo))uplo_=uplo
  lda=max(1,size(A,1))
  n  =size(A,2)
  allocate(ipvt(n))
  !Test hermiticity:
  bool=.false.
  testH: do i=1,size(A,1)
     do j=1,size(A,2)
        c=A(i,j)-conjg(A(j,i))
        if(c/=cmplx(0.d0,0.d0,8))bool=.true.
        exit testH
     enddo
  enddo testH
  if(bool)stop "Error MATRIX/Z_mat_invertHER: A not Hermitian"
  !
  call zhetrf(uplo_,n,A,lda,ipvt,lwork_guess,-1,info)
  if(info/=0)stop "Error MATRIX/Z_mat_invertHER: 1st call zhetrf"
  lwork=lwork_guess(1) ; allocate(work(lwork))
  call zhetrf(uplo_,n,A,lda,ipvt,work,lwork,info)
  if(info/=0)stop "Error MATRIX/Z_mat_invertHERE: 2nd call zhetrf"
  call zhetri(uplo_,n,A,lda,ipvt,work,info)
  if(info/=0)stop "Error MATRIX/Z_mat_invertHERE: zhetri"
  deallocate(ipvt,work)
  if(uplo_=="U")then
     forall(i=1:size(A,1),j=1:size(A,2),i>j)A(i,j)=conjg(A(j,i))
  elseif(uplo_=="L")then
     forall(i=1:size(A,1),j=1:size(A,2),i<j)A(i,j)=conjg(A(j,i))
  endif
end subroutine Zinv_her
