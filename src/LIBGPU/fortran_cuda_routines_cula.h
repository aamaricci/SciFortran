!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
 
      !------------------------------------------!
 subroutine cula__cgetrf_(n,a,ipiv)
 implicit none
 integer    :: n
 complex(4) :: a(n,n)
 integer    :: ipiv(n),status
 integer,external :: CULA_CGETRF
   ipiv=0
#ifdef _CULA
   status= CULA_CGETRF(n,n,a,n,ipiv)
#endif
 end subroutine
       !------------------------------------------!
 subroutine cula__zgetrf_(n,a,ipiv)
 implicit none
 integer(4)       :: n
 complex(8)       :: a(n,n)
 integer(4)       :: ipiv(n),INFO,status
 integer,external :: CULA_ZGETRF,CULA_INITIALIZE
   ipiv=0
  !status=CULA_ZGETRF(n,n,a,n,ipiv)
   call ZGETRF(n,n,a,n,ipiv,status)
 end subroutine
       !------------------------------------------!
 subroutine cula__dgetrf_(n,a,ipiv)
 implicit none
 integer    :: n,status
 real(8)    :: a(n,n)
 integer(4) :: ipiv(n)
 integer,external :: CULA_DGETRF
   ipiv=0
#ifdef _CULA
   status=CULA_DGETRF(n,n,a,n,ipiv)
#endif
 end subroutine
       !------------------------------------------!
 subroutine cula__sgetrf_(n,a,ipiv)
 implicit none
 integer    :: n,status
 real(4)    :: a(n,n)
 integer(4) :: ipiv(n)
 integer,external :: CULA_SGETRF
    ipiv=0
#ifdef _CULA
    status=CULA_SGETRF(n,n,a,n,ipiv)
#endif
 end subroutine
       !------------------------------------------!
 subroutine cula__cgetri_(n,a,ipiv)
 implicit none
 integer    :: n,status
 complex(4) :: a(n,n)
 integer(4)  :: ipiv(n)
 integer,external :: CULA_CGETRI
#ifdef _CULA
    status=CULA_CGETRI(n,a,n,ipiv)
#endif
 end subroutine
       !------------------------------------------!
 subroutine cula__zgetri_(n,a,ipiv)
 implicit none
 integer          :: n,status
 complex(8)       :: a(n,n),work(2*n)
 integer(4)       :: ipiv(n)
 integer,external :: CULA_ZGETRI,CULA_INITIALIZE
    if(n/=size(a,1).or.n/=size(a,2)) stop 'error cula_zgetri_'
  ! status=CULA_ZGETRI(n,a,ipiv,status)
    call ZGETRI(n,a,n,ipiv,work,2*n,status)
 end subroutine
       !------------------------------------------!
 subroutine cula__dgetri_(n,a,ipiv)
 implicit none
 integer    :: n,status
 real(8)    :: a(n,n)
 integer    :: ipiv(n)
 integer,external :: CULA_DGETRI
#ifdef _CULA
    status=CULA_DGETRI(n,a,n,ipiv)
#endif
 end subroutine
       !------------------------------------------!
 subroutine cula__sgetri_(n,a,ipiv)
 implicit none
 integer    :: n,status
 real(4)    :: a(n,n)
 integer    :: ipiv(n)
 integer,external :: CULA_SGETRI
#ifdef _CULA
    status=CULA_SGETRI(n,a,n,ipiv)
#endif
 end subroutine
       !------------------------------------------!

 subroutine cula_svd_cs(M,N,A,values,U,V)
 implicit none
 integer          :: M,N
 complex(4)       :: A(M,N),U(M,M),V(N,N)
 real(4)          :: values(N)
 integer          :: status,i
 integer,external :: cula_cgesvd
   if(M<N) stop 'error cula_svd leading dimension should be rows'
#ifdef _CULA
   status = cula_cgesvd('A','A',M,N,A,M,values,U,M,V,N)
#endif
 end subroutine

       !------------------------------------------!

 subroutine cula_svd_c(M,N,A,values,U,V)
 implicit none
 integer          :: M,N
 complex(8)       :: A(M,N),U(M,M),V(N,N)
 real(8)          :: values(N)
 integer          :: status,i
 integer,external :: cula_zgesvd
   if(M<N) stop 'error cula_svd leading dimension should be rows'
#ifdef _CULA
   status = cula_zgesvd('A','A',M,N,A,M,values,U,M,V,N);
#endif
 end subroutine

       !------------------------------------------!

 subroutine cula_svd_s(M,N,A,values,U,V)
 implicit none
 integer          :: M,N
 real(4)          :: A(M,N),U(M,M),V(N,N)
 real(4)          :: values(N)
 integer          :: status,i
 integer,external :: cula_sgesvd
   if(M<N) stop 'error cula_svd leading dimension should be rows'
#ifdef _CULA
   status = cula_sgesvd('A','A',M,N,A,M,values,U,M,V,N);
#endif
 end subroutine

       !------------------------------------------!

 subroutine cula_svd_d(M,N,A,values,U,V)
 implicit none
 integer          :: M,N
 real(8)          :: A(M,N),U(M,M),V(N,N)
 real(8)          :: values(N)
 integer          :: status,i
 integer,external :: cula_dgesvd
   if(M<N) stop 'error cula_svd leading dimension should be rows'
#ifdef _CULA
   status = cula_dgesvd('A','A',M,N,A,M,values,U,M,V,N);
#endif
 end subroutine

       !------------------------------------------!

 subroutine cula_eigenvector_square_cs(lsize,A,values,copyA)
 implicit none
 integer          :: lsize
 complex(4)       :: A(lsize,lsize),U(size(A,1),size(A,1)),V(size(A,1),size(A,1))
 real(4)          :: values(lsize)
 integer          :: m,n,status,i
 complex(4)       :: copyA(lsize,lsize)
 integer,external :: cula_cgesvd
   m=size(A,1) ; n=size(A,2);
   if(size(values)/=m) stop 'error cula not vector with eigenvalues is wrong'
   if(m/=n) stop 'error in cula get eigenvector for square matrix, not square'
   do i=1,size(A,1)
    copyA(i,i)=copyA(i,i)+shift_diag_positive_def
   enddo
#ifdef _CULA
   status = cula_cgesvd('N','A',m,m,copyA,m,values,U,m,A,m);
#endif
   if(minval(values)<1.0) stop 'error in cula_eigenvector, matrix not positive definite (neg eigenvalues)'
   do i=1,size(A,1)
     copyA(i,i)=copyA(i,i)  -shift_diag_positive_def
    values(i)=values(i) -shift_diag_positive_def
   enddo
 end subroutine

       !------------------------------------------!

 subroutine cula_eigenvector_square_c(lsize,A,values,copyA)
 implicit none
 integer          :: lsize
 complex(8)       :: A(lsize,lsize),U(size(A,1),size(A,1)),V(size(A,1),size(A,1))
 real(8)          :: values(lsize)
 integer          :: m,n,status,i
 complex(8)       :: copyA(lsize,lsize)
 integer,external :: cula_zgesvd
   m=size(A,1) ; n=size(A,2);
   if(size(values)/=m) stop 'error cula not vector with eigenvalues is wrong'
   if(m/=n) stop 'error in cula get eigenvector for square matrix, not square'
   do i=1,size(A,1)
    copyA(i,i)=copyA(i,i)+shift_diag_positive_def
   enddo
#ifdef _CULA
   status = cula_zgesvd('N','A',m,m,copyA,m,values,U,m,A,m);
#endif
   if(minval(values)<1.0) stop 'error in cula_eigenvector, matrix not positive definite (neg eigenvalues)'
   do i=1,size(A,1)
     copyA(i,i)=copyA(i,i)  -shift_diag_positive_def
    values(i)=values(i) -shift_diag_positive_def
   enddo
 end subroutine

       !------------------------------------------!

 subroutine cula_eigenvector_square_d(lsize,A,values,copyA)
 implicit none
 integer          :: lsize
 real(8)          :: A(lsize,lsize),U(size(A,1),size(A,1)),V(size(A,1),size(A,1)),values(lsize)
 integer          :: m,n,status,i
 real(8)          :: copyA(lsize,lsize)
 integer,external :: cula_dgesvd
   m=size(A,1) ; n=size(A,2); 
   if(size(values)/=m) stop 'error cula not vector with eigenvalues is wrong'
   if(m/=n) stop 'error in cula get eigenvector for square matrix, not square'
   do i=1,size(A,1)
    copyA(i,i)=copyA(i,i)+shift_diag_positive_def
   enddo
#ifdef _CULA
   status = cula_dgesvd('N','A',m,m,copyA,m,values,U,m,A,m);
#endif
   if(minval(values)<1.0) stop 'error in cula_eigenvector, matrix not positive definite (neg eigenvalues)'
   do i=1,size(A,1)
     copyA(i,i)=copyA(i,i)  -shift_diag_positive_def
    values(i)=values(i) -shift_diag_positive_def
   enddo
 end subroutine

       !------------------------------------------!

 subroutine cula_eigenvector_square_r(lsize,A,values,copyA)
 implicit none
 integer          :: lsize
 real(4)          :: A(lsize,lsize),U(size(A,1),size(A,1)),V(size(A,1),size(A,1)),values(lsize)
 integer          :: m,n,status,i
 real(4)          :: copyA(lsize,lsize)
 integer,external :: cula_sgesvd
   m=size(A,1) ; n=size(A,2);
   if(size(values)/=m) stop 'error cula not vector with eigenvalues is wrong'
   if(m/=n) stop 'error in cula get eigenvector for square matrix, not square'
   do i=1,size(A,1)
    copyA(i,i)=copyA(i,i)+shift_diag_positive_def
   enddo
#ifdef _CULA
   status = cula_sgesvd('N','A',m,m,copyA,m,values,U,m,A,m);
#endif
   CALL cula_CHECK_STATUS(STATUS)
   if(minval(values)<1.0) stop 'error in cula_eigenvector, matrix not positive definite (neg eigenvalues)'
   do i=1,size(A,1)
     copyA(i,i)=copyA(i,i)  -shift_diag_positive_def
    values(i)=values(i) -shift_diag_positive_def
   enddo
 end subroutine

       !------------------------------------------!

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  subroutine cula_sgerf_test
    EXTERNAL   :: CULA_INITIALIZE
    EXTERNAL   :: CULA_SGEQRF
    EXTERNAL   :: CULA_SHUTDOWN
    INTEGER    :: CULA_INITIALIZE
    INTEGER    :: CULA_SGEQRF
    INTRINSIC  :: MAX,MIN
    INTEGER    :: M,N,K,STATUS
    PARAMETER(M=8192,N=8192,K=8192)
    REAL       :: A(M, N)
    REAL       :: TAU(K)
#ifdef _CULA
        WRITE(*,*) 'Initializing CULA'
        STATUS = CULA_INITIALIZE()
        CALL cula_CHECK_STATUS(STATUS)
        WRITE(*,*) 'Calling CULA_SGEQRF'
        STATUS = CULA_SGEQRF(M, N, A, M, TAU)
        CALL cula_CHECK_STATUS(STATUS)
        WRITE(*,*) 'Shutting down CULA'
        CALL CULA_SHUTDOWN()
#endif
        stop
  end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  subroutine cula_CHECK_STATUS(STATUS)
   INTEGER STATUS
   INTEGER INFO
   INTEGER CULA_GETERRORINFO
#ifdef _CULA
         IF (STATUS .NE. 0) THEN
            IF (STATUS .EQ. 6) THEN
               INFO = CULA_GETERRORINFO()
               WRITE(*,*) 'Invalid value for parameter ', INFO
            ELSE IF (STATUS .EQ. 9) THEN
               INFO = CULA_GETERRORINFO()
               WRITE(*,*) 'Runtime error (', INFO ,')'
            ELSE
               call CULA_GETSTATUSSTRING(STATUS)
            ENDIF
            STOP 1
         ENDIF
#endif
  END subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
