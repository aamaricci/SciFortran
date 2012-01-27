!###############################################################
! PROGRAM  : UTILS
! TYPE     : Module
! PURPOSE  : 
!###############################################################
MODULE BROYDN_UTILS
  implicit none
  INTEGER, PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8
  INTERFACE get_diag
     MODULE PROCEDURE get_diag_rv, get_diag_dv
  END INTERFACE get_diag
  INTERFACE outerdiff
     MODULE PROCEDURE outerdiff_r,outerdiff_d,outerdiff_i
  END INTERFACE outerdiff
  INTERFACE outerprod
     MODULE PROCEDURE outerprod_r,outerprod_d
  END INTERFACE outerprod
  INTERFACE put_diag
     MODULE PROCEDURE put_diag_rv, put_diag_r
  END INTERFACE put_diag
  INTERFACE assert_eq
     MODULE PROCEDURE assert_eq2,assert_eq3,assert_eq4,assert_eqn
  END INTERFACE assert_eq
contains

  FUNCTION get_diag_rv(mat)
    REAL(4),DIMENSION(:,:),INTENT(IN) :: mat
    REAL(4), DIMENSION(size(mat,1)) :: get_diag_rv
    INTEGER :: j
    j=assert_eq2(size(mat,1),size(mat,2),'get_diag_rv')
    do j=1,size(mat,1)
       get_diag_rv(j)=mat(j,j)
    end do
  END FUNCTION get_diag_rv
  !BL
  FUNCTION get_diag_dv(mat)
    REAL(8),DIMENSION(:,:),INTENT(IN) :: mat
    REAL(8), DIMENSION(size(mat,1)) :: get_diag_dv
    INTEGER :: j
    j=assert_eq2(size(mat,1),size(mat,2),'get_diag_dv')
    do j=1,size(mat,1)
       get_diag_dv(j)=mat(j,j)
    end do
  END FUNCTION get_diag_dv
  !
  FUNCTION assert_eq2(n1,n2,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2
    INTEGER :: assert_eq2
    if (n1 == n2) then
       assert_eq2=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq2'
    end if
  END FUNCTION assert_eq2
  !BL
  FUNCTION assert_eq3(n1,n2,n3,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2,n3
    INTEGER :: assert_eq3
    if (n1 == n2 .and. n2 == n3) then
       assert_eq3=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq3'
    end if
  END FUNCTION assert_eq3
  !BL
  FUNCTION assert_eq4(n1,n2,n3,n4,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2,n3,n4
    INTEGER :: assert_eq4
    if (n1 == n2 .and. n2 == n3 .and. n3 == n4) then
       assert_eq4=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq4'
    end if
  END FUNCTION assert_eq4
  !BL
  FUNCTION assert_eqn(nn,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, DIMENSION(:), INTENT(IN) :: nn
    INTEGER :: assert_eqn
    if (all(nn(2:) == nn(1))) then
       assert_eqn=nn(1)
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eqn'
    end if
  END FUNCTION assert_eqn


  FUNCTION lower_triangle(j,k,extra)
    INTEGER, INTENT(IN) :: j,k
    INTEGER, OPTIONAL, INTENT(IN) :: extra
    LOGICAL, DIMENSION(j,k) :: lower_triangle
    INTEGER :: n
    n=0;if (present(extra)) n=extra
    lower_triangle=(outerdiff(arth_i(1,1,j),arth_i(1,1,k)) > -n)
  contains
    function arth_i(first,increment,n)
      INTEGER, INTENT(IN) :: first,increment,n
      INTEGER, DIMENSION(n) :: arth_i
      INTEGER :: k,k2,temp
      if (n > 0) arth_i(1)=first
      if (n <= NPAR_ARTH) then
         do k=2,n
            arth_i(k)=arth_i(k-1)+increment
         end do
      else
         do k=2,NPAR2_ARTH
            arth_i(k)=arth_i(k-1)+increment
         end do
         temp=increment*NPAR2_ARTH
         k=NPAR2_ARTH
         do
            if (k >= n) exit
            k2=k+k
            arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
            temp=temp+temp
            k=k2
         end do
      end if
    END FUNCTION arth_i
  END FUNCTION lower_triangle


  FUNCTION outerdiff_r(a,b)
    REAL(4), DIMENSION(:), INTENT(IN) :: a,b
    REAL(4), DIMENSION(size(a),size(b)) :: outerdiff_r
    outerdiff_r = spread(a,dim=2,ncopies=size(b)) - spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerdiff_r
  !BL
  FUNCTION outerdiff_d(a,b)
    REAL(8), DIMENSION(:), INTENT(IN) :: a,b
    REAL(8), DIMENSION(size(a),size(b)) :: outerdiff_d
    outerdiff_d = spread(a,dim=2,ncopies=size(b)) - spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerdiff_d
  !BL
  FUNCTION outerdiff_i(a,b)
    INTEGER, DIMENSION(:), INTENT(IN) :: a,b
    INTEGER, DIMENSION(size(a),size(b)) :: outerdiff_i
    outerdiff_i = spread(a,dim=2,ncopies=size(b)) - spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerdiff_i

  FUNCTION outerprod_r(a,b)
    REAL(4), DIMENSION(:), INTENT(IN) :: a,b
    REAL(4), DIMENSION(size(a),size(b)) :: outerprod_r
    outerprod_r = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerprod_r
  !BL
  FUNCTION outerprod_d(a,b)
    REAL(8), DIMENSION(:), INTENT(IN) :: a,b
    REAL(8), DIMENSION(size(a),size(b)) :: outerprod_d
    outerprod_d = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerprod_d

  SUBROUTINE put_diag_rv(diagv,mat)
    REAL(8), DIMENSION(:), INTENT(IN) :: diagv
    REAL(8), DIMENSION(:,:), INTENT(INOUT) :: mat
    INTEGER :: j,n
    n=assert_eq2(size(diagv),min(size(mat,1),size(mat,2)),'put_diag_rv')
    do j=1,n
       mat(j,j)=diagv(j)
    end do
  END SUBROUTINE put_diag_rv
  !BL
  SUBROUTINE put_diag_r(scal,mat)
    REAL(8), INTENT(IN) :: scal
    REAL(8), DIMENSION(:,:), INTENT(INOUT) :: mat
    INTEGER :: j,n
    n = min(size(mat,1),size(mat,2))
    do j=1,n
       mat(j,j)=scal
    end do
  END SUBROUTINE put_diag_r

  SUBROUTINE unit_matrix(mat)
    REAL(8), DIMENSION(:,:), INTENT(OUT) :: mat
    INTEGER :: i,n
    n=min(size(mat,1),size(mat,2))
    mat(:,:)=0.0
    do i=1,n
       mat(i,i)=1.0
    end do
  END SUBROUTINE unit_matrix

  FUNCTION vabs(v)
    REAL(8), DIMENSION(:), INTENT(IN) :: v
    REAL(8) :: vabs
    vabs=sqrt(dot_product(v,v))
  END FUNCTION vabs

  FUNCTION ifirstloc(mask)
    LOGICAL, DIMENSION(:), INTENT(IN) :: mask
    INTEGER               :: ifirstloc
    INTEGER, DIMENSION(1) :: loc
    loc=maxloc(merge(1,0,mask))
    ifirstloc=loc(1)
    if (.not. mask(ifirstloc)) ifirstloc=size(mask)+1
  END FUNCTION ifirstloc
END MODULE BROYDN_UTILS

