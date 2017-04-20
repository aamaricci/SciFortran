subroutine d_jacobi(a,d,v,nrot)
  real(8),dimension(:,:),intent(inout)               :: a
  real(8),dimension(size(a,1)),intent(out)           :: d
  real(8),dimension(size(a,1),size(a,2)),intent(out) :: v
  integer,intent(out)                                :: nrot
  integer                                            :: i,ip,iq,n
  real(8)                                            :: c,g,h,s,sm,t,tau,theta,tresh
  real(8),dimension(size(d))                         :: b,z
  !
  n=size(a,1);if(size(a,2)/=n)stop "Error in Jacobi: size(a)!=n**2 - a not a square matrix"
  !
  v = eye(n)
  forall(i=1:n)b(i)=a(i,i)
  !
  d=b
  z=0d0
  nrot=0
  !
  do i=1,50
     sm=sum(abs(a),mask=upper_triangle(n,n))
     if (sm == 0.0) return
     tresh=merge(0.2d0*sm/n**2, 0d0, i<4 )
     do ip=1,n-1
        do iq=ip+1,n
           g=100d0*abs(a(ip,iq))
           if ((i > 4) .and. (abs(d(ip))+g == abs(d(ip))) &
                .and. (abs(d(iq))+g == abs(d(iq)))) then
              a(ip,iq)=0.0
           else if (abs(a(ip,iq)) > tresh) then
              h=d(iq)-d(ip)
              if (abs(h)+g == abs(h)) then
                 t=a(ip,iq)/h
              else
                 theta=0.5d0*h/a(ip,iq)
                 t=1.0d0/(abs(theta)+sqrt(1.0d0+theta**2))
                 if (theta < 0.0) t=-t
              end if
              c=1d0/sqrt(1+t**2)
              s=t*c
              tau=s/(1d0+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq)=0d0
              call jrotate(a(1:ip-1,ip),a(1:ip-1,iq))
              call jrotate(a(ip,ip+1:iq-1),a(ip+1:iq-1,iq))
              call jrotate(a(ip,iq+1:n),a(iq,iq+1:n))
              call jrotate(v(:,ip),v(:,iq))
              nrot=nrot+1
           end if
        end do
     end do
     b(:)=b(:)+z(:)
     d(:)=b(:)
     z(:)=0d0
  end do
  stop "Too many iteration in Jacobi"
contains
  subroutine jrotate(a1,a2)
    real(8),dimension(:), intent(inout) :: a1,a2
    real(8),dimension(size(a1))         :: wk1
    wk1(:)=a1(:)
    a1(:)=a1(:)-s*(a2(:)+a1(:)*tau)
    a2(:)=a2(:)+s*(wk1(:)-a2(:)*tau)
  end subroutine jrotate
end subroutine d_jacobi

subroutine c_jacobi(A,D,U,sweep)
  implicit none
  complex(8),dimension(:,:),intent(inout)               :: a
  real(8),dimension(size(a,1)),intent(out)              :: d
  complex(8),dimension(size(a,1),size(a,2)),intent(out) :: U
  integer                                               :: n
  integer                                               :: p, q, j
  real(8)                                               :: red, off, thresh
  real(8)                                               :: t, delta, invc, s
  complex(8)                                            :: x, y, Apq
  real(8)                                               :: ev(2,size(a,1))
  integer                                               :: sweep
  real(8)                                               :: SYM_EPS=tiny(1d0)
  !
  n=size(a,1);if(size(a,2)/=n)stop "Error in Jacobi: size(a)!=n**2 - a not a square matrix"
  !
  do p = 1, n
     ev(1,p) = 0d0
     ev(2,p) = dble(A(p,p))
     d(p) = ev(2,p)
  enddo
  U = eye(n)
  red = 0.04d0/n**4
  do sweep = 1, 50
     off = 0
     do q = 2, n
        do p = 1, q - 1
           off = off + abs(A(p,q))**2
        enddo
     enddo
     if( .not. off .gt. SYM_EPS ) return !goto 10
     thresh = 0
     if( sweep .lt. 4 ) thresh = off*red
     do q = 2, n
        do p = 1, q - 1
           Apq = A(p,q)
           off = abs(Apq)**2
           if( sweep .gt. 4 .and. off .lt. SYM_EPS*(ev(2,p)**2 + ev(2,q)**2) ) then
              A(p,q) = 0
           else if( off .gt. thresh ) then
              t = .5D0*(ev(2,p) - ev(2,q))
              t = 1/(t + sign(sqrt(t**2 + off), t))
              delta = t*off
              ev(1,p) = ev(1,p) + delta
              ev(2,p) = d(p) + ev(1,p)
              ev(1,q) = ev(1,q) - delta
              ev(2,q) = d(q) + ev(1,q)
              invc = sqrt(delta*t + 1)
              s = t/invc
              t = delta/(invc + 1)
              do j = 1, p - 1
                 x = A(j,p)
                 y = A(j,q)
                 A(j,p) = x + s*(conjg(Apq)*y - t*x)
                 A(j,q) = y - s*(Apq*x + t*y)
              enddo
              do j = p + 1, q - 1
                 x = A(p,j)
                 y = A(j,q)
                 A(p,j) = x + s*(Apq*conjg(y) - t*x)
                 A(j,q) = y - s*(Apq*conjg(x) + t*y)
              enddo
              do j = q + 1, n
                 x = A(p,j)
                 y = A(q,j)
                 A(p,j) = x + s*(Apq*y - t*x)
                 A(q,j) = y - s*(conjg(Apq)*x + t*y)
              enddo
              A(p,q) = 0d0
              do j = 1, n
                 x = U(p,j)
                 y = U(q,j)
                 ! #if UCOLS
                 ! U(p,j) = x + s*(conjg(Apq)*y - t*x)
                 ! U(q,j) = y - s*(Apq*x + t*y)
                 ! #else
                 U(p,j) = x + s*(Apq*y - t*x)
                 U(q,j) = y - s*(conjg(Apq)*x + t*y)
                 ! #endif
              enddo
           endif
        enddo
     enddo
     do p = 1, n
        ev(1,p) = 0
        d(p) = ev(2,p)
     enddo
  enddo
  stop "Too many iteration in Jacobi"
  ! 10  continue
  !     ! if( sort .eq. 0 ) return
  !     do p = 1, n - 1
  !        j = p
  !        t = d(p)
  !        do q = p + 1, n
  !           if( 1*(t - d(q)) .gt. 0 ) then
  !              j = q
  !              t = d(q)
  !           endif
  !        enddo
  !        if( j .ne. p ) then
  !           d(j) = d(p)
  !           d(p) = t
  !           do q = 1, n
  !              x = U(p,q)
  !              U(p,q) = U(j,q)
  !              U(j,q) = x
  !           enddo
  !        endif
  !     enddo
end subroutine C_jacobi
