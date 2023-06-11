!+-----------------------------------------------------------------+
!PURPOSE  : Build the BETHE Lattice structure of the problem
!+-----------------------------------------------------------------+
subroutine bethe_lattice(dos,ome,Lk,D)
  real(8)          :: dos(Lk),ome(Lk)
  integer          :: Lk
  real(8),optional :: D
  integer          :: ie
  real(8)          :: de,e,D_
  complex(8)       :: gf,zeta
  D_=1.d0;if(present(D))D_=D
  de= 2.d0*D/dble(Lk-1)
  write(*,"(A,I8,A)")"Bethe Lattice with:",Lk," e-points"
  open(10,file="DOSbethe.lattice")
  do ie=1,Lk
     e=-D_ + dble(ie-1)*de
     dos(ie)=dens_bethe(e,D_)*de
     ome(ie)=e
     write(10,*)e,dos(ie)/de
  enddo
  close(10)
end subroutine bethe_lattice



!+-------------------------------------------------------------------+
!purpose  : calculate the non-interacting dos for BETHE lattice 
!+-------------------------------------------------------------------+
elemental function dens_bethe(x,D)
  real(8),intent(in)          :: x
  real(8),intent(in),optional :: d
  real(8)                     :: dens_bethe,d_
  complex(8)                  :: root,d2
  d_=1.d0;if(present(d))d_=d
  d2=dcmplx(d_,0.d0)
  root=dcmplx((1.d0-1.d0*((x/d_))**2),0.d0)
  root=sqrt(root)
  dens_bethe=(2.d0/(3.141592653589793238d0*d_))*root
end function dens_bethe




!+------------------------------------------------------------------+
!purpose  : get the hilber transfom of a given "zeta" with bethe dos
!+------------------------------------------------------------------+
elemental function gfbethe(w,zeta,d)
  real(8),intent(in)    :: w,d
  complex(8),intent(in) :: zeta
  complex(8)            :: gfbethe,sqroot
  real(8)               :: sq,sig
  sqroot=sqrt(zeta**2-d**2)
  sq=dimag(sqroot)
  sig=w*sq/abs(w*sq)
  gfbethe=2.d0/(zeta+sig*sqroot)
  return
end function gfbethe


!+------------------------------------------------------------------+
!purpose  : get the hilber transfom of a given "zeta" with bethe dos
!+------------------------------------------------------------------+
function gfbether(w,zeta,d)
  real(8)               :: w,d
  complex(8)            :: zeta
  complex(8)            :: gfbether,sqroot
  real(8)               :: sig
  if(dreal(zeta)==0.d0)zeta=dcmplx(1.d-8,dimag(zeta))
  sqroot=sqrt(zeta**2-d**2)
  sig=dreal(zeta)/abs(dreal(zeta))
  gfbether=2.d0/(zeta+sig*sqroot)
end function gfbether




!+-----------------------------------------------------------------+
!purpose  : build the bethe lattice structure of the problem
!+-----------------------------------------------------------------+
subroutine bethe_guess_g0(g0,d,beta,hloc) 
  complex(8),dimension(:) :: g0
  real(8)                 :: d,beta,hloc
  integer                 :: i
  real(8)                 :: wm
  complex(8)              :: zeta
  do i=1,size(g0)
     wm    = 3.141592653589793238d0/beta*(2*i-1)
     zeta  = dcmplx(0.d0,wm) - hloc
     g0(i) = gfbethe(wm,zeta,d)
  enddo
end subroutine bethe_guess_g0
