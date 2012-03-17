
!+-----------------------------------------------------------------+
!PROGRAM  : BuildLattice
!TYPE     : Subroutine
!PURPOSE  : Build the BETHE Lattice structure of the problem
!+-----------------------------------------------------------------+
subroutine bethe_lattice(dos,ome,Lk,D)
  real(8)          :: dos(Lk),ome(Lk)
  integer          :: Lk
  real(8),optional :: D
  !
  integer          :: ie
  real(8)          :: de,e,D_
  complex(8)       :: gf,zeta
  !
  D_=1.d0;if(present(D))D_=D
  de= 2.d0*D/dble(Lk-1)
  write(*,"(A,I8,A)")"Bethe Lattice with:",Lk," e-points"
  call system("if [ ! -d LATTICEinfo ]; then mkdir LATTICEinfo; fi")
  open(10,file="LATTICEinfo/DOSbethe.lattice")
  do ie=1,Lk
     e=-D + dble(ie-1)*de
     dos(ie)=dens_bethe(e,D)*de
     ome(ie)=e
     ! zeta=cmplx(e,eps_)
     ! gf=gfbether(e,zeta,D)
     ! dos(ie)=-aimag(gf)/pi*de
     write(10,*)e,dos(ie)/de
  enddo
  close(10)
end subroutine bethe_lattice



!***************************************************************
!***************************************************************
!***************************************************************



!+-------------------------------------------------------------------+
!purpose  : calculate the non-interacting dos for BETHE lattice 
!+-------------------------------------------------------------------+
ELEMENTAL FUNCTION dens_bethe(x,D)
  REAL(8),intent(in) :: x
  REAL(8),intent(in),optional :: D
  REAL(8) :: dens_bethe,D_
  complex(8):: root,D2
  D_=1.d0;if(present(D))D_=D
  D2=cmplx(D_,0.d0,8)
  root=cmplx((1.d0-1.d0*((x/D))**2),0.d0,8)
  root=sqrt(root)
  dens_bethe=(2.d0/(3.141592653589793238d0*D))*root
END FUNCTION dens_bethe



!*******************************************************************
!*******************************************************************
!*******************************************************************




!+------------------------------------------------------------------+
!PURPOSE  : get the Hilber transfom of a given "zeta" with Bethe DOS
!+------------------------------------------------------------------+
ELEMENTAL FUNCTION gfbethe(w,zeta,d)
  real(8),intent(in) :: w,d
  complex(8),intent(in) :: zeta
  complex(8) :: gfbethe,sqroot
  real(8)    :: sq,sig
  sqroot=cdsqrt(zeta**2-d**2)
  sq=aimag(sqroot)
  sig=w*sq/abs(w*sq)
  gfbethe=2.d0/(zeta+sig*sqroot)
  return
end FUNCTION gfbethe


!*******************************************************************
!*******************************************************************
!*******************************************************************


!+------------------------------------------------------------------+
!PURPOSE  : get the Hilber transfom of a given "zeta" with Bethe DOS
!+------------------------------------------------------------------+
ELEMENTAL FUNCTION gfbether(w,zeta,d)
  real(8),intent(in) :: w,d
  complex(8),intent(in) :: zeta
  complex(8) :: gfbether,sqroot
  real(8)    :: sig
  sqroot=cdsqrt(zeta**2-d**2)
  sig=real(zeta,8)/abs(real(zeta,8))
  gfbether=2.d0/(zeta+sig*sqroot)
end FUNCTION gfbether
! ELEMENTAL FUNCTION gfbether(w,zeta,d)
!   real(8),intent(in) :: w,d
!   complex(8),intent(in) :: zeta
!   complex(8) :: gfbether,sqroot
!   real(8)    :: sq,sig,x
!   x=w;if(w==0.d0)x=w+1.d-5
!   sqroot=cdsqrt(zeta**2-d**2)
!   sq=aimag(sqroot)
!   sig=w*sq/abs(w*sq)
!   gfbether=2.d0/(zeta+sig*sqroot)
!   if(aimag(gfbether) >= 0.d0)&
!        gfbether=2.d0/(zeta-sig*sqroot)
!   return
! end FUNCTION gfbether
!*******************************************************************
!*******************************************************************
!*******************************************************************
