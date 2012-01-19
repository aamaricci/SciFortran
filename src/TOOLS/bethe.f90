
!+-----------------------------------------------------------------+
!PROGRAM  : BuildLattice
!TYPE     : Subroutine
!PURPOSE  : Build the BETHE Lattice structure of the problem
!+-----------------------------------------------------------------+
subroutine bethe_lattice(wt_,epsik_,Lk,D_,eps_)
  integer    :: N,ie,Nk,Lk
  real(8)    :: de,e,pi,D,eps
  complex(8) :: gf,zeta
  real(8)    :: wt_(Lk),epsik_(Lk)
  real(8),optional    :: D_,eps_
  D=1.d0;if(present(D_))D=D_
  eps=1.d-4 ;if(present(eps_))eps=eps_
  write(*,"(A,I8,A)")"Bethe Lattice with:",Lk," e-points"
  pi=acos(-1.d0)
  de= 2.d0*D/dble(Lk-1) !energy step
  do ie=1,Lk
     e=-D + dble(ie-1)*de
     zeta=cmplx(e,eps_)
     gf=gfbether(e,zeta,D)
     wt_(ie)=-aimag(gf)/pi*de
     epsik_(ie)=e
  enddo
  call get_free_dos(epsik_,wt_,file='DOSfree.lattice',&
       wmin=2.d0*minval(epsik_),wmax=2.d0*maxval(epsik_),eps=eps_)
  call system("if [ ! -d LATTICEinfo ]; then mkdir LATTICEinfo; fi")
  call system("mv *.lattice LATTICEinfo/ 2>/dev/null")
end subroutine bethe_lattice



!***************************************************************
!***************************************************************
!***************************************************************



!+-------------------------------------------------------------------+
!purpose  : calculate the non-interacting dos for BETHE lattice 
!+-------------------------------------------------------------------+
ELEMENTAL FUNCTION dens_bethe(x,ts)
  REAL(8),intent(in) :: x
  REAL(8),intent(in),optional :: ts
  REAL(8) :: dens_bethe,D
  complex(8):: root,D2
  D=1.d0;if(present(ts))D=2.d0*ts
  D2=D*cmplx(1.d0,0.d0)
  root=cmplx((1.d0-1.d0*((x/D))**2),0.d0)
  root=sqrt(root)
  dens_bethe=(2.d0/(3.141592653589793238d0*D))*root
  return
END FUNCTION dens_bethe



!*******************************************************************
!*******************************************************************
!*******************************************************************




!+------------------------------------------------------------------+
!PROGRAM  : BETHE 
!TYPE     : function
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
!PROGRAM  : BETHE R
!TYPE     : function
!PURPOSE  : get the Hilber transfom of a given "zeta" with Bethe DOS
!+------------------------------------------------------------------+
ELEMENTAL FUNCTION gfbether(wr,zeta,d)
  real(8),intent(in) :: wr,d
  complex(8),intent(in) :: zeta
  complex(8) :: gfbether,sqroot
  real(8)    :: sq,sig,w
  w=wr;if(wr==0.d0)w=wr+1.d-5
  sqroot=sqrt(zeta**2-d**2)
  sq=aimag(sqroot)
  sig=w*sq/dabs(w*sq)
  gfbether=2.d0/(zeta+sig*sqroot)
  if(aimag(gfbether) >= 0.d0)&
       gfbether=2.d0/(zeta-sig*sqroot)
  return
end FUNCTION gfbether
!*******************************************************************
!*******************************************************************
!*******************************************************************
