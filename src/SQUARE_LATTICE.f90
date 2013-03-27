  !###############################################################
  !     PURPOSE  : Setup a 2D square lattice structure. 
  !     COMMENTS : with some linear algebra can be generalized to
  ! any lattice structure. 
  !     AUTHORS  : Adriano Amaricci
  !###############################################################
  include "VECTORS.f90"
  module SQUARE_LATTICE
    USE VECTORS
    USE COMMON_VARS
    implicit none
    private

    !parameters:
    !=========================================================
    real(8),parameter,public                       :: alat=1.d0 !lattice constant

    !k-grid:
    !=========================================================
    type(vect2D),dimension(:,:),allocatable,public :: kgrid
    integer,allocatable,dimension(:,:),public      :: kindex
    integer,allocatable,dimension(:),public        :: ik2ix,ik2iy


    !MGXM path in k-grid:
    !=========================================================
    type(vect2D),dimension(:),allocatable,public   :: kgrid_MGXMpath
    real(8),allocatable,dimension(:),public        :: MGXMpath
    integer,allocatable,dimension(:),public        :: ik2ix_MGXMpath,ik2iy_MGXMpath


    !local variables: lattice weight & hopping parameters:
    !=========================================================
    real(8),save                                   :: t_hopping,t_prime


    public :: square_lattice_dispersion
    public :: square_lattice_velocity
    public :: square_lattice_dimension
    public :: square_lattice_structure
    public :: square_lattice_dispersion_array
    public :: square_lattice_MGXMpath_dimension
    public :: square_lattice_MGXMpath_structure
    public :: square_lattice_reduxGrid_dimension
    public :: square_lattice_reduxGrid_index
    public :: square_lattice_reduxGrid_dispersion_array



  contains


    !+-----------------------------------------------------------------+
    !PURPOSE  : 
    !+-----------------------------------------------------------------+
    function square_lattice_dispersion(ka,ts,tsp) result(epsk)
      type(vect2D),intent(in)      :: ka
      real(8),intent(in),optional  :: ts,tsp
      real(8)                      :: epsk
      if(present(ts))t_hopping=ts
      if(present(tsp))t_prime=tsp
      if(t_hopping==0.d0)then
         print*,"warning in +square_lattice_dispersion: t_hopping=0!"
         call sleep(1)
      endif
      epsk=-2.d0*t_hopping*(cos(ka%x) + cos(ka%y)) &
           -4.0*t_prime*cos(ka%x)*cos(ka%y)
    end function square_lattice_dispersion



    !***************************************************************
    !***************************************************************
    !***************************************************************




    !+-------------------------------------------------------------------+
    !PROGRAM  : 
    !TYPE     : Function
    !PURPOSE  : 
    !+-------------------------------------------------------------------+
    function square_lattice_velocity(kt,ts,tsp) result(velocity)
      type(vect2D),intent(in) :: kt
      real(8),optional        :: ts,tsp
      type(vect2D)            :: velocity
      if(present(ts))t_hopping=ts
      if(present(tsp))t_prime=tsp
      if(t_hopping==0.d0)then
         print*,"warning in +square_lattice_velocity: t_hopping=0!"
         call sleep(1)
      endif
      velocity%x = 2.0*t_hopping*sin(kt%x) + 4.0*t_prime*sin(kt%x)*cos(kt%y)
      velocity%y = 2.0*t_hopping*sin(kt%y) + 4.0*t_prime*cos(kt%x)*sin(kt%y)
    end function Square_lattice_velocity




    !***************************************************************
    !***************************************************************
    !***************************************************************




    !+-----------------------------------------------------------------+
    !PURPOSE  : 
    !+-----------------------------------------------------------------+
    integer function square_lattice_dimension(Nx,Ny)
      integer          :: Nx,Lk,Nk
      integer,optional :: Ny 
      Nk=Nx/2+1
      Lk=Nk*(Nk+1)/2 ; if(present(Ny))Lk=(Nx+1)*(Ny+1)
      square_lattice_dimension = Lk
    end function square_lattice_dimension



    !***************************************************************
    !***************************************************************
    !***************************************************************




    !+-----------------------------------------------------------------+
    !PURPOSE  : Build the Lattice structure for the 2D square lattice
    !+-----------------------------------------------------------------+
    function square_lattice_structure(Lk,Nx,Ny) result(wt)
      integer              :: Nx,Lk
      integer,optional     :: Ny
      real(8),dimension(Lk):: wt
      integer              :: ik,ix,iy,Nk,Lk_
      real(8)              :: peso 
      real(8)              :: Kx,Ky
      type(vect2D)         :: ai,aj,bi,bj

      Nk=Nx/2+1 ; Lk_=Nk*(Nk+1)/2 ; if(present(Ny))Lk_=(Nx+1)*(Ny+1)
      if(Lk_ /= Lk)call error("LATTICE: the +input1.dimension in +build_2dsquare_lattice is wrong!")

      if(present(Ny))then
         write(*,"(A,I8,A)")"Full BZ:",Lk," k-points"
      else
         write(*,"(A,I8,A)")"Red. BZ:",Lk," k-points"         
      endif

      ai=alat*Xver       ; aj=alat*Yver
      bi=(pi2/alat)*Xver ; bj=(pi2/alat)*Yver

      if(present(Ny))then
         call build_full_reciprocal_lattice()
      else
         call build_reduced_reciprocal_lattice() 
      endif

      call plot_reciprocal_lattice(bi,bj,kgrid)
      call system("if [ ! -d LATTICEinfo ]; then mkdir LATTICEinfo; fi")
      call system("mv *.lattice LATTICEinfo/ 2>/dev/null")

    contains

      !-------------------------
      subroutine build_full_reciprocal_lattice()
        allocate(kindex(0:Nx,0:Ny),kgrid(0:Nx,0:Ny)) 
        allocate(ik2ix(Lk),ik2iy(Lk))
        ik2ix=0 ; ik2iy=0 ; wt=0.d0
        kgrid=Vzero
        ik=0
        do ix=0,Nx
           do iy=0,Ny
              ik=ik+1
              Kx=dble(ix)/dble(Nx+1)
              Ky=dble(iy)/dble(Ny+1)

              ik2ix(ik)=ix
              ik2iy(ik)=iy
              kindex(ix,iy)=ik

              kgrid(ix,iy)=Kx*bi + Ky*bj  - pi*Vone 
           enddo
        enddo
        wt=1.d0/real(Lk,8)
      end subroutine build_full_reciprocal_lattice
      !-------------------------

      !-------------------------
      subroutine build_reduced_reciprocal_lattice()
        allocate(kindex(0:Nk,0:Nk),kgrid(0:Nk,0:Nk))
        allocate(ik2ix(Lk),ik2iy(Lk))
        ik2ix=0 ; ik2iy=0 ; wt=0.d0
        kgrid=Vzero
        ik=0
        do ix=0,Nx/2
           do iy=0,ix          
              ik=ik+1
              Kx=dble(ix)/dble(Nx)
              Ky=dble(iy)/dble(Nx)
              ik2ix(ik)=ix
              ik2iy(ik)=iy
              kgrid(ix,iy)=Kx*bi + Ky*bj - pi*Vone 
              kindex(ix,iy)=ik
              if (ix==0) then
                 peso=1.d0       !center
              elseif(ix==Nx/2) then
                 if (iy==0) then
                    peso=2.d0    ! point (pi,0)
                 elseif(iy==Nx/2) then 
                    peso=1.d0    ! corner 
                 else
                    peso=4.d0    ! border
                 endif
              else
                 if (iy==ix) then
                    peso=4.d0    ! diagonal
                 elseif (iy==0) then
                    peso=4.d0    ! x-axis
                 else
                    peso=8.d0    ! all other points
                 endif
              endif
              wt(ik)=peso/dble(Nx**2)
           enddo
        enddo
      end subroutine build_reduced_reciprocal_lattice

      subroutine plot_reciprocal_lattice(bi,bj,grid)
        integer :: i,j
        type(vect2D) :: bi,bj
        type(vect2D),dimension(0:,0:) :: grid
        open(40,file='Reciprocal.lattice')
        do i=0,size(grid,1)-1
           do j=0,size(grid,2)-1
              write(40,*) grid(i,j)%x,grid(i,j)%y
           enddo
        enddo
        write(40,*)''
        write(40,*)0.d0,0.d0
        write(40,*)''
        write(40,*)bi
        write(40,*)''
        write(40,*)bj
        write(40,*)''
        write(40,*)bi+bj
        write(40,*)''
        write(40,*)(-1.d0)*bi
        write(40,*)''
        write(40,*)(-1.d0)*bj
        write(40,*)''
        write(40,*)(-1.d0)*bi+(-1.d0)*bj
        close(40)
      end subroutine plot_reciprocal_lattice

    end function square_lattice_structure




    !***************************************************************
    !***************************************************************
    !***************************************************************




    !+-----------------------------------------------------------------+
    !PURPOSE  :   Build dispersion relation arrays:\epsilon(\ka)=epsk(i)
    !+-----------------------------------------------------------------+
    function square_lattice_dispersion_array(Lk,ts,tsp) result(epsik)!,sorted_epsik,sorted_ik)
      integer                     :: Lk
      real(8),dimension(Lk)       :: epsik
      integer                     :: ik,ix,iy
      real(8),intent(in)          :: ts
      real(8),intent(in),optional :: tsp
      ! real(8),dimension(size(epsik)),optional    :: sorted_epsik
      ! integer,dimension(size(epsik)),optional    :: sorted_ik
      t_hopping=ts
      t_prime  =0.d0 ; if(present(tsp))t_prime=tsp
      do ik=1,Lk
         ix=ik2ix(ik)
         iy=ik2iy(ik)
         epsik(ik)=square_lattice_dispersion(kgrid(ix,iy))
      enddo
      ! if(present(sorted_epsik).AND.present(sorted_ik))then
      !    sorted_epsik=epsik;call sort(sorted_epsik,sorted_ik)
      ! elseif(present(sorted_epsik))then
      !    sorted_epsik=epsik;call sort(sorted_epsik)
      ! endif
    end function square_lattice_dispersion_array



    !***************************************************************
    !***************************************************************
    !***************************************************************





    !+-----------------------------------------------------------------+
    !PURPOSE  : 
    !+-----------------------------------------------------------------+
    integer function square_lattice_MGXMpath_dimension(Nx_,Ny_)
      integer :: Nx_,Ny_
      square_lattice_MGXMpath_dimension = 2*(Nx_/2+1) + (Ny_/2+1)
    end function square_lattice_MGXMpath_dimension



    !***************************************************************
    !***************************************************************
    !***************************************************************




    !+-----------------------------------------------------------------+
    !PROGRAM  : 
    !TYPE     : Subroutine
    !PURPOSE  : 
    !+-----------------------------------------------------------------+
    subroutine square_lattice_MGXMpath_structure(Nx,Ny,epsik_,ipoint_)
      integer :: ik,ix,iy,Lpath,Nx,Ny
      real(8),dimension(:),intent(inout) :: epsik_
      integer,dimension(4),intent(inout) :: ipoint_

      !a - start from M=(pi,pi) and go to \Gamma=(0,0)
      !b - go from \Gamma(0,0) to X=(pi,0)
      !c - go from X=(pi,0) to M=(pi,pi)
      Lpath=size(epsik_)
      allocate(kgrid_MGXMpath(Lpath))
      allocate(ik2ix_MGXMpath(Lpath),ik2iy_MGXMpath(Lpath))
      allocate(MGXMpath(Lpath))
      ik=0
      ipoint_(1)=0
      do ix=Nx,Nx/2,-1
         iy=ix    !Diagonal condition
         ik=ik+1 ; kgrid_MGXMpath(ik)=kgrid(ix,iy)
         ik2ix_MGXMpath(ik)=ix ; ik2iy_MGXMpath(ik)=iy
      enddo
      ipoint_(2)=ik
      do ix=Nx/2,Nx
         iy=Ny/2   !x-axis condition
         ik=ik+1 ; kgrid_MGXMpath(ik)=kgrid(ix,iy)
         ik2ix_MGXMpath(ik)=ix ; ik2iy_MGXMpath(ik)=iy
      enddo
      ipoint_(3)=ik
      do iy=Ny/2,Ny
         ix=Nx     !y-axis condition
         ik=ik+1 ; kgrid_MGXMpath(ik)=kgrid(ix,iy)
         ik2ix_MGXMpath(ik)=ix ; ik2iy_MGXMpath(ik)=iy
      enddo
      ipoint_(4)=Lpath

      do ik=1,Lpath
         MGXMpath(ik)=dble(ik)
         epsik_(ik)=square_lattice_dispersion(kgrid_MGXMpath(ik))
      enddo
      return
    end subroutine square_lattice_MGXMpath_structure



    !***************************************************************
    !***************************************************************
    !***************************************************************




    !+-----------------------------------------------------------------+
    !PURPOSE  : 
    !+-----------------------------------------------------------------+
    subroutine square_lattice_reduxGrid_dimension(Lk_,step_,Lkredux_)
      integer               :: count,ik
      integer,intent(inout) :: Lkredux_
      integer,intent(in)    :: Lk_
      integer,intent(in)    :: step_
      if(Lk_ > Lkredux_)then
         count=0
         do ik=1,Lk_,step_
            count=count+1
         enddo
         Lkredux_=count
      else
         Lkredux_=Lk_
      endif
      return
    end subroutine square_lattice_reduxGrid_dimension



    !***************************************************************
    !***************************************************************
    !***************************************************************







    !+-----------------------------------------------------------------+
    !PURPOSE  : 
    !+-----------------------------------------------------------------+
    subroutine square_lattice_reduxGrid_index(Lk_,step_,reduced_index)
      integer,intent(in)                  :: Lk_
      integer,intent(in)                  :: step_
      integer                             :: ik,count
      integer,dimension(:),intent(inout)  :: reduced_index
      count=0
      do ik=1,Lk_,step_
         count=count+1
         reduced_index(count) = ik
      enddo
      return
    end subroutine square_lattice_reduxGrid_index




    !***************************************************************
    !***************************************************************
    !***************************************************************






    !+-----------------------------------------------------------------+
    !PURPOSE  : 
    !+-----------------------------------------------------------------+
    subroutine square_lattice_reduxGrid_dispersion_array(epsik_,reduced_index,reduced_epsik)
      real(8),dimension(:),intent(in)        :: epsik_
      integer,dimension(:),intent(in)        :: reduced_index
      real(8),dimension(size(reduced_index)) :: reduced_epsik
      integer                                :: ik,Lkreduced
      Lkreduced=size(reduced_index)
      forall(ik=1:Lkreduced)reduced_epsik(ik) = epsik_(reduced_index(ik))
      return
    end subroutine square_lattice_reduxGrid_dispersion_array



    !***************************************************************
    !***************************************************************
    !***************************************************************


    ! include "sortEPSIK.f90"


  end module SQUARE_LATTICE
