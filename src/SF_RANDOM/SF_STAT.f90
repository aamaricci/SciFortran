MODULE SF_STAT
   USE SF_ARRAYS, only: linspace
   USE SF_INTEGRATE, only: simps
   USE SF_IOTOOLS, only: free_unit,splot3d
   USE SF_LINALG, only: det,inv
   implicit none
   private

   type,public :: histogram
      integer                       :: n=0
      real(8),dimension(:),pointer  :: range
      real(8),dimension(:),pointer  :: bin
   end type histogram

   type,public :: pdf_kernel
      integer                          :: N=0
      real(8)                          :: xmin
      real(8)                          :: xmax
      real(8)                          :: dx
      integer                          :: Ndata=0
      real(8),dimension(:),allocatable :: x
      real(8),dimension(:),allocatable :: pdf
      real(8)                          :: sigma
      logical                          :: status=.false.
      logical                          :: variance=.false.
      logical                          :: rescale=.false.
   end type pdf_kernel

   type,public :: pdf_kernel_2d
      integer,dimension(2)               :: N
      real(8),dimension(2)               :: xmin
      real(8),dimension(2)               :: xmax
      real(8),dimension(2)               :: dx
      integer                            :: Ndata=0
      real(8),dimension(:),allocatable   :: x,y
      real(8),dimension(:,:),allocatable :: pdf
      real(8),dimension(2,2)             :: Sigma
      logical                            :: status=.false.
      logical                            :: variance=.false.
      logical                            :: rescale=.false.
   end type pdf_kernel_2d


   interface pdf_allocate
      module procedure :: pdf_allocate_1d
      module procedure :: pdf_allocate_2d
   end interface pdf_allocate

   interface pdf_deallocate
      module procedure :: pdf_deallocate_1d
      module procedure :: pdf_deallocate_2d
   end interface pdf_deallocate


   interface pdf_save
      module procedure :: pdf_save_1d
      module procedure :: pdf_save_2d
   end interface pdf_save


   interface pdf_read
      module procedure :: pdf_read_1d
      module procedure :: pdf_read_2d
   end interface pdf_read


   interface pdf_set_range
      module procedure :: pdf_set_range_1d
      module procedure :: pdf_set_range_2d
   end interface pdf_set_range


   interface pdf_push_sigma
      module procedure :: pdf_push_sigma_1d
      module procedure :: pdf_push_sigma_2d
   end interface pdf_push_sigma


   interface pdf_get_sigma
      module procedure :: pdf_get_sigma_1d
      module procedure :: pdf_get_sigma_2d
   end interface pdf_get_sigma


   interface pdf_sigma
      module procedure :: pdf_sigma_data_1d
      module procedure :: pdf_sigma_sdev_1d
      module procedure :: pdf_sigma_data_2d
      module procedure :: pdf_sigma_sdev_2d
   end interface pdf_sigma


   interface pdf_accumulate
      module procedure :: pdf_accumulate_s_1d
      module procedure :: pdf_accumulate_v_1d
      module procedure :: pdf_accumulate_s_2d
   end interface pdf_accumulate



   interface pdf_normalize
      module procedure :: pdf_normalize_1d
      module procedure :: pdf_normalize_2d
   end interface pdf_normalize


   interface pdf_print
      module procedure :: pdf_print_pfile_1d
      module procedure :: pdf_print_pfile_2d
   end interface pdf_print

   interface pdf_write
      module procedure :: pdf_print_pfile_1d
      module procedure :: pdf_print_pfile_2d
   end interface pdf_write



   interface pdf_mean
      module procedure :: pdf_mean_1d
   end interface pdf_mean

   interface pdf_var
      module procedure :: pdf_var_1d
   end interface pdf_var

   interface pdf_sdev
      module procedure :: pdf_sdev_1d
   end interface pdf_sdev

   interface pdf_moment
      module procedure :: pdf_moment_1d
   end interface pdf_moment

   interface pdf_skew
      module procedure :: pdf_skew_1d
   end interface pdf_skew

   interface pdf_curt
      module procedure :: pdf_curt_1d
   end interface pdf_curt

   interface pdf_print_moments
      module procedure :: pdf_print_moments_pfile_1d
   end interface pdf_print_moments



   public :: get_moments
   public :: get_mean,get_sd,get_var,get_skew,get_curt
   public :: get_covariance

   public :: pdf_allocate
   public :: pdf_deallocate
   public :: pdf_set_range
   public :: pdf_sigma
   public :: pdf_push_sigma
   public :: pdf_accumulate
   public :: pdf_normalize
   public :: pdf_save
   public :: pdf_read
   public :: pdf_print
   public :: pdf_print_moments
   public :: pdf_mean
   public :: pdf_var
   public :: pdf_sdev
   public :: pdf_skew
   public :: pdf_curt
   public :: pdf_moment


   public :: histogram_allocate
   public :: histogram_deallocate
   public :: histogram_set_range_uniform
   public :: histogram_accumulate
   public :: histogram_get_range
   public :: histogram_get_value
   public :: histogram_print
   public :: histogram_reset


contains


   subroutine get_moments(data,ave,sdev,var,skew,curt)
      real(8), intent(out),optional     :: ave,sdev,var,skew,curt
      real(8)                           :: ave_,sdev_,var_,skew_,curt_
      real(8), dimension(:), intent(in) :: data
      integer                           :: n
      real(8)                           :: ep,adev
      real(8), dimension(size(data))    :: p,s
      n=size(data)
      if (n <= 1) then
         print*,'data.size must be at least 2'
         stop
      endif
      ave_=sum(data(:))/n
      s(:)=data(:)-ave_
      ep=sum(s(:))
      adev=sum(abs(s(:)))/n
      p(:)=s(:)*s(:)
      var_=sum(p(:))
      p(:)=p(:)*s(:)
      skew_=sum(p(:))
      p(:)=p(:)*s(:)
      curt_=sum(p(:))
      var_=(var_-ep**2/n)/(n-1)
      sdev_=sqrt(var_)
      if (var_ /= 0.0) then
         skew_=skew_/(n*sdev_**3)
         curt_=curt_/(n*var_**2)-3.0d0
      else
         skew_=0.d0
         curt_=0.d0
      endif
      if(present(ave))ave=ave_
      if(present(sdev))sdev=sdev_
      if(present(var))var=var_
      if(present(skew))skew=skew_
      if(present(curt))curt=curt_
   end subroutine get_moments


   function get_mean(data) result(mean)
      real(8), dimension(:), intent(in) :: data
      real(8)                           :: mean
      call get_moments(data,ave=mean)
   end function get_mean


   function get_sd(data) result(sd)
      real(8), dimension(:), intent(in) :: data
      real(8)                           :: sd
      call get_moments(data,sdev=sd)
   end function get_sd


   function get_var(data) result(var)
      real(8), dimension(:), intent(in) :: data
      real(8)                           :: var
      call get_moments(data,var=var)
   end function get_var


   function get_skew(data) result(skew)
      real(8), dimension(:), intent(in) :: data
      real(8)                           :: skew
      call get_moments(data,skew=skew)
   end function get_skew


   function get_curt(data) result(curt)
      real(8), dimension(:), intent(in) :: data
      real(8)                           :: curt
      call get_moments(data,curt=curt)
   end function get_curt


   function get_covariance(data,mean) result(covariance)
      real(8),dimension(:,:),intent(in)            :: data
      real(8),dimension(size(data,1)),intent(in)   :: mean
      real(8),dimension(size(data,1),size(data,1)) :: covariance
      real(8),dimension(size(data,2))              :: Xi,Xj
      integer                                      :: i,j,ncol,L
      ncol=size(data,1)
      L   =size(data,2)
      covariance=0.d0
      do i=1,ncol
         Xi(:)=data(i,:)-mean(i)
         do j=1,ncol
            Xj(:)=data(j,:)-mean(j)
            covariance(i,j) = sum(Xi(:)*Xj(:))/real(L-1,8)
         enddo
      enddo
   end function get_covariance




   !!! 1D PDF
   subroutine pdf_allocate_1d(self,N)
      type(pdf_kernel) :: self
      integer          :: N
      if(self%status)stop "PDF_ALLOCATE: PDF already allocated"
      self%N    = N
      self%Ndata= 0
      allocate(self%x(N))
      allocate(self%pdf(N))
      self%status=.true.
      self%pdf = 0d0
      self%sigma = 0d0
   end subroutine pdf_allocate_1d
   !
   subroutine pdf_deallocate_1d(self)
      type(pdf_kernel) :: self
      if(.not.self%status)stop "PDF_DEALLOCATE: PDF not allocated"
      self%xmin=0d0
      self%xmax=0d0
      self%dx  =0d0
      self%sigma = 0d0
      self%N=0
      self%Ndata= 0
      deallocate(self%x)
      deallocate(self%pdf)
      self%status=.false.
   end subroutine pdf_deallocate_1d
   !
   subroutine pdf_save_1d(self,pfile)
      type(pdf_kernel) :: self
      character(len=*) :: pfile
      integer          :: i,unit
      if(.not.self%status)stop "PDF_SAVE: PDF not allocated"
      open(free_unit(unit),file=trim(pfile))
      write(unit,*)self%N
      write(unit,*)self%xmin
      write(unit,*)self%xmax
      write(unit,*)self%dx
      write(unit,*)self%Ndata
      write(unit,*)self%x
      write(unit,*)self%pdf
      write(unit,*)self%sigma
      write(unit,*)self%status
      write(unit,*)self%variance
      write(unit,*)self%rescale
      close(unit)
   end subroutine pdf_save_1d
   !
   subroutine pdf_read_1d(self,pfile)
      type(pdf_kernel) :: self
      character(len=*) :: pfile
      integer          :: i,N,unit
      if(.not.self%status)write(*,"(A)")"PDF_READ: PDF not allocated"
      open(free_unit(unit),file=trim(pfile))
      read(unit,*)N
      call pdf_allocate(self,N)
      read(unit,*)self%xmin
      read(unit,*)self%xmax
      read(unit,*)self%dx
      read(unit,*)self%Ndata
      read(unit,*)self%x
      read(unit,*)self%pdf
      read(unit,*)self%sigma
      read(unit,*)self%status
      read(unit,*)self%variance
      read(unit,*)self%rescale
      close(unit)
   end subroutine pdf_read_1d
   !
   subroutine pdf_set_range_1d(self,a,b)
      type(pdf_kernel) :: self
      real(8)          :: a
      real(8)          :: b
      if(.not.self%status)stop "PDF_SET_RANGE: PDF not allocated"
      self%xmin=a
      self%xmax=b
      self%x = linspace(a,b,self%N,mesh=self%dx)
   end subroutine pdf_set_range_1d
   !
   subroutine pdf_push_sigma_1d(self,sigma)
      type(pdf_kernel) :: self
      real(8)          :: sigma
      if(.not.self%status)stop "PDF_SET_SIGMA: PDF not allocated"
      self%sigma = sigma
      self%variance=.true.
   end subroutine pdf_push_sigma_1d
   !
   subroutine pdf_get_sigma_1d(self,sigma)
      type(pdf_kernel) :: self
      real(8)          :: sigma
      if(.not.self%status)stop "PDF_SET_SIGMA: PDF not allocated"
      sigma = self%sigma
   end subroutine pdf_get_sigma_1d
   !
   subroutine pdf_sigma_data_1d(self,data,h)
      type(pdf_kernel)     :: self
      real(8),dimension(:) :: data
      real(8)              :: h
      if(.not.self%status)stop "PDF_SET_SIGMA: PDF not allocated"
      h = get_sd(data)
      h = (4d0/3d0/size(data))**(1d0/5d0)*h !Silverman's rule of thumb.
   end subroutine pdf_sigma_data_1d
   !
   subroutine pdf_sigma_sdev_1d(self,sdev,N,h)
      type(pdf_kernel) :: self
      real(8)          :: sdev
      integer          :: N
      real(8)          :: h
      if(.not.self%status)stop "PDF_SET_SIGMA: PDF not allocated"
      h = (4d0/3d0/N)**(1d0/5d0)*sdev !Silverman's rule of thumb.
   end subroutine pdf_sigma_sdev_1d
   !
   subroutine pdf_accumulate_s_1d(self,data,sigma)
      type(pdf_kernel) :: self
      real(8)          :: data,sigma_
      real(8),optional :: sigma
      if(.not.self%status)stop "PDF_ACCUMULATE: PDF not allocated"
      if(self%variance)then
         sigma_ = self%sigma
      else
         if(present(sigma))then
            sigma_ = sigma
         else
            stop "PDF_ACCUMULATE: PDF sigma not set or passed"
         endif
      endif
      self%pdf   = self%pdf + gaussian_kernel_1d(self%x,data,sigma_)
      self%Ndata = self%Ndata + 1
   end subroutine pdf_accumulate_s_1d
   !
   subroutine pdf_accumulate_v_1d(self,data,sigma)
      type(pdf_kernel) :: self
      real(8)          :: data(:),sigma_
      real(8),optional :: sigma
      integer          :: i
      if(.not.self%status)stop "PDF_ACCUMULATE: PDF not allocated"
      if(self%variance)then
         sigma_ = self%sigma
      else
         if(present(sigma))then
            sigma_ = sigma
         else
            stop "PDF_ACCUMULATE: PDF sigma not set or passed"
         endif
      endif
      !
      do i=1,size(data)
         self%pdf   = self%pdf + gaussian_kernel_1d(self%x,data(i),sigma_)
         self%Ndata = self%Ndata + 1
      enddo
   end subroutine pdf_accumulate_v_1d
   !
   elemental function gaussian_kernel_1d(x,mean,sigma) result(gaussian_kernel)
      real(8),intent(in) :: x,mean,sigma
      real(8)            :: gaussian_kernel
      real(8),parameter  :: pi2=2d0*acos(-1d0)
      gaussian_kernel=exp(-0.5d0*((x-mean)/sigma)**2)/sqrt(pi2)/sigma
   end function gaussian_kernel_1d
   !
   subroutine pdf_normalize_1d(self)
      type(pdf_kernel) :: self
      if(.not.self%status)stop "PDF_NORMALIZE: PDF not allocated"
      if(.not.self%rescale)self%pdf = self%pdf/self%Ndata
      self%rescale = .true.
      self%pdf = self%pdf/sum(self%pdf)/self%dx
   end subroutine pdf_normalize_1d
   !
   subroutine pdf_print_pfile_1d(self,pfile,normalize)
      type(pdf_kernel) :: self
      character(len=*) :: pfile
      logical,optional :: normalize
      logical          :: normalize_
      integer          :: i,unit
      if(.not.self%status)stop "PDF_WRITE: PDF not allocated"
      normalize_ = .true.; if(present(normalize))normalize_=normalize
      if(normalize_)call pdf_normalize_1d(self)
      open(free_unit(unit),file=trim(pfile),access='append')
      do i=1,self%N
         write(unit,*)self%x(i),self%pdf(i)
      enddo
      write(unit,*)
      close(unit)
   end subroutine pdf_print_pfile_1d
   !
   function pdf_mean_1d(self) result(mean)
      type(pdf_kernel) :: self
      real(8)          :: mean
      mean = simps(self%x*self%pdf,self%xmin,self%xmax)
   end function pdf_mean_1d
   !
   function pdf_var_1d(self) result(var)
      type(pdf_kernel) :: self
      real(8)          :: var,mean
      mean = pdf_mean_1d(self)
      var  = simps((self%x-mean)**2*self%pdf,self%xmin,self%xmax)
   end function pdf_var_1d
   !
   function pdf_sdev_1d(self) result(sdev)
      type(pdf_kernel) :: self
      real(8)          :: var,sdev
      var  = pdf_var_1d(self)
      sdev = sqrt(var)
   end function pdf_sdev_1d
   !
   function pdf_moment_1d(self,n,mu) result(mom)
      type(pdf_kernel) :: self
      integer          :: n
      real(8),optional :: mu
      real(8)          :: mom,mu_
      mu_ = 0d0;if(present(mu))mu_=mu
      mom = simps((self%x-mu_)**n*self%pdf,self%xmin,self%xmax)
   end function pdf_moment_1d
   !
   function pdf_skew_1d(self) result(skew)
      type(pdf_kernel) :: self
      real(8)          :: mean,sdev,skew
      mean = pdf_mean_1d(self)
      sdev = pdf_sdev_1d(self)
      skew = pdf_moment_1d(self,3)
      skew = skew - 3*mean*sdev**2 - mean**3
      skew = skew/sdev**3
   end function pdf_skew_1d
   !
   function pdf_curt_1d(self) result(curt)
      type(pdf_kernel) :: self
      real(8)          :: mean,var,curt
      mean = pdf_mean_1d(self)
      var  = pdf_var_1d(self)
      curt = pdf_moment_1d(self,4,mean)
      curt = curt/var**2
   end function pdf_curt_1d
   !
   subroutine pdf_print_moments_pfile_1d(self,pfile)
      type(pdf_kernel) :: self
      character(len=*) :: pfile
      integer          :: i,unit
      if(.not.self%status)stop "PDF_WRITE: PDF not allocated"
      open(free_unit(unit),file=trim(pfile))
      write(unit,*)pdf_mean_1d(self),pdf_sdev_1d(self),pdf_skew_1d(self),pdf_curt_1d(self)
      close(unit)
   end subroutine pdf_print_moments_pfile_1d

   !!! 2D PDF
   subroutine pdf_allocate_2d(self,Nvec)
      type(pdf_kernel_2d) :: self
      integer             :: Nvec(2)
      if(self%status)stop "PDF_ALLOCATE: PDF already allocated"
      self%N     = Nvec
      self%Ndata = 0
      allocate(self%x(self%N(1)))
      allocate(self%y(self%N(2)))
      allocate(self%pdf(self%N(1),self%N(2)))
      self%status=.true.
      self%pdf   = 0d0
      self%sigma = 0d0
   end subroutine pdf_allocate_2d
   !
   subroutine pdf_deallocate_2d(self)
      type(pdf_kernel_2d) :: self
      if(.not.self%status)stop "PDF_DEALLOCATE: PDF not allocated"
      self%xmin = 0d0
      self%xmax = 0d0
      self%dx   = 0d0
      self%Sigma= 0d0
      self%N    = 0
      self%Ndata= 0
      deallocate(self%x)
      deallocate(self%y)
      deallocate(self%pdf)
      self%status=.false.
   end subroutine pdf_deallocate_2d
   !
   subroutine pdf_save_2d(self,pfile)
      type(pdf_kernel_2d) :: self
      character(len=*) :: pfile
      integer          :: i,unit
      if(.not.self%status)stop "PDF_SAVE: PDF not allocated"
      open(free_unit(unit),file=trim(pfile))
      write(unit,*)self%N
      write(unit,*)self%xmin
      write(unit,*)self%xmax
      write(unit,*)self%dx
      write(unit,*)self%Ndata
      write(unit,*)self%x
      write(unit,*)self%y
      write(unit,*)self%pdf
      write(unit,*)self%sigma
      write(unit,*)self%status
      write(unit,*)self%variance
      write(unit,*)self%rescale
      close(unit)
   end subroutine pdf_save_2d
   !
   subroutine pdf_read_2d(self,pfile)
      type(pdf_kernel_2d) :: self
      character(len=*) :: pfile
      integer          :: i,N1,N2,unit
      if(.not.self%status)then
         print*,"PDF_READ: PDF not allocated"
      endif
      open(free_unit(unit),file=trim(pfile))
      read(unit,*)N1,N2
      call pdf_allocate(self,[N1,N2])
      read(unit,*)self%xmin
      read(unit,*)self%xmax
      read(unit,*)self%dx
      read(unit,*)self%Ndata
      read(unit,*)self%x
      read(unit,*)self%y
      read(unit,*)self%pdf
      read(unit,*)self%sigma
      read(unit,*)self%status
      read(unit,*)self%variance
      read(unit,*)self%rescale
      close(unit)
   end subroutine pdf_read_2d
   !
   subroutine pdf_set_range_2d(self,a,b)
      type(pdf_kernel_2d) :: self
      real(8)             :: a(2)
      real(8)             :: b(2)
      if(.not.self%status)stop "PDF_SET_RANGE: PDF not allocated"
      self%xmin=a
      self%xmax=b
      self%x = linspace(a(1),b(1),self%N(1),mesh=self%dx(1))
      self%y = linspace(a(2),b(2),self%N(2),mesh=self%dx(2))
   end subroutine pdf_set_range_2d
   !
   subroutine pdf_push_sigma_2d(self,sigma)
      type(pdf_kernel_2d) :: self
      real(8)             :: sigma(2,2)
      if(.not.self%status)stop "PDF_PUSH_SIGMA: PDF not allocated"
      self%sigma = sigma
      self%variance=.true.
   end subroutine pdf_push_sigma_2d
   !
   subroutine pdf_get_sigma_2d(self,sigma)
      type(pdf_kernel_2d) :: self
      real(8)             :: sigma(2,2)
      if(.not.self%status)stop "PDF_GET_SIGMA: PDF not allocated"
      sigma = self%sigma
   end subroutine pdf_get_sigma_2d
   !
   subroutine pdf_sigma_data_2d(self,data,h)
      type(pdf_kernel_2d)    :: self
      real(8),dimension(:,:) :: data
      real(8)                :: h(2,2)
      integer                :: L
      if(.not.self%status)stop "PDF_=SIGMA: PDF not allocated"
      L = size(data,2)
      if(any(shape(data)/=[2,L]))stop "PDF_SET_SIGMA: DATA wrong dimensions"
      h = 0d0
      h(1,1) = get_var(data(1,:))
      h(2,2) = get_var(data(2,:))
      h(1,1) = (1d0/size(data,2))**(2/6d0)*h(1,1) !Silverman's rule of thumb.
      h(2,2) = (1d0/size(data,2))**(2/6d0)*h(2,2) !Silverman's rule of thumb.
   end subroutine pdf_sigma_data_2d
   !
   subroutine pdf_sigma_sdev_2d(self,sdev,Nvec,h)
      type(pdf_kernel_2d) :: self
      real(8)             :: sdev(2)
      integer             :: Nvec(2)
      real(8)             :: h(2,2)
      if(.not.self%status)stop "PDF_SET_SIGMA: PDF not allocated"
      h      = 0d0
      h(1,1) = (1d0/Nvec(1))**(2/6d0)*sdev(1)**2 !Silverman's rule of thumb.
      h(2,2) = (1d0/Nvec(2))**(2/6d0)*sdev(2)**2 !Silverman's rule of thumb.
   end subroutine pdf_sigma_sdev_2d
   !
   subroutine pdf_accumulate_s_2d(self,data,sigma)
      type(pdf_kernel_2d) :: self
      real(8)             :: data(2)
      real(8),optional    :: sigma(2,2)
      real(8)             :: sigma_(2,2)
      integer             :: i,j
      if(.not.self%status)stop "PDF_ACCUMULATE: PDF not allocated"
      if(self%variance)then
         sigma_ = self%sigma
      else
         if(present(sigma))then
            sigma_ = sigma
         else
            stop "PDF_ACCUMULATE: PDF sigma not set or passed"
         endif
      endif
      !
      self%pdf   = self%pdf + gaussian_kernel_2d(self%x,self%y,data,sigma_)
      self%Ndata = self%Ndata + 1
   end subroutine pdf_accumulate_s_2d
   !
   function gaussian_kernel_2d(x,y,mean,sigma) result(gaussian_kernel)
      real(8),intent(in) :: x(:),y(:)
      real(8),intent(in) :: mean(2)
      real(8),intent(in) :: sigma(2,2)
      real(8)            :: gaussian_kernel(size(x),size(y))
      real(8),parameter  :: pi2=2d0*acos(-1d0)
      real(8)            :: detSigma
      real(8)            :: InvSigma(2,2)
      real(8)            :: gaussian_factor
      real(8)            :: arg,xvec(2)
      integer            :: i,j
      !
      detSigma = det(sigma)
      gaussian_factor=1/pi2/sqrt(detSigma)
      !
      InvSigma = Sigma
      call inv(InvSigma)
      do i=1,size(x)
         do j=1,size(y)
            xvec = [x(i),y(j)] - mean
            arg  = dot_product(xvec,matmul(InvSigma,xvec))
            !
            gaussian_kernel(i,j)=gaussian_factor*exp(-0.5d0*arg)
            !
         enddo
      enddo
   end function gaussian_kernel_2d
   !
   subroutine pdf_normalize_2d(self)
      type(pdf_kernel_2d) :: self
      if(.not.self%status)stop "PDF_NORMALIZE: PDF not allocated"
      if(.not.self%rescale)self%pdf = self%pdf/self%Ndata
      self%rescale = .true.
      self%pdf = self%pdf/sum(self%pdf)/product(self%dx)
   end subroutine pdf_normalize_2d
   !
   subroutine pdf_print_pfile_2d(self,pfile,normalize)
      type(pdf_kernel_2d) :: self
      character(len=*)    :: pfile
      logical,optional    :: normalize
      logical             :: normalize_
      integer             :: i,unit
      if(.not.self%status)stop "PDF_WRITE: PDF not allocated"
      normalize_ = .true.; if(present(normalize))normalize_=normalize
      if(normalize_)call pdf_normalize(self)
      call splot3d(pfile,self%x,self%y,self%pdf)
   end subroutine pdf_print_pfile_2d





   !!! HISTOGRAMS
   function histogram_allocate(n) result(h)
      integer,intent(in) :: n
      type(histogram)    :: h
      if(n<=0)then
         print*,"histogram length must be positive integer. n=",n
         stop
      endif
      allocate(h%range(0:n));h%range=0.d0
      allocate(h%bin(0:n))    ;h%bin=0.d0
      h%n=n
   end function histogram_allocate
   !
   subroutine histogram_deallocate(h)
      type(histogram)    :: h
      deallocate(h%range)
      deallocate(h%bin)
      h%n=0
   end subroutine  histogram_deallocate
   !
   subroutine histogram_reset(h)
      type(histogram)    :: h
      h%range=0.d0
      h%bin=0.d0
   end subroutine histogram_reset
   !
   subroutine histogram_set_range_uniform(h,xmin,xmax)
      type(histogram),intent(inout) :: h
      real(8),intent(in)            :: xmin,xmax
      integer                       :: i,n
      real(8)                       :: f1,f2
      if(xmin>=xmax)then
         print*,"histogram%range: xmin must be less than xmax:",xmin,xmax
         stop
      endif
      n=h%n
      do i=0,n
         f1= real(n-i,8)/real(n,8)
         f2= real(i,8)/real(n,8)
         h%range(i) = f1*xmin + f2*xmax
      enddo
      h%bin=0.d0
   end subroutine histogram_set_range_uniform
   !
   subroutine histogram_accumulate(h,x,w)
      type(histogram),intent(inout) :: h
      real(8),intent(in)            :: x,w
      integer                       :: i,index
      index=0
      call find_index(h%n,h%range,x,index)
      if(index>=h%n)then
         print*,"index lies outside valid range of 0 .. n - 1"
         stop
      endif
      h%bin(index)=h%bin(index)+w
   end subroutine histogram_accumulate
   !
   subroutine find_index(n,range,x,index)
      integer,intent(in)                :: n
      real(8),dimension(0:n),intent(in) :: range
      real(8),intent(in)                :: x
      integer,intent(out)               :: index
      integer                           :: i,upper,lower,mid
      if((x<range(0)) .OR. (x>range(n)))then
         print*,"X out of range!"
         return
      endif
      upper=n
      lower=0
      do while((upper-lower>1))
         mid = (upper+lower)/2    !int(dble(upper + lower)/2.d0)
         if( x >= range(mid))then
            lower=mid
         else
            upper=mid
         endif
      enddo
      index=lower
      if(x<range(lower) .OR. x>range(lower+1))then
         print*,"error: x not found within range!"
         stop
      endif
   end subroutine find_index
   !
   subroutine histogram_get_range(h,index,lower,upper)
      type(histogram),intent(in) :: h
      integer,intent(in)         :: index
      real(8),intent(out)        :: lower,upper
      if(index>=h%n)then
         print*,"error: *i lies outside valid range=0...n-1"
         stop
      endif
      lower=h%range(index)
      upper=h%range(index+1)
   end subroutine histogram_get_range
   !
   function histogram_get_value(h,index) result(value)
      type(histogram),intent(in) :: h
      integer,intent(in)         :: index
      real(8)                    :: value
      if(index>=h%n)then
         print*,"error: *index lies outside valid range=0...n-1"
         stop
      endif
      value=h%bin(index)
   end function histogram_get_value
   !
   subroutine histogram_print(h,unit)
      type(histogram),intent(in) :: h
      integer,intent(in)         :: unit
      integer                    :: i,n
      real(8)                    :: lower,upper,bin_value
      n=h%n
      call histogram_get_range(h,0,lower,upper)
      write(unit,"(2F12.7)")lower,0.d0
      do i=0,n-1
         call histogram_get_range(h,i,lower,upper)
         bin_value = histogram_get_value(h,i)
         write(unit,"(2F12.7)")lower,bin_value
         write(unit,"(2F12.7)")upper,bin_value
      enddo
      write(unit,"(2F12.7)")upper,0.d0
   end subroutine histogram_print



END MODULE SF_STAT

