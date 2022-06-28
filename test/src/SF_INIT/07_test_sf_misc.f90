program test_SF_MISC
  USE SF_MISC
  USE ASSERTING
  implicit none


  call test_uniq()

  call test_sort()
  
contains


  subroutine test_uniq()
    integer,parameter                :: L=10
    integer,dimension(L)             :: i_array
    integer,dimension(L)             :: i_array_save
    integer,dimension(:),allocatable :: i_array_uniqd1(:)
    integer,dimension(:),allocatable :: i_array_uniqd2(:)
    !
    real(8),dimension(L)             :: d_array
    real(8),dimension(L)             :: d_array_save
    real(8),dimension(:),allocatable :: d_array_uniqd1(:)
    real(8),dimension(:),allocatable :: d_array_uniqd2(:)
    !
    logical,dimension(:),allocatable :: mask
    integer                          :: ndim
    i_array     = [6,3,6,6,6,9,1,6,10,6]
    i_array_save= i_array
    write(*,"(10I4)")i_array
    call uniq_array(i_array,i_array_uniqd1,mask)
    ndim = size(mask)
    allocate(i_array_uniqd2(ndim))
    i_array_uniqd2 = pack(i_array_save,mask)
    write(*,"(5I4)")i_array_uniqd1
    write(*,"(5I4)")i_array_uniqd2
    call assert(i_array_uniqd1,[6,3,9,1,10],"UNIQ INT")
    !
    deallocate(mask)
    !
    d_array     = 1d0*[1,2,3,2,2,9,1,6,7,6]
    d_array_save= d_array
    write(*,"(10F5.1)")d_array
    call uniq_array(d_array,d_array_uniqd1,mask)
    ndim = size(mask)
    allocate(d_array_uniqd2(ndim))
    d_array_uniqd2 = pack(d_array_save,mask)
    write(*,"(6F5.1)")d_array_uniqd1
    write(*,"(6F5.1)")d_array_uniqd2
    call assert(d_array_uniqd2,[1d0,2d0,3d0,9d0,6d0,7d0],"UNIQ DBLE")
  end subroutine test_uniq




  subroutine test_sort()
    USE SF_RANDOM
    USE SF_IOTOOLS, only:str
    integer,parameter :: L=50000
    real(8)           :: t_start,t_stop
    real(8)           :: array1(L),time1
    real(8)           :: array2(L),time2
    real(8)           :: array3(L),time3
    real(8)           :: array_init(L),test_array(10)
    integer           :: indx(L)
    integer           :: i

    test_array = [&
         0.00002894d0,&
         0.00005877d0,&
         0.00006432d0,&
         0.00015180d0,&
         0.00015700d0,&
         0.99981638d0,&
         0.99992071d0,&
         0.99993130d0,&
         0.99995289d0,&
         0.99997527d0]

    call mt_init(123456)
    do i=1,L
       array_init(i) = mersenne()
    enddo
    array1 = array_init
    array2 = array_init
    array3 = array_init

    indx  = 0
    call cpu_time(t_start)
    call sort_insertion(array1,indx)
    call cpu_time(t_stop)
    time1=t_stop-t_start
    call assert([array1(1:5),array1(L-4:L)],test_array,"SORT INSERTION",tol=1d-6)
    write(*,"(A,A)")"Sort Insertion  time: ",str(time1)
    write(*,"(5F12.8,A3,5F12.8)")array1(1:5),"...",array1(L-4:L)


    indx  = 0
    call cpu_time(t_start)
    call sort_qsort(array2,indx)
    call cpu_time(t_stop)
    time2=t_stop-t_start
    call assert([array2(1:5),array2(L-4:L)],test_array,"SORT QUICKSORT 1",tol=1d-6)
    write(*,"(A,A)")"Sort Quicksort1 time: ",str(time2)
    write(*,"(5F12.8,A3,5F12.8)")array2(1:5),"...",array2(L-4:L)




    indx  = 0
    call cpu_time(t_start)
    call sort_quicksort(array3,indx)
    call cpu_time(t_stop)
    time3=t_stop-t_start
    call assert([array2(1:5),array2(L-4:L)],test_array,"SORT QUICKSORT 2",tol=1d-6)
    write(*,"(A,A)")"Sort Quicksort2 time: ",str(time3)
    write(*,"(5F12.8,A3,5F12.8)")array3(1:5),"...",array3(L-4:L)




  end subroutine test_sort



end program
