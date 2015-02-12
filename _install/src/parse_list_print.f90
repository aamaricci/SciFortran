subroutine help_input_node(c)
  type(input_node)   :: c
  character(len=255) :: name
  integer            :: clen
  integer            :: unit,i
  name=c%name
  call s_blank_delete(name)
  p_buffer=trim(name)//" : "//trim(c%comment)
  write(*,"(1x,A)")trim(p_buffer)
end subroutine help_input_node



!---------------------------------------------------------------------
!PURPOSE:
!---------------------------------------------------------------------
subroutine print_input_node(c,file)
  type(input_node) :: c
  character(len=*)    :: file
  select case(c%type)
  case('i')
     if(size(c%i)==1)then
        call i_print_to_file(c%i(1),c%name,file,c%comment)
     else
        call iv_print_to_file(c%i(:),c%name,file,c%comment)
     endif
  case('d')
     if(size(c%d)==1)then
        call d_print_to_file(c%d(1),c%name,file,c%comment)
     else
        call dv_print_to_file(c%d(:),c%name,file,c%comment)
     endif
  case('ch')
     if(size(c%ch)==1)then
        call ch_print_to_file(c%ch(1),c%name,file,c%comment)
     else
        call chv_print_to_file(c%ch(:),c%name,file,c%comment)
     endif
  case('l')
     if(size(c%l)==1)then
        call l_print_to_file(c%l(1),c%name,file,c%comment)
     else
        call lv_print_to_file(c%l(:),c%name,file,c%comment)
     endif
  end select
end subroutine print_input_node

subroutine i_print_to_file(variable,name,file,comment)
  integer          :: variable
  character(len=*) :: name,comment
  character(len=*) :: file
  character(len=255) :: blank=""
  integer          :: clen
  integer          :: unit,i
  call s_blank_delete(name)
  p_buffer=trim(name)//"="//txtfy(variable)
  call s_blank_delete(p_buffer)
  clen=pos_comment-len(trim(p_buffer))
  if(clen<=0)clen=1
  p_buffer=trim(p_buffer)//blank(1:clen)//"!"//trim(comment)
  unit=free_unit()
  open(unit,file="used."//file,position='append',status=trim(file_status));file_status='old'
  write(unit,"(1x,A)")trim(p_buffer)
  write(*,"(1x,A)")trim(p_buffer)
  close(unit)
end subroutine i_print_to_file

subroutine d_print_to_file(variable,name,file,comment)
  real(8)            :: variable
  character(len=*)   :: name,comment
  character(len=*)   :: file
  integer            :: unit
  character(len=255) :: blank=""
  integer            :: clen
  call s_blank_delete(name)
  p_buffer=trim(name)//"="//txtfy(variable)
  call s_blank_delete(p_buffer)
  clen=pos_comment-len(trim(p_buffer))
  if(clen<=0)clen=1
  p_buffer=trim(p_buffer)//blank(1:clen)//"!"//trim(comment)
  unit=free_unit()
  open(unit,file="used."//file,position='append',status=trim(file_status));file_status='old'
  write(unit,"(1x,A)")trim(p_buffer)
  write(*,"(1x,A)")trim(p_buffer)
  close(unit)
end subroutine d_print_to_file

subroutine ch_print_to_file(variable,name,file,comment)
  character(len=*) :: variable
  character(len=*) :: name,comment
  character(len=*) :: file
  integer          :: unit
  character(len=255) :: blank=""
  integer          :: clen
  call s_blank_delete(name)
  p_buffer=trim(name)//"="//trim(variable)
  call s_blank_delete(p_buffer)
  clen=pos_comment-len(trim(p_buffer))
  if(clen<=0)clen=1
  p_buffer=trim(p_buffer)//blank(1:clen)//"!"//trim(comment)
  unit=free_unit()
  open(unit,file="used."//file,position='append',status=trim(file_status));file_status='old'
  write(unit,"(1x,A)")trim(p_buffer)
  write(*,"(1x,A)")trim(p_buffer)
  close(unit)
end subroutine ch_print_to_file

subroutine l_print_to_file(variable,name,file,comment)
  logical          :: variable
  character(len=*) :: name,comment
  character(len=*) :: file
  integer          :: unit
  character(len=255) :: blank=""
  integer          :: clen
  call s_blank_delete(name)
  p_buffer=trim(name)//"="//txtfy(variable)
  call s_blank_delete(p_buffer)
  clen=pos_comment-len(trim(p_buffer))
  if(clen<=0)clen=1
  p_buffer=trim(p_buffer)//blank(1:clen)//"!"//trim(comment)
  unit=free_unit()
  open(unit,file="used."//file,position='append',status=trim(file_status));file_status='old'
  write(unit,"(1x,A)")trim(p_buffer)
  write(*,"(1x,A)")trim(p_buffer)
  close(unit)
end subroutine l_print_to_file

!==========================1-dimension==================================

subroutine iv_print_to_file(variable,name,file,comment)
  integer,dimension(:) :: variable
  character(len=*)     :: name,comment
  character(len=*)     :: file
  integer              :: unit,i
  character(len=255) :: blank=""
  integer          :: clen
  call s_blank_delete(name)
  p_buffer=trim(name)//"="
  do i=1,size(variable)-1
     p_buffer=trim(p_buffer)//trim(txtfy(variable(i)))//","
  end do
  p_buffer=trim(p_buffer)//trim(txtfy(variable(size(variable))))
  call s_blank_delete(p_buffer)
  clen=pos_comment-len(trim(p_buffer))
  if(clen<=0)clen=1
  p_buffer=trim(p_buffer)//blank(1:clen)//"!"//trim(comment)
  unit=free_unit()
  open(unit,file="used."//file,position='append',status=trim(file_status));file_status='old'
  write(unit,"(1x,A)")trim(p_buffer)
  write(*,"(1x,A)")trim(p_buffer)
  close(unit)
end subroutine iv_print_to_file

subroutine dv_print_to_file(variable,name,file,comment)
  real(8),dimension(:) :: variable
  character(len=*)     :: name,comment
  character(len=*)     :: file
  integer              :: unit,i
  character(len=255) :: blank=""
  integer          :: clen
  call s_blank_delete(name)
  p_buffer=trim(name)//"="
  do i=1,size(variable)-1
     p_buffer=trim(p_buffer)//trim(txtfy(variable(i)))//","
  end do
  p_buffer=trim(p_buffer)//trim(txtfy(variable(size(variable))))
  call s_blank_delete(p_buffer)
  clen=pos_comment-len(trim(p_buffer))
  if(clen<=0)clen=1
  p_buffer=trim(p_buffer)//blank(1:clen)//"!"//trim(comment)
  unit=free_unit()
  open(unit,file="used."//file,position='append',status=trim(file_status));file_status='old'
  write(unit,"(1x,A)")trim(p_buffer)
  write(*,"(1x,A)")trim(p_buffer)
  close(unit)
end subroutine dv_print_to_file


subroutine chv_print_to_file(variable,name,file,comment)
  character(len=*),dimension(:) :: variable
  character(len=*)     :: name,comment
  character(len=*)     :: file
  integer              :: unit,i
  character(len=255) :: blank=""
  integer          :: clen
  call s_blank_delete(name)
  p_buffer=trim(name)//"="
  do i=1,size(variable)-1
     p_buffer=trim(p_buffer)//trim(variable(i))//","
  end do
  p_buffer=trim(p_buffer)//trim(variable(size(variable)))
  call s_blank_delete(p_buffer)
  clen=pos_comment-len(trim(p_buffer))
  if(clen<=0)clen=1
  p_buffer=trim(p_buffer)//blank(1:clen)//"!"//trim(comment)
  unit=free_unit()
  open(unit,file="used."//file,position='append',status=trim(file_status));file_status='old'
  write(unit,"(1x,A)")trim(p_buffer)
  write(*,"(1x,A)")trim(p_buffer)
  close(unit)
end subroutine chv_print_to_file

subroutine lv_print_to_file(variable,name,file,comment)
  logical,dimension(:) :: variable
  character(len=*)     :: name,comment
  character(len=*)     :: file
  integer              :: unit,i
  character(len=255) :: blank=""
  integer          :: clen
  call s_blank_delete(name)
  p_buffer=trim(name)//"="
  do i=1,size(variable)-1
     p_buffer=trim(p_buffer)//trim(txtfy(variable(i)))//","
  end do
  p_buffer=trim(p_buffer)//trim(txtfy(variable(size(variable))))
  call s_blank_delete(p_buffer)
  clen=pos_comment-len(trim(p_buffer))
  if(clen<=0)clen=1
  p_buffer=trim(p_buffer)//blank(1:clen)//"!"//trim(comment)
  unit=free_unit()
  open(unit,file="used."//file,position='append',status=trim(file_status));file_status='old'
  write(unit,"(1x,A)")trim(p_buffer)
  write(*,"(1x,A)")trim(p_buffer)
  close(unit)
end subroutine lv_print_to_file


! subroutine im_print_to_file(variable,name,file,comment)
!   integer,dimension(:,:) :: variable
!   character(len=*)     :: name,comment
!   character(len=*)     :: file
!   integer              :: unit,i,j
! character(len=255) :: blank=""
! integer          :: clen
!   call s_blank_delete(name)
!   p_buffer=trim(name)//"="
!   do i=1,size(variable,1)-1
!      do j=1,size(variable,2)
!         p_buffer=trim(p_buffer)//txtfy(variable(i,j))//","
!      end do
!   end do
!   i=size(variable,1)
!   do j=1,size(variable,2)-1
!      p_buffer=trim(p_buffer)//txtfy(variable(i,j))//","
!   end do
!   j=size(variable,2)
!   p_buffer=trim(p_buffer)//txtfy(variable(i,j))
!   call s_blank_delete(p_buffer)
! clen=pos_comment-len(trim(p_buffer))
! if(clen<=0)clen=1
! p_buffer=trim(p_buffer)//blank(1:clen)//"!"//trim(comment)
!   unit=free_unit()
!   open(unit,file="used."//file,position='append',status=trim(file_status));file_status='old'
!   write(unit,"(1x,A)")trim(p_buffer)
!   write(*,"(1x,A)")trim(p_buffer)
!   close(unit)
! end subroutine im_print_to_file

! subroutine dm_print_to_file(variable,name,file,comment)
!   real(8),dimension(:,:) :: variable
!   character(len=*)     :: name,comment
!   character(len=*)     :: file
!   integer              :: unit,i,j
! character(len=255) :: blank=""
! integer          :: clen
!   call s_blank_delete(name)
!   p_buffer=trim(name)//"="
!   do i=1,size(variable,1)-1
!      do j=1,size(variable,2)
!         p_buffer=trim(p_buffer)//txtfy(variable(i,j))//","
!      end do
!   end do
!   i=size(variable,1)
!   do j=1,size(variable,2)-1
!      p_buffer=trim(p_buffer)//txtfy(variable(i,j))//","
!   end do
!   j=size(variable,2)
!   p_buffer=trim(p_buffer)//txtfy(variable(i,j))
!   call s_blank_delete(p_buffer)
! clen=pos_comment-len(trim(p_buffer))
! if(clen<=0)clen=1
! p_buffer=trim(p_buffer)//blank(1:clen)//"!"//trim(comment)
!   unit=free_unit()
!   open(unit,file="used."//file,position='append',status=trim(file_status));file_status='old'
!   write(unit,"(1x,A)")trim(p_buffer)
!   write(*,"(1x,A)")trim(p_buffer)
!   close(unit)
! end subroutine dm_print_to_file
