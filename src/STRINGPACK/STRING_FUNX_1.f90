! -----------------------------------------------
!by David Frank  dave_frank@hotmail.com
! http://home.earthlink.net/~dave_gemini/strings.f90
! -----------------------------------------------
MODULE STRING_FUNX
  implicit none            

  ! Copy (generic) char array to string or string to char array
  ! Clen           returns same as LEN      unless last non-blank char = null
  ! Clen_trim      returns same as LEN_TRIM    "              "
  ! Ctrim          returns same as TRIM        "              "
  ! Count_Items    in string that are blank or comma separated
  ! Reduce_Blanks  in string to 1 blank between items, last char not blank
  ! Replace_Text   in all occurances in string with replacement string
  ! Spack          pack string's chars == extract string's chars
  ! Tally          occurances in string of text arg
  ! Translate      text arg via indexed code table
  ! Upper/Lower    case the text arg

  interface copy
     module procedure copy_a2s, copy_s2a
  end interface copy

contains

  pure function copy_a2s(a)  result (s)    ! copy char array to string
    character,intent(in) :: a(:)
    character(size(a)) :: s
    integer :: i
    do i = 1,size(a)
       s(i:i) = a(i)
    end do
  end function copy_a2s
  ! ------------------------
  pure function copy_s2a(s)  result (a)   ! copy s(1:clen(s)) to char array
    character(*),intent(in) :: s
    character :: a(len(s))
    integer :: i
    do i = 1,len(s)
       a(i) = s(i:i)
    end do
  end function copy_s2a


  !******************************************
  !******************************************
  !******************************************


  pure integer function clen(s)      ! returns same result as len unless:
    character(*),intent(in) :: s       ! last non-blank char is null
    integer :: i
    clen = len(s)
    i = len_trim(s)
    if (s(i:i) == char(0)) clen = i-1  ! len of c string
  end function clen


  !******************************************
  !******************************************
  !******************************************



  pure integer function clen_trim(s) ! returns same result as len_trim unless:
    character(*),intent(in) :: s       ! last char non-blank is null, if true:
    integer :: i                       ! then len of c string is returned, note:
    ! ctrim is only user of this function
    i = len_trim(s) ; clen_trim = i
    if (s(i:i) == char(0)) clen_trim = clen(s)   ! len of c string
  end function clen_trim


  !******************************************
  !******************************************
  !******************************************



  function ctrim(s1)  result(s2)     ! returns same result as trim unless:
    character(*),intent(in)  :: s1     ! last non-blank char is null in which
    character(clen_trim(s1)) :: s2     ! case trailing blanks prior to null
    s2 = s1                            ! are output
  end function ctrim


  !******************************************
  !******************************************
  !******************************************


  integer function count_items(s1)  ! in string or c string that are blank or comma separated
    character(*) :: s1
    character(clen(s1)) :: s
    integer :: i, k
    s = s1                            ! remove possible last char null
    k = 0  ; if (s /= ' ') k = 1      ! string has at least 1 item
    do i = 1,len_trim(s)-1
       if (s(i:i) /= ' '.and.s(i:i) /= ',' &
            .and.s(i+1:i+1) == ' '.or.s(i+1:i+1) == ',') k = k+1
    end do
    count_items = k
  end function count_items


  !******************************************
  !******************************************
  !******************************************


  function reduce_blanks(s)  result (outs)
    character(*)      :: s
    character(len_trim(s)) :: outs
    integer           :: i, k, n
    n = 0  ; k = len_trim(s)          ! k=index last non-blank (may be null)
    do i = 1,k-1                      ! dont process last char yet
       n = n+1 ; outs(n:n) = s(i:i)
       if (s(i:i+1) == '  ') n = n-1  ! backup/discard consecutive output blank
    end do
    n = n+1  ; outs(n:n)  = s(k:k)    ! last non-blank char output (may be null)
    if (n < k) outs(n+1:) = ' '       ! pad trailing blanks
  end function reduce_blanks


  !******************************************
  !******************************************
  !******************************************


  function replace_text (s,text,rep)  result(outs)
    character(*)        :: s,text,rep
    character(len(s)+100) :: outs     ! provide outs with extra 100 char len
    integer             :: i, nt, nr
    outs = s ; nt = len_trim(text) ; nr = len_trim(rep)
    do
       i = index(outs,text(:nt)) ; if (i == 0) exit
       outs = outs(:i-1) // rep(:nr) // outs(i+nt:)
    end do
  end function replace_text


  !******************************************
  !******************************************
  !******************************************


  function spack (s,sex)  result (outs)
    character(*) :: s,sex
    character(len(s)) :: outs
    character :: asex(len(sex))   ! array of sex chars to extract
    integer   :: i, n
    n = 0  ;  asex = copy(sex)
    do i = 1,len(s)
       if (.not.any(s(i:i) == asex)) cycle   ! dont pack char
       n = n+1 ; outs(n:n) = s(i:i)
    end do
    outs(n+1:) = ' '     ! pad with trailing blanks
  end function spack


  integer function tally (s,text)
    character(*) :: s, text
    integer :: i, nt
    tally = 0 ; nt = len_trim(text)
    do i = 1,len(s)-nt+1
       if (s(i:i+nt-1) == text(:nt)) tally = tally+1
    end do
  end function tally


  !******************************************
  !******************************************
  !******************************************


  function translate(s1,codes)  result (s2)
    character(*)       :: s1, codes(2)
    character(len(s1)) :: s2
    character          :: ch
    integer            :: i, j

    do i = 1,len(s1)
       ch = s1(i:i)
       j = index(codes(1),ch) ; if (j > 0) ch = codes(2)(j:j)
       s2(i:i) = ch
    end do
  end function translate


  !******************************************
  !******************************************
  !******************************************


  function upper(s1)  result (s2)
    character(*)       :: s1
    character(len(s1)) :: s2
    character          :: ch
    integer,parameter  :: duc = ichar('a') - ichar('a')
    integer            :: i
    do i = 1,len(s1)
       ch = s1(i:i)
       if (ch >= 'a'.and.ch <= 'z') ch = char(ichar(ch)+duc)
       s2(i:i) = ch
    end do
  end function upper


  !******************************************
  !******************************************
  !******************************************


  function lower(s1)  result (s2)
    character(*)       :: s1
    character(len(s1)) :: s2
    character          :: ch
    integer,parameter  :: duc = ichar('a') - ichar('a')
    integer            :: i
    do i = 1,len(s1)
       ch = s1(i:i)
       if (ch >= 'a'.and.ch <= 'z') ch = char(ichar(ch)-duc)
       s2(i:i) = ch
    end do
  end function lower

END MODULE STRING_FUNX
