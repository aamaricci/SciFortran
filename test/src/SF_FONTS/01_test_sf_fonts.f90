program test_SF_FONTS
  USE SF_FONTS
  USE ASSERTING
  implicit none

  character(len=27) :: text="Lorem ipsum dolor sit amet."

  write(*,"(A)")bold(text)
  write(*,"(A)")underline(text)
  write(*,"(A)")highlight(text)
  write(*,"(A)")erased(text)
  write(*,"(A)")font_red(text)
  write(*,"(A)")font_green(text)
  write(*,"(A)")font_yellow(text)
  write(*,"(A)")font_blue(text)
  write(*,"(A)")bold_red(text)
  write(*,"(A)")bold_green(text)
  write(*,"(A)")bold_yellow(text)
  write(*,"(A)")bold_blue(text)
  write(*,"(A)")bg_red(text)
  write(*,"(A)")bg_green(text)
  write(*,"(A)")bg_yellow(text)
  write(*,"(A)")bg_blue(text)

end program test_SF_FONTS
