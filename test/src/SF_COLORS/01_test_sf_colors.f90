program test_SF_COLORS
  USE SF_COLORS
  USE ASSERTING
  implicit none
  type(rgb_color) :: color1,color2,color3,my_yellow,mix_color
  integer :: rgb_yellow=16776960

  color1 = red
  color2 = green
  color3 = blue

  my_yellow = red + green
  call assert(rgb(my_yellow),rgb_yellow,"RGB YELLOW")
  write(*,*)my_yellow

  mix_color = 0.6274509803921569d0*red + 0.12549019607843137d0*green + 0.9411764705882353d0*blue
  call assert(rgb(mix_color),10494192,"RGB PURPLE")
  write(*,*)purple
  write(*,*)mix_color
end program test_SF_COLORS
