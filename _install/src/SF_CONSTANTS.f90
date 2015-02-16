MODULE SF_MPI_VARS
  implicit none
  integer          :: MPIID=0
  integer          :: MPISIZE=1
  integer          :: MPIERR
  character(len=3) :: MPICHAR
  ! public :: init_mpi
  ! public :: finalize_mpi
  ! contains
  !   subroutine init_mpi
  !     call MPI_INIT(mpiERR)
  !     call MPI_COMM_RANK(MPI_COMM_WORLD,mpiID,mpiERR)
  !     call MPI_COMM_SIZE(MPI_COMM_WORLD,mpiSIZE,mpiERR)
  !     write(*,"(A,I4,A,I4,A)")'Processor ',mpiID,' of ',mpiSIZE,' is alive'
  !     call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  !   end subroutine init_mpi
  !   subroutine finalize_mpi
  !     call MPI_FINALIZE(mpiERR)
  !   end subroutine finalize_mpi
END MODULE SF_MPI_VARS

MODULE SF_OMP_VARS
  implicit none
  integer :: omp_num_threads
  integer :: omp_id
  integer :: omp_size
END MODULE SF_OMP_VARS



MODULE SF_COLORS
  implicit none

  type rgb_color
     integer :: r,g,b
  end type rgb_color
  type(rgb_color),parameter ::  black               =rgb_color(0,0,0)
  type(rgb_color),parameter ::  red                 =rgb_color(255,0,0)
  type(rgb_color),parameter ::  green               =rgb_color(0, 255, 0)
  type(rgb_color),parameter ::  orange              =rgb_color(255,193,0)
  type(rgb_color),parameter ::  blue                =rgb_color(0, 0, 255)
  type(rgb_color),parameter ::  yellow              =rgb_color(255,255,0)
  type(rgb_color),parameter ::  cyan                =rgb_color(0,255,255)
  type(rgb_color),parameter ::  magenta             =rgb_color(159, 0, 159)
  !X11/rgb.txt color codes:
  type(rgb_color),parameter :: snow                 =rgb_color(255,250,250) 
  type(rgb_color),parameter :: ghostwhite           =rgb_color(248,248,255) 
  type(rgb_color),parameter :: whitesmoke           =rgb_color(245,245,245) 
  type(rgb_color),parameter :: gainsboro            =rgb_color(220,220,220) 
  type(rgb_color),parameter :: floralwhite          =rgb_color(255,250,240) 
  type(rgb_color),parameter :: oldlace              =rgb_color(253,245,230) 
  type(rgb_color),parameter :: linen                =rgb_color(250,240,230) 
  type(rgb_color),parameter :: antiquewhite         =rgb_color(250,235,215) 
  type(rgb_color),parameter :: papayawhip           =rgb_color(255,239,213) 
  type(rgb_color),parameter :: blanchedalmond       =rgb_color(255,235,205) 
  type(rgb_color),parameter :: bisque               =rgb_color(255,228,196) 
  type(rgb_color),parameter :: peachpuff            =rgb_color(255,218,185) 
  type(rgb_color),parameter :: navajowhite          =rgb_color(255,222,173) 
  type(rgb_color),parameter :: moccasin             =rgb_color(255,228,181) 
  type(rgb_color),parameter :: cornsilk             =rgb_color(255,248,220) 
  type(rgb_color),parameter :: ivory                =rgb_color(255,255,240) 
  type(rgb_color),parameter :: lemonchiffon         =rgb_color(255,250,205) 
  type(rgb_color),parameter :: seashell             =rgb_color(255,245,238) 
  type(rgb_color),parameter :: honeydew             =rgb_color(240,255,240) 
  type(rgb_color),parameter :: mintcream            =rgb_color(245,255,250) 
  type(rgb_color),parameter :: azure                =rgb_color(240,255,255) 
  type(rgb_color),parameter :: aliceblue            =rgb_color(240,248,255) 
  type(rgb_color),parameter :: lavender             =rgb_color(230,230,250) 
  type(rgb_color),parameter :: lavenderblush        =rgb_color(255,240,245) 
  type(rgb_color),parameter :: mistyrose            =rgb_color(255,228,225) 
  type(rgb_color),parameter :: white                =rgb_color(255,255,255) 
  type(rgb_color),parameter :: darkslategray        =rgb_color(47,79,79) 
  type(rgb_color),parameter :: darkslategrey        =rgb_color(47,79,79) 
  type(rgb_color),parameter :: dimgray              =rgb_color(105,105,105) 
  type(rgb_color),parameter :: dimgrey              =rgb_color(105,105,105) 
  type(rgb_color),parameter :: slategray            =rgb_color(112,128,144) 
  type(rgb_color),parameter :: slategrey            =rgb_color(112,128,144) 
  type(rgb_color),parameter :: lightslategray       =rgb_color(119,136,153) 
  type(rgb_color),parameter :: lightslategrey       =rgb_color(119,136,153) 
  type(rgb_color),parameter :: gray                 =rgb_color(190,190,190) 
  type(rgb_color),parameter :: grey                 =rgb_color(190,190,190) 
  type(rgb_color),parameter :: lightgrey            =rgb_color(211,211,211) 
  type(rgb_color),parameter :: lightgray            =rgb_color(211,211,211) 
  type(rgb_color),parameter :: midnightblue         =rgb_color(25,25,112) 
  type(rgb_color),parameter :: navy                 =rgb_color(0,0,128) 
  type(rgb_color),parameter :: navyblue             =rgb_color(0,0,128) 
  type(rgb_color),parameter :: cornflowerblue       =rgb_color(100,149,237) 
  type(rgb_color),parameter :: darkslateblue        =rgb_color(72,61,139) 
  type(rgb_color),parameter :: slateblue            =rgb_color(106,90,205) 
  type(rgb_color),parameter :: mediumslateblue      =rgb_color(123,104,238) 
  type(rgb_color),parameter :: lightslateblue       =rgb_color(132,112,255) 
  type(rgb_color),parameter :: mediumblue           =rgb_color(0,0,205) 
  type(rgb_color),parameter :: royalblue            =rgb_color(65,105,225) 
  type(rgb_color),parameter :: dodgerblue           =rgb_color(30,144,255) 
  type(rgb_color),parameter :: deepskyblue          =rgb_color(0,191,255) 
  type(rgb_color),parameter :: skyblue              =rgb_color(135,206,235) 
  type(rgb_color),parameter :: lightskyblue         =rgb_color(135,206,250) 
  type(rgb_color),parameter :: steelblue            =rgb_color(70,130,180) 
  type(rgb_color),parameter :: lightsteelblue       =rgb_color(176,196,222) 
  type(rgb_color),parameter :: lightblue            =rgb_color(173,216,230) 
  type(rgb_color),parameter :: powderblue           =rgb_color(176,224,230) 
  type(rgb_color),parameter :: paleturquoise        =rgb_color(175,238,238) 
  type(rgb_color),parameter :: darkturquoise        =rgb_color(0,206,209) 
  type(rgb_color),parameter :: mediumturquoise      =rgb_color(72,209,204) 
  type(rgb_color),parameter :: turquoise            =rgb_color(64,224,208) 
  type(rgb_color),parameter :: lightcyan            =rgb_color(224,255,255) 
  type(rgb_color),parameter :: cadetblue            =rgb_color(95,158,160) 
  type(rgb_color),parameter :: mediumaquamarine     =rgb_color(102,205,170) 
  type(rgb_color),parameter :: aquamarine           =rgb_color(127,255,212) 
  type(rgb_color),parameter :: darkgreen            =rgb_color(0,100,0) 
  type(rgb_color),parameter :: darkolivegreen       =rgb_color(85,107,47) 
  type(rgb_color),parameter :: darkseagreen         =rgb_color(143,188,143) 
  type(rgb_color),parameter :: seagreen             =rgb_color(46,139,87) 
  type(rgb_color),parameter :: mediumseagreen       =rgb_color(60,179,113) 
  type(rgb_color),parameter :: lightseagreen        =rgb_color(32,178,170) 
  type(rgb_color),parameter :: palegreen            =rgb_color(152,251,152) 
  type(rgb_color),parameter :: springgreen          =rgb_color(0,255,127) 
  type(rgb_color),parameter :: lawngreen            =rgb_color(124,252,0) 
  type(rgb_color),parameter :: chartreuse           =rgb_color(127,255,0) 
  type(rgb_color),parameter :: mediumspringgreen    =rgb_color(0,250,154) 
  type(rgb_color),parameter :: greenyellow          =rgb_color(173,255,47) 
  type(rgb_color),parameter :: limegreen            =rgb_color(50,205,50) 
  type(rgb_color),parameter :: yellowgreen          =rgb_color(154,205,50) 
  type(rgb_color),parameter :: forestgreen          =rgb_color(34,139,34) 
  type(rgb_color),parameter :: olivedrab            =rgb_color(107,142,35) 
  type(rgb_color),parameter :: darkkhaki            =rgb_color(189,183,107) 
  type(rgb_color),parameter :: khaki                =rgb_color(240,230,140) 
  type(rgb_color),parameter :: palegoldenrod        =rgb_color(238,232,170) 
  type(rgb_color),parameter :: lightgoldenrodyellow =rgb_color(250,250,210) 
  type(rgb_color),parameter :: lightyellow          =rgb_color(255,255,224) 
  type(rgb_color),parameter :: gold                 =rgb_color(255,215,0) 
  type(rgb_color),parameter :: lightgoldenrod       =rgb_color(238,221,130) 
  type(rgb_color),parameter :: goldenrod            =rgb_color(218,165,32) 
  type(rgb_color),parameter :: darkgoldenrod        =rgb_color(184,134,11) 
  type(rgb_color),parameter :: rosybrown            =rgb_color(188,143,143) 
  type(rgb_color),parameter :: indianred            =rgb_color(205,92,92) 
  type(rgb_color),parameter :: saddlebrown          =rgb_color(139,69,19) 
  type(rgb_color),parameter :: sienna               =rgb_color(160,82,45) 
  type(rgb_color),parameter :: peru                 =rgb_color(205,133,63) 
  type(rgb_color),parameter :: burlywood            =rgb_color(222,184,135) 
  type(rgb_color),parameter :: beige                =rgb_color(245,245,220) 
  type(rgb_color),parameter :: wheat                =rgb_color(245,222,179) 
  type(rgb_color),parameter :: sandybrown           =rgb_color(244,164,96) 
  type(rgb_color),parameter :: chocolate            =rgb_color(210,105,30) 
  type(rgb_color),parameter :: firebrick            =rgb_color(178,34,34) 
  type(rgb_color),parameter :: brown                =rgb_color(165,42,42) 
  type(rgb_color),parameter :: darksalmon           =rgb_color(233,150,122) 
  type(rgb_color),parameter :: salmon               =rgb_color(250,128,114) 
  type(rgb_color),parameter :: lightsalmon          =rgb_color(255,160,122) 
  type(rgb_color),parameter :: darkorange           =rgb_color(255,140,0) 
  type(rgb_color),parameter :: coral                =rgb_color(255,127,80) 
  type(rgb_color),parameter :: lightcoral           =rgb_color(240,128,128) 
  type(rgb_color),parameter :: tomato               =rgb_color(255,99,71) 
  type(rgb_color),parameter :: orangered            =rgb_color(255,69,0) 
  type(rgb_color),parameter :: hotpink              =rgb_color(255,105,180) 
  type(rgb_color),parameter :: deeppink             =rgb_color(255,20,147) 
  type(rgb_color),parameter :: pink                 =rgb_color(255,192,203) 
  type(rgb_color),parameter :: lightpink            =rgb_color(255,182,193) 
  type(rgb_color),parameter :: palevioletred        =rgb_color(219,112,147) 
  type(rgb_color),parameter :: maroon               =rgb_color(176,48,96) 
  type(rgb_color),parameter :: mediumvioletred      =rgb_color(199,21,133) 
  type(rgb_color),parameter :: violetred            =rgb_color(208,32,144) 
  type(rgb_color),parameter :: violet               =rgb_color(238,130,238) 
  type(rgb_color),parameter :: plum                 =rgb_color(221,160,221) 
  type(rgb_color),parameter :: orchid               =rgb_color(218,112,214) 
  type(rgb_color),parameter :: mediumorchid         =rgb_color(186,85,211) 
  type(rgb_color),parameter :: darkorchid           =rgb_color(153,50,204) 
  type(rgb_color),parameter :: darkviolet           =rgb_color(148,0,211) 
  type(rgb_color),parameter :: blueviolet           =rgb_color(138,43,226) 
  type(rgb_color),parameter :: purple               =rgb_color(160,32,240) 
  type(rgb_color),parameter :: mediumpurple         =rgb_color(147,112,219) 
  type(rgb_color),parameter :: thistle              =rgb_color(216,191,216) 
  type(rgb_color),parameter :: snow1                =rgb_color(255,250,250) 
  type(rgb_color),parameter :: snow2                =rgb_color(238,233,233) 
  type(rgb_color),parameter :: snow3                =rgb_color(205,201,201) 
  type(rgb_color),parameter :: snow4                =rgb_color(139,137,137) 
  type(rgb_color),parameter :: seashell1            =rgb_color(255,245,238) 
  type(rgb_color),parameter :: seashell2            =rgb_color(238,229,222) 
  type(rgb_color),parameter :: seashell3            =rgb_color(205,197,191) 
  type(rgb_color),parameter :: seashell4            =rgb_color(139,134,130) 
  type(rgb_color),parameter :: antiquewhite1        =rgb_color(255,239,219) 
  type(rgb_color),parameter :: antiquewhite2        =rgb_color(238,223,204) 
  type(rgb_color),parameter :: antiquewhite3        =rgb_color(205,192,176) 
  type(rgb_color),parameter :: antiquewhite4        =rgb_color(139,131,120) 
  type(rgb_color),parameter :: bisque1              =rgb_color(255,228,196) 
  type(rgb_color),parameter :: bisque2              =rgb_color(238,213,183) 
  type(rgb_color),parameter :: bisque3              =rgb_color(205,183,158) 
  type(rgb_color),parameter :: bisque4              =rgb_color(139,125,107) 
  type(rgb_color),parameter :: peachpuff1           =rgb_color(255,218,185) 
  type(rgb_color),parameter :: peachpuff2           =rgb_color(238,203,173) 
  type(rgb_color),parameter :: peachpuff3           =rgb_color(205,175,149) 
  type(rgb_color),parameter :: peachpuff4           =rgb_color(139,119,101) 
  type(rgb_color),parameter :: navajowhite1         =rgb_color(255,222,173) 
  type(rgb_color),parameter :: navajowhite2         =rgb_color(238,207,161) 
  type(rgb_color),parameter :: navajowhite3         =rgb_color(205,179,139) 
  type(rgb_color),parameter :: navajowhite4         =rgb_color(139,121,94) 
  type(rgb_color),parameter :: lemonchiffon1        =rgb_color(255,250,205) 
  type(rgb_color),parameter :: lemonchiffon2        =rgb_color(238,233,191) 
  type(rgb_color),parameter :: lemonchiffon3        =rgb_color(205,201,165) 
  type(rgb_color),parameter :: lemonchiffon4        =rgb_color(139,137,112) 
  type(rgb_color),parameter :: cornsilk1            =rgb_color(255,248,220) 
  type(rgb_color),parameter :: cornsilk2            =rgb_color(238,232,205) 
  type(rgb_color),parameter :: cornsilk3            =rgb_color(205,200,177) 
  type(rgb_color),parameter :: cornsilk4            =rgb_color(139,136,120) 
  type(rgb_color),parameter :: ivory1               =rgb_color(255,255,240) 
  type(rgb_color),parameter :: ivory2               =rgb_color(238,238,224) 
  type(rgb_color),parameter :: ivory3               =rgb_color(205,205,193) 
  type(rgb_color),parameter :: ivory4               =rgb_color(139,139,131) 
  type(rgb_color),parameter :: honeydew1            =rgb_color(240,255,240) 
  type(rgb_color),parameter :: honeydew2            =rgb_color(224,238,224) 
  type(rgb_color),parameter :: honeydew3            =rgb_color(193,205,193) 
  type(rgb_color),parameter :: honeydew4            =rgb_color(131,139,131) 
  type(rgb_color),parameter :: lavenderblush1       =rgb_color(255,240,245) 
  type(rgb_color),parameter :: lavenderblush2       =rgb_color(238,224,229) 
  type(rgb_color),parameter :: lavenderblush3       =rgb_color(205,193,197) 
  type(rgb_color),parameter :: lavenderblush4       =rgb_color(139,131,134) 
  type(rgb_color),parameter :: mistyrose1           =rgb_color(255,228,225) 
  type(rgb_color),parameter :: mistyrose2           =rgb_color(238,213,210) 
  type(rgb_color),parameter :: mistyrose3           =rgb_color(205,183,181) 
  type(rgb_color),parameter :: mistyrose4           =rgb_color(139,125,123) 
  type(rgb_color),parameter :: azure1               =rgb_color(240,255,255) 
  type(rgb_color),parameter :: azure2               =rgb_color(224,238,238) 
  type(rgb_color),parameter :: azure3               =rgb_color(193,205,205) 
  type(rgb_color),parameter :: azure4               =rgb_color(131,139,139) 
  type(rgb_color),parameter :: slateblue1           =rgb_color(131,111,255) 
  type(rgb_color),parameter :: slateblue2           =rgb_color(122,103,238) 
  type(rgb_color),parameter :: slateblue3           =rgb_color(105,89,205) 
  type(rgb_color),parameter :: slateblue4           =rgb_color(71,60,139) 
  type(rgb_color),parameter :: royalblue1           =rgb_color(72,118,255) 
  type(rgb_color),parameter :: royalblue2           =rgb_color(67,110,238) 
  type(rgb_color),parameter :: royalblue3           =rgb_color(58,95,205) 
  type(rgb_color),parameter :: royalblue4           =rgb_color(39,64,139) 
  type(rgb_color),parameter :: blue1                =rgb_color(0,0,255) 
  type(rgb_color),parameter :: blue2                =rgb_color(0,0,238) 
  type(rgb_color),parameter :: blue3                =rgb_color(0,0,205) 
  type(rgb_color),parameter :: blue4                =rgb_color(0,0,139) 
  type(rgb_color),parameter :: dodgerblue1          =rgb_color(30,144,255) 
  type(rgb_color),parameter :: dodgerblue2          =rgb_color(28,134,238) 
  type(rgb_color),parameter :: dodgerblue3          =rgb_color(24,116,205) 
  type(rgb_color),parameter :: dodgerblue4          =rgb_color(16,78,139) 
  type(rgb_color),parameter :: steelblue1           =rgb_color(99,184,255) 
  type(rgb_color),parameter :: steelblue2           =rgb_color(92,172,238) 
  type(rgb_color),parameter :: steelblue3           =rgb_color(79,148,205) 
  type(rgb_color),parameter :: steelblue4           =rgb_color(54,100,139) 
  type(rgb_color),parameter :: deepskyblue1         =rgb_color(0,191,255) 
  type(rgb_color),parameter :: deepskyblue2         =rgb_color(0,178,238) 
  type(rgb_color),parameter :: deepskyblue3         =rgb_color(0,154,205) 
  type(rgb_color),parameter :: deepskyblue4         =rgb_color(0,104,139) 
  type(rgb_color),parameter :: skyblue1             =rgb_color(135,206,255) 
  type(rgb_color),parameter :: skyblue2             =rgb_color(126,192,238) 
  type(rgb_color),parameter :: skyblue3             =rgb_color(108,166,205) 
  type(rgb_color),parameter :: skyblue4             =rgb_color(74,112,139) 
  type(rgb_color),parameter :: lightskyblue1        =rgb_color(176,226,255) 
  type(rgb_color),parameter :: lightskyblue2        =rgb_color(164,211,238) 
  type(rgb_color),parameter :: lightskyblue3        =rgb_color(141,182,205) 
  type(rgb_color),parameter :: lightskyblue4        =rgb_color(96,123,139) 
  type(rgb_color),parameter :: slategray1           =rgb_color(198,226,255) 
  type(rgb_color),parameter :: slategray2           =rgb_color(185,211,238) 
  type(rgb_color),parameter :: slategray3           =rgb_color(159,182,205) 
  type(rgb_color),parameter :: slategray4           =rgb_color(108,123,139) 
  type(rgb_color),parameter :: lightsteelblue1      =rgb_color(202,225,255) 
  type(rgb_color),parameter :: lightsteelblue2      =rgb_color(188,210,238) 
  type(rgb_color),parameter :: lightsteelblue3      =rgb_color(162,181,205) 
  type(rgb_color),parameter :: lightsteelblue4      =rgb_color(110,123,139) 
  type(rgb_color),parameter :: lightblue1           =rgb_color(191,239,255) 
  type(rgb_color),parameter :: lightblue2           =rgb_color(178,223,238) 
  type(rgb_color),parameter :: lightblue3           =rgb_color(154,192,205) 
  type(rgb_color),parameter :: lightblue4           =rgb_color(104,131,139) 
  type(rgb_color),parameter :: lightcyan1           =rgb_color(224,255,255) 
  type(rgb_color),parameter :: lightcyan2           =rgb_color(209,238,238) 
  type(rgb_color),parameter :: lightcyan3           =rgb_color(180,205,205) 
  type(rgb_color),parameter :: lightcyan4           =rgb_color(122,139,139) 
  type(rgb_color),parameter :: paleturquoise1       =rgb_color(187,255,255) 
  type(rgb_color),parameter :: paleturquoise2       =rgb_color(174,238,238) 
  type(rgb_color),parameter :: paleturquoise3       =rgb_color(150,205,205) 
  type(rgb_color),parameter :: paleturquoise4       =rgb_color(102,139,139) 
  type(rgb_color),parameter :: cadetblue1           =rgb_color(152,245,255) 
  type(rgb_color),parameter :: cadetblue2           =rgb_color(142,229,238) 
  type(rgb_color),parameter :: cadetblue3           =rgb_color(122,197,205) 
  type(rgb_color),parameter :: cadetblue4           =rgb_color(83,134,139) 
  type(rgb_color),parameter :: turquoise1           =rgb_color(0,245,255) 
  type(rgb_color),parameter :: turquoise2           =rgb_color(0,229,238) 
  type(rgb_color),parameter :: turquoise3           =rgb_color(0,197,205) 
  type(rgb_color),parameter :: turquoise4           =rgb_color(0,134,139) 
  type(rgb_color),parameter :: cyan1                =rgb_color(0,255,255) 
  type(rgb_color),parameter :: cyan2                =rgb_color(0,238,238) 
  type(rgb_color),parameter :: cyan3                =rgb_color(0,205,205) 
  type(rgb_color),parameter :: cyan4                =rgb_color(0,139,139) 
  type(rgb_color),parameter :: darkslategray1       =rgb_color(151,255,255) 
  type(rgb_color),parameter :: darkslategray2       =rgb_color(141,238,238) 
  type(rgb_color),parameter :: darkslategray3       =rgb_color(121,205,205) 
  type(rgb_color),parameter :: darkslategray4       =rgb_color(82,139,139) 
  type(rgb_color),parameter :: aquamarine1          =rgb_color(127,255,212) 
  type(rgb_color),parameter :: aquamarine2          =rgb_color(118,238,198) 
  type(rgb_color),parameter :: aquamarine3          =rgb_color(102,205,170) 
  type(rgb_color),parameter :: aquamarine4          =rgb_color(69,139,116) 
  type(rgb_color),parameter :: darkseagreen1        =rgb_color(193,255,193) 
  type(rgb_color),parameter :: darkseagreen2        =rgb_color(180,238,180) 
  type(rgb_color),parameter :: darkseagreen3        =rgb_color(155,205,155) 
  type(rgb_color),parameter :: darkseagreen4        =rgb_color(105,139,105) 
  type(rgb_color),parameter :: seagreen1            =rgb_color(84,255,159) 
  type(rgb_color),parameter :: seagreen2            =rgb_color(78,238,148) 
  type(rgb_color),parameter :: seagreen3            =rgb_color(67,205,128) 
  type(rgb_color),parameter :: seagreen4            =rgb_color(46,139,87) 
  type(rgb_color),parameter :: palegreen1           =rgb_color(154,255,154) 
  type(rgb_color),parameter :: palegreen2           =rgb_color(144,238,144) 
  type(rgb_color),parameter :: palegreen3           =rgb_color(124,205,124) 
  type(rgb_color),parameter :: palegreen4           =rgb_color(84,139,84) 
  type(rgb_color),parameter :: springgreen1         =rgb_color(0,255,127) 
  type(rgb_color),parameter :: springgreen2         =rgb_color(0,238,118) 
  type(rgb_color),parameter :: springgreen3         =rgb_color(0,205,102) 
  type(rgb_color),parameter :: springgreen4         =rgb_color(0,139,69) 
  type(rgb_color),parameter :: green1               =rgb_color(0,255,0) 
  type(rgb_color),parameter :: green2               =rgb_color(0,238,0) 
  type(rgb_color),parameter :: green3               =rgb_color(0,205,0) 
  type(rgb_color),parameter :: green4               =rgb_color(0,139,0) 
  type(rgb_color),parameter :: chartreuse1          =rgb_color(127,255,0) 
  type(rgb_color),parameter :: chartreuse2          =rgb_color(118,238,0) 
  type(rgb_color),parameter :: chartreuse3          =rgb_color(102,205,0) 
  type(rgb_color),parameter :: chartreuse4          =rgb_color(69,139,0) 
  type(rgb_color),parameter :: olivedrab1           =rgb_color(192,255,62) 
  type(rgb_color),parameter :: olivedrab2           =rgb_color(179,238,58) 
  type(rgb_color),parameter :: olivedrab3           =rgb_color(154,205,50) 
  type(rgb_color),parameter :: olivedrab4           =rgb_color(105,139,34) 
  type(rgb_color),parameter :: darkolivegreen1      =rgb_color(202,255,112) 
  type(rgb_color),parameter :: darkolivegreen2      =rgb_color(188,238,104) 
  type(rgb_color),parameter :: darkolivegreen3      =rgb_color(162,205,90) 
  type(rgb_color),parameter :: darkolivegreen4      =rgb_color(110,139,61) 
  type(rgb_color),parameter :: khaki1               =rgb_color(255,246,143) 
  type(rgb_color),parameter :: khaki2               =rgb_color(238,230,133) 
  type(rgb_color),parameter :: khaki3               =rgb_color(205,198,115) 
  type(rgb_color),parameter :: khaki4               =rgb_color(139,134,78) 
  type(rgb_color),parameter :: lightgoldenrod1      =rgb_color(255,236,139) 
  type(rgb_color),parameter :: lightgoldenrod2      =rgb_color(238,220,130) 
  type(rgb_color),parameter :: lightgoldenrod3      =rgb_color(205,190,112) 
  type(rgb_color),parameter :: lightgoldenrod4      =rgb_color(139,129,76) 
  type(rgb_color),parameter :: lightyellow1         =rgb_color(255,255,224) 
  type(rgb_color),parameter :: lightyellow2         =rgb_color(238,238,209) 
  type(rgb_color),parameter :: lightyellow3         =rgb_color(205,205,180) 
  type(rgb_color),parameter :: lightyellow4         =rgb_color(139,139,122) 
  type(rgb_color),parameter :: yellow1              =rgb_color(255,255,0) 
  type(rgb_color),parameter :: yellow2              =rgb_color(238,238,0) 
  type(rgb_color),parameter :: yellow3              =rgb_color(205,205,0) 
  type(rgb_color),parameter :: yellow4              =rgb_color(139,139,0) 
  type(rgb_color),parameter :: gold1                =rgb_color(255,215,0) 
  type(rgb_color),parameter :: gold2                =rgb_color(238,201,0) 
  type(rgb_color),parameter :: gold3                =rgb_color(205,173,0) 
  type(rgb_color),parameter :: gold4                =rgb_color(139,117,0) 
  type(rgb_color),parameter :: goldenrod1           =rgb_color(255,193,37) 
  type(rgb_color),parameter :: goldenrod2           =rgb_color(238,180,34) 
  type(rgb_color),parameter :: goldenrod3           =rgb_color(205,155,29) 
  type(rgb_color),parameter :: goldenrod4           =rgb_color(139,105,20) 
  type(rgb_color),parameter :: darkgoldenrod1       =rgb_color(255,185,15) 
  type(rgb_color),parameter :: darkgoldenrod2       =rgb_color(238,173,14) 
  type(rgb_color),parameter :: darkgoldenrod3       =rgb_color(205,149,12) 
  type(rgb_color),parameter :: darkgoldenrod4       =rgb_color(139,101,8) 
  type(rgb_color),parameter :: rosybrown1           =rgb_color(255,193,193) 
  type(rgb_color),parameter :: rosybrown2           =rgb_color(238,180,180) 
  type(rgb_color),parameter :: rosybrown3           =rgb_color(205,155,155) 
  type(rgb_color),parameter :: rosybrown4           =rgb_color(139,105,105) 
  type(rgb_color),parameter :: indianred1           =rgb_color(255,106,106) 
  type(rgb_color),parameter :: indianred2           =rgb_color(238,99,99) 
  type(rgb_color),parameter :: indianred3           =rgb_color(205,85,85) 
  type(rgb_color),parameter :: indianred4           =rgb_color(139,58,58) 
  type(rgb_color),parameter :: sienna1              =rgb_color(255,130,71) 
  type(rgb_color),parameter :: sienna2              =rgb_color(238,121,66) 
  type(rgb_color),parameter :: sienna3              =rgb_color(205,104,57) 
  type(rgb_color),parameter :: sienna4              =rgb_color(139,71,38) 
  type(rgb_color),parameter :: burlywood1           =rgb_color(255,211,155) 
  type(rgb_color),parameter :: burlywood2           =rgb_color(238,197,145) 
  type(rgb_color),parameter :: burlywood3           =rgb_color(205,170,125) 
  type(rgb_color),parameter :: burlywood4           =rgb_color(139,115,85) 
  type(rgb_color),parameter :: wheat1               =rgb_color(255,231,186) 
  type(rgb_color),parameter :: wheat2               =rgb_color(238,216,174) 
  type(rgb_color),parameter :: wheat3               =rgb_color(205,186,150) 
  type(rgb_color),parameter :: wheat4               =rgb_color(139,126,102) 
  type(rgb_color),parameter :: tan1                 =rgb_color(255,165,79) 
  type(rgb_color),parameter :: tan2                 =rgb_color(238,154,73) 
  type(rgb_color),parameter :: tan3                 =rgb_color(205,133,63) 
  type(rgb_color),parameter :: tan4                 =rgb_color(139,90,43) 
  type(rgb_color),parameter :: chocolate1           =rgb_color(255,127,36) 
  type(rgb_color),parameter :: chocolate2           =rgb_color(238,118,33) 
  type(rgb_color),parameter :: chocolate3           =rgb_color(205,102,29) 
  type(rgb_color),parameter :: chocolate4           =rgb_color(139,69,19) 
  type(rgb_color),parameter :: firebrick1           =rgb_color(255,48,48) 
  type(rgb_color),parameter :: firebrick2           =rgb_color(238,44,44) 
  type(rgb_color),parameter :: firebrick3           =rgb_color(205,38,38) 
  type(rgb_color),parameter :: firebrick4           =rgb_color(139,26,26) 
  type(rgb_color),parameter :: brown1               =rgb_color(255,64,64) 
  type(rgb_color),parameter :: brown2               =rgb_color(238,59,59) 
  type(rgb_color),parameter :: brown3               =rgb_color(205,51,51) 
  type(rgb_color),parameter :: brown4               =rgb_color(139,35,35) 
  type(rgb_color),parameter :: salmon1              =rgb_color(255,140,105) 
  type(rgb_color),parameter :: salmon2              =rgb_color(238,130,98) 
  type(rgb_color),parameter :: salmon3              =rgb_color(205,112,84) 
  type(rgb_color),parameter :: salmon4              =rgb_color(139,76,57) 
  type(rgb_color),parameter :: lightsalmon1         =rgb_color(255,160,122) 
  type(rgb_color),parameter :: lightsalmon2         =rgb_color(238,149,114) 
  type(rgb_color),parameter :: lightsalmon3         =rgb_color(205,129,98) 
  type(rgb_color),parameter :: lightsalmon4         =rgb_color(139,87,66) 
  type(rgb_color),parameter :: orange1              =rgb_color(255,165,0) 
  type(rgb_color),parameter :: orange2              =rgb_color(238,154,0) 
  type(rgb_color),parameter :: orange3              =rgb_color(205,133,0) 
  type(rgb_color),parameter :: orange4              =rgb_color(139,90,0) 
  type(rgb_color),parameter :: darkorange1          =rgb_color(255,127,0) 
  type(rgb_color),parameter :: darkorange2          =rgb_color(238,118,0) 
  type(rgb_color),parameter :: darkorange3          =rgb_color(205,102,0) 
  type(rgb_color),parameter :: darkorange4          =rgb_color(139,69,0) 
  type(rgb_color),parameter :: coral1               =rgb_color(255,114,86) 
  type(rgb_color),parameter :: coral2               =rgb_color(238,106,80) 
  type(rgb_color),parameter :: coral3               =rgb_color(205,91,69) 
  type(rgb_color),parameter :: coral4               =rgb_color(139,62,47) 
  type(rgb_color),parameter :: tomato1              =rgb_color(255,99,71) 
  type(rgb_color),parameter :: tomato2              =rgb_color(238,92,66) 
  type(rgb_color),parameter :: tomato3              =rgb_color(205,79,57) 
  type(rgb_color),parameter :: tomato4              =rgb_color(139,54,38) 
  type(rgb_color),parameter :: orangered1           =rgb_color(255,69,0) 
  type(rgb_color),parameter :: orangered2           =rgb_color(238,64,0) 
  type(rgb_color),parameter :: orangered3           =rgb_color(205,55,0) 
  type(rgb_color),parameter :: orangered4           =rgb_color(139,37,0) 
  type(rgb_color),parameter :: red1                 =rgb_color(255,0,0) 
  type(rgb_color),parameter :: red2                 =rgb_color(238,0,0) 
  type(rgb_color),parameter :: red3                 =rgb_color(205,0,0) 
  type(rgb_color),parameter :: red4                 =rgb_color(139,0,0) 
  type(rgb_color),parameter :: debianred            =rgb_color(215,7,81) 
  type(rgb_color),parameter :: deeppink1            =rgb_color(255,20,147) 
  type(rgb_color),parameter :: deeppink2            =rgb_color(238,18,137) 
  type(rgb_color),parameter :: deeppink3            =rgb_color(205,16,118) 
  type(rgb_color),parameter :: deeppink4            =rgb_color(139,10,80) 
  type(rgb_color),parameter :: hotpink1             =rgb_color(255,110,180) 
  type(rgb_color),parameter :: hotpink2             =rgb_color(238,106,167) 
  type(rgb_color),parameter :: hotpink3             =rgb_color(205,96,144) 
  type(rgb_color),parameter :: hotpink4             =rgb_color(139,58,98) 
  type(rgb_color),parameter :: pink1                =rgb_color(255,181,197) 
  type(rgb_color),parameter :: pink2                =rgb_color(238,169,184) 
  type(rgb_color),parameter :: pink3                =rgb_color(205,145,158) 
  type(rgb_color),parameter :: pink4                =rgb_color(139,99,108) 
  type(rgb_color),parameter :: lightpink1           =rgb_color(255,174,185) 
  type(rgb_color),parameter :: lightpink2           =rgb_color(238,162,173) 
  type(rgb_color),parameter :: lightpink3           =rgb_color(205,140,149) 
  type(rgb_color),parameter :: lightpink4           =rgb_color(139,95,101) 
  type(rgb_color),parameter :: palevioletred1       =rgb_color(255,130,171) 
  type(rgb_color),parameter :: palevioletred2       =rgb_color(238,121,159) 
  type(rgb_color),parameter :: palevioletred3       =rgb_color(205,104,137) 
  type(rgb_color),parameter :: palevioletred4       =rgb_color(139,71,93) 
  type(rgb_color),parameter :: maroon1              =rgb_color(255,52,179) 
  type(rgb_color),parameter :: maroon2              =rgb_color(238,48,167) 
  type(rgb_color),parameter :: maroon3              =rgb_color(205,41,144) 
  type(rgb_color),parameter :: maroon4              =rgb_color(139,28,98) 
  type(rgb_color),parameter :: violetred1           =rgb_color(255,62,150) 
  type(rgb_color),parameter :: violetred2           =rgb_color(238,58,140) 
  type(rgb_color),parameter :: violetred3           =rgb_color(205,50,120) 
  type(rgb_color),parameter :: violetred4           =rgb_color(139,34,82) 
  type(rgb_color),parameter :: magenta1             =rgb_color(255,0,255) 
  type(rgb_color),parameter :: magenta2             =rgb_color(238,0,238) 
  type(rgb_color),parameter :: magenta3             =rgb_color(205,0,205) 
  type(rgb_color),parameter :: magenta4             =rgb_color(139,0,139) 
  type(rgb_color),parameter :: orchid1              =rgb_color(255,131,250) 
  type(rgb_color),parameter :: orchid2              =rgb_color(238,122,233) 
  type(rgb_color),parameter :: orchid3              =rgb_color(205,105,201) 
  type(rgb_color),parameter :: orchid4              =rgb_color(139,71,137) 
  type(rgb_color),parameter :: plum1                =rgb_color(255,187,255) 
  type(rgb_color),parameter :: plum2                =rgb_color(238,174,238) 
  type(rgb_color),parameter :: plum3                =rgb_color(205,150,205) 
  type(rgb_color),parameter :: plum4                =rgb_color(139,102,139) 
  type(rgb_color),parameter :: mediumorchid1        =rgb_color(224,102,255) 
  type(rgb_color),parameter :: mediumorchid2        =rgb_color(209,95,238) 
  type(rgb_color),parameter :: mediumorchid3        =rgb_color(180,82,205) 
  type(rgb_color),parameter :: mediumorchid4        =rgb_color(122,55,139) 
  type(rgb_color),parameter :: darkorchid1          =rgb_color(191,62,255) 
  type(rgb_color),parameter :: darkorchid2          =rgb_color(178,58,238) 
  type(rgb_color),parameter :: darkorchid3          =rgb_color(154,50,205) 
  type(rgb_color),parameter :: darkorchid4          =rgb_color(104,34,139) 
  type(rgb_color),parameter :: purple1              =rgb_color(155,48,255) 
  type(rgb_color),parameter :: purple2              =rgb_color(145,44,238) 
  type(rgb_color),parameter :: purple3              =rgb_color(125,38,205) 
  type(rgb_color),parameter :: purple4              =rgb_color(85,26,139) 
  type(rgb_color),parameter :: mediumpurple1        =rgb_color(171,130,255) 
  type(rgb_color),parameter :: mediumpurple2        =rgb_color(159,121,238) 
  type(rgb_color),parameter :: mediumpurple3        =rgb_color(137,104,205) 
  type(rgb_color),parameter :: mediumpurple4        =rgb_color(93,71,139) 
  type(rgb_color),parameter :: thistle1             =rgb_color(255,225,255) 
  type(rgb_color),parameter :: thistle2             =rgb_color(238,210,238) 
  type(rgb_color),parameter :: thistle3             =rgb_color(205,181,205) 
  type(rgb_color),parameter :: thistle4             =rgb_color(139,123,139) 
  type(rgb_color),parameter :: gray0                =rgb_color(0,0,0) 
  type(rgb_color),parameter :: grey0                =rgb_color(0,0,0) 
  type(rgb_color),parameter :: gray1                =rgb_color(3,3,3) 
  type(rgb_color),parameter :: grey1                =rgb_color(3,3,3) 
  type(rgb_color),parameter :: gray2                =rgb_color(5,5,5) 
  type(rgb_color),parameter :: grey2                =rgb_color(5,5,5) 
  type(rgb_color),parameter :: gray3                =rgb_color(8,8,8) 
  type(rgb_color),parameter :: grey3                =rgb_color(8,8,8) 
  type(rgb_color),parameter :: gray4                =rgb_color(10,10,10) 
  type(rgb_color),parameter :: grey4                =rgb_color(10,10,10) 
  type(rgb_color),parameter :: gray5                =rgb_color(13,13,13) 
  type(rgb_color),parameter :: grey5                =rgb_color(13,13,13) 
  type(rgb_color),parameter :: gray6                =rgb_color(15,15,15) 
  type(rgb_color),parameter :: grey6                =rgb_color(15,15,15) 
  type(rgb_color),parameter :: gray7                =rgb_color(18,18,18) 
  type(rgb_color),parameter :: grey7                =rgb_color(18,18,18) 
  type(rgb_color),parameter :: gray8                =rgb_color(20,20,20) 
  type(rgb_color),parameter :: grey8                =rgb_color(20,20,20) 
  type(rgb_color),parameter :: gray9                =rgb_color(23,23,23) 
  type(rgb_color),parameter :: grey9                =rgb_color(23,23,23) 
  type(rgb_color),parameter :: gray10               =rgb_color(26,26,26) 
  type(rgb_color),parameter :: grey10               =rgb_color(26,26,26) 
  type(rgb_color),parameter :: gray11               =rgb_color(28,28,28) 
  type(rgb_color),parameter :: grey11               =rgb_color(28,28,28) 
  type(rgb_color),parameter :: gray12               =rgb_color(31,31,31) 
  type(rgb_color),parameter :: grey12               =rgb_color(31,31,31) 
  type(rgb_color),parameter :: gray13               =rgb_color(33,33,33) 
  type(rgb_color),parameter :: grey13               =rgb_color(33,33,33) 
  type(rgb_color),parameter :: gray14               =rgb_color(36,36,36) 
  type(rgb_color),parameter :: grey14               =rgb_color(36,36,36) 
  type(rgb_color),parameter :: gray15               =rgb_color(38,38,38) 
  type(rgb_color),parameter :: grey15               =rgb_color(38,38,38) 
  type(rgb_color),parameter :: gray16               =rgb_color(41,41,41) 
  type(rgb_color),parameter :: grey16               =rgb_color(41,41,41) 
  type(rgb_color),parameter :: gray17               =rgb_color(43,43,43) 
  type(rgb_color),parameter :: grey17               =rgb_color(43,43,43) 
  type(rgb_color),parameter :: gray18               =rgb_color(46,46,46) 
  type(rgb_color),parameter :: grey18               =rgb_color(46,46,46) 
  type(rgb_color),parameter :: gray19               =rgb_color(48,48,48) 
  type(rgb_color),parameter :: grey19               =rgb_color(48,48,48) 
  type(rgb_color),parameter :: gray20               =rgb_color(51,51,51) 
  type(rgb_color),parameter :: grey20               =rgb_color(51,51,51) 
  type(rgb_color),parameter :: gray21               =rgb_color(54,54,54) 
  type(rgb_color),parameter :: grey21               =rgb_color(54,54,54) 
  type(rgb_color),parameter :: gray22               =rgb_color(56,56,56) 
  type(rgb_color),parameter :: grey22               =rgb_color(56,56,56) 
  type(rgb_color),parameter :: gray23               =rgb_color(59,59,59) 
  type(rgb_color),parameter :: grey23               =rgb_color(59,59,59) 
  type(rgb_color),parameter :: gray24               =rgb_color(61,61,61) 
  type(rgb_color),parameter :: grey24               =rgb_color(61,61,61) 
  type(rgb_color),parameter :: gray25               =rgb_color(64,64,64) 
  type(rgb_color),parameter :: grey25               =rgb_color(64,64,64) 
  type(rgb_color),parameter :: gray26               =rgb_color(66,66,66) 
  type(rgb_color),parameter :: grey26               =rgb_color(66,66,66) 
  type(rgb_color),parameter :: gray27               =rgb_color(69,69,69) 
  type(rgb_color),parameter :: grey27               =rgb_color(69,69,69) 
  type(rgb_color),parameter :: gray28               =rgb_color(71,71,71) 
  type(rgb_color),parameter :: grey28               =rgb_color(71,71,71) 
  type(rgb_color),parameter :: gray29               =rgb_color(74,74,74) 
  type(rgb_color),parameter :: grey29               =rgb_color(74,74,74) 
  type(rgb_color),parameter :: gray30               =rgb_color(77,77,77) 
  type(rgb_color),parameter :: grey30               =rgb_color(77,77,77) 
  type(rgb_color),parameter :: gray31               =rgb_color(79,79,79) 
  type(rgb_color),parameter :: grey31               =rgb_color(79,79,79) 
  type(rgb_color),parameter :: gray32               =rgb_color(82,82,82) 
  type(rgb_color),parameter :: grey32               =rgb_color(82,82,82) 
  type(rgb_color),parameter :: gray33               =rgb_color(84,84,84) 
  type(rgb_color),parameter :: grey33               =rgb_color(84,84,84) 
  type(rgb_color),parameter :: gray34               =rgb_color(87,87,87) 
  type(rgb_color),parameter :: grey34               =rgb_color(87,87,87) 
  type(rgb_color),parameter :: gray35               =rgb_color(89,89,89) 
  type(rgb_color),parameter :: grey35               =rgb_color(89,89,89) 
  type(rgb_color),parameter :: gray36               =rgb_color(92,92,92) 
  type(rgb_color),parameter :: grey36               =rgb_color(92,92,92) 
  type(rgb_color),parameter :: gray37               =rgb_color(94,94,94) 
  type(rgb_color),parameter :: grey37               =rgb_color(94,94,94) 
  type(rgb_color),parameter :: gray38               =rgb_color(97,97,97) 
  type(rgb_color),parameter :: grey38               =rgb_color(97,97,97) 
  type(rgb_color),parameter :: gray39               =rgb_color(99,99,99) 
  type(rgb_color),parameter :: grey39               =rgb_color(99,99,99) 
  type(rgb_color),parameter :: gray40               =rgb_color(102,102,102) 
  type(rgb_color),parameter :: grey40               =rgb_color(102,102,102) 
  type(rgb_color),parameter :: gray41               =rgb_color(105,105,105) 
  type(rgb_color),parameter :: grey41               =rgb_color(105,105,105) 
  type(rgb_color),parameter :: gray42               =rgb_color(107,107,107) 
  type(rgb_color),parameter :: grey42               =rgb_color(107,107,107) 
  type(rgb_color),parameter :: gray43               =rgb_color(110,110,110) 
  type(rgb_color),parameter :: grey43               =rgb_color(110,110,110) 
  type(rgb_color),parameter :: gray44               =rgb_color(112,112,112) 
  type(rgb_color),parameter :: grey44               =rgb_color(112,112,112) 
  type(rgb_color),parameter :: gray45               =rgb_color(115,115,115) 
  type(rgb_color),parameter :: grey45               =rgb_color(115,115,115) 
  type(rgb_color),parameter :: gray46               =rgb_color(117,117,117) 
  type(rgb_color),parameter :: grey46               =rgb_color(117,117,117) 
  type(rgb_color),parameter :: gray47               =rgb_color(120,120,120) 
  type(rgb_color),parameter :: grey47               =rgb_color(120,120,120) 
  type(rgb_color),parameter :: gray48               =rgb_color(122,122,122) 
  type(rgb_color),parameter :: grey48               =rgb_color(122,122,122) 
  type(rgb_color),parameter :: gray49               =rgb_color(125,125,125) 
  type(rgb_color),parameter :: grey49               =rgb_color(125,125,125) 
  type(rgb_color),parameter :: gray50               =rgb_color(127,127,127) 
  type(rgb_color),parameter :: grey50               =rgb_color(127,127,127) 
  type(rgb_color),parameter :: gray51               =rgb_color(130,130,130) 
  type(rgb_color),parameter :: grey51               =rgb_color(130,130,130) 
  type(rgb_color),parameter :: gray52               =rgb_color(133,133,133) 
  type(rgb_color),parameter :: grey52               =rgb_color(133,133,133) 
  type(rgb_color),parameter :: gray53               =rgb_color(135,135,135) 
  type(rgb_color),parameter :: grey53               =rgb_color(135,135,135) 
  type(rgb_color),parameter :: gray54               =rgb_color(138,138,138) 
  type(rgb_color),parameter :: grey54               =rgb_color(138,138,138) 
  type(rgb_color),parameter :: gray55               =rgb_color(140,140,140) 
  type(rgb_color),parameter :: grey55               =rgb_color(140,140,140) 
  type(rgb_color),parameter :: gray56               =rgb_color(143,143,143) 
  type(rgb_color),parameter :: grey56               =rgb_color(143,143,143) 
  type(rgb_color),parameter :: gray57               =rgb_color(145,145,145) 
  type(rgb_color),parameter :: grey57               =rgb_color(145,145,145) 
  type(rgb_color),parameter :: gray58               =rgb_color(148,148,148) 
  type(rgb_color),parameter :: grey58               =rgb_color(148,148,148) 
  type(rgb_color),parameter :: gray59               =rgb_color(150,150,150) 
  type(rgb_color),parameter :: grey59               =rgb_color(150,150,150) 
  type(rgb_color),parameter :: gray60               =rgb_color(153,153,153) 
  type(rgb_color),parameter :: grey60               =rgb_color(153,153,153) 
  type(rgb_color),parameter :: gray61               =rgb_color(156,156,156) 
  type(rgb_color),parameter :: grey61               =rgb_color(156,156,156) 
  type(rgb_color),parameter :: gray62               =rgb_color(158,158,158) 
  type(rgb_color),parameter :: grey62               =rgb_color(158,158,158) 
  type(rgb_color),parameter :: gray63               =rgb_color(161,161,161) 
  type(rgb_color),parameter :: grey63               =rgb_color(161,161,161) 
  type(rgb_color),parameter :: gray64               =rgb_color(163,163,163) 
  type(rgb_color),parameter :: grey64               =rgb_color(163,163,163) 
  type(rgb_color),parameter :: gray65               =rgb_color(166,166,166) 
  type(rgb_color),parameter :: grey65               =rgb_color(166,166,166) 
  type(rgb_color),parameter :: gray66               =rgb_color(168,168,168) 
  type(rgb_color),parameter :: grey66               =rgb_color(168,168,168) 
  type(rgb_color),parameter :: gray67               =rgb_color(171,171,171) 
  type(rgb_color),parameter :: grey67               =rgb_color(171,171,171) 
  type(rgb_color),parameter :: gray68               =rgb_color(173,173,173) 
  type(rgb_color),parameter :: grey68               =rgb_color(173,173,173) 
  type(rgb_color),parameter :: gray69               =rgb_color(176,176,176) 
  type(rgb_color),parameter :: grey69               =rgb_color(176,176,176) 
  type(rgb_color),parameter :: gray70               =rgb_color(179,179,179) 
  type(rgb_color),parameter :: grey70               =rgb_color(179,179,179) 
  type(rgb_color),parameter :: gray71               =rgb_color(181,181,181) 
  type(rgb_color),parameter :: grey71               =rgb_color(181,181,181) 
  type(rgb_color),parameter :: gray72               =rgb_color(184,184,184) 
  type(rgb_color),parameter :: grey72               =rgb_color(184,184,184) 
  type(rgb_color),parameter :: gray73               =rgb_color(186,186,186) 
  type(rgb_color),parameter :: grey73               =rgb_color(186,186,186) 
  type(rgb_color),parameter :: gray74               =rgb_color(189,189,189) 
  type(rgb_color),parameter :: grey74               =rgb_color(189,189,189) 
  type(rgb_color),parameter :: gray75               =rgb_color(191,191,191) 
  type(rgb_color),parameter :: grey75               =rgb_color(191,191,191) 
  type(rgb_color),parameter :: gray76               =rgb_color(194,194,194) 
  type(rgb_color),parameter :: grey76               =rgb_color(194,194,194) 
  type(rgb_color),parameter :: gray77               =rgb_color(196,196,196) 
  type(rgb_color),parameter :: grey77               =rgb_color(196,196,196) 
  type(rgb_color),parameter :: gray78               =rgb_color(199,199,199) 
  type(rgb_color),parameter :: grey78               =rgb_color(199,199,199) 
  type(rgb_color),parameter :: gray79               =rgb_color(201,201,201) 
  type(rgb_color),parameter :: grey79               =rgb_color(201,201,201) 
  type(rgb_color),parameter :: gray80               =rgb_color(204,204,204) 
  type(rgb_color),parameter :: grey80               =rgb_color(204,204,204) 
  type(rgb_color),parameter :: gray81               =rgb_color(207,207,207) 
  type(rgb_color),parameter :: grey81               =rgb_color(207,207,207) 
  type(rgb_color),parameter :: gray82               =rgb_color(209,209,209) 
  type(rgb_color),parameter :: grey82               =rgb_color(209,209,209) 
  type(rgb_color),parameter :: gray83               =rgb_color(212,212,212) 
  type(rgb_color),parameter :: grey83               =rgb_color(212,212,212) 
  type(rgb_color),parameter :: gray84               =rgb_color(214,214,214) 
  type(rgb_color),parameter :: grey84               =rgb_color(214,214,214) 
  type(rgb_color),parameter :: gray85               =rgb_color(217,217,217) 
  type(rgb_color),parameter :: grey85               =rgb_color(217,217,217) 
  type(rgb_color),parameter :: gray86               =rgb_color(219,219,219) 
  type(rgb_color),parameter :: grey86               =rgb_color(219,219,219) 
  type(rgb_color),parameter :: gray87               =rgb_color(222,222,222) 
  type(rgb_color),parameter :: grey87               =rgb_color(222,222,222) 
  type(rgb_color),parameter :: gray88               =rgb_color(224,224,224) 
  type(rgb_color),parameter :: grey88               =rgb_color(224,224,224) 
  type(rgb_color),parameter :: gray89               =rgb_color(227,227,227) 
  type(rgb_color),parameter :: grey89               =rgb_color(227,227,227) 
  type(rgb_color),parameter :: gray90               =rgb_color(229,229,229) 
  type(rgb_color),parameter :: grey90               =rgb_color(229,229,229) 
  type(rgb_color),parameter :: gray91               =rgb_color(232,232,232) 
  type(rgb_color),parameter :: grey91               =rgb_color(232,232,232) 
  type(rgb_color),parameter :: gray92               =rgb_color(235,235,235) 
  type(rgb_color),parameter :: grey92               =rgb_color(235,235,235) 
  type(rgb_color),parameter :: gray93               =rgb_color(237,237,237) 
  type(rgb_color),parameter :: grey93               =rgb_color(237,237,237) 
  type(rgb_color),parameter :: gray94               =rgb_color(240,240,240) 
  type(rgb_color),parameter :: grey94               =rgb_color(240,240,240) 
  type(rgb_color),parameter :: gray95               =rgb_color(242,242,242) 
  type(rgb_color),parameter :: grey95               =rgb_color(242,242,242) 
  type(rgb_color),parameter :: gray96               =rgb_color(245,245,245) 
  type(rgb_color),parameter :: grey96               =rgb_color(245,245,245) 
  type(rgb_color),parameter :: gray97               =rgb_color(247,247,247) 
  type(rgb_color),parameter :: grey97               =rgb_color(247,247,247) 
  type(rgb_color),parameter :: gray98               =rgb_color(250,250,250) 
  type(rgb_color),parameter :: grey98               =rgb_color(250,250,250) 
  type(rgb_color),parameter :: gray99               =rgb_color(252,252,252) 
  type(rgb_color),parameter :: grey99               =rgb_color(252,252,252) 
  type(rgb_color),parameter :: gray100              =rgb_color(255,255,255) 
  type(rgb_color),parameter :: grey100              =rgb_color(255,255,255) 
  type(rgb_color),parameter :: darkgrey             =rgb_color(169,169,169) 
  type(rgb_color),parameter :: darkgray             =rgb_color(169,169,169) 
  type(rgb_color),parameter :: darkblue             =rgb_color(0,0,139) 
  type(rgb_color),parameter :: darkcyan             =rgb_color(0,139,139) 
  type(rgb_color),parameter :: darkmagenta          =rgb_color(139,0,139) 
  type(rgb_color),parameter :: darkred              =rgb_color(139,0,0) 
  type(rgb_color),parameter :: lightgreen           =rgb_color(144,238,144) 

  interface operator(+)
     module procedure add_colors
  end interface operator(+)
  interface assignment(=)
     module procedure equal_colors
  end interface assignment(=)
  interface operator(-)
     module procedure subtract_colors
  end interface operator(-)
  interface operator(*)
     module procedure scalar_left_color
  end interface operator(*)
  interface operator(.dot.)
     module procedure dot_scalar_colors
  end interface operator(.dot.)

contains

  function rgb(c) result(num)
    type(rgb_color),intent(in) :: c
    integer :: num
    num = int(c%r)*65536 + int(c%g)*256 + int(c%b)
  end function rgb

  elemental subroutine equal_colors(C1,C2)
    type(rgb_color),intent(in)    :: C2
    type(rgb_color),intent(inout) :: C1
    C1%r = C2%r
    C1%g = C2%g
    C1%b = C2%b
  end subroutine equal_colors

  elemental function add_colors(c1,c2) result(c)
    type(rgb_color),intent(in) :: c1,c2
    type(rgb_color)            :: c
    c%r = c1%r + c2%r
    c%g = c1%g + c2%g
    c%b = c1%b + c2%b
  end function add_colors

  elemental function subtract_colors(c1,c2) result(c)
    type(rgb_color),intent(in) :: c1,c2
    type(rgb_color)            :: c
    c%r = c1%r - c2%r
    c%g = c1%g - c2%g
    c%b = c1%b - c2%b
  end function subtract_colors

  elemental function scalar_left_color(k,cin) result(cout)
    real(8),intent(in)         :: k
    type(rgb_color),intent(in) :: cin
    type(rgb_color)            :: cout
    cout%r = k*cin%r
    cout%g = k*cin%g
    cout%b = k*cin%b
  end function scalar_left_color

  elemental function scalar_right_color(k,cin) result(cout)
    real(8),intent(in)         :: k
    type(rgb_color),intent(in) :: cin
    type(rgb_color)            :: cout
    cout%r = cin%r*k
    cout%g = cin%g*k
    cout%b = cin%b*k
  end function scalar_right_color

  function dot_scalar_colors(v,cin) result(cout)
    real(8),dimension(:),intent(in)               :: v
    type(rgb_color),dimension(size(v)),intent(in) :: cin
    type(rgb_color)                               :: cout
    integer :: i
    cout=rgb_color(0,0,0)
    do i=1,size(v)
       cout%r = cout%r + v(i)*cin(i)%r
       cout%g = cout%g + v(i)*cin(i)%g
       cout%b = cout%b + v(i)*cin(i)%b
    enddo
  end function dot_scalar_colors

  function pick_color(string) result(crgb)
    character(len=*) :: string
    type(rgb_color)  :: crgb
    select case(string)
    case("black")
       crgb=black
    case("red")
       crgb=red
    case("green")
       crgb=green
    case("orange")
       crgb=orange
    case("blue")
       crgb=blue
    case("yellow")
       crgb=yellow
    case("cyan")
       crgb=cyan
    case("magenta")
       crgb=magenta
    case default
       stop "pick_color: color name does not exist"
    end select
  end function pick_color
END MODULE SF_COLORS





MODULE SF_PAULI
  implicit none
  private

  complex(8),parameter :: zero=(0.d0,0.d0)
  complex(8),parameter :: xi=(0.d0,1.d0)
  complex(8),parameter :: one=(1.d0,0.d0)

  complex(8),dimension(2,2),parameter :: pauli_0=reshape([one,zero,zero,one],[2,2])
  complex(8),dimension(2,2),parameter :: pauli_x=reshape([zero,one,one,zero],[2,2])
  complex(8),dimension(2,2),parameter :: pauli_y=reshape([zero,xi,-xi,zero],[2,2])
  complex(8),dimension(2,2),parameter :: pauli_z=reshape([one,zero,zero,-one],[2,2])
  !
  complex(8),dimension(2,2),parameter :: pauli_1=pauli_x
  complex(8),dimension(2,2),parameter :: pauli_2=pauli_y
  complex(8),dimension(2,2),parameter :: pauli_3=pauli_z
  !
  complex(8),dimension(2,2),parameter :: pauli_tau_0=pauli_0
  complex(8),dimension(2,2),parameter :: pauli_tau_x=pauli_x
  complex(8),dimension(2,2),parameter :: pauli_tau_y=pauli_y
  complex(8),dimension(2,2),parameter :: pauli_tau_z=pauli_z
  !
  complex(8),dimension(2,2),parameter :: pauli_tau_1=pauli_tau_x
  complex(8),dimension(2,2),parameter :: pauli_tau_2=pauli_tau_y
  complex(8),dimension(2,2),parameter :: pauli_tau_3=pauli_tau_z
  !
  complex(8),dimension(2,2),parameter :: pauli_sigma_0=pauli_tau_0
  complex(8),dimension(2,2),parameter :: pauli_sigma_x=pauli_tau_x
  complex(8),dimension(2,2),parameter :: pauli_sigma_y=pauli_tau_y
  complex(8),dimension(2,2),parameter :: pauli_sigma_z=pauli_tau_z
  !
  complex(8),dimension(2,2),parameter :: pauli_sigma_1=pauli_tau_x
  complex(8),dimension(2,2),parameter :: pauli_sigma_2=pauli_tau_y
  complex(8),dimension(2,2),parameter :: pauli_sigma_3=pauli_tau_z


  public :: pauli_0
  public :: pauli_x
  public :: pauli_y
  public :: pauli_z
  !
  public :: pauli_1
  public :: pauli_2
  public :: pauli_3
  !
  public :: pauli_tau_0
  public :: pauli_tau_x
  public :: pauli_tau_y
  public :: pauli_tau_z
  !
  public :: pauli_tau_1
  public :: pauli_tau_2
  public :: pauli_tau_3
  !
  public :: pauli_sigma_0
  public :: pauli_sigma_x
  public :: pauli_sigma_y
  public :: pauli_sigma_z
  !
  public :: pauli_sigma_1
  public :: pauli_sigma_2
  public :: pauli_sigma_3


  public :: kronecker_product_pauli_matrices
  public :: kronecker_product_pauli_vector

contains


  !---------------------------------------------------------------------
  !PURPOSE: return the Kronecker's product of 2 Pauli's matrices. 
  !---------------------------------------------------------------------
  function kronecker_product_pauli_matrices(sigma1,sigma2) result(gamma)
    complex(8),dimension(2,2),intent(in) :: sigma1,sigma2
    complex(8),dimension(4,4)            :: gamma
    gamma = c_kronecker_product(sigma1,2,2,sigma2,2,2)
  end function kronecker_product_pauli_matrices



  !---------------------------------------------------------------------
  !PURPOSE: return the Kronecker's product of n Pauli's matrices. 
  ! The especification of the order of the matrices is given as input on the 
  ! vector vec_ord_pm, that has dimension npm.
  !---------------------------------------------------------------------
  function kronecker_product_pauli_vector(vec_ord_pm, npm) result(kron_prod_n_pauli_mat)
    integer, intent(in)     :: npm
    integer, intent(in)     :: vec_ord_pm(npm)
    complex(8)              :: kron_prod_n_pauli_mat(2**npm,2**npm)
    integer                 :: d2
    complex(8)              :: M2(2,2)
    complex(8), allocatable :: M1(:,:), M1_kp_M2(:,:)
    integer                 :: d1, i
    d2=2
    do i=1,npm-1
       select case(vec_ord_pm(i+1))
       case (0)
          M2 = pauli_sigma_0
       case (1)
          M2 = pauli_sigma_1
       case (2)
          M2 = pauli_sigma_2
       case (3)
          M2 = pauli_sigma_3
       end select
       d1 = 2**i
       if(i==1) then
          allocate(M1(d1,d1))
          select case(vec_ord_pm(i))
          case (0) 
             M1 = pauli_sigma_0
          case (1) 
             M1 = pauli_sigma_1
          case (2) 
             M1 = pauli_sigma_2
          case (3) 
             M1 = pauli_sigma_3
          end select
       endif
       allocate(M1_kp_M2(d1*d2,d1*d2))
       M1_kp_M2 = c_kronecker_product(M1,d1,d1,M2,d2,d2)  
       deallocate(M1)
       allocate(M1(1:d1*d2,1:d1*d2))
       M1 = M1_kp_M2
       deallocate(M1_kp_M2)
    end do
    kron_prod_n_pauli_mat = M1
    deallocate(M1)
  end function kronecker_product_pauli_vector




  ! !---------------------------------------------------------------------
  ! !PURPOSE: Function to compute the tensor product (M1_kp_M2) of 
  ! ! two complex matrices M1 and M2. nr1(nr2) and nc1(nc2) are 
  ! ! the number of rows and columns of the Matrix M1 and M2
  ! !---------------------------------------------------------------------
  function c_kronecker_product(M1, nr1, nc1, M2, nr2, nc2) result(M1_kp_M2)
    integer               :: i, j
    integer,intent(in)    :: nr1,nc1,nr2,nc2
    complex(8),intent(in) :: M1(nr1,nc1), M2(nr2,nc2)
    complex(8)            :: M1_kp_M2(nr1*nr2,nc1*nc2)
    M1_kp_M2 = zero
    forall(i =1:nr1,j=1:nc1)
       M1_kp_M2(nr2*(i-1)+1 : nr2*i , nc2*(j-1)+1 : nc2*j)  =  M1(i,j)*M2
    end forall
  end function c_kronecker_product

END MODULE SF_PAULI







module SF_CONSTANTS
  USE SF_MPI_VARS
  USE SF_OMP_VARS
  USE SF_COLORS
  USE SF_PAULI
  implicit none

  !PARAMETERS
  !===============================================================
  complex(8),parameter,public :: zero=(0.d0,0.d0)
  complex(8),parameter,public :: xi=(0.d0,1.d0)
  complex(8),parameter,public :: one=(1.d0,0.d0)
  real(8),parameter,public    :: sqrt2 = 1.41421356237309504880169d0
  real(8),parameter,public    :: sqrt3 = 1.73205080756887729352745d0
  real(8),parameter,public    :: sqrt6 = 2.44948974278317809819728d0
  real(8),parameter,public    :: pi    = 3.14159265358979323846264338327950288419716939937510d0
  real(8),parameter,public    :: pi2   = 6.28318530717959d0
  real(8),parameter,public    :: gamma_euler = 0.57721566490153286060d0  !euler s constant
  real(8),parameter,public    :: euler= 2.7182818284590452353602874713526624977572470936999595749669676277240766303535d0
  integer,parameter,public    :: max_int  = huge(1) 
  real(8),parameter,public    :: max_real = huge(1.d0)
  real(8),parameter,public    :: epsilonr=epsilon(1.d0),epsilonq=1.d-30
  integer,parameter,public    :: dbl=8,dp=8        ! "double" precision
  integer,parameter,public    :: ddp=16            ! "quad"   precision
  integer,parameter,public    :: sp = kind(1.0)    ! "single" precision


  real(8),parameter,public ::                                 Avogadro_constant=  0.602214129000D+24
  real(8),parameter,public ::                                     Bohr_magneton=  0.927400968000D-23
  real(8),parameter,public ::                             Bohr_magneton_in_eVoT=  0.578838180660D-04
  real(8),parameter,public ::                             Bohr_magneton_in_HzoT=  0.139962455500D+11
  real(8),parameter,public ::         Bohr_magneton_in_inverse_meters_per_tesla=     46.686449800000
  real(8),parameter,public ::                              Bohr_magneton_in_KoT=      0.671713880000
  real(8),parameter,public ::                                       Bohr_radius=  0.529177210920D-10
  real(8),parameter,public ::                                Boltzmann_constant=  0.138064880000D-22
  real(8),parameter,public ::                        Boltzmann_constant_in_eVoK=  0.861733240000D-04
  real(8),parameter,public ::                        Boltzmann_constant_in_HzoK=  0.208366180000D+11
  real(8),parameter,public ::   Boltzmann_constant_in_inverse_meters_per_kelvin=     69.503476000000
  real(8),parameter,public ::                                Compton_wavelength=  0.242631023890D-11
  real(8),parameter,public ::                      Compton_wavelength_over_2_pi=  0.386159268000D-12
  real(8),parameter,public ::                                 electric_constant=  0.885418781700D-11
  real(8),parameter,public ::                  electron_charge_to_mass_quotient= -0.175882008800D+12
  real(8),parameter,public ::                                 electron_g_factor= -0.200231930436D+01
  real(8),parameter,public ::                           electron_gyromag__ratio=  0.176085970800D+12
  real(8),parameter,public ::                 electron_gyromag__ratio_over_2_pi=  0.280249526600D+05
  real(8),parameter,public ::                                electron_mag__mom_= -0.928476430000D-23
  real(8),parameter,public ::         electron_mag__mom__to_Bohr_magneton_ratio= -0.100115965218D+01
  real(8),parameter,public ::                                     electron_mass=  0.910938291000D-30
  real(8),parameter,public ::                   electron_mass_energy_equivalent=  0.818710506000D-13
  real(8),parameter,public ::            electron_mass_energy_equivalent_in_MeV=      0.510998928000
  real(8),parameter,public ::                                     electron_volt=  0.160217656500D-18
  real(8),parameter,public ::       electron_volt_atomic_mass_unit_relationship=  0.107354415000D-08
  real(8),parameter,public ::                electron_volt_hartree_relationship=      0.036749323790
  real(8),parameter,public ::                  electron_volt_hertz_relationship=  0.241798934800D+15
  real(8),parameter,public ::          electron_volt_inverse_meter_relationship=  0.806554429000D+06
  real(8),parameter,public ::                  electron_volt_joule_relationship=  0.160217656500D-18
  real(8),parameter,public ::                 electron_volt_kelvin_relationship=  0.116045190000D+05
  real(8),parameter,public ::               electron_volt_kilogram_relationship=  0.178266184500D-35
  real(8),parameter,public ::                                 elementary_charge=  0.160217656500D-18
  real(8),parameter,public ::                          elementary_charge_over_h=  0.241798934800D+15
  real(8),parameter,public ::                                  Faraday_constant=  0.964853365000D+05
  real(8),parameter,public ::Faraday_constant_for_conventional_electric_current=  0.964853321000D+05
  real(8),parameter,public ::                           fine_structure_constant=  0.729735256980D-02
  real(8),parameter,public ::                                Josephson_constant=  0.483597870000D+15
  real(8),parameter,public ::                  joule_electron_volt_relationship=  0.624150934000D+19
  real(8),parameter,public ::                          joule_hertz_relationship=  0.150919031100D+34
  real(8),parameter,public ::                  joule_inverse_meter_relationship=  0.503411701000D+25
  real(8),parameter,public ::                         joule_kelvin_relationship=  0.724297160000D+23
  real(8),parameter,public ::                       joule_kilogram_relationship=  0.111265005600D-16
  real(8),parameter,public ::              kelvin_atomic_mass_unit_relationship=  0.925108680000D-13
  real(8),parameter,public ::                 kelvin_electron_volt_relationship=  0.861733240000D-04
  real(8),parameter,public ::                       kelvin_hartree_relationship=  0.316681140000D-05
  real(8),parameter,public ::                         kelvin_hertz_relationship=  0.208366180000D+11
  real(8),parameter,public ::                 kelvin_inverse_meter_relationship=     69.503476000000
  real(8),parameter,public ::                         kelvin_joule_relationship=  0.138064880000D-22
  real(8),parameter,public ::                      kelvin_kilogram_relationship=  0.153617900000D-39
  real(8),parameter,public ::            kilogram_atomic_mass_unit_relationship=  0.602214129000D+27
  real(8),parameter,public ::               kilogram_electron_volt_relationship=  0.560958885000D+36
  real(8),parameter,public ::                     kilogram_hartree_relationship=  0.206148596800D+35
  real(8),parameter,public ::                       kilogram_hertz_relationship=  0.135639260800D+50
  real(8),parameter,public ::               kilogram_inverse_meter_relationship=  0.452443873000D+42
  real(8),parameter,public ::                       kilogram_joule_relationship=  0.898755178700D+17
  real(8),parameter,public ::                      kilogram_kelvin_relationship=  0.650965820000D+40
  real(8),parameter,public ::                      lattice_parameter_of_silicon=  0.543102050400D-09
  real(8),parameter,public ::                            natural_unit_of_action=  0.105457172600D-33
  real(8),parameter,public ::                    natural_unit_of_action_in_eV_s=  0.658211928000D-15
  real(8),parameter,public ::                            natural_unit_of_energy=  0.818710506000D-13
  real(8),parameter,public ::                     natural_unit_of_energy_in_MeV=      0.510998928000
  real(8),parameter,public ::                            natural_unit_of_length=  0.386159268000D-12
  real(8),parameter,public ::                              natural_unit_of_mass=  0.910938291000D-30
  real(8),parameter,public ::                            natural_unit_of_mom_um=  0.273092429000D-21
  real(8),parameter,public ::                   natural_unit_of_mom_um_in_MeVoc=      0.510998928000
  real(8),parameter,public ::                              natural_unit_of_time=  0.128808866833D-20
  real(8),parameter,public ::                          natural_unit_of_velocity=  0.299792458000D+09
  real(8),parameter,public ::                 Newtonian_constant_of_gravitation=  0.667384000000D-10
  real(8),parameter,public ::                                   Planck_constant=  0.662606957000D-33
  real(8),parameter,public ::                           Planck_constant_in_eV_s=  0.413566751600D-14
  real(8),parameter,public ::                         Planck_constant_over_2_pi=  0.105457172600D-33
  real(8),parameter,public ::                                  Rydberg_constant=  0.109737315685D+08
  real(8),parameter,public ::                    Rydberg_constant_times_c_in_Hz=  0.328984196036D+16
  real(8),parameter,public ::                   Rydberg_constant_times_hc_in_eV=     13.605692530000
  real(8),parameter,public ::                    Rydberg_constant_times_hc_in_J=  0.217987217100D-17
  real(8),parameter,public ::                          speed_of_light_in_vacuum=  0.299792458000D+09
  real(8),parameter,public ::                  standard_acceleration_of_gravity=      9.806650000000
  real(8),parameter,public ::                         Stefan_Boltzmann_constant=  0.567037300000D-07





  public :: timestamp



contains



  !+-------------------------------------------------------------------+
  !PURPOSE  : prints the current YMDHMS date as a time stamp.
  ! Example: 31 May 2001   9:45:54.872 AM
  !+-------------------------------------------------------------------+
  subroutine timestamp(unit)
    integer,optional        :: unit
    integer                 :: unit_
    integer(4),dimension(8) :: data
    unit_=6;if(present(unit))unit_=unit
    call date_and_time(values=data)
    call print_date(data,unit_)
  end subroutine timestamp



  !+-------------------------------------------------------------------+
  !PURPOSE  : print actual date
  !+-------------------------------------------------------------------+
  subroutine print_date(dummy,unit)
    integer(4),dimension(8) :: dummy
    integer                 :: unit
    integer(4)                          :: year
    integer(4)                          :: mese
    integer(4)                          :: day
    integer(4)                          :: h
    integer(4)                          :: m
    integer(4)                          :: s
    integer(4)                          :: ms
    character(len=9),parameter,dimension(12) :: month = (/ &
         'January  ', 'February ', 'March    ', 'April    ', &
         'May      ', 'June     ', 'July     ', 'August   ', &
         'September', 'October  ', 'November ', 'December ' /)
    year = dummy(1)
    mese = dummy(2)
    day  = dummy(3)
    h    = dummy(5)
    m    = dummy(6)
    s    = dummy(7)
    ms   = dummy(8)
    write(unit,"(A,i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3)")&
         "Timestamp: +",day,trim(month(mese)),year, h,':',m,':',s,'.',ms
    write(unit,*)""
  end subroutine print_date

END MODULE SF_CONSTANTS


