module CONSTANTS
  implicit none
  private

  type,public:: constant
     character(len=64) :: name
     real(8)           :: val
     real(8)           :: uncert
     character(len=16) :: units
  end type constant

  integer,parameter,public :: nmax_constant=334
  !constants updated to 2013
  !check this website:
  !http://physics.nist.gov/cuu/Constants/index.html
  type(constant),dimension(nmax_constant),public :: fundamental_constant=(/& 
       constant('lattice spacing of silicon                             ',    192.015571400000D-12,      3.200000000000D-18,'    m               '),&
       constant('alpha particle-electron mass ratio                     ',      7.294299536100D+03,      2.900000000000D-06,'                    '),&
       constant('alpha particle mass                                    ',      6.644656750000D-27,    290.000000000000D-36,'    kg              '),&
       constant('alpha particle mass energy equivalent                  ',    597.191967000000D-12,     26.000000000000D-18,'    J               '),&
       constant('alpha particle mass energy equivalent in MeV           ',      3.727379240000D+03,     82.000000000000D-06,'    MeV             '),&
       constant('alpha particle mass in u                               ',      4.001506179125D+00,     62.000000000000D-12,'    u               '),&
       constant('alpha particle molar mass                              ',      4.001506179125D-03,     62.000000000000D-15,'    kg mol^-1       '),&
       constant('alpha particle-proton mass ratio                       ',      3.972599689330D+00,    360.000000000000D-12,'                    '),&
       constant('Angstrom star                                          ',    100.001495000000D-12,     90.000000000000D-18,'    m               '),&
       constant('atomic mass constant                                   ',      1.660538921000D-27,     73.000000000000D-36,'    kg              '),&
       constant('atomic mass constant energy equivalent                 ',    149.241795400000D-12,      6.600000000000D-18,'    J               '),&
       constant('atomic mass constant energy equivalent in MeV          ',    931.494061000000D+00,     21.000000000000D-06,'    MeV             '),&
       constant('atomic mass unit-electron volt relationship            ',    931.494061000000D+06,     21.000000000000D+00,'    eV              '),&
       constant('atomic mass unit-hartree relationship                  ',     34.231776845000D+06,     24.000000000000D-03,'    E_h             '),&
       constant('atomic mass unit-hertz relationship                    ',    225.234271680000D+21,    160.000000000000D+12,'    Hz              '),&
       constant('atomic mass unit-inverse meter relationship            ',    751.300660420000D+12,    530.000000000000D+03,'    m^-1            '),&
       constant('atomic mass unit-joule relationship                    ',    149.241795400000D-12,      6.600000000000D-18,'    J               '),&
       constant('atomic mass unit-kelvin relationship                   ',     10.809540800000D+12,      9.800000000000D+06,'    K               '),&
       constant('atomic mass unit-kilogram relationship                 ',      1.660538921000D-27,     73.000000000000D-36,'    kg              '),&
       constant('atomic unit of 1st hyperpolarizability                 ',     32.063614490000D-54,    710.000000000000D-63,'    C^3 m^3 J^-2    '),&
       constant('atomic unit of 2nd hyperpolarizability                 ',     62.353805400000D-66,      2.800000000000D-72,'    C^4 m^4 J^-3    '),&
       constant('atomic unit of action                                  ',    105.457172600000D-36,      4.700000000000D-42,'    J s             '),&
       constant('atomic unit of charge                                  ',    160.217656500000D-21,      3.500000000000D-27,'    C               '),&
       constant('atomic unit of charge density                          ',      1.081202338000D+12,     24.000000000000D+03,'    C m^-3          '),&
       constant('atomic unit of current                                 ',      6.623617950000D-03,    150.000000000000D-12,'    A               '),&
       constant('atomic unit of electric dipole mom.                    ',      8.478353260000D-30,    190.000000000000D-39,'    C m             '),&
       constant('atomic unit of electric field                          ',    514.220652000000D+09,     11.000000000000D+03,'    V m^-1          '),&
       constant('atomic unit of electric field gradient                 ',      9.717362000000D+21,    210.000000000000D+12,'    V m^-2          '),&
       constant('atomic unit of electric polarizability                 ',     16.487772754000D-42,     16.000000000000D-51,'    C^2 m^2 J^-1    '),&
       constant('atomic unit of electric potential                      ',     27.211385050000D+00,    600.000000000000D-09,'    V               '),&
       constant('atomic unit of electric quadrupole mom.                ',    448.655133100000D-42,      9.900000000000D-48,'    C m^2           '),&
       constant('atomic unit of energy                                  ',      4.359744340000D-18,    190.000000000000D-27,'    J               '),&
       constant('atomic unit of force                                   ',     82.387227800000D-09,      3.600000000000D-15,'    N               '),&
       constant('atomic unit of length                                  ',     52.917721092000D-12,     17.000000000000D-21,'    m               '),&
       constant('atomic unit of mag. dipole mom.                        ',     18.548019360000D-24,    410.000000000000D-33,'    J T^-1          '),&
       constant('atomic unit of mag. flux density                       ',    235.051746400000D+03,      5.200000000000D-03,'    T               '),&
       constant('atomic unit of magnetizability                         ',     78.910366070000D-30,    130.000000000000D-39,'    J T^-2          '),&
       constant('atomic unit of mass                                    ',    910.938291000000D-33,     40.000000000000D-39,'    kg              '),&
       constant('atomic unit of mom.um                                  ',      1.992851740000D-24,     88.000000000000D-33,'    kg m s^-1       '),&
       constant('atomic unit of permittivity                            ',    111.265005600000D-12,      0.000000000000D+00,'    F m^-1          '),&
       constant('atomic unit of time                                    ',     24.188843265020D-18,    120.000000000000D-30,'    s               '),&
       constant('atomic unit of velocity                                ',      2.187691263790D+06,    710.000000000000D-06,'    m s^-1          '),&
       constant('Avogadro constant                                      ',    602.214129000000D+21,     27.000000000000D+15,'    mol^-1          '),&
       constant('Bohr magneton                                          ',      9.274009680000D-24,    200.000000000000D-33,'    J T^-1          '),&
       constant('Bohr magneton in eV/T                                  ',     57.883818066000D-06,     38.000000000000D-15,'    eV T^-1         '),&
       constant('Bohr magneton in Hz/T                                  ',     13.996245550000D+09,    310.000000000000D+00,'    Hz T^-1         '),&
       constant('Bohr magneton in inverse meters per tesla              ',     46.686449800000D+00,      1.000000000000D-06,'    m^-1 T^-1       '),&
       constant('Bohr magneton in K/T                                   ',    671.713880000000D-03,    610.000000000000D-09,'    K T^-1          '),&
       constant('Bohr radius                                            ',     52.917721092000D-12,     17.000000000000D-21,'    m               '),&
       constant('Boltzmann constant                                     ',     13.806488000000D-24,     13.000000000000D-30,'    J K^-1          '),&
       constant('Boltzmann constant in eV/K                             ',     86.173324000000D-06,     78.000000000000D-12,'    eV K^-1         '),&
       constant('Boltzmann constant in Hz/K                             ',     20.836618000000D+09,     19.000000000000D+03,'    Hz K^-1         '),&
       constant('Boltzmann constant in inverse meters per kelvin        ',     69.503476000000D+00,     63.000000000000D-06,'    m^-1 K^-1       '),&
       constant('characteristic impedance of vacuum                     ',    376.730313461000D+00,      0.000000000000D+00,'    ohm             '),&
       constant('classical electron radius                              ',      2.817940326700D-15,      2.700000000000D-24,'    m               '),&
       constant('Compton wavelength                                     ',      2.426310238900D-12,      1.600000000000D-21,'    m               '),&
       constant('Compton wavelength over 2 pi                           ',    386.159268000000D-15,    250.000000000000D-24,'    m               '),&
       constant('conductance quantum                                    ',     77.480917346000D-06,     25.000000000000D-15,'    S               '),&
       constant('conventional value of Josephson constant               ',    483.597900000000D+12,      0.000000000000D+00,'    Hz V^-1         '),&
       constant('conventional value of von Klitzing constant            ',     25.812807000000D+03,      0.000000000000D+00,'    ohm             '),&
       constant('Cu x unit                                              ',    100.207697000000D-15,     28.000000000000D-21,'    m               '),&
       constant('deuteron-electron mag. mom. ratio                      ',   -466.434553700000D-06,      3.900000000000D-12,'                    '),&
       constant('deuteron-electron mass ratio                           ',      3.670482965200D+03,      1.500000000000D-06,'                    '),&
       constant('deuteron g factor                                      ',    857.438230800000D-03,      7.200000000000D-09,'                    '),&
       constant('deuteron mag. mom.                                     ',      4.330734890000D-27,    100.000000000000D-36,'    J T^-1          '),&
       constant('deuteron mag. mom. to Bohr magneton ratio              ',    466.975455600000D-06,      3.900000000000D-12,'                    '),&
       constant('deuteron mag. mom. to nuclear magneton ratio           ',    857.438230800000D-03,      7.200000000000D-09,'                    '),&
       constant('deuteron mass                                          ',      3.343583480000D-27,    150.000000000000D-36,'    kg              '),&
       constant('deuteron mass energy equivalent                        ',    300.506297000000D-12,     13.000000000000D-18,'    J               '),&
       constant('deuteron mass energy equivalent in MeV                 ',      1.875612859000D+03,     41.000000000000D-06,'    MeV             '),&
       constant('deuteron mass in u                                     ',      2.013553212712D+00,     77.000000000000D-12,'    u               '),&
       constant('deuteron molar mass                                    ',      2.013553212712D-03,     77.000000000000D-15,'    kg mol^-1       '),&
       constant('deuteron-neutron mag. mom. ratio                       ',   -448.206520000000D-03,    110.000000000000D-09,'                    '),&
       constant('deuteron-proton mag. mom. ratio                        ',    307.012207000000D-03,      2.400000000000D-09,'                    '),&
       constant('deuteron-proton mass ratio                             ',      1.999007500970D+00,    180.000000000000D-12,'                    '),&
       constant('deuteron rms charge radius                             ',      2.142400000000D-15,      2.100000000000D-18,'    m               '),&
       constant('electric constant                                      ',      8.854187817000D-12,      0.000000000000D+00,'    F m^-1          '),&
       constant('electron charge to mass quotient                       ',   -175.882008800000D+09,      3.900000000000D+03,'    C kg^-1         '),&
       constant('electron-deuteron mag. mom. ratio                      ',     -2.143923498000D+03,     18.000000000000D-06,'                    '),&
       constant('electron-deuteron mass ratio                           ',    272.443710950000D-06,    110.000000000000D-15,'                    '),&
       constant('electron g factor                                      ',     -2.002319304362D+00,    530.000000000000D-15,'                    '),&
       constant('electron gyromag. ratio                                ',    176.085970800000D+09,      3.900000000000D+03,'    s^-1 T^-1       '),&
       constant('electron gyromag. ratio over 2 pi                      ',     28.024952660000D+03,    620.000000000000D-06,'    MHz T^-1        '),&
       constant('electron-helion mass ratio                             ',    181.954307610000D-06,    170.000000000000D-15,'                    '),&
       constant('electron mag. mom.                                     ',     -9.284764300000D-24,    210.000000000000D-33,'    J T^-1          '),&
       constant('electron mag. mom. anomaly                             ',      1.159652180760D-03,    270.000000000000D-15,'                    '),&
       constant('electron mag. mom. to Bohr magneton ratio              ',     -1.001159652181D+00,    270.000000000000D-15,'                    '),&
       constant('electron mag. mom. to nuclear magneton ratio           ',     -1.838281970900D+03,    750.000000000000D-09,'                    '),&
       constant('electron mass                                          ',    910.938291000000D-33,     40.000000000000D-39,'    kg              '),&
       constant('electron mass energy equivalent                        ',     81.871050600000D-15,      3.600000000000D-21,'    J               '),&
       constant('electron mass energy equivalent in MeV                 ',    510.998928000000D-03,     11.000000000000D-09,'    MeV             '),&
       constant('electron mass in u                                     ',    548.579909460000D-06,    220.000000000000D-15,'    u               '),&
       constant('electron molar mass                                    ',    548.579909460000D-09,    220.000000000000D-18,'    kg mol^-1       '),&
       constant('electron-muon mag. mom. ratio                          ',    206.766989600000D+00,      5.200000000000D-06,'                    '),&
       constant('electron-muon mass ratio                               ',      4.836331660000D-03,    120.000000000000D-12,'                    '),&
       constant('electron-neutron mag. mom. ratio                       ',    960.920500000000D+00,    230.000000000000D-06,'                    '),&
       constant('electron-neutron mass ratio                            ',    543.867344610000D-06,    320.000000000000D-15,'                    '),&
       constant('electron-proton mag. mom. ratio                        ',   -658.210684800000D+00,      5.400000000000D-06,'                    '),&
       constant('electron-proton mass ratio                             ',    544.617021780000D-06,    220.000000000000D-15,'                    '),&
       constant('electron-tau mass ratio                                ',    287.592000000000D-06,     26.000000000000D-09,'                    '),&
       constant('electron to alpha particle mass ratio                  ',    137.093355578000D-06,     55.000000000000D-15,'                    '),&
       constant('electron to shielded helion mag. mom. ratio            ',    864.058257000000D+00,     10.000000000000D-06,'                    '),&
       constant('electron to shielded proton mag. mom. ratio            ',   -658.227597100000D+00,      7.200000000000D-06,'                    '),&
       constant('electron-triton mass ratio                             ',    181.920006530000D-06,    170.000000000000D-15,'                    '),&
       constant('electron volt                                          ',    160.217656500000D-21,      3.500000000000D-27,'    J               '),&
       constant('electron volt-atomic mass unit relationship            ',      1.073544150000D-09,     24.000000000000D-18,'    u               '),&
       constant('electron volt-hartree relationship                     ',     36.749323790000D-03,    810.000000000000D-12,'    E_h             '),&
       constant('electron volt-hertz relationship                       ',    241.798934800000D+12,      5.300000000000D+06,'    Hz              '),&
       constant('electron volt-inverse meter relationship               ',    806.554429000000D+03,     18.000000000000D-03,'    m^-1            '),&
       constant('electron volt-joule relationship                       ',    160.217656500000D-21,      3.500000000000D-27,'    J               '),&
       constant('electron volt-kelvin relationship                      ',     11.604519000000D+03,     11.000000000000D-03,'    K               '),&
       constant('electron volt-kilogram relationship                    ',      1.782661845000D-36,     39.000000000000D-45,'    kg              '),&
       constant('elementary charge                                      ',    160.217656500000D-21,      3.500000000000D-27,'    C               '),&
       constant('elementary charge over h                               ',    241.798934800000D+12,      5.300000000000D+06,'    A J^-1          '),&
       constant('Faraday constant                                       ',     96.485336500000D+03,      2.100000000000D-03,'    C mol^-1        '),&
       constant('Faraday constant for conventional electric current     ',     96.485332100000D+03,      4.300000000000D-03,'    C_90 mol^-1     '),&
       constant('Fermi coupling constant                                ',     11.663640000000D-06,     50.000000000000D-12,'    GeV^-2          '),&
       constant('fine-structure constant                                ',      7.297352569800D-03,      2.400000000000D-12,'                    '),&
       constant('first radiation constant                               ',    374.177153000000D-18,     17.000000000000D-24,'    W m^2           '),&
       constant('first radiation constant for spectral radiance         ',    119.104286900000D-18,      5.300000000000D-24,'    W m^2 sr^-1     '),&
       constant('hartree-atomic mass unit relationship                  ',     29.212623246000D-09,     21.000000000000D-18,'    u               '),&
       constant('hartree-electron volt relationship                     ',     27.211385050000D+00,    600.000000000000D-09,'    eV              '),&
       constant('Hartree energy                                         ',      4.359744340000D-18,    190.000000000000D-27,'    J               '),&
       constant('Hartree energy in eV                                   ',     27.211385050000D+00,    600.000000000000D-09,'    eV              '),&
       constant('hartree-hertz relationship                             ',      6.579683920729D+15,     33.000000000000D+03,'    Hz              '),&
       constant('hartree-inverse meter relationship                     ',     21.947463137080D+06,    110.000000000000D-06,'    m^-1            '),&
       constant('hartree-joule relationship                             ',      4.359744340000D-18,    190.000000000000D-27,'    J               '),&
       constant('hartree-kelvin relationship                            ',    315.775040000000D+03,    290.000000000000D-03,'    K               '),&
       constant('hartree-kilogram relationship                          ',     48.508697900000D-36,      2.100000000000D-42,'    kg              '),&
       constant('helion-electron mass ratio                             ',      5.495885275400D+03,      5.000000000000D-06,'                    '),&
       constant('helion g factor                                        ',     -4.255250613000D+00,     50.000000000000D-09,'                    '),&
       constant('helion mag. mom.                                       ',    -10.746174860000D-27,    270.000000000000D-36,'    J T^-1          '),&
       constant('helion mag. mom. to Bohr magneton ratio                ',     -1.158740958000D-03,     14.000000000000D-12,'                    '),&
       constant('helion mag. mom. to nuclear magneton ratio             ',     -2.127625306000D+00,     25.000000000000D-09,'                    '),&
       constant('helion mass                                            ',      5.006412340000D-27,    220.000000000000D-36,'    kg              '),&
       constant('helion mass energy equivalent                          ',    449.953902000000D-12,     20.000000000000D-18,'    J               '),&
       constant('helion mass energy equivalent in MeV                   ',      2.808391482000D+03,     62.000000000000D-06,'    MeV             '),&
       constant('helion mass in u                                       ',      3.014932246800D+00,      2.500000000000D-09,'    u               '),&
       constant('helion molar mass                                      ',      3.014932246800D-03,      2.500000000000D-12,'    kg mol^-1       '),&
       constant('helion-proton mass ratio                               ',      2.993152670700D+00,      2.500000000000D-09,'                    '),&
       constant('hertz-atomic mass unit relationship                    ',      4.439821668900D-24,      3.100000000000D-33,'    u               '),&
       constant('hertz-electron volt relationship                       ',      4.135667516000D-15,     91.000000000000D-24,'    eV              '),&
       constant('hertz-hartree relationship                             ',    151.982984600450D-18,    760.000000000000D-30,'    E_h             '),&
       constant('hertz-inverse meter relationship                       ',      3.335640951000D-09,      0.000000000000D+00,'    m^-1            '),&
       constant('hertz-joule relationship                               ',    662.606957000000D-36,     29.000000000000D-42,'    J               '),&
       constant('hertz-kelvin relationship                              ',     47.992434000000D-12,     44.000000000000D-18,'    K               '),&
       constant('hertz-kilogram relationship                            ',      7.372496680000D-51,    330.000000000000D-60,'    kg              '),&
       constant('inverse fine-structure constant                        ',    137.035999074000D+00,     44.000000000000D-09,'                    '),&
       constant('inverse meter-atomic mass unit relationship            ',      1.331025051200D-15,    940.000000000000D-27,'    u               '),&
       constant('inverse meter-electron volt relationship               ',      1.239841930000D-06,     27.000000000000D-15,'    eV              '),&
       constant('inverse meter-hartree relationship                     ',     45.563352527550D-09,    230.000000000000D-21,'    E_h             '),&
       constant('inverse meter-hertz relationship                       ',    299.792458000000D+06,      0.000000000000D+00,'    Hz              '),&
       constant('inverse meter-joule relationship                       ',    198.644568400000D-27,      8.800000000000D-33,'    J               '),&
       constant('inverse meter-kelvin relationship                      ',     14.387770000000D-03,     13.000000000000D-09,'    K               '),&
       constant('inverse meter-kilogram relationship                    ',      2.210218902000D-42,     98.000000000000D-51,'    kg              '),&
       constant('inverse of conductance quantum                         ',     12.906403721700D+03,      4.200000000000D-06,'    ohm             '),&
       constant('Josephson constant                                     ',    483.597870000000D+12,     11.000000000000D+06,'    Hz V^-1         '),&
       constant('joule-atomic mass unit relationship                    ',      6.700535850000D+09,    300.000000000000D+00,'    u               '),&
       constant('joule-electron volt relationship                       ',      6.241509340000D+18,    140.000000000000D+09,'    eV              '),&
       constant('joule-hartree relationship                             ',    229.371248000000D+15,     10.000000000000D+09,'    E_h             '),&
       constant('joule-hertz relationship                               ',      1.509190311000D+33,     67.000000000000D+24,'    Hz              '),&
       constant('joule-inverse meter relationship                       ',      5.034117010000D+24,    220.000000000000D+15,'    m^-1            '),&
       constant('joule-kelvin relationship                              ',     72.429716000000D+21,     66.000000000000D+15,'    K               '),&
       constant('joule-kilogram relationship                            ',     11.126500560000D-18,      0.000000000000D+00,'    kg              '),&
       constant('kelvin-atomic mass unit relationship                   ',     92.510868000000D-15,     84.000000000000D-21,'    u               '),&
       constant('kelvin-electron volt relationship                      ',     86.173324000000D-06,     78.000000000000D-12,'    eV              '),&
       constant('kelvin-hartree relationship                            ',      3.166811400000D-06,      2.900000000000D-12,'    E_h             '),&
       constant('kelvin-hertz relationship                              ',     20.836618000000D+09,     19.000000000000D+03,'    Hz              '),&
       constant('kelvin-inverse meter relationship                      ',     69.503476000000D+00,     63.000000000000D-06,'    m^-1            '),&
       constant('kelvin-joule relationship                              ',     13.806488000000D-24,     13.000000000000D-30,'    J               '),&
       constant('kelvin-kilogram relationship                           ',    153.617900000000D-42,    140.000000000000D-48,'    kg              '),&
       constant('kilogram-atomic mass unit relationship                 ',    602.214129000000D+24,     27.000000000000D+18,'    u               '),&
       constant('kilogram-electron volt relationship                    ',    560.958885000000D+33,     12.000000000000D+27,'    eV              '),&
       constant('kilogram-hartree relationship                          ',     20.614859680000D+33,    910.000000000000D+24,'    E_h             '),&
       constant('kilogram-hertz relationship                            ',    135.639260800000D+48,      6.000000000000D+42,'    Hz              '),&
       constant('kilogram-inverse meter relationship                    ',    452.443873000000D+39,     20.000000000000D+33,'    m^-1            '),&
       constant('kilogram-joule relationship                            ',     89.875517870000D+15,      0.000000000000D+00,'    J               '),&
       constant('kilogram-kelvin relationship                           ',      6.509658200000D+39,      5.900000000000D+33,'    K               '),&
       constant('lattice parameter of silicon                           ',    543.102050400000D-12,      8.900000000000D-18,'    m               '),&
       constant('Loschmidt constant (273.15 K, 100 kPa)                 ',     26.516462000000D+24,     24.000000000000D+18,'    m^-3            '),&
       constant('Loschmidt constant (273.15 K, 101.325 kPa)             ',     26.867805000000D+24,     24.000000000000D+18,'    m^-3            '),&
       constant('mag. constant                                          ',      1.256637061400D-06,      0.000000000000D+00,'    N A^-2          '),&
       constant('mag. flux quantum                                      ',      2.067833758000D-15,     46.000000000000D-24,'    Wb              '),&
       constant('molar gas constant                                     ',      8.314462100000D+00,      7.500000000000D-06,'    J mol^-1 K^-1   '),&
       constant('molar mass constant                                    ',      1.000000000000D-03,      0.000000000000D+00,'    kg mol^-1       '),&
       constant('molar mass of carbon-12                                ',     12.000000000000D-03,      0.000000000000D+00,'    kg mol^-1       '),&
       constant('molar Planck constant                                  ',    399.031271760000D-12,    280.000000000000D-21,'    J s mol^-1      '),&
       constant('molar Planck constant times c                          ',    119.626565779000D-03,     84.000000000000D-12,'    J m mol^-1      '),&
       constant('molar volume of ideal gas (273.15 K, 100 kPa)          ',     22.710953000000D-03,     21.000000000000D-09,'    m^3 mol^-1      '),&
       constant('molar volume of ideal gas (273.15 K, 101.325 kPa)      ',     22.413968000000D-03,     20.000000000000D-09,'    m^3 mol^-1      '),&
       constant('molar volume of silicon                                ',     12.058833010000D-06,    800.000000000000D-15,'    m^3 mol^-1      '),&
       constant('Mo x unit                                              ',    100.209952000000D-15,     53.000000000000D-21,'    m               '),&
       constant('muon Compton wavelength                                ',     11.734441030000D-15,    300.000000000000D-24,'    m               '),&
       constant('muon Compton wavelength over 2 pi                      ',      1.867594294000D-15,     47.000000000000D-24,'    m               '),&
       constant('muon-electron mass ratio                               ',    206.768284300000D+00,      5.200000000000D-06,'                    '),&
       constant('muon g factor                                          ',     -2.002331841800D+00,      1.300000000000D-09,'                    '),&
       constant('muon mag. mom.                                         ',    -44.904480700000D-27,      1.500000000000D-33,'    J T^-1          '),&
       constant('muon mag. mom. anomaly                                 ',      1.165920910000D-03,    630.000000000000D-12,'                    '),&
       constant('muon mag. mom. to Bohr magneton ratio                  ',     -4.841970440000D-03,    120.000000000000D-12,'                    '),&
       constant('muon mag. mom. to nuclear magneton ratio               ',     -8.890596970000D+00,    220.000000000000D-09,'                    '),&
       constant('muon mass                                              ',    188.353147500000D-30,      9.600000000000D-36,'    kg              '),&
       constant('muon mass energy equivalent                            ',     16.928336670000D-12,    860.000000000000D-21,'    J               '),&
       constant('muon mass energy equivalent in MeV                     ',    105.658371500000D+00,      3.500000000000D-06,'    MeV             '),&
       constant('muon mass in u                                         ',    113.428926700000D-03,      2.900000000000D-09,'    u               '),&
       constant('muon molar mass                                        ',    113.428926700000D-06,      2.900000000000D-12,'    kg mol^-1       '),&
       constant('muon-neutron mass ratio                                ',    112.454517700000D-03,      2.800000000000D-09,'                    '),&
       constant('muon-proton mag. mom. ratio                            ',     -3.183345107000D+00,     84.000000000000D-09,'                    '),&
       constant('muon-proton mass ratio                                 ',    112.609527200000D-03,      2.800000000000D-09,'                    '),&
       constant('muon-tau mass ratio                                    ',     59.464900000000D-03,      5.400000000000D-06,'                    '),&
       constant('natural unit of action                                 ',    105.457172600000D-36,      4.700000000000D-42,'    J s             '),&
       constant('natural unit of action in eV s                         ',    658.211928000000D-18,     15.000000000000D-24,'    eV s            '),&
       constant('natural unit of energy                                 ',     81.871050600000D-15,      3.600000000000D-21,'    J               '),&
       constant('natural unit of energy in MeV                          ',    510.998928000000D-03,     11.000000000000D-09,'    MeV             '),&
       constant('natural unit of length                                 ',    386.159268000000D-15,    250.000000000000D-24,'    m               '),&
       constant('natural unit of mass                                   ',    910.938291000000D-33,     40.000000000000D-39,'    kg              '),&
       constant('natural unit of mom.um                                 ',    273.092429000000D-24,     12.000000000000D-30,'    kg m s^-1       '),&
       constant('natural unit of mom.um in MeV/c                        ',    510.998928000000D-03,     11.000000000000D-09,'    MeV/c           '),&
       constant('natural unit of time                                   ',      1.288088668330D-21,    830.000000000000D-33,'    s               '),&
       constant('natural unit of velocity                               ',    299.792458000000D+06,      0.000000000000D+00,'    m s^-1          '),&
       constant('neutron Compton wavelength                             ',      1.319590906800D-15,      1.100000000000D-24,'    m               '),&
       constant('neutron Compton wavelength over 2 pi                   ',    210.019415680000D-18,    170.000000000000D-27,'    m               '),&
       constant('neutron-electron mag. mom. ratio                       ',      1.040668820000D-03,    250.000000000000D-12,'                    '),&
       constant('neutron-electron mass ratio                            ',      1.838683660500D+03,      1.100000000000D-06,'                    '),&
       constant('neutron g factor                                       ',     -3.826085450000D+00,    900.000000000000D-09,'                    '),&
       constant('neutron gyromag. ratio                                 ',    183.247179000000D+06,     43.000000000000D+00,'    s^-1 T^-1       '),&
       constant('neutron gyromag. ratio over 2 pi                       ',     29.164694300000D+00,      6.900000000000D-06,'    MHz T^-1        '),&
       constant('neutron mag. mom.                                      ',     -9.662364700000D-27,      2.300000000000D-33,'    J T^-1          '),&
       constant('neutron mag. mom. to Bohr magneton ratio               ',     -1.041875630000D-03,    250.000000000000D-12,'                    '),&
       constant('neutron mag. mom. to nuclear magneton ratio            ',     -1.913042720000D+00,    450.000000000000D-09,'                    '),&
       constant('neutron mass                                           ',      1.674927351000D-27,     74.000000000000D-36,'    kg              '),&
       constant('neutron mass energy equivalent                         ',    150.534963100000D-12,      6.600000000000D-18,'    J               '),&
       constant('neutron mass energy equivalent in MeV                  ',    939.565379000000D+00,     21.000000000000D-06,'    MeV             '),&
       constant('neutron mass in u                                      ',      1.008664916000D+00,    430.000000000000D-12,'    u               '),&
       constant('neutron molar mass                                     ',      1.008664916000D-03,    430.000000000000D-15,'    kg mol^-1       '),&
       constant('neutron-muon mass ratio                                ',      8.892484000000D+00,    220.000000000000D-09,'                    '),&
       constant('neutron-proton mag. mom. ratio                         ',   -684.979340000000D-03,    160.000000000000D-09,'                    '),&
       constant('neutron-proton mass difference                         ',      2.305573920000D-30,    760.000000000000D-39,'                    '),&
       constant('neutron-proton mass difference energy equivalent       ',    207.214650000000D-15,     68.000000000000D-21,'                    '),&
       constant('neutron-proton mass difference energy equivalent in MeV',      1.293332170000D+00,    420.000000000000D-09,'                    '),&
       constant('neutron-proton mass difference in u                    ',      1.388449190000D-03,    450.000000000000D-12,'                    '),&
       constant('neutron-proton mass ratio                              ',      1.001378419170D+00,    450.000000000000D-12,'                    '),&
       constant('neutron-tau mass ratio                                 ',    528.790000000000D-03,     48.000000000000D-06,'                    '),&
       constant('neutron to shielded proton mag. mom. ratio             ',   -684.996940000000D-03,    160.000000000000D-09,'                    '),&
       constant('Newtonian constant of gravitation                      ',     66.738400000000D-12,      8.000000000000D-15,'    m^3 kg^-1 s^-2  '),&
       constant('Newtonian constant of gravitation over h-bar c         ',      6.708370000000D-39,    800.000000000000D-45,'    (GeV/c^2)^-2    '),&
       constant('nuclear magneton                                       ',      5.050783530000D-27,    110.000000000000D-36,'    J T^-1          '),&
       constant('nuclear magneton in eV/T                               ',     31.524512605000D-09,     22.000000000000D-18,'    eV T^-1         '),&
       constant('nuclear magneton in inverse meters per tesla           ',     25.426235270000D-03,    560.000000000000D-12,'    m^-1 T^-1       '),&
       constant('nuclear magneton in K/T                                ',    365.826820000000D-06,    330.000000000000D-12,'    K T^-1          '),&
       constant('nuclear magneton in MHz/T                              ',      7.622593570000D+00,    170.000000000000D-09,'    MHz T^-1        '),&
       constant('Planck constant                                        ',    662.606957000000D-36,     29.000000000000D-42,'    J s             '),&
       constant('Planck constant in eV s                                ',      4.135667516000D-15,     91.000000000000D-24,'    eV s            '),&
       constant('Planck constant over 2 pi                              ',    105.457172600000D-36,      4.700000000000D-42,'    J s             '),&
       constant('Planck constant over 2 pi in eV s                      ',    658.211928000000D-18,     15.000000000000D-24,'    eV s            '),&
       constant('Planck constant over 2 pi times c in MeV fm            ',    197.326971800000D+00,      4.400000000000D-06,'    MeV fm          '),&
       constant('Planck length                                          ',     16.161990000000D-36,    970.000000000000D-42,'    m               '),&
       constant('Planck mass                                            ',     21.765100000000D-09,      1.300000000000D-12,'    kg              '),&
       constant('Planck mass energy equivalent in GeV                   ',     12.209320000000D+18,    730.000000000000D+12,'    GeV             '),&
       constant('Planck temperature                                     ',    141.683300000000D+30,      8.500000000000D+27,'    K               '),&
       constant('Planck time                                            ',     53.910600000000D-45,      3.200000000000D-48,'    s               '),&
       constant('proton charge to mass quotient                         ',     95.788335800000D+06,      2.100000000000D+00,'    C kg^-1         '),&
       constant('proton Compton wavelength                              ',      1.321409856230D-15,    940.000000000000D-27,'    m               '),&
       constant('proton Compton wavelength over 2 pi                    ',    210.308910470000D-18,    150.000000000000D-27,'    m               '),&
       constant('proton-electron mass ratio                             ',      1.836152672450D+03,    750.000000000000D-09,'                    '),&
       constant('proton g factor                                        ',      5.585694713000D+00,     46.000000000000D-09,'                    '),&
       constant('proton gyromag. ratio                                  ',    267.522200500000D+06,      6.300000000000D+00,'    s^-1 T^-1       '),&
       constant('proton gyromag. ratio over 2 pi                        ',     42.577480600000D+00,      1.000000000000D-06,'    MHz T^-1        '),&
       constant('proton mag. mom.                                       ',     14.106067430000D-27,    330.000000000000D-36,'    J T^-1          '),&
       constant('proton mag. mom. to Bohr magneton ratio                ',      1.521032210000D-03,     12.000000000000D-12,'                    '),&
       constant('proton mag. mom. to nuclear magneton ratio             ',      2.792847356000D+00,     23.000000000000D-09,'                    '),&
       constant('proton mag. shielding correction                       ',     25.694000000000D-06,     14.000000000000D-09,'                    '),&
       constant('proton mass                                            ',      1.672621777000D-27,     74.000000000000D-36,'    kg              '),&
       constant('proton mass energy equivalent                          ',    150.327748400000D-12,      6.600000000000D-18,'    J               '),&
       constant('proton mass energy equivalent in MeV                   ',    938.272046000000D+00,     21.000000000000D-06,'    MeV             '),&
       constant('proton mass in u                                       ',      1.007276466812D+00,     90.000000000000D-12,'    u               '),&
       constant('proton molar mass                                      ',      1.007276466812D-03,     90.000000000000D-15,'    kg mol^-1       '),&
       constant('proton-muon mass ratio                                 ',      8.880243310000D+00,    220.000000000000D-09,'                    '),&
       constant('proton-neutron mag. mom. ratio                         ',     -1.459898060000D+00,    340.000000000000D-09,'                    '),&
       constant('proton-neutron mass ratio                              ',    998.623478260000D-03,    450.000000000000D-12,'                    '),&
       constant('proton rms charge radius                               ',    877.500000000000D-18,      5.100000000000D-18,'    m               '),&
       constant('proton-tau mass ratio                                  ',    528.063000000000D-03,     48.000000000000D-06,'                    '),&
       constant('quantum of circulation                                 ',    363.694755200000D-06,    240.000000000000D-15,'    m^2 s^-1        '),&
       constant('quantum of circulation times 2                         ',    727.389510400000D-06,    470.000000000000D-15,'    m^2 s^-1        '),&
       constant('Rydberg constant                                       ',     10.973731568539D+06,     55.000000000000D-06,'    m^-1            '),&
       constant('Rydberg constant times c in Hz                         ',      3.289841960364D+15,     17.000000000000D+03,'    Hz              '),&
       constant('Rydberg constant times hc in eV                        ',     13.605692530000D+00,    300.000000000000D-09,'    eV              '),&
       constant('Rydberg constant times hc in J                         ',      2.179872171000D-18,     96.000000000000D-27,'    J               '),&
       constant('Sackur-Tetrode constant (1 K, 100 kPa)                 ',     -1.151707800000D+00,      2.300000000000D-06,'                    '),&
       constant('Sackur-Tetrode constant (1 K, 101.325 kPa)             ',     -1.164870800000D+00,      2.300000000000D-06,'                    '),&
       constant('second radiation constant                              ',     14.387770000000D-03,     13.000000000000D-09,'    m K             '),&
       constant('shielded helion gyromag. ratio                         ',    203.789465900000D+06,      5.100000000000D+00,'    s^-1 T^-1       '),&
       constant('shielded helion gyromag. ratio over 2 pi               ',     32.434100840000D+00,    810.000000000000D-09,'    MHz T^-1        '),&
       constant('shielded helion mag. mom.                              ',    -10.745530440000D-27,    270.000000000000D-36,'    J T^-1          '),&
       constant('shielded helion mag. mom. to Bohr magneton ratio       ',     -1.158671471000D-03,     14.000000000000D-12,'                    '),&
       constant('shielded helion mag. mom. to nuclear magneton ratio    ',     -2.127497718000D+00,     25.000000000000D-09,'                    '),&
       constant('shielded helion to proton mag. mom. ratio              ',   -761.766558000000D-03,     11.000000000000D-09,'                    '),&
       constant('shielded helion to shielded proton mag. mom. ratio     ',   -761.786131300000D-03,      3.300000000000D-09,'                    '),&
       constant('shielded proton gyromag. ratio                         ',    267.515326800000D+06,      6.600000000000D+00,'    s^-1 T^-1       '),&
       constant('shielded proton gyromag. ratio over 2 pi               ',     42.576386600000D+00,      1.000000000000D-06,'    MHz T^-1        '),&
       constant('shielded proton mag. mom.                              ',     14.105704990000D-27,    350.000000000000D-36,'    J T^-1          '),&
       constant('shielded proton mag. mom. to Bohr magneton ratio       ',      1.520993128000D-03,     17.000000000000D-12,'                    '),&
       constant('shielded proton mag. mom. to nuclear magneton ratio    ',      2.792775598000D+00,     30.000000000000D-09,'                    '),&
       constant('speed of light in vacuum                               ',    299.792458000000D+06,      0.000000000000D+00,'    m s^-1          '),&
       constant('standard acceleration of gravity                       ',      9.806650000000D+00,      0.000000000000D+00,'    m s^-2          '),&
       constant('standard atmosphere                                    ',    101.325000000000D+03,      0.000000000000D+00,'    Pa              '),&
       constant('standard-state pressure                                ',    100.000000000000D+03,      0.000000000000D+00,'    Pa              '),&
       constant('Stefan-Boltzmann constant                              ',     56.703730000000D-09,    210.000000000000D-15,'    W m^-2 K^-4     '),&
       constant('tau Compton wavelength                                 ',    697.787000000000D-18,     63.000000000000D-21,'    m               '),&
       constant('tau Compton wavelength over 2 pi                       ',    111.056000000000D-18,     10.000000000000D-21,'    m               '),&
       constant('tau-electron mass ratio                                ',      3.477150000000D+03,    310.000000000000D-03,'                    '),&
       constant('tau mass                                               ',      3.167470000000D-27,    290.000000000000D-33,'    kg              '),&
       constant('tau mass energy equivalent                             ',    284.678000000000D-12,     26.000000000000D-15,'    J               '),&
       constant('tau mass energy equivalent in MeV                      ',      1.776820000000D+03,    160.000000000000D-03,'    MeV             '),&
       constant('tau mass in u                                          ',      1.907490000000D+00,    170.000000000000D-06,'    u               '),&
       constant('tau molar mass                                         ',      1.907490000000D-03,    170.000000000000D-09,'    kg mol^-1       '),&
       constant('tau-muon mass ratio                                    ',     16.816700000000D+00,      1.500000000000D-03,'                    '),&
       constant('tau-neutron mass ratio                                 ',      1.891110000000D+00,    170.000000000000D-06,'                    '),&
       constant('tau-proton mass ratio                                  ',      1.893720000000D+00,    170.000000000000D-06,'                    '),&
       constant('Thomson cross section                                  ',     66.524587340000D-30,    130.000000000000D-39,'    m^2             '),&
       constant('triton-electron mass ratio                             ',      5.496921526700D+03,      5.000000000000D-06,'                    '),&
       constant('triton g factor                                        ',      5.957924896000D+00,     76.000000000000D-09,'                    '),&
       constant('triton mag. mom.                                       ',     15.046094470000D-27,    380.000000000000D-36,'    J T^-1          '),&
       constant('triton mag. mom. to Bohr magneton ratio                ',      1.622393657000D-03,     21.000000000000D-12,'                    '),&
       constant('triton mag. mom. to nuclear magneton ratio             ',      2.978962448000D+00,     38.000000000000D-09,'                    '),&
       constant('triton mass                                            ',      5.007356300000D-27,    220.000000000000D-36,'    kg              '),&
       constant('triton mass energy equivalent                          ',    450.038741000000D-12,     20.000000000000D-18,'    J               '),&
       constant('triton mass energy equivalent in MeV                   ',      2.808921005000D+03,     62.000000000000D-06,'    MeV             '),&
       constant('triton mass in u                                       ',      3.015500713400D+00,      2.500000000000D-09,'    u               '),&
       constant('triton molar mass                                      ',      3.015500713400D-03,      2.500000000000D-12,'    kg mol^-1       '),&
       constant('triton-proton mass ratio                               ',      2.993717030800D+00,      2.500000000000D-09,'                    '),&
       constant('unified atomic mass unit                               ',      1.660538921000D-27,     73.000000000000D-36,'    kg              '),&
       constant('von Klitzing constant                                  ',     25.812807443400D+03,      8.400000000000D-06,'    ohm             '),&
       constant('weak mixing angle                                      ',    222.300000000000D-03,      2.100000000000D-03,'                    '),&
       constant('Wien frequency displacement law constant               ',     58.789254000000D+09,     53.000000000000D+03,'    Hz K^-1         ')&
       /)

  public :: value_constant
  public :: name_constant
  public :: units_constant
  public :: parse_constant
  public :: print_constant
contains

  function value_constant(key) result(value)
    integer        :: key    
    type(constant) :: const
    real(8)        :: value
    if(key>nmax_constant)then
       print*,"No constant associated to this number. nmax_constant=",nmax_constant
       value=0.d0
       return
    endif
    value=fundamental_constant(key)%val
  end function value_constant

  function name_constant(key) result(name)
    integer          :: key    
    type(constant)   :: const
    character(len=64) :: name
    if(key>nmax_constant)then
       print*,"No constant associated to this number. nmax_constant=",nmax_constant
       name=""
       return
    endif
    name=fundamental_constant(key)%name
  end function name_constant

  function units_constant(key) result(name)
    integer          :: key    
    type(constant)   :: const
    character(len=16) :: name
    if(key>nmax_constant)then
       print*,"No constant associated to this number. nmax_constant=",nmax_constant
       name=""
       return
    endif
    name=fundamental_constant(key)%units
  end function units_constant

  subroutine print_constant(key)
    integer        :: key    
    type(constant) :: const
    if(key>nmax_constant)then
       print*,"No constant associated to this number. nmax_constant=",nmax_constant
       return
    endif
    const=fundamental_constant(key)
    write(*,*)"Name     = ",trim(adjustl(trim(const%name)))
    write(*,*)"Value    =",const%val," err",const%uncert
    write(*,*)"Units    =",const%units
  end subroutine print_constant


  function parse_constant(buffer) result(const)
    type(constant)    :: const
    character(len=*)  :: buffer
    character(len=32) :: cvalue
    integer :: irep
    !GET NAME:
    const%name = adjustl(buffer(1:59))
    !GET VALUE:
    cvalue     = buffer(59:84)
    call s_s_delete2(cvalue,"...",irep)
    call s_blank_delete(cvalue)
    read(cvalue,*)const%val
    !GET UNCERTAINTY
    cvalue     = buffer(84:109)
    call s_blank_delete(cvalue)
    if(trim(cvalue)=="(exact)")then
       const%uncert=0.d0
    else
       read(cvalue,*)const%uncert
    endif
    !GET UNITS
    const%units = adjustl(buffer(109:))
  end function parse_constant


  subroutine get_constant_from_file(fileIN,Nmax,fileOUT)
    character(len=*)   :: fileIN
    integer            :: i,Nmax
    character(len=*)   :: fileOUT
    type(constant)     :: const
    character(len=160) ::string
    open(10,file=fileIN)
    open(11,file=fileOUT)
    do i=1,nmax
       read(10,"(A128)")string
       const =  parse_constant(string)       
       write(*,"(A20,A55,A2,EN24.12,A1,EN24.12,A2,A20,A4)")"constant('",const%name,"',",const%val,",",const%uncert,",'",const%units,"'),&"
       write(11,"(A20,A55,A2,EN24.12,A1,EN24.12,A2,A20,A4)")"constant('",const%name,"',",const%val,",",const%uncert,",'",const%units,"'),&"
    enddo
  end subroutine get_constant_from_file


  subroutine s_blank_delete ( s )
    !! S_BLANK_DELETE removes blanks from a string, left justifying the remainder.
    !    All TAB characters are also removed.
    !  Parameters:
    !    Input/output, character ( len = * ) S, the string to be transformed.
    character              ch
    integer   ( kind = 4 ) get
    integer   ( kind = 4 ) put
    character ( len = * )  s
    integer   ( kind = 4 ) s_length
    character, parameter :: tab = achar ( 9 )
    put = 0
    s_length = len_trim ( s )
    do get = 1, s_length
       ch = s(get:get)
       if ( ch /= ' ' .and. ch /= tab ) then
          put = put + 1
          s(put:put) = ch
       end if
    end do
    s(put+1:s_length) = ' '
    return
  end subroutine s_blank_delete

  subroutine s_s_delete2 ( s, sub, irep )
    !! S_S_DELETE2 recursively removes a substring from a string.
    !    The remainder is left justified and padded with blanks.
    !    The substitution is recursive, so
    !    that, for example, removing all occurrences of "ab" from
    !    "aaaaabbbbbQ" results in "Q".
    !  Parameters:
    !    Input/output, character ( len = * ) S, the string to be transformed.
    !    Input, character ( len = * ) SUB, the substring to be removed.
    !    Output, integer ( kind = 4 ) IREP, the number of occurrences of
    !    the substring.
    integer   ( kind = 4 ) ihi
    integer   ( kind = 4 ) irep
    integer   ( kind = 4 ) loc
    integer   ( kind = 4 ) nsub
    character ( len = * )  s
    integer   ( kind = 4 ) s_length
    character ( len = * )  sub
    s_length = len ( s )
    nsub = len ( sub )
    irep = 0
    ihi = s_length
    do while ( 0 < ihi )
       loc = index ( s(1:ihi), sub )
       if ( loc == 0 ) then
          return
       end if
       irep = irep + 1
       call s_chop ( s, loc, loc+nsub-1 )
       ihi = ihi - nsub
    end do
    return
  end subroutine s_s_delete2

  subroutine s_chop ( s, ilo, ihi )
    !! S_CHOP "chops out" a portion of a string, and closes up the hole.
    !  Example:
    !    S = 'Fred is not a jerk!'
    !    call s_chop ( S, 9, 12 )
    !    S = 'Fred is a jerk!    '
    !  Parameters:
    !    Input/output, character ( len = * ) S, the string to be transformed.
    !    Input, integer ( kind = 4 ) ILO, IHI, the locations of the first and last
    !    characters to be removed.
    integer   ( kind = 4 ) ihi
    integer   ( kind = 4 ) ihi2
    integer   ( kind = 4 ) ilo
    integer   ( kind = 4 ) ilo2
    character ( len = * )  s
    integer   ( kind = 4 ) s_length
    s_length = len ( s )
    ilo2 = max ( ilo, 1 )
    ihi2 = min ( ihi, s_length )
    if ( ihi2 < ilo2 ) then
       return
    end if
    s(ilo2:s_length+ilo2-ihi2-1) = s(ihi2+1:s_length)
    s(s_length+ilo2-ihi2:s_length) = ' '
  end subroutine s_chop

end module CONSTANTS
