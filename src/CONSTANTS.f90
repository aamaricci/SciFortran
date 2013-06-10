!include "contants_text.h"
module CONSTANTS
  implicit none
  private

  type,public:: constant
     character(len=128):: name
     real(8)           :: val
     real(8)           :: uncert
     character(len=21) :: units
  end type constant

  integer,parameter,public :: nmax_constant=334
  !constants updated to 2013
  !check this website:
  !http://physics.nist.gov/cuu/Constants/index.html
  type(constant),dimension(nmax_constant),public  :: fundamental_constant=[& 
       constant('lattice spacing of silicon',	  1.920155714000000D-010 ,	  3.200000000000000D-018 ,	 'm'),&
       constant('alpha particle-electron mass ratio',	   7294.29953610000      ,	  2.900000000000000D-006 ,	 ''),&
       constant('alpha particle mass',	  6.644656750000000D-027 ,	  2.900000000000000D-034 ,	 'kg'),&
       constant('alpha particle mass energy equivalent',	  5.971919670000000D-010 ,	  2.600000000000000D-017 ,	 'J'),&
       constant('alpha particle mass energy equivalent in MeV',	   3727.37924000000      ,	  8.200000000000000D-005 ,	 'MeV'),&
       constant('alpha particle mass in u',	   4.00150617912500      ,	  6.200000000000001D-011 ,	 'u'),&
       constant('alpha particle molar mass',	  4.001506179125000D-003 ,	  6.200000000000000D-014 ,	 'kg mol^-1'),&
       constant('alpha particle-proton mass ratio',	   3.97259968933000      ,	  3.600000000000000D-010 ,	 ''),&
       constant('Angstrom star',	  1.000014950000000D-010 ,	  9.000000000000000D-017 ,	 'm'),&
       constant('atomic mass constant',	  1.660538921000000D-027 ,	  7.300000000000000D-035 ,	 'kg'),&
       constant('atomic mass constant energy equivalent',	  1.492417954000000D-010 ,	  6.600000000000000D-018 ,	 'J'),&
       constant('atomic mass constant energy equivalent in MeV',	   931.494061000000      ,	  2.100000000000000D-005 ,	 'MeV'),&
       constant('atomic mass unit-electron volt relationship',	   931494061.000000      ,	   21.0000000000000      ,	 'eV'),&
       constant('atomic mass unit-hartree relationship',	   34231776.8450000      ,	  2.400000000000000D-002 ,	 'E_h'),&
       constant('atomic mass unit-hertz relationship',	  2.252342716800000D+023 ,	   160000000000000.      ,	 'Hz'),&
       constant('atomic mass unit-inverse meter relationship',	   751300660420000.      ,	   530000.000000000      ,	 'm^-1'),&
       constant('atomic mass unit-joule relationship',	  1.492417954000000D-010 ,	  6.600000000000000D-018 ,	 'J'),&
       constant('atomic mass unit-kelvin relationship',	   10809540800000.0      ,	   9800000.00000000      ,	 'K'),&
       constant('atomic mass unit-kilogram relationship',	  1.660538921000000D-027 ,	  7.300000000000000D-035 ,	 'kg'),&
       constant('atomic unit of 1st hyperpolarizability',	  3.206361449000000D-053 ,	  7.100000000000000D-061 ,	 'C^3 m^3 J^-2'),&
       constant('atomic unit of 2nd hyperpolarizability',	  6.235380540000000D-065 ,	  2.800000000000000D-072 ,	 'C^4 m^4 J^-3'),&
       constant('atomic unit of action',	  1.054571726000000D-034 ,	  4.700000000000000D-042 ,	 'J s'),&
       constant('atomic unit of charge',	  1.602176565000000D-019 ,	  3.500000000000000D-027 ,	 'C'),&
       constant('atomic unit of charge density',	   1081202338000.00      ,	   24000.0000000000      ,	 'C m^-3'),&
       constant('atomic unit of current',	  6.623617950000000D-003 ,	  1.500000000000000D-010 ,	 'A'),&
       constant('atomic unit of electric dipole mom.',	  8.478353260000001D-030 ,	  1.900000000000000D-037 ,	 'C m'),&
       constant('atomic unit of electric field',	   514220652000.000      ,	   11000.0000000000      ,	 'V m^-1'),&
       constant('atomic unit of electric field gradient',	  9.717362000000000D+021 ,	   210000000000000.      ,	 'V m^-2'),&
       constant('atomic unit of electric polarizability',	  1.648777275400000D-041 ,	  1.600000000000000D-050 ,	 'C^2 m^2 J^-1'),&
       constant('atomic unit of electric potential',	   27.2113850500000      ,	  6.000000000000000D-007 ,	 'V'),&
       constant('atomic unit of electric quadrupole mom.',	  4.486551331000000D-040 ,	  9.900000000000000D-048 ,	 'C m^2'),&
       constant('atomic unit of energy',	  4.359744340000000D-018 ,	  1.900000000000000D-025 ,	 'J'),&
       constant('atomic unit of force',	  8.238722780000000D-008 ,	  3.600000000000000D-015 ,	 'N'),&
       constant('atomic unit of length',	  5.291772109200000D-011 ,	  1.700000000000000D-020 ,	 'm'),&
       constant('atomic unit of mag. dipole mom.',	  1.854801936000000D-023 ,	  4.100000000000000D-031 ,	 'J T^-1'),&
       constant('atomic unit of mag. flux density',	   235051.746400000      ,	  5.200000000000000D-003 ,	 'T'),&
       constant('atomic unit of magnetizability',	  7.891036607000000D-029 ,	  1.300000000000000D-037 ,	 'J T^-2'),&
       constant('atomic unit of mass',	  9.109382910000000D-031 ,	  4.000000000000000D-038 ,	 'kg'),&
       constant('atomic unit of mom.um',	  1.992851740000000D-024 ,	  8.800000000000000D-032 ,	 'kg m s^-1'),&
       constant('atomic unit of permittivity',	  1.112650056000000D-010 ,	  0.000000000000000D+000 ,	 'F m^-1'),&
       constant('atomic unit of time',	  2.418884326502000D-017 ,	  1.200000000000000D-028 ,	 's'),&
       constant('atomic unit of velocity',	   2187691.26379000      ,	  7.100000000000000D-004 ,	 'm s^-1'),&
       constant('Avogadro constant',	  6.022141290000000D+023 ,	  2.700000000000000D+016 ,	 'mol^-1'),&
       constant('Bohr magneton',	  9.274009680000000D-024 ,	  2.000000000000000D-031 ,	 'J T^-1'),&
       constant('Bohr magneton in eV/T',	  5.788381806600000D-005 ,	  3.800000000000000D-014 ,	 'eV T^-1'),&
       constant('Bohr magneton in Hz/T',	   13996245550.0000      ,	   310.000000000000      ,	 'Hz T^-1'),&
       constant('Bohr magneton in inverse meters per tesla',	   46.6864498000000      ,	  1.000000000000000D-006 ,	 'm^-1 T^-1'),&
       constant('Bohr magneton in K/T',	  0.671713880000000      ,	  6.100000000000000D-007 ,	 'K T^-1'),&
       constant('Bohr radius',	  5.291772109200000D-011 ,	  1.700000000000000D-020 ,	 'm'),&
       constant('Boltzmann constant',	  1.380648800000000D-023 ,	  1.300000000000000D-029 ,	 'J K^-1'),&
       constant('Boltzmann constant in eV/K',	  8.617332400000000D-005 ,	  7.800000000000000D-011 ,	 'eV K^-1'),&
       constant('Boltzmann constant in Hz/K',	   20836618000.0000      ,	   19000.0000000000      ,	 'Hz K^-1'),&
       constant('Boltzmann constant in inverse meters per kelvin',69.5034760000000,6.300000000000000D-005 ,'m^-1 K^-1'),&
       constant('characteristic impedance of vacuum',	   376.730313461000      ,	  0.000000000000000D+000 ,	 'ohm'),&
       constant('classical electron radius',	  2.817940326700000D-015 ,	  2.700000000000000D-024 ,	 'm'),&
       constant('Compton wavelength',	  2.426310238900000D-012 ,	  1.600000000000000D-021 ,	 'm'),&
       constant('Compton wavelength over 2 pi',	  3.861592680000000D-013 ,	  2.500000000000000D-022 ,	 'm'),&
       constant('conductance quantum',	  7.748091734599999D-005 ,	  2.500000000000000D-014 ,	 'S'),&
       constant('conventional value of Josephson constant',	   483597900000000.      ,	  0.000000000000000D+000 ,	 'Hz V^-1'),&
       constant('conventional value of von Klitzing constant',	   25812.8070000000      ,	  0.000000000000000D+000 ,	 'ohm'),&
       constant('Cu x unit',	  1.002076970000000D-013 ,	  2.800000000000000D-020 ,	 'm'),&
       constant('deuteron-electron mag. mom. ratio',	 -4.664345537000000D-004 ,	  3.900000000000000D-012 ,	 ''),&
       constant('deuteron-electron mass ratio',	   3670.48296520000      ,	  1.500000000000000D-006 ,	 ''),&
       constant('deuteron g factor',	  0.857438230800000      ,	  7.200000000000000D-009 ,	 ''),&
       constant('deuteron mag. mom.',	  4.330734890000000D-027 ,	  9.999999999999999D-035 ,	 'J T^-1'),&
       constant('deuteron mag. mom. to Bohr magneton ratio',	  4.669754556000000D-004 ,	  3.900000000000000D-012 ,	 ''),&
       constant('deuteron mag. mom. to nuclear magneton ratio',	  0.857438230800000      ,	  7.200000000000000D-009 ,	 ''),&
       constant('deuteron mass',	  3.343583480000000D-027 ,	  1.500000000000000D-034 ,	 'kg'),&
       constant('deuteron mass energy equivalent',	  3.005062970000000D-010 ,	  1.300000000000000D-017 ,	 'J'),&
       constant('deuteron mass energy equivalent in MeV',	   1875.61285900000      ,	  4.100000000000000D-005 ,	 'MeV'),&
       constant('deuteron mass in u',	   2.01355321271200      ,	  7.700000000000001D-011 ,	 'u'),&
       constant('deuteron molar mass',	  2.013553212712000D-003 ,	  7.700000000000000D-014 ,	 'kg mol^-1'),&
       constant('deuteron-neutron mag. mom. ratio',	 -0.448206520000000      ,	  1.100000000000000D-007 ,	 ''),&
       constant('deuteron-proton mag. mom. ratio',	  0.307012207000000      ,	  2.400000000000000D-009 ,	 ''),&
       constant('deuteron-proton mass ratio',	   1.99900750097000      ,	  1.800000000000000D-010 ,	 ''),&
       constant('deuteron rms charge radius',	  2.142400000000000D-015 ,	  2.100000000000000D-018 ,	 'm'),&
       constant('electric constant',	  8.854187817000000D-012 ,	  0.000000000000000D+000 ,	 'F m^-1'),&
       constant('electron charge to mass quotient',	  -175882008800.000      ,	   3900.00000000000      ,	 'C kg^-1'),&
       constant('electron-deuteron mag. mom. ratio',	  -2143.92349800000      ,	  1.800000000000000D-005 ,	 ''),&
       constant('electron-deuteron mass ratio',	  2.724437109500000D-004 ,	  1.100000000000000D-013 ,	 ''),&
       constant('electron g factor',	  -2.00231930436200      ,	  5.300000000000000D-013 ,	 ''),&
       constant('electron gyromag. ratio',	   176085970800.000      ,	   3900.00000000000      ,	 's^-1 T^-1'),&
       constant('electron gyromag. ratio over 2 pi',	   28024.9526600000      ,	  6.200000000000000D-004 ,	 'MHz T^-1'),&
       constant('electron-helion mass ratio',	  1.819543076100000D-004 ,	  1.700000000000000D-013 ,	 ''),&
       constant('electron mag. mom.',	 -9.284764300000000D-024 ,	  2.100000000000000D-031 ,	 'J T^-1'),&
       constant('electron mag. mom. anomaly',	  1.159652180760000D-003 ,	  2.700000000000000D-013 ,	 ''),&
       constant('electron mag. mom. to Bohr magneton ratio',	  -1.00115965218100      ,	  2.700000000000000D-013 ,	 ''),&
       constant('electron mag. mom. to nuclear magneton ratio',	  -1838.28197090000      ,	  7.500000000000000D-007 ,	 ''),&
       constant('electron mass',	  9.109382910000000D-031 ,	  4.000000000000000D-038 ,	 'kg'),&
       constant('electron mass energy equivalent',	  8.187105060000000D-014 ,	  3.600000000000000D-021 ,	 'J'),&
       constant('electron mass energy equivalent in MeV',	  0.510998928000000      ,	  1.100000000000000D-008 ,	 'MeV'),&
       constant('electron mass in u',	  5.485799094600000D-004 ,	  2.200000000000000D-013 ,	 'u'),&
       constant('electron molar mass',	  5.485799094600000D-007 ,	  2.200000000000000D-016 ,	 'kg mol^-1'),&
       constant('electron-muon mag. mom. ratio',	   206.766989600000      ,	  5.200000000000000D-006 ,	 ''),&
       constant('electron-muon mass ratio',	  4.836331660000000D-003 ,	  1.200000000000000D-010 ,	 ''),&
       constant('electron-neutron mag. mom. ratio',	   960.920500000000      ,	  2.300000000000000D-004 ,	 ''),&
       constant('electron-neutron mass ratio',	  5.438673446100000D-004 ,	  3.200000000000000D-013 ,	 ''),&
       constant('electron-proton mag. mom. ratio',	  -658.210684800000      ,	  5.400000000000000D-006 ,	 ''),&
       constant('electron-proton mass ratio',	  5.446170217800000D-004 ,	  2.200000000000000D-013 ,	 ''),&
       constant('electron-tau mass ratio',	  2.875920000000000D-004 ,	  2.600000000000000D-008 ,	 ''),&
       constant('electron to alpha particle mass ratio',	  1.370933555780000D-004 ,	  5.500000000000000D-014 ,	 ''),&
       constant('electron to shielded helion mag. mom. ratio',	   864.058257000000      ,	  1.000000000000000D-005 ,	 ''),&
       constant('electron to shielded proton mag. mom. ratio',	  -658.227597100000      ,	  7.200000000000000D-006 ,	 ''),&
       constant('electron-triton mass ratio',	  1.819200065300000D-004 ,	  1.700000000000000D-013 ,	 ''),&
       constant('electron volt',	  1.602176565000000D-019 ,	  3.500000000000000D-027 ,	 'J'),&
       constant('electron volt-atomic mass unit relationship',	  1.073544150000000D-009 ,	  2.400000000000000D-017 ,	 'u'),&
       constant('electron volt-hartree relationship',	  3.674932379000000D-002 ,	  8.100000000000000D-010 ,	 'E_h'),&
       constant('electron volt-hertz relationship',	   241798934800000.      ,	   5300000.00000000      ,	 'Hz'),&
       constant('electron volt-inverse meter relationship',	   806554.429000000      ,	  1.800000000000000D-002 ,	 'm^-1'),&
       constant('electron volt-joule relationship',	  1.602176565000000D-019 ,	  3.500000000000000D-027 ,	 'J'),&
       constant('electron volt-kelvin relationship',	   11604.5190000000      ,	  1.100000000000000D-002 ,	 'K'),&
       constant('electron volt-kilogram relationship',	  1.782661845000000D-036 ,	  3.900000000000000D-044 ,	 'kg'),&
       constant('elementary charge',	  1.602176565000000D-019 ,	  3.500000000000000D-027 ,	 'C'),&
       constant('elementary charge over h',	   241798934800000.      ,	   5300000.00000000      ,	 'A J^-1'),&
       constant('Faraday constant',	   96485.3365000000      ,	  2.100000000000000D-003 ,	 'C mol^-1'),&
       constant('Faraday constant for conventional electric current',  96485.3321000000   ,4.300000000000000D-003 ,'C_90 mol^-1'),&
       constant('Fermi coupling constant',	  1.166364000000000D-005 ,	  5.000000000000000D-011 ,	 'GeV^-2'),&
       constant('fine-structure constant',	  7.297352569800000D-003 ,	  2.400000000000000D-012 ,	 ''),&
       constant('first radiation constant',	  3.741771530000000D-016 ,	  1.700000000000000D-023 ,	 'W m^2'),&
       constant('first radiation constant for spectral radiance', 1.191042869000000D-016 , 5.300000000000000D-024 ,'W m^2 sr^-1'),&
       constant('hartree-atomic mass unit relationship',	  2.921262324600000D-008 ,  2.100000000000000D-017 ,	 'u'),&
       constant('hartree-electron volt relationship',	   27.2113850500000      ,6.000000000000000D-007 ,	 'eV'),&
       constant('Hartree energy',	  4.359744340000000D-018 , 1.900000000000000D-025 ,	 'J'),&
       constant('Hartree energy in eV',	   27.2113850500000      , 6.000000000000000D-007 ,	 'eV'),&
       constant('hartree-hertz relationship',	  6.579683920729000D+015 ,  33000.0000000000      ,	 'Hz'),&
       constant('hartree-inverse meter relationship',	   21947463.1370800     ,1.100000000000000D-004 ,	 'm^-1'),&
       constant('hartree-joule relationship',	  4.359744340000000D-018 ,	  1.900000000000000D-025 ,	 'J'),&
       constant('hartree-kelvin relationship',	   315775.040000000      ,	  0.290000000000000      ,	 'K'),&
       constant('hartree-kilogram relationship',	  4.850869790000000D-035 ,	  2.100000000000000D-042 ,	 'kg'),&
       constant('helion-electron mass ratio',	   5495.88527540000      ,	  5.000000000000000D-006 ,	 ''),&
       constant('helion g factor',	  -4.25525061300000      ,	  5.000000000000000D-008 ,	 ''),&
       constant('helion mag. mom.',	 -1.074617486000000D-026 ,	  2.700000000000000D-034 ,	 'J T^-1'),&
       constant('helion mag. mom. to Bohr magneton ratio',	 -1.158740958000000D-003 ,	  1.400000000000000D-011 ,	 ''),&
       constant('helion mag. mom. to nuclear magneton ratio',	  -2.12762530600000      ,	  2.500000000000000D-008 ,	 ''),&
       constant('helion mass',	  5.006412340000000D-027 ,	  2.200000000000000D-034 ,	 'kg'),&
       constant('helion mass energy equivalent',	  4.499539020000000D-010 ,	  2.000000000000000D-017 ,	 'J'),&
       constant('helion mass energy equivalent in MeV',	   2808.39148200000      ,	  6.200000000000000D-005 ,	 'MeV'),&
       constant('helion mass in u',	   3.01493224680000      ,	  2.500000000000000D-009 ,	 'u'),&
       constant('helion molar mass',	  3.014932246800000D-003 ,	  2.500000000000000D-012 ,	 'kg mol^-1'),&
       constant('helion-proton mass ratio',	   2.99315267070000      ,	  2.500000000000000D-009 ,	 ''),&
       constant('hertz-atomic mass unit relationship',	  4.439821668900000D-024 ,	  3.100000000000000D-033 ,	 'u'),&
       constant('hertz-electron volt relationship',	  4.135667516000000D-015 ,	  9.100000000000000D-023 ,	 'eV'),&
       constant('hertz-hartree relationship',	  1.519829846004500D-016 ,	  7.600000000000000D-028 ,	 'E_h'),&
       constant('hertz-inverse meter relationship',	  3.335640951000000D-009 ,	  0.000000000000000D+000 ,	 'm^-1'),&
       constant('hertz-joule relationship',	  6.626069570000000D-034 ,	  2.900000000000000D-041 ,	 'J'),&
       constant('hertz-kelvin relationship',	  4.799243400000000D-011 ,	  4.400000000000000D-017 ,	 'K'),&
       constant('hertz-kilogram relationship',	  7.372496680000000D-051 ,	  3.300000000000000D-058 ,	 'kg'),&
       constant('inverse fine-structure constant',	   137.035999074000      ,	  4.400000000000000D-008 ,	 ''),&
       constant('inverse meter-atomic mass unit relationship',	  1.331025051200000D-015 ,	  9.399999999999999D-025 ,	 'u'),&
       constant('inverse meter-electron volt relationship',	  1.239841930000000D-006 ,	  2.700000000000000D-014 ,	 'eV'),&
       constant('inverse meter-hartree relationship',	  4.556335252755000D-008 ,	  2.300000000000000D-019 ,	 'E_h'),&
       constant('inverse meter-hertz relationship',	   299792458.000000      ,	  0.000000000000000D+000 ,	 'Hz'),&
       constant('inverse meter-joule relationship',	  1.986445684000000D-025 ,	  8.800000000000000D-033 ,	 'J'),&
       constant('inverse meter-kelvin relationship',	  1.438777000000000D-002 ,	  1.300000000000000D-008 ,	 'K'),&
       constant('inverse meter-kilogram relationship',	  2.210218902000000D-042 ,	  9.800000000000000D-050 ,	 'kg'),&
       constant('inverse of conductance quantum',	   12906.4037217000      ,	  4.200000000000000D-006 ,	 'ohm'),&
       constant('Josephson constant',	   483597870000000.      ,	   11000000.0000000      ,	 'Hz V^-1'),&
       constant('joule-atomic mass unit relationship',	   6700535850.00000      ,	   300.000000000000      ,	 'u'),&
       constant('joule-electron volt relationship',	  6.241509340000000D+018 ,	   140000000000.000      ,	 'eV'),&
       constant('joule-hartree relationship',	  2.293712480000000D+017 ,	   10000000000.0000      ,	 'E_h'),&
       constant('joule-hertz relationship',	  1.509190311000000D+033 ,	  6.700000000000000D+025 ,	 'Hz'),&
       constant('joule-inverse meter relationship',	  5.034117010000000D+024 ,	  2.200000000000000D+017 ,	 'm^-1'),&
       constant('joule-kelvin relationship',	  7.242971600000000D+022 ,	  6.600000000000000D+016 ,	 'K'),&
       constant('joule-kilogram relationship',	  1.112650056000000D-017 ,	  0.000000000000000D+000 ,	 'kg'),&
       constant('kelvin-atomic mass unit relationship',	  9.251086800000000D-014 ,	  8.400000000000000D-020 ,	 'u'),&
       constant('kelvin-electron volt relationship',	  8.617332400000000D-005 ,	  7.800000000000000D-011 ,	 'eV'),&
       constant('kelvin-hartree relationship',	  3.166811400000000D-006 ,	  2.900000000000000D-012 ,	 'E_h'),&
       constant('kelvin-hertz relationship',	   20836618000.0000      ,	   19000.0000000000      ,	 'Hz'),&
       constant('kelvin-inverse meter relationship',	   69.5034760000000      ,	  6.300000000000000D-005 ,	 'm^-1'),&
       constant('kelvin-joule relationship',	  1.380648800000000D-023 ,	  1.300000000000000D-029 ,	 'J'),&
       constant('kelvin-kilogram relationship',	  1.536179000000000D-040 ,	  1.400000000000000D-046 ,	 'kg'),&
       constant('kilogram-atomic mass unit relationship',	  6.022141290000000D+026 ,	  2.700000000000000D+019 ,	 'u'),&
       constant('kilogram-electron volt relationship',	  5.609588850000000D+035 ,	  1.200000000000000D+028 ,	 'eV'),&
       constant('kilogram-hartree relationship',	  2.061485968000000D+034 ,	  9.100000000000001D+026 ,	 'E_h'),&
       constant('kilogram-hertz relationship', 1356392608.000000D+040 , 6.000000000000001D+042 ,'Hz'),&
       constant('kilogram-inverse meter relationship',	  4.524438730000000D+041 ,	  2.000000000000000D+034 ,	 'm^-1'),&
       constant('kilogram-joule relationship',	  8.987551787000000D+016 ,	  0.000000000000000D+000 ,	 'J'),&
       constant('kilogram-kelvin relationship',	  6.509658200000001D+039 ,	  5.900000000000000D+033 ,	 'K'),&
       constant('lattice parameter of silicon',	  5.431020504000000D-010 ,	  8.900000000000000D-018 ,	 'm'),&
       constant('Loschmidt constant (273.15 K, 100 kPa)',	  2.651646200000000D+025 ,	  2.400000000000000D+019 ,	 'm^-3'),&
       constant('Loschmidt constant (273.15 K, 101.325 kPa)',	  2.686780500000000D+025 ,	  2.400000000000000D+019 ,	 'm^-3'),&
       constant('mag. constant',	  1.256637061400000D-006 ,	  0.000000000000000D+000 ,	 'N A^-2'),&
       constant('mag. flux quantum',	  2.067833758000000D-015 ,	  4.600000000000000D-023 ,	 'Wb'),&
       constant('molar gas constant',	   8.31446210000000      ,	  7.500000000000000D-006 ,	 'J mol^-1 K^-1'),&
       constant('molar mass constant',	  1.000000000000000D-003 ,	  0.000000000000000D+000 ,	 'kg mol^-1'),&
       constant('molar mass of carbon-12',	  1.200000000000000D-002 ,	  0.000000000000000D+000 ,	 'kg mol^-1'),&
       constant('molar Planck constant',	  3.990312717600000D-010 ,	  2.800000000000000D-019 ,	 'J s mol^-1'),&
       constant('molar Planck constant times c',	  0.119626565779000      ,	  8.399999999999999D-011 , 'J m mol^-1'),&
       constant('molar volume of ideal gas (273.15 K, 100 kPa)',  2.271095300000000D-002 ,  2.100000000000000D-008 , 'm^3 mol^-1'),&
       constant('molar volume of ideal gas (273.15 K, 101.325 kPa)', 2.241396800000000D-002 ,2.000000000000000D-008,'m^3 mol^-1'),&
       constant('molar volume of silicon',  1.205883301000000D-005 ,8.000000000000000D-013 ,'m^3 mol^-1'),&
       constant('Mo x unit',	  1.002099520000000D-013 ,	  5.300000000000000D-020 ,	 'm'),&
       constant('muon Compton wavelength',	  1.173444103000000D-014 ,	  3.000000000000000D-022 ,	 'm'),&
       constant('muon Compton wavelength over 2 pi',	  1.867594294000000D-015 ,	  4.700000000000000D-023 ,	 'm'),&
       constant('muon-electron mass ratio',	   206.768284300000      ,	  5.200000000000000D-006 ,	 ''),&
       constant('muon g factor',	  -2.00233184180000      ,	  1.300000000000000D-009 ,	 ''),&
       constant('muon mag. mom.',	 -4.490448070000000D-026 ,	  1.500000000000000D-033 ,	 'J T^-1'),&
       constant('muon mag. mom. anomaly',	  1.165920910000000D-003 ,	  6.300000000000000D-010 ,	 ''),&
       constant('muon mag. mom. to Bohr magneton ratio',	 -4.841970440000000D-003 ,	  1.200000000000000D-010 ,	 ''),&
       constant('muon mag. mom. to nuclear magneton ratio',	  -8.89059697000000      ,	  2.200000000000000D-007 ,	 ''),&
       constant('muon mass',	  1.883531475000000D-028 ,	  9.600000000000000D-036 ,	 'kg'),&
       constant('muon mass energy equivalent',	  1.692833667000000D-011 ,	  8.600000000000000D-019 ,	 'J'),&
       constant('muon mass energy equivalent in MeV',	   105.658371500000      ,	  3.500000000000000D-006 ,	 'MeV'),&
       constant('muon mass in u',	  0.113428926700000      ,	  2.900000000000000D-009 ,	 'u'),&
       constant('muon molar mass',	  1.134289267000000D-004 ,	  2.900000000000000D-012 ,	 'kg mol^-1'),&
       constant('muon-neutron mass ratio',	  0.112454517700000      ,	  2.800000000000000D-009 ,	 ''),&
       constant('muon-proton mag. mom. ratio',	  -3.18334510700000      ,	  8.400000000000000D-008 ,	 ''),&
       constant('muon-proton mass ratio',	  0.112609527200000      ,	  2.800000000000000D-009 ,	 ''),&
       constant('muon-tau mass ratio',	  5.946490000000000D-002 ,	  5.400000000000000D-006 ,	 ''),&
       constant('natural unit of action',	  1.054571726000000D-034 ,	  4.700000000000000D-042 ,	 'J s'),&
       constant('natural unit of action in eV s',	  6.582119280000000D-016 ,	  1.500000000000000D-023 ,	 'eV s'),&
       constant('natural unit of energy',	  8.187105060000000D-014 ,	  3.600000000000000D-021 ,	 'J'),&
       constant('natural unit of energy in MeV',	  0.510998928000000      ,	  1.100000000000000D-008 ,	 'MeV'),&
       constant('natural unit of length',	  3.861592680000000D-013 ,	  2.500000000000000D-022 ,	 'm'),&
       constant('natural unit of mass',	  9.109382910000000D-031 ,	  4.000000000000000D-038 ,	 'kg'),&
       constant('natural unit of mom.um',	  2.730924290000000D-022 ,	  1.200000000000000D-029 ,	 'kg m s^-1'),&
       constant('natural unit of mom.um in MeV/c',	  0.510998928000000      ,	  1.100000000000000D-008 ,	 'MeV/c'),&
       constant('natural unit of time',	  1.288088668330000D-021 ,	  8.300000000000000D-031 ,	 's'),&
       constant('natural unit of velocity',	   299792458.000000      ,	  0.000000000000000D+000 ,	 'm s^-1'),&
       constant('neutron Compton wavelength',	  1.319590906800000D-015 ,	  1.100000000000000D-024 ,	 'm'),&
       constant('neutron Compton wavelength over 2 pi',	  2.100194156800000D-016 ,	  1.700000000000000D-025 ,	 'm'),&
       constant('neutron-electron mag. mom. ratio',	  1.040668820000000D-003 ,	  2.500000000000000D-010 ,	 ''),&
       constant('neutron-electron mass ratio',	   1838.68366050000      ,	  1.100000000000000D-006 ,	 ''),&
       constant('neutron g factor',	  -3.82608545000000      ,	  9.000000000000000D-007 ,	 ''),&
       constant('neutron gyromag. ratio',	   183247179.000000      ,	   43.0000000000000      ,	 's^-1 T^-1'),&
       constant('neutron gyromag. ratio over 2 pi',	   29.1646943000000      ,	  6.900000000000000D-006 ,	 'MHz T^-1'),&
       constant('neutron mag. mom.',	 -9.662364700000001D-027 ,	  2.300000000000000D-033 ,	 'J T^-1'),&
       constant('neutron mag. mom. to Bohr magneton ratio',	 -1.041875630000000D-003 ,	  2.500000000000000D-010 ,	 ''),&
       constant('neutron mag. mom. to nuclear magneton ratio',	  -1.91304272000000      ,	  4.500000000000000D-007 ,	 ''),&
       constant('neutron mass',	  1.674927351000000D-027 ,	  7.400000000000000D-035 ,	 'kg'),&
       constant('neutron mass energy equivalent',	  1.505349631000000D-010 ,	  6.600000000000000D-018 ,	 'J'),&
       constant('neutron mass energy equivalent in MeV',	   939.565379000000      ,	  2.100000000000000D-005 ,	 'MeV'),&
       constant('neutron mass in u',	   1.00866491600000      ,	  4.300000000000000D-010 ,	 'u'),&
       constant('neutron molar mass',	  1.008664916000000D-003 ,	  4.300000000000000D-013 ,	 'kg mol^-1'),&
       constant('neutron-muon mass ratio',	   8.89248400000000      ,	  2.200000000000000D-007 ,	 ''),&
       constant('neutron-proton mag. mom. ratio',	 -0.684979340000000      ,	  1.600000000000000D-007 ,	 ''),&
       constant('neutron-proton mass difference',	  2.305573920000000D-030 ,	  7.600000000000001D-037 ,	 ''),&
       constant('neutron-proton mass difference energy equivalent',	  2.072146500000000D-013 ,	  6.799999999999999D-020 ,	 ''),&
       constant('neutron-proton mass difference energy equivalent in MeV',1.29333217000000,4.200000000000000D-007 ,''),&
       constant('neutron-proton mass difference in u',1.388449190000000D-003,4.500000000000000D-010,''),&
       constant('neutron-proton mass ratio',	   1.00137841917000      ,	  4.500000000000000D-010 ,	 ''),&
       constant('neutron-tau mass ratio',	  0.528790000000000      ,	  4.800000000000000D-005 ,	 ''),&
       constant('neutron to shielded proton mag. mom. ratio',	 -0.684996940000000      ,	  1.600000000000000D-007 ,	 ''),&
       constant('Newtonian constant of gravitation',6.673840000000000D-011 , 8.000000000000001D-015 ,'m^3 kg^-1 s^-2'),&
       constant('Newtonian constant of gravitation over h-bar c',6.708369999999999D-039 , 8.000000000000001D-043 ,'(GeV/c^2)^-2'),&
       constant('nuclear magneton',5.050783530000000D-027 ,1.100000000000000D-034 ,'J T^-1'),&
       constant('nuclear magneton in eV/T',  3.152451260500000D-008 ,	  2.200000000000000D-017 ,	 'eV T^-1'),&
       constant('nuclear magneton in inverse meters per tesla',  2.542623527000000D-002 ,	  5.600000000000000D-010 ,	 'm^-1 T^-1'),&
       constant('nuclear magneton in K/T',  3.658268200000000D-004 ,	  3.300000000000000D-010 ,	 'K T^-1'),&
       constant('nuclear magneton in MHz/T',   7.62259357000000      ,	  1.700000000000000D-007 ,	 'MHz T^-1'),&
       constant('Planck constant',	  6.626069570000000D-034 ,	  2.900000000000000D-041 ,	 'J s'),&
       constant('Planck constant in eV s',	  4.135667516000000D-015 ,	  9.100000000000000D-023 ,	 'eV s'),&
       constant('Planck constant over 2 pi',	  1.054571726000000D-034 ,	  4.700000000000000D-042 ,	 'J s'),&
       constant('Planck constant over 2 pi in eV s',	  6.582119280000000D-016 ,	  1.500000000000000D-023 ,	 'eV s'),&
       constant('Planck constant over 2 pi times c in MeV fm',	   197.326971800000      ,	  4.400000000000000D-006 ,	 'MeV fm'),&
       constant('Planck length',	  1.616199000000000D-035 ,	  9.700000000000000D-040 ,	 'm'),&
       constant('Planck mass',	  2.176510000000000D-008 ,	  1.300000000000000D-012 ,	 'kg'),&
       constant('Planck mass energy equivalent in GeV',	  1.220932000000000D+019 ,	   730000000000000.      ,	 'GeV'),&
       constant('Planck temperature',	  1.416833000000000D+032 ,	  8.500000000000000D+027 ,	 'K'),&
       constant('Planck time',	  5.391060000000000D-044 ,	  3.200000000000000D-048 ,	 's'),&
       constant('proton charge to mass quotient',	   95788335.8000000      ,	   2.10000000000000      ,	 'C kg^-1'),&
       constant('proton Compton wavelength',	  1.321409856230000D-015 ,	  9.399999999999999D-025 ,	 'm'),&
       constant('proton Compton wavelength over 2 pi',	  2.103089104700000D-016 ,	  1.500000000000000D-025 ,	 'm'),&
       constant('proton-electron mass ratio',	   1836.15267245000      ,	  7.500000000000000D-007 ,	 ''),&
       constant('proton g factor',	   5.58569471300000      ,	  4.600000000000000D-008 ,	 ''),&
       constant('proton gyromag. ratio',	   267522200.500000      ,	   6.30000000000000      ,	 's^-1 T^-1'),&
       constant('proton gyromag. ratio over 2 pi',	   42.5774806000000      ,	  1.000000000000000D-006 ,	 'MHz T^-1'),&
       constant('proton mag. mom.',	  1.410606743000000D-026 ,	  3.300000000000000D-034 ,	 'J T^-1'),&
       constant('proton mag. mom. to Bohr magneton ratio',	  1.521032210000000D-003 ,	  1.200000000000000D-011 ,	 ''),&
       constant('proton mag. mom. to nuclear magneton ratio',	   2.79284735600000      ,	  2.300000000000000D-008 ,	 ''),&
       constant('proton mag. shielding correction',	  2.569400000000000D-005 ,	  1.400000000000000D-008 ,	 ''),&
       constant('proton mass',	  1.672621777000000D-027 ,	  7.400000000000000D-035 ,	 'kg'),&
       constant('proton mass energy equivalent',	  1.503277484000000D-010 ,	  6.600000000000000D-018 ,	 'J'),&
       constant('proton mass energy equivalent in MeV',	   938.272046000000      ,	  2.100000000000000D-005 ,	 'MeV'),&
       constant('proton mass in u',	   1.00727646681200      ,	  9.000000000000000D-011 ,	 'u'),&
       constant('proton molar mass',	  1.007276466812000D-003 ,	  9.000000000000000D-014 ,	 'kg mol^-1'),&
       constant('proton-muon mass ratio',	   8.88024331000000      ,	  2.200000000000000D-007 ,	 ''),&
       constant('proton-neutron mag. mom. ratio',	  -1.45989806000000      ,	  3.400000000000000D-007 ,	 ''),&
       constant('proton-neutron mass ratio',	  0.998623478260000      ,	  4.500000000000000D-010 ,	 ''),&
       constant('proton rms charge radius',	  8.775000000000000D-016 ,	  5.100000000000000D-018 ,	 'm'),&
       constant('proton-tau mass ratio',	  0.528063000000000      ,	  4.800000000000000D-005 ,	 ''),&
       constant('quantum of circulation',	  3.636947552000000D-004 ,	  2.400000000000000D-013 ,	 'm^2 s^-1'),&
       constant('quantum of circulation times 2',	  7.273895104000000D-004 ,	  4.700000000000000D-013 ,	 'm^2 s^-1'),&
       constant('Rydberg constant',	   10973731.5685390      ,	  5.500000000000000D-005 ,	 'm^-1'),&
       constant('Rydberg constant times c in Hz',	  3.289841960364000D+015 ,	   17000.0000000000      ,	 'Hz'),&
       constant('Rydberg constant times hc in eV',	   13.6056925300000      ,	  3.000000000000000D-007 ,	 'eV'),&
       constant('Rydberg constant times hc in J',	  2.179872171000000D-018 ,	  9.600000000000000D-026 ,	 'J'),&
       constant('Sackur-Tetrode constant (1 K, 100 kPa)',	  -1.15170780000000      ,	  2.300000000000000D-006 ,	 ''),&
       constant('Sackur-Tetrode constant (1 K, 101.325 kPa)',	  -1.16487080000000      ,	  2.300000000000000D-006 ,	 ''),&
       constant('second radiation constant',	  1.438777000000000D-002 ,	  1.300000000000000D-008 ,	 'm K'),&
       constant('shielded helion gyromag. ratio',	   203789465.900000      ,	   5.10000000000000      ,	 's^-1 T^-1'),&
       constant('shielded helion gyromag. ratio over 2 pi',	   32.4341008400000      ,	  8.100000000000000D-007 ,	 'MHz T^-1'),&
       constant('shielded helion mag. mom.',	 -1.074553044000000D-026 ,	  2.700000000000000D-034 ,	 'J T^-1'),&
       constant('shielded helion mag. mom. to Bohr magneton ratio',	 -1.158671471000000D-003 ,	  1.400000000000000D-011 ,	 ''),&
       constant('shielded helion mag. mom. to nuclear magneton ratio',	  -2.12749771800000      ,	  2.500000000000000D-008 ,	 ''),&
       constant('shielded helion to proton mag. mom. ratio',	 -0.761766558000000      ,	  1.100000000000000D-008 ,	 ''),&
       constant('shielded helion to shielded proton mag. mom. ratio',	 -0.761786131300000      ,	  3.300000000000000D-009 ,	 ''),&
       constant('shielded proton gyromag. ratio',	   267515326.800000      ,	   6.60000000000000      ,	 's^-1 T^-1'),&
       constant('shielded proton gyromag. ratio over 2 pi',	   42.5763866000000      ,	  1.000000000000000D-006 ,	 'MHz T^-1'),&
       constant('shielded proton mag. mom.',	  1.410570499000000D-026 ,	  3.500000000000000D-034 ,	 'J T^-1'),&
       constant('shielded proton mag. mom. to Bohr magneton ratio',	  1.520993128000000D-003 ,	  1.700000000000000D-011 ,	 ''),&
       constant('shielded proton mag. mom. to nuclear magneton ratio',	   2.79277559800000      ,	  3.000000000000000D-008 ,	 ''),&
       constant('speed of light in vacuum',	   299792458.000000      ,	  0.000000000000000D+000 ,	 'm s^-1'),&
       constant('standard acceleration of gravity',	   9.80665000000000      ,	  0.000000000000000D+000 ,	 'm s^-2'),&
       constant('standard atmosphere',	   101325.000000000      ,	  0.000000000000000D+000 ,	 'Pa'),&
       constant('standard-state pressure',	   100000.000000000      ,	  0.000000000000000D+000 ,	 'Pa'),&
       constant('Stefan-Boltzmann constant',	  5.670373000000000D-008 ,	  2.100000000000000D-013 ,	 'W m^-2 K^-4'),&
       constant('tau Compton wavelength',	  6.977870000000000D-016 ,	  6.300000000000000D-020 ,	 'm'),&
       constant('tau Compton wavelength over 2 pi',	  1.110560000000000D-016 ,	  9.999999999999999D-021 ,	 'm'),&
       constant('tau-electron mass ratio',	   3477.15000000000      ,	  0.310000000000000      ,	 ''),&
       constant('tau mass',	  3.167470000000000D-027 ,	  2.900000000000000D-031 ,	 'kg'),&
       constant('tau mass energy equivalent',	  2.846780000000000D-010 ,	  2.600000000000000D-014 ,	 'J'),&
       constant('tau mass energy equivalent in MeV',	   1776.82000000000      ,	  0.160000000000000      ,	 'MeV'),&
       constant('tau mass in u',	   1.90749000000000      ,	  1.700000000000000D-004 ,	 'u'),&
       constant('tau molar mass',	  1.907490000000000D-003 ,	  1.700000000000000D-007 ,	 'kg mol^-1'),&
       constant('tau-muon mass ratio',	   16.8167000000000      ,	  1.500000000000000D-003 ,	 ''),&
       constant('tau-neutron mass ratio',	   1.89111000000000      ,	  1.700000000000000D-004 ,	 ''),&
       constant('tau-proton mass ratio',	   1.89372000000000      ,	  1.700000000000000D-004 ,	 ''),&
       constant('Thomson cross section',	  6.652458734000000D-029 ,	  1.300000000000000D-037 ,	 'm^2'),&
       constant('triton-electron mass ratio',	   5496.92152670000      ,	  5.000000000000000D-006 ,	 ''),&
       constant('triton g factor',	   5.95792489600000      ,	  7.600000000000001D-008 ,	 ''),&
       constant('triton mag. mom.',	  1.504609447000000D-026 ,	  3.800000000000000D-034 ,	 'J T^-1'),&
       constant('triton mag. mom. to Bohr magneton ratio',	  1.622393657000000D-003 ,	  2.100000000000000D-011 ,	 ''),&
       constant('triton mag. mom. to nuclear magneton ratio',	   2.97896244800000      ,	  3.800000000000000D-008 ,	 ''),&
       constant('triton mass',	  5.007356300000000D-027 ,	  2.200000000000000D-034 ,	 'kg'),&
       constant('triton mass energy equivalent',	  4.500387410000000D-010 ,	  2.000000000000000D-017 ,	 'J'),&
       constant('triton mass energy equivalent in MeV',	   2808.92100500000      ,	  6.200000000000000D-005 ,	 'MeV'),&
       constant('triton mass in u',	   3.01550071340000      ,	  2.500000000000000D-009 ,	 'u'),&
       constant('triton molar mass',	  3.015500713400000D-003 ,	  2.500000000000000D-012 ,	 'kg mol^-1'),&
       constant('triton-proton mass ratio',	   2.99371703080000      ,	  2.500000000000000D-009 ,	 ''),&
       constant('unified atomic mass unit',	  1.660538921000000D-027 ,	  7.300000000000000D-035 ,	 'kg'),&
       constant('von Klitzing constant',	   25812.8074434000      ,	  8.399999999999999D-006 ,	 'ohm'),&
       constant('weak mixing angle',	  0.222300000000000      ,	  2.100000000000000D-003 ,	 ''),&
       constant('Wien frequency displacement law constant',	   58789254000.0000      ,	   53000.0000000000      ,	 'Hz K^-1')&
       ]



  public :: value_constant
  public :: name_constant
  public :: units_constant
  public :: parse_constant
  public :: print_constant
contains

  function value_constant(key) result(value_)
    integer        :: key    
    type(constant) :: const
    real(8)        :: value_
    if(key>nmax_constant)then
       print*,"No constant associated to this number. nmax_constant=",nmax_constant
       value_=0.d0
       return
    endif
    value_=fundamental_constant(key)%val
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


  ! subroutine get_constant_from_file(fileIN,Nmax,fileOUT)
  !   character(len=*)   :: fileIN
  !   integer            :: i,Nmax
  !   character(len=*)   :: fileOUT
  !   type(constant)     :: const
  !   character(len=160) ::string
  !   open(10,file=fileIN)
  !   open(11,file=fileOUT)
  !   do i=1,nmax
  !      read(10,"(A128)")string
  !      const =  parse_constant(string)
  !      write(*,"(A20,A55,A2,EN24.12,A1,EN24.12,A2,A20,A4)")"constant("",const%name,"",",const%val",",const%uncert,","",const%units,","),& "
  !      write(11,"(A20,A55,A2,EN24.12,A1,EN24.12,A2,A20,A4)")"constant("",const%name,"",",const%val,",",const%uncert,","",const%units,","),& "
  !   enddo
  ! end subroutine get_constant_from_file


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
       if ( ch /= " " .and. ch /= tab ) then
          put = put + 1
          s(put:put) = ch
       end if
    end do
    s(put+1:s_length) = " "
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
    !    S = "Fred is not a jerk!"
    !    call s_chop ( S, 9, 12 )
    !    S = "Fred is a jerk!    "
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
    s(s_length+ilo2-ihi2:s_length) = " "
  end subroutine s_chop

end module CONSTANTS
