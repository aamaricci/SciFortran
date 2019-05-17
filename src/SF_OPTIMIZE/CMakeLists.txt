# SF_OPTIMIZE library
file(GLOB SF_OPTIMIZE_LOCAL_SRC
  optimize_broyden_routines.f90
  minimize.f
  cgfam.f
  cgsearch.f
  optimize_cgfit_routines.f90
  timer.f
  linpack.f
  lbfgsb.f
  ROOT_FINDING.f90
  OPTIMIZATION.f90
  SF_OPTIMIZE.f90
  )
SET(SF_OPTIMIZE_src ${SF_OPTIMIZE_LOCAL_SRC})
ADD_LIBRARY(SF_OPTIMIZELIB OBJECT ${SF_OPTIMIZE_src})
