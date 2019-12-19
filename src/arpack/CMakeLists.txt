SET(UTIL_OBJS
  util/icnteq.f
  util/icopy.f
  util/iset.f
  util/iswap.f
  util/ivout.f
  util/second.f)

SET(UTIL_SOBJ
  util/svout.f
  util/smout.f)
SET(UTIL_DOBJ
  util/dvout.f
  util/dmout.f)
SET(UTIL_COBJ
  util/cvout.f
  util/cmout.f)
SET(UTIL_ZOBJ
  util/zvout.f
  util/zmout.f)


SET(P_UTIL_OBJS
  parpack/src/pcontext.f
  parpack/util/pivout.f)
SET(P_UTIL_SOBJ
  parpack/util/psvout.f
  parpack/util/psmout.f)
SET(P_UTIL_DOBJ
  parpack/util/pdvout.f
  parpack/util/pdmout.f)
SET(P_UTIL_COBJ
  parpack/util/pcvout.f
  parpack/util/pcmout.f)
SET(P_UTIL_ZOBJ
  parpack/util/pzvout.f
  parpack/util/pzmout.f)


SET(SRC_SOBJ
  src/sgetv0.f
  ##src/slaqrb.f
  src/sstqrb.f
  src/ssortc.f
  src/ssortr.f
  src/sstatn.f
  src/sstats.f 
  src/snaitr.f
  src/snapps.f
  src/snaup2.f
  src/snaupd.f
  src/snconv.f
  src/sneigh.f
  src/sngets.f 
  src/ssaitr.f
  src/ssapps.f
  src/ssaup2.f
  src/ssaupd.f
  src/ssconv.f
  src/sseigt.f
  src/ssgets.f 
  src/sneupd.f
  src/sseupd.f
  src/ssesrt.f)

SET(SRC_DOBJ
  src/dgetv0.f
  ##src/dlaqrb.f
  src/dstqrb.f
  src/dsortc.f
  src/dsortr.f
  src/dstatn.f
  src/dstats.f
  src/dnaitr.f
  src/dnapps.f
  src/dnaup2.f
  src/dnaupd.f
  src/dnconv.f
  src/dneigh.f
  src/dngets.f 
  src/dsaitr.f
  src/dsapps.f
  src/dsaup2.f
  src/dsaupd.f
  src/dsconv.f
  src/dseigt.f
  src/dsgets.f 
  src/dneupd.f
  src/dseupd.f
  src/dsesrt.f)

SET(SRC_COBJ
  src/cnaitr.f
  src/cnapps.f
  src/cnaup2.f
  src/cnaupd.f
  src/cneigh.f
  src/cneupd.f
  src/cngets.f
  src/cgetv0.f
  src/csortc.f
  src/cstatn.f)

SET(SRC_ZOBJ
  src/znaitr.f
  src/znapps.f
  src/znaup2.f
  src/znaupd.f
  src/zneigh.f
  src/zneupd.f
  src/zngets.f
  src/zgetv0.f
  src/zsortc.f
  src/zstatn.f)


SET(P_SRC_SOBJ
  parpack/src/psgetv0.f  
  parpack/src/psnaitr.f
  parpack/src/psnapps.f
  parpack/src/psnaup2.f
  parpack/src/psnaupd.f
  parpack/src/psneigh.f
  parpack/src/psngets.f 
  parpack/src/pssaitr.f
  parpack/src/pssapps.f
  parpack/src/pssaup2.f
  parpack/src/pssaupd.f
  parpack/src/psseigt.f
  parpack/src/pssgets.f 
  parpack/src/psneupd.f
  parpack/src/psseupd.f
  parpack/src/pslarnv.f
  parpack/src/pslamch10.f
  parpack/src/psnorm2.f)

SET(P_SRC_DOBJ
  parpack/src/pdgetv0.f 
  parpack/src/pdnaitr.f
  parpack/src/pdnapps.f
  parpack/src/pdnaup2.f
  parpack/src/pdnaupd.f
  parpack/src/pdneigh.f
  parpack/src/pdngets.f 
  parpack/src/pdsaitr.f
  parpack/src/pdsapps.f
  parpack/src/pdsaup2.f
  parpack/src/pdsaupd.f
  parpack/src/pdseigt.f
  parpack/src/pdsgets.f 
  parpack/src/pdneupd.f
  parpack/src/pdseupd.f
  parpack/src/pdlarnv.f
  parpack/src/pdlamch10.f
  parpack/src/pdnorm2.f)

SET(P_SRC_COBJ
  parpack/src/pcnaitr.f
  parpack/src/pcnapps.f
  parpack/src/pcnaup2.f
  parpack/src/pcnaupd.f
  parpack/src/pcneigh.f 
  parpack/src/pcneupd.f
  parpack/src/pcngets.f
  parpack/src/pcgetv0.f
  parpack/src/pscnorm2.f
  parpack/src/pclarnv.f)

SET(P_SRC_ZOBJ
  parpack/src/pznaitr.f
  parpack/src/pznapps.f
  parpack/src/pznaup2.f
  parpack/src/pznaupd.f
  parpack/src/pzneigh.f 
  parpack/src/pzneupd.f
  parpack/src/pzngets.f
  parpack/src/pzgetv0.f
  parpack/src/pdznorm2.f
  parpack/src/pzlarnv.f)


SET(SINGLE    ${UTIL_OBJS} ${UTIL_SOBJ}  ${SRC_SOBJ})
SET(DOUBLE    ${UTIL_OBJS} ${UTIL_DOBJ}  ${SRC_DOBJ})
SET(COMPLEX   ${UTIL_OBJS} ${UTIL_COBJ}  ${SRC_COBJ})
SET(COMPLEX16 ${UTIL_OBJS} ${UTIL_ZOBJ}  ${SRC_ZOBJ})

IF(USE_MPI)
  SET(SINGLE    ${SINGLE}    ${P_UTIL_OBJS} ${P_UTIL_SOBJ}  ${P_SRC_SOBJ})
  SET(DOUBLE    ${DOUBLE}    ${P_UTIL_OBJS} ${P_UTIL_DOBJ}  ${P_SRC_DOBJ})
  SET(COMPLEX   ${COMPLEX}   ${P_UTIL_OBJS} ${P_UTIL_COBJ}  ${P_SRC_COBJ})
  SET(COMPLEX16 ${COMPLEX16} ${P_UTIL_OBJS} ${P_UTIL_ZOBJ}  ${P_SRC_ZOBJ})  
ENDIF(USE_MPI)

SET(PATCH "zdotc.patch")
MESSAGE(STATUS "Testing P-ARPACK patch: ${PATCH}")
EXECUTE_PROCESS(
  COMMAND patch -Rsfp0 --dry-run 
  INPUT_FILE "${PATCH}"
  WORKING_DIRECTORY "${SF_SRC}/arpack"
  OUTPUT_VARIABLE OUTPUT
  RESULT_VARIABLE RESULT)
IF(RESULT EQUAL 0)
  MESSAGE(STATUS "Patch applied: ${PATCH}")
ELSE()
  MESSAGE(STATUS "Patch not applied. Applying: ${PATCH}")
  EXECUTE_PROCESS(
    COMMAND patch -s -p0 
    INPUT_FILE "${PATCH}"
    WORKING_DIRECTORY "${SF_SRC}/arpack"
    OUTPUT_VARIABLE OUTPUT
    RESULT_VARIABLE RESULT2)
  IF(RESULT2 EQUAL 0)
    MESSAGE(STATUS "Patch successfully applied: ${PATCH}")
  ELSE()
    MESSAGE(FATAL_ERROR "Error applying patch ${PATCH}")
  ENDIF()
ENDIF()


# SET(SCRIPT "zdotc.sh")
# EXECUTE_PROCESS(
#   COMMAND bash "${SCRIPT}"
#   WORKING_DIRECTORY "${SF_SRC}/arpack"
#   OUTPUT_VARIABLE OUTPUT
#   RESULT_VARIABLE RESULT)
# IF(RESULT EQUAL 0)
#   MESSAGE(STATUS "Patch successfully applied: ${PATCH}")
# ELSE()
#   MESSAGE(FATAL_ERROR "Error applying patch ${PATCH}")
# ENDIF()

set(ARPACK_LOCAL_SRC ${DOUBLE} ${COMPLEX16})

SET(ARPACK_src ${ARPACK_LOCAL_SRC})

ADD_LIBRARY(ARPACKLIB OBJECT ${ARPACK_src})
