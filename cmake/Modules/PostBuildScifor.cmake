FUNCTION(BUILD_ENV_MODULE FILE)
  MESSAGE(STATUS "${BoldYellow}Writing tmp env. module to ${FILE}${ColourReset}")
  FILE(WRITE  ${FILE}  "#%Modules\n")
  FILE(APPEND ${FILE} "set project ${PROJECT_NAME}\n")
  FILE(APPEND ${FILE} "set root ${PREFIX}\n")
  FILE(APPEND ${FILE} "set plat ${FC_PLAT}\n")
  FILE(APPEND ${FILE} "set version \"${VERSION} (${FC_PLAT})\"\n")
  FILE(APPEND ${FILE} "set compiler ${CMAKE_Fortran_COMPILER}\n")
  FILE(READ   ${SF_ENV}/module CONTENTS)
  FILE(APPEND ${FILE} "${CONTENTS}")
ENDFUNCTION()



FUNCTION(BUILD_CONFIGVARS FILE)
  MESSAGE(STATUS "${BoldYellow}Setup configvars.sh${ColourReset}")
  FILE(READ   ${SF_BIN}/configvars.sh CONTENTS)
  FILE(WRITE  ${FILE} "${CONTENTS}")
  FILE(APPEND ${FILE} "add_library_to_system ${CMAKE_INSTALL_PREFIX}\n")
ENDFUNCTION()


FUNCTION(BUILD_PKCONFIG FILE)
  MESSAGE(STATUS "${BoldYellow}Writing tmp pk-config file to ${FILE}${ColourReset}")
  FILE(WRITE ${FILE} "prefix=${CMAKE_INSTALL_PREFIX}\n")
  FILE(READ   ${SF_ETC}/${PROJECT_NAME}.pc CONTENTS)
  FILE(APPEND ${FILE} "${CONTENTS}")
  FILE(APPEND ${FILE} "Version:${VERSION}\n")
ENDFUNCTION()
