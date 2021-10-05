#Build the Library module for environment modules
FILE(REMOVE_RECURSE ${LIB_TMP_ETC}/modules)

#Create default .version file if following conditions are verified:
#1. master branch
#2. no debug
#3. default prefix
IF( USE_DEFAULT_MODULE_NAME AND (NOT "${CMAKE_BUILD_TYPE}" MATCHES "DEBUG") AND (GIT_BRANCH MATCHES "master"))
  SET(TMP_VER_MODULE_FILE ${LIB_TMP_ETC}/modules/${VERSION_PATH}/.version)
  CONFIGURE_FILE(${LIB_ENV}/version.in ${TMP_VER_MODULE_FILE}  @ONLY)
  MESSAGE(STATUS "${Red}Version file${ColourReset}: ${TMP_VER_MODULE_FILE}")
ENDIF()

#Build the library module for environment modules
SET(TMP_MODULE_NAME "${MODULE_NAME}" CACHE PATH "Prefix prepended to install directories")
SET(TMP_ENV_MODULE_FILE ${LIB_TMP_ETC}/modules/${TMP_MODULE_NAME})
CONFIGURE_FILE(${LIB_ENV}/module.in ${TMP_ENV_MODULE_FILE} @ONLY)
MESSAGE(STATUS "${Red}Module file${ColourReset}: ${MODULE_NAME}")


#Build the user CONFIG scripts (sourced in user shell config file, i.e. .bashrc)
SET(USER_CONFIG_SCRIPT ${PROJECT_NAME}_config_user.sh)
SET(TMP_CONFIGVARS_USER_FILE ${LIB_TMP_ETC}/${USER_CONFIG_SCRIPT})
CONFIGURE_FILE(${LIB_ETC}/${PROJECT_NAME}_config_user.sh.in ${TMP_CONFIGVARS_USER_FILE} @ONLY)

SET(GLOBAL_CONFIG_SCRIPT ${PROJECT_NAME}_config_global.sh)
SET(TMP_CONFIGVARS_GLOBAL_FILE ${LIB_TMP_ETC}/${GLOBAL_CONFIG_SCRIPT})
CONFIGURE_FILE(${LIB_ETC}/${PROJECT_NAME}_config_global.sh.in ${TMP_CONFIGVARS_GLOBAL_FILE} @ONLY)


FILE(WRITE  ${LIB_TMP_VER}  "${VERSION}\n")

INSTALL(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION ${LIB_TARGET_INC})

INSTALL(CODE "execute_process(COMMAND \"${CMAKE_COMMAND}\" -E rm -f ${LIB_VERSION_FILE})")

INSTALL(TARGETS ${PROJECT_NAME} DESTINATION ${LIB_TARGET_LIB})

INSTALL(DIRECTORY ${LIB_TMP_ETC}/ DESTINATION ${LIB_TARGET_ETC})

INSTALL(FILES ${TMP_CONFIGVARS_USER_FILE} DESTINATION ${LIB_TARGET_BIN}/
  PERMISSIONS ${PERMISSION_777} SETUID)

INSTALL(FILES ${TMP_CONFIGVARS_GLOBAL_FILE} DESTINATION ${LIB_TARGET_BIN}/
  PERMISSIONS ${PERMISSION_777} SETUID)

INSTALL(FILES ${LIB_TMP_VER} DESTINATION ${LIB_TARGET_DIR} 
  PERMISSIONS ${PERMISSION_777} SETUID)

INSTALL(FILES ${LIB_TARGET_ETC}/${PROJECT_NAME}.pc DESTINATION $ENV{HOME}/.pkgconfig.d/
  PERMISSIONS ${PERMISSION_777} SETUID)

INSTALL(DIRECTORY ${LIB_TARGET_ETC}/modules/ DESTINATION $ENV{HOME}/.modules.d)


MESSAGE( STATUS "${Red}Library version:${ColourReset} ${VERSION}")
MESSAGE( STATUS "${Red}Library will be installed in:${ColourReset} ${CMAKE_INSTALL_PREFIX}")
MESSAGE( STATUS "
>> ${Red}TO CONCLUDE INSTALLATION${ColourReset} <<
Compile with:
$ make
Install with:
$ make install

Uninstall with:
$ make uninstall
")

INSTALL(CODE "MESSAGE(
\"
ADD LIBRARY TO YOUR SYSTEM: 
Pick ONE method below [or add it in your bash profile, e.g. ~/.bashrc]:
${Yellow}Method 1: use the provided ${PROJECT_NAME} environment module${ColourReset}:
   $ module use $HOME/.modules.d
   $ module load ${PROJECT_NAME}/${FC_PLAT}

${Yellow}Method 2: source the config script${ColourReset}:
   $ source ${LIB_TARGET_BIN}/${USER_CONFIG_FILE}

${Yellow}Method 3: use pkg-config with the provided ${PROJECT_NAME}.pc${ColourReset}:
   $ export PKG_CONFIG_PATH=${LIB_TARGET_ETC}/:$PKG_CONFIG_PATH
   $ pkg-config --cflags --libs ${PROJECT_NAME}

${Yellow}Method ADMIN: Add this line to the system shell configuration file, e.g. /etc/bash.bashrc${ColourReset}
   $ source ${LIB_TARGET_BIN}/${GLOBAL_CONFIG_FILE}
\")
")



# Add a distclean target to the Makefile
ADD_CUSTOM_TARGET(distclean 
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_MODULE_PATH}/DistClean.cmake
)


# Uninstall target
if(NOT TARGET uninstall)
  CONFIGURE_FILE(
    "${CMAKE_MODULE_PATH}/ConfigUninstall.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

  ADD_CUSTOM_TARGET(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake
    "${PROJECT_NAME}/${TMP_MODULE_NAME}" "${PROJECT_NAME}.pc" )
  
ENDIF()
