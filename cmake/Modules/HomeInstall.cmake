SET(CMAKE_SOURCE_DIR ${CMAKE_ARGV3})
SET(CMAKE_INSTALL_PREFIX ${CMAKE_ARGV4})
SET(PROJECT_NAME ${CMAKE_ARGV5})

SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules")
SET(SF_TARGET_LIB ${CMAKE_INSTALL_PREFIX}/lib)
SET(SF_TARGET_INC ${CMAKE_INSTALL_PREFIX}/include)
SET(SF_TARGET_ETC ${CMAKE_INSTALL_PREFIX}/etc)
SET(SF_TARGET_BIN ${CMAKE_INSTALL_PREFIX}/bin)
SET(PKCONFIG_FILE ${SF_TARGET_ETC}/${PROJECT_NAME}.pc)


INCLUDE(${CMAKE_MODULE_PATH}/ColorsMsg.cmake)

FILE(INSTALL ${SF_TARGET_ETC}/modules/ DESTINATION $ENV{HOME}/.modules.d
  USE_SOURCE_PERMISSIONS)

FILE(INSTALL ${PKCONFIG_FILE} DESTINATION $ENV{HOME}/.pkgconfig.d
  USE_SOURCE_PERMISSIONS)


MESSAGE( STATUS "After installation usage options:")
MESSAGE( STATUS
  "
${Red} GLOBAL INSTALLATION${ColourReset}:
${Yellow} Add this line to the system shell configuration file (i.e. /etc/bash.bashrc)${ColourReset}:
   source ${DT_TARGET_BIN}/scifor_config_global.sh

${Red} USER INSTALLATION${ColourReset}:
Pick ONE choice or add in your bash profile (i.e. ~/.bashrc):
${Yellow}A. source the config script: ${SF_TARGET_BIN}/dmft_tools_config_user.sh${ColourReset}:
   source ${DT_TARGET_BIN}/scifor_config_user.sh

${Yellow}B. use the provided ${PROJECT_NAME} environment module in ${SF_TARGET_ETC}${ColourReset}:
   module use $HOME/.modules.d
   module load ${PROJECT_NAME}/${FC_PLAT}

${Yellow}C. use pkg-config with the provided ${PROJECT_NAME}.pc in ${SF_TARGET_ETC}${ColourReset}:
   export PKG_CONFIG_PATH=${DT_TARGET_ETC}/:$PKG_CONFIG_PATH
   pkg-config --cflags --libs ${PROJECT_NAME} (to get lib info)
")
