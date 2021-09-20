##################################################
# SCIFOR PROJECT
##################################################
CMAKE_MINIMUM_REQUIRED(VERSION 3.0.0)
PROJECT(scifor Fortran)
SET(VERSION 4.7.3)


MESSAGE(STATUS "OS: ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_VERSION}")
SET(BUILD_TYPE "RELEASE" CACHE STRING "Build type")
STRING(TOUPPER "${BUILD_TYPE}" BUILD_TYPE)
MESSAGE(STATUS "BUILD_TYPE=${BUILD_TYPE}")
SET(CMAKE_BUILD_TYPE ${BUILD_TYPE} CACHE INTERNAL "Build type" FORCE)

# Add our local CMAKE modules to the module path
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

INCLUDE(${CMAKE_MODULE_PATH}/ColorsMsg.cmake)

# Check that Fortran 90 is supported
IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
   MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
ENDIF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)



IF( (${CMAKE_Fortran_COMPILER_ID} MATCHES Intel) OR (${CMAKE_Fortran_COMPILER_ID} MATCHES GNU))
  MESSAGE(STATUS "Fortran Compiler id   = ${CMAKE_Fortran_COMPILER_ID}")
  MESSAGE(STATUS "Fortran Compiler ver. = ${CMAKE_Fortran_COMPILER_VERSION}")
ELSEIF()
  MESSAGE(FATAL_ERROR "Unsupported Fortran compiler (use Intel or GNU). Try export FC=<your FC compiler> ")
ENDIF()  

############################################################
# Set Fortran options based on BUILD_TYPE and FC ID
############################################################
# -mcmodel=large  this is to remove the 2Gb limit of virtual memory allocation
if(CMAKE_Fortran_COMPILER_ID MATCHES GNU) # this is gfortran
  SET(CMAKE_Fortran_MODDIR_FLAG   "-J")
  SET(CMAKE_Fortran_FLAGS         "-cpp -ffree-line-length-none -fPIC -w ")
  IF(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 10.0.0)
    SET(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -Wno-argument-mismatch")
  ELSE()
    SET(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
  ENDIF()
  SET(CMAKE_Fortran_FLAGS_TESTING "-O2 -funroll-loops")
  SET(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -p -g -Wsurprising -Waliasing -fwhole-file -fcheck=all -fbacktrace -fbounds-check")
  SET(CMAKE_Fortran_FLAGS_RELEASE "-O3   -funroll-loops")   
elseif(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  SET(CMAKE_Fortran_MODDIR_FLAG   "-module ") #remember the ending white space here 
  SET(CMAKE_Fortran_FLAGS         "-fpp")
  SET(CMAKE_Fortran_FLAGS_TESTING "-O2 -ftz")
  SET(CMAKE_Fortran_FLAGS_DEBUG   "-p -O0 -g -fpe1 -warn -debug extended -traceback -check all,noarg_temp_created")
  SET(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ftz")

elseif(CMAKE_Fortran_COMPILER_ID MATCHES G95)
  SET(CMAKE_Fortran_MODDIR_FLAG   "-fmod=")
  SET(CMAKE_Fortran_FLAGS         "-cpp")
  SET(CMAKE_Fortran_FLAGS_TESTING "-O1  -fsloppy-char")
  SET(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -ftrace=full")
  SET(CMAKE_Fortran_FLAGS_RELEASE "-O3")
  
elseif(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
  SET(CMAKE_Fortran_MODDIR_FLAG   "-module ")
  SET(CMAKE_Fortran_FLAGS         "-")
  SET(CMAKE_Fortran_FLAGS         "")
  SET(CMAKE_Fortran_FLAGS_DEBUG   "-g -O0 -Mframe")
  SET(CMAKE_Fortran_FLAGS_RELEASE "-O3 -mcmodel=medium -fast -Munroll")
endif()

IF( "${BUILD_TYPE}" MATCHES "DEBUG")
  MESSAGE(STATUS "Fortran Compiler options = ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_DEBUG}")
ELSEIF("${BUILD_TYPE}" MATCHES "TESTING")
  MESSAGE(STATUS "Fortran Compiler options = ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_TESTING}")
ELSEIF("${BUILD_TYPE}" MATCHES "RELEASE")
  MESSAGE(STATUS "Fortran Compiler options = ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELEASE}")
ENDIF()



# Set some options the user may choose
OPTION(USE_MPI "Use the MPI library for parallelization" ON)
OPTION(VERBOSE "Build a verbose Makefile" OFF)
OPTION(WITH_BLAS_LAPACK "Build using internal Blas/Lapack" OFF)
OPTION(WITH_SCALAPACK "Build ScaLapack support" ON)
IF(VERBOSE)
  SET(CMAKE_VERBOSE_MAKEFILE ON)
ELSE(VERBOSE)
  SET(CMAKE_VERBOSE_MAKEFILE OFF)
ENDIF(VERBOSE)



############################################################
# Look for MPI if needed
############################################################
# reSet compiler to MPI is required 
IF(USE_MPI)
  FIND_PACKAGE(MPI REQUIRED)
  IF(MPI_Fortran_FOUND)
    SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
    MESSAGE(STATUS "${Yellow}Set Fortran compiler FC to ${ColourReset}${CMAKE_Fortran_COMPILER}, ID=${CMAKE_Fortran_COMPILER_ID}")
    SET(MPI_CPP "MPI")		#pre-processing option
  ELSE()
    MESSAGE(FATAL_ERROR "MPI Found but No MPI-Fortran compiler can be determined.")    
  ENDIF()
ELSE(USE_MPI)
  SET(MPI_CPP "")
ENDIF(USE_MPI)








############################################################
# Check for .git and extract the current SHA1 revision.
############################################################
SET(GIT_SHA1 0)
SET(GIT_BRANCH "master")
SET(GIT_TAG ${VERSION})
IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)
  SET(GIT_REPO TRUE)
  FIND_PACKAGE(Git)
  IF(GIT_FOUND)
    EXECUTE_PROCESS(
      COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      OUTPUT_VARIABLE GIT_SHA1
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    MESSAGE( STATUS "${Yellow}Git SHA1: ${ColourReset}${GIT_SHA1}" )
    #
    EXECUTE_PROCESS(
      COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      OUTPUT_VARIABLE GIT_BRANCH
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    MESSAGE( STATUS "${Yellow}Git branch: ${ColourReset}${GIT_BRANCH}" )
    #
    EXECUTE_PROCESS(
      COMMAND ${GIT_EXECUTABLE} describe --tags
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      OUTPUT_VARIABLE GIT_TAG
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    MESSAGE(STATUS "${Yellow}Git tag: ${ColourReset}${GIT_TAG}")
  ENDIF()
ELSE()
  SET(GIT_REPO FALSE)
ENDIF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)
SET(VERSION ${GIT_TAG})    	#Update version


############################################################
# Check for Blas/Lapack/ScaLapack from different sources
############################################################
# Try to find MKL: IF found use it ELSE try finding Blas/Lapack
# in the OS, ELSE using internal Blas/Lapack
IF(NOT WITH_BLAS_LAPACK)
  FIND_PACKAGE(MKL QUIET)
  IF(MKL_FOUND)
    MESSAGE(STATUS "MKL found at: ${MKL_ROOT_DIR}")
    SET(BLA_LIBRARIES ${MKL_LIBRARIES})    
  ELSE()
    FIND_PACKAGE(BLAS QUIET)
    FIND_PACKAGE(LAPACK QUIET)
    IF(BLAS_FOUND AND LAPACK_FOUND)
      MESSAGE(STATUS "BLAS found as: ${BLAS_LIBRARIES}")
      MESSAGE(STATUS "LAPACK found as: ${LAPACK_LIBRARIES}")
      SET(BLA_LIBRARIES ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    ELSE()
      SET(WITH_BLAS_LAPACK "ON")
    ENDIF()
  ENDIF()
ENDIF()

SET(SCALAPACK_FOUND FALSE)
IF(WITH_SCALAPACK)
  IF(MKL_SCALAPACK_FOUND)
    MESSAGE(STATUS "MKL supports Scalapack+Blacs")
    SET(SCALAPACK_FOUND TRUE)
  ELSE()
    FIND_PACKAGE(SCALAPACK QUIET)	#TRY TO FIND SCALAPACK IN THE SYSTEM
    IF(SCALAPACK_FOUND)
      MESSAGE(STATUS "Scalapack found: ${SCALAPACK_LIBRARIES}")
    ENDIF()
  ENDIF()  
ENDIF()



#set MPI preprocessing defition:
IF(USE_MPI)
  ADD_DEFINITIONS(-D_MPI)
ELSE(USE_MPI)
  ADD_DEFINITIONS(-D_)
ENDIF(USE_MPI)

IF((WITH_SCALAPACK) AND (SCALAPACK_FOUND))
  ADD_DEFINITIONS(-D_SCALAPACK)
ENDIF()


IF( "${BUILD_TYPE}" MATCHES "DEBUG")
  ADD_DEFINITIONS(-D_DEBUG)
ENDIF()




############################################################
# Define the actual files and folders that make up the build
############################################################
SET(USER_HOME $ENV{HOME})
SET(USER $ENV{USER})

STRING(TOLOWER "${CMAKE_Fortran_COMPILER_ID}" FC_ID)
STRING(TOLOWER "${CMAKE_BUILD_TYPE}" BUILD_TYPE)

#default prefix is $HOME/opt/<libname>/<fc_id>/[<git_branch>[/<debug>]]/<version>
SET(PREFIX_DEF_LOC "$ENV{HOME}/opt")
SET(PREFIX_PROJ "${PROJECT_NAME}")
SET(PREFIX_PATH "${FC_ID}")

#user can change the default location of the library but not the remaining naming conventions 
SET(PREFIX  "${PREFIX_DEF_LOC}" CACHE PATH "Prefix prepended to install directories")

#if prefix has not be changed use default module name conventions
IF(PREFIX STREQUAL "${PREFIX_DEF_LOC}")
  SET(USE_DEFAULT_MODULE_NAME TRUE)
ENDIF()

#if not master branch, include simplified branch name
IF( (NOT GIT_BRANCH MATCHES "master") )
    SET(PREFIX_PATH  "${PREFIX_PATH}/${GIT_BRANCH}")
ENDIF()

#If DEBUG, add /debug to 
IF("${CMAKE_BUILD_TYPE}" MATCHES "DEBUG")
  SET(PREFIX_PATH  "${PREFIX_PATH}/${BUILD_TYPE}")
ENDIF()

#set default prefix:
SET(FULL_VER "${PREFIX_PATH}/${VERSION}")
SET(VERSION_PATH "${PREFIX_PROJ}/${PREFIX_PATH}")
SET(DEFAULT_PREFIX "${PREFIX}/${PREFIX_PROJ}/${PREFIX_PATH}/${VERSION}")
SET(CMAKE_INSTALL_PREFIX "${DEFAULT_PREFIX}" CACHE INTERNAL "Prefix prepended to install directories" FORCE)

#Now set module name corresponding to given prefix:
#if prefix is default nothing to be done.
#if user defined prefix prepend
IF(USE_DEFAULT_MODULE_NAME)
  SET(MODULE_NAME "${PREFIX_PROJ}/${PREFIX_PATH}/${VERSION}")
ELSE()
  SET(MODULE_NAME "${PREFIX_PROJ}/user_prefix/${PREFIX_PATH}/${VERSION}")
ENDIF()


# Define the library name
SET(SCIFORLIB scifor)

# Define some directories
SET(SF_SRC ${CMAKE_SOURCE_DIR}/src)
SET(SF_ETC ${CMAKE_SOURCE_DIR}/etc)
SET(SF_ENV ${CMAKE_SOURCE_DIR}/etc/environment_modules)
SET(SF_BIN ${CMAKE_SOURCE_DIR}/bin)
SET(SF_TMP_INC ${CMAKE_BINARY_DIR}/include)
SET(SF_TMP_ETC ${CMAKE_BINARY_DIR}/etc)
SET(SF_TMP_BIN ${CMAKE_BINARY_DIR}/bin)
SET(SF_TMP_VER ${CMAKE_BINARY_DIR}/version)
SET(SF_TARGET_LIB ${CMAKE_INSTALL_PREFIX}/lib)
SET(SF_TARGET_INC ${CMAKE_INSTALL_PREFIX}/include)
SET(SF_TARGET_ETC ${CMAKE_INSTALL_PREFIX}/etc)
SET(SF_TARGET_BIN ${CMAKE_INSTALL_PREFIX}/bin)
SET(SF_TARGET_DIR ${CMAKE_INSTALL_PREFIX})

# Have the .mod files placed in the INC folder
SET(CMAKE_Fortran_MODULE_DIRECTORY ${SF_TMP_INC})
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_Fortran_MODULE_DIRECTORY})


IF(WITH_BLAS_LAPACK)
  SET(SRC_DIR_BLAS ${SF_SRC}/blas)
  ADD_SUBDIRECTORY(${SRC_DIR_BLAS})

  SET(SRC_DIR_LAPACK ${SF_SRC}/lapack)
  ADD_SUBDIRECTORY(${SRC_DIR_LAPACK})
ENDIF()


SET(SRC_DIR_ARPACK ${SF_SRC}/arpack)
ADD_SUBDIRECTORY(${SRC_DIR_ARPACK})


SET(SRC_DIR_FFTPACK ${SF_SRC}/fftpack)
ADD_SUBDIRECTORY(${SRC_DIR_FFTPACK})


SET(SRC_DIR_MINPACK ${SF_SRC}/minpack)
ADD_SUBDIRECTORY(${SRC_DIR_MINPACK})


SET(SRC_DIR_QUADPACK ${SF_SRC}/quadpack)
ADD_SUBDIRECTORY(${SRC_DIR_QUADPACK})


SET(SRC_DIR_SF_INIT ${SF_SRC}/SF_INIT)
ADD_SUBDIRECTORY(${SRC_DIR_SF_INIT})


SET(SRC_DIR_SF_MPI ${SF_SRC}/SF_MPI)
ADD_SUBDIRECTORY(${SRC_DIR_SF_MPI})


SET(SRC_DIR_SF_IOTOOLS ${SF_SRC}/SF_IOTOOLS)
ADD_SUBDIRECTORY(${SRC_DIR_SF_IOTOOLS})


SET(SRC_DIR_SF_DERIVATE ${SF_SRC}/SF_DERIVATE)
ADD_SUBDIRECTORY(${SRC_DIR_SF_DERIVATE})


SET(SRC_DIR_SF_INTEGRATE ${SF_SRC}/SF_INTEGRATE)
ADD_SUBDIRECTORY(${SRC_DIR_SF_INTEGRATE})


SET(SRC_DIR_SF_INTERPOLATE ${SF_SRC}/SF_INTERPOLATE)
ADD_SUBDIRECTORY(${SRC_DIR_SF_INTERPOLATE})


SET(SRC_DIR_SF_LINALG ${SF_SRC}/SF_LINALG)
ADD_SUBDIRECTORY(${SRC_DIR_SF_LINALG})


SET(SRC_DIR_SF_FFT ${SF_SRC}/SF_FFT)
ADD_SUBDIRECTORY(${SRC_DIR_SF_FFT})


SET(SRC_DIR_SF_OPTIMIZE ${SF_SRC}/SF_OPTIMIZE)
ADD_SUBDIRECTORY(${SRC_DIR_SF_OPTIMIZE})


SET(SRC_DIR_SF_SPECIAL ${SF_SRC}/SF_SPECIAL)
ADD_SUBDIRECTORY(${SRC_DIR_SF_SPECIAL})


SET(SRC_DIR_SF_RANDOM ${SF_SRC}/SF_RANDOM)
ADD_SUBDIRECTORY(${SRC_DIR_SF_RANDOM})


SET(SRC_DIR_SF_SP_LINALG ${SF_SRC}/SF_SP_LINALG)
ADD_SUBDIRECTORY(${SRC_DIR_SF_SP_LINALG})


SET(SRC_DIR_SF_PARSE_INPUT ${SF_SRC}/SF_PARSE_INPUT)
ADD_SUBDIRECTORY(${SRC_DIR_SF_PARSE_INPUT})

ADD_SUBDIRECTORY(${SF_SRC})

# Build the scifor library and have it placed in the target lib folder
IF(WITH_BLAS_LAPACK)
  ADD_LIBRARY(scifor 
    $<TARGET_OBJECTS:BLASLIB>
    $<TARGET_OBJECTS:LAPACKLIB>
    $<TARGET_OBJECTS:ARPACKLIB>
    $<TARGET_OBJECTS:FFTPACKLIB>
    $<TARGET_OBJECTS:MINPACKLIB>
    $<TARGET_OBJECTS:QUADPACKLIB>
    $<TARGET_OBJECTS:SF_INITLIB>
    $<TARGET_OBJECTS:SF_MPILIB>
    $<TARGET_OBJECTS:SF_IOTOOLSLIB>
    $<TARGET_OBJECTS:SF_LINALGLIB>
    $<TARGET_OBJECTS:SF_DERIVATELIB>
    $<TARGET_OBJECTS:SF_OPTIMIZELIB>
    $<TARGET_OBJECTS:SF_INTEGRATELIB>
    $<TARGET_OBJECTS:SF_SPECIALLIB>
    $<TARGET_OBJECTS:SF_FFTLIB>
    $<TARGET_OBJECTS:SF_INTERPOLATELIB>
    $<TARGET_OBJECTS:SF_RANDOMLIB>
    $<TARGET_OBJECTS:SF_SP_LINALGLIB>
    $<TARGET_OBJECTS:SF_PARSE_INPUTLIB>
    $<TARGET_OBJECTS:SCIFORLIB>
    )
ELSE()
  ADD_LIBRARY(scifor 
    $<TARGET_OBJECTS:ARPACKLIB>
    $<TARGET_OBJECTS:FFTPACKLIB>
    $<TARGET_OBJECTS:MINPACKLIB>
    $<TARGET_OBJECTS:QUADPACKLIB>
    $<TARGET_OBJECTS:SF_INITLIB>
    $<TARGET_OBJECTS:SF_MPILIB>
    $<TARGET_OBJECTS:SF_IOTOOLSLIB>
    $<TARGET_OBJECTS:SF_LINALGLIB>
    $<TARGET_OBJECTS:SF_DERIVATELIB>
    $<TARGET_OBJECTS:SF_OPTIMIZELIB>
    $<TARGET_OBJECTS:SF_INTEGRATELIB>
    $<TARGET_OBJECTS:SF_SPECIALLIB>
    $<TARGET_OBJECTS:SF_FFTLIB>
    $<TARGET_OBJECTS:SF_INTERPOLATELIB>
    $<TARGET_OBJECTS:SF_RANDOMLIB>
    $<TARGET_OBJECTS:SF_SP_LINALGLIB>
    $<TARGET_OBJECTS:SF_PARSE_INPUTLIB>
    $<TARGET_OBJECTS:SCIFORLIB>
    )
ENDIF()



SET(SF_VERSION_FILE ${SF_SRC}/${PROJECT_NAME}_version_current_sha1.inc)
ADD_CUSTOM_COMMAND(OUTPUT ${SF_VERSION_FILE}
  COMMAND ${CMAKE_COMMAND} -E echo "character(len=41),parameter,public :: ${PROJECT_NAME}_version_sha1 = \"${GIT_SHA1}\"" > ${SF_VERSION_FILE}
  COMMENT "${Yellow}Creating ${SF_VERSION_FILE}: ${ColourReset}"
  VERBATIM)
ADD_CUSTOM_TARGET(SF_VERSION_SHA1 DEPENDS ${SF_VERSION_FILE})



ADD_DEPENDENCIES(SCIFORLIB
  SF_VERSION_SHA1
  SF_INITLIB
  SF_MPILIB
  SF_IOTOOLSLIB
  SF_FFTLIB
  SF_DERIVATELIB
  SF_LINALGLIB
  SF_OPTIMIZELIB
  SF_SPECIALLIB 
  SF_INTEGRATELIB
  SF_INTERPOLATELIB
  SF_RANDOMLIB
  SF_SP_LINALGLIB
  SF_PARSE_INPUTLIB
  )




SET(PERMISSION_777
  OWNER_WRITE OWNER_READ OWNER_EXECUTE
  GROUP_WRITE GROUP_READ GROUP_EXECUTE
  WORLD_WRITE WORLD_READ WORLD_EXECUTE)



#Build the SCIFOR module for environment modules
FILE(REMOVE_RECURSE ${SF_TMP_ETC}/modules/)
#Create default .version file if following conditions are verified:
#1. master branch
#2. no debug
#3. default prefix
IF( (GIT_BRANCH MATCHES "master")
    AND (NOT "${CMAKE_BUILD_TYPE}" MATCHES "DEBUG")
    AND (PREFIX STREQUAL "${DEFAULT_PREFIX}")
    ) 
  SET(TMP_VER_MODULE_FILE ${SF_TMP_ETC}/modules/${VERSION_PATH}/.version)
  CONFIGURE_FILE(${SF_ENV}/version.in ${TMP_VER_MODULE_FILE}  @ONLY)
  MESSAGE(STATUS "${Red}Version file${ColourReset}: ${TMP_VER_MODULE_FILE}")
ENDIF()

#Build the SCIFOR module for environment modules
SET(TMP_MODULE_NAME "${MODULE_NAME}" CACHE PATH "Prefix prepended to install directories")
SET(TMP_ENV_MODULE_FILE ${SF_TMP_ETC}/modules/${TMP_MODULE_NAME})
CONFIGURE_FILE(${SF_ENV}/module.in ${TMP_ENV_MODULE_FILE} @ONLY)
MESSAGE(STATUS "${Red}Module file${ColourReset}: ${MODULE_NAME}")

#Build the PKG-CONFIG file
INCLUDE(${CMAKE_MODULE_PATH}/BuildPkgConfigFile.cmake)
SET(TMP_PKCONFIG_FILE ${SF_TMP_ETC}/${PROJECT_NAME}.pc)
BUILD_PKCONFIG(${TMP_PKCONFIG_FILE})


#Build the user CONFIG scripts (sourced in user shell config file, i.e. .bashrc)
SET(USER_CONFIG_SCRIPT scifor_config_user.sh)
SET(TMP_CONFIGVARS_USER_FILE ${SF_TMP_ETC}/${USER_CONFIG_SCRIPT})
CONFIGURE_FILE(${SF_ETC}/scifor_config_user.sh.in ${TMP_CONFIGVARS_USER_FILE} @ONLY)

SET(GLOBAL_CONFIG_SCRIPT scifor_config_global.sh)
SET(TMP_CONFIGVARS_GLOBAL_FILE ${SF_TMP_ETC}/${GLOBAL_CONFIG_SCRIPT})
CONFIGURE_FILE(${SF_ETC}/scifor_config_global.sh.in ${TMP_CONFIGVARS_GLOBAL_FILE} @ONLY)



FILE(WRITE  ${SF_TMP_VER}  "${VERSION}\n")

INSTALL(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION ${SF_TARGET_INC})

INSTALL(CODE "execute_process(COMMAND \"${CMAKE_COMMAND}\" -E remove ${SF_VERSION_FILE})")

INSTALL(TARGETS scifor DESTINATION ${SF_TARGET_LIB})

INSTALL(DIRECTORY ${SF_TMP_ETC}/ DESTINATION ${SF_TARGET_ETC})

INSTALL(FILES ${SF_BIN}/scifor_completion.sh DESTINATION ${SF_TARGET_BIN}/
  PERMISSIONS ${PERMISSION_777} SETUID)

INSTALL(FILES ${TMP_CONFIGVARS_USER_FILE} DESTINATION ${SF_TARGET_BIN}/
  PERMISSIONS ${PERMISSION_777} SETUID)

INSTALL(FILES ${TMP_CONFIGVARS_GLOBAL_FILE} DESTINATION ${SF_TARGET_BIN}/
  PERMISSIONS ${PERMISSION_777} SETUID)

INSTALL(FILES ${SF_TMP_VER} DESTINATION ${SF_TARGET_DIR} 
  PERMISSIONS ${PERMISSION_777} SETUID)

INSTALL(FILES ${SF_TARGET_ETC}/${PROJECT_NAME}.pc DESTINATION $ENV{HOME}/.pkgconfig.d/
  PERMISSIONS ${PERMISSION_777} SETUID)

INSTALL(DIRECTORY ${SF_TARGET_ETC}/modules/ DESTINATION $ENV{HOME}/.modules.d)


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
   $ source ${SF_TARGET_BIN}/${USER_CONFIG_FILE}

${Yellow}Method 3: use pkg-config with the provided ${PROJECT_NAME}.pc${ColourReset}:
   $ export PKG_CONFIG_PATH=${SF_TARGET_ETC}/:$PKG_CONFIG_PATH
   $ pkg-config --cflags --libs ${PROJECT_NAME}

${Yellow}Method ADMIN: Add this line to the system shell configuration file, e.g. /etc/bash.bashrc${ColourReset}
   $ source ${SF_TARGET_BIN}/${GLOBAL_CONFIG_FILE}
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


