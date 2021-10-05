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
IF("${BUILD_TYPE}" MATCHES "DEBUG")
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

