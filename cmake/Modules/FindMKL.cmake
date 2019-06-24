# Find a Intel® Math Kernel Library (Intel® MKL) installation and provide
# all necessary variables and macros to compile software for it.
#
# MKLROOT is required in your system
#
# we use mkl_link_tool to get the library needed depending on variables
#
# The following are set after the configuration is done:
#  MKL_FOUND        -  system has MKL
#  MKL_ROOT_DIR     -  path to the MKL base directory
#  MKL_INCLUDE      -  the MKL include directory
#  MKL_LIBRARIES    -  MKL libraries
#
#
#
# Sample usage:
#    find_package(MKL REQUIRED)
#    if (MKL_FOUND)
#        include_directories(${MKL_INCLUDE})
#        # and for each of your dependent executable/library targets:
#        target_link_libraries(<YourTarget> ${MKL_LIBRARIES})
#    endif()
#
#-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -lm -ldl
#-lmkl_gf_lp64 -lmkl_gf_thread -lmkl_core -lpthread -lm -ldl
#
# AUTHOR
# Adriano Amaricci (adriano.amaricci.AT.sissa.it)
# essentially based on previous work of. Simplified version for SciFortran
# Joan MASSICH (joan.massich-vall.AT.inria.fr).
# Alexandre GRAMFORT (Alexandre.Gramfort.AT.inria.fr)
# Théodore PAPADOPOULO (papadop.AT.inria.fr)



SET(CMAKE_FIND_DEBUG_MODE 0)

# Find MKL ROOT
FIND_PATH(MKL_ROOT_DIR NAMES include/mkl.h include/mkl.fi  PATHS $ENV{MKLROOT})

# Convert symlinks to real paths

GET_FILENAME_COMPONENT(MKL_ROOT_DIR ${MKL_ROOT_DIR} REALPATH)

IF (NOT MKL_ROOT_DIR)
  IF (MKL_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "Could not find MKL: please set environment variable {MKLROOT}")
  ELSE()
    UNSET(MKL_ROOT_DIR CACHE)
  ENDIF()
  RETURN()
ENDIF()
  

SET(MKL_INCLUDE_DIR ${MKL_ROOT_DIR}/include)
  
# set arguments to call the MKL provided tool for linking
SET(MKL_LINK_TOOL ${MKL_ROOT_DIR}/tools/mkl_link_tool)

IF(EXISTS "${MKL_LINK_TOOL}")
  IF(APPLE)			#do something specific for Apple
    IF(${CMAKE_Fortran_COMPILER_ID} MATCHES GNU)
      set(MKL_LINK_TOOL_LIBS ${MKL_LINK_TOOL} -check_mkl_presence -libs -l static -p no)
      set(MKL_LINK_TOOL_INCS ${MKL_LINK_TOOL} -check_mkl_presence -opts )
    ELSE()
      set(MKL_LINK_TOOL_LIBS ${MKL_LINK_TOOL} -check_mkl_presence -c intel_f -libs -l static -p no)
      set(MKL_LINK_TOOL_INCS ${MKL_LINK_TOOL} -check_mkl_presence -c intel_f -opts)
    ENDIF()    
  ELSE()			#Linux only system    
    IF(${CMAKE_Fortran_COMPILER_ID} MATCHES GNU)
      set(MKL_LINK_TOOL_LIBS ${MKL_LINK_TOOL} -check_mkl_presence -c gnu_f -libs -l static -p no)
      set(MKL_LINK_TOOL_INCS ${MKL_LINK_TOOL} -check_mkl_presence -c gnu_f -opts)
    ELSE()
      set(MKL_LINK_TOOL_LIBS ${MKL_LINK_TOOL} -check_mkl_presence -c intel_f -libs -l static -p no)
      set(MKL_LINK_TOOL_INCS ${MKL_LINK_TOOL} -check_mkl_presence -c intel_f -opts)
    ENDIF()    
  ENDIF()
  
  EXECUTE_PROCESS(COMMAND  ${MKL_LINK_TOOL_LIBS}
    OUTPUT_VARIABLE MKL_LIBRARIES
    RESULT_VARIABLE COMMAND_WORKED
    TIMEOUT 2 ERROR_QUIET)
  if ( MKL_FIND_REQUIRED AND (NOT ${COMMAND_WORKED} EQUAL 0))
    MESSAGE(FATAL_ERROR "Cannot find MKL libraries. The mkl_link_tool command executed was:\n ${MKL_LINK_TOOL_LIBS}.")
  endif()
  
  EXECUTE_PROCESS(COMMAND ${MKL_LINK_TOOL_INCS}
    OUTPUT_VARIABLE MKL_INCLUDE
    RESULT_VARIABLE COMMAND_WORKED
    TIMEOUT 2 ERROR_QUIET)
  if ( MKL_FIND_REQUIRED AND (NOT ${COMMAND_WORKED} EQUAL 0))
    MESSAGE(FATAL_ERROR "Cannot find MKL libraries. The mkl_link_tool command executed was:\n ${MKL_LINK_TOOL_INCS}.")
  endif()

  SET(MKL_LIBRARIES_ ${MKL_LIBRARIES})
  SET(MKL_INCLUDE_ ${MKL_INCLUDE})
  STRING(STRIP ${MKL_LIBRARIES_} MKL_LIBRARIES)
  STRING(STRIP ${MKL_INCLUDE_} MKL_INCLUDE)
  

ELSE()
  
  #ON OSX THERE IS NO SUPPORT FOR GNU_F SO THERE IS NO DISTINCTION GNU/INTEL
  IF(APPLE)
    SET(LP64_LIB       "libmkl_intel_lp64.a")
    SET(SEQUENTIAL_LIB "libmkl_sequential.a")
    SET(THREAD_LIB     "libmkl_intel_thread.a")
    SET(CORE_LIB       "libmkl_core.a")  
  ELSE()
    IF(${CMAKE_Fortran_COMPILER_ID} MATCHES GNU)
      SET(LP64_LIB       "libmkl_gf_lp64.a")
      SET(SEQUENTIAL_LIB "libmkl_sequential.a")
      SET(THREAD_LIB     "libmkl_gnu_thread.a")
      SET(CORE_LIB       "libmkl_core.a")  
    ELSEIF(${CMAKE_Fortran_COMPILER_ID} MATCHES INTEL)
      SET(LP64_LIB       "libmkl_intel_lp64.a")
      SET(SEQUENTIAL_LIB "libmkl_sequential.a")
      SET(THREAD_LIB     "libmkl_intel_thread.a")
      SET(CORE_LIB       "libmkl_core.a")
    ELSE()
      MESSAGE(FATAL_ERROR "${CMAKE_Fortran_COMPILER_ID} not supported in MKL")
    ENDIF()
  ENDIF()

  FIND_LIBRARY(MKL_LP64_LIBRARY
    NAMES ${LP64_LIB}
    PATHS ${MKL_ROOT_DIR}/lib
    ${MKL_ROOT_DIR}/lib/intel64
    $ENV{INTEL}/mkl/lib/intel64
    NO_DEFAULT_PATH)
  
  FIND_LIBRARY(MKL_SEQUENTIAL_LIBRARY
    NAMES ${SEQUENTIAL_LIB}
    PATHS ${MKL_ROOT_DIR}/lib
    ${MKL_ROOT_DIR}/lib/intel64
    $ENV{INTEL}/mkl/lib/intel64
    NO_DEFAULT_PATH)
  
  FIND_LIBRARY(MKL_THREAD_LIBRARY
    NAMES ${THREAD_LIB}
    PATHS ${MKL_ROOT_DIR}/lib
    ${MKL_ROOT_DIR}/lib/intel64
    $ENV{INTEL}/mkl/lib/intel64
    NO_DEFAULT_PATH)
  
  FIND_LIBRARY(MKL_CORE_LIBRARY
    NAMES ${CORE_LIB}
    PATHS ${MKL_ROOT_DIR}/lib
    ${MKL_ROOT_DIR}/lib/intel64
    $ENV{INTEL}/mkl/lib/intel64
    NO_DEFAULT_PATH)
  
  SET(MKL_INCLUDE ${MKL_INCLUDE_DIR})
  SET(MKL_LIBRARIES "${MKL_LP64_LIBRARY} ${MKL_SEQUENTIAL_LIBRARY} ${MKL_THREAD_LIBRARY} ${MKL_CORE_LIBRARY}")
  
ENDIF()




INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(MKL DEFAULT_MSG MKL_ROOT_DIR MKL_LIBRARIES MKL_INCLUDE )

MARK_AS_ADVANCED(MKL_INCLUDE MKL_LIBRARIES MKL_ROOT_DIR)

if (CMAKE_FIND_DEBUG_MODE)
  MESSAGE(STATUS "Found MKL_LIBRARIES:${MKL_LIBRARIES}")
  MESSAGE(STATUS "Found MKL_INCLUDE:${MKL_INCLUDE}")
endif()

