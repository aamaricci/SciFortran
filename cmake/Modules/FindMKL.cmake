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
# AUTHORS
# Adriano Amaricci (adriano.amaricci.AT.sissa.it)
# Lorenzo Crippa (crippa.AT.sissa.it)
#
# essentially based on previous work of: 
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

EXECUTE_PROCESS(COMMAND  ${MKL_LINK_TOOL}
  OUTPUT_QUIET  
  RESULT_VARIABLE RET)

IF(NOT RET EQUAL "0")
  UNSET(MKL_LINK_TOOL)
ENDIF()



IF(EXISTS "${MKL_LINK_TOOL}")

  #CHECK IF SCALAPACK IS ACTUALLY SUPPORTED
  FIND_LIBRARY(MKL_SCALAPACK_LIBRARY
    NAMES mkl_scalapack_lp64 mkl_scalapack_core mkl_scalapack
    PATHS ${MKL_ROOT_DIR}/lib
    ${MKL_ROOT_DIR}/lib/intel64
    $ENV{INTEL}/mkl/lib/intel64
    NO_DEFAULT_PATH)
  
  FIND_LIBRARY(MKL_BLACS_LIBRARY
    NAMES mkl_blacs mkl_blacs_lp64 mkl_blacs_mpich_lp64 mkl_blacs_intelmpi_lp64 mkl_blacs_openmpi_lp64
    PATHS ${MKL_ROOT_DIR}/lib
    ${MKL_ROOT_DIR}/lib/intel64
    $ENV{INTEL}/mkl/lib/intel64
    NO_DEFAULT_PATH)

  IF(MKL_SCALAPACK_LIBRARY AND MKL_BLACS_LIBRARY)
    MESSAGE(STATUS "Found MKL Scalapack support")
    SET(MKL_SCALAPACK_FOUND TRUE)
  ENDIF()

  #These options can be made available to the user in the end
  SET(MKL_PARALLEL  "--parallel=no")
  #SET(MKL_LINKING   "--linking=static")
  IF(MKL_SCALAPACK_FOUND)
    SET(MKL_SCALAPACK "--cluster_library=scalapack")
  ENDIF()
  
  IF(APPLE)			#do something specific for Apple
    IF(${CMAKE_Fortran_COMPILER_ID} MATCHES GNU)
      set(MKL_LINK_TOOL_LIBS ${MKL_LINK_TOOL} -check_mkl_presence -libs ${MKL_PARALLEL} ${MKL_LINKING} ${MKL_SCALAPACK})
      set(MKL_LINK_TOOL_INCS ${MKL_LINK_TOOL} -check_mkl_presence -opts )
    ELSE()
      set(MKL_LINK_TOOL_LIBS ${MKL_LINK_TOOL} -check_mkl_presence -c intel_f -libs ${MKL_PARALLEL} ${MKL_LINKING} ${MKL_SCALAPACK})
      set(MKL_LINK_TOOL_INCS ${MKL_LINK_TOOL} -check_mkl_presence -c intel_f -opts)
    ENDIF()    
  ELSE()			#Linux only system    
    IF(${CMAKE_Fortran_COMPILER_ID} MATCHES GNU)
      set(MKL_LINK_TOOL_LIBS ${MKL_LINK_TOOL} -check_mkl_presence -c gnu_f -libs ${MKL_PARALLEL} ${MKL_LINKING} ${MKL_SCALAPACK})
      set(MKL_LINK_TOOL_INCS ${MKL_LINK_TOOL} -check_mkl_presence -c gnu_f -opts)
    ELSE()
      set(MKL_LINK_TOOL_LIBS ${MKL_LINK_TOOL} -check_mkl_presence -c intel_f -libs ${MKL_PARALLEL} ${MKL_LINKING} ${MKL_SCALAPACK})
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

  #Take care of the way mkl_link_tool handles the variable $(MKLROOT) in the linking flag. Replace it with actual
  #value of the MKLROOT
  STRING(REPLACE "$(MKLROOT)" ${MKL_ROOT_DIR} MKL_LIBRARIES ${MKL_LIBRARIES}) 
  # STRING(REGEX REPLACE "\n$" "" MKL_LIBRARIES "${MKL_LIBRARIES}")
  
  SET(MKL_LIBRARIES_ ${MKL_LIBRARIES})
  SET(MKL_INCLUDE_ ${MKL_INCLUDE})
  IF(APPLE)
    STRING(STRIP ${MKL_LIBRARIES_} ${MKL_LIBRARIES})
    STRING(STRIP ${MKL_INCLUDE_} ${MKL_INCLUDE})
  ELSE()
    STRING(STRIP $MKL_LIBRARIES_ ${MKL_LIBRARIES})
    STRING(STRIP $MKL_INCLUDE_ ${MKL_INCLUDE})
  ENDIF()

#NO MKL_LINK_TOOL: DO OUR OWN
ELSE()
  
  SET(FC_ID ${CMAKE_Fortran_COMPILER_ID})
  STRING(TOUPPER "${FC_ID}" FC_ID)	

  #ON OSX THERE IS NO SUPPORT FOR GNU_F SO THERE IS NO DISTINCTION GNU/INTEL
  IF(APPLE)
    SET(LP64_LIB       "libmkl_intel_lp64.a")
    SET(SEQUENTIAL_LIB "libmkl_sequential.a")
    SET(THREAD_LIB     "libmkl_intel_thread.a")
    SET(CORE_LIB       "libmkl_core.a")
    SET(PTHREAD_LIB    "libpthread.a")
    SET(MATH_LIB       "libm.a")
    SET(DL_LIB         "libdl.a")
  ELSE()
    IF(${FC_ID} MATCHES GNU)
      SET(LP64_LIB       "libmkl_gf_lp64.a")
      SET(SEQUENTIAL_LIB "libmkl_sequential.a")
      SET(THREAD_LIB     "libmkl_gnu_thread.a")
      SET(CORE_LIB       "libmkl_core.a")
      SET(PTHREAD_LIB    "libpthread.a")
      SET(MATH_LIB       "libm.a")
      SET(DL_LIB         "libdl.a")
    ELSEIF(${FC_ID} MATCHES INTEL)
      SET(LP64_LIB       "libmkl_intel_lp64.a")
      SET(SEQUENTIAL_LIB "libmkl_sequential.a")
      SET(THREAD_LIB     "libmkl_intel_thread.a")
      SET(CORE_LIB       "libmkl_core.a")
      SET(PTHREAD_LIB    "libpthread.a")
      SET(MATH_LIB       "libm.a")
      SET(DL_LIB         "libdl.a")
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
  IF(NOT MKL_LP64_LIBRARY)
    MESSAGE(FATAL_ERROR "Can not find ${LP64_LIB} in MKL/ OS")
  ENDIF()
    
  FIND_LIBRARY(MKL_SEQUENTIAL_LIBRARY
    NAMES ${SEQUENTIAL_LIB}
    PATHS ${MKL_ROOT_DIR}/lib
    ${MKL_ROOT_DIR}/lib/intel64
    $ENV{INTEL}/mkl/lib/intel64
    NO_DEFAULT_PATH)
  IF(NOT MKL_SEQUENTIAL_LIBRARY)
    MESSAGE(FATAL_ERROR "Can not find ${SEQUENTIAL_LIB} in MKL/ OS")
  ENDIF()
  
  FIND_LIBRARY(MKL_THREAD_LIBRARY
    NAMES ${THREAD_LIB}
    PATHS ${MKL_ROOT_DIR}/lib
    ${MKL_ROOT_DIR}/lib/intel64
    $ENV{INTEL}/mkl/lib/intel64
    NO_DEFAULT_PATH)
  IF(NOT MKL_THREAD_LIBRARY)
    MESSAGE(FATAL_ERROR "Can not find ${THREAD_LIB} in MKL/ OS")
  ENDIF()
  
  FIND_LIBRARY(MKL_CORE_LIBRARY
    NAMES ${CORE_LIB}
    PATHS ${MKL_ROOT_DIR}/lib
    ${MKL_ROOT_DIR}/lib/intel64
    $ENV{INTEL}/mkl/lib/intel64
    NO_DEFAULT_PATH)
  IF(NOT MKL_CORE_LIBRARY)
    MESSAGE(FATAL_ERROR "Can not find ${CORE_LIB} in MKL/ OS")
  ENDIF()
  
  FIND_LIBRARY(PTHREAD_LIBRARY
    NAMES ${PTHREAD_LIB}
    NO_DEFAULT_PATH)
  IF(NOT PTHREAD_LIBRARY)
    MESSAGE(FATAL_ERROR "Can not find ${PTHREAD_LIB} in OS")
  ENDIF()
  
  FIND_LIBRARY(MATH_LIBRARY
    NAMES ${MATH_LIB}
    NO_DEFAULT_PATH)
  IF(NOT MATH_LIBRARY)
    MESSAGE(FATAL_ERROR "Can not find ${MATH_LIB} in OS")
  ENDIF()
  
  FIND_LIBRARY(DL_LIBRARY
    NAMES ${DL_LIB}
    NO_DEFAULT_PATH)
  IF(NOT DL_LIBRARY)
    MESSAGE(FATAL_ERROR "Can not find ${DL_LIB} in OS")
  ENDIF()


  # #TRY TO FIND SCALAPACK SUPPORT
  FIND_LIBRARY(MKL_SCALAPACK_LIBRARY
    NAMES mkl_scalapack_lp64 mkl_scalapack_core
    PATHS ${MKL_ROOT_DIR}/lib
    ${MKL_ROOT_DIR}/lib/intel64
    $ENV{INTEL}/mkl/lib/intel64
    NO_DEFAULT_PATH)
  
  FIND_LIBRARY(MKL_BLACS_LIBRARY
    NAMES mkl_blacs_lp64 mkl_blacs_mpich_lp64 mkl_blacs_intelmpi_lp64 mkl_blacs_openmpi_lp64
    PATHS ${MKL_ROOT_DIR}/lib
    ${MKL_ROOT_DIR}/lib/intel64
    $ENV{INTEL}/mkl/lib/intel64
    NO_DEFAULT_PATH)

  IF(MKL_SCALAPACK_LIBRARY AND MKL_BLACS_LIBRARY)
    MESSAGE(STATUS "Found MKL Scalapack support")
    SET(MKL_SCALAPACK_FOUND TRUE)
  ENDIF()

  SET(MKL_INCLUDE ${MKL_INCLUDE_DIR})
  IF(MKL_SCALAPACK_FOUND)
    SET(MKL_LIBRARIES "${MKL_LP64_LIBRARY} ${MKL_SEQUENTIAL_LIBRARY} 
${MKL_THREAD_LIBRARY} ${MKL_CORE_LIBRARY} ${MKL_SCALAPACK_LIBRARY} ${MKL_BLACS_LIBRARY} 
${PTHREAD_LIBRARY} ${MATH_LIBRARY} ${DL_LIBRARY}")
  ELSE()
    SET(MKL_LIBRARIES "${MKL_LP64_LIBRARY} ${MKL_SEQUENTIAL_LIBRARY} 
${MKL_THREAD_LIBRARY} ${MKL_CORE_LIBRARY} ${PTHREAD_LIBRARY} ${MATH_LIBRARY} ${DL_LIBRARY}")
  ENDIF()
  
  
ENDIF()




INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(MKL DEFAULT_MSG MKL_ROOT_DIR MKL_LIBRARIES MKL_INCLUDE )

MARK_AS_ADVANCED(MKL_INCLUDE MKL_LIBRARIES MKL_ROOT_DIR)

if (CMAKE_FIND_DEBUG_MODE)
  MESSAGE(STATUS "Found MKL_LIBRARIES:${MKL_LIBRARIES}")
  MESSAGE(STATUS "Found MKL_INCLUDE:${MKL_INCLUDE}")
endif()

