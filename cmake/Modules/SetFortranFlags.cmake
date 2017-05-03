# Determine and set the Fortran compiler flags we want 

# Make sure that the default build type is RELEASE if not specified.
INCLUDE(${CMAKE_MODULE_PATH}/SetCompileFlag.cmake)

# Make sure the build type is uppercase
STRING(TOUPPER "${CMAKE_BUILD_TYPE}" BT)

IF(BT STREQUAL "RELEASE")
    SET(CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are DEBUG, RELEASE, or TESTING."
      FORCE)
ELSEIF(BT STREQUAL "DEBUG")
    SET (CMAKE_BUILD_TYPE DEBUG CACHE STRING
      "Choose the type of build, options are DEBUG, RELEASE, or TESTING."
      FORCE)
ELSEIF(BT STREQUAL "TESTING")
    SET (CMAKE_BUILD_TYPE TESTING CACHE STRING
      "Choose the type of build, options are DEBUG, RELEASE, or TESTING."
      FORCE)
ELSEIF(NOT BT)
    SET(CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are DEBUG, RELEASE, or TESTING."
      FORCE)
    MESSAGE(STATUS "CMAKE_BUILD_TYPE not given, defaulting to RELEASE")
ELSE()
    MESSAGE(FATAL_ERROR "CMAKE_BUILD_TYPE not valid, choices are DEBUG, RELEASE, or TESTING")
ENDIF(BT STREQUAL "RELEASE")


# If the compiler flags have already been set, return now
IF(CMAKE_Fortran_FLAGS_RELEASE AND CMAKE_Fortran_FLAGS_TESTING AND CMAKE_Fortran_FLAGS_DEBUG)
    RETURN ()
ENDIF(CMAKE_Fortran_FLAGS_RELEASE AND CMAKE_Fortran_FLAGS_TESTING AND CMAKE_Fortran_FLAGS_DEBUG)


#####################
### GENERAL FLAGS ###
#####################
# Set preprocessing directive:
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}" Fortran REQUIRED 
  "-cpp"        # GNU and possibly all compilers
  "-fpp"       # Intel
  "/fpp"       # Intel Windows
  "/cpp"       # GNU Windows
  )

# There is some bug where -march=native doesn't work on Mac
IF(APPLE)
  SET(GNUNATIVE "-mtune=native")
ELSE()
  SET(GNUNATIVE "-march=native")
ENDIF()

IF(${CMAKE_Fortran_COMPILER_ID} MATCHES Intel)
  SET(CMAKE_Fortran_MODDIR_FLAG "-module ")
ELSEIF(${CMAKE_Fortran_COMPILER_ID} MATCHES GNU)
  SET(CMAKE_Fortran_MODDIR_FLAG "-J")
ELSEIF(${CMAKE_Fortran_COMPILER_ID} MATCHES G95)
  SET(CMAKE_Fortran_MODDIR_FLAG "-fmod=")
ELSEIF(${CMAKE_Fortran_COMPILER_ID} MATCHES PGI)
  SET(CMAKE_Fortran_MODDIR_FLAG "-module ")
ENDIF()

 
# Optimize for the host's architecture
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}" Fortran
  ${GNUNATIVE}    # GNU
  "-xHost"        # Intel
  "/QxHost"       # Intel Windows
  "-ta=host"      # Portland Group
  )

# # Don't add underscores in symbols for C-compatability
# SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}" Fortran "-fno-underscoring")




###################
### DEBUG FLAGS ###
###################
# NOTE: debugging symbols (-g or /debug:full) are already on by default

# Disable optimizations
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}" Fortran REQUIRED
  "-O0" # All compilers not on Windows
  "/Od" # Intel Windows
  )

# Turn on all warnings 
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}" Fortran
  "-warn all" # Intel
  "/warn:all" # Intel Windows
  "-Wall"     # GNU
  # Portland Group (on by default)
  )

# Traceback
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}" Fortran
  "-traceback"   # Intel/Portland Group
  "/traceback"   # Intel Windows
  "-fbacktrace"  # GNU (gfortran)
  "-ftrace=full" # GNU (g95)
  )

# Check array bounds
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}" Fortran
  "-check bounds"  # Intel
  "/check:bounds"  # Intel Windows
  "-fcheck=bounds" # GNU (New style)
  "-fbounds-check" # GNU (Old style)
  "-Mbounds"       # Portland Group
  )



	      
#####################
### TESTING FLAGS ###
#####################
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_TESTING "${CMAKE_Fortran_FLAGS_TESTING}" Fortran  REQUIRED
  "-O1" # All compilers not on Windows
  "/O1" # Intel Windows
  )


SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_TESTING "${CMAKE_Fortran_FLAGS_TESTING}" Fortran
  "-static"			# GNU
  "-static-intel"		# Intel
  )





#####################
### RELEASE FLAGS ###
#####################
# NOTE: agressive optimizations (-O3) are already turned on by default

# Unroll loops
SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}" Fortran
  "-funroll-loops" # GNU
  "-unroll"        # Intel
  "/unroll"        # Intel Windows
  "-Munroll"       # Portland Group
  )


# # Inline functions
# SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}" Fortran
#   "-inline"            # Intel
#   "/Qinline"           # Intel Windows
#   "-finline-functions" # GNU
#   "-Minline"           # Portland Group
#   )

# # Interprocedural (link-time) optimizations
# SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}" Fortran
#   "-ipo"     # Intel
#   "/Qipo"    # Intel Windows
#   "-flto"    # GNU
#   "-Mipa"    # Portland Group
#   )

# # Single-file optimizations
# SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}" Fortran
#   "-ip"  # Intel
#   "/Qip" # Intel Windows
#   )


# # Vectorize code
# SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}" Fortran
#   "-vec-report0"  # Intel
#   "/Qvec-report0" # Intel Windows
#   "-Mvect"        # Portland Group
#   )







# if (NOT DEFINED DEFAULT_Fortran_FLAGS_SET)
# if(CMAKE_Fortran_COMPILER_ID MATCHES GNU) # this is gfortran
#     set(CMAKE_Fortran_FLAGS         "-cpp")
#     set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -fbacktrace")
#     set(CMAKE_Fortran_FLAGS_RELEASE "-w -O3 -ffast-math -funroll-loops -ftree-vectorize")
#     if(ENABLE_BOUNDS_CHECK)
#         set(CMAKE_Fortran_FLAGS
#             "${CMAKE_Fortran_FLAGS} -fbounds-check"
#             )
#     endif()

# elseif(CMAKE_Fortran_COMPILER_ID MATCHES G95)
#     set(CMAKE_Fortran_FLAGS         "-Wno=155 -fno-second-underscore -fsloppy-char")
#     set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -ftrace=full")
#     set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
#     if(ENABLE_BOUNDS_CHECK)
#         set(CMAKE_Fortran_FLAGS
#             "${CMAKE_Fortran_FLAGS} -Wall -fbounds-check"
#             )
#     endif()

# elseif(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
#     set(CMAKE_Fortran_FLAGS         "-w -fpp -assume byterecl -traceback")
#     set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
#     set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -xW -ip")
#     if(ENABLE_BOUNDS_CHECK)
#         set(CMAKE_Fortran_FLAGS
#             "${CMAKE_Fortran_FLAGS} -check bounds -fpstkchk -check pointers -check uninit -check output_conversion"
#             )
#     endif()
# elseif(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
#     set(CMAKE_Fortran_FLAGS         "")
#     set(CMAKE_Fortran_FLAGS_DEBUG   "-g -O0 -Mframe")
#     set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -mcmodel=medium -fast -Munroll")
#     if(ENABLE_BOUNDS_CHECK)
#         set(CMAKE_Fortran_FLAGS
#             "${CMAKE_Fortran_FLAGS} "
#             )
#     endif()

# elseif(CMAKE_Fortran_COMPILER_ID MATCHES XL)
#     set(CMAKE_Fortran_FLAGS         "-qzerosize -qextname")
#     set(CMAKE_Fortran_FLAGS_DEBUG   "-g")
#     set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
#     set_source_files_properties(${FREE_FORTRAN_SOURCES}
#         PROPERTIES COMPILE_FLAGS
#         "-qfree"
#         )
#     set_source_files_properties(${FIXED_FORTRAN_SOURCES}
#         PROPERTIES COMPILE_FLAGS
#         "-qfixed"
#         )
# endif()
# save_compiler_flags(Fortran)
# endif ()
