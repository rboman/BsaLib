cmake_minimum_required( VERSION 3.15 )

set( project BsaLib )

# list( APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
# include(defaultFlags)

# get version from FILE
file(STRINGS "${CMAKE_SOURCE_DIR}/VERSION" PROJECT_VERSION)
string(REPLACE "." ";" VERSION_LIST ${PROJECT_VERSION})
list(GET VERSION_LIST 0 VERSION_MAJOR)
list(GET VERSION_LIST 1 VERSION_MINOR)
list(GET VERSION_LIST 2 VERSION_PATCH)
set( version ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH} )
message( STATUS "Version ${version}" )
unset(VERSION_LIST)


project( ${project} 
   LANGUAGES Fortran 
   VERSION ${version} 
   DESCRIPTION "Library for Bispectral analysis of systems under non-Gaussian random actions"
)
include(CheckFortranSourceCompiles)
include(CheckFortranSourceRuns)


# # DEBUG: check default configuration flags
# foreach(config ${CMAKE_CONFIGURATION_TYPES})
#    string( TOUPPER ${config} config )
#    message( STATUS "Default flags (${config}):  ${CMAKE_Fortran_FLAGS_${config}}" )
# endforeach()



# Project configurable options
option(enable-sequential-code  "Disable any code parallelisation"                   OFF)
option(enable-single           "Enable single precision real kind"                  OFF)
option(enable-heap-arrays      "Allocatable arrays are put on the heap"             ON)
option(enable-shared-runtime   "Link to shared Fortran runtime libraries"           OFF)
option(enable-openmp           "Enable OpenMP parallelisation"                      ON)
option(enable-mkl              "Enable Intel MKL libraries (Intel compilers only)"  ON)
option(enable-sym-ev-routine   "Enable usage of 'dsyev' instead of 'dgesvd'"        ON)
option(enable-par-compilation  "Enable parallelisation of sources compilation"      ON)
option(enable-gpu-code         "Enable GPU kernel code with OpenCL (defaulr) or CUDA backends"  OFF)
option(enable-cuda             "IF GPU code enabled, use CUDA instead of OpenCL (default)"      OFF)


# manage option dependencies
if (enable-sequential-code)
   message( STATUS "Sequential code execution is enabled. Disabling parallelisation.")
   set( enable-openmp   OFF )
   set( enable-mkl      OFF )
   set( enable-gpu-code OFF )
   message( STATUS "OMP: ${enable-openmp};  MKL: ${enable-mkl};  GPU: ${enable-gpu-code}")
endif()


if (WIN32)
   set( CMAKE_DEBUG_POSTFIX "d" )
endif()


# check if intel compiler or not
if (CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
   message( STATUS "NOTE: Using Intel compiler")
   set(intel-compiler ON)
else()
   message( STATUS "NOTE: Compiler ID:  ${CMAKE_Fortran_COMPILER_ID}")
   set(intel-compiler OFF)
endif()


# Flags common to all configurations
if (WIN32)
   add_compile_options(
      $<$<COMPILE_LANGUAGE:Fortran>:-warn:all>
      $<$<COMPILE_LANGUAGE:Fortran>:-stand:f18>
   )
else()
   add_compile_options(
      $<$<COMPILE_LANGUAGE:Fortran>:-warn all>
      $<$<COMPILE_LANGUAGE:Fortran>:-stand f18>
      )
endif()


# Configuration specific options
if (intel-compiler AND enable-heap-arrays)
   add_compile_options( $<$<COMPILE_LANGUAGE:Fortran>:-heap-arrays0> )
endif()

if (intel-compiler)
   add_compile_options(
      $<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANGUAGE:Fortran>>:-O3>
      $<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANGUAGE:Fortran>>:-Ob1>
   )
else()
   add_compile_options(
      $<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANGUAGE:Fortran>>:-Ofast>
   )
endif()


add_compile_definitions( $<$<CONFIG:DEBUG>:_BSA_DEBUG> )

if (enable-gpu-code)
   add_subdirectory( BsaCL )
endif()

add_subdirectory( src )

# reset default flags, if ever changed (hopefully backed up!)
foreach(config ${CMAKE_CONFIGURATION_TYPES})
   string( TOUPPER ${config} config )
   if (DEFINED CMAKE_Fortran_FLAGS_${config}-default)
      set( CMAKE_Fortran_FLAGS_${config} ${CMAKE_Fortran_FLAGS_${config}-default} )
   endif()
endforeach()

add_subdirectory( app )

