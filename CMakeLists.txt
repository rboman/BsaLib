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


# Project configurable options
option(enable-par-compilation "Enable parallelisation of sources compilation"      ON)
option(enable-single          "Enable single precision real kind"                  OFF)
option(enable-openmp          "Enable OpenMP parallelisation"                      ON)
option(enable-heap-arrays     "Allocatable arrays are put on the heap"             ON)
option(enable-mkl             "Enable Intel MKL libraries (Intel compilers only)"  ON)
option(enable-shared-runtime  "Link to shared Fortran runtime libraries"           OFF)
option(enable-sym-ev-routine  "Enable usage of 'dsyev' instead of 'dgesvd'"        ON)


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
foreach(config ${CMAKE_CONFIGURATION_TYPES})
   string( TOUPPER ${config} config )
   if (intel-compiler)
      add_compile_options(
         $<$<CONFIG:${config}>:-warn:all>
         $<$<CONFIG:${config}>:-stand:f18>
      )
   else()
      add_compile_options(
         $<$<CONFIG:${config}>:-warn all>
         $<$<CONFIG:${config}>:-stand f18>
      )
   endif()
endforeach()


# Configuration specific options
if (intel-compiler AND enable-heap-arrays)
   add_compile_options( "-heap-arrays0" )
endif()

if (intel-compiler)
   add_compile_options(
      $<$<CONFIG:RELEASE>:-O3>
      $<$<CONFIG:RELEASE>:-Ob1>
   )
else()
   add_compile_options(
      $<$<CONFIG:RELEASE>:-Ofast>
   )
endif()


add_subdirectory( src )
add_subdirectory( app )

