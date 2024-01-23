#
# Copyright (c), Michele Esposito Marzino
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#      * Redistributions of source code must retain the above copyright notice, this
#        list of conditions and the following disclaimer.
#
#      * Redistributions in binary form must reproduce the above copyright notice,
#        this list of conditions and the following disclaimer in the documentation
#        and/or other materials provided with the distribution.
#
#      * Neither the name of the copyright holder nor the names of its
#        contributors may be used to endorse or promote products derived from
#        this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

cmake_minimum_required( VERSION 3.15 )


###########################################################
# Main project name
###########################################################
file(STRINGS "${CMAKE_SOURCE_DIR}/VERSION" PROJECT_VERSION)
set( project BsaLib )
project( ${project} 
   LANGUAGES Fortran 
   VERSION ${PROJECT_VERSION} 
   DESCRIPTION "Library for Bispectral analysis of systems under non-Gaussian wind action"
)
include(CheckFortranSourceCompiles)
include(CheckFortranSourceRuns)


###########################################################
# Options
###########################################################
option(enable-sequential-code  "Disable any code parallelisation"                   OFF)
option(enable-single           "Enable single floating precision"                   OFF)
option(enable-heap-arrays      "Allocatable arrays are put on the heap"             ON)
option(enable-shared-runtime   "Link to shared Fortran runtime libraries"           OFF)
option(enable-openmp           "Enable OpenMP parallelisation"                      ON)
option(enable-mkl              "Enable Intel MKL libraries (Intel compilers only)"  ON)
option(enable-sym-ev-routine   "Enable usage of 'dsyev' instead of 'dgesvd'"        ON)
option(enable-par-compilation  "Enable parallelisation of sources compilation"      ON)

option(enable-gpu-code         "Enable GPU kernel code with OpenCL (defaulr) or CUDA backends"  OFF)
option(enable-const-ptr-def    "In GPU related code, declare internal pointers `const`"         ON)
option(enable-cuda             "IF GPU code enabled, use CUDA instead of OpenCL (default)"      OFF)


if (WIN32)
   set( CMAKE_DEBUG_POSTFIX "d" )
endif()


###########################################################
# General section
###########################################################
if (enable-sequential-code)
   message( STATUS "Sequential code execution is enabled. Disabling parallelisation.")
   set( enable-openmp   OFF )
   set( enable-mkl      OFF )
   set( enable-gpu-code OFF )
   message( STATUS "OMP: ${enable-openmp};  MKL: ${enable-mkl};  GPU: ${enable-gpu-code}")
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

if (enable-single)
   message( STATUS "NOTE: enabling single precision floating point arithmetic")
   add_compile_definitions( BSA_SINGLE_FLOATING_PRECISION )
endif()


###########################################################
# Targets compilation
###########################################################
if (enable-gpu-code)
   if (enable-cuda AND enable-ocl-ptx-code-gen)
      set( enable-ocl-ptx-code-gen OFF )
   endif()
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

