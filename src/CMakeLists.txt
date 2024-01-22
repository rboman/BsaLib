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

set(SRC
   ${CMAKE_CURRENT_SOURCE_DIR}/BsaLib.f90
   ${CMAKE_CURRENT_SOURCE_DIR}/CONSTANTS/constants.f90
   ${CMAKE_CURRENT_SOURCE_DIR}/IO/Logging.f90
   ${CMAKE_CURRENT_SOURCE_DIR}/IO/LoggingImpl.f90
   ${CMAKE_CURRENT_SOURCE_DIR}/IO/io.f90
   ${CMAKE_CURRENT_SOURCE_DIR}/bsa/BsaLibImpl.f90
   ${CMAKE_CURRENT_SOURCE_DIR}/bsa/classic/BsaClassicImpl.f90
   ${CMAKE_CURRENT_SOURCE_DIR}/bsa/data/BsaLibData.f90
   ${CMAKE_CURRENT_SOURCE_DIR}/bsa/data/BsaLibDataImpl.f90
   ${CMAKE_CURRENT_SOURCE_DIR}/bsa/functions/functions.f90
   ${CMAKE_CURRENT_SOURCE_DIR}/bsa/functions/functionsImpl.f90
   ${CMAKE_CURRENT_SOURCE_DIR}/bsa/meshing/BsaMesherImpl.f90
   ${CMAKE_CURRENT_SOURCE_DIR}/bsa/meshing/point/MPoint.f90
   ${CMAKE_CURRENT_SOURCE_DIR}/bsa/meshing/policy/MPolicy.f90
   ${CMAKE_CURRENT_SOURCE_DIR}/bsa/meshing/zones/M2DPolygZone.f90
   ${CMAKE_CURRENT_SOURCE_DIR}/bsa/meshing/zones/MRectZone.f90
   ${CMAKE_CURRENT_SOURCE_DIR}/bsa/meshing/zones/MRectZoneImpl.f90
   ${CMAKE_CURRENT_SOURCE_DIR}/bsa/meshing/zones/MTriangZone.f90
   ${CMAKE_CURRENT_SOURCE_DIR}/bsa/meshing/zones/MTriangZoneImpl.f90
   ${CMAKE_CURRENT_SOURCE_DIR}/bsa/meshing/zones/MZone.f90
   ${CMAKE_CURRENT_SOURCE_DIR}/settings/Settings.f90
   ${CMAKE_CURRENT_SOURCE_DIR}/settings/SettingsImpl.f90
   ${CMAKE_CURRENT_SOURCE_DIR}/structure/StructureData.f90
   ${CMAKE_CURRENT_SOURCE_DIR}/structure/StructureImpl.f90
   ${CMAKE_CURRENT_SOURCE_DIR}/timing/Timer.f90
   ${CMAKE_CURRENT_SOURCE_DIR}/utils/utility.f90
   ${CMAKE_CURRENT_SOURCE_DIR}/wind/WindPSDImpl.f90
   ${CMAKE_CURRENT_SOURCE_DIR}/wind/WindSetImpl.f90
   ${CMAKE_CURRENT_SOURCE_DIR}/wind/WindType.f90
)


# TODO: at some point, enable shared library generation
add_library( ${project} STATIC ${SRC} )
set_target_properties( ${project} 
   PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX}
)


# TODO: see if we can manage this better
# NOTE: disabling anoying -O2 default release flag...
set( CMAKE_Fortran_FLAGS_RELEASE-default ${CMAKE_Fortran_FLAGS_RELEASE} PARENT_SCOPE )
set( CMAKE_Fortran_FLAGS_RELEASE "" )


# optional compile flags
if (enable-par-compilation)
   message( STATUS "NOTE: Enabling parallel sources compilation" )
   if (intel-compiler)
      target_compile_options( ${project} PRIVATE "-MP" )
   else()
      message( WARNING "Enable parallel source compilation for non MSVC!" )
   endif()
endif()

if (intel-compiler AND enable-mkl)
   message( STATUS "NOTE: Using Intel MKL libraries." )
   target_compile_options( ${project} PRIVATE "-Qmkl:parallel" )
else()
   find_package( LAPACK REQUIRED )
   if (NOT LAPACK_FOUND )
      message( SEND_ERROR "Could not find LAPACK." )
   endif()
   target_link_libraries( ${project} ${LAPACK_LIBRARIES} )
   target_link_options( ${project} AFTER PUBLIC ${LAPACK_LINKER_FLAGS} )
endif()

if (enable-openmp)
   if (intel-compiler)
      target_compile_options( ${project} PRIVATE "-Qopenmp" )
   else()
      target_compile_options( ${project} PRIVATE "-fopenmp" )
   endif()
endif()

if (NOT enable-shared-runtime)
   if (intel-compiler)
      target_compile_options( ${project} PUBLIC "-libs:static")
   else()
      message( WARNING "Check how to link to static runtime libraries for non MSVC!" )
   endif()
else()
   message( STATUS "NOTE: Using SHARED runtime libraries!")
endif()

if (NOT enable-sym-ev-routine)
   message( STATUS "NOTE: Using 'dgesvd' method")
   target_compile_definitions( ${project} PRIVATE "-D_BSA_USE_SVD_METHOD" )
else()
   message( STATUS "NOTE: Using 'dsyev'  method")
endif()

if (enable-gpu-code)
   target_compile_definitions( ${project} PRIVATE _BSA_USE_GPU )
   if (enable-cuda)
      target_compile_definitions( ${project} PRIVATE _BSA_USE_CUDA )
   endif()
   target_link_libraries( ${project} BsaCL_Fortran )
endif()


# Set module directory
set(LIB_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/modules)
if(NOT EXISTS "${LIB_MOD_DIR}")
  make_directory("${LIB_MOD_DIR}")
endif()
set_target_properties( ${project} PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_MOD_DIR} )

# target_include_directories(${project} 
#    PUBLIC
#    $<BUILD_INTERFACE:${LIB_MOD_DIR}>
#    $<INSTALL_INTERFACE:${CMAKE_INSTALL_MODULEDIR}>
# )

