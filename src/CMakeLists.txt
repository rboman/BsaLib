
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
   target_compile_definitions( ${project} PRIVATE "-DBSA_USE_SVD_METHOD__" )
else()
   message( STATUS "NOTE: Using 'dsyev'  method")
endif()


# Set module directory
set(LIB_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/modules)
if(NOT EXISTS "${LIB_MOD_DIR}")
  make_directory("${LIB_MOD_DIR}")
endif()
set_target_properties( ${project}
   PROPERTIES
   Fortran_MODULE_DIRECTORY ${LIB_MOD_DIR}
)

# target_include_directories(${project} 
#    PUBLIC
#    $<BUILD_INTERFACE:${LIB_MOD_DIR}>
#    $<INSTALL_INTERFACE:${CMAKE_INSTALL_MODULEDIR}>
# )


