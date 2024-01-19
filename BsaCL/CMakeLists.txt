cmake_minimum_required( VERSION 3.15 )

project( BsaCL LANGUAGES C Fortran )

if (enable-cuda)
   message( FATAL_ERROR "ERROR: CUDA support not yet implemented in CMake build. Use OpenCL instead.")
else()
   find_package( OpenCL REQUIRED )
   if (NOT OpenCL_FOUND)
      message( SEND_ERROR "COUld not find OpenCL library.")
   endif()
endif()

add_subdirectory( src )
add_subdirectory( BsaCL_Fortran )

