enable_language( Fortran )
add_executable( bsa bsa.f90 )

target_link_libraries( bsa BsaLib )

if (intel-compiler)
   target_compile_options( bsa PRIVATE -w )
endif()

if (enable-gpu-code)
   if (enable-cuda)
      set_target_properties(bsa PROPERTIES OUTPUT_NAME "bsa_cuda")
      set_property(TARGET bsa PROPERTY CUDA_SEPARABLE_COMPILATION ON)
      if(APPLE)
         # We need to add the path to the driver (libcuda.dylib) as an rpath,
         # so that the static cuda runtime can find it at runtime.
         set_property(TARGET bsa 
            PROPERTY
            BUILD_RPATH ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})
      endif()
   else()
      set_target_properties(bsa PROPERTIES OUTPUT_NAME "bsa_cl")
   endif()
endif()

set(LIB_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/modules)
if(NOT EXISTS "${LIB_MOD_DIR}")
  make_directory("${LIB_MOD_DIR}")
endif()
set_target_properties( bsa
   PROPERTIES
   Fortran_MODULE_DIRECTORY ${LIB_MOD_DIR}
)
