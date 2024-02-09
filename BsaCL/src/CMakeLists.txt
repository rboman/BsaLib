# TODO: activate generation of shared library as well
add_library( BsaCL STATIC bsacl.c )

if (enable-gpu-double)
   target_compile_definitions( BsaCL PRIVATE BSACL_USE_DOUBLE_PRECISION )
endif()


if (enable-ocl-ptx-code-gen)
   message( STATUS "NOTE: Compiling PTX code ON")
   target_compile_definitions( BsaCL PRIVATE BSACL_GENERATE_PTX_BINARY )
endif()

if (enable-const-ptr-def)
   target_compile_definitions( BsaCL PRIVATE BSACL_USE_CONST_EXT_POINTERS )
endif()

if (enable-cuda)
   set_target_properties( BsaCL PROPERTIES CUDA_ARCHITECTURES native )
   set_source_files_properties( bsacl.c PROPERTIES LANGUAGE CUDA )
   target_compile_options( BsaCL PRIVATE -diag-suppress 546 )
   target_compile_options( BsaCL PRIVATE
      $<$<CONFIG:Release>: --use_fast_math --extra-device-vectorization>
   )
   target_compile_definitions( BsaCL PRIVATE BSACL_USE_CUDA__ )
   set_target_properties( BsaCL PROPERTIES POSITION_INDEPENDENT_CODE ON )
else()
   target_compile_definitions( BsaCL PRIVATE "BASE_DIRECTORY=\"${CMAKE_CURRENT_SOURCE_DIR}/\"" )
   target_link_libraries( BsaCL OpenCL::OpenCL )
endif()

target_include_directories( BsaCL BEFORE PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../include" )

if (MSVC AND NOT enable-shared-runtime)
   set_property(TARGET BsaCL PROPERTY
      MSVC_RUNTIME_LIBRARY "MultiThreadedDebug")
endif()
