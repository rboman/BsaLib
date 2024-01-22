# TODO: activate generation of shared library as well
add_library( BsaCL STATIC bsacl.c )

if (enable-ocl-ptx-code-gen)
   message( STATUS "NOTE: Compiling PTX code ON")
   target_compile_definitions( BsaCL PRIVATE BSACL_GENERATE_PTX_BINARY )
endif()

if (enable-const-ptr-def)
   target_compile_definitions( BsaCL PRIVATE BSACL_USE_CONST_EXT_POINTERS )
endif()

target_compile_definitions( BsaCL PRIVATE "BASE_DIRECTORY=\"${CMAKE_CURRENT_SOURCE_DIR}/\"" )

target_link_libraries( BsaCL OpenCL::OpenCL )
target_include_directories( BsaCL BEFORE PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../include" )
if (MSVC AND NOT enable-shared-runtime)
   set_property(TARGET BsaCL PROPERTY
      MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()
