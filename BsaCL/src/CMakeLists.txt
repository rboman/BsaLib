
# TODO: activate generation of shared library as well
add_library( BsaCL STATIC bsacl.c )
target_link_libraries( BsaCL OpenCL::OpenCL )
target_include_directories( BsaCL BEFORE PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../include" )
if (MSVC)
   set_property(TARGET BsaCL PROPERTY
      MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()
