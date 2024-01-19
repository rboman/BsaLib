
project( BsaCL_Fortran LANGUAGES Fortran )

add_library( BsaCL_Fortran 
   STATIC
   bsacl.f90
   bsacl_c.f90
   bsacl_impl.f90
)

target_link_libraries( BsaCL_Fortran BsaCL )

set(cur_mod_dir ${CMAKE_CURRENT_BINARY_DIR}/modules)
if(NOT EXISTS "${cur_mod_dir}")
  make_directory("${cur_mod_dir}")
endif()
set_target_properties( ${PROJECT_NAME} PROPERTIES Fortran_MODULE_DIRECTORY ${cur_mod_dir} )
target_include_directories( ${PROJECT_NAME} 
   PUBLIC $<BUILD_INTERFACE:${cur_mod_dir}>
)
unset( cur_mod_dir )
