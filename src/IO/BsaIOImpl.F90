submodule(BsaLib_IO) BsaLib_IO_Impl

   use BsaLib_Data, only: bsa_Abort, settings, struct_data, wd
   implicit none

contains

   module subroutine allocKOMsg(name_, istat, emsg)
      character(len = *), intent(in) :: name_, emsg
      integer, intent(in) :: istat

      write(unit_debug_, fmt='(3a)') &
         '[ERROR] variable  "', name_, '"  could not be allocated.'
      write(unit_debug_, fmt='(15x, a, i0, 2a)') &
         'Exit code  ', istat, '. Error message:  ', emsg(1 : len_trim(emsg))
      call bsa_Abort()
   end subroutine



   module subroutine deallocKOMsg(name_, istat, emsg)
      character(len = *), intent(in) :: name_, emsg
      integer, intent(in) :: istat

      write(unit_debug_, fmt='(3a)') &
         '[ERROR] variable  "', name_, '"  could not be de-allocated.'
      write(unit_debug_, fmt='(15x, a, i0, 2a)') &
         'Exit code  ', istat, '. Error message:  ', emsg(1 : len_trim(emsg))

      call bsa_Abort()
   end subroutine



   module subroutine io_printUserData()
      character(len=64) :: fmt
      character(len=64) :: fmt2
      integer(int32) :: i

      write(fmt, '(a)') '("    - ", a, i10)'

      write(unit_debug_, *) '### settings:'
      write(unit_debug_, fmt) 'SUBANALYSIS TYPE   = ',  settings%i_suban_type_
      write(unit_debug_, fmt) 'VERSION            = ',  settings%i_vers_
      write(unit_debug_, fmt) 'CONVENTION USED    = ',  settings%i_def_scaling_
      write(unit_debug_, fmt) 'COMPUTE PSDs       = ',  settings%i_compute_psd_
      write(unit_debug_, fmt) 'COMPUTE BISP       = ',  settings%i_compute_bisp_
      write(unit_debug_, fmt) 'ONLY DIAG ELEMENTS = ',  settings%i_only_diag_
      write(unit_debug_, fmt) 'TESTING MODE       = ',  settings%i_test_mode_
      write(unit_debug_, fmt) 'DUMP MODAL INFO    = ',  settings%i_dump_modal_
      write(unit_debug_, fmt) 'USE BISP      SYM  = ',  settings%i_bisp_sym_
      write(unit_debug_, fmt) 'USE 3D MATRIX SYM  = ',  settings%i_spctr_sym_
      write(unit_debug_, fmt) 'N. OF FREQUENCIES  = ',  settings%nfreqs_
      write(unit_debug_, '("    - ", a, g10.5)') 'DELTA FREQ         = ',  settings%df_
      write(unit_debug_, fmt) 'USE "SVD" DECOMP   = ',  settings%i_use_svd_
      write(unit_debug_, fmt) 'BKG_BASE_RFMT      = ',  settings%bkg_base_rfmnt_
      write(unit_debug_, '("    - ", a, g10.5)') 'BKG_AERA_EXT       = ',  settings%bkg_area_ext_
      write(unit_debug_, '("    - ", a, g10.5)') 'GEN_PEAK_AREA_EXT  = ',  settings%peak_area_ext_
      write(unit_debug_, '("    - ", a, g10.5)') 'MAX AREA EXTENSION = ',  settings%max_area_ext_
      write(unit_debug_, fmt) 'DO FULL COVERAGE   = ',  settings%i_full_coverage_


      write(unit_debug_, *) '### structure:'
      write(unit_debug_, fmt) 'NLIBS         = ', struct_data%nlibs_
      write(unit_debug_, fmt) 'NNODES        = ', struct_data%nn_
      write(unit_debug_, fmt) 'NDOFS         = ', struct_data%ndofs_
      write(unit_debug_, fmt) 'NLIBS LOADED  = ', struct_data%nlibs_load_
      write(unit_debug_, '(*(i5))') struct_data%libs_load_
      write(unit_debug_, fmt) 'NODES LOADED  = ', struct_data%nn_load_
      write(unit_debug_, '(*(10i5))') struct_data%n_load_
      write(fmt2, '(a)') '( "    - ", a)'             ! BUG: only actual loaded nodes are saved !!
      write(unit_debug_, fmt2) 'NODAL COORDS  = '
      write(fmt2, '(a)') '( "n.", i5, ":", 3(2x, g10.4) )'
      do i = 1, struct_data%nn_
         write(unit_debug_, fmt2) i, struct_data%coords_(:, i)
      enddo
      write(unit_debug_, fmt) 'N. MODES      = ', struct_data%modal_%nm_
      write(unit_debug_, fmt) 'N. MODES EFF  = ', struct_data%modal_%nm_eff_
      write(fmt2, '(a)') '( "    - ", a, /, *(g10.4) )'
      write(unit_debug_, fmt2) 'MODES KEPT    = ', struct_data%modal_%modes_
      write(unit_debug_, fmt2) 'NAT. FREQS    = ', struct_data%modal_%nat_freqs_
      write(fmt2, '(a, i5, a)') ' ( "    - ", a, /,  *(',  struct_data%modal_%nm_, '(2x, g12.6), /) )'
      ! WARNING: creates a temporary array..
      write(unit_debug_, fmt2) 'MOD. MAT    = ', (struct_data%modal_%phi_(i, :), i = 1, struct_data%ndofs_)
      write(unit_debug_, fmt2) 'M*          = ', struct_data%modal_%Mm_
      write(unit_debug_, fmt2) 'K*          = ', struct_data%modal_%Km_
      ! WARNING: creates a temporary array..
      write(unit_debug_, fmt2) 'C*          = ', (struct_data%modal_%Cm_(i, :), i = 1, struct_data%modal_%nm_)
      write(unit_debug_, fmt2) 'XSI         = ', struct_data%modal_%xsi_



      write(unit_debug_, fmt) 'WIND ZONES  = ', wd%nz_
      write(unit_debug_, fmt) 'PSD TYPE    = ', wd%i_psd_type_
      write(unit_debug_, fmt) 'EQ. NOD. VEL= ', wd%i_eq_nod_wind_speed_
      write(unit_debug_, fmt) 'WIND PROF   = ', wd%i_wind_prof_
      write(unit_debug_, fmt) 'WIND iVERT  = ', wd%i_vert_
      write(unit_debug_, fmt) 'TURB COMP   = ', wd%i_ntc_
      write(unit_debug_, '(*(i5))') wd%tc_
      write(unit_debug_, fmt) 'WIND DIRS   = ', wd%i_ndirs_
      write(unit_debug_, '(*(i5))') wd%dirs_
      write(unit_debug_, fmt) 'WIND SPEEDS = ', size(wd%u_node_)
      write(unit_debug_, '(*(10("  ", g10.4), /))') wd%u_node_
      write(unit_debug_, fmt) 'NODE W. ZONE= '
      write(unit_debug_, '(*(10("  ", g10.4), /))') wd%wz_node_
      write(unit_debug_, fmt) 'NODE W. ALT = '
      write(unit_debug_, '(*(10("  ", g10.4), /))') wd%wAlt_node_

      ! WARNING: creates temporary because of non-contiguous memory..
      write(unit_debug_, fmt) 'NODAL COHER = ', size(wd%nod_corr_, 1)
      write(unit_debug_, '(*(10("  ", g10.4), /))') wd%nod_corr_(:, 1)

      write(unit_debug_, fmt) 'W. F. C.    = '
      write(fmt2, '( 2(a, i5), a)') ' ( ', size(wd%wfc_, 2), '(', size(wd%wfc_, 1), '(2x, g10.4), /) )'
      do i = 1, size(wd%wfc_, 3)
         write(unit_debug_, fmt2) wd%wfc_(:, :, i)
      enddo

      write(fmt2, '(a, i2, a)') '( ', wd%nz_, '( 3g10.4, 1x, / ) )'
      write(unit_debug_, fmt) 'WZ. STD. DEV= '
      write(unit_debug_, fmt2) wd%sigmaUVW_wz_
      write(fmt2, '(a, i2, a)') '( ', wd%nz_, '(g10.4, 1x), / )'
      write(unit_debug_, fmt) 'WZ. ZREF    = '
      write(unit_debug_, fmt2) wd%Zref_wz_
      write(unit_debug_, fmt) 'WZ. UB ZREF = '
      write(unit_debug_, fmt2) wd%u_mean_ref_wz_
      write(unit_debug_, fmt) 'WZ. INC ANG = '
      write(unit_debug_, fmt2) wd%incAng_wz_
      write(fmt2, '(a, i2, a)') '( ', wd%nz_ + 1, '(g10.4, 1x), / )'
      write(unit_debug_, fmt) 'WZ. LIMITS  = '
      write(unit_debug_, fmt2) wd%limits_wz_


      write(fmt2, '(a, i2, a)') '( ', wd%nz_, '( 3(3g12.4, 1x, /) ) )'
      write(unit_debug_, fmt) 'WZ. LXYZ    = '
      write(unit_debug_, fmt2) wd%turb_scales_wz_
      write(unit_debug_, fmt) 'WZ. CORR COF= '
      write(unit_debug_, fmt2) wd%corrCoeffs_wz_
      write(unit_debug_, fmt) 'WZ. CORR EXP= '
      write(unit_debug_, fmt2) wd%corrExp_wz_
      write(unit_debug_, fmt) 'WZ. lROTW2G = '
      write(unit_debug_, fmt2) wd%rot_LW2G_wz_
   end subroutine

end submodule
