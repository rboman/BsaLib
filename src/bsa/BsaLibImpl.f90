!! This file is part of BSA Library.
!! Copyright (C) 2023  Michele Esposito Marzino 
!!
!! BSA Library is free software: you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation, either version 3 of the License, or
!! (at your option) any later version.
!!
!! BSA Library is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with BSA Library.  If not, see <https://www.gnu.org/licenses/>.
submodule(BsaLib) BsaLib_Impl

   use BsaLib_Data
   use BsaLib_Utility
   use BsaLib_IO
   implicit none (type, external)

   character(len = :), allocatable :: out_dir_ ! output directory

   logical :: only_diag_elems_ = .false.
   logical :: is_only_msh_     = .false.
   logical :: header_called_   = .false.


contains


   module subroutine bsa_printBSAHeader()
      print *
      print *
      print *, '              ____________________________________________ '
      print *, '             |     _____         ____                     |'
      print *, '             |      /   \       /              /\         |'
      print *, '             |     /____/       \___          /  \        |'
      print *, '             |    /    \            \        /____\       |'
      print *, '             |   /_____/  .    _____/  .   _/      \_  .  |'
      print *, '             |____________________________________________|'
      print *
      print *

      header_called_ = .true.
   end subroutine bsa_printBSAHeader





   module subroutine bsa_enableGPU()
      is_gpu_enabled_ = .true.
   end subroutine




   module subroutine bsa_Init()
      integer(int32) :: istat
      character(len = 256) :: emsg

      if (.not. header_called_) call bsa_printBSAHeader()

      ! if (.not. allocated(out_dir_)) out_dir_ = BSA_OUT_DIRNAME_DEFAULT
      ! istat = util_createDirIfNotExist(out_dir_)

      call bsa_openFileHandles_()

      if (.not. allocated(settings)) then
         allocate(settings, stat=istat, errmsg=emsg)
         if (istat /= 0) call allocKOMsg('settings', istat, emsg)
      endif

      if (.not. allocated(wd)) then
         allocate(wd, stat=istat, errmsg=emsg)
         if (istat /= 0) call allocKOMsg('wd', istat, emsg)
      endif

      if (.not. allocated(struct_data)) then
         allocate(struct_data, stat=istat, errmsg=emsg)
         if (istat /= 0) call allocKOMsg('struct_data', istat, emsg)
      endif

      if (.not. allocated(timer)) then
         allocate(timer, stat=istat, errmsg=emsg)
         if (istat /= 0) call allocKOMsg('timer', istat, emsg)
      endif
   end subroutine



   module subroutine bsa_forceBsaClsExecution(bool)
      logical, intent(in) :: bool
      force_cls_execution_ = bool
   end subroutine



   module subroutine bsa_setMaxBkgPeakRestriction(bool)
      logical, intent(in) :: bool
      do_restrict_bkgpeak_ = bool
   end subroutine




   module subroutine bsa_setPODTruncationThreshold(rval)
      real(real64), intent(in) :: rval

      if (rval > 0.0_real64) then
         do_trunc_POD_ = .true.
         POD_trunc_lim_ = rval / 100.0_real64
      endif
   end subroutine


   module subroutine bsa_setPODNOfModesKept(nmodes)
      integer(int32), intent(in) :: nmodes

      nPODmodes_set_ = .true.
      nmodes_POD_    = nmodes
   end subroutine




   module subroutine bsa_Run(m2mf_cls, m2mr_cls, m2o2mr_cls, m3mf_msh, m3mr_msh, m3mf_cls, m3mr_cls)
      use BsaLib_Functions
      real(bsa_real_t), target, allocatable, dimension(:) :: &
         m2mf_cls, m2mr_cls, m2o2mr_cls, m3mf_msh, m3mr_msh, m3mf_cls, m3mr_cls

      ! user asked for nothing ??
      if (settings%i_compute_psd_ == 0 .and. settings%i_compute_bisp_ == 0) then
         print '(1x, 2a)', &
            WARNMSG, 'Both  PSD  and  BISP  computation are disabled.'
         return
      endif


      print *
      block
         integer(int32) :: itmp

         call io_setExportSpecifiers()

         call validateAll_() ! check before doing some bad things..

         ! NOTE: recall this function since validateAll_() might have changed some..
         call setBsaFunctionLocalVars()

         call io_printUserData()

! BUG: force GPU usage for the moment. Testing only
#ifdef BSA_USE_GPU
         is_gpu_enabled_ = .true.
#endif


         if (is_gpu_enabled_) then
#ifdef BSA_USE_GPU
# ifdef BSA_USE_CUDA
            call bsacl_AcquirePSDId(wd%i_psd_type_)
# endif
            call bsacl_AcquireStructModMat(struct_data%modal_%phi_, struct_data%modal_%nat_freqs_)
            call bsacl_AcquireModalMatrices(struct_data%modal_%Mm_, struct_data%modal_%Cm_, struct_data%modal_%Km_)
            call bsacl_AcquireLoadedNodesList(struct_data%n_load_)
            call bsacl_AcquireTotalNOfNodes(struct_data%nn_)
            call bsacl_AcquireUsedModesList(struct_data%modal_%modes_)
            call bsacl_AcquireWindCoeffs(wd%wfc_)
            call bsacl_AcquirePhiTimesCMat(PHItimesC_local_)
            call bsacl_AcquireTurbComponentsList(wd%tc_)

# ifdef BSACL_ENABLE_EVALFCT_PTR
            call bsacl_AcquireEvaluationFunc(evaluatePSD)
# endif

            call bsacl_AcquireNodalCorrelation(wd%nod_corr_)

            call bsacl_AcquireWindNodalVelocities(wd%u_node_)
            call bsacl_AcquireWindNodalWindZones(wd%wz_node_)
            call bsacl_AcquireWindTurbScales(wd%turb_scales_wz_, wd%nz_)
            call bsacl_AcquireWindTurbStd(wd%sigmaUVW_wz_, wd%nz_)

            call bsacl_SetDeviceType(BSACL_DEVICE_TYPE_GPU)

            ierr_cl_ = bsacl_Init(1)
            if (ierr_cl_ /= 0) then
               print '(1x, a, a)', ERRMSG, 'Failed to initialise BsaCL.'
               goto 998
            endif

            ierr_cl_ = bsacl_InitDeviceMemory()
            if (ierr_cl_ /= 0) then
               print '(1x, a, a)', ERRMSG, 'Failed to initialise BSACL device memory.'
               goto 998
            endif
#else
            print '(1x, 2a)', &
               WARNMSG, 'This BSA version does not support GPU offloading. Using CPU only.'
            is_gpu_enabled_ = .false.
#endif ! BSA_USE_GPU
         endif ! gpu enabled


         ! OK. check done.
         if (settings%i_only_diag_ == 1) then

            only_diag_elems_ = .true.

            itmp        = struct_data%nn_load_ * struct_data%nlibs_load_
            dimNf_psd_  = itmp
            dimNf_bisp_ = itmp
            dimM_psd_   = struct_data%modal_%nm_eff_
            dimM_bisp_  = dimM_psd_
            dimNr_psd_  = struct_data%ndofs_
            dimNr_bisp_ = struct_data%ndofs_


            ! classic fct pointers
            getBFM_vect_cls => getFM_diag_tnlm_vect_cls_
            getBRM_vect_cls => getRM_diag_vect_cls_

            getBFM_scalar_cls => getFM_diag_tnlm_scalar_cls_
            getBRM_scalar_cls => getRM_diag_scalar_cls_


            ! mesher fct pointers
            getBFM_msh => getFM_diag_tnm_scalar_msh_
            getBRM_msh => getRM_diag_scalar_msh_

         else ! == 0, full

            itmp        = struct_data%nn_load_ * struct_data%nlibs_load_
            dimNf_psd_  = itmp * itmp ! itmp^2
            dimNf_bisp_ = dimNf_psd_ * itmp ! itmp^3

            dimM_psd_   = struct_data%modal_%nm_eff_**2
            dimM_bisp_  = dimM_psd_ * struct_data%modal_%nm_eff_  ! nm^3

            itmp        = struct_data%ndofs_
            dimNr_psd_  = itmp ! * itmp
            dimNr_bisp_ = itmp ! * dimNr_psd_


            ! classic fct pointers
            getBFM_vect_cls => getFM_full_tnm_vect_cls_
            getBRM_vect_cls => getRM_full_vect_cls_

            getBFM_scalar_cls => getFM_full_tnm_scalar_cls_
            getBRM_scalar_cls => getRM_full_scalar_cls_


            ! mesher fct pointers
            if (settings%i_use_svd_ == 0) then
               getBFM_msh => getFM_full_tnm_scalar_msh_
            else
               getBFM_msh => getFM_full_tm_scalar_msh_POD_
            endif
            getBRM_msh => getRM_full_scalar_msh_
         end if


         if (.not. do_export_brm_ .and. .not. associated(write_brm_fptr_)) &
            write_brm_fptr_ => exportBRM_void_internal_

#ifdef _BSA_CHECK_NOD_COH_SVD
         settings%i_suban_type_   = 2  ! force MSH execution
#endif


#ifdef BSA_USE_GPU
         ! BUG: this has to be removed. Use GPU also with mesher
         if (is_gpu_enabled_) then
            settings%i_suban_type_   = 1  ! force CLS execution
            settings%i_compute_bisp_ = 1
            settings%i_compute_psd_  = 0
         endif
#endif

         is_only_msh_ = settings%i_suban_type_ == 2
         if (is_only_msh_ .or. settings%i_suban_type_ == 3) &
            call mainMesher_(m3mf_msh, m3mr_msh)

#ifdef _BSA_CHECK_NOD_COH_SVD
         goto 998
#endif

         ! NOTE: in case we cannot have 2nd order moments, force it here
         if (.not. is_only_msh_ .or. force_cls_execution_) then
            if (is_only_msh_) then ! only 2nd order stats
               settings%i_compute_bisp_ = 0
               settings%i_compute_psd_  = 1
            endif
            call mainClassic_(m2mf_cls, m2mr_cls, m2o2mr_cls, m3mf_cls, m3mr_cls)
         endif
      end block

      998 continue

#ifdef BSA_USE_GPU
      if (is_gpu_enabled_) then
         call bsacl_Finalise()
         if (ierr_cl_ == BSACL_SUCCESS) then
            print '(1x, 2a)', INFOMSG, "BSACL returned correctly."
         else
            call bsa_Abort(" BSACL returned with error.")
         endif
      endif
#endif
   end subroutine bsa_Run






   subroutine validateAll_()

      call setExportPathPrefix_()

      if (do_trunc_POD_) then
         if (POD_trunc_lim_ == 0.0_real64 .or. POD_trunc_lim_ == 1.0_real64)   do_trunc_POD_  = .false.
      endif
      if (nPODmodes_set_) then
         if (nmodes_POD_ <= 0_int32 .or. nmodes_POD_ >= struct_data%nn_load_) nPODmodes_set_ = .false.
         if (nPODmodes_set_) &
            print '(/ 1x, a, a, i0, a /)', INFOMSG, 'Using  ', nmodes_POD_, '  POD modes.'
      endif
      if (I_BKG_PEAK_DELTAF_BFM_REFMT_FCT_ <= 0) I_BKG_PEAK_DELTAF_BFM_REFMT_FCT_ = 2
      if (I_BKG_PEAK_DELTAF_BFM_REFMT_FCT_ <= 0) I_BKG_PEAK_DELTAF_BFM_REFMT_FCT_ = 3

      associate(ibispsym => settings%i_bisp_sym_)
         if (.not. ( ibispsym == BSA_SPATIAL_SYM_NONE .or. &
                     ibispsym == BSA_SPATIAL_SYM_HALF .or. &
                     ibispsym == BSA_SPATIAL_SYM_FOUR )) then
            print '(1x, a, a, i0)', WARNMSG, 'Unsupported value  "iBispSym"= ', ibispsym
            print '(1x, a, a)',     MSGCONT, 'Valid values are:  0 (FULL, default), 2 (HALF), 4 (FOURTH).'
            print '(1x, a, a)',     MSGCONT, 'Setting default value.'
            ibispsym = BSA_SPATIAL_SYM_NONE
         endif

         if (ibispsym == BSA_SPATIAL_SYM_FOUR .and. settings%i_3d_sym_ == 1) then
            print '(1x, a, a)', WARNMSG, 'Cannot use 3D matrix symmetry if computing only 1/4 in space.'
            print '(1x, a, a)', MSGCONT, 'Disabling it..'
            settings%i_3d_sym_ = 0
         endif
      end associate


      if (struct_data%ndofs_ == 0) struct_data%ndofs_ = struct_data%nn_ * struct_data%nlibs_
      if (struct_data%ndofs_ /= size(struct_data%modal_%phi_, 1)) call bsa_Abort("NDOFs does not match modal matrix size.")
      if (do_validate_modal_) call validateModalInfo_()

      if (.not. allocated(struct_data%bkg_peak_width_)) then
         block
            real(bsa_real_t) :: vtmp(3, 3, wd%nz_)
            real(bsa_real_t) :: vtmp2(3, 3)
            integer(int32)   :: itmp

            do itmp = 1, wd%nz_
               vtmp(:, :, itmp) = wd%turb_scales_wz_(:, :, itmp) / wd%u_mean_ref_wz_(itmp)
            enddo

            vtmp2 = maxval(vtmp, dim=3)

            call struct_data%computeBKGPeakWidths(vtmp2)
         end block
      endif


      ! BUG: check correctly
      if (.not. (allocated(wd%tc_) .and. allocated(wd%dirs_))) &
         call wd%SetTurbCompsAndDirsDefault()

      if (.not. associated(wd%phi_times_A_ndegw_)) then
         print '(1x, a, a)', &
            WARNMSG, 'Using local PHItimesC instance.'
         print '(1x, a, a /)', &
            MSGCONT, 'Consider using external for less memory usage.'
         call setPhitimesCLocalInstance_()
      endif
   end subroutine validateAll_




   !> Modal info validation step.
   !> Mainly to avoid having mode shapes non 1-normalised
   !> (i.e. modes in torsion, etc..)
   !> NOTE: here is where NMODES_EFF is actually set.
   !>       Better not to give user the chance to do it.
   subroutine validateModalInfo_()
      real(bsa_real_t), dimension(struct_data%modal_%nm_) :: maxvals
      integer(int32) :: i, j, nskip, ilocmax(1), ilib
      integer(int32) :: nmk, istat
      integer(int32), allocatable :: modesk(:)
      character(len = 256) :: emsg

      maxvals = maxval(abs(struct_data%modal_%phi_), dim=1)
      nskip   = 0
      do i = 1, struct_data%modal_%nm_

         if (maxvals(i) == 1._bsa_real_t) cycle

         ! ! NOTE: might also be greater than 1..
         ! if (maxvals(i) < 1._RDP .or. maxvals(i) > 1._RDP + incr) then

         print '(1x, a, a, i3, a)', &
            WARNMSG, 'Mode  ', i, '  is not 1-normalised. Ignoring it..'

         nskip = nskip + 1

         ! check at which lib happens the MAX
         ilocmax = maxloc(abs(struct_data%modal_%phi_(:, i)), kind = int32)
         ilib    = mod(ilocmax(1), struct_data%nlibs_)
         if (ilib == 0) ilib = struct_data%nlibs_
         print '(1x, a, a, i0 /)', MSGCONT, 'Its local max abs value is for LIB= ', ilib
         ! endif
      enddo

      if (nskip == 0) then
         struct_data%modal_%nm_eff_ = struct_data%modal_%nm_
         if (.not. allocated(struct_data%modal_%modes_)) call struct_data%SetKeptModesDefault()
         return
      endif

      ! we have found some mode shapes not 1-normalised
      nmk = struct_data%modal_%nm_ - nskip

      allocate(modesk(nmk), stat=istat, errmsg=emsg)
      if (istat /= 0) call allocKOMsg('modesk', istat, emsg)

      j = 1
      do i = 1, struct_data%modal_%nm_

         if (maxvals(i) == 1._bsa_real_t) then
            modesk(j) = i
            j = j + 1
         endif
      enddo

      struct_data%modal_%nm_eff_ = nmk
      struct_data%modal_%modes_  = int(modesk, kind=bsa_int_t)
      deallocate(modesk, stat=istat)
   end subroutine validateModalInfo_






   !> This routine sets WindData PHItimesC variable
   !> to point to a locally computed (and allocated)
   !> instance, instead of externally acquired one.
   !> This is to avoid memory error if NMODES_EFF < NMODES
   !> but externally PHItimesC is allocated using NMODES.
   subroutine setPhitimesCLocalInstance_()
      integer(int32) :: ndegw, nlib_l, nnodes_l, nmodes, ndofs
      integer(int32) :: id, in, im, n, m, skip
      integer(int32) :: istat
      character(len = 256) :: emsg

      if (.not. associated(wd%wfc_)) &
         call bsa_Abort('Wind force coefficients were not acquired. Aborting.')

      nlib_l  = size(wd%wfc_, 1)
      ndegw   = size(wd%wfc_, 2)
      nnodes_l= size(wd%wfc_, 3)

      if (nnodes_l /= struct_data%nn_load_) call bsa_Abort('N. of loaded nodes does not match.')

      if (nlib_l /= struct_data%nlibs_load_) call bsa_Abort('N. of loaded libs  does not match.')

      if (.not. associated(struct_data%modal_%phi_)) call bsa_Abort('Modal matrix was not acquired. Aborting.')

      ndofs  = size(struct_data%modal_%phi_, 1)
      nmodes = size(struct_data%modal_%phi_, 2)

      if (ndofs /= struct_data%nn_ * struct_data%nlibs_) call bsa_Abort('Incorrect value for  ndofs.')

      if (nmodes /= struct_data%modal_%nm_) &
         call bsa_Abort('N. of modes does not match between local PHI instance and nm_ value.')

      if (struct_data%modal_%nm_eff_ == 0) &
         call bsa_Abort('Effective number of modes to be kept is 0. Validate modal info before this step.')

      if (struct_data%nn_load_ == 0) call bsa_Abort('No nodes are loaded. Aborting.')

      allocate(PHItimesC_local_(struct_data%modal_%nm_eff_, nnodes_l, ndegw), &
         stat=istat, errmsg=emsg)
      if (istat /= 0) call allocKOMsg('PHItimesC_local_', istat, emsg)
      PHItimesC_local_ = 0._bsa_real_t


      do in = 1, nnodes_l

         n    = struct_data%n_load_(in)
         skip = (n - 1) * struct_data%nlibs_

         do id = 1, ndegw
            do im = 1, struct_data%modal_%nm_eff_

               m = struct_data%modal_%modes_(im)

               PHItimesC_local_(im, in, id) = PHItimesC_local_(im, in, id) + &
                  sum(wd%wfc_(:, id, in) * struct_data%modal_%phi_(skip + struct_data%libs_load_, m))
            enddo
         enddo
      enddo

      call wd%SetPhitimesC(PHItimesC_local_)

#ifdef BSA_DEBUG
      write(unit_debug_, '(1x, a, a)') &
         INFOMSG, '@BsaLibImpl::setPhitimesCLocalInstance_() : local PhiTimesC instance computed -- ok.'
#endif
   end subroutine setPhitimesCLocalInstance_





   logical pure module function bsa_isCleaned()
      bsa_isCleaned = is_data_cleaned_
   end function



   module subroutine bsa_Finalise()
      call cleanBSAData_()
   end subroutine bsa_Finalise







!=========================================================================================
!=========================================================================================
!=========================================================================================
!
!  SETTERS  section
!
!=========================================================================================
!=========================================================================================
!=========================================================================================




   module subroutine bsa_setOutputDirectory(dirname)
      character(len=*), intent(in) :: dirname

      out_dir_ = io_appendFilesep(dirname)
   end subroutine


   module subroutine bsa_setExportDirectory(dirname)
      !! Sets export directory to a different path than outdir.
      character(len = *), intent(in) :: dirname

      call io_setExportDirectory(dirname)
   end subroutine


   module subroutine bsa_setExportInCurrDir()
      call io_setExportInCurrDir()
   end subroutine


   module subroutine bsa_setOutUnit(iunit)
      integer(bsa_int_t), intent(in) :: iunit

      unit_debug_ = iunit
   end subroutine

   module subroutine bsa_setOutFileName(fname)
      character(len=*), intent(in) :: fname

      undebug_fname_ = fname
   end subroutine



   subroutine setExportPathPrefix_()
      if (export_in_cwd_) then
         exp_dir_ = ''
         out_dir_ = ''
      else
         if (allocated(exp_dir_)) return
         if (.not. allocated(out_dir_)) out_dir_ = BSA_OUT_DIRNAME_DEFAULT
         exp_dir_ = out_dir_
      endif
   end subroutine 



   module subroutine bsa_setSpatialSymmetry(isym)
      integer(bsa_int_t), intent(in) :: isym

      select case (isym)
         case (BSA_SPATIAL_SYM_NONE)
            settings%i_bisp_sym_ = isym
         case (BSA_SPATIAL_SYM_FOUR)
            settings%i_bisp_sym_ = isym
         case default
            if (.not. isym == BSA_SPATIAL_SYM_HALF) then
               print '(1x, 2a, i0, a)', &
                  WARNMSG, 'Invalid   ', isym, '  spatial symmetry value.'
               print '(1x, 2a)', &
                  MSGCONT, 'Setting default (HALF).'
            endif
            settings%i_bisp_sym_ = BSA_SPATIAL_SYM_HALF
      end select
   end subroutine


   module subroutine bsa_setBfmMLR(bool)
      logical, intent(in) :: bool

      test_no_bfm_mlr_ = .not. bool
   end subroutine


   module subroutine bsa_setPremeshType(itype)
      integer(bsa_int_t), intent(in) :: itype

      select case (itype)
         case (BSA_PREMESH_TYPE_DIAG_CREST_YES)
            ipre_mesh_type = itype
         case default
            if (.not. itype == BSA_PREMESH_TYPE_DIAG_CREST_NO) then
               print '(1x, 2a, i0, a)', &
                  WARNMSG, 'Invalid   ', itype, '  pre-meshing type value.'
               print '(1x, 2a)', &
                  MSGCONT, 'Setting default (NO DIAG).'
            endif
      end select
   end subroutine


   module subroutine bsa_setPremeshMode(imode)
      integer(bsa_int_t), intent(in) :: imode

      select case (imode)
         case (BSA_PREMESH_MODE_BASE)
            ipre_mesh_mode = imode
         case default
            if (.not. imode == BSA_PREMESH_MODE_ZONE_REFINED) then
               print '(1x, 2a, i0, a)', &
                  WARNMSG, 'Invalid   ', imode, '  pre-meshing mode value.'
               print '(1x, 2a)', &
                  MSGCONT, 'Setting default (ZONE REFINED).'
            endif
      end select
   end subroutine



   module subroutine bsa_doValidateModalData(bool)
      logical, intent(in) :: bool

      do_validate_modal_ = bool
   end subroutine



   ! module subroutine bsa_doValidateZoneDeltas(bool)
   !    logical, intent(in) :: bool

   !    do_validate_deltas_ = bool
   ! end subroutine


   module subroutine bsa_setValidateDeltasPolicy(id)
      integer, intent(in) :: id

      select case (id)
      case (BSA_VALIDATE_DELTAS_POLICY_NONE)
         do_validate_deltas_ = .false.
         return
      case (BSA_VALIDATE_DELTAS_POLICY_LIGHT)
         I_BKG_PEAK_DELTAF_BFM_REFMT_FCT_ = 2
         I_BKG_PEAK_DELTAF_BFM_REFMT_FCT_ = 2
      case (BSA_VALIDATE_DELTAS_POLICY_MEDIUM)
         I_BKG_PEAK_DELTAF_BFM_REFMT_FCT_ = 3
         I_BKG_PEAK_DELTAF_BFM_REFMT_FCT_ = 3
      case (BSA_VALIDATE_DELTAS_POLICY_HIGH)
         I_BKG_PEAK_DELTAF_BFM_REFMT_FCT_ = 4
         I_BKG_PEAK_DELTAF_BFM_REFMT_FCT_ = 4
      case (BSA_VALIDATE_DELTAS_POLICY_STRICT)
         I_BKG_PEAK_DELTAF_BFM_REFMT_FCT_ = 5
         I_BKG_PEAK_DELTAF_BFM_REFMT_FCT_ = 5
      case default
         I_BKG_PEAK_DELTAF_BFM_REFMT_FCT_ = 2
         I_BKG_PEAK_DELTAF_BFM_REFMT_FCT_ = 3
      end select

      if (.not.do_validate_deltas_) &
         do_validate_deltas_ = .true.
   end subroutine


   module subroutine bsa_setValidateDeltasValues(ibkg, ires)
      integer, intent(in) :: ibkg, ires

      I_BKG_PEAK_DELTAF_BFM_REFMT_FCT_ = ibkg
      I_RES_PEAK_DELTAF_BFM_REFMT_FCT_ = ires
   end subroutine


   module subroutine bsa_closeUnitsAtEnd()
      close_deb_unit_ = .true.
   end subroutine



   module subroutine bsa_setExportFileFormat(iform)
      integer(bsa_int_t), intent(in) :: iform

      call io_setExportFileFormat(iform)
   end subroutine


   module subroutine bsa_setExportAppendMode(imode)
      integer(bsa_int_t), intent(in) :: imode

      call io_setExportAppendMode(imode)
   end subroutine







   ! TODO: settings might be set via direct assignment ??


   !=====================================
   !   SETTINGS
   !=====================================


   module elemental function bsa_isFullComp() result(bool)
      logical :: bool

      bool = .not. only_diag_elems_
   end function


   module subroutine bsa_setSubanType(isuban)
      integer(bsa_int_t), intent(in) :: isuban

      call settings%SetSubanType(isuban)
   end subroutine


   module subroutine bsa_setVersion(ivers)
      integer(bsa_int_t), intent(in) :: ivers

      call settings%SetVersion(ivers)
   end subroutine


   module subroutine bsa_setScalingConv(iconv)
      integer(bsa_int_t), intent(in) :: iconv

      call settings%SetScalingType(iconv)
   end subroutine


   module subroutine bsa_setSpectraComputation(ipsd, ibisp)
      integer(bsa_int_t), intent(in), optional :: ipsd, ibisp

      call settings%ActivateSpectraComputation(ipsd, ibisp)
   end subroutine



   module subroutine bsa_setSpectraExtension(ionlydiag)
      integer(bsa_int_t), intent(in) :: ionlydiag

      call settings%SetExtension(ionlydiag)
   end subroutine

   module subroutine bsa_setTestMode(itest)
      integer(bsa_int_t), intent(in) :: itest

      call settings%TestMode(itest)
   end subroutine


   module subroutine bsa_setSymmetries(ibispsym, i3dsym)
      integer(bsa_int_t), intent(in) :: ibispsym, i3dsym

      call settings%setSymmetries(ibispsym, i3dsym)
   end subroutine


   module subroutine bsa_setupClassic(nfreqs, df)
      integer(bsa_int_t), intent(in) :: nfreqs
      real(bsa_real_t), intent(in)  :: df

      call settings%setClsSettings(nfreqs, df)
   end subroutine

   module subroutine bsa_setupMesher(isvd, bkgrfmt, bkgaext, genpaext, maxaext, ifcov, idumpmod)
      integer(bsa_int_t), intent(in) :: isvd, bkgrfmt, maxaext
      integer(bsa_int_t), intent(in) :: bkgaext, genpaext, ifcov, idumpmod

      call settings%SetMshrSetts(isvd, bkgrfmt, bkgaext, genpaext, maxaext, ifcov, idumpmod)
   end subroutine






   !=====================================
   !   WIND DATA
   !=====================================

   module subroutine bsa_setWindDirections(dirs, ndirs)
      integer(bsa_int_t), intent(in) :: dirs(:)
      integer(bsa_int_t), intent(in), optional :: ndirs

      call wd%setWindDirections(dirs, ndirs)
   end subroutine


   module subroutine bsa_setWindTurbComps(tc, ntc)
      integer(bsa_int_t), intent(in) :: tc(:)
      integer(bsa_int_t), intent(in), optional :: ntc

      call wd%setTurbComps(tc, ntc)
   end subroutine



   module subroutine bsa_setWindVertProf(iwprof)
      integer(bsa_int_t), intent(in) :: iwprof
      call wd%SetWindvertProf(iwprof)
   end subroutine


   module subroutine bsa_setPSDType(ipsd)
      integer(bsa_int_t), intent(in) :: ipsd

      call wd%SetPSDType(ipsd)
   end subroutine


   module subroutine bsa_setWindAltDir(ivert)
      integer(bsa_int_t), intent(in) :: ivert

      call wd%SetMainVertDir(ivert)
   end subroutine


   module subroutine bsa_setWindZoneLimits(lim, ilim)
      real(bsa_real_t), intent(in) :: lim(..)
      integer(bsa_int_t), intent(in), optional :: ilim(..)

      call wd%SetWindZoneLimits(lim, ilim)
   end subroutine


   module subroutine bsa_setAirDensity(aird)
      real(bsa_real_t), intent(in) :: aird

      call wd%SetAirDensity(aird)
   end subroutine


   module subroutine bsa_setGlobalRotMatW2G(rotW2G)
      real(bsa_real_t), intent(in) :: rotW2G(3, 3)

      call wd%SetGlobalW2G(rotW2G)
   end subroutine


   module subroutine bsa_setWZMeanWindVel(mat)
      real(bsa_real_t), target, intent(in) :: mat(:)

      call wd%SetWZMeanWindVel(mat)
   end subroutine


   module subroutine bsa_setWZRefAlt(Zref)
      real(bsa_real_t), target, intent(in) :: Zref(:)

      call wd%SetWZRefAlt(Zref)
   end subroutine


   module subroutine bsa_setTurbWindScales(L)
      real(bsa_real_t), target, intent(in) :: L(3, 3, *)

      call wd%SetTurbWindScales(L)
   end subroutine


   module subroutine bsa_setTurbWindSDT(sigma)
      real(bsa_real_t), target, intent(in) :: sigma(3, *)

      call wd%SetTurbWindSDT(sigma)
   end subroutine


   module subroutine bsa_setWindCorrCoeffs(ccoeffs)
      real(bsa_real_t), target, intent(in) :: ccoeffs(3, 3, *)

      call wd%SetWindCorrCoeffs(ccoeffs)
   end subroutine


   module subroutine bsa_setWindCorrExpnts(cexpn)
      real(bsa_real_t), target, intent(in) :: cexpn(3, 3, *)

      call wd%SetWindCorrExpnts(cexpn)
   end subroutine


   module subroutine bsa_setIncidenceAngles(incang)
      real(bsa_real_t), target, intent(in) :: incang(:)

      call wd%SetIncidenceAngles(incang)
   end subroutine


   module subroutine bsa_setWZRotMatW2G(rotW2G_L)
      real(bsa_real_t), target, intent(in) :: rotW2G_L(3, 3, *)

      call wd%SetLocalRotMatW2G(rotW2G_L)
   end subroutine



   module subroutine bsa_setNodalVel(Unod)
      real(bsa_real_t), target, intent(in) :: Unod(:)

      call wd%SetNodalVel(Unod)
   end subroutine


   module subroutine bsa_setNodalWindZones(NodWZ)
      integer(bsa_int_t), target, intent(in) :: NodWZ(:)

      call wd%SetNodalWindZones(NodWZ)
   end subroutine


   module subroutine bsa_setNodalWindAltitudes(WnodAlt)
      real(bsa_real_t), target, intent(in) :: WnodAlt(:)

      call wd%SetNodalWindAltitudes(WnodAlt)
   end subroutine


   module subroutine bsa_setSpatialNodalCorr(nodCorr)
      real(bsa_real_t), target, intent(in) :: nodCorr(:, :)

      call wd%SetSpatialNodalCorr(nodCorr)
   end subroutine



   module subroutine bsa_setWindFCoeffs(wfc)
      real(bsa_real_t), target, intent(in) :: wfc(:, :, :)

      call wd%SetWindFCoeffs(wfc)
   end subroutine



   module subroutine bsa_setPhitimesC(phiTc)
      real(bsa_real_t), target, intent(in) :: phiTc(:, :, :)

      call wd%SetPhitimesC(phiTc)
   end subroutine bsa_setPhitimesC




   !=====================================
   !   NODAL DATA
   !=====================================

   module subroutine bsa_setNodalCoords(nn, coords)
      integer(bsa_int_t), intent(in)  :: nn
      real(bsa_real_t), target, allocatable :: coords(:, :)

      call struct_data%SetNodalCoords(nn, coords)
   end subroutine



   module subroutine bsa_setNodalNOfDOFs(nlibs)
      integer(bsa_int_t), intent(in) :: nlibs

      call struct_data%SetNOfNodalDOFs(nlibs)
   end subroutine



   module subroutine bsa_setTotalNOfNodes(nn)
      integer(bsa_int_t), intent(in) :: nn

      call struct_data%SetTotalNOfNodes(nn)
   end subroutine



   module subroutine bsa_setLoadedNodalDOFs(libs_l, nlibs_l)
      integer(bsa_int_t), intent(in), target, allocatable :: libs_l(:)
      integer(bsa_int_t), intent(in), optional :: nlibs_l
      integer(bsa_int_t) :: siz

      if (.not. allocated(libs_l)) return

      if (.not. present(nlibs_l)) then
         siz = size(libs_l)
      else
         if (.not. (nlibs_l == size(libs_l))) &
            call bsa_Abort('Passed number of loaded LIBs does not match size of array.')
         siz = nlibs_l
      endif

      call struct_data%SetLoadedNodalDOFs(siz, libs_l)
   end subroutine



   module subroutine bsa_setLoadedNodes(nodes_l, nn_l)
      integer(bsa_int_t), intent(in), target, allocatable :: nodes_l(:)
      integer(bsa_int_t), intent(in), optional :: nn_l
      integer(bsa_int_t) :: siz

      if (.not. allocated(nodes_l)) return

      if (.not. present(nn_l)) then
         siz = size(nodes_l)
      else
         if (.not. (nn_l == size(nodes_l))) &
            call bsa_Abort('Passed number of loaded LIBs does not match size of array.')
         siz = nn_l
      endif

      call struct_data%SetLoadedNodes(siz, nodes_l)
   end subroutine









   !=====================================
   !   MODAL DATA
   !=====================================

   module subroutine bsa_setModalInfo(ndofs, nm, Phi, natf)
      integer(bsa_int_t), intent(in) :: ndofs, nm
      real(bsa_real_t), intent(in), target :: Phi(ndofs, nm), natf(nm)

      call struct_data%SetModalInfo(ndofs, nm, Phi, natf)
   end subroutine



   module subroutine bsa_setModalMatrices(nm, Mgen, Kgen, Cgen)
      integer(bsa_int_t), intent(in) :: nm
      real(bsa_real_t), intent(in), target, dimension(nm) :: Mgen, Kgen
      real(bsa_real_t), intent(in), target :: Cgen(nm, nm)

      call struct_data%SetModalMatrices(nm, Mgen, Kgen, Cgen)
   end subroutine


   module subroutine bsa_setKeptModalShapes(modes)
      integer(bsa_int_t), intent(in) :: modes(:)

      call struct_data%SetKeptModes(modes)
   end subroutine



   module subroutine bsa_setTotDamping(xsi)
      real(bsa_real_t), target, intent(in) :: xsi(:)

      call struct_data%SetTotDamping(xsi)
   end subroutine


   module pure function bsa_getUsedModeShapes() result(modes)
      integer(bsa_int_t), allocatable :: modes(:)

      modes = struct_data%modal_%modes_
   end function



!=========================================================================================
!=========================================================================================
!=========================================================================================
!
!   COMPUTING   SECTION
!
!=========================================================================================
!=========================================================================================
!=========================================================================================




   module subroutine bsa_computeBRdecomp(m2mf, bkg, res)
      use BsaLib_Functions, only: getBR_SFm_val_
      real(bsa_real_t), intent(in)  :: m2mf(:)
      real(bsa_real_t), allocatable, intent(out) :: bkg(:), res(:)

      integer(int32) :: im, m

      associate(nm => struct_data%modal_%nm_eff_, modes => struct_data%modal_%modes_, &
         Km => struct_data%modal_%Km_, f => struct_data%modal_%nat_freqs_)


         block
            integer(int32) :: istat
            allocate(bkg(nm), stat=istat)
            if (istat /= 0) then
               print '(1x, a, a)', &
                  ERRMSG, 'Cannot allocate resources for BR decomposition computation. Skipping.'
               return
            endif
            allocate(res(nm), stat=istat)
            if (istat /= 0) then
               print '(1x, a, a)', &
                  ERRMSG, 'Cannot allocate resources for BR decomposition computation. Skipping.'
               return
            endif
         end block

         block
            integer   :: idim2, ipsd, ibisp, dimPSD, dimBSP, id
            integer   :: iun, idxi, idxe, itc_, idir_
            real(bsa_real_t) :: fnat, SFm_fnat, rtmp(1), Km_loc2_
            real(bsa_real_t), allocatable :: S_uvw(:, :), S_pad(:)


            ! NOTE: backup this data, for later reset to right values
            !       We want ONLY PSDs here..
            ipsd       = settings%i_compute_psd_
            ibisp      = settings%i_compute_bisp_
            dimPSD     = dimM_psd_
            dimBSP     = dimM_bisp_
            ! NOTE: force these values before calling "getBFM_scalar_cls" fct pointer.
            dimM_psd_  = 1
            dimM_bisp_ = 1
            settings%i_compute_psd_  = 1
            settings%i_compute_bisp_ = 0



            ! BUG: avoid code copy-paste !!
            idim2 = struct_data%nn_load_ * wd%i_ndirs_ * wd%i_ntc_
            allocate(S_uvw(nm, idim2))
            allocate(S_pad(idim2))

            idxi = 1
            idxe = struct_data%nn_load_
            do itc_ = 1, wd%i_ntc_

               do idir_ = 1, wd%i_ndirs_

                  ! BUG:  difference between idir and itc ???
                  ! NOTE: done for all the loaded nodes at once !
                  S_uvw(:, idxi:idxe) = wd%evalPSD(nm, f(modes), &
                     struct_data%nn_load_, struct_data%n_load_, wd%dirs_(idir_), wd%tc_(itc_))

                  idxi = idxe + 1
                  idxe = idxe + struct_data%nn_load_
               enddo ! i direction
            enddo ! i turb comp


            ! do concurrent (im = 1 : nm) local(fnat, SFm_fnat, m, Km_loc2_) &
            !       shared(bkg, res, Km, f, modes, S_pad, S_uvw, rtmp)

            !    m        = modes(im)
            !    Km_loc2_ = Km(m)
            !    Km_loc2_ = Km_loc2_ * Km_loc2_

            !    if (settings%i_only_diag_ == 1) then
            !       bkg(im) = m2mf(im) / Km_loc2_
            !    else
            !       id       = (im-1)*nm + 1
            !       bkg(im) = m2mf(id) / Km_loc2_
            !    endif

            !    fnat = f(m)

            !    ! BUG: adapt back to use already existing procedures..
            !    ! call getBFM_scalar_cls(1, 1, fnat, 0.0_RDP, S_uvw, S_pad, SFm_fnat, rtmp)
            !    call getBR_SFm_val_(nm, S_uvw, fnat, im, m, SFm_fnat)
            !    if (settings%i_def_scaling_ == 1) SFm_fnat = SFm_fnat / CST_PIt4

            !    res(im) = CST_PIGREC * CST_PIt2 * fnat *  SFm_fnat &
            !                / (2 * struct_data%modal_%xsi_(m) * Km_loc2_)
            ! enddo

            do im = 1 , nm

               m        = modes(im)
               Km_loc2_ = Km(m)
               Km_loc2_ = Km_loc2_ * Km_loc2_

               if (settings%i_only_diag_ == 1) then
                  bkg(im) = m2mf(im) / Km_loc2_
               else
                  id      = (im-1)*nm + im
                  bkg(im) = m2mf(id) / Km_loc2_
               endif


               fnat = f(m)
               ! call getBFM_scalar_cls(im, 1, fnat, 0.0_RDP, S_uvw, S_pad, SFm_fnat, rtmp)
               call getBR_SFm_val_(nm, S_uvw, fnat, im, m, SFm_fnat)
               ! if (settings%i_def_scaling_ == 1) SFm_fnat = SFm_fnat / CST_PIt4

               res(im) = CST_PIGREC * CST_PIt2 * fnat *  SFm_fnat &
                           / (2 * struct_data%modal_%xsi_(m) * Km_loc2_)
            enddo

            ! reset old (right) values
            settings%i_compute_psd_  = ipsd
            settings%i_compute_bisp_ = ibisp
            dimM_psd_  = dimPSD
            dimM_bisp_ = dimBSP

         end block
      end associate

   end subroutine bsa_computeBRdecomp







   module subroutine bsa_computePeakFactors(&
         m2, m2o2, obs_time, peak_g, sk, peak_ng_pos, peak_ng_neg)
      real(bsa_real_t), intent(in)  :: m2(:), m2o2(:)
      real(bsa_real_t), intent(in)  :: obs_time
      real(bsa_real_t), allocatable, intent(inout) :: peak_g(:)
      real(bsa_real_t), intent(in), allocatable    :: sk(:)
      real(bsa_real_t), allocatable, intent(inout) :: peak_ng_pos(:)
      real(bsa_real_t), allocatable, intent(inout), optional :: peak_ng_neg(:)

      !> Euler's constant
      real(bsa_real_t), parameter   :: gamma_ = 0.5772d0
      real(bsa_real_t), allocatable :: beta(:)

      if (all(m2o2 == 0)) then
         print '(/ 1x, a, a)', &
            WARNMSG, '"m2_ord2"  is null. Cannot compute extremes. Skipping.'
         return
      endif

      beta = sqrt(m2o2 / m2) / CST_PIt2 * obs_time
      if (any(beta < 1.d0)) then
         print '(1x, a, f7.2, a)', &
            WARNMSG, 'Observation time of ', obs_time, ' sec.  is too short.'
         return
      endif
      beta = sqrt(2 * log(beta))

      peak_g = gamma_ / beta
      peak_g = peak_g + beta

      if (allocated(sk)) then

         block
            real(real64), parameter :: PI2 = CST_PIGREC * CST_PIGREC
            real(bsa_real_t), allocatable :: rtmp(:), g4(:), h3(:), h40(:), h4(:)
            real(bsa_real_t), allocatable :: k(:), beta2(:), beta3(:), pk_ng_neg_(:)

            ! NOTE: excess kurtosis evaluated empirically
            !       based on the "parabolic" relationship with
            !       skewness coefficient, quite acceptable in Wind Engineering
            !       (from conducted experiments).
            !       However, in future, maybe estimate numerically g4 as well!
            rtmp = sk * sk
            g4   = CST_3d2 * rtmp

            ! NOTE: formulation based on Kwon-Kareem-2014-revisited paper (Eq. 1)
            h40 = sqrt(1.d0 + 1.5d0 * g4)
            h3  = sk / (4.d0 +  2.d0 *h40)
            h4  = (h40 - 1.d0) / 18.0d0

            ! ! NOTE: Revised Hermite Model improved formulation (Eqs. 7)
            ! h40 = ((1 + 1.25d0 * g4)**(1.d0 / 3.d0) - 1.d0) / 10.d0
            ! h3  = sk / 6.d0 * (1.d0 - 0.015d0*abs(sk) + 0.3d0*rtmp) / (1.d0 + 0.2d0*g4)
            ! h4  = h40 * (1.d0 - 1.43d0*rtmp/g4)**(1.d0 - 0.1d0 * g4**(0.8d0))

            k     = 1.d0 / (sqrt(1 + 2.d0*h3*h3 + 6.d0*h4*h4))
            beta2 = beta  * beta
            beta3 = beta2 * beta

            ! third term, multipliying h4
            peak_ng_pos = 5.44d0 / (beta3)
            peak_ng_pos = peak_ng_pos + (3.d0 / beta * (PI2 / 6.d0 - gamma_ + (gamma_**2)))
            peak_ng_pos = peak_ng_pos + (beta3 + 3*beta*(gamma_ - 1))
            peak_ng_pos = peak_ng_pos * h4

            pk_ng_neg_  = peak_ng_pos
            rtmp        = h3 * (beta2 + 2.d0 * gamma_ - 1. + 1.98d0 / beta2)
            peak_ng_pos = peak_ng_pos + rtmp
            pk_ng_neg_  = pk_ng_neg_  - rtmp

            peak_ng_pos = peak_ng_pos + peak_g
            peak_ng_pos = peak_ng_pos * k

            if (present(peak_ng_neg)) then
               peak_ng_neg = pk_ng_neg_ + peak_g
               peak_ng_neg = peak_ng_neg * k
            endif
         end block
      endif

   end subroutine









!=========================================================================================
!=========================================================================================
!=========================================================================================
!
!  OUTPUT / EXPORTING
!
!=========================================================================================
!=========================================================================================
!=========================================================================================


   subroutine bsa_openFileHandles_()
      !! BUG: not really adapted to logic...

      ! DEBUG unit
      if (.not. allocated(undebug_fname_)) undebug_fname_ = 'bsadebug.bsa'
      call io_getVerifiedFile(unit_debug_, undebug_fname_)
      open(unit=unit_debug_          & 
         , file=undebug_fname_       &
         , status=IO_STATUS_REPLACE  &
         , form=IO_FORM_FORMATTED    &
         , action=IO_ACTION_WRITE )

#ifdef _BSA_EXPORT_POD_TRUNC_INFO
      open(unit=iun_POD_trunc_       &
         , file=iun_POD_trunc_fname_ &
         , status=IO_STATUS_REPLACE  &
         , form=IO_FORM_UNFORMATTED  &
         , access=IO_ACCESS_STREAM   &
         , action=IO_ACTION_WRITE )
#endif
   end subroutine




   module subroutine bsa_exportBR_nocompute_(fname, bkg, res, xsi)
      !! BUG: adapt to a more general XSI management..
      character(len = *), intent(in) :: fname
      real(bsa_real_t), intent(in)   :: bkg(:), res(:), xsi(:)
      integer(int32) :: iun, im
      ! integer(int32) :: s2, j

      ! s2  = size(bkg, 2)
      iun = io_openExportFileByName(fname)
      if (iun == 0) call bsa_Abort()
      write(iun, *) struct_data%modal_%nm_eff_
      write(iun, *) struct_data%modal_%modes_
      ! write(iun, *) s2
      ! do j = 1, s2
         write(iun, *) xsi(:)
         do im = 1, struct_data%modal_%nm_eff_
            write(iun, *) bkg(im), res(im)
         enddo
      ! enddo
      close(iun)
   end subroutine





   module subroutine bsa_exportMomentToFile(fname, vec)
      character(len = *), intent(in) :: fname
      real(bsa_real_t), intent(in)   :: vec(:)
      integer(int32) :: iun, i, dim

      iun = io_openExportFileByName(exp_dir_ // fname)
      if (iun == 0) call bsa_Abort()
      dim = size(vec)
      write(iun, *) dim
      do i = 1, dim
         write(iun, *) vec(i)
      enddo
      close(iun)
   end subroutine




   module subroutine bsa_exportSkewness_nocompute_(fname, sk)
      character(len = *), intent(in) :: fname
      real(bsa_real_t), intent(in)   :: sk(:)

      call exportSkewness_(fname, sk)
   end subroutine



   module subroutine bsa_exportSkewness_compute_(fname, dim, m2, m3)
      character(len = *), intent(in) :: fname
      integer(bsa_int_t), intent(in) :: dim
      real(bsa_real_t), intent(in)   :: m2(:), m3(:)
      real(bsa_real_t), allocatable  :: sk(:)

      sk = computeSkewness_(dim, m2, m3, only_diag_elems_)
      call exportSkewness_(fname, sk)
   end subroutine




   function computeSkewness_(dim, m2, m3, only_diag) result(sk)
#ifdef BSA_DEBUG
      use, intrinsic :: ieee_arithmetic
#endif
      integer(bsa_int_t), intent(in) :: dim
      real(bsa_real_t), intent(in)   :: m2(:), m3(:)
      logical, intent(in)            :: only_diag

      real(real64), parameter :: cst3d2 = 3._real64 / 2._real64
      real(bsa_real_t), allocatable :: sk(:)

      if (only_diag) then
         sk = m3 / (m2)**(cst3d2)
         return
      endif

      block
         integer(int32) :: szm2, szm3
         integer(int32) :: pm3, pm2
         integer(int32) ::  k,  j,  i, l
         integer(int32) :: ik, ij, ii
         integer(int32) :: s2

         real(bsa_real_t), allocatable :: sigm(:)
         real(bsa_real_t) :: denK, denJ

         szm2 = size(m2, 1)
         szm3 = size(m3, 1)

         allocate(sk(szm3))
         sk = 0._bsa_real_t

         sigm = sqrt(m2)  ! std

         pm3 = 1
         ik  = 1
         do k = 1, dim

            denK = sigm(ik)

            ij = 1
            do j = 1, dim

               denJ = denK * sigm(ij)

               ii = 1
               do i = 1, dim

                  sk(pm3) = m3(pm3) / (denJ * sigm(ii))

#ifdef BSA_DEBUG
                  if (ieee_is_nan(sk(pm3))) then
                     print '(1x, a, a, 2i6)', &
                        ERRMSG, 'SK is NaN at indexes   ', pm3, l
                     goto 99 ! exit loop
                  endif
#endif

                  pm3 = pm3 + 1
                  ii  = i * dim + i + 1
               enddo ! dim i

               ij = j * dim + j + 1
            enddo ! dim j

            ik = k * dim + k + 1
         enddo ! dim k

         99 continue
      end block
   end function computeSkewness_





   subroutine exportSkewness_(fname, vec)
      use BsaLib_Data, only: struct_data
      character(len = *), intent(in) :: fname
      real(bsa_real_t), intent(in)   :: vec(:)
      integer(int32) :: iun, i, dim
      logical :: is_modal
      iun = io_openExportFileByName(exp_dir_ // fname)
      if (iun == 0) call bsa_Abort()

      dim = size(vec)
      is_modal = dim == dimM_bisp_

      if (is_modal) then
         ! modal header
         write(iun, *) struct_data%modal_%nm_eff_
         write(iun, *) struct_data%modal_%modes_
      endif

      ! actual data
      write(iun, *) dim
      do i = 1, dim
         write(iun, *) vec(i)
      enddo
      close(iun)
   end subroutine






   module subroutine bsa_exportPeakOrExtremesToFile(fname, rvar)
      character(len = *), intent(in) :: fname
      real(bsa_real_t), intent(in)   :: rvar(:)
      integer(int32) :: ndofs, iun, i

      iun   = io_openExportFileByName(fname)
      if (iun == 0) call bsa_Abort()

      ndofs = size(rvar)
      write(iun, *) ndofs
      do i = 1, ndofs
         write(iun, *) rvar(i)
      enddo
      close(iun)
   end subroutine




   ! BUG: should also providfe a way to pass pointer to user defined exporting data 
   !      structure that has to be finally dereferenced in actual exporting routine!
   module subroutine bsa_setBRMExportFunction(fptr)
      procedure(exportBRMinterf_vect_), pointer, intent(in) :: fptr

      write_brm_fptr_ => fptr

      ! if user provides its own function, make sure it does not get overridden
      i_brmexport_mode_ = BSA_EXPORT_BRM_MODE_USR
   end subroutine





   subroutine exportBRM_void_internal_(f1, f2, brm, pdata)
      real(bsa_real_t), intent(in)  :: f1(:), f2(:), brm(:, :)
      class(*), pointer, intent(in) :: pdata

      ! do nothing
   end subroutine



   ! BUG: this should be called via a function pointer
   subroutine exportBRM_baseHeaderWriter_internal_(pdata)
      type(BrmExportBaseData_t), pointer, intent(in) :: pdata

      ! BUG: maybe general header written only ONCE in a separate procedure..
      ! do print general header
      if (pdata%i_doNotPrintGenHeader_ == 0) then
         write(unit_dump_brm_) pdata%nm_
         write(unit_dump_brm_) pdata%modes_
         write(unit_dump_brm_) pdata%ncomb_
         write(unit_dump_brm_) pdata%ispsym_
         write(unit_dump_brm_) pdata%nzones_
      endif

      ! do print zone info header
      if (pdata%i_doNotPrintZonHeader_ == 0) then
         write(unit_dump_brm_) pdata%idZone_
         write(unit_dump_brm_) pdata%nI_
         write(unit_dump_brm_) pdata%nJ_

! #ifdef _OPENMP
!          print *, ' Dumping zone with id, ni, nj  =  ', &
!             pdata%idZone_, pdata%nI_, pdata%nJ_
! #endif
      endif
   end subroutine



   subroutine exportBRM_base_internal_(fi, fj, brm, pdata)
      real(bsa_real_t), intent(in)  :: fi(:), fj(:), brm(:, :)
      class(*), pointer, intent(in) :: pdata

      ! Need to verify if to print headers
      if (associated(pdata)) then
         select type (pdata)
            type is (BrmExportBaseData_t)
               call exportBRM_baseHeaderWriter_internal_(pdata)
            class default
               call bsa_Abort("Expecting pdata to be of type  ""BrmExportBaseData_t"".")
         end select
      endif

      block
         integer(int32) :: i, siz

         siz = size(fi)
         do i= 1, siz
            write(unit_dump_brm_) real(fi(i), kind=real32), real(fj(i), kind=real32), real(brm(:, i), kind=real32)
         enddo
      endblock
   end subroutine




   module subroutine bsa_setBRMExportDefaultMode(imode)
      integer(int32), intent(in) :: imode
      integer(int32) :: iost

      select case (imode)

         case (BSA_EXPORT_BRM_MODE_NONE)
            return

         case default  ! includes  (BSA_EXPORT_BRM_MODE_BASE)
            if (.not. imode == BSA_EXPORT_BRM_MODE_BASE) &
               print '(1x, a, a)', WARNMSG, 'Unknown BRM export mode. Setting default  (base).'

            do_export_brm_ = .true.
            write_brm_fptr_  => exportBRM_base_internal_
      end select

      open(&
         unit=unit_dump_brm_,        &
         file=brm_export_file_name_, &
         form=IO_FORM_UNFORMATTED,   &
         access=IO_ACCESS_STREAM,    &
         status=IO_STATUS_REPLACE,   &
         iostat=iost)

      if (iost /= 0) then
         print '(1x, 4a)', &
            WARNMSG, 'Error while opening BRM export file  "', brm_export_file_name_, '".'
         print '(1x, 2a)', MSGCONT, 'Disabling exporting.'
         do_export_brm_  = .false.
         write_brm_fptr_ => null()
      endif
   end subroutine




   module subroutine bsa_saveCoordinatesToFile(fname, coords)
      character(len = *), intent(in) :: fname
      real(bsa_real_t), intent(in), target, optional :: coords(:, :)
      real(bsa_real_t), pointer :: coords_(:, :)
      integer(int32)  :: iun, istat, i, nn_

      if (.not. present(coords) .and. .not. associated(struct_data%coords_)) then
         print '(1x, a, a)', &
            WARNMSG, 'Cannot save coordinates to file. Data not provided. Skipping.'
         return
      endif

      ! TODO: adapt to output in BSA folder
      open(unit=iun, file=fname, form='formatted', action='write', status='replace', &
         iostat=istat)
      if (istat /= 0) return

      if (present(coords)) then
         coords_ => coords
      else
         coords_ => struct_data%coords_
      endif

      nn_ = size(coords_, 2)
      write(iun, *) nn_
      do i = 1, nn_
         write(iun, *) coords_(:, i) ! dims (3, nn)
      enddo
      close(iun)
   end subroutine







   module subroutine bsa_exportPSDToFile(fname, psd, f)
      character(len = *), intent(in) :: fname
      real(bsa_real_t),   intent(in) :: psd(:, :)
      real(bsa_real_t),   intent(in), optional :: f(:)

      real(bsa_real_t), allocatable :: tmp(:)
      integer(int32) :: s1, s2, iun, j

      iun = io_openExportFileByName(exp_dir_ // fname)
      if (iun == 0) call bsa_Abort()

      s1 = size(psd, 1)
      s2 = size(psd, 2)
      write(iun, *) s1
      write(iun, *) s2

      ! NOTE: different from writing bisp
      ! BUG: need this dummy variable in order to avoid 
      !      I/O runtime warning...
      allocate(tmp(s2))
      if (present(f)) then
         do j = 1, s1
            tmp = psd(j, :)
            write(iun, '(*(g0, 1x))') f(j), tmp
         enddo
      else
         do j = 1, s1
            tmp = psd(j, :)
            write(iun, '(*(g0, 1x))') tmp
         enddo
      endif
      close(iun)
      deallocate(tmp)
   end subroutine




   module subroutine bsa_exportBispToFile(fname, bisp)
      character(len = *), intent(in) :: fname
      real(bsa_real_t), intent(in)   :: bisp(:, :, :)

      integer(int32) :: s1, s2, s3, iun, i, j

      iun = io_openExportFileByName(exp_dir_ // fname)
      if (iun == 0) call bsa_Abort()

      s1 = size(bisp, 1)
      s2 = size(bisp, 2)
      if (s2 /= s1) call bsa_Abort('First two dimensions of bisp do not match.')
      s3 = size(bisp, 3)

      write(iun, *) s1
      write(iun, *) s2
      write(iun, *) s3
      do j = 1, s3
         do i = 1, s2
            write(iun, '(*(g0, 1x))') bisp(:, i, j)
         enddo
      enddo
      close(iun)
   end subroutine


end submodule
