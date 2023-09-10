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
module BsaLib_Data
   
#include "../../precisions"
   
   use Logging
   use BsaLib_Timing
   use BsaLib_CONSTANTS
   use BsaLib_Settings
   use BsaLib_Structure
   use BsaLib_WindData
   !$ use omp_lib
#ifdef __BSA_CL
   use BsaCL
#endif
   implicit none
   public

   
   !====================================
   ! MODULE VARIABLES
   !====================================

   type(settings_t),      allocatable, target :: settings
   type(WindData_t),      allocatable, target :: wd
   type(StructureData_t), allocatable, target :: struct_data
   type(timer_t),         allocatable, target :: timer
   type(logger_t),        allocatable, target :: logger_debug

   logical :: is_data_cleaned_ = .false.
   logical :: close_deb_unit_  = .true.

   logical :: do_validate_modal_ = .true.

   integer(kind = 4) :: dimNf_psd_ = 0, dimNf_bisp_ = 0
   integer(kind = 4) :: dimNr_psd_ = 0, dimNr_bisp_ = 0
   integer(kind = 4) :: dimM_psd_  = 0, dimM_bisp_  = 0
   
   real(RDP), allocatable, target :: PHItimesC_local_(:, :, :)

   real(RDP), allocatable :: peak_exts_(:)
   logical :: do_restrict_bkgpeak_ = .false.

   logical :: do_export_brm_ = .false.
   integer(kind = 4) :: i_brmexport_mode_ = BSA_EXPORT_BRM_MODE_BASE
   character(len = *), parameter :: brm_export_file_name_ = 'bsaexport.brm'
#ifdef __BSA_OMP
   procedure(exportBRMinterf_vect_all_), pointer :: write_brm_fptr_  => null()
#else
   procedure(exportBRMinterf_scalar_),   pointer :: write_brm_fptr_  => null()
#endif
   type, public :: BrmExportBaseData_t
      
      integer(kind = 4) :: i_doNotPrintGenHeader_ = 0   ! == 0  means DO PRINT !!
      integer(kind = 4) :: nm_     = 0
      integer(kind = 4) :: ncomb_  = 0
      integer(kind = 4) :: ispsym_ = 0
      integer(kind = 4) :: nzones_ = 0
      integer(kind = 4), pointer :: modes_(:) => null()

      integer(kind = 4) :: i_doNotPrintZonHeader_ = 0
      integer(kind = 4) :: idZone_ = 0
      integer(kind = 4) :: nI_     = 0
      integer(kind = 4) :: nJ_     = 0
   end type


#ifdef __BSA_CL
   integer, target :: ierr_cl_
#endif


#ifdef __BSA_CHECK_NOD_COH_SVD
   real(RDP), allocatable :: nod_corr_full_(:, :)
   real(RDP), allocatable :: nod_corr_EVLs_(:), nod_corr_EVTs_(:, :)
#endif


   ! =========================
   ! ==== classic related
   !
   logical :: force_cls_execution_ = .false.
   integer(kind = 4), parameter :: MAX_VECT_ALLOC_ELEMS = 1000000000 ! 1B -> almost 8Gb
   integer(kind = 4) :: ifr = 0, jfr = 0
   ! real(RDP), pointer :: m2mf_cls_ptr_(:), m2mr_cls_ptr_(:)   ! 2nd order moments
   ! real(RDP), pointer :: m3mf_cls_ptr_(:), m3mr_cls_ptr_(:)   ! 3rd order moments
   
   procedure(getBFMClsVect), pointer :: getBFM_vect_cls => null()
   procedure(getBRMClsVect), pointer :: getBRM_vect_cls => null()
   abstract interface
      subroutine getBFMClsVect(f, Suvw, psd, bisp)
         import :: RDP
         import :: settings, struct_data, wd
         real(RDP), intent(in) :: f(settings%nfreqs_)
         real(RDP), intent(in) :: Suvw(settings%nfreqs_, struct_data%nn_load_ * wd%i_ndirs_ * wd%i_ntc_)
         real(RDP), allocatable, intent(inout) :: psd(:, :), bisp(:, :, :)
      end subroutine

      subroutine getBRMClsVect(f, psd, bisp)
         import :: RDP
         import :: settings
         real(RDP), intent(in)                 :: f(settings%nfreqs_)
         real(RDP), allocatable, intent(inout) :: psd(:, :), bisp(:, :, :)
      end subroutine
   end interface


   procedure(getBFMClsScalar), pointer :: getBFM_scalar_cls => null()
   procedure(getBRMClsScalar), pointer :: getBRM_scalar_cls => null()
   abstract interface
      pure subroutine getBFMClsScalar(ii, ij, fi, fj, Suvw, Suvw_pad, psd, bisp)
         import :: RDP, dimM_psd_, dimM_bisp_
         import :: settings, struct_data, wd
         integer, intent(in)   :: ii, ij
         real(RDP), intent(in) :: fi, fj
         real(RDP), intent(in) :: Suvw(settings%nfreqs_, struct_data%nn_load_ * wd%i_ndirs_ * wd%i_ntc_)
         real(RDP), intent(in)    :: Suvw_pad(struct_data%nn_load_ * wd%i_ndirs_ * wd%i_ntc_)
         real(RDP), intent(inout) :: psd(dimM_psd_), bisp(dimM_bisp_)
      end subroutine

      subroutine getBRMClsScalar(ii, ij, fi, fj, psdin, psdout, bispin, bispout)
         import :: RDP, dimM_psd_, dimM_bisp_
         integer, intent(in)    :: ii, ij  ! freqs indexes
         real(RDP), intent(in)  :: fi, fj
         real(RDP), intent(in)  :: psdin(dimM_psd_), bispin(dimM_bisp_)
         real(RDP), intent(out) :: psdout(dimM_psd_), bispout(dimM_bisp_)
      end subroutine
   end interface




   ! =========================
   ! ==== mesher related
   !
   real(RDP), pointer :: m3mf_msh_ptr_(:) => null(), m3mr_msh_ptr_(:) => null()

#ifndef __BSA_OMP
   !> Shared instance of undumped BFM.
   !> It holds the max N. of points of all the dumped zones, so that no overflows occur.
   real(RDP), allocatable :: bfm_undump(:, :)
#endif

   integer(kind = 4), public :: ipre_mesh_type = BSA_PREMESH_TYPE_DIAG_CREST_NO
   integer(kind = 4), public :: ipre_mesh_mode = BSA_PREMESH_MODE_ZONE_REFINED
   integer(kind = 4), public :: msh_iZone

   !> Controls if checking zone's deltas or not.
   logical :: do_validate_deltas_ = .true.
   
   integer(kind = 4), public :: I_BKG_PEAK_DELTAF_BFM_REFMT_FCT_ = 2
   integer(kind = 4), public :: I_RES_PEAK_DELTAF_BFM_REFMT_FCT_ = 3
   
   ! Total pre-mesh/post-mesh phase points
   integer(kind = 4), public :: msh_bfmpts_pre_
   integer(kind = 4), public :: msh_bfmpts_post_
   integer(kind = 4), public :: msh_brmpts_post_

   !> Code "-1" means that interest modes have to be read from the NEXT (right)
   !> limit only, since this is the first limit in the list
   integer(kind = 4), public, parameter :: CODE_PRE_PEAK_OK = -1
   
   !> Code "-2" means that interest modes have to be read from 
   !> NEXT (right) limit only, since this is the first limit in the list.
   !> However, info is MISSING from the left side (BKG) peak
   !> in which some modes fall within, but we cannot know which
   !> one is close to this zone's limit, in order to determine if
   !> to be added to its interest modes.
   !> Hence, it is a sort of "ALARM". Info will not be accurate in this case.
   integer(kind = 4), public, parameter :: CODE_PRE_PEAK_KO = -2

   !> Limit zones interest modes indexes
   integer(kind = 4), public, allocatable :: msh_ZoneLimsInterestModes(:)
   
   !> Tot n. of zones counter.
   integer(kind = 4), public, target :: msh_NZones = 0

   !> Controls whether employing new BFM MLR method or not
   logical :: test_no_bfm_mlr_ = .false.

   !> Controls whether to perform modal truncation or not
   logical        :: do_trunc_POD_  = .false.
   real(kind = 8) :: POD_trunc_lim_ = 0.d0

   !> Stores width of background peak
   real(RDP) :: bkg_peakw_ = 0._RDP

   ! Mesher function pointer (pre/post meshing)
   procedure(getMshBFM), pointer :: getBFM_msh => null()
   procedure(getMshBRM), pointer :: getBRM_msh => null()
   abstract interface
      function getMshBFM(fi, fj) result(vals)
         import RDP, dimM_bisp_
         real(RDP), intent(in) :: fi, fj
         real(RDP) :: vals(dimM_bisp_)
      end function

      function getMshBRM(bfm, fi, fj) result(vals)
         import RDP, dimM_bisp_
         real(RDP), intent(in) :: bfm(dimM_bisp_)
         real(RDP), intent(in) :: fi, fj
         real(RDP) :: vals(dimM_bisp_)
      end function
   end interface


   interface
      module function evaluatePSD(f, nf, itc) result(PSD)
         integer(kind = 4), intent(in) :: nf, itc
         real(kind = 8), intent(in)  :: f(nf)
         real(kind = 8), allocatable, target :: PSD(:, :)
      end function

      module subroutine cleanBSAData_()
      end subroutine

      module subroutine bsa_Abort(emsg)
         character(len = *), intent(in), optional :: emsg
      end subroutine
   end interface



end module BsaLib_Data