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
module BsaLib_Settings

#include "../precisions"
   
   implicit none
   private

   !> Minimum rounding precision to guarantee.
   integer, public :: i_min_round_prec_ = 10

   type, public :: settings_t
   
      ! GENERAL settings


      !> Subanalysis type
      !> ==1, classic
      !> ==2, mesher
      !> ==3, BOTH (comparison) 
      integer(kind = 4) :: i_suban_type_ = 2

      !> NOTE: only for "classic" suban type.
      !> Manage which version to use (specially for dev testing).
      !> ==1 uses old adapted spectra
      !> ==2, uses new spectra, sistematic
      integer(kind = 4) :: i_vers_ = 2

      !> Turbulence PSDs scaling convention
      !> ==1, pulsations  (-infty, +infty)
      !> ==2, frequencies (0, +infty)
      integer(kind = 4) :: i_def_scaling_ = 1

      !> If ==0, do not compute PSDs
      integer(kind = 4) :: i_compute_psd_ = 1

      !> If ==0, do not compute BISPs
      integer(kind = 4) :: i_compute_bisp_ = 1

      !> Whether computing FULL matrices or not.
      !> True if ==1.
      integer(kind = 4) :: i_only_diag_ = 0

      !> Activate for testing some new features.
      integer(kind = 4) :: i_test_mode_ = 0



      ! CLASSIC settings


      !> If suban=="classic", number of sistematic freqs.
      integer(kind = 4) :: nfreqs_ = 0

      !> If suban=="classic", constant delta frequency.
      real(RDP) :: df_ = 0._RDP

      !> If ==0, using VECTORISED functions version.
      !> Otherwise (==1), using SCALAR versions.
      integer(kind = 4) :: i_scalar_vers_ = 0

      !> Bisp symmetry case.
      !> 0 = full
      !> 2 = half
      !> 4 = fourth
      integer(kind = 4) :: i_bisp_sym_ = 0

      !> 3D bisp matrix symmetry exploitation.
      !> 0 = no
      !> 1 = yes
      !> NOTE: if i_bisp_sym_==4, automatically 0
      integer(kind = 4) :: i_3d_sym_ = 0




      ! MESHER settings


      !> If ==1, using SVD to S_uvw matrices
      integer(kind = 4) :: i_use_svd_ = 1

      !> How many points (per side) for meshing 
      !> main central BKG peak zone.
      integer(kind = 4) :: bkg_base_rfmnt_ = 20

      !> 
      integer(kind = 4) :: max_area_extension_ = 2

      !> How much to extend BKG peak area influence.
      integer(kind = 4) :: bkg_area_extension_ = 2

      !> How much to extend general peak area influence.
      integer(kind = 4) :: gen_peak_area_extension_ = 3

      !> If true, we get up to 2*max_freq.
      integer(kind = 4) :: i_full_coverage_ = 1

      !> Controls wheter to include modal info when
      !> writing to dump file.
      !> Unactive by default!
      !> Warn if gets activated.
      integer(kind = 4) :: i_dump_modal_ = 0


   contains

      procedure, public, pass :: SetSubanType
      procedure, public, pass :: SetVersion
      procedure, public, pass :: SetScalingType
      procedure, public, pass :: ActivateSpectraComputation
      procedure, public, pass :: SetExtension
      procedure, public, pass :: TestMode
      procedure, public, pass :: setSymmetries
      procedure, public, pass :: setClsSettings
      procedure, public, pass :: SetMshrSetts
   end type settings_t



   interface

      !> Sets sub analysis type
      module subroutine SetSubanType(this, isuban)
         class(settings_t), intent(inout) :: this
         integer(kind = 4), intent(in) :: isuban
      end subroutine

      !> Set version
      module subroutine SetVersion(this, ivers)
         class(settings_t), intent(inout) :: this
         integer(kind = 4), intent(in) :: ivers
      end subroutine

      !> Sets PSDs scaling convention.
      !>   ==1, default, pulsation (-infty, +infty)
      !>   ==2, frequencies [WARNING] (0, +infty)
      module subroutine SetScalingType(this, idefsc)
         class(settings_t), intent(inout) :: this
         integer(kind = 4), intent(in) :: idefsc
      end subroutine


      !> Allows to control spectra comptation.
      !> Pass 0 to deactivate.
      module subroutine ActivateSpectraComputation(this, ipsd, ibisp)
         class(settings_t), intent(inout) :: this
         integer(kind = 4), intent(in), optional    :: ipsd, ibisp
      end subroutine


      !> Allows to specify whether full 2d/3d matrices are to be computed.
      !> If 0, only main diagonal elements are computed (uncorrelated case).
      module subroutine SetExtension(this, ionlydiag)
         class(settings_t), intent(inout) :: this
         integer(kind = 4), intent(in) :: ionlydiag
      end subroutine


      !> Controls whether testing mode is active.
      module subroutine TestMode(this, itest)
         class(settings_t), intent(inout) :: this
         integer(kind = 4), intent(in) :: itest
      end subroutine


      module subroutine setSymmetries(this, ibispsym, i3dsym)
         class(settings_t), intent(inout) :: this
         integer(kind = 4), intent(in)    :: ibispsym, i3dsym
      end subroutine


      !> Sets main Classic suban settings.
      module subroutine setClsSettings(this, nfreqs, df)
         class(settings_t), intent(inout) :: this
         integer(kind = 4), intent(in)    :: nfreqs
         real(RDP), intent(in) :: df
      end subroutine


      !> Sets main Mesher suban settings.
      module subroutine SetMshrSetts(this, isvd, bkgrfmt, bkgaext, genpaext, maxaext, ifcov, idumpmod)
         class(settings_t), intent(inout) :: this
         integer(kind = 4), intent(in) :: isvd, bkgrfmt, bkgaext, genpaext, maxaext, ifcov, idumpmod
      end subroutine

   end interface


end module