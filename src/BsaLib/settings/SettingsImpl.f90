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
submodule(BsaLib_Settings) BsaLib_SettingsImpl

#include "../precisions"

   use BsaLib_IO, only: INFOMSG, WARNMSG, ERRMSG, MSGCONT, DBGMSG &
                        , BSA_SETTS_DATA_DUMPFILE, unit_debug_
   use BsaLib_Data, only: bsa_Abort
   use BsaLib_CONSTANTS, only: BSA_SPATIAL_SYM_NONE
   implicit none

contains



   module subroutine SetSubanType(this, isuban)
      class(settings_t), intent(inout) :: this
      integer(kind = 4), intent(in) :: isuban

      if (isuban < 0 .or. isuban > 3) call bsa_Abort('Invalid "sub-an" value.')
      this%i_suban_type_ = isuban
   end subroutine SetSubanType


   module subroutine SetVersion(this, ivers)
      class(settings_t), intent(inout) :: this
      integer(kind = 4), intent(in) :: ivers

      if (ivers < 0 .or. ivers > 2) call bsa_Abort('Invalid "ivers" value.')
      this%i_vers_ = ivers
   end subroutine SetVersion


   module subroutine SetScalingType(this, idefsc)
      class(settings_t), intent(inout) :: this
      integer(kind = 4), intent(in) :: idefsc

      if (idefsc < 0 .or. idefsc > 2) call bsa_Abort('Invalid "idefsc" value.')
      this%i_def_scaling_ = idefsc
   end subroutine SetScalingType


   module subroutine ActivateSpectraComputation(this, ipsd, ibisp)
      class(settings_t), intent(inout) :: this
      integer(kind = 4), intent(in), optional    :: ipsd, ibisp

      if (present(ipsd)) then
         if (ipsd < 0 .or. ipsd > 1) call bsa_Abort('Invalid "ipsd" value.')
         this%i_compute_psd_ = ipsd
      endif

      if (present(ibisp)) then
         if (ibisp < 0 .or. ibisp > 1) call bsa_Abort('Invalid "ibisp" value.')
         this%i_compute_bisp_ = ibisp
      endif
   end subroutine ActivateSpectraComputation



   module subroutine SetExtension(this, ionlydiag)
      class(settings_t), intent(inout) :: this
      integer(kind = 4), intent(in) :: ionlydiag

      if (ionlydiag < 0 .or. ionlydiag > 1) call bsa_Abort('Invalid "ionlydiag" value.')
      this%i_only_diag_ = ionlydiag
   end subroutine SetExtension



   module subroutine TestMode(this, itest)
      class(settings_t), intent(inout) :: this
      integer(kind = 4), intent(in) :: itest

      if (itest < 0 .or. itest > 1) call bsa_Abort('Invalid "itest" value.')
      this%i_test_mode_ = itest
   end subroutine TestMode


   module subroutine setSymmetries(this, ibispsym, i3dsym)
      class(settings_t), intent(inout) :: this
      integer(kind = 4), intent(in)    :: ibispsym, i3dsym

      if (ibispsym == 0) then
         this%i_bisp_sym_ =  BSA_SPATIAL_SYM_NONE
      elseif (ibispsym == 2 .or. ibispsym == 4) then
         this%i_bisp_sym_ = ibispsym
      else
         call bsa_Abort('Invalid "ibispsym" value.')
      endif

      if (i3dsym < 0 .or. i3dsym > 1) call bsa_Abort('Invalid "i3dsym" value.')
      this%i_3d_sym_ = i3dsym
   end subroutine



   module subroutine setClsSettings(this, nfreqs, df)
      class(settings_t), intent(inout) :: this
      integer(kind = 4), intent(in)    :: nfreqs
      real(RDP), intent(in) :: df

      if (nfreqs <= 0) call bsa_Abort('Invalid "nfreqs" value.')
      this%nfreqs_ = nfreqs

      if (df <= 0._RDP) call bsa_Abort('Invalid "df" value.')
      this%df_ = df

#ifdef __BSA_DEBUG
      write(unit_debug_, '(1x, a, a)') &
         INFOMSG, '@SettingsImpl::setClsSettings() : setting Bsa Classic settings -- ok.'
#endif
   end subroutine setClsSettings





   module subroutine SetMshrSetts(this, isvd, bkgrfmt, bkgaext, genpaext, maxaext, ifcov, idumpmod)
      class(settings_t), intent(inout) :: this
      integer(kind = 4), intent(in) :: isvd, bkgrfmt, bkgaext, genpaext, maxaext, ifcov, idumpmod

! #ifdef __BSA_DEBUG
!       write(unit_debug_, *) ' @SettingsImpl::SetMshrSetts() : init setting Bsa Mesher settings...'
! #endif

      this%i_use_svd_               = isvd
      this%bkg_base_rfmnt_          = bkgrfmt
      this%bkg_area_extension_      = bkgaext
      this%gen_peak_area_extension_ = genpaext
      this%max_area_extension_      = maxaext
      this%i_full_coverage_         = ifcov
      this%i_dump_modal_            = idumpmod

#ifdef __BSA_DEBUG
      write(unit_debug_, '(1x, a, a)') &
         INFOMSG, '@SettingsImpl::SetMshrSetts() : setting Bsa Mesher settings -- ok.'
#endif
   end subroutine SetMshrSetts

end submodule