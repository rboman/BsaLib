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
module BsaLib_Functions

#include "../../precisions"
   
   use BsaLib_Data, only: wd, struct_data, settings, dimM_bisp_, dimM_psd_
   implicit none
   public


   private :: wd, struct_data, settings, dimM_bisp_, dimM_psd_


   ! make a local internal copy
   integer(kind = 4) :: NFREQS, NNODES, NNODESL, NLIBS, NLIBSL
   integer(kind = 4) :: NMODES, NMODES_EFF
   integer(kind = 4), allocatable :: MODES(:)
   integer(kind = 4) :: NPSDEL, NTCOMPS, NDIRS = 1
   integer(kind = 4), allocatable :: TCOMPS(:), DIRS(:)

   integer(kind = 8)              :: MSHR_SVD_LWORK = - 1
   integer(kind = 8), allocatable :: MSHR_SVD_INFO
   double precision, allocatable  :: MSHR_SVD_WORK(:)
   

   interface

      module subroutine setBsaFunctionLocalVars()
      end subroutine


      module subroutine prefetchSVDWorkDim_()
      end subroutine

      module subroutine cleanSVDWorkInfo_()
      end subroutine



      module function getFM_full_tnm_scalar_msh_(fi, fj) result(bfm)
         real(RDP), intent(in) :: fi, fj
         real(RDP) :: bfm(dimM_bisp_)
      end function


      module function getFM_full_tm_scalar_msh_POD_(fi, fj) result(bfm)
         real(RDP), intent(in) :: fi, fj
         real(RDP) :: bfm(dimM_bisp_)
      end function


      module function getRM_full_scalar_msh_(bfm, fi, fj) result(brm)
         real(RDP), intent(in) :: bfm(dimM_bisp_), fi, fj
         real(RDP) :: brm(dimM_bisp_)
      end function



      module function getFM_diag_tnm_scalar_msh_(fi, fj) result(bfm)
         real(RDP), intent(in) :: fi, fj
         real(RDP) :: bfm(dimM_bisp_)
      end function


      module function getRM_diag_scalar_msh_(bfm, fi, fj) result(brm)
         real(RDP), intent(in) :: bfm(dimM_bisp_), fi, fj
         real(RDP) :: brm(dimM_bisp_)
      end function









      !> BUG: this routine is adapted to the case where we use
      !>      convention on PULSATION.
      !>      Please, adpapt it to the case of convention over FREQUENCIES.
      module subroutine getFM_full_tnlm_vect_cls_(f, Suvw, psd, bisp)
         real(RDP), intent(in)         :: f(NFREQS)
         real(RDP), intent(in)         :: Suvw(NFREQS, NPSDEL)
         real(RDP), allocatable, intent(inout) :: psd(:, :), bisp(:, :, :)
      end subroutine


      module subroutine getFM_full_tnm_vect_cls_(f, Suvw, psd, bisp)
         real(RDP), intent(in) :: f(NFREQS)
         real(RDP), intent(in) :: Suvw(NFREQS, NPSDEL)
         real(RDP), allocatable, intent(inout) :: psd(:, :), bisp(:, :, :)
      end subroutine



      module subroutine getRM_full_vect_cls_(f, psd, bisp)
         real(RDP), intent(in)                 :: f(NFREQS)
         real(RDP), allocatable, intent(inout) :: psd(:, :), bisp(:, :, :)
      end subroutine



      module subroutine getFM_diag_tnlm_vect_cls_(f, Suvw, psd, bisp)
         real(RDP), intent(in) :: f(NFREQS)
         real(RDP), intent(in) :: Suvw(NFREQS, NPSDEL)
         real(RDP), intent(inout), allocatable :: psd(:, :), bisp(:, :, :)
      end subroutine   



      module subroutine getRM_diag_vect_cls_(f, psd, bisp)
         real(RDP), intent(in)                 :: f(NFREQS)
         real(RDP), allocatable, intent(inout) :: psd(:, :), bisp(:, :, :)
      end subroutine







      !> BUG: this routine is adapted to the case where we use
      !>      convention on PULSATION.
      !>      Please, adapt it to the case of convention over FREQUENCIES.
      module pure subroutine getFM_full_tnlm_scalar_cls_(ii, ij, fi, fj, Suvw, Suvw_pad, psd, bisp)
         integer, intent(in)   :: ii, ij
         real(RDP), intent(in) :: fi, fj
         real(RDP), intent(in) :: Suvw(NFREQS, NPSDEL)
         real(RDP), intent(in) :: Suvw_pad(NPSDEL)
         real(RDP), intent(inout) :: psd(dimM_psd_), bisp(dimM_bisp_)
      end subroutine



      module pure subroutine getFM_full_tnm_scalar_cls_(ii, ij, fi, fj, Suvw, Suvw_pad, psd, bisp)
         integer, intent(in)      :: ii, ij
         real(RDP), intent(in)    :: fi, fj
         real(RDP), intent(in)    :: Suvw(NFREQS, NPSDEL)
         real(RDP), intent(in)    :: Suvw_pad(NPSDEL)
         real(RDP), intent(inout) :: psd(dimM_psd_), bisp(dimM_bisp_)
      end subroutine



      module subroutine getRM_full_scalar_cls_(ii, ij, fi, fj, psdin, psdout, bispin, bispout)
         integer, intent(in)    :: ii, ij
         real(RDP), intent(in)  :: fi, fj
         real(RDP), intent(in)  :: psdin(dimM_psd_), bispin(dimM_bisp_)
         real(RDP), intent(out) :: psdout(dimM_psd_), bispout(dimM_bisp_)
      end subroutine


      !> BUG: this routine is adapted to the case where we use
      !>      convention on PULSATION.
      !>      Please, adapt it to the case of convention over FREQUENCIES.
      module pure subroutine getFM_diag_tnlm_scalar_cls_(ii, ij, fi, fj, Suvw, Suvw_pad, psd, bisp)
         integer, intent(in)   :: ii, ij
         real(RDP), intent(in) :: fi, fj
         real(RDP), intent(in) :: Suvw(NFREQS, NPSDEL)
         real(RDP), intent(in) :: Suvw_pad(NPSDEL)
         real(RDP), intent(inout) :: psd(dimM_psd_), bisp(dimM_bisp_)
      end subroutine



      module subroutine getRM_diag_scalar_cls_(ii, ij, fi, fj, psdin, psdout, bispin, bispout)
         integer, intent(in)    :: ii, ij  ! freqs indexes
         real(RDP), intent(in)  :: fi, fj
         real(RDP), intent(in)  :: psdin(dimM_psd_), bispin(dimM_bisp_)
         real(RDP), intent(out) :: psdout(dimM_psd_), bispout(dimM_bisp_)
      end subroutine

      


      module pure subroutine getBR_SFm_val_(nm, Suvw, fnat, im, m, psd)
         integer, intent(in)      :: im, m, nm
         real(RDP), intent(in)    :: Suvw(nm, NPSDEL), fnat
         real(RDP), intent(inout) :: psd
      end subroutine

   end interface

end module BsaLib_Functions