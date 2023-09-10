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
submodule(BsaLib_WindData) BsaLib_WindDataImpl

#include "../precisions"

   use BsaLib_IO, only: INFOMSG, WARNMSG, ERRMSG, MSGCONT, DBGMSG &
                        , BSA_WIND_DATA_DUMPFILE
   use BsaLib_Data, only: bsa_Abort
   implicit none

contains



   module subroutine SetWindvertProf(this, ivaru)
      class(WindData_t), intent(inout) :: this
      integer(kind = 4), intent(in) :: ivaru

      if (ivaru < 1 .or. ivaru > 3) call bsa_Abort('Invalid "ivaru" value.')
      
      this%i_wind_prof_ = ivaru

#ifdef __BSA_DEBUG
      write(unit_debug_, '(1x, a, a, a)') &
         INFOMSG, '@WindImpl::SetWindvertProf() : wind vert prof set to ', &
         trim(CST_WIND_V_PROFILES(this%i_wind_prof_))
#endif
   end subroutine SetWindvertProf




   module subroutine SetMainvertDir(this, ivert)
      class(WindData_t), intent(inout) :: this
      integer(kind = 4), intent(in)    :: ivert

      if (ivert < 1 .or. ivert > 3) call bsa_Abort('Invalid "ivert" value.')
      this%i_vert_ = ivert

#ifdef __BSA_DEBUG
      write(unit_debug_, '(1x, a, a, i0)') INFOMSG, '@WindImpl::SetMainvertDir() : wind vert direction set to  ', this%i_vert_
#endif
   end subroutine SetMainvertDir




   module subroutine SetWindZoneLimits(this, lim, ilim)
      class(WindData_t), intent(inout) :: this
      real(RDP), intent(in), target    :: lim(..)
      integer(kind = 4), intent(in), optional :: ilim(..)


      select rank (lim)
      rank (0)

         if (.not. present(ilim)) call bsa_Abort('Please provide "ilim" value.')

         select rank (ilim)
         rank (0)
            if (ilim < 1 .or. ilim > this%nz_ + 1) call bsa_Abort('Invalid "ilim" value.')
            this%limits_wz_(ilim) = lim
         rank default
            call bsa_Abort('"ilim" must be an integer scalar value.')
         endselect

      rank (1)

         if (present(ilim)) then

            select rank (ilim)
            rank (1)
               if (.not. (size(lim(:)) == size(ilim(:)))) &
                  call bsa_Abort('sizes of "lim" and "ilim" do not match.')
                  
               this%limits_wz_(ilim(:)) = lim(:)
            rank default
               call bsa_Abort('expecting "ilim" to be a 1-rank array.')
            endselect

         else

            if (.not. this%nz_ == 0) then
               if (.not. (size(lim(:)) == this%nz_ + 1)) & 
                  call bsa_Abort('size of "lim" does not match number of wind zones.')
            else
               this%nz_ = size(lim(:)) - 1
            endif
            this%limits_wz_ => lim(:)
         endif ! present ilim
         
         
      rank default
         call bsa_Abort('Expeting "lim" either to be a scalar or a 1-rank array.')
      endselect
   end subroutine SetWindZoneLimits





   module subroutine SetAirDensity(aird)
      real(RDP), intent(in) :: aird

      if (aird < 0._RDP) call bsa_Abort('Air density has a negative value.')
      air_dens_ = aird

#ifdef __BSA_DEBUG
      write(unit_debug_, '(1x, a, a)') INFOMSG, '@WindImpl::SetAirDensity() : air density set -- ok.'
#endif
   end subroutine SetAirDensity




   module subroutine SetGlobalW2G(mat)
      real(RDP), intent(in) :: mat(3, 3)
      integer(kind = 4)    :: istat
      character(len = 256) :: emsg

      if (.not. allocated(rot_W2G_)) then
         allocate(rot_W2G_(3, 3), stat=istat, errmsg=emsg)
         if (istat == 0) then
#ifdef __BSA_ALLOC_DEBUG
            call allocOKMsg('rot_W2G_', [3, 3], loc(rot_W2G_), sizeof(rot_W2G_))
#endif
         else
            call allocKOMsg('rot_W2G_', istat, emsg)
         endif
      endif
      rot_W2G_ = mat

#ifdef __BSA_DEBUG
      write(unit_debug_, '(1x, a, a)') INFOMSG, '@WindImpl::SetGlobalW2G() : global rotation matrix WIND-GLOB set -- ok.'
#endif
   end subroutine SetGlobalW2G




   module subroutine SetWZMeanWindVel(this, UBref)
      class(WindData_t), intent(inout) :: this
      real(RDP), target, intent(in)    :: UBref(this%nz_)

      this%u_mean_ref_wz_ => UBref

#ifdef __BSA_DEBUG
      write(unit_debug_, '(1x, a, a)') &
         INFOMSG, &
         '@WindImpl::SetWZMeanWindVel() : wind zone mean wind speeds (at reference altitude) set -- ok.'
#endif
   end subroutine SetWZMeanWindVel




   module subroutine SetWZRefAlt(this, Zref)
      class(WindData_t), intent(inout) :: this
      real(RDP), target, intent(in)    :: Zref(this%nz_)

      this%Zref_wz_ => Zref

#ifdef __BSA_DEBUG
      write(unit_debug_, '(1x, a, a)') INFOMSG, '@WindImpl::SetWZRefAlt() : wind zones reference altitudes set -- ok.'
#endif
   end subroutine SetWZRefAlt




   module subroutine SetTurbWindScales(this, L)
      class(WindData_t), intent(inout) :: this
      real(RDP), target, intent(in)    :: L(3, 3, this%nz_)

      this%turb_scales_wz_ => L
   end subroutine SetTurbWindScales




   module subroutine SetTurbWindSDT(this, wtstd)
      class(WindData_t), intent(inout) :: this
      real(RDP), target, intent(in)    :: wtstd(3, this%nz_)

      this%sigmaUVW_wz_ => wtstd
   end subroutine SetTurbWindSDT




   module subroutine SetWindCorrCoeffs(this, corrcoeff)
      class(WindData_t), intent(inout) :: this
      real(RDP), target, intent(in)    :: corrcoeff(3, 3, this%nz_)

      this%corrCoeffs_wz_ => corrcoeff
   end subroutine SetWindCorrCoeffs




   module subroutine SetWindCorrExpnts(this, correxpn)
      class(WindData_t), intent(inout) :: this
      real(RDP), target, intent(in)    :: correxpn(3, 3, this%nz_)

      this%corrExp_wz_ => correxpn
   end subroutine SetWindCorrExpnts





   module subroutine SetIncidenceAngles(this, incang)
      class(WindData_t), intent(inout) :: this
      real(RDP), target, intent(in)    :: incang(this%nz_)

      this%incAng_wz_ => incang
   end subroutine SetIncidenceAngles





   module subroutine SetLocalRotMatW2G(this, rotW2G_L)
      class(WindData_t), intent(inout) :: this
      real(RDP), target, intent(in)    :: rotW2G_L(3, 3, this%nz_)

      this%rot_LW2G_wz_ => rotW2G_L
   end subroutine SetLocalRotMatW2G




   module subroutine setWindDirections(this, dirs, ndirs)
      class(WindData_t) :: this
      integer(kind = 4), intent(in) :: dirs(:)
      integer(kind = 4), intent(in), optional :: ndirs
      integer(kind = 4) :: itmp

      itmp = size(dirs)

      if (present(ndirs)) then
         if (itmp /= ndirs) then
            print '(1x, a, a/)', &
               ERRMSG, 'Size mismatch in setting spatial directions.'
            call bsa_Abort()
         endif
      endif

      this%i_ndirs_ = itmp
      this%dirs_    = dirs(1 : itmp)
   end subroutine


   
   module subroutine setTurbComps(this, tc, ntc)
      class(WindData_t) :: this
      integer(kind = 4), intent(in) :: tc(:)
      integer(kind = 4), intent(in), optional :: ntc
      integer(kind = 4) :: itmp

      itmp = size(tc)

      if (present(ntc)) then
         if (itmp /= ntc) then
            print '(1x, a, a/)', &
               ERRMSG, 'Size mismatch in setting turbulent components.'
            call bsa_Abort()
         endif
      endif

      this%i_ntc_ = itmp
      this%tc_    = tc(1 : itmp)
   end subroutine





   module subroutine SetTurbCompsAndDirsDefault(this)
      class(WindData_t), intent(inout) :: this
      integer :: itmp
      integer(kind = 4)    :: istat
      character(len = 256) :: emsg

      this%i_ntc_ = 1
      if (.not. allocated(this%tc_)) then
         allocate(this%tc_(1), stat=istat, errmsg=emsg)
         if (istat == 0) then
#ifdef __BSA_ALLOC_DEBUG
            call allocOKMsg('this % tc_', 1, loc(this%tc_), sizeof(this%tc_))
#endif
         else
            call allocKOMsg('this % tc_', istat, emsg)
         endif
      else
         itmp = size(this%tc_)
         if (itmp /= 1) then
            deallocate(this%tc_)
            allocate(this%tc_(1), stat=istat, errmsg=emsg)
            if (istat == 0) then
#ifdef __BSA_ALLOC_DEBUG
               call allocOKMsg('this % tc_', 1, loc(this%tc_), sizeof(this%tc_))
#endif
            else
               call allocKOMsg('this % tc_', istat, emsg)
            endif
         endif
      endif
      this%tc_(1) = 1

      this%i_ndirs_ = 1
      if (.not. allocated(this%dirs_)) then
         allocate(this%dirs_(1))
      else
         itmp = size(this%dirs_)
         if (itmp /= 1) then
            deallocate(this%dirs_)
            allocate(this%dirs_(1))
         endif
      endif
      this%dirs_(1) = 1

#ifdef __BSA_DEBUG
      write(unit_debug_, '(1x, a, a)') &
         INFOMSG, &
         '@WindImpl::SetTurbCompsAndDirsDefault() : default direction and turbulent components set -- ok.'
#endif
   end subroutine SetTurbCompsAndDirsDefault





   module subroutine SetNodalVel(this, Unod)
      class(WindData_t), intent(inout) :: this
      real(RDP), target, intent(in) :: Unod(:)

! #ifdef __BSA_DEBUG
!       write(unit_debug_, '(1x, a, a)') INFOMSG, '@WindImpl::SetNodalVel() : setting nodal velocities...'
! #endif

      this%u_node_ => Unod

#ifdef __BSA_DEBUG
      write(unit_debug_, '(1x, a, a)') INFOMSG, '@WindImpl::SetNodalVel() : setting nodal velocities -- ok.'
#endif
   end subroutine SetNodalVel




   module subroutine SetNodalWindZones(this, NodWZ)
      class(WindData_t), intent(inout) :: this
      integer(kind = 4), target, intent(in) :: NodWZ(:)

! #ifdef __BSA_DEBUG
!       write(unit_debug_, '(1x, a, a)') INFOMSG, '@WindImpl::SetNodalWindZones() : setting nodal wind zones...'
! #endif

      this%wz_node_ => NodWZ

#ifdef __BSA_DEBUG
      write(unit_debug_, '(1x, a, a)') INFOMSG, '@WindImpl::SetNodalWindZones() : setting nodal wind zones -- ok.'
#endif
   end subroutine SetNodalWindZones





   module subroutine SetNodalWindAltitudes(this, WnodAlt)
      class(WindData_t), intent(inout) :: this
      real(RDP), target, intent(in) :: WnodAlt(:)

! #ifdef __BSA_DEBUG
!       write(unit_debug_, '(1x, a, a)') INFOMSG, '@WindImpl::SetNodalWindAltitudes() : setting nodal wind altitudes...'
! #endif

      this%wAlt_node_ => WnodAlt

#ifdef __BSA_DEBUG
      write(unit_debug_, '(1x, a, a)') INFOMSG, '@WindImpl::SetNodalWindAltitudes() : setting nodal wind altitudes -- ok.'
#endif
   end subroutine SetNodalWindAltitudes






   module subroutine SetSpatialNodalCorr(this, nodCorr)
      class(WindData_t), intent(inout) :: this
      real(RDP), target, intent(in) :: nodCorr(:, :)

! #ifdef __BSA_DEBUG
!       write(unit_debug_, '(1x, a, a)') INFOMSG, '@WindImpl::SetSpatialNodalCorr() : setting spatial nodal correlation...'
! #endif

      this%nod_corr_ => nodCorr

#ifdef __BSA_DEBUG
      write(unit_debug_, '(1x, a, a)') INFOMSG, '@WindImpl::SetSpatialNodalCorr() : setting spatial nodal correlation -- ok.'
#endif
   end subroutine SetSpatialNodalCorr



   module subroutine getFull2DNodCorrMat(this, nn, nodcorr2d)
      use BsaLib_Utility, only: util_getCorrVectIndex
      class(WindData_t), intent(in) :: this
      integer(kind = 4), intent(in) :: nn
      real(RDP), allocatable, intent(inout) :: nodcorr2d(:, :)
      integer :: i_ = 0, j_, id_

      if (.not. associated(this%nod_corr_)) return

      if (.not. allocated(nodcorr2d)) allocate(nodcorr2d(nn, nn), stat=i_)
      if (i_ /= 0) return
      do j_ = 1, nn
         do i_ = 1, nn
            id_ = util_getCorrVectIndex(i_, j_, nn)
            nodcorr2d(i_, j_) = this%nod_corr_(id_, 1)
         enddo
      enddo
   end subroutine




   module subroutine SetWindFCoeffs(this, wfc)
      class(WindData_t), intent(inout) :: this
      real(RDP), target, intent(in) :: wfc(:, :, :)

! #ifdef __BSA_DEBUG
!       write(unit_debug_, '(1x, a, a)') INFOMSG, '@WindImpl::SetWindFCoeffs() : setting wind forces coefficients...'
! #endif

      this%wfc_ => wfc

#ifdef __BSA_DEBUG
      write(unit_debug_, '(1x, a, a)') INFOMSG, '@WindImpl::SetWindFCoeffs() : setting wind forces coefficients -- ok.'
#endif
   end subroutine SetWindFCoeffs





   module subroutine SetPhitimesC(this, phiTc)
      class(WindData_t), intent(inout) :: this
      real(RDP), target, intent(in)    :: phiTc(:, :, :)


      this%phi_times_A_ndegw_ => phiTc

#ifdef __BSA_DEBUG
      write(unit_debug_, '(1x, a, a)') INFOMSG, '@WindImpl::SetPhitimesC() : setting projected wind forces coefficients -- ok.'
#endif
   end subroutine SetPhitimesC







   module subroutine clean(this)
      class(WindData_t) :: this

#ifdef __BSA_DEBUG
      write(unit_debug_, '(1x, a, a)') INFOMSG, '@WindImpl::clean() : cleaning up...'
#endif

      if (allocated(rot_W2G_)) deallocate(rot_W2G_)

      if (associated(this%u_node_))            nullify(this%u_node_)
      if (associated(this%wz_node_))           nullify(this%wz_node_)
      if (associated(this%nod_corr_))          nullify(this%nod_corr_)
      if (associated(this%wfc_))               nullify(this%wfc_)
      if (associated(this%phi_times_A_ndegw_)) nullify(this%phi_times_A_ndegw_)

      if (associated(this%sigmaUVW_wz_))    nullify(this%sigmaUVW_wz_)
      if (associated(this%turb_scales_wz_)) nullify(this%turb_scales_wz_)
      if (associated(this%corrCoeffs_wz_))  nullify(this%corrCoeffs_wz_)
      if (associated(this%corrExp_wz_))     nullify(this%corrExp_wz_)
      if (associated(this%u_mean_ref_wz_))  nullify(this%u_mean_ref_wz_)
      if (associated(this%Zref_wz_))        nullify(this%Zref_wz_)
      if (associated(this%rot_LW2G_wz_))    nullify(this%rot_LW2G_wz_)
      if (associated(this%limits_wz_))      nullify(this%limits_wz_)
      if (associated(this%incAng_wz_))      nullify(this%incAng_wz_)

      if (allocated(this%psd_)) deallocate(this%psd_)

#ifdef __BSA_DEBUG
      write(unit_debug_, '(1x, a, a)') INFOMSG, '@WindImpl::clean() : cleaning up -- ok.'
#endif
   end subroutine clean






end submodule