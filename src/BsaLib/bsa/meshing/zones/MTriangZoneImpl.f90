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
submodule(BsaLib_MTriangZone) BsaLib_MTriangZoneImpl

! #ifndef BSA_M3MF_ONLY_PREMESH_
! # define BSA_M3MF_ONLY_PREMESH_ 0
! #else
! # if (BSA_M3MF_ONLY_PREMESH_ != 0 && BSA_M3MF_ONLY_PREMESH_ != 1)
! #  undef BSA_M3MF_ONLY_PREMESH_
! #  define BSA_M3MF_ONLY_PREMESH_ 0
! # endif
! #endif

   use BsaLib_Data, only: bsa_Abort
   use BsaLib_IO, only: INFOMSG, WARNMSG, ERRMSG, MSGCONT, DBGMSG, NOTEMSG
   implicit none


contains



   !> Default MTriangZone basic constructor
   module function MTriangZone_constr_def(name) result(this)
      character(len = *), intent(in), optional :: name
      type(MTriangZone_t) :: this

      if (present(name)) call this%zoneName(name)
   end function




   !> Gets rect base along I-dir
   module elemental function baseI_triang(this) result(res)
      class(MTriangZone_t), intent(in) :: this
      real(bsa_real_t) :: res

      res = getPointsDistance(this%Cpt_, this%Bpt_)
   end function


   !> Gets rect base along J-dir
   module elemental function baseJ_triang(this) result(res)
      class(MTriangZone_t), intent(in) :: this
      real(bsa_real_t) :: res

      res = getPointsDistance(this%Cpt_, this%Apt_)
   end function




   module subroutine deduceDeltas_triang(this)
      class(MTriangZone_t), intent(inout) :: this

      ! NOTE: do nothing
      !       DONE JUST BECAUSE WE NEED TO IMPLEMENT THIS
      !       METHOD SINCE IT IS AN ABSTRACT COMING FROM 
      !       PARENT CLASS.
   end subroutine





   !> Returns total N of zone's meshing points.
   !> NOTE: for the moment, ONLY triangles rectangles supported
   !> NOTE: overrides default inherited from parent class.
   module pure function zoneTotNPts_triang(this) result(npt)
      class(MTriangZone_t), intent(in) :: this
      integer :: npt

      npt = this%ni_ * this%ni_

      ! BUG: might cast to real and back to int ??
      npt = (npt + this%ni_ ) / 2
   end function





   !> Computes PAB angle, using "cosine" rule.
   !> NOTE: this only works if Cpt is the point at which
   !>       the wider angle is located!
   module subroutine setPABangle(this)
      class(MTriangZone_t), intent(inout) :: this

      real(bsa_real_t) :: CA, AB, CB, cang

      CA = getPointsDistance(this%Cpt_, this%Apt_)
      AB = getPointsDistance(this%Apt_, this%Bpt_)
      CB = getPointsDistance(this%Cpt_, this%Bpt_)

      cang = (CA*CA + AB*AB - CB*CB) / (2._bsa_real_t * CA * AB)
      this%PABang_ = acos(cang)
   end subroutine




   !> Used if we want some specific deltas
   !> along I and J directions
   module subroutine adaptToDeltas(this, dfi, dfj)
      class(MTriangZone_t), intent(inout) :: this
      real(bsa_real_t), value :: dfi, dfj

      real(bsa_real_t) :: CA, CB
      integer   :: ni, nj

      CA = this%baseJ()
      CB = this%baseI()

      ! n. of segments
      ni = floor(CB / dfi)
      nj = floor(CA / dfj)

      ! BUG: is this necessary??
      ! readapt deltas to fit in FIXED sides (bases)
      dfi = CB / ni
      dfj = CA / nj

      ! now get actual n. of points (refinements)
      this%ni_ = ni + 1
      this%nj_ = nj + 1
   end subroutine





   !> BUG: this routine might hide a bug.
   !>
   !> Deduce triangle zone rotation, knowing the 
   !> coordinates of the 3 triang points.
   !>
   !> Still, triangle rotation is given by the rotation
   !> needed to be done to reach the CA side, starting from
   !> Y positive axis of the GRS.
   module subroutine deduceRotation(this)
      class(MTriangZone_t), intent(inout) :: this

      real(bsa_real_t) :: dx, dy, ang

      ! BUG: not always 0 when supposed to.
      !      machine precision + finite floating point repr instabilities
      ! NOTE: do not call distance methods 
      !       since WE WANT TO KEEP SIGNS. 
      dx = this%Apt_%freqI() - this%Cpt_%freqI()
      dy = this%Apt_%freqJ() - this%Cpt_%freqJ()

      ang = atan(abs(dx) / abs(dy))

      if (.not. (dx >= 0._bsa_real_t .and. dy >= 0._bsa_real_t)) then ! NOT 1st quadrant

         if (dy < 0._bsa_real_t) then

            if (dx >= 0._bsa_real_t) then ! 2nd quadrant

               ang = CST_PIGREC - ang
            else ! 3rd quadrant

               ang = CST_PIGREC + ang
            endif

         else ! 4th quadrant

            ! BUG: precision instability.
            !      Tru to mitigate it, vary bad!!
            if (ang > 1e-14_bsa_real_t) then ! fixed tolerance acceptance

               ang = CST_PIt2 - ang
            else ! assume it 0

               ang = 0._bsa_real_t
            endif
         endif
      endif

      this%rot_ = ang
   end subroutine deduceRotation





   !> Gets actua points' delta increments based on this zone's
   !> rotation w.r.t GRS, when moving along the mesh.
   !> If NOT inverted:
   !>   - df_var : total delta along CA side (J-dir in LRS)
   !>   - df_cst : total delta along CB side (I-dir in LRS)
   module subroutine getRotatedUnaryDF(this, df_I_var, df_J_var, df_I_cst, df_J_cst)
      class(MTriangZone_t), intent(inout) :: this
      real(bsa_real_t), intent(out) :: df_I_var, df_J_var
      ! logical, intent(in), optional    :: invert
      real(bsa_real_t), intent(out), optional :: df_I_cst, df_J_cst

      logical :: do_invert = .false.
      real(bsa_real_t) :: ang

      
      if (present(df_I_cst) .and. present(df_J_cst)) do_invert = .true.


      if (this%rot_ <= CST_PId2) then ! 1st quadrant

         df_I_var = sin(this%rot_)
         df_J_var = cos(this%rot_)

      elseif (this%rot_ <= CST_PIGREC) then ! 2nd quadrant

         ang = this%rot_ - CST_PId2
         df_I_var =   cos(ang)
         df_J_var = - sin(ang)

      elseif (this%rot_ <= CST_PIt3d2) then ! 3rd quadrant

         ang = this%rot_ - CST_PIGREC
         df_I_var = - sin(ang)
         df_J_var = - cos(ang)

      elseif (this%rot_ <= CST_PIt2) then ! 4th quadrant

         ang = this%rot_ - CST_PIt3d2
         df_I_var = - cos(ang)
         df_J_var =   sin(ang)

      endif

      ! BUG: needed to avoid machine rounding precision
      if (abs(df_I_var) < MACHINE_PRECISION) df_I_var = 0._bsa_real_t
      if (abs(df_J_var) < MACHINE_PRECISION) df_J_var = 0._bsa_real_t


      if (do_invert) then ! BUG: ?????????

         df_I_cst =   df_J_var
         df_J_cst = - df_I_var

         ! BUG: needed to avoid machine rounding precision
         if (abs(df_I_cst) < MACHINE_PRECISION) df_I_cst = 0._bsa_real_t
         if (abs(df_J_cst) < MACHINE_PRECISION) df_J_cst = 0._bsa_real_t
      endif

   end subroutine getRotatedUnaryDF









   !> Defines a triang zone from 3 triang points. 
   !> No need to specify which one is A or B point, this is 
   !> automatically deduced.
   !> Automatically choses refinements definition.
   module subroutine defineFromPts_norm(this, Cp, P1, P2)
      class(MTriangZone_t), intent(inout) :: this
      class(MPoint_t), intent(in)         :: Cp, P1, P2

      if (this%ni_ == 0 .or. this%nj_ == 0) & 
         call bsa_Abort('Missing refinement information. Aborting.')


      block
         integer(int32) :: ni, nj
         character(len = *), parameter :: msg_segm = &
            ERRMSG//'Triangle collapsed into a segment. Aborting.'

         ni = this%ni_
         nj = this%nj_


         ! Try to deduce which one is A and B points.
         ! NOTE: this relies on the fact that we know
         !       a priori which is the corner point.
         if (P1%freqJ() > P2%freqJ() .or. P1%freqJ() < P2%freqJ()) then ! NE, E, SE || SW, W, NW

            this%Apt_ = P1
            this%Bpt_ = P2

         else ! same FJ

            if (P1%freqJ() > Cp%freqJ()) then ! w2 > 0

               if (P1%freqI() < P2%freqI()) then

                  this%Apt_ = P1
                  this%Bpt_ = P2
               
               elseif (P2%freqI() < P1%freqI()) then

                  this%Apt_ = P2
                  this%Bpt_ = P1
               else
                  call bsa_Abort(msg_segm)
               endif

            elseif (P1%freqJ() < Cp%freqJ()) then ! w2 < 0

               if (P1%freqI() > P2%freqI()) then

                  this%Apt_ = P1
                  this%Bpt_ = P2

               elseif (P2%freqI() > P1%freqI()) then

                  this%Apt_ = P2
                  this%Bpt_ = P1
               else
                  call bsa_Abort(msg_segm)
               endif

            else
               call bsa_Abort('Three points are aligned on the same line. Aborting.')
            endif
         endif

      end block


      ! NOTE: make a copy
      this%Cpt_ = MPoint(Cp)

      call this%setPABangle()
      call this%deduceRotation()

   end subroutine defineFromPts_norm









   !> Defines a triang zone from 3 triang points. 
   !> No need to specify which one is A or B point, this is 
   !> automatically deduced.
   !> Also, possible to specify either a value.
   !> If "is_refinement", then is either refinement along I or J direction.
   !> NOTE: this is the defult if no optional is passed.
   !>
   !> If "is_refinement" is false, but still a value is passed, 
   !> it is interpreted as a delta value, still along either I or J direction.
   module subroutine defineFromPts_refmtsORdeltas(&
         this, Cp, P1, P2, is_refinement, val_types, val1, val2)
      class(MTriangZone_t), intent(inout) :: this
      class(MPoint_t), intent(in)    :: Cp, P1, P2
      logical, intent(in)            :: is_refinement
      character(len = *), intent(in), optional :: val_types
      real(bsa_real_t), intent(in), optional   :: val1, val2

      integer(int32)   :: ni, nj
      real(bsa_real_t) :: dfi, dfj


      if (is_refinement) then

         if (present(val_types) .and. present(val1) .and. present(val2)) then ! is present, acquire

            block
               integer :: slen, i

               slen = len(val_types)
               if (slen == 0 .or. slen > 2) &
                  call bsa_Abort('"val_types" is either empty or too long.')

               do i = 1, slen

                  if (val_types(i:i) == 'i') then

                     if (this%ni_ == 0) then ! we need to set it
                        if (i == 1) then
                           if (.not.(present(val1))) &
                              call bsa_Abort('Missing base information along I-dir')
                           ni = val1
                        else ! == 2
                           if (.not.(present(val2))) &
                              call bsa_Abort('Missing base information along I-dir')
                           ni = val2
                        endif
                        this%ni_ = ni
                     endif

                  elseif (val_types(i:i) == 'j') then

                     if (this%nj_ == 0) then ! we need to set it
                        if (i == 1) then
                           if (.not.(present(val1))) &
                              call bsa_Abort('Missing base information along J-dir')
                           nj = val1
                        else ! == 2
                           if (.not.(present(val2))) &
                              call bsa_Abort('Missing base information along J-dir')
                           nj = val2
                        endif
                        this%nj_ = nj
                     endif

                  endif
               enddo

            end block

         endif ! optional vals present         


      else ! use deltas



         if (.not. (present(val_types) .and. present(val1) .and. present(val2))) &
            call bsa_Abort('Please provide both deltas.')
            
         
         block
            integer :: slen, i

            slen = len(val_types)
            if (slen == 0 .or. slen > 2) &
               call bsa_Abort('"val_types" is either empty or too long.')

            do i = 1, slen

               if (val_types(i:i) == 'i') then

                  if (i == 1) then
                     if (.not.(present(val1))) &
                        call bsa_Abort('Missing delta information along I-dir')
                     dfi = val1
                  else ! == 2
                     if (.not.(present(val2))) &
                        call bsa_Abort('Missing delta information along I-dir')
                     dfi = val2
                  endif

               elseif (val_types(i:i) == 'j') then

                  if (i == 1) then
                     if (.not.(present(val1))) &
                        call bsa_Abort('Missing delta information along J-dir')
                     dfj = val1
                  else ! == 2
                     if (.not.(present(val2))) &
                        call bsa_Abort('Missing delta information along J-dir')
                     dfj = val2
                  endif

               endif
            enddo
         end block

      endif ! use refinements / deltas

      ! Common definition.
      call this%defineFromPts_norm(Cp, P1, P2)

      if (.not. is_refinement) call this%adaptToDeltas(dfi, dfj)
   end subroutine defineFromPts_refmtsORdeltas








   !> Computes a completely defined triang zone
   module subroutine compute(this)
      class(MTriangZone_t), intent(inout) :: this

      real(bsa_real_t) :: df_CA, df_CB


      ! get (full) deltas in straight directions
      df_CA = getPointsDistance(this%Cpt_, this%Apt_) ! base J
      df_CB = getPointsDistance(this%Cpt_, this%Bpt_) ! base I


      ! BUG: Allow for different triangle sides (NOT isosceles)

      ! BUG: it might be that they're not equal 
      !      when they should. So, no direct comparison
      !      instead use MACHINE PRECISION
      if (abs(df_CA - df_CB) <= MACHINE_PRECISION) then
         call computeISOTriangle(this)
      else
         call bsa_Abort('Triangle must be "isosceles". Aborting.')
      endif
   end subroutine compute




   subroutine computeISOTriangle(this)
      use BsaLib_Data, only: &
         dimM_bisp_, getBFM_msh, settings  &
         , m3mf_msh_ptr_, msh_NZones, msh_bfmpts_pre_
      class(MTriangZone_t), intent(inout) :: this

      if (this%ni_ /= this%nj_) &
         call bsa_Abort('Different sides'' refinements not yet allowed. Aborting.')

      block
         integer(int32)   :: Np, Np_m1, Nj_m1
         real(bsa_real_t) :: df_CA, df_CST, dfMAJi, dfMAJj, dfMINi, dfMINj

         integer(int32)   :: id, i, j
         real(bsa_real_t) :: base_fi, base_fj, fi, fj

         integer(int32) :: tot
         real(bsa_real_t), allocatable :: bfm(:, :)

#ifdef BSA_M3MF_ONLY_PREMESH_
         real(bsa_real_t) :: dw, ctr_infl, brd_infl, vtx_infl_rect, vtx_infl_triang
         real(bsa_real_t), allocatable :: intg(:)
#endif

         type(MPoint_t) :: pt

         ! caching
         Np    = this%ni_

         ! tot = (Np*Np + Np) / 2
         tot   = Np*Np
         tot   = tot + Np
         tot   = tot / 2
         
         Np_m1 = Np - 1
         Nj_m1 = Np_m1


         ! get unary delta freq increments in both I and J directions
         ! get (full) deltas in straight directions
         df_CA  = getPointsDistance(this%Cpt_, this%Apt_) ! base J
         df_CST = df_CA / Np_m1
         call this%getRotatedUnaryDF(dfMAJi, dfMAJj, dfMINi, dfMINj)

         ! get actualised values
         dfMAJi = dfMAJi * df_CST
         dfMAJj = dfMAJj * df_CST
         dfMINi = dfMINi * df_CST
         dfMINj = dfMINj * df_CST

#ifdef BSA_M3MF_ONLY_PREMESH_
         ! influences for integration
         dw = df_CST * CST_PIt2
         ctr_infl = dw * dw
         brd_infl = ctr_infl / 2._bsa_real_t
         vtx_infl_rect   = brd_infl / 2._bsa_real_t
         vtx_infl_triang = vtx_infl_rect / 2._bsa_real_t
         
         allocate(intg(dimM_bisp_))
#endif
         allocate(bfm(dimM_bisp_, tot))

         
         !=======================================================
         !
         ! FIRST COLUMN (CA side)
         !

         ! track base points
         base_fi = this%Cpt_%freqI()
         base_fj = this%Cpt_%freqJ()

         ! first line
         fi = base_fi
         fj = base_fj
         bfm(:, 1) = getBFM_msh(fi, fj)
#ifdef BSA_M3MF_ONLY_PREMESH_
         intg = bfm(:, 1) * vtx_infl_rect
#endif

         ! internal lines
         ! NOTE: since is beginning, we can use j directly
         !       to index into bfm.
         do j = 2, Nj_m1

            fi = fi + dfMAJi
            fj = fj + dfMAJj
            bfm(:, j) = getBFM_msh(fi, fj)
#ifdef BSA_M3MF_ONLY_PREMESH_
            intg = intg + bfm(:, j) * brd_infl
#endif
         enddo

         ! last line
         fi = fi + dfMAJi
         fj = fj + dfMAJj
         pt = MPoint(fi, fj)

         if (.not. pt == this%Apt_) &
            call bsa_Abort('First end point does not match triang A point. Aborting.')

         id = j
         bfm(:, id) = getBFM_msh(fi, fj)
#ifdef BSA_M3MF_ONLY_PREMESH_
         intg = intg + bfm(:, id) * vtx_infl_triang
#endif

         
         
         !=======================================================
         !
         ! INTERNAL COLUMNS (Parallel to CA)
         !
         do i = 2, Np_m1

            Nj_m1 = Nj_m1 - 1
            id    = id + 1

            base_fi = base_fi + dfMINi
            base_fj = base_fj + dfMINj

            ! first line
            fi = base_fi
            fj = base_fj
            bfm(:, id) = getBFM_msh(fi, fj)
#ifdef BSA_M3MF_ONLY_PREMESH_
            intg = intg + bfm(:, id) * brd_infl
#endif

            ! internal lines
            do j = 2, Nj_m1

               id = id + 1
               fi = fi + dfMAJi
               fj = fj + dfMAJj
               bfm(:, id) = getBFM_msh(fi, fj)
#ifdef BSA_M3MF_ONLY_PREMESH_
               intg = intg + bfm(:, id) * ctr_infl
#endif
            enddo

            ! last line
            id = id + 1
            fi = fi + dfMAJi
            fj = fj + dfMAJj
            bfm(:, id) = getBFM_msh(fi, fj)
#ifdef BSA_M3MF_ONLY_PREMESH_
            intg = intg + bfm(:, id) * brd_infl
#endif
         enddo ! i


         !=======================================================
         !
         ! LAST COLUMN (B point !!)
         !
         base_fi = base_fi + dfMINi
         base_fj = base_fj + dfMINj

         pt = MPoint(base_fi, base_fj)

         if (.not. pt == this%Bpt_) &
            call bsa_Abort('End point does not match B triang point. Aborting.')

         id = id + 1
         bfm(:, id) = getBFM_msh(base_fi, base_fj)
#ifdef BSA_M3MF_ONLY_PREMESH_
         intg = intg + bfm(:, id) * vtx_infl_triang
#endif


! #ifdef __BSA_DEBUG
         if (.not. id == tot) &
            call bsa_Abort('"id" does not equal tot N of Triang zone''s points.')
! #endif

         !$omp critical
#ifdef BSA_M3MF_ONLY_PREMESH_
         m3mf_msh_ptr_   = m3mf_msh_ptr_ + (intg * settings%i_bisp_sym_) ! update main integral
#endif
         msh_NZones      = msh_NZones + 1                ! update n. of zones count
         msh_bfmpts_pre_ = msh_bfmpts_pre_ + tot         ! update tot num of meshing points

         ! eventually, update zone with max N of points
         if (tot > msh_max_zone_NPts) msh_max_zone_NPts = tot

         ! dump zone info
         call DumpZone(this, bfm)
         !$omp end critical

      end block

! #ifdef __BSA_DEBUG
!       write(unit_debug_, *) ' @MTriangZoneImpl::computeISOTriangle() : init -- ok.'
! #endif
   end subroutine computeISOTriangle







   !> Dumps a triang zone data for later reconstruction
   module subroutine dumpTZ(this)
      class(MTriangZone_t), intent(in) :: this

      write(unit_dump_bfm_) MZone_ID%TRIANGLE

      ! 3 pts
      write(unit_dump_bfm_) this%Cpt_%freqI(), this%Cpt_%freqJ()
      write(unit_dump_bfm_) this%Apt_%freqI(), this%Apt_%freqJ()
      write(unit_dump_bfm_) this%Bpt_%freqI(), this%Bpt_%freqJ()
      
      ! NOTE: useless, since rot might reconstructed from points
      write(unit_dump_bfm_) this%rot_
      write(unit_dump_bfm_) this%ni_, this%nj_

#ifdef __BSA_ZONE_DEBUG
      write(unit=4533, fmt=*) &
         'Refms at  TZ=', trim(this%name_), this%ni_, this%nj_, &
         'thread id= ', omp_get_thread_num()
#endif
   end subroutine dumpTZ




   !> Undumps a triang zone data for later reconstruction
   module subroutine undumpTZ(this)
      class(MTriangZone_t), intent(inout) :: this
      real(bsa_real_t) :: rval1, rval2

      ! 3 pts
      read(unit_dump_bfm_) rval1, rval2
      call this%Cpt_%setFreqs(rval1, rval2)

      read(unit_dump_bfm_) rval1, rval2
      call this%Apt_%setFreqs(rval1, rval2)

      read(unit_dump_bfm_) rval1, rval2
      call this%Bpt_%setFreqs(rval1, rval2)


      read(unit_dump_bfm_) this%rot_
      read(unit_dump_bfm_) this%ni_, this%nj_
   end subroutine undumpTZ




   !> Implementation of triang zone interpolation methods
   module subroutine interpolateTZ( this &
#ifdef __BSA_OMP
      , bfm, pdata &
#endif
      & )
      class(MTriangZone_t), intent(inout) :: this
#ifdef __BSA_OMP
      real(bsa_real_t), intent(in)  :: bfm(:, :)
      class(*), pointer, intent(in) :: pdata

      ! NOTE: for the moment only supporting HTPC method
      call interpolateTZ_HTPC_v3(this, bfm, pdata)
#else
      call interpolateTZ_HTPC_v3(this)
#endif
   end subroutine





   !> Queries Triang zone n of points if it had ni/nj refinements
   elemental function getTriangZoneEquivNPts(ni, nj) result(npt)
      integer(int32), intent(in) :: ni, nj
      integer(int32) :: npt

      ! NOTE: for the moment ni==nj
      npt = ni * nj

      ! BUG: do we really need to divide by integer??
      npt = int((npt + ni) / 2._real32, kind=int32)
   end function




   !> Implementation of HTPC interpolation method for a rectangle zone,
   !> including MultiLevel-Refinement for BFM data.
   subroutine interpolateTZ_HTPC_v3( this &
#ifdef __BSA_OMP
      , bfm_undump, pdata &
#endif
      & )
      use BsaLib_Data, only: &
#ifndef __BSA_OMP
         bfm_undump, &
#endif
         dimM_bisp_, getBFM_msh, getBRM_msh, m3mf_msh_ptr_, m3mr_msh_ptr_, settings  &
         , msh_bfmpts_post_, msh_brmpts_post_  &
         , write_brm_fptr_, do_export_brm_, BrmExportBaseData_t

      class(MTriangZone_t), intent(inout) :: this
#ifdef __BSA_OMP
      real(bsa_real_t), intent(in)  :: bfm_undump(:, :)
      class(*), pointer, intent(in) :: pdata
#endif

      integer(int32)   :: interp_fact, njOld, njNew_piprev, njNew_picurr, njNew_tmp, njtmp
      integer(int32)   :: ni, nj, ni_bfm_ref_, nj_bfm_ref_
      integer(int32)   :: nipI, nipJ, ipos, nPtsPost
      integer(int32)   :: i, ist, n_segs_bfm_ref_i_, n_segs_bfm_ref_j_
      integer(int32)   :: n_pts_bfm_ref_i_, n_pts_bfm_ref_j_
      real(bsa_real_t) :: dfIi_bfm_lv0_, dfIj_bfm_lv0_, dfJi_bfm_lv0_, dfJj_bfm_lv0_
      real(bsa_real_t) :: dfIi_bfm_ref_, dfIj_bfm_ref_, dfJi_bfm_ref_, dfJj_bfm_ref_
      real(bsa_real_t) :: dfI_bfm_lv0_, dfJ_bfm_lv0_, dfI_bfm_ref_, dfJ_bfm_ref_, df_cst
      real(bsa_real_t) :: dfIi_brm_interp, dfIj_brm_interp, dfJi_brm_interp, dfJj_brm_interp
      real(bsa_real_t) :: dfI_brm_interp_, dfJ_brm_interp_, dw_cst, df_diag_old, df_diag_bfm_ref, df_diag_brm_interp
      real(bsa_real_t) :: vtx_infl_r, brd_infl_r, ctr_infl_r, brd_infl_t, vtx_infl_t

      ! HTPC indexes
      integer(int32) :: pIcurr, pIprev, pJhead, pJtail

      ! Pos in general BFM undumped data
      integer(int32) :: i_bfm_old, i_bfm_ref_i, i_bfm_ref_j
      integer(int32) :: i_brm, i_brm_shift, i_brm_offsetJ
#ifndef __BSA_OMP
      integer(int32) :: i_brm_write_
#endif
      integer(int32) :: i_bfm_interpJ

      ! freqs
      real(bsa_real_t) :: fi_baseptI, fj_baseptI, fi, fj, fi_baseptJ, fj_baseptJ
#ifdef __BSA_OMP
      real(bsa_real_t), allocatable, dimension(:) :: fi_v_, fj_v_
#endif

      real(bsa_real_t), allocatable :: bfm_new_left(:, :), bfm_new_right(:, :)
      real(bsa_real_t), allocatable :: bfm_interp(:, :)

      real(bsa_real_t) :: dfJtail, dfJhead, dfIcurr, dfIprev, dfJ_oldtmp
      real(bsa_real_t) :: bfmtail(dimM_bisp_), bfmhead(dimM_bisp_)

      real(bsa_real_t), allocatable :: brm(:, :)
#ifndef BSA_M3MF_ONLY_PREMESH_
      real(bsa_real_t) :: vtx_infl_r_bfm, brd_infl_r_bfm, ctr_infl_r_bfm, brd_infl_t_bfm, vtx_infl_t_bfm
      real(bsa_real_t) :: intg_bfm(dimM_bisp_)
#endif
      real(bsa_real_t) :: intg(dimM_bisp_)


      ! BUG: for the moment, we take the max
      !      such to ensure having the SAME n. of points
      !      along both sides.
      !      Later, consider supporting more general approach.
      interp_fact = max(this%policy_%interp_I_fct_, this%policy_%interp_J_fct_)

      ! get unary delta freqs increments in I and J directions
      ! (old BFM values, to be interpolated)
      njOld  = this%nj_ - 1
      df_cst = getPointsDistance(this%Cpt_, this%Apt_) / njOld

      ! NOTE: first two refer to MAJOUR direction (J), then MINOR (I)
      call this%getRotatedUnaryDF(dfJi_bfm_lv0_, dfJj_bfm_lv0_, dfIi_bfm_lv0_, dfIj_bfm_lv0_)

      ! actualise them (scaled w.r.t actual sides' length)
      dfIi_bfm_lv0_ = dfIi_bfm_lv0_ * df_cst
      dfIj_bfm_lv0_ = dfIj_bfm_lv0_ * df_cst
      dfJi_bfm_lv0_ = dfJi_bfm_lv0_ * df_cst
      dfJj_bfm_lv0_ = dfJj_bfm_lv0_ * df_cst

      dfI_bfm_lv0_ = df_cst
      dfJ_bfm_lv0_ = df_cst

      ! get n. of BFM refinement segments (between two old ones)
      n_segs_bfm_ref_i_ = 1 ! original
      n_segs_bfm_ref_j_ = 1
      do i = 1, this%policy_%n_interp_bfm_lvs_
         n_segs_bfm_ref_i_ = n_segs_bfm_ref_i_ * this%policy_%interp_bfm_I_fct_
         n_segs_bfm_ref_j_ = n_segs_bfm_ref_j_ * this%policy_%interp_bfm_J_fct_
      enddo
      n_pts_bfm_ref_i_ = n_segs_bfm_ref_i_ - 1
      n_pts_bfm_ref_j_ = n_segs_bfm_ref_j_ - 1
      
      ! get refined (BFM) deltas (GRS)
      dfIi_bfm_ref_ = dfIi_bfm_lv0_ / n_segs_bfm_ref_i_
      dfIj_bfm_ref_ = dfIj_bfm_lv0_ / n_segs_bfm_ref_i_
      dfJi_bfm_ref_ = dfJi_bfm_lv0_ / n_segs_bfm_ref_j_
      dfJj_bfm_ref_ = dfJj_bfm_lv0_ / n_segs_bfm_ref_j_


      ! get BRM interpolated deltas (GRS) (taken from refined BFM deltas this time)
      dfIi_brm_interp = dfIi_bfm_ref_ / interp_fact
      dfIj_brm_interp = dfIj_bfm_ref_ / interp_fact
      dfJi_brm_interp = dfJi_bfm_ref_ / interp_fact
      dfJj_brm_interp = dfJj_bfm_ref_ / interp_fact


      ! get absolute deltas (in LRS), along I and J directions
      dfI_bfm_ref_    = dfI_bfm_lv0_ / n_segs_bfm_ref_i_
      dfI_brm_interp_ = dfI_bfm_ref_ / interp_fact

      dfJ_bfm_ref_    = dfJ_bfm_lv0_ / n_segs_bfm_ref_j_
      dfJ_brm_interp_ = dfJ_bfm_ref_ / interp_fact


#ifndef BSA_M3MF_ONLY_PREMESH_
      ! influence areas for BRM integration
      dw_cst         = dfI_bfm_ref_ * CST_PIt2
      ctr_infl_r_bfm = dw_cst * dw_cst
      brd_infl_r_bfm = ctr_infl_r_bfm / 2._bsa_real_t
      vtx_infl_r_bfm = brd_infl_r_bfm / 2._bsa_real_t
      brd_infl_t_bfm = brd_infl_r_bfm
      vtx_infl_t_bfm = vtx_infl_r_bfm / 2._bsa_real_t
#endif

      ! influence areas for BRM integration
      dw_cst     = dfI_brm_interp_ * CST_PIt2
      ctr_infl_r = dw_cst * dw_cst
      brd_infl_r = ctr_infl_r / 2._bsa_real_t
      vtx_infl_r = brd_infl_r / 2._bsa_real_t
      brd_infl_t = brd_infl_r
      vtx_infl_t = vtx_infl_r / 2._bsa_real_t

      ! get actualised BFM-refined and BRM-interp  refinements (along borders)
      ni_bfm_ref_ = (this%ni_ - 1)
      nj_bfm_ref_ = (this%nj_ - 1)

      ! take a backup since this will be decremented by one each time we move pi_curr
      njOld = nj_bfm_ref_

      ni          = ni_bfm_ref_ * (n_segs_bfm_ref_i_ * interp_fact) + 1
      nj          = nj_bfm_ref_ * (n_segs_bfm_ref_j_ * interp_fact) + 1
      ni_bfm_ref_ = ni_bfm_ref_ * n_segs_bfm_ref_i_ + 1
      nj_bfm_ref_ = nj_bfm_ref_ * n_segs_bfm_ref_j_ + 1
      
      njNew_picurr = nj

      ! number of BRM points to interpolate (insert)
      ! between two know BFM (refined) points' direction lines.
      nipI = interp_fact - 1
      nipJ = nipI

      i_brm_offsetJ = nipI * nj

      ! TODO: deltas along the hypotenuse
      !       Last section when moving pj_head
      df_diag_old        = sqrt(dfI_bfm_lv0_**2 + dfJ_bfm_lv0_**2)
      df_diag_bfm_ref    = df_diag_old / n_segs_bfm_ref_i_
      df_diag_brm_interp = df_diag_bfm_ref / interp_fact

      ! allocate data      
      nPtsPost = getTriangZoneEquivNPts(ni, nj)
      allocate(brm(dimM_bisp_, nPtsPost), stat=ist)
      if (ist /= 0) call bsa_Abort("Error allocating ""brm"" in interpolating RZ.")
      brm = 0._bsa_real_t

      allocate(bfm_new_left(dimM_bisp_, nj), stat=ist)
      if (ist /= 0) call bsa_Abort("Error allocating ""bfm_new_left"" in interpolating RZ.")
      bfm_new_left  = 0._bsa_real_t

      allocate(bfm_new_right(dimM_bisp_, nj - nipJ - 1), stat=ist)
      if (ist /= 0) call bsa_Abort("Error allocating ""bfm_new_right"" in interpolating RZ.")
      bfm_new_right = 0._bsa_real_t
      
      allocate(bfm_interp(dimM_bisp_, nj), stat=ist)
      if (ist /= 0) call bsa_Abort("Error allocating ""bfm_interp"" in interpolating RZ.")
      bfm_interp = 0._bsa_real_t


#ifdef __BSA_OMP
      allocate(fi_v_(nPtsPost), stat=ist)
      if (ist /= 0) call bsa_Abort("Error allocating ""fi"" in interpolating RZ.")
      fi = 0._bsa_real_t

      allocate(fj_v_(nPtsPost), stat=ist)
      if (ist /= 0) call bsa_Abort("Error allocating ""fj"" in interpolating RZ.")
      fj = 0._bsa_real_t

# define __FREQ_I_  fi_v_(i_brm) = fi
# define __FREQ_J_  fj_v_(i_brm) = fj
# define __FREQ_I_shift_  fi_v_(i_brm_shift) = fi
# define __FREQ_J_shift_  fj_v_(i_brm_shift) = fj
# define __PDATA ,pdata
      if (do_export_brm_ .and. associated(pdata)) then
         select type (pdata)
            class is (BrmExportBaseData_t)
               pdata%nI_ = ni
               pdata%nJ_ = nj
         end select
      endif

#else

! # define __FREQ_I_
! # define __FREQ_J_
# define __FREQ_I_shift_
# define __FREQ_J_shift_
# define __PDATA ,null()

#endif


      ! before starting, interpolate very first column
      ! along J dir, to get new mesh from old
      ! Then for all the others, it will be done inside
      ! the loop over the NI (old) points of the old mesh.
      ! NOTE: integrate as well.
      
      ! init point vertex
      fi_baseptI = this%Cpt_%freqI()
      fj_baseptI = this%Cpt_%freqJ()
      fi         = fi_baseptI
      fj         = fj_baseptI
      fi_baseptJ = fi
      fj_baseptJ = fj

      bfmtail            = bfm_undump(:, 1)
#ifndef BSA_M3MF_ONLY_PREMESH_
      intg_bfm           = bfmtail * vtx_infl_r_bfm
#endif
      bfm_new_left(:, 1) = bfmtail
      
      brm(:, 1) = getBRM_msh(bfmtail, fi, fj)
#ifdef __BSA_OMP
      __FREQ_I_
      __FREQ_J_
#else
      call write_brm_fptr_(fi, fj, brm(:, 1)  __PDATA)
#endif
      intg      = brm(:, 1) * vtx_infl_r

      i_brm = 1
      do pJhead = 2, this%nj_ ! loop on all OLD BFM saved points (J-dir)

         do i_bfm_ref_j = 1, n_pts_bfm_ref_j_ ! loop on all REF BFM pts between 2 old.
            
            ! compute head
            fi_baseptJ = fi_baseptJ + dfJi_bfm_ref_
            fj_baseptJ = fj_baseptJ + dfJj_bfm_ref_
            bfmhead    = getBFM_msh(fi_baseptJ, fj_baseptJ)
#ifndef BSA_M3MF_ONLY_PREMESH_
            intg_bfm   = intg_bfm + bfmhead * brd_infl_r_bfm
#endif

            dfJhead = dfJ_bfm_ref_
            dfJtail = 0._bsa_real_t
            do pJtail = 1, nipJ ! interp (J-dir) between tail-head

               fi = fi + dfJi_brm_interp
               fj = fj + dfJj_brm_interp

               ! update actual distances head/tail
               dfJtail = dfJtail + dfJ_brm_interp_
               dfJhead = dfJhead - dfJ_brm_interp_

               ! interpolation along J dir (between HEAD-TAIL)
               ! NOTE: save BFM for later use.
               i_brm                  = i_brm + 1
               bfm_new_left(:, i_brm) = (bfmhead * dfJtail + bfmtail * dfJhead) / dfJ_bfm_ref_
               brm(:, i_brm)          = getBRM_msh(bfm_new_left(:, i_brm), fi, fj)
#ifdef __BSA_OMP
               __FREQ_I_
               __FREQ_J_
#else
               call write_brm_fptr_(fi, fj, brm(:, i_brm), null())
#endif
               intg = intg + brm(:, i_brm) * brd_infl_r ! NOTE: it is a border point
            enddo

            ! treat head (new BFM refined point)
            fi = fi + dfJi_brm_interp
            fj = fj + dfJj_brm_interp

#ifdef __BSA_DEBUG
            ! DEBUG: they should equate fi/fj_baseptI!!
            if (abs(fi - fi_baseptJ) > MACHINE_PRECISION .or. &
               abs(fj - fj_baseptJ) > MACHINE_PRECISION) &
                  call bsa_Abort(&
                     'BAD (1): fi or fj at the end of a BFM ref segment does not coincide..')
#endif

            i_brm = i_brm + 1
            bfm_new_left(:, i_brm) = bfmhead
            brm(:, i_brm) = getBRM_msh(bfm_new_left(:, i_brm), fi, fj)
#ifdef __BSA_OMP
            __FREQ_I_
            __FREQ_J_
#else
            call write_brm_fptr_(fi, fj, brm(:, i_brm), null())
#endif
            ! NOTE: it is a border point, except for the very last one (VERTEX)
            intg = intg + brm(:, i_brm) * brd_infl_r

            bfmtail = bfmhead
         enddo ! n. of (exact) ref points for BFM

         !
         ! NOTE: now new head is OLD BFM (next) point!
         !
         bfmhead  = bfm_undump(:, pJhead) ! OK because we stored it NJ majour.
#ifndef BSA_M3MF_ONLY_PREMESH_
         intg_bfm = intg_bfm + bfmhead * brd_infl_r_bfm
#endif
         
#ifdef __BSA_DEBUG
         ! DEBUG: 
         if (n_pts_bfm_ref_j_ > 0) then 
            if (abs((fi_baseptI + (dfJi_bfm_lv0_*(pJhead-1))) - (fi_baseptJ + dfJi_bfm_ref_)) > MACHINE_PRECISION .or. &
                  abs((fj_baseptI + (dfJj_bfm_lv0_*(pJhead-1))) - (fj_baseptJ + dfJj_bfm_ref_)) > MACHINE_PRECISION) &
                     call bsa_Abort(&
                        'BAD (2-t): fi or fj at the end of a BFM ref segment does not coincide..')
         endif
#endif
         
         dfJhead = dfJ_bfm_ref_
         dfJtail = 0._bsa_real_t
         do pJtail = 1, nipJ ! interp (J-dir) between tail-head

            fi = fi + dfJi_brm_interp
            fj = fj + dfJj_brm_interp

            ! update actual distances head/tail
            dfJtail = dfJtail + dfJ_brm_interp_
            dfJhead = dfJhead - dfJ_brm_interp_

            ! interpolation along J dir (between HEAD-TAIL)
            ! NOTE: save it for later use.
            i_brm                  = i_brm + 1
            bfm_new_left(:, i_brm) = (bfmhead * dfJtail + bfmtail * dfJhead) / dfJ_bfm_ref_
            brm(:, i_brm)          = getBRM_msh(bfm_new_left(:, i_brm), fi, fj)
#ifdef __BSA_OMP
            __FREQ_I_
            __FREQ_J_
#else
            call write_brm_fptr_(fi, fj, brm(:, i_brm), null())
#endif
            intg = intg + brm(:, i_brm) * brd_infl_r ! NOTE: it is a border point
         enddo ! pJtail = 1, nipJ

         ! here, treat head, TAIL==HEAD (head - old mesh)
         fi    = fi + dfJi_brm_interp
         fj    = fj + dfJj_brm_interp
         i_brm                  = i_brm + 1
         bfm_new_left(:, i_brm) = bfmhead
         brm(:, i_brm)          = getBRM_msh(bfm_new_left(:, i_brm), fi, fj)
#ifdef __BSA_OMP
         __FREQ_I_
         __FREQ_J_
#else
         call write_brm_fptr_(fi, fj, brm(:, i_brm), null())
#endif

         ! NOTE: it is a border point, except for the very last one (VERTEX)
         intg = intg + brm(:, i_brm) * brd_infl_r

         ! old head (old mesh) becomes new tail
         ! TODO: we could change bfmhead here, so that
         !       it might be already ready for next loop.
         bfmtail = bfmhead

         ! NOTE: set new bfm ref base freqs as (current) head (old-mesh) point
         fi_baseptJ = fi
         fj_baseptJ = fj
      enddo ! pJhead = 2, this%nj_

      ! NOTE: removing excess contribution for last HEAD (VERTEX)
      intg = intg + brm(:, i_brm) * (vtx_infl_t - brd_infl_r)
#ifndef BSA_M3MF_ONLY_PREMESH_
      intg_bfm = intg_bfm - bfmtail * vtx_infl_t_bfm
#endif



      !
      i_bfm_old = this%nj_

      do pIcurr = 2, this%ni_ ! loop on all OLD BFM infl lines (I-dir)

         ! before doing any computation,
         ! we need to interpolate BFM along J
         ! at new CURRENT (I) infl line (including ref infl lines)
         ! NOTE: once we go through, integrate as well.

         do i_bfm_ref_i = 1, n_pts_bfm_ref_i_ ! loop on all REF BFM pts between 2 old (I-dir)

            ! computing BRM offset from pi_prev and pi_curr infl lines.
            i_brm_shift = i_brm
            njNew_tmp   = njNew_picurr
            do pJhead = 1, nipI
               njNew_tmp   = njNew_tmp - 1
               i_brm_shift = i_brm_shift + njNew_tmp
            enddo
#ifndef __BSA_OMP
            i_brm_write_ = i_brm_shift
#endif

            ! reset freqs to point to new base -> prev base moved by ref deltas!
            ! (now pi_prev-pj_tail)
            ! NOTE: still keep prev base in memory here since they might serve later.
            fi_baseptJ = fi_baseptI + dfIi_bfm_ref_  ! reset J bases to match next I
            fj_baseptJ = fj_baseptI + dfIj_bfm_ref_
            fi         = fi_baseptJ
            fj         = fj_baseptJ

            bfmtail             = getBFM_msh(fi, fj)
#ifndef BSA_M3MF_ONLY_PREMESH_
            intg_bfm            = intg_bfm + bfmtail * brd_infl_r_bfm
#endif
            bfm_new_right(:, 1) = bfmtail

            i_brm_shift         = i_brm_shift + 1
            brm(:, i_brm_shift) = getBRM_msh(bfm_new_right(:, 1), fi, fj)
            ! call write_brm_fptr_(fi, fj, brm(:, i_brm_shift), null())
            __FREQ_I_shift_
            __FREQ_J_shift_

            intg = intg + brm(:, i_brm_shift) * brd_infl_r

            !
            ! here is where we compute right infl-line (for I-dir interpolation)
            !
            i_bfm_interpJ = 1
            do pJhead = 2, njOld ! loop on all OLD BFM saved points (J-dir)

               do i_bfm_ref_j = 1, n_pts_bfm_ref_j_ ! loop on all REF BFM pts between 2 old.

                  fi_baseptJ = fi_baseptJ + dfJi_bfm_ref_
                  fj_baseptJ = fj_baseptJ + dfJj_bfm_ref_
                  bfmhead    = getBFM_msh(fi_baseptJ, fj_baseptJ)
#ifndef BSA_M3MF_ONLY_PREMESH_
                  intg_bfm   = intg_bfm + bfmhead * ctr_infl_r_bfm
#endif
         
                  ! once we moved head, restore init, distances from head/tail
                  dfJhead = dfJ_bfm_ref_
                  dfJtail = 0._bsa_real_t
                  do pJtail = 1, nipJ ! interp (J-dir) between tail-head
         
                     fi = fi + dfJi_brm_interp
                     fj = fj + dfJj_brm_interp
         
                     ! update actual distances head/tail
                     dfJtail = dfJtail + dfJ_brm_interp_
                     dfJhead = dfJhead - dfJ_brm_interp_
         
                     ! interpolation along J dir (between HEAD-TAIL)
                     i_bfm_interpJ = i_bfm_interpJ + 1
                     bfm_new_right(:, i_bfm_interpJ) = &
                        (bfmhead * dfJtail + bfmtail * dfJhead) / dfJ_bfm_ref_

                     i_brm_shift         = i_brm_shift + 1
                     brm(:, i_brm_shift) = getBRM_msh(bfm_new_right(:, i_bfm_interpJ), fi, fj)
                     ! call write_brm_fptr_(fi, fj, brm(:, i_brm_shift), null())
                     __FREQ_I_shift_
                     __FREQ_J_shift_

                     ! NOTE: it is a center point, except for very last row -> BORDER
                     intg  = intg + brm(:, i_brm_shift) * ctr_infl_r
                  enddo ! pJtail = 1, nipJ

                  ! tail in now head (new BFM refined point)
                  fi = fi + dfJi_brm_interp
                  fj = fj + dfJj_brm_interp

#ifdef __BSA_DEBUG
                  ! DEBUG:
                  if (abs(fi - fi_baseptJ) > MACHINE_PRECISION .or. &
                     abs(fj - fj_baseptJ) > MACHINE_PRECISION) &
                        call bsa_Abort(&
                           'BAD (4-t): fi or fj at the end of a BFM ref segment does not coincide..')
#endif

                  i_bfm_interpJ                   = i_bfm_interpJ + 1
                  bfm_new_right(:, i_bfm_interpJ) = bfmhead

                  i_brm_shift         = i_brm_shift + 1
                  brm(:, i_brm_shift) = getBRM_msh(bfm_new_right(:, i_bfm_interpJ), fi, fj)
                  ! call write_brm_fptr_(fi, fj, brm(:, i_brm_shift), null())
                  __FREQ_I_shift_
                  __FREQ_J_shift_

                  ! NOTE: it is a center point, except for the very last one (BORDER)
                  intg = intg + brm(:, i_brm_shift) * ctr_infl_r

                  bfmtail = bfmhead
               enddo ! n. of (exact) ref points for BFM
      
               ! here, treat head, TAIL==HEAD (old mesh, J-dir)
               ! NOTE: it is a center point, except for the very last one (BORDER)
               fi_baseptJ = fi_baseptJ + dfJi_bfm_ref_
               fj_baseptJ = fj_baseptJ + dfJj_bfm_ref_
               bfmhead    = getBFM_msh(fi_baseptJ, fj_baseptJ)
#ifndef BSA_M3MF_ONLY_PREMESH_
               intg_bfm   = intg_bfm + bfmhead * ctr_infl_r_bfm
#endif

               ! once we moved head, restore init, distances from head/tail
               dfJhead = dfJ_bfm_ref_
               dfJtail = 0._bsa_real_t
               do pJtail = 1, nipJ
      
                  fi = fi + dfJi_brm_interp
                  fj = fj + dfJj_brm_interp
      
                  ! update actual distances head/tail
                  dfJtail = dfJtail + dfJ_brm_interp_
                  dfJhead = dfJhead - dfJ_brm_interp_
      
                  ! interpolation along J dir (between HEAD-TAIL)
                  i_bfm_interpJ = i_bfm_interpJ + 1
                  bfm_new_right(:, i_bfm_interpJ) = &
                     (bfmhead * dfJtail + bfmtail * dfJhead) / dfJ_bfm_ref_

                  i_brm_shift         = i_brm_shift + 1
                  brm(:, i_brm_shift) = getBRM_msh(bfm_new_right(:, i_bfm_interpJ), fi, fj)
                  ! call write_brm_fptr_(fi, fj, brm(:, i_brm_shift), null())
                  __FREQ_I_shift_
                  __FREQ_J_shift_
      
                  ! NOTE: it is a center point, except for very last row -> BORDER
                  intg  = intg + brm(:, i_brm_shift) * ctr_infl_r
               enddo ! pJtail = 1, nipJ

               ! here treat this new head
               fi = fi + dfJi_brm_interp
               fj = fj + dfJj_brm_interp
               i_bfm_interpJ                   = i_bfm_interpJ + 1
               bfm_new_right(:, i_bfm_interpJ) = bfmhead

               i_brm_shift         = i_brm_shift + 1
               brm(:, i_brm_shift) = getBRM_msh(bfm_new_right(:, i_bfm_interpJ), fi, fj)
               ! call write_brm_fptr_(fi, fj, brm(:, i_brm_shift), null())
               __FREQ_I_shift_
               __FREQ_J_shift_

               intg = intg + brm(:, i_brm_shift) * ctr_infl_r

               bfmtail = bfmhead ! old head becomes new tail
            enddo ! pJhead = 2, njOld

            !
            ! Last segments
            ! I.e. next BFM head (computed) lies on the hypotenuse
            !
            do i_bfm_ref_j = 1, n_pts_bfm_ref_j_ - (i_bfm_ref_i - 1)

               fi_baseptJ = fi_baseptJ + dfJi_bfm_ref_
               fj_baseptJ = fj_baseptJ + dfJj_bfm_ref_
               bfmhead    = getBFM_msh(fi_baseptJ, fj_baseptJ)
#ifndef BSA_M3MF_ONLY_PREMESH_
               intg_bfm   = intg_bfm + bfmhead * ctr_infl_r_bfm
#endif

               dfJhead = dfJ_bfm_ref_
               dfJtail = 0._bsa_real_t
               do pJtail = 1, nipJ

                  fi = fi + dfJi_brm_interp
                  fj = fj + dfJj_brm_interp

                  dfJtail = dfJtail + dfJ_brm_interp_
                  dfJhead = dfJhead - dfJ_brm_interp_

                  i_bfm_interpJ                   = i_bfm_interpJ + 1
                  bfm_new_right(:, i_bfm_interpJ) = &
                        (bfmhead * dfJtail + bfmtail * dfJhead) / dfJ_bfm_ref_

                  i_brm_shift         = i_brm_shift + 1
                  brm(:, i_brm_shift) = getBRM_msh(bfm_new_right(:, i_bfm_interpJ), fi, fj)
                  ! call write_brm_fptr_(fi, fj, brm(:, i_brm_shift), null())
                  __FREQ_I_shift_
                  __FREQ_J_shift_

                  intg = intg + brm(:, i_brm_shift) * ctr_infl_r
               enddo

               ! here treat actual HEAD (last on diag line -> remove integral surplus)
               fi = fi + dfJi_brm_interp
               fj = fj + dfJj_brm_interp
               i_bfm_interpJ                   = i_bfm_interpJ + 1
               bfm_new_right(:, i_bfm_interpJ) = bfmhead

               i_brm_shift         = i_brm_shift + 1
               brm(:, i_brm_shift) = getBRM_msh(bfm_new_right(:, i_bfm_interpJ), fi, fj)
               ! call write_brm_fptr_(fi, fj, brm(:, i_brm_shift), null())
               __FREQ_I_shift_
               __FREQ_J_shift_

               intg = intg + brm(:, i_brm_shift) * ctr_infl_r

               bfmtail = bfmhead
            enddo

            ! removing excess of very last HEAD, accounted as center, it is BORDER.
            ! NOTE: even worse for very last HEAD which happens to be End point.
            !       there, it is a VERTEX point.
            !       However, it's the very last element in brm, we can remove it after.
            intg = intg - brm(:, i_brm_shift) * brd_infl_t
#ifndef BSA_M3MF_ONLY_PREMESH_
            intg_bfm = intg_bfm - bfmtail * brd_infl_t_bfm
#endif


            ! backup NJ new intrp n. of points for next iteration
            ! NOTE: also, refers to NJ point at pi_curr J inlf line.
            njNew_piprev = njNew_picurr
            njNew_picurr = i_bfm_interpJ ! NOTE: at last iter, should be 1


            !
            ! Now INTERPOLATE along I-dir between left and right BFM infl lines.
            !
            dfIcurr = dfI_bfm_ref_ ! reset I-dir CURR-PREV distances
            dfIprev = 0._bsa_real_t

            dfJ_oldtmp = dfJ_bfm_lv0_
            njtmp      = nipJ
            do pIprev = 1, nipI ! interp (I-dir) between prev-curr

               ! bulk I-dir interpolation until pj_head section level.
               ! Then, after treat that triang shaped zone separately.

               dfIprev = dfIprev + dfI_brm_interp_
               dfIcurr = dfIcurr - dfI_brm_interp_
               
               bfm_interp(:, 1 : njNew_picurr) = &
                  (  bfm_new_left (:, 1 : njNew_picurr) * dfIcurr + &
                     bfm_new_right(:, 1 : njNew_picurr) * dfIprev ) / dfI_bfm_ref_

               ! once we have the values, go through them to integrate
               ! NOTE: reset base freqs pointers, this time moving them along INTERP mesh
               fi = fi_baseptI + (dfIi_brm_interp * pIprev)
               fj = fj_baseptI + (dfIj_brm_interp * pIprev)

               i_brm         = i_brm + 1
               brm(:, i_brm) = getBRM_msh(bfm_interp(:, 1), fi, fj)
#ifdef __BSA_OMP
               __FREQ_I_
               __FREQ_J_
#else
               call write_brm_fptr_(fi, fj, brm(:, i_brm), null())
#endif
               intg = intg + brm(:, i_brm) * brd_infl_r
               
               do pJtail = 2, njNew_picurr

                  fi = fi + dfJi_brm_interp
                  fj = fj + dfJj_brm_interp
                  i_brm         = i_brm + 1
                  brm(:, i_brm) = getBRM_msh(bfm_interp(:, pJtail), fi, fj)
#ifdef __BSA_OMP
                  __FREQ_I_
                  __FREQ_J_
#else
                  call write_brm_fptr_(fi, fj, brm(:, i_brm), null())
#endif
                  intg = intg + brm(:, i_brm) * ctr_infl_r
               enddo

               !==============================================
               ! here we have to treat triang shaped zone.
               ! NOTE: tmp head is interpolated along hypotenuse!!
               bfmhead = (&
                  bfm_new_left (:, njNew_piprev) * (df_diag_brm_interp * pIprev)  + &
                  bfm_new_right(:, njNew_picurr) * (df_diag_brm_interp * (nipI + 1 - pIprev)) ) / df_diag_bfm_ref

               ! once we have the head, continue interpolating along J between tail and tmp head
               ! NOTE: everytime we move pi_prev, actual Nj points reduce by 1
               njtmp      = njtmp - 1
               dfJ_oldtmp = dfJ_oldtmp - dfJ_brm_interp_
               dfJhead    = dfJ_oldtmp
               dfJtail    = 0._bsa_real_t
               ipos       = njNew_picurr ! tail index position
               do pJtail = 1, njtmp

                  fi = fi + dfJi_brm_interp
                  fj = fj + dfJj_brm_interp

                  dfJtail = dfJtail + dfJ_brm_interp_
                  dfJhead = dfJhead - dfJ_brm_interp_

                  ! NOTE: tail is stored in bfm_interp(:, njNew_picurr)
                  ipos = ipos + 1
                  bfm_interp(:, ipos) = &
                     (bfmhead * dfJtail + dfJhead * bfm_interp(:, njNew_picurr)) / dfJ_oldtmp

                  i_brm         = i_brm + 1
                  brm(:, i_brm) = getBRM_msh(bfm_interp(:, ipos), fi, fj)
#ifdef __BSA_OMP
                  __FREQ_I_
                  __FREQ_J_
#else
                  call write_brm_fptr_(fi, fj, brm(:, i_brm), null())
#endif
                  intg = intg + brm(:, i_brm) * ctr_infl_r
               enddo

               ! treat HEAD (on hypotenuse) separately
               fi = fi + dfJi_brm_interp
               fj = fj + dfJj_brm_interp
               i_brm         = i_brm + 1
               brm(:, i_brm) = getBRM_msh(bfmhead, fi, fj)
#ifdef __BSA_OMP
               __FREQ_I_
               __FREQ_J_
#else
               call write_brm_fptr_(fi, fj, brm(:, i_brm), null())
#endif
               intg = intg + brm(:, i_brm) * vtx_infl_t
               !==============================================

            enddo ! pIprev = 1, nipI (interp points along I dir)

            ! now update bases along I (CURR now, PREV next iteration!)
            fi_baseptI = fi_baseptI + dfIi_bfm_ref_
            fj_baseptI = fj_baseptI + dfIj_bfm_ref_


#ifndef __BSA_OMP
            ! ! BUG: check this condition
            ! if (njNew_picurr /= ((njOld - 1) * (nipJ + 1) + 1)) &
            !    call bsa_Abort("""njNew_picurr"" does not match computed value.")
            
            fi = fi_baseptI
            fj = fj_baseptI
            do pIprev = 1, njNew_picurr
               i_brm_write_ = i_brm_write_ + 1
               call write_brm_fptr_(fi, fj, brm(:, i_brm_write_), null())
               fi = fi + dfJi_brm_interp
               fj = fj + dfJj_brm_interp
            enddo
#endif

            bfm_new_left(:, 1 : njNew_picurr) = bfm_new_right(:, 1 : njNew_picurr)

            ! once finished with this section (CURR-PREV), since we skip J column
            ! at pi_curr infl line, we reset general BRM index to point to
            ! previously shifted one.
            i_brm = i_brm_shift

            ! njOld = njOld - 1

         enddo ! BFM ref points along I dir



         !
         ! Here, right BFM infl line is one where we have old BFM mesh points !
         !

         ! computing BRM offset from pi_prev and pi_curr J infl lines.
         i_brm_shift = i_brm
         njNew_tmp   = njNew_picurr
         do pJhead = 1, nipI
            njNew_tmp   = njNew_tmp - 1
            i_brm_shift = i_brm_shift + njNew_tmp
         enddo
#ifndef __BSA_OMP
         i_brm_write_ = i_brm_shift
#endif

         ! again, I bases refer to PREV infl line.
         fi_baseptJ = fi_baseptI + dfIi_bfm_ref_
         fj_baseptJ = fj_baseptI + dfIj_bfm_ref_
         fi         = fi_baseptJ
         fj         = fj_baseptJ

         i_bfm_old           = i_bfm_old + 1
         bfmtail             = bfm_undump(:, i_bfm_old)
#ifndef BSA_M3MF_ONLY_PREMESH_
         intg_bfm            = intg_bfm + bfmtail * brd_infl_r_bfm
#endif
         bfm_new_right(:, 1) = bfmtail
         
         i_brm_shift         = i_brm_shift + 1
         brm(:, i_brm_shift) = getBRM_msh(bfm_new_right(:, 1), fi, fj)
         ! call write_brm_fptr_(fi, fj, brm(:, i_brm_shift), null())
         __FREQ_I_shift_
         __FREQ_J_shift_

         intg = intg + brm(:, i_brm_shift) * brd_infl_r

         ! compute right infl-line
         i_bfm_interpJ = 1
         do pJhead = 2, njOld

            do i_bfm_ref_j = 1, n_pts_bfm_ref_j_ ! loop on all REF BFM pts between 2 old.

               fi_baseptJ = fi_baseptJ + dfJi_bfm_ref_
               fj_baseptJ = fj_baseptJ + dfJj_bfm_ref_
               bfmhead    = getBFM_msh(fi_baseptJ, fj_baseptJ)
#ifndef BSA_M3MF_ONLY_PREMESH_
               intg_bfm   = intg_bfm + bfmhead * ctr_infl_r_bfm
#endif
      
               ! once we moved head, restore init, distances from head/tail
               dfJhead = dfJ_bfm_ref_
               dfJtail = 0._bsa_real_t
               do pJtail = 1, nipJ
      
                  fi = fi + dfJi_brm_interp
                  fj = fj + dfJj_brm_interp
      
                  ! update actual distances head/tail
                  dfJtail = dfJtail + dfJ_brm_interp_
                  dfJhead = dfJhead - dfJ_brm_interp_
      
                  ! interpolation along J dir (between HEAD-TAIL)
                  i_bfm_interpJ = i_bfm_interpJ + 1
                  bfm_new_right(:, i_bfm_interpJ) = &
                     (bfmhead * dfJtail + bfmtail * dfJhead) / dfJ_bfm_ref_

                  i_brm_shift         = i_brm_shift + 1
                  brm(:, i_brm_shift) = getBRM_msh(bfm_new_right(:, i_bfm_interpJ), fi, fj)
                  ! call write_brm_fptr_(fi, fj, brm(:, i_brm_shift), null())
                  __FREQ_I_shift_
                  __FREQ_J_shift_
      
                  ! NOTE: it is a center point, except for very last row -> BORDER
                  intg  = intg + brm(:, i_brm_shift) * ctr_infl_r
               enddo ! pJtail = 1, nipJ

               ! tail in now head (new BFM refined point)
               fi = fi + dfJi_brm_interp
               fj = fj + dfJj_brm_interp

               i_bfm_interpJ                   = i_bfm_interpJ + 1
               bfm_new_right(:, i_bfm_interpJ) = bfmhead

               i_brm_shift         = i_brm_shift + 1
               brm(:, i_brm_shift) = getBRM_msh(bfm_new_right(:, i_bfm_interpJ), fi, fj)
               ! call write_brm_fptr_(fi, fj, brm(:, i_brm_shift), null())
               __FREQ_I_shift_
               __FREQ_J_shift_

               ! NOTE: it is a center point, except for the very last one (BORDER)
               intg = intg + brm(:, i_brm_shift) * ctr_infl_r

               bfmtail = bfmhead
            enddo ! i_bfm_ref_j = 1, n_pts_bfm_ref_j_

            ! next head is OLD BFM point
            fi = fi_baseptJ
            fj = fj_baseptJ
            fi_baseptJ = fi_baseptJ + dfJi_bfm_ref_
            fj_baseptJ = fj_baseptJ + dfJj_bfm_ref_
            
            i_bfm_old = i_bfm_old + 1
            bfmhead   = bfm_undump(:, i_bfm_old)
#ifndef BSA_M3MF_ONLY_PREMESH_
            intg_bfm  = intg_bfm + bfmhead * ctr_infl_r_bfm
#endif
   
            ! once we moved head, restore init distances from head/tail
            dfJhead = dfJ_bfm_ref_
            dfJtail = 0._bsa_real_t
            do pJtail = 1, nipJ
   
               fi = fi + dfJi_brm_interp
               fj = fj + dfJj_brm_interp
   
               ! update actual distances head/tail
               dfJtail = dfJtail + dfJ_brm_interp_
               dfJhead = dfJhead - dfJ_brm_interp_
   
               ! interpolation along J dir (between HEAD-TAIL)
               ! NOTE: save it for later use.
               i_bfm_interpJ                   = i_bfm_interpJ + 1
               bfm_new_right(:, i_bfm_interpJ) = &
                  (bfmhead * dfJtail + bfmtail * dfJhead) / dfJ_bfm_ref_

               i_brm_shift         = i_brm_shift + 1
               brm(:, i_brm_shift) = getBRM_msh(bfm_new_right(:, i_bfm_interpJ), fi, fj)
               ! call write_brm_fptr_(fi, fj, brm(:, i_brm_shift), null())
               __FREQ_I_shift_
               __FREQ_J_shift_
   
               ! NOTE: it is a center point, except for very last row -> BORDER
               intg  = intg + brm(:, i_brm_shift) * ctr_infl_r
            enddo ! pJtail = 1, nipJ
   
            ! here, treat head, TAIL==HEAD
            fi = fi + dfJi_brm_interp
            fj = fj + dfJj_brm_interp
            i_bfm_interpJ                   = i_bfm_interpJ + 1
            bfm_new_right(:, i_bfm_interpJ) = bfmhead
            
            i_brm_shift         = i_brm_shift + 1
            brm(:, i_brm_shift) = getBRM_msh(bfm_new_right(:, i_bfm_interpJ), fi, fj)
            ! call write_brm_fptr_(fi, fj, brm(:, i_brm_shift), null())
            __FREQ_I_shift_
            __FREQ_J_shift_
   
            ! NOTE: it is a center point, except for the very last one (BORDER)
            intg = intg + brm(:, i_brm_shift) * ctr_infl_r
   
            ! old head becomes new tail
            bfmtail = bfmhead
         enddo ! pJhead
         
         ! removing excess of very last HEAD
         ! accounted as center, it is BORDER
         ! NOTE: even worse for very last HEAD which happens to be End point.
         !       there, it is a VERTEX point.
         !       However, it's the very last element in brm, we can remove it after.
         intg = intg - brm(:, i_brm_shift) * brd_infl_t
#ifndef BSA_M3MF_ONLY_PREMESH_
         intg_bfm = intg_bfm - bfmtail * brd_infl_t_bfm
#endif


         ! backup NJ new intrp n. of points for next iteration
         ! NOTE: also, refers to NJ point at pi_curr J inlf line.
         njNew_piprev = njNew_picurr
         njNew_picurr = i_bfm_interpJ ! NOTE: at last iter, should be 1




         ! ok, here we now have BFM values (interpolated along J)
         ! at CURR and PREV (I) index pointers.
         ! We have to interpolate along I between CURR and PREV, i.e.
         ! prev has to start moving toward CURR.

         dfIcurr    = dfI_bfm_ref_ ! reset I-dir CURR-PREV distances
         dfIprev    = 0._bsa_real_t
         dfJ_oldtmp = dfJ_bfm_ref_
         njtmp      = nipJ
         do pIprev = 1, nipI ! interpolate along I

            ! bulk I-dir interpolation until pj_head section level.
            ! Then, after treat that triang shaped zone separately.

            dfIprev = dfIprev + dfI_brm_interp_
            dfIcurr = dfIcurr - dfI_brm_interp_
            
            bfm_interp(:, 1 : njNew_picurr) = &
               (  bfm_new_left (:, 1 : njNew_picurr) * dfIcurr + &
                  bfm_new_right(:, 1 : njNew_picurr) * dfIprev ) / dfI_bfm_ref_

            ! once we have the values, go through them to integrate
            ! NOTE: reset base freqs pointers, this time moving them along INTERP mesh
            fi = fi_baseptI + (dfIi_brm_interp * pIprev)
            fj = fj_baseptI + (dfIj_brm_interp * pIprev)

            i_brm         = i_brm + 1
            brm(:, i_brm) = getBRM_msh(bfm_interp(:, 1), fi, fj)
#ifdef __BSA_OMP
            __FREQ_I_
            __FREQ_J_
#else
            call write_brm_fptr_(fi, fj, brm(:, i_brm), null())
#endif
            intg = intg + brm(:, i_brm) * brd_infl_r
            
            do pJtail = 2, njNew_picurr

               fi = fi + dfJi_brm_interp
               fj = fj + dfJj_brm_interp
               
               i_brm         = i_brm + 1
               brm(:, i_brm) = getBRM_msh(bfm_interp(:, pJtail), fi, fj)
#ifdef __BSA_OMP
               __FREQ_I_
               __FREQ_J_
#else
               call write_brm_fptr_(fi, fj, brm(:, i_brm), null())
#endif
               intg = intg + brm(:, i_brm) * ctr_infl_r
            enddo

            !==============================================
            ! here we have to treat triang shaped zone.
            ! NOTE: tmp head is interpolated along hypotenuse!!
            bfmhead = (&
               bfm_new_left (:, njNew_piprev) * (df_diag_brm_interp * pIprev)  + &
               bfm_new_right(:, njNew_picurr) * (df_diag_brm_interp * (nipI + 1 - pIprev)) ) / df_diag_bfm_ref

            ! once we have the head, continue interpolating along J between tail and tmp head
            ! NOTE: everytime we move pi_prev, actual Nj points reduce by 1
            njtmp      = njtmp - 1
            dfJ_oldtmp = dfJ_oldtmp - dfJ_brm_interp_
            dfJhead    = dfJ_oldtmp
            dfJtail    = 0._bsa_real_t
            ipos       = njNew_picurr ! tail index position
            do pJtail = 1, njtmp

               fi = fi + dfJi_brm_interp
               fj = fj + dfJj_brm_interp

               dfJtail = dfJtail + dfJ_brm_interp_
               dfJhead = dfJhead - dfJ_brm_interp_

               ! NOTE: tail is stored in bfm_interp(:, njNew_picurr)
               ipos = ipos + 1
               bfm_interp(:, ipos) = &
                  (bfmhead * dfJtail + dfJhead * bfm_interp(:, njNew_picurr)) / dfJ_oldtmp

               i_brm         = i_brm + 1
               brm(:, i_brm) = getBRM_msh(bfm_interp(:, ipos), fi, fj)
#ifdef __BSA_OMP
               __FREQ_I_
               __FREQ_J_
#else
               call write_brm_fptr_(fi, fj, brm(:, i_brm), null())
#endif
               intg = intg + brm(:, i_brm) * ctr_infl_r
            enddo

            ! treat HEAD (on hypotenuse) separately
            fi = fi + dfJi_brm_interp
            fj = fj + dfJj_brm_interp
            i_brm         = i_brm + 1
            brm(:, i_brm) = getBRM_msh(bfmhead, fi, fj)
#ifdef __BSA_OMP
            __FREQ_I_
            __FREQ_J_
#else
            call write_brm_fptr_(fi, fj, brm(:, i_brm), null())
#endif
            intg = intg + brm(:, i_brm) * vtx_infl_t
            !==============================================

         enddo ! pIprev = 1, nipI


         ! now we can update (PREV) base freqs to match CURR
         ! NOTE: CURR J infl line has already been computed!
         fi_baseptI = fi_baseptI + dfIi_bfm_ref_
         fj_baseptI = fj_baseptI + dfIj_bfm_ref_


#ifndef __BSA_OMP
         ! BUG: check this condition
         if (njNew_picurr /= ((njOld - 1) * (nipJ + 1) + 1)) &
            call bsa_Abort("""njNew_picurr"" does not match computed value.")
         
         fi = fi_baseptI
         fj = fj_baseptI
         do pIprev = 1, njNew_picurr
            i_brm_write_ = i_brm_write_ + 1
            call write_brm_fptr_(fi, fj, brm(:, i_brm_write_), null())
            fi = fi + dfJi_brm_interp
            fj = fj + dfJj_brm_interp
         enddo
#endif


         bfm_new_left(:, 1 : njNew_picurr) = bfm_new_right(:, 1 : njNew_picurr)


         ! once finished with this section (CURR-PREV), since we skip J column
         ! at pi_curr infl line, we reset general BRM index to point to
         ! previously shifted one.
         i_brm = i_brm_shift

         
         njOld = njOld - 1

      enddo ! pIcurr = 2, this%ni_


! #ifdef __BSA_DEBUG
      if (i_brm /= nPtsPost) &
         call bsa_Abort('"i_bfm_old" does not equal rect zone''s n. of (interpolated) points.')
! #endif


      ! removing overestimation for B point
      ! NOTE: should be stored at very last index!!
      intg = intg - brm(:, i_brm) * vtx_infl_t
#ifndef BSA_M3MF_ONLY_PREMESH_
      intg_bfm = intg_bfm - bfmtail * vtx_infl_t_bfm
#endif

      !$omp critical
#ifdef __BSA_OMP
      call write_brm_fptr_(fi_v_, fj_v_, brm, pdata)
#endif

      msh_bfmpts_post_ = msh_bfmpts_post_ + getTriangZoneEquivNPts(ni_bfm_ref_, nj_bfm_ref_)
      msh_brmpts_post_ = msh_brmpts_post_ + nPtsPost

#ifndef BSA_M3MF_ONLY_PREMESH_
      m3mf_msh_ptr_ = m3mf_msh_ptr_ + (intg_bfm * settings%i_bisp_sym_)   ! update main BFM integral
#endif
      m3mr_msh_ptr_ = m3mr_msh_ptr_ + (intg     * settings%i_bisp_sym_)   ! update main BRM integral
      !$omp end critical

   end subroutine interpolateTZ_HTPC_v3















!    subroutine interpolateTZ_HTPC_v2(this)
! #ifdef __BSA_OMP
!       use BsaLib_Data, only: dimM_bisp_, getBRM_msh, m3mr_msh_ptr_
! #else
!       use BsaLib_Data, only: dimM_bisp_, getBRM_msh, bfm_undump, m3mr_msh_ptr_
! #endif
!       class(MTriangZone_t), intent(inout) :: this

!       integer   :: interp_fact, ni, nj, njOld, njNew_piprev, njNew_picurr, njNew_tmp, njtmp
!       integer   :: nipI, nipJ, zNintrpPts, ipos
!       real(bsa_real_t) :: df_cst, dfIi_old, dfIj_old, dfJi_old, dfJj_old
!       real(bsa_real_t) :: dfIi_interp, dfIj_interp, dfJi_interp, dfJj_interp
!       real(bsa_real_t) :: dfI_old, dfJ_old, dfI_interp, dfJ_interp
!       real(bsa_real_t) :: df_diag_old, df_diag_interp
!       real(bsa_real_t) :: dw_cst
!       real(bsa_real_t) :: vtx_infl_r, brd_infl_r, ctr_infl_r, brd_infl_t, vtx_infl_t

!       ! HTPC indexes
!       integer :: picurr, piprev, pjhead, pjtail

!       !> Pos in general BFM undumped data
!       integer :: i_bfm

!       ! BRM current position index tracker (general)
!       integer :: i_brm, i_brm_shift

!       integer :: i_bfm_interpJ


!       ! freqs
!       real(bsa_real_t) :: base_fi, base_fj, fi, fj


!       real(bsa_real_t), allocatable :: bfm_new_left(:, :), bfm_new_right(:, :)
!       real(bsa_real_t), allocatable :: bfm_interp(:, :)

!       real(bsa_real_t) :: dfJtail, dfJhead, dfIcurr, dfIprev, dfJ_oldtmp
!       real(bsa_real_t) :: bfmtail(dimM_bisp_), bfmhead(dimM_bisp_)


!       real(bsa_real_t), allocatable :: brm(:, :)
!       real(bsa_real_t) :: intg(dimM_bisp_)


!       ! BUG: for the moment, we take the max
!       !      such to ensure having the SAME n. of points
!       !      along both sides.
!       !      Later, consider supporting more general approach.
!       interp_fact = max(this%policy_%interp_I_fct_, this%policy_%interp_J_fct_)

! 		! take a backup since this will be decremented by one
!       ! each time we move pi_curr
!       njOld = this%nj_ - 1

      
!       ! get unary delta freqs increments in I and J directions
!       ! (old values, to be interpolated)
!       df_cst = getPointsDistance(this%Cpt_, this%Apt_) / njOld

!       ! NOTE: first two refer to MAJOUR direction (J), then MINOR (I)
!       call this%getRotatedUnaryDF(dfJi_old, dfJj_old, dfIi_old, dfIj_old)

!       ! actualise them
!       dfIi_old = dfIi_old * df_cst
!       dfIj_old = dfIj_old * df_cst
!       dfJi_old = dfJi_old * df_cst
!       dfJj_old = dfJj_old * df_cst


!       ! get interpolated deltas (in GRS)
!       dfIi_interp = dfIi_old / interp_fact
!       dfIj_interp = dfIj_old / interp_fact
!       dfJi_interp = dfJi_old / interp_fact
!       dfJj_interp = dfJj_old / interp_fact



!       ! get absolute deltas (in LRS) along I and J directions
!       ! NOTE: they refer to interpolation
!       dfI_old = df_cst
!       dfJ_old = df_cst
!       dfI_interp = df_cst / interp_fact
!       dfJ_interp = df_cst / interp_fact


!       ! influence areas for integration
!       dw_cst = dfI_interp * CST_PIt2
!       ctr_infl_r = dw_cst*dw_cst
!       brd_infl_r = ctr_infl_r / 2._RDP
!       vtx_infl_r = brd_infl_r / 2._RDP
!       brd_infl_t = brd_infl_r
!       vtx_infl_t = vtx_infl_r / 2._RDP

!       ! get actualised refinements (along borders)
!       ni = (this%ni_ - 1) * interp_fact + 1
!       nj = njOld * interp_fact + 1
!       njNew_picurr = nj


!       ! n of points to interpolate between two known points' direction lines
!       ! Refers to the number of points (or equally the number of segments)
!       ! that are contained between two known pre-meshed points
!       ! TODO: use max ref to have same n of points along both sides
!       nipI = interp_fact - 1
!       nipJ = nipI

!       ! TODO: deltas along the hypotenuse
!       !       Last section when moving pj_head
!       df_diag_old    = sqrt(dfI_old**2 + dfJ_old**2)
!       df_diag_interp = df_diag_old / interp_fact


!       ! allocate data
!       zNintrpPts = getTriangZoneEquivNPts(ni, nj)
!       allocate(brm(dimM_bisp_, zNintrpPts))
!       brm = 0._RDP


!       ! NOTE: allocate intil before-last section.
!       !       Keep upper triang shaped zone out
!       !       treating it seperately.
!       !       This means that we have to take one interp section
!       !       out from allocation.
!       ! NOTE: however, or left, needed to store nj values, 
!       !       since at first iteration along J, we have 
!       !       nj point. Not optimal..
!       ! NOTE: 0 assignment is like memset() in C.
!       allocate(bfm_new_left(dimM_bisp_, nj))
!       bfm_new_left  = 0._RDP
!       allocate(bfm_new_right(dimM_bisp_, nj - nipJ - 1))
!       bfm_new_right = 0._RDP

!       allocate(bfm_interp, source=bfm_new_left)
!       bfm_interp = 0._RDP




!       ! before starting, interpolate very first column
!       ! along J dir, to get new mesh from old
!       ! Then for all the others, it will be done inside
!       ! the loop over the NI (old) points of the old mesh.
!       ! NOTE: integrate as well.
      
!       ! init point vertex
!       base_fi = this%Cpt_%freqI()
!       base_fj = this%Cpt_%freqJ()
!       fi = base_fi
!       fj = base_fj
!       bfmtail = bfm_undump(:, 1)
!       bfm_new_left(:, 1) = bfmtail
!       brm(:, 1) = getBRM_msh(bfm_new_left(:, 1), base_fi, base_fj)
!       intg      = brm(:, 1) * vtx_infl_r

!       i_brm = 1
!       do pjhead = 2, this%nj_

!          ! OK because we stored irt NJ majour.
!          bfmhead = bfm_undump(:, pjhead)

!          ! once we moved head, restore init
!          ! distances from head/tail
!          dfJhead = dfJ_old
!          dfJtail = 0._RDP

!          do pjtail = 1, nipJ

!             fi = fi + dfJi_interp
!             fj = fj + dfJj_interp

!             ! update actual distances head/tail
!             dfJtail = dfJtail + dfJ_interp
!             dfJhead = dfJhead - dfJ_interp

!             ! interpolation along J dir (between HEAD-TAIL)
!             ! NOTE: save it for later use.
!             i_brm = i_brm + 1
!             bfm_new_left(:, i_brm) = (bfmhead * dfJtail + bfmtail * dfJhead) / dfJ_old
!             brm(:, i_brm) = getBRM_msh(bfm_new_left(:, i_brm), fi, fj)

!             ! NOTE: it is a border point
!             intg = intg + brm(:, i_brm) * brd_infl_r
!          enddo ! pjtail = 1, nipJ

!          ! here, treat head, TAIL==HEAD
!          fi = fi + dfJi_interp
!          fj = fj + dfJj_interp

!          i_brm = i_brm + 1
!          bfm_new_left(:, i_brm) = bfmhead
!          brm(:, i_brm) = getBRM_msh(bfm_new_left(:, i_brm), fi, fj)

!          ! NOTE: it is a border point, except for the very last one (VERTEX_t)
!          intg = intg + brm(:, i_brm) * brd_infl_r

!          ! old head becomes new tail
!          ! TODO: we could change bfmhead here, so that
!          !       it might be already ready for next loop.
!          bfmtail = bfmhead
!       enddo

!       ! NOTE: removing excess contribution for last HEAD (VERTEX_t)
!       intg  = intg + brm(:, i_brm) * (vtx_infl_t - brd_infl_r)

!       i_bfm  = this%nj_


!       ! NOTE: starting from 2 because we have to take
!       !       infl lines in consequent groups of two.
!       do picurr = 2, this%ni_

!          ! before doing any computation,
!          ! we need to interpolate BFM along J
!          ! at new CURRENT (I) infl line
!          ! NOTE: once we go through, integrate as well.


!          ! computing BRM offset.
!          ! NOTE: njNew_picurr refers at NJ points in the new mesh
!          !       at pi_prev infl line (along J).
!          i_brm_shift = i_brm
!          njNew_tmp   = njNew_picurr
!          do pjhead = 1, nipI
!             njNew_tmp   = njNew_tmp - 1
!             i_brm_shift = i_brm_shift + njNew_tmp
!          enddo


!          ! reset freqs to point to new base -> prev base moved by old deltas!
!          ! (now pi_prev-pj_tail)
!          ! NOTE: still keep prev base in memory here since they might serve later.
!          fi = base_fi + dfIi_old
!          fj = base_fj + dfIj_old

!          i_bfm   = i_bfm + 1
!          bfmtail = bfm_undump(:, i_bfm)
!          bfm_new_right(:, 1) = bfmtail
         
!          i_brm_shift         = i_brm_shift + 1
!          ! brm(:, i_brm_shift) = getBRM_msh(bfm_new_right(:, 1), base_fi, base_fj)
!          brm(:, i_brm_shift) = getBRM_msh(bfm_new_right(:, 1), fi, fj)
!          intg = intg + brm(:, i_brm_shift) * brd_infl_r

!          i_bfm_interpJ = 1
!          do pjhead = 2, njOld

!             i_bfm   = i_bfm + 1
!             bfmhead = bfm_undump(:, i_bfm)
   
!             ! once we moved head, restore init
!             ! distances from head/tail
!             dfJhead = dfJ_old
!             dfJtail = 0._RDP
   
   
!             do pjtail = 1, nipJ
   
!                fi = fi + dfJi_interp
!                fj = fj + dfJj_interp
   
!                ! update actual distances head/tail
!                dfJtail = dfJtail + dfJ_interp
!                dfJhead = dfJhead - dfJ_interp
   
!                ! interpolation along J dir (between HEAD-TAIL)
!                ! NOTE: save it for later use.
!                i_bfm_interpJ = i_bfm_interpJ + 1
!                bfm_new_right(:, i_bfm_interpJ) = &
!                   (bfmhead * dfJtail + bfmtail * dfJhead) / dfJ_old

!                i_brm_shift = i_brm_shift + 1
!                brm(:, i_brm_shift) = getBRM_msh(bfm_new_right(:, i_bfm_interpJ), fi, fj)
   
!                ! NOTE: it is a border point
!                intg  = intg + brm(:, i_brm_shift) * ctr_infl_r
!             enddo ! pjtail = 1, nipJ
   
!             ! here, treat head, TAIL==HEAD
!             fi = fi + dfJi_interp
!             fj = fj + dfJj_interp
!             i_bfm_interpJ = i_bfm_interpJ + 1
!             bfm_new_right(:, i_bfm_interpJ) = bfmhead
!             i_brm_shift = i_brm_shift + 1
!             brm(:, i_brm_shift) = getBRM_msh(bfm_new_right(:, i_bfm_interpJ), fi, fj)
   
!             ! NOTE: it is a center point, except for the very last one (BORDER_t)
!             intg = intg + brm(:, i_brm_shift) * ctr_infl_r
   
!             ! old head becomes new tail
!             bfmtail = bfmhead
!          enddo ! pjhead = 2, njOld


!          ! removing excess of very last HEAD
!          ! accounted as center, it is BORDER_t
!          ! NOTE: even worse for very last HEAD which happens to be B point.
!          !       there, it is a VERTEX_t point.
!          !       However, it's the very last element in brm , we can remove it after.
!          intg = intg + brm(:, i_brm_shift) * (ctr_infl_r - brd_infl_t)


!          ! backup NJ new intrp n. of points for next iteration
!          ! NOTE: also, refers to NJ point at pi_curr J inlf line.
!          njNew_piprev = njNew_picurr
!          njNew_picurr = i_bfm_interpJ ! NOTE: at last iter, should be 1




!          ! ok, here we now have BFM values (interp along J)
!          ! at CURR and PREV (I) index pointers.
!          ! We have to interpolate along I between CURR and PREV, i.e.
!          ! prev has to start moving toward CURR.


!          ! reset I-dir CURR-PREV distances
!          dfIcurr    = dfI_old
!          dfIprev    = 0._RDP
         
!          dfJ_oldtmp = dfJ_old
!          njtmp      = nipJ

!          ! interpolate along I
!          do piprev = 1, nipI

!             ! bulk I-dir interpolation until pj_head section level.
!             ! Then, after treat that triang shaped zone separately.
            
!             dfIprev = dfIprev + dfI_interp
!             dfIcurr = dfIcurr - dfI_interp

!             bfm_interp(:, 1 : njNew_picurr) = &
!                (  bfm_new_left (:, 1 : njNew_picurr) * dfIcurr + &
!                   bfm_new_right(:, 1 : njNew_picurr) * dfIprev ) / dfI_old

!             ! once we have the values, go through them to integrate
!             ! NOTE: reset base freqs pointers, this time moving them along INTERP mesh
!             fi = base_fi + (dfIi_interp * piprev)
!             fj = base_fj + (dfIj_interp * piprev)

!             i_brm = i_brm + 1
!             brm(:, i_brm) = getBRM_msh(bfm_interp(:, 1), fi, fj)
!             intg = intg + brm(:, i_brm) * brd_infl_r

            
!             do pjtail = 2, njNew_picurr

!                fi = fi + dfJi_interp
!                fj = fj + dfJj_interp
               
!                i_brm = i_brm + 1
!                brm(:, i_brm) = getBRM_msh(bfm_interp(:, pjtail), fi, fj)
!                intg = intg + brm(:, i_brm) * ctr_infl_r
!             enddo


!             !==============================================
!             ! here we have to treat triang shaped zone.

!             ! NOTE: tmp head is interpolated along hypotenuse!!
!             bfmhead = (&
!                bfm_new_left (:, njNew_piprev) * (df_diag_interp * piprev)  + &
!                bfm_new_right(:, njNew_picurr) * (df_diag_interp * (nipI+1-piprev)) ) / df_diag_old


!             ! once we have the head, continue interpolating along J between tail
!             ! and tmp head
!             ! NOTE: everytime we move pi_prev, actual Nj points reduce by 1
!             njtmp      = njtmp - 1
!             dfJ_oldtmp = dfJ_oldtmp - dfJ_interp
!             dfJhead    = dfJ_oldtmp
!             dfJtail    = 0._RDP
!             ipos       = njNew_picurr ! tail index position
!             do pjtail = 1, njtmp

!                fi = fi + dfJi_interp
!                fj = fj + dfJj_interp

!                dfJtail = dfJtail + dfJ_interp
!                dfJhead = dfJhead - dfJ_interp

!                ! NOTE: tail is stored in bfm_interp(:, njNew_picurr)
!                ipos = ipos + 1
!                bfm_interp(:, ipos) = &
!                   (bfmhead * dfJtail + dfJhead * bfm_interp(:, njNew_picurr)) / dfJ_oldtmp

!                i_brm = i_brm + 1
!                brm(:, i_brm) = getBRM_msh(bfm_interp(:, ipos), fi, fj)
!                intg = intg + brm(:, i_brm) * ctr_infl_r
!             enddo

!             ! treat HEAD (on hypotenuse) separately
!             fi = fi + dfJi_interp
!             fj = fj + dfJj_interp

!             i_brm = i_brm + 1
!             brm(:, i_brm) = getBRM_msh(bfmhead, fi, fj)
!             intg = intg + brm(:, i_brm) * vtx_infl_t

!          enddo ! piprev = 1, nipI


!          ! now we can update (PREV) base freqs.
!          base_fi = base_fi + dfIi_old
!          base_fj = base_fj + dfIj_old


!          bfm_new_left(: , 1:njNew_picurr) = bfm_new_right(:, 1:njNew_picurr)


!          ! once finished with this section (CURR-PREV), since we skip J column
!          ! at pi_curr infl line, we reset general BRM index to point to
!          ! previously shifted one.
!          i_brm = i_brm_shift

         
!          ! NOTE: at last iter, should be 1!
!          ! NOTE: this is true just because ni==nj.
!          njOld = njOld - 1

!       enddo ! picurr = 2, this%ni_


! #ifdef __BSA_DEBUG
!       if (i_brm /= zNintrpPts) call bsa_Abort('"i_bfm" does not equal triang zone''s n. of points.')
! #endif


!       ! removing overestimation for B point
!       ! NOTE: should be stored at very last index!!
!       intg = intg - brm(:, i_brm) * vtx_infl_t

!       ! updating main integral
!       m3mr_msh_ptr_ = m3mr_msh_ptr_ + intg

!    end subroutine interpolateTZ_HTPC_v2

















end submodule