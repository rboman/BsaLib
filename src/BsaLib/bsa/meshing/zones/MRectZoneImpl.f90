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
submodule(BsaLib_MRectZone) BsaLib_MRectZoneImpl

#include "../../../precisions"

! #ifndef BSA_M3MF_ONLY_PREMESH_
! # define BSA_M3MF_ONLY_PREMESH_ 0
! #else
! # if (BSA_M3MF_ONLY_PREMESH_ != 0 && BSA_M3MF_ONLY_PREMESH_ != 1)
! #  undef BSA_M3MF_ONLY_PREMESH_
! #  define BSA_M3MF_ONLY_PREMESH_ 0
! # endif
! #endif

   use BsaLib_CONSTANTS
   use BsaLib_Data, only: bsa_Abort
   use BsaLib_IO, only: INFOMSG, WARNMSG, ERRMSG, MSGCONT, DBGMSG, NOTEMSG
   implicit none


contains


   module function MRectZone_t_custom_constructor(rot, name) result(this)
      real(RDP), intent(in), optional :: rot
      character(len=*), intent(in), optional :: name
      ! NOTE: compiler uses default initialisation here (built-in)
      type(MRectZone_t) :: this

      ! From Base type
      call DefaultInitBaseZone(this)

      if (present(rot))  this%rot_  = rot

      ! NOTE: here no need to call setName() method.
      if (present(name)) call this%zoneName(name)
   end function




   !> Gets rect base along I-dir
   module elemental function baseI_rct(this) result(res)
      class(MRectZone_t), intent(in) :: this
      real(RDP) :: res

      res = this%base_I_
   end function


   !> Gets rect base along J-dir
   module elemental function baseJ_rct(this) result(res)
      class(MRectZone_t), intent(in) :: this
      real(RDP) :: res

      res = this%base_J_
   end function






   !> Get A point.
   !> A point is the point defined, starting from I point, 
   !> along the J-dir (Y-axis) parallel side.
   module pure function getAPoint(this) result(pt)
      class(MRectZone_t), intent(in) :: this
      type(MPoint_t) :: pt

      pt = this%Ipt_%getNewPointFromDistAndRot(this%base_J_, this%rot_)
   end function



   !> Get B point.
   !> B point is the point defined, starting from I point, 
   !> along the I-dir (X-axis) parallel side.
   module pure function getBPoint(this) result(pt)
      class(MRectZone_t), intent(in) :: this
      type(MPoint_t) :: pt
      real(RDP) :: ang

      ! in 4-th quadrant, decrement by 3/2*pi -> first quadrant
      if (this%rot_ > CST_PIt3d2) then

         ang = this%rot_ - CST_PIt3d2

      else ! otherwise, increment by 1/2*pi
         
         ang = this%rot_ + CST_PId2
      endif
      pt = this%Ipt_%getNewPointFromDistAndRot(this%base_I_, ang)
   end function



   module subroutine deduceDeltas_rect(this)
      class(MRectZone_t), intent(inout) :: this
      integer :: nsegi, nsegj

      !
      this%refmts_set_ = .true.

      ! n. of segments along I-dir
      nsegi = this%ni_ - 1

      ! n. of segments along J-dir
      nsegj = this%nj_ - 1

      this%deltaf_I_ = this%base_I_ / nsegi
      this%deltaf_J_ = this%base_J_ / nsegj
   end subroutine






   module subroutine setDeltas(this, dfi, dfj, adapt)
      class(MRectZone_t), intent(inout) :: this
      real(RDP), intent(in)             :: dfi, dfj
      logical, intent(in), optional     :: adapt
      logical :: do_adapt = .false.

      if (present(adapt) .and. adapt) do_adapt = .true.

      this%deltaf_I_   = dfi
      this%deltaf_J_   = dfj
      this%deltas_set_ = .true.
      call this%deduceRefinements(do_adapt)
   end subroutine



   module subroutine deduceRefinements(this, adapt)
      class(MRectZone_t), intent(inout) :: this
      logical, intent(in) :: adapt
      integer :: ni, nj

      ! NOTE: as such, they mean N. OF SEGMENTS
      ni = floor(this%base_I_ / this%deltaf_I_)
      nj = floor(this%base_J_ / this%deltaf_J_)

      ! NOTE: we want even N. of segments -> odd N. of points
      if (.not. mod(ni, 2) == 0) ni = ni + 1
      if (.not. mod(nj, 2) == 0) nj = nj + 1


      ! BUG: if we imposed deltas greater than bases -> at least 1 segment!
      if (ni == 0) ni = 1
      if (nj == 0) nj = 1


      if (adapt) then ! we readapt deltas, to fit in base

         this%deltaf_I_ = this%base_I_ / ni
         this%deltaf_J_ = this%base_J_ / nj

      else ! recompute bases, such to preserve desired deltas.

         this%base_I_ = this%deltaf_I_ * ni
         this%base_J_ = this%deltaf_J_ * nj
      endif

      ! Now, get actual number of points along sides (odd)
      this%ni_ = ni + 1
      this%nj_ = nj + 1
   end subroutine







   !> Defines zone bases (sides) from given Deltas, specifying refinement.
   !> NOTE: Init and End zone points MUST be known.
   module subroutine defineFromDeltas_refinements(this, pt, loc, dfi, dfj, ni, nj)
      class(MRectZone_t), intent(inout) :: this

      !> Point specification at location 'loc'
      class(MPoint_t), intent(in) :: pt

      !> Location of the specified point.
      !> Can be:
      !>   'i' -> init point
      !>   'c' -> centre point
      !>   'e' -> end point
      character(len=1) :: loc

      !> Delta values
      real(RDP), intent(in) :: dfi, dfj

      !> Refinements
      integer, value :: ni, nj
!DIR$ ATTRIBUTES VALUE :: ni, nj


      if (.not. this%isGRSAligned()) call bsa_Abort('Rect zone is not GRS aligned.')


! #ifdef __BSA_DEBUG
!       write(unit_debug_, *) ' @MRectZoneImpl::defineFromDeltas_refinements() : init...'
! #endif


      ! force them odd
      if (mod(ni, 2) == 0) ni = ni + 1
      if (mod(nj, 2) == 0) nj = nj + 1

      this%ni_ = ni
      this%nj_ = nj

      ! here we refer to number of segments (even)
      ni = ni - 1
      nj = nj - 1

      ! deduce bases
      this%base_I_ = ni * dfi
      this%base_J_ = nj * dfj


      this%refmts_set_ = .true.
      this%deltaf_I_   = dfi
      this%deltaf_J_   = dfj


      ! NOTE: loc(x) returns the integer address of variable x !!
      !
      call this%define(pt, loc=loc)


! #ifdef __BSA_DEBUG
!       write(unit_debug_, *) ' @MRectZoneImpl::defineFromDeltas_refinements() : init -- ok.'
! #endif
   end subroutine






   !> Defines zone bases (sides) from given Deltas, specifying max deltas values.
   !> NOTE: Init and End zone points MUST be known.
   !> If force is present and true, forces deltas to reach max values specified
   !> by 'valI' and 'valJ'. In this case, deltas are then slightly adjusted to fit.
   module subroutine defineFromDeltas_maxvalues(this, pt, loc, dfi, dfj, maxF_i, maxF_j, force, exceed)
      class(MRectZone_t), intent(inout) :: this

      !> Point specification, which location is specified by 'loc'
      class(MPoint_t), intent(in) :: pt

      !> Point location:
      !>   'i' -> Init location
      !>   'c' -> centre (NOT YET IMPLEMENTED)
      !>   'e' -> end    (NOT YET IMPLEMENTED)
      character(len=1), intent(in) :: loc

      !> Delta values
      real(RDP), value :: dfi, dfj
!DIR$ ATTRIBUTES VALUE :: dfi, dfj

      !> Max deltas values
      real(RDP) :: maxF_i, maxF_j

      !> adjusts deltas to max values specified.
      logical, intent(in), optional :: force


      logical, intent(in), optional :: exceed


      if (.not. this%isGRSAligned()) &
         call bsa_Abort('Zone is not aligned with GRS.')

      if (.not. loc == 'i') &
         call bsa_Abort(&
            'Cannot define deltas from max values if given point location is not "i" (Init).')


! #ifdef __BSA_DEBUG
!       write(unit_debug_, *) ' @MRectZoneImpl::defineFromDeltas_maxvalues() : init...'
! #endif


      block
         ! local
         character(len = *), parameter :: invalid_max_vals = &
            ERRMSG//'Invalid max values (out of allowed bounds).'
         logical   :: do_force
         real(RDP) :: fi, fj, bi, bj
         integer   :: ni, nj


         if (present(force) .and. force) then
            do_force = .true.
         else
            do_force = .false.
         endif

         ! initialise
         bi = 0._RDP
         bj = 0._RDP
         fi = 0._RDP
         fj = 0._RDP
         ni = 0
         nj = 0

         if (do_force) then

            fi = pt%freqI()
            fj = pt%freqJ()

            if (this%rot_ == 0._RDP) then
               if (maxF_i <= fi .or. maxF_j <= fj) call bsa_Abort(invalid_max_vals)

               bi = maxF_i - fi
               bj = maxF_j - fj

            elseif (this%rot_ == CST_PId2) then ! 1/2 * pi
               if (maxF_i <= fi .or. maxF_j >= fj) call bsa_Abort(invalid_max_vals)

               bj = maxF_i - fi
               bi = fj - maxF_j

            elseif (this%rot_ == CST_PIGREC) then
               if (maxF_i >= fi .or. maxF_j >= fj) call bsa_Abort(invalid_max_vals)
               
               bi = fi - maxF_i
               bj = fj - maxF_j

            elseif (this%rot_ == CST_PIt3d2) then ! 3/2 * pi
               if (maxF_i >= fi .or. maxF_j <= fj) call bsa_Abort(invalid_max_vals)
               
               bj = fi - maxF_i
               bi = maxF_j - fj

            end if


            ! invert deltas if needed
            if (.not. (this%rot_ == 0 .or. this%rot_ == CST_PIGREC)) then
               block
                  real(RDP) :: tmp

                  tmp = dfi
                  dfi = dfj
                  dfj = tmp
               end block
            endif

            ! get approx n. of segments
            ni = ceiling(bi / dfi)
            nj = ceiling(bj / dfj)

            ! get refactored deltas
            dfi = bi / ni
            dfj = bj / nj

            ! actual n. of points (odd)
            ni = ni + 1
            nj = nj + 1



         else ! DO NOT force (keep those deltas)


            block
               real(RDP) :: di, dj
               logical   :: do_exceed = .false.

               
               ! defaults
               if (present(exceed) .and. exceed) do_exceed= .true.
               di = dfi
               dj = dfj
               ni = 1
               nj = 1


               if (this%rot_ == 0._RDP) then

                  fi = pt%freqI() + di
                  fj = pt%freqJ() + dj

                  if (maxF_i < fi .or. maxF_j < fj) call bsa_Abort(invalid_max_vals)

                  do while (fi <= maxF_i)
                     ni = ni + 1
                     fi = fi + di
                  enddo
                  do while (fj <= maxF_j)
                     nj = nj + 1
                     fj = fj + dj
                  enddo


               elseif (this%rot_ == CST_PId2) then ! 1/2 * pi
                  fi = pt%freqI() + dj
                  fj = pt%freqJ() - di
                  
                  if (maxF_i < fi .or. maxF_j > fj) call bsa_Abort(invalid_max_vals)

                  do while (fi <= maxF_i)
                     nj = nj + 1
                     fi = fi + dj
                  enddo
                  do while (fj >= maxF_j)
                     ni = ni + 1
                     fj = fj - di
                  enddo


               elseif (this%rot_ == CST_PIGREC) then
                  di = - di
                  dj = - dj
                  fi = pt%freqI() + di
                  fj = pt%freqJ() + dj

                  if (maxF_i > fi .or. maxF_j > fj) call bsa_Abort(invalid_max_vals)

                  do while (fi >= maxF_i)
                     ni = ni + 1
                     fi = fi + di
                  enddo
                  do while (fj >= maxF_j)
                     nj = nj + 1
                     fj = fj + dj
                  enddo

                  di = abs(di)
                  dj = abs(dj)


               ! TODO: check this
               elseif (this%rot_ == CST_PIt3d2) then ! 3/2 * pi
                  fi = pt%freqI() - dj
                  fj = pt%freqJ() + di

                  if (maxF_i > fi .or. maxF_j < fj) call bsa_Abort(invalid_max_vals)

                  do while (fi >= maxF_i)
                     nj = nj + 1
                     fi = fi + dj
                  enddo
                  do while (fj <= maxF_j)
                     ni = ni + 1
                     fj = fj - di
                  enddo

               end if ! rot


               if (ni == 1 .or. nj == 1) call bsa_Abort('At least one max value is too small.')


               if (this%rot_ == 0._RDP) then
                  bi = (ni - 1) * di
                  bj = (nj - 1) * dj

                  if (do_exceed) then
                     if (pt%freqI() + bi < maxF_i) then
                        bi = ni + di
                        ni = ni + 1
                     endif
                     if (pt%freqJ() + bj < maxF_j) then
                        bj = nj * dj
                        nj = nj + 1
                     endif
                  endif


               elseif (this%rot_ == CST_PId2) then
                  bi = (ni - 1) * dj
                  bj = (nj - 1) * di

                  if (do_exceed) then
                     if (pt%freqI() + bj < maxF_i) then
                        bj = nj * di
                        nj = nj + 1
                     endif
                     if (pt%freqJ() - bi > maxF_j) then
                        bi = ni * dj
                        ni = ni + 1
                     endif
                  endif


               elseif (this%rot_ == CST_PIGREC) then
                  bi = (ni - 1) * di
                  bj = (nj - 1) * dj

                  if (do_exceed) then
                     if (pt%freqI() - bi > maxF_i) then
                        bi = ni * di
                        ni = ni + 1
                     endif
                     if (pt%freqJ() - bj > maxF_j) then
                        bj = nj * dj
                        nj = nj + 1
                     endif
                  endif


               elseif (this%rot_ == CST_PIt3d2) then
                  bi = (ni - 1) * dj
                  bj = (nj - 1) * di

                  if (do_exceed) then
                     if (pt%freqI() - bj > maxF_i) then
                        bj = nj * di
                        nj = nj + 1
                     endif
                     if (pt%freqJ() + bi < maxF_j) then
                        bi = ni * dj
                        ni = ni + 1
                     endif
                  endif

               end if ! rot

            end block


         endif ! force

         
         ! backup actual n. of points
         this%ni_ = ni
         this%nj_ = nj


         ! BUG: this is a code copy
         this%refmts_set_ = .true.
         this%deltaf_I_   = dfi
         this%deltaf_J_   = dfj


         call this%define(pt, loc, bi, bj)

      end block


! #ifdef __BSA_DEBUG
!       write(unit_debug_, *) ' @MRectZoneImpl::defineFromDeltas_maxvalues() : init -- ok.'
! #endif
   end subroutine






   module subroutine defineFromEndPtCoordAndBase_norm(&
      this, Pi, coord_val, coord_ty_ch, baseval, base_dir, called)

      class(MRectZone_t), intent(inout) :: this
      class(MPoint_t), intent(in) :: Pi
      real(RDP), intent(in) :: coord_val
      character(len = 1), intent(in) :: coord_ty_ch
      real(RDP), intent(in) :: baseval
      character(len = 1), intent(in)  :: base_dir
      logical, intent(in)             :: called

      real(RDP)      :: ang
      type(MPoint_t) :: pt

      if (.not. (coord_ty_ch == 'i' .or. coord_ty_ch == 'j')) &
         call bsa_Abort('Unvalid coordinate type identifier. Must be one of "i"/"j".')

      if (.not. (base_dir == 'i' .or. base_dir == 'j')) &
         call bsa_Abort('Unvalid base direction identifier. Must be one of "i"/"j".')


! #ifdef __BSA_DEBUG
!       write(unit_debug_, *) ' @MRectZoneImpl::defineFromEndPtCoordAndBase_norm() : init...'
! #endif

      ! backup init point
      this%Ipt_ = MPoint(Pi)

      
      ! get A or B point
      if (base_dir == 'i') then ! B point

         this%base_I_ = baseval

         ! if 4th quadrant, go back to 1st one
         if (this%rot_ >= CST_PIt3d2) then
            ang = this%rot_ - CST_PIt3d2
         else
            ang = this%rot_ + CST_PId2
         endif

      else ! A pt

         this%base_J_ = baseval
         ang = this%rot_
      endif


      ! NOTE: here pt is either A or B, depending on base passed
      pt = this%Ipt_%getNewPointFromDistAndRot(baseval, ang)

      ! deduce missing end point coordinate
      ! From it, compute missing rect base for complete definition
      call this%getOtherBase(pt, base_dir, coord_ty_ch, coord_val)

      if (called) return

      ! once we have both bases, redefine deltas
      ! NOTE: this assumes refinements have been already set
      call this%deduceDeltas()

! #ifdef __BSA_DEBUG
!       write(unit_debug_, *) ' @MRectZoneImpl::defineFromEndPtCoordAndBase_norm() : init -- ok.'
! #endif
   end subroutine defineFromEndPtCoordAndBase_norm




   module subroutine defineFromEndPtCoordAndBase_forceDeltas(&
      this, Pi, coord_val, coord_ty_ch, baseval, base_dir, dfi, dfj)

      class(MRectZone_t), intent(inout) :: this
      class(MPoint_t), intent(in) :: Pi
      real(RDP), intent(in) :: coord_val
      character(len = 1), intent(in) :: coord_ty_ch
      real(RDP), intent(in) :: baseval
      character(len = 1), intent(in)  :: base_dir
      real(RDP), intent(in)           :: dfi, dfj


! #ifdef __BSA_DEBUG
!       write(unit_debug_, *) ' @MRectZoneImpl::defineFromEndPtCoordAndBase_forceDeltas() : init...'
! #endif

      call this%defineFromEndPtCoordAndBase_norm(&
         Pi, coord_val, coord_ty_ch, baseval, base_dir, .true.)

      ! forcing deltas
      call this%setDeltas(dfi, dfj, .true.)

! #ifdef __BSA_DEBUG
!       write(unit_debug_, *) ' @MRectZoneImpl::defineFromEndPtCoordAndBase_forceDeltas() : init -- ok.'
! #endif
   end subroutine defineFromEndPtCoordAndBase_forceDeltas










   module subroutine define(this, pt, loc, base_i, base_j)
      class(MRectZone_t), intent(inout) :: this
      class(MPoint_t), intent(in)       :: pt
      character(len=1), optional, intent(in) :: loc
      real(RDP), intent(in), optional        :: base_i, base_j

      character(len=1) :: location = 'i'

      if (present(loc)) location = loc

      if (location == 'c') then
         
         block
            real(RDP) :: bid2, bjd2

            if (present(base_i)) then
               this%base_I_ = base_i
               bid2         = base_i / 2
            else
               bid2         = this%base_I_ / 2
            endif
      
            if (present(base_j)) then
               this%base_J_ = base_j
               bjd2         = base_j / 2
            else
               bjd2         = this%base_J_ / 2
            endif

            this%Ipt_ = MPoint(pt%freqI() - bid2, pt%freqJ() - bjd2)
            this%Ept_ = MPoint(pt%freqI() + bid2, pt%freqJ() + bjd2)
         end block

      else ! loc == 'i' or 'e'

         ! NOTE: don't forget to update bases if passed!
         if (present(base_i)) this%base_I_ = base_i
         if (present(base_j)) this%base_J_ = base_j
         
         call this%setIEpts(pt, loc)
      endif
   end subroutine define






   module subroutine setIEpts(this, pt, loc)
      class(MRectZone_t), intent(inout) :: this
      class(MPoint_t), intent(in)       :: pt
      character(len=1), intent(in)      :: loc

      real(RDP) :: c, s, ang, di, dj

      if (.not. (loc == 'i' .or. loc == 'e')) &
         call bsa_Abort('Invalid point location.')


      if (this%rot_ < CST_PId2) then ! FIRST quadrant

         c = cos(this%rot_)
         s = sin(this%rot_)

         di = this%base_I_ * c + this%base_J_ * s
         dj = this%base_J_ * c - this%base_I_ * s

      elseif (this%rot_ < CST_PIGREC) then ! SECOND quadrant

         ang = this%rot_ - CST_PId2
         c   = cos(ang)
         s   = sin(ang)

         di = this%base_J_ * c - this%base_I_ * s
         dj = - (this%base_J_ * s + this%base_I_ * c)

      elseif (this%rot_ < CST_PIt3d2) then ! THIRD quadrant

         ang = this%rot_ - CST_PIGREC
         c   = cos(ang)
         s   = sin(ang)

         di = - (this%base_J_ * s + this%base_I_ * c)
         dj = - this%base_J_ * c  + this%base_I_ * s


      elseif (this%rot_ < CST_PIt2) then ! FOURTH quadrant

         ang = this%rot_ - CST_PIt3d2
         c   = cos(ang)
         s   = sin(ang)

         di = - this%base_J_ * c + this%base_I_ * s
         dj =   this%base_J_ * s + this%base_I_ * c
      endif

      

      select case (loc)

      case ('i')
         this%Ipt_ = MPoint(pt)

         this%Ept_ = MPoint(pt%freqI() + di, pt%freqJ() + dj)

      case ('e')
         this%Ept_ = MPoint(pt)

         this%Ipt_ = MPoint(pt%freqI() - di, pt%freqJ() - dj)
      end select
   end subroutine





   !> Avoid setting a delta smaller than given limit
   module elemental impure subroutine validateDeltas(this, lval)
      class(MRectZone_t), intent(inout) :: this
      real(RDP), intent(in) :: lval

      real(RDP) :: dfi, dfj
      logical   :: coarsen = .false.

      dfi = this%deltaf_I_
      dfj = this%deltaf_J_
      if (dfi < lval .or. dfi > lval) then
         dfi     = lval
         coarsen = .true.
      endif
      if (dfj < lval .or. dfj > lval) then
         dfj     = lval
         coarsen = .true.
      endif
      if (coarsen) then
         print '( 1x, 2a, g10.5, " [Hz])." )', &
            WARNMSG, 'Detected at least one zone deltas too  small/big. Setting to optimal value  (', &
            lval
         call this%setDeltas(dfi, dfj, .true.)
      endif
   end subroutine




   ! !> Returns euqivalent rotation in the FIRST quadrant.
   ! elemental function getFirstQuadEquivRot(rot) result(rot1)
   !    real(RDP), intent(in) :: rot
   !    real(RDP) :: rot1

   !    if (rot < CST_PId2) then ! FIRST quad, keep it
   !       rot1 = rot
   !    elseif (rot < CST_PIGREC) then ! SECOND quad
   !       rot1 = rot - CST_PId2
   !    elseif (rot < CST_3d2) then ! THIRD quad
   !       rot1 = rot - CST_PIGREC
   !    elseif (rot < CST_PIt2) then ! FOURTH quad
   !       rot1 = rot - CST_3d2
   !    endif
   ! end function





   !> Automatically computes the second remaining (unknown) rect base
   !> based on the point's coordinates and the base that we already defined.
   !> BUG: cannot see when a base is negative (i.e. when END pt is "behind" INIT one).
   module subroutine getOtherBase(this, pt, base_dir, known_coord, coord_val)
      class(MRectZone_t), intent(inout) :: this
      class(MPoint_t), intent(in)  :: pt
      character(len=1), intent(in) :: base_dir, known_coord
      real(RDP), intent(in) :: coord_val

      real(RDP) :: kd, rot, cd, fi, fj
      type(MPoint_t) :: Pe

      
      ! NOTE: preinitialise to avoid errors !!!!!!!!!!!!
      cd = 0._RDP


      if (base_dir == 'i') then ! we passed B point, we know side along I-dir

         if (known_coord == 'i') then ! we search DJ

            ! BUG: not always 0 when it should be
            kd = abs(coord_val - pt%freqI())

            ! BUG: forcing it to zero if below some precision
            if (kd < MACHINE_PRECISION) then

#ifdef __BSA_DEBUG
               write(unit_debug_, '(a, a)') &
                  WARNMSG, '(1) kd < machine precision. Assuming kd == 0.'
#endif
            else

               if (this%rot_ < CST_PId2) then
                  rot = this%rot_
                  cd  = kd / tan(rot)
               elseif (this%rot_ < CST_PIGREC) then
                  rot  = this%rot_ - CST_PId2
                  cd   = - kd * tan(rot)
               elseif (this%rot_ < CST_PIt3d2) then
                  rot  = this%rot_ - CST_PIGREC
                  cd   = - kd / tan(rot)
               elseif (this%rot_ < CST_PIt2) then
                  rot = this%rot_ - CST_PIt3d2
                  cd  = kd * tan(rot)
               endif
            endif

            fj = pt%freqJ() + cd
            Pe = MPoint(coord_val, fj)


         elseif (known_coord == 'j') then ! we search DI

            kd = abs(coord_val - pt%freqJ())

            ! BUG: forcing it to zero if below some precision
            if (kd < MACHINE_PRECISION) then
               
#ifdef __BSA_DEBUG
               write(unit_debug_, '(a, a)') &
                  WARNMSG, '(2) kd < machine precision. Assuming kd == 0.'
#endif
            else

               if (this%rot_ < CST_PId2) then
                  rot = this%rot_
                  cd  = kd * tan(rot)
               elseif (this%rot_ < CST_PIGREC) then
                  rot  = this%rot_ - CST_PId2
                  cd   = kd / tan(rot)
               elseif (this%rot_ < CST_PIt3d2) then
                  rot  = this%rot_ - CST_PIGREC
                  cd   = - kd * tan(rot)
               elseif (this%rot_ < CST_PIt2) then
                  rot = this%rot_ - CST_PIt3d2
                  cd  = - kd / tan(rot)
               endif
            endif

            fi = pt%freqI() + cd
            Pe = MPoint(fi, coord_val)

         endif ! known_coord

         this%base_J_ = getPointsDistance(pt, Pe)



      elseif (base_dir == 'j') then ! we pass A pt, we know base along J-dir


         if (known_coord == 'i') then

            kd = abs(coord_val - pt%freqI())

            ! BUG: forcing it to zero if below some precision
            if (kd < MACHINE_PRECISION) then
               
#ifdef __BSA_DEBUG
               write(unit_debug_, '(a, a)') &
                  WARNMSG, '(3) kd < machine precision. Assuming kd == 0.'
#endif
            else
               
               if (this%rot_ < CST_PId2) then
                  rot = this%rot_
                  cd  = kd * tan(rot)
               elseif (this%rot_ < CST_PIGREC) then
                  rot  = this%rot_ - CST_PId2
                  cd   = - kd / tan(rot)
               elseif (this%rot_ < CST_PIt3d2) then
                  rot  = this%rot_ - CST_PIGREC
                  cd   = - kd * tan(rot)
               elseif (this%rot_ < CST_PIt2) then
                  rot = this%rot_ - CST_PIt3d2
                  cd  = kd / tan(rot)
               endif
            endif

            fj = pt%freqJ() + cd
            Pe = MPoint(coord_val, fj)


         elseif (known_coord == 'j') then


            kd = abs(coord_val - pt%freqJ())

            ! BUG: forcing it to zero if below some precision
            if (kd < MACHINE_PRECISION) then
               
#ifdef __BSA_DEBUG
               write(unit_debug_, '(a, a)') &
                  WARNMSG, '(4) kd < machine precision. Assuming kd == 0.'
#endif
            else
               
               if (this%rot_ < CST_PId2) then
                  rot = this%rot_
                  cd  = kd / tan(rot)
               elseif (this%rot_ < CST_PIGREC) then
                  rot  = this%rot_ - CST_PId2
                  cd   = - kd * tan(rot)
               elseif (this%rot_ < CST_PIt3d2) then
                  rot  = this%rot_ - CST_PIGREC
                  cd   = - kd / tan(rot)
               elseif (this%rot_ < CST_PIt2) then
                  rot = this%rot_ - CST_PIt3d2
                  cd  = kd * tan(rot)
               endif
            endif

            fi = pt%freqI() + cd
            Pe = MPoint(fi, coord_val)

         endif ! known_coord


         this%base_I_ = getPointsDistance(pt, Pe)


      endif ! base_dir

      ! Now once defined save End point
      this%Ept_ = Pe
   end subroutine getOtherBase





   !> Gets actualised frequency deltas along two main
   !> sides directions (I, J), actualised based on *this
   !> zone rotation w.r.t. GRS.
   module subroutine getIJfsteps(this, dfIx, dfIy, dfJx, dfJy)
      class(MRectZone_t), intent(in) :: this
      real(RDP), intent(out) :: dfIx, dfIy, dfJx, dfJy

      real(RDP) :: c, s, ang

      if (this%rot_ < CST_PId2) then ! FIRST quadrant

         c = cos(this%rot_)
         s = sin(this%rot_)

         dfIx = this%deltaf_I_ * c
         dfIy = - this%deltaf_I_ * s

         dfJx = this%deltaf_J_ * s
         dfJy = this%deltaf_J_ * c

      elseif (this%rot_ < CST_PIGREC) then ! SECOND quadrant

         ang = this%rot_ - CST_PId2
         c   = cos(ang)
         s   = sin(ang)

         dfIx = - this%deltaf_I_ * s
         dfIy = - this%deltaf_I_ * c

         dfJx = this%deltaf_J_ * c
         dfJy = - this%deltaf_J_ * s

      elseif (this%rot_ < CST_PIt3d2) then ! THIRD quadrant

         ang = this%rot_ - CST_PIGREC
         c   = cos(ang)
         s   = sin(ang)

         dfIx = - this%deltaf_I_ * c
         dfIy = this%deltaf_I_ * s

         dfJx = - this%deltaf_J_ * s
         dfJy = - this%deltaf_J_ * c

      elseif (this%rot_ < CST_PIt2) then ! FOURTH quadrant

         ang = this%rot_ - CST_PIt3d2
         c   = cos(ang)
         s   = sin(ang)

         dfIx = this%deltaf_I_ * s
         dfIy = this%deltaf_I_ * c

         dfJx = - this%deltaf_J_ * c
         dfJy = this%deltaf_J_ * s

      endif
   end subroutine







   module function reconstructZoneBaseMesh(this) result(msh)
      class(MRectZone_t), intent(in) :: this
      !> BUG: might be 2-rank array instead of 3!
      real(RDP) :: msh(2, this%nj_, this%ni_)

      real(RDP) :: dfIi, dfIj, dfJi, dfJj
      real(RDP) :: base_fi, base_fj, fi, fj
      integer   :: i, j

      call this%getIJfsteps(dfIi, dfIj, dfJi, dfJj)

      fi = this%Ipt_%freqI()
      fj = this%Ipt_%freqJ()
      base_fi = fi
      base_fj = fj

      msh(:, 1, 1) = [fj, fi]

      ! internal lines
      do j = 2, this%nj_
         fi = fi + dfJi
         fj = fj + dfJj
         msh(:, j, 1) = [fj, fi]
      enddo ! pj_head

      
      ! internal columns
      do i = 2, this%ni_

         base_fi = base_fi + dfIi
         base_fj = base_fj + dfIj

         fi = base_fi
         fj = base_fj

         msh(:, 1, i) = [fj, fi]

         do j = 2, this%nj_
            fi = fi + dfJi
            fj = fj + dfJj
            msh(:, j, i) = [fj, fi]
         enddo ! pj_head
      enddo ! pi_head
   end function 






   !> Actual zone comutation (pre phase).
   module subroutine compute_s(this)
      use BsaLib_Data, only: &
         dimM_bisp_, getBFM_msh, settings  &
         , m3mf_msh_ptr_, msh_NZones, msh_bfmpts_pre_
      class(MRectZone_t), intent(inout) :: this


      if (.not. (this%refmts_set_ .or. this%deltas_set_)) &
         call bsa_Abort('Either deltas or refinements must be set before computing a zone.')

      
      if (this%base_I_ <= MACHINE_PRECISION .or. this%base_J_ <= MACHINE_PRECISION) then
         print '(1x, 4a)', &
            WARNMSG, 'One sub-zone at   ', &
               this%name_(1 : len_trim(this%name_)), '   is empty.'
         goto 998
      endif

      if (this%deltaf_I_ <= MACHINE_PRECISION .or. this%deltaf_J_ <= MACHINE_PRECISION) &
         call bsa_Abort("At least one delta freq is zero.")
      

      block
         real(RDP) :: dfIi, dfIj, dfJi, dfJj
         real(RDP) :: base_fi, base_fj, fi, fj

         integer :: niM1, njM1
         integer :: i, j, idbfm, zNp

         real(RDP), allocatable :: bfm(:, :)

#ifdef BSA_M3MF_ONLY_PREMESH_
         real(RDP) :: dwI, dwJ
         real(RDP) :: ctr_infl, brd_infl, vtx_infl
         real(RDP), allocatable :: intg(:)
#endif

         call this%getIJfsteps(dfIi, dfIj, dfJi, dfJj)

      
#ifdef BSA_M3MF_ONLY_PREMESH_
         ! deltas in [rad/s] (to compute influence areas)
         dwI = this%deltaf_I_ * CST_PIt2
         dwJ = this%deltaf_J_ * CST_PIt2
         ctr_infl = dwI * dwJ
         brd_infl = ctr_infl / 2
         vtx_infl = brd_infl / 2

         allocate(intg(dimM_bisp_))
#endif


         ! get before last refmts indexes (along I and J dirs)
         niM1 = this%ni_ - 1
         njM1 = this%nj_ - 1


         ! allocate memory
         ! NOTE: but dimBISP as 1st dimensions, since
         !       we are gonna calling getBFM for each 
         !       couple of nodal indexes, returning a 
         !       dimBISP 1D vector, so better memory access.
         !       However, later maybe better invert.
         zNp = this%nj_ * this%ni_
         allocate(bfm(dimM_bisp_, zNp))


         !=========================================================
         ! FIRST COLUMN (along J-dir, from I to A)
         !
         base_fi = this%Ipt_%freqI()
         base_fj = this%Ipt_%freqJ()

         fi = base_fi
         fj = base_fj

         bfm(:, 1) = getBFM_msh(fi, fj)
#ifdef BSA_M3MF_ONLY_PREMESH_
         intg(:)   = bfm(:, 1) * vtx_infl
#endif

#ifdef __BSA_CHECK_NOD_COH_SVD
         return
#endif

         ! internal lines
         do j = 2, njM1

            fi = fi + dfJi
            fj = fj + dfJj
            bfm(:, j) = getBFM_msh(fi, fj)
#ifdef BSA_M3MF_ONLY_PREMESH_
            intg(:)   = intg(:) + bfm(:, j) * brd_infl
#endif
         enddo

         ! BUG: handle in case 2 > njM1 ???
         if (njM1 == 1 .and. (.not. j==2)) j = 2

         fi = fi + dfJi
         fj = fj + dfJj
         bfm(:, j) = getBFM_msh(fi, fj)
#ifdef BSA_M3MF_ONLY_PREMESH_
         intg(:)   = intg(:) + bfm(:, j) * vtx_infl
#endif
         idbfm     = j + 1



         !=========================================================
         ! INTERNAL COLUMNS
         !
         do i = 2, niM1

            ! update base freqs moving along I local direction (X)
            base_fi = base_fi + dfIi
            base_fj = base_fj + dfIj
            fi = base_fi
            fj = base_fj

            bfm(:, idbfm) = getBFM_msh(fi, fj)
#ifdef BSA_M3MF_ONLY_PREMESH_
            intg(:)       = intg(:) + bfm(:, idbfm) * brd_infl
#endif
            idbfm         = idbfm + 1

            ! internal lines
            do j = 2, njM1

               fi = fi + dfJi
               fj = fj + dfJj
               bfm(:, idbfm) = getBFM_msh(fi, fj)
#ifdef BSA_M3MF_ONLY_PREMESH_
               intg(:)       = intg(:) + bfm(:, idbfm) * ctr_infl
#endif
               idbfm         = idbfm + 1
            enddo

            ! last line
            ! BUG: handle in case 2 > njM1 ???
            if (njM1 == 1 .and. (.not. j==2)) j = 2

            fi = fi + dfJi
            fj = fj + dfJj
            bfm(:, idbfm) = getBFM_msh(fi, fj)
#ifdef BSA_M3MF_ONLY_PREMESH_
            intg(:)       = intg(:) + bfm(:, idbfm) * brd_infl
#endif
            idbfm         = idbfm + 1
         enddo


         !=========================================================
         ! LAST COLUMN (from B to E)
         !
         ! BUG: handle in case 2 > niM1 ???
         if (niM1 == 1 .and. (.not. i==2)) i = 2

         base_fi = base_fi + dfIi
         base_fj = base_fj + dfIj
         fi = base_fi
         fj = base_fj
         ! first line
         bfm(:, idbfm) = getBFM_msh(fi, fj)
#ifdef BSA_M3MF_ONLY_PREMESH_
         intg(:)       = intg(:) + bfm(:, idbfm) * vtx_infl
#endif
         idbfm         = idbfm + 1

         ! internal lines
         do j = 2, njM1
            fi = fi + dfJi
            fj = fj + dfJj
            bfm(:, idbfm) = getBFM_msh(fi, fj)
#ifdef BSA_M3MF_ONLY_PREMESH_
            intg(:)       = intg(:) + bfm(:, idbfm) * brd_infl
#endif
            idbfm         = idbfm + 1
         enddo

         ! last line
         ! BUG: handle in case 2 > njM1 ???
         if (njM1 == 1 .and. (.not. j==2)) j = 2

         fi = fi + dfJi
         fj = fj + dfJj
         bfm(:, idbfm) = getBFM_msh(fi, fj)
#ifdef BSA_M3MF_ONLY_PREMESH_
         intg(:)       = intg(:) + bfm(:, idbfm) * vtx_infl
#endif

! #ifdef __BSA_DEBUG
         if (idbfm /= zNp) then
            print *, 'idbfm , zNp  =  ', idbfm, zNp
            call bsa_Abort('"idbfm" does not equally tot N of Rect zone''s points.')
         endif
! #endif

         !$omp critical
#ifdef BSA_M3MF_ONLY_PREMESH_
         m3mf_msh_ptr_   = m3mf_msh_ptr_ + (intg * settings%i_bisp_sym_) ! update main integral
#endif
         msh_NZones      = msh_NZones + 1          ! update n. of zones count
         msh_bfmpts_pre_ = msh_bfmpts_pre_ + zNp   ! update tot num of meshing points
         
         ! eventually, update zone with max N of points
         if (zNp > msh_max_zone_NPts) msh_max_zone_NPts = zNp

         call DumpZone(this, bfm)   ! dump zone info
         !$omp end critical

      end block

      ! NOTE: reset them to 0 for ensuring next zone correct setup
      998 continue
      this%refmts_set_ = .false.
      this%deltas_set_ = .false.

! #ifdef __BSA_DEBUG
!       write(unit_debug_, *) ' @MRectZoneImpl::compute_s() : init -- ok.'
! #endif
   end subroutine compute_s






   !> Gets vertex point pt coordinates in Nth quadrant
   !> w.r.t. Center point, in zone's LRS.
   module pure function getNthQuadVtx(this, iquad) result(pt)
      class(MRectZone_t), intent(in) :: this
      integer, intent(in) :: iquad
      type(MPoint_t) :: pt
      real(RDP) :: c, s, rot, di, dj

      if (iquad < 1 .or. iquad > 4) return

      if (this%rot_ < CST_PId2) then

         c = cos(this%rot_)
         s = sin(this%rot_)

         select case (iquad)
         case (1)
            di = this%base_I_ * c + this%base_J_ * s
            dj = this%base_J_ * c - this%base_I_ * s

         case (2)
            di =   this%base_I_ * c
            dj = - this%base_I_ * s

         case (3)
            pt = this%Ipt_
            return

         case (4)
            di = this%base_J_ * s
            dj = this%base_J_ * c

         end select


      elseif (this%rot_ < CST_PIGREC) then

         rot = this%rot_ - CST_PId2
         c   = cos(rot)
         s   = sin(rot)

         select case (iquad)
         case (1)
            di =   this%base_J_ * c
            dj = - this%base_J_ * s

         case (2)
            di =    this%base_J_ * c - this%base_I_ * s
            dj = - (this%base_J_ * s + this%base_I_ * c)

         case (3)
            di = - this%base_I_ * s
            dj = - this%base_I_ * c

         case (4)
            pt = this%Ipt_
            return

         end select

      elseif (this%rot_ < CST_PIt3d2) then

         rot = this%rot_ - CST_PIGREC
         c   = cos(rot)
         s   = sin(rot)

         select case (iquad)
         case (1)
            pt = this%Ipt_
            return

         case (2)
            di = - this%base_J_ * s
            dj = - this%base_J_ * c

         case (3)
            di = - (this%base_J_ * s + this%base_I_ * c)
            dj = - this%base_J_ * c  + this%base_I_ * s

         case (4)
            di = - this%base_I_ * c
            dj =   this%base_I_ * s
         end select

      elseif (this%rot_ < CST_PIt2) then

         rot = this%rot_ - CST_PIt3d2
         c   = cos(rot)
         s   = sin(rot)

         select case (iquad)
         case (1)
            di = this%base_I_ * s
            dj = this%base_I_ * c

         case (2)
            pt = this%Ipt_
            return

         case (3)
            di = - this%base_J_ * c
            dj =   this%base_J_ * s

         case (4)
            di = - this%base_J_ * c + this%base_I_ * s
            dj =   this%base_J_ * s + this%base_I_ * c
         end select

      endif

      pt = MPoint(this%Ipt_%freqI() + di, this%Ipt_%freqJ() + dj)
   end function








   !> Dumps a RECTANGULAR zone.
   !>
   !> NOTE: Each specific zone dumping method is called 
   !>       from the STATIC MSaver_t procedure dump().
   !>       This makes no need for specific saver to have
   !>       their own dump() implementation, since this
   !>       current imlementation is what we are looking for.
   module subroutine dumpRZ(this)
      class(MRectZone_t), intent(in) :: this

      write(unit_dump_bfm_) MZone_ID%RECTANGLE

      ! init pt
      write(unit_dump_bfm_) this%Ipt_%freqI(), this%Ipt_%freqJ()

      write(unit_dump_bfm_) this%rot_
      write(unit_dump_bfm_) this%base_I_, this%base_J_
      
      ! NOTE: maybe useless ?
      write(unit_dump_bfm_) this%ni_, this%nj_

#ifdef __BSA_ZONE_DEBUG
      write(unit=4533, fmt=*) &
         'Refms at  RZ=', trim(this%name_), this%ni_, this%nj_, &
			'thread id= ', omp_get_thread_num()
#endif
   end subroutine dumpRZ





   !> Undumps a RECTANGULAR zone.
   !>
   !> NOTE: Each specific zone dumping method is called 
   !>       from the STATIC procedure UndumpZone() in MZone Module.
   module subroutine undumpRZ(this)
      class(MRectZone_t), intent(inout) :: this
      
      real(RDP)      :: rval1, rval2
      integer        :: ival1, ival2

      ! init point
      read(unit_dump_bfm_) rval1, rval2
      call this%Ipt_%setfreqs(rval1, rval2)

      ! rotation
      read(unit_dump_bfm_) rval1
      this%rot_ = rval1

      ! sides
      read(unit_dump_bfm_) rval1, rval2
      this%base_I_ = rval1
      this%base_J_ = rval2
      
      ! refinements
      read(unit_dump_bfm_) ival1, ival2
      call this%setRefinements(ival1, ival2, .true.)
   end subroutine undumpRZ







   !> Implementation of rect zone interpolation methods
   module subroutine interpolateRZ( this &
#ifdef __BSA_OMP
      , bfm, pdata &
#endif
      & )
      class(MRectZone_t), intent(inout) :: this
#ifdef __BSA_OMP
      real(RDP), intent(in) :: bfm(:, :)
      class(*), pointer, intent(in) :: pdata

      ! NOTE: for the moment only supporting HTPC method
      call interpolateRZ_HTPC_v3(this, bfm, pdata)
#else
      call interpolateRZ_HTPC_v3(this)
#endif
   end subroutine interpolateRZ









   !> Implementation of HTPC interpolation method for a rectangle zone,
   !> including MultiLevel-Refinement for BFM data.
   subroutine interpolateRZ_HTPC_v3( this &
#ifdef __BSA_OMP
      , bfm_undump, pdata &
#endif
      & )
      use BsaLib_Data, only: &
#ifndef __BSA_OMP
         bfm_undump, &
#endif
         dimM_bisp_, getBFM_msh, getBRM_msh, m3mf_msh_ptr_, m3mr_msh_ptr_, settings  &
         , msh_bfmpts_post_, msh_brmpts_post_, do_validate_deltas_ &
         , msh_ZoneLimsInterestModes, peak_exts_ &
         , write_brm_fptr_, do_export_brm_, BrmExportBaseData_t  &
         , I_BKG_PEAK_DELTAF_BFM_REFMT_FCT_, I_RES_PEAK_DELTAF_BFM_REFMT_FCT_ &
         , CODE_PRE_PEAK_OK, CODE_PRE_PEAK_KO &
         , bkg_peakw_

      use BsaLib_MPolicy, only: MPolicy_t
      !
      class(MRectZone_t), intent(inout) :: this
#ifdef __BSA_OMP
      real(RDP), intent(in) :: bfm_undump(:, :)
      class(*), pointer, intent(in) :: pdata
#endif

      type(MPolicy_t) :: pol
      integer   :: ni, nj, ni_bfm_ref_, nj_bfm_ref_
      integer   :: nipI, nipJ, nPtsPost, n_im_, im_idx_
      integer   :: i, ist, n_segs_bfm_ref_i_, n_segs_bfm_ref_j_
      integer   :: n_pts_bfm_ref_i_, n_pts_bfm_ref_j_
      real(RDP) :: dfIi_bfm_lv0_, dfIj_bfm_lv0_, dfJi_bfm_lv0_, dfJj_bfm_lv0_
      real(RDP) :: dfIi_bfm_ref_, dfIj_bfm_ref_, dfJi_bfm_ref_, dfJj_bfm_ref_
      real(RDP) :: dfI_bfm_lv0_, dfJ_bfm_lv0_, dfI_bfm_ref_, dfJ_bfm_ref_
      real(RDP) :: dfIi_brm_interp, dfIj_brm_interp, dfJi_brm_interp, dfJj_brm_interp
      real(RDP) :: dfI_brm_interp_, dfJ_brm_interp_, dwI, dwJ
      real(RDP) :: vtx_infl, brd_infl, ctr_infl
      integer(kind = 4), allocatable :: inter_modes_(:)

      ! HTPC indexes
      integer :: pIcurr, pIprev, pJhead, pJtail

      ! Pos in general BFM undumped data
      integer :: i_bfm_old, i_bfm_ref_i, i_bfm_ref_j
      integer :: i_brm, i_brm_shift, i_brm_write_, i_brm_offsetJ
      integer :: i_bfm_interpJ, i_ftc

      ! freqs
      real(RDP) :: fi_baseptI, fj_baseptI, fi_baseptJ, fj_baseptJ
      real(RDP) :: fi, fj
#ifdef __BSA_OMP
      real(RDP), allocatable, dimension(:) :: fi_v_, fj_v_
#endif


      real(RDP), allocatable :: bfm_new_left(:, :), bfm_new_right(:, :)
      real(RDP), allocatable :: bfm_interp(:, :)

      real(RDP) :: dfJtail, dfJhead, dfIcurr, dfIprev
      real(RDP) :: bfmtail(dimM_bisp_), bfmhead(dimM_bisp_)

      real(RDP), allocatable :: brm(:, :)
      real(RDP) :: intg(dimM_bisp_)

#ifndef BSA_M3MF_ONLY_PREMESH_
      real(RDP) :: vtx_infl_bfm, brd_infl_bfm, ctr_infl_bfm
      real(RDP) :: intg_bfm(dimM_bisp_)
#endif


      pol = this%policy() ! get zone's policy

      ! get original (pre-meshing) deltas (LEVEL 0)
      ! NOTE: keep them in memory unchanged since they might serve later on.
      call this%getIJfsteps(dfIi_bfm_lv0_, dfIj_bfm_lv0_, dfJi_bfm_lv0_, dfJj_bfm_lv0_)


      ! get n. of BFM refinement segments (between two old ones)
      n_segs_bfm_ref_i_ = 1 ! original
      n_segs_bfm_ref_j_ = 1
      do i = 1, pol%n_interp_bfm_lvs_
         n_segs_bfm_ref_i_ = n_segs_bfm_ref_i_ * pol%interp_bfm_I_fct_
         n_segs_bfm_ref_j_ = n_segs_bfm_ref_j_ * pol%interp_bfm_J_fct_
      enddo
      dfI_bfm_lv0_ = this%deltaf_I_
      dfI_bfm_ref_ = dfI_bfm_lv0_ / n_segs_bfm_ref_i_
      dfJ_bfm_lv0_ = this%deltaf_J_
      dfJ_bfm_ref_ = dfJ_bfm_lv0_ / n_segs_bfm_ref_j_


      ! Validate BFM_ref deltas
      ! NOTE: for the moment, check only if too coarse. 
      !       Later on, check also if too fine!
      if (do_validate_deltas_) then

         ! NOTE: 0 denotes that interest modes are to be inferenced from index 1.
         !       In fact, there are 3 scenarios.
         !          1.  next zone is pre-peak, and next peak interest modes' start from 1.
         !              BKG does not include any resonant peak.
         !          2.  next zone is pre-peak, and next peak interest modes' DO NOT start from 1.
         !              BKG does include some resonant peaks (from 1-less-index or next peak zone)
         !          3.  next zone is peak.
         !              BKG does include this, plus all previous resonant peaks.

         im_idx_ = this%id_im_
         if (im_idx_ == 0) im_idx_ = 1
         n_im_   = msh_ZoneLimsInterestModes(im_idx_)

         ! NOTE: this is allowed since in Pre-Mesh we have already +1 incremented pointer index
         !       for all pre-peak zones. So, only negative index is possible for very first 
         !       pre-peak zone right after BKG, for which we had set pointer index to 0!
         if (n_im_ == CODE_PRE_PEAK_OK) then
            i_ftc = I_BKG_PEAK_DELTAF_BFM_REFMT_FCT_
            dwI   = bkg_peakw_
         else
            if (n_im_ < 0) then
               im_idx_ = im_idx_ + 1
               n_im_   = msh_ZoneLimsInterestModes(im_idx_)
            endif

            i_ftc        = I_RES_PEAK_DELTAF_BFM_REFMT_FCT_
            inter_modes_ = msh_ZoneLimsInterestModes(im_idx_ + 1 : im_idx_ + n_im_)
            
            ! BUG: introduce I and J peak widths!
            dwI = minval(peak_exts_(inter_modes_))  ! base MIN deltaf
         endif
         
         if (dfI_bfm_ref_ > dwI / i_ftc) then
            do while (dfI_bfm_ref_ > dwI / i_ftc)
               n_segs_bfm_ref_i_ = n_segs_bfm_ref_i_ + 1
               dfI_bfm_ref_      = dfI_bfm_lv0_ / n_segs_bfm_ref_i_
            enddo
         else
            do while (dfI_bfm_ref_ < dwI / i_ftc .and. n_segs_bfm_ref_i_ > 1)
               n_segs_bfm_ref_i_ = n_segs_bfm_ref_i_ - 1
               dfI_bfm_ref_      = dfI_bfm_lv0_ / n_segs_bfm_ref_i_
            enddo
         endif

         if (dfJ_bfm_ref_ > dwI / i_ftc) then
            do while (dfJ_bfm_ref_ > dwI / i_ftc)
               n_segs_bfm_ref_j_ = n_segs_bfm_ref_j_ + 1
               dfJ_bfm_ref_      = dfJ_bfm_lv0_ / n_segs_bfm_ref_j_
            enddo
         else
            do while (dfJ_bfm_ref_ < dwI / i_ftc .and. n_segs_bfm_ref_j_ > 1)
               n_segs_bfm_ref_j_ = n_segs_bfm_ref_j_ - 1
               dfJ_bfm_ref_      = dfJ_bfm_lv0_ / n_segs_bfm_ref_j_
            enddo
         endif
      endif
      n_pts_bfm_ref_i_ = n_segs_bfm_ref_i_ - 1
      n_pts_bfm_ref_j_ = n_segs_bfm_ref_j_ - 1


      ! get refined (BFM) deltas (GRS)
      dfIi_bfm_ref_ = dfIi_bfm_lv0_ / n_segs_bfm_ref_i_
      dfIj_bfm_ref_ = dfIj_bfm_lv0_ / n_segs_bfm_ref_i_
      dfJi_bfm_ref_ = dfJi_bfm_lv0_ / n_segs_bfm_ref_j_
      dfJj_bfm_ref_ = dfJj_bfm_lv0_ / n_segs_bfm_ref_j_

      ! get BRM interpolated deltas (GRS) (taken from refined BFM deltas this time)
      dfIi_brm_interp = dfIi_bfm_ref_ / pol%interp_I_fct_
      dfIj_brm_interp = dfIj_bfm_ref_ / pol%interp_I_fct_
      dfJi_brm_interp = dfJi_bfm_ref_ / pol%interp_J_fct_
      dfJj_brm_interp = dfJj_bfm_ref_ / pol%interp_J_fct_

      ! get absolute deltas (in LRS), along I and J directions
      dfI_brm_interp_ = dfI_bfm_ref_ / pol%interp_I_fct_
      dfJ_brm_interp_ = dfJ_bfm_ref_ / pol%interp_J_fct_


#ifndef BSA_M3MF_ONLY_PREMESH_
      ! compute BFM (refined) influence areas for integration
      dwI = dfI_bfm_ref_ * CST_PIt2
      dwJ = dfJ_bfm_ref_ * CST_PIt2
      ctr_infl_bfm = dwI * dwJ
      brd_infl_bfm = ctr_infl / 2._RDP
      vtx_infl_bfm = brd_infl / 2._RDP
#endif

      ! compute BRM influence areas for integration
      dwI = dfI_brm_interp_ * CST_PIt2
      dwJ = dfJ_brm_interp_ * CST_PIt2
      ctr_infl = dwI * dwJ
      brd_infl = ctr_infl / 2._RDP
      vtx_infl = brd_infl / 2._RDP

      ! get actualised BFM-refined and BRM-interp  refinements (along borders)
      ni_bfm_ref_ = (this%ni_ - 1)
      nj_bfm_ref_ = (this%nj_ - 1)
      ni          = ni_bfm_ref_ * (n_segs_bfm_ref_i_ * pol%interp_I_fct_) + 1
      nj          = nj_bfm_ref_ * (n_segs_bfm_ref_j_ * pol%interp_J_fct_) + 1
      ni_bfm_ref_ = ni_bfm_ref_ * n_segs_bfm_ref_i_ + 1
      nj_bfm_ref_ = nj_bfm_ref_ * n_segs_bfm_ref_j_ + 1
      
      ! number of BRM points to interpolate (insert)
      ! between two know BFM (refined) points' direction lines.
      nipI = pol%interp_I_fct_ - 1
      nipJ = pol%interp_J_fct_ - 1
      i_brm_offsetJ = nipI * nj

      ! allocate data
      nPtsPost = ni * nj
      allocate(brm(dimM_bisp_, nPtsPost), stat=ist)
      if (ist /= 0) call bsa_Abort("Error allocating ""brm"" in interpolating RZ.")
      brm = 0._RDP

      allocate(bfm_new_left(dimM_bisp_, nj), stat=ist)
      if (ist /= 0) call bsa_Abort("Error allocating ""bfm_new_left"" in interpolating RZ.")
      bfm_new_left  = 0._RDP

      allocate(bfm_new_right(dimM_bisp_, nj), stat=ist)
      if (ist /= 0) call bsa_Abort("Error allocating ""bfm_new_right"" in interpolating RZ.")
      bfm_new_right = 0._RDP
      
      allocate(bfm_interp(dimM_bisp_, nj), stat=ist)
      if (ist /= 0) call bsa_Abort("Error allocating ""bfm_interp"" in interpolating RZ.")
      bfm_interp = 0._RDP


#ifdef __BSA_OMP
      allocate(fi_v_(nPtsPost), stat=ist)
      if (ist /= 0) call bsa_Abort("Error allocating ""fi"" in interpolating RZ.")
      fi = 0._RDP

      allocate(fj_v_(nPtsPost), stat=ist)
      if (ist /= 0) call bsa_Abort("Error allocating ""fj"" in interpolating RZ.")
      fj = 0._RDP

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
            ! class default
            !    brm_export_data_ => null()  ! NOTE: produces "error #8201: Associate name cannot be a pointer."
         end select
      endif

#else

! # define __FREQ_I_
! # define __FREQ_J_
# define __FREQ_I_shift_
# define __FREQ_J_shift_
# define __PDATA ,null()

#endif


! #ifdef __BSA_DEBUG
!       print *, nPtsPost, ni_bfm_ref_ * nj_bfm_ref_, this%ni_ * this%nj_
! #endif


      ! before starting, interpolate very first column
      ! along J dir, to get new mesh from old
      ! Then for all the others, it will be done inside
      ! the loop over the NI (old) points of the old mesh.
      ! NOTE: integrate as well.
      

      ! init point vertex
      i_brm = 1
      fi_baseptI = this%Ipt_%freqI()
      fj_baseptI = this%Ipt_%freqJ()
      fi         = fi_baseptI
      fj         = fj_baseptI
      fi_baseptJ = fi
      fj_baseptJ = fj

      bfmtail            = bfm_undump(:, 1)
#ifndef BSA_M3MF_ONLY_PREMESH_
      intg_bfm           = bfmtail * vtx_infl_bfm
#endif
      bfm_new_left(:, 1) = bfmtail
      
      brm(:, 1) = getBRM_msh(bfmtail, fi, fj)
#ifdef __BSA_OMP
      __FREQ_I_
      __FREQ_J_
#else
      call write_brm_fptr_(fi, fj, brm(:, 1)  __PDATA)
#endif
      intg = brm(:, 1) * vtx_infl

      do pJhead = 2, this%nj_ ! loop on all OLD BFM saved points (J-dir)

         do i_bfm_ref_j = 1, n_pts_bfm_ref_j_ ! loop on all REF BFM pts between 2 old.
            
            ! compute head
            fi_baseptJ = fi_baseptJ + dfJi_bfm_ref_
            fj_baseptJ = fj_baseptJ + dfJj_bfm_ref_
            bfmhead    = getBFM_msh(fi_baseptJ, fj_baseptJ)
#ifndef BSA_M3MF_ONLY_PREMESH_
            intg_bfm   = intg_bfm + bfmhead * brd_infl_bfm
#endif

            dfJhead = dfJ_bfm_ref_
            dfJtail = 0._RDP
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
               intg = intg + brm(:, i_brm) * brd_infl ! NOTE: it is a border point
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

            i_brm                  = i_brm + 1
            bfm_new_left(:, i_brm) = bfmhead
            brm(:, i_brm) = getBRM_msh(bfm_new_left(:, i_brm), fi, fj)
#ifdef __BSA_OMP
            __FREQ_I_
            __FREQ_J_
#else
            call write_brm_fptr_(fi, fj, brm(:, i_brm), null())
#endif
            ! NOTE: it is a border point, except for the very last one (VERTEX)
            intg = intg + brm(:, i_brm) * brd_infl

            bfmtail = bfmhead
         enddo ! n. of (exact) ref points for BFM

         !
         ! NOTE: now new head is OLD BFM (next) point!
         !
         bfmhead  = bfm_undump(:, pJhead) ! OK because we stored it NJ majour.
#ifndef BSA_M3MF_ONLY_PREMESH_
         intg_bfm = intg_bfm + bfmhead * brd_infl_bfm
#endif
         

#ifdef __BSA_DEBUG
         if (n_pts_bfm_ref_j_ > 0) then 
            if (abs((fi_baseptI + (dfJi_bfm_lv0_*(pJhead-1))) - (fi_baseptJ + dfJi_bfm_ref_)) > MACHINE_PRECISION .or. &
                  abs((fj_baseptI + (dfJj_bfm_lv0_*(pJhead-1))) - (fj_baseptJ + dfJj_bfm_ref_)) > MACHINE_PRECISION) &
                     call bsa_Abort(&
                        'BAD (2): fi or fj at the end of a BFM ref segment does not coincide..')
         endif
#endif

         
         dfJhead = dfJ_bfm_ref_
         dfJtail = 0._RDP
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
            intg = intg + brm(:, i_brm) * brd_infl ! NOTE: it is a border point
         enddo ! pJtail = 1, nipJ

         ! here, treat head, TAIL==HEAD (head - old mesh)
         fi = fi + dfJi_brm_interp
         fj = fj + dfJj_brm_interp
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
         intg = intg + brm(:, i_brm) * brd_infl

         ! old head (old mesh) becomes new tail
         ! TODO: we could change bfmhead here, so that
         !       it might be already ready for next loop.
         bfmtail = bfmhead

         ! NOTE: set new bfm ref base freqs as (current) head (old-mesh) point
         fi_baseptJ = fi
         fj_baseptJ = fj
      enddo

      ! NOTE: removing excess contribution for last HEAD (VERTEX)
      intg     = intg - brm(:, i_brm) * vtx_infl
#ifndef BSA_M3MF_ONLY_PREMESH_
      intg_bfm = intg_bfm - bfmtail * vtx_infl_bfm
#endif


      !
      i_bfm_old = this%nj_
      do pIcurr = 2, this%ni_ ! loop on all OLD BFM infl lines (I-dir)

         ! before doing any computation,
         ! we need to interpolate BFM along J
         ! at new CURRENT (I) infl line (including ref infl lines)
         ! NOTE: once we go through, integrate as well.

         do i_bfm_ref_i = 1, n_pts_bfm_ref_i_ ! loop on all REF BFM pts between 2 old (I-dir)

            ! computing BRM offset from pi_prev and pi_curr J infl lines.
            i_brm_shift  = i_brm + i_brm_offsetJ
            i_brm_write_ = i_brm_shift

            ! reset base freqs to point to new base -> prev base moved by ref BFM deltas (I-dir)
            ! NOTE: still keep prev base in memory here since they might serve later.
            fi_baseptJ = fi_baseptI + dfIi_bfm_ref_  ! reset J bases to match next I
            fj_baseptJ = fj_baseptI + dfIj_bfm_ref_
            fi         = fi_baseptJ
            fj         = fj_baseptJ

            bfmtail             = getBFM_msh(fi, fj)
#ifndef BSA_M3MF_ONLY_PREMESH_
            intg_bfm            = intg_bfm + bfmtail * brd_infl_bfm
#endif
            bfm_new_right(:, 1) = bfmtail

            i_brm_shift         = i_brm_shift + 1
            brm(:, i_brm_shift) = getBRM_msh(bfm_new_right(:, 1), fi, fj)

            ! call write_brm_fptr_(fi, fj, brm(:, i_brm_shift), null())
            __FREQ_I_shift_
            __FREQ_J_shift_

            intg = intg + brm(:, i_brm_shift) * brd_infl

            i_bfm_interpJ = 1
            do pJhead = 2, this%nj_ ! loop on all OLD BFM saved points (J-dir)

               do i_bfm_ref_j = 1, n_pts_bfm_ref_j_ ! loop on all REF BFM pts between 2 old.

                  fi_baseptJ = fi_baseptJ + dfJi_bfm_ref_
                  fj_baseptJ = fj_baseptJ + dfJj_bfm_ref_
                  bfmhead    = getBFM_msh(fi_baseptJ, fj_baseptJ)
#ifndef BSA_M3MF_ONLY_PREMESH_
                  intg_bfm   = intg_bfm + bfmhead * ctr_infl_bfm
#endif
                  
                  ! once we moved head, restore init, distances from head/tail
                  dfJhead = dfJ_bfm_ref_
                  dfJtail = 0._RDP
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
                     intg  = intg + brm(:, i_brm_shift) * ctr_infl
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
                  intg = intg + brm(:, i_brm_shift) * ctr_infl

                  bfmtail = bfmhead
               enddo ! n. of (exact) ref points for BFM (J-dir)
      
               ! here, next head is special (lies on an OLD BFM I-dir infl line).
               ! NOTE: don't forget to interpolate between this head and tail!!
               ! NOTE: it is a center point, except for the very last one (BORDER)
               fi_baseptJ = fi_baseptJ + dfJi_bfm_ref_
               fj_baseptJ = fj_baseptJ + dfJj_bfm_ref_
               bfmhead    = getBFM_msh(fi_baseptJ, fj_baseptJ)
#ifndef BSA_M3MF_ONLY_PREMESH_
               intg_bfm   = intg_bfm + bfmhead * vtx_infl_bfm
#endif

               ! once we moved head, restore init, distances from head/tail
               dfJhead = dfJ_bfm_ref_
               dfJtail = 0._RDP
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
                  intg  = intg + brm(:, i_brm_shift) * ctr_infl
               enddo ! pJtail = 1, nipJ

               ! here treat this new head
               fi = fi + dfJi_brm_interp  ! they should equate fi_baseptJ
               fj = fj + dfJj_brm_interp
               i_bfm_interpJ                   = i_bfm_interpJ + 1
               bfm_new_right(:, i_bfm_interpJ) = bfmhead
               i_brm_shift         = i_brm_shift + 1
               brm(:, i_brm_shift) = getBRM_msh(bfm_new_right(:, i_bfm_interpJ), fi, fj)
               ! call write_brm_fptr_(fi, fj, brm(:, i_brm_shift), null())
               __FREQ_I_shift_
               __FREQ_J_shift_

               intg = intg + brm(:, i_brm_shift) * ctr_infl

               bfmtail = bfmhead ! old head becomes new tail
            enddo ! pJhead = 2, this%nj_

            ! removing excess of very last HEAD, accounted as center, it is BORDER.
            ! NOTE: even worse for very last HEAD which happens to be End point.
            !       there, it is a VERTEX point.
            !       However, it's the very last element in brm, we can remove it after.
            intg     = intg - brm(:, i_brm_shift) * brd_infl
#ifndef BSA_M3MF_ONLY_PREMESH_
            intg_bfm = intg_bfm - bfmtail * brd_infl_bfm
#endif

            ! Now INTERPOLATE along I-dir between left and right BFM infl lines.
            !
            dfIcurr = dfI_bfm_ref_ ! reset I-dir CURR-PREV distances
            dfIprev = 0._RDP
            do pIprev = 1, nipI ! interp (I-dir) between prev-curr

               ! bulk I-dir interpolation until pj_head section level.
               ! Then, after treat that triang shaped zone separately.

               dfIprev = dfIprev + dfI_brm_interp_
               dfIcurr = dfIcurr - dfI_brm_interp_
               
               bfm_interp = &
                  (  bfm_new_left  * dfIcurr + &
                     bfm_new_right * dfIprev ) / dfI_bfm_ref_

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
               intg = intg + brm(:, i_brm) * brd_infl
               
               do pJtail = 2, nj

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
                  intg = intg + brm(:, i_brm) * ctr_infl
               enddo

               ! removing excess from having accounted values at last iter HEAD as
               ! center points (they are BORDER)
               intg = intg - brm(:, i_brm) * brd_infl

            enddo ! pIprev = 1, nipI


            ! now update bases along I (CURR now, PREV next iteration!)
            fi_baseptI = fi_baseptI + dfIi_bfm_ref_
            fj_baseptI = fj_baseptI + dfIj_bfm_ref_


#ifndef __BSA_OMP
            ! Now, we can write actual new BFM I-dir infl line.
            fi = fi_baseptI
            fj = fj_baseptI
            do pIprev = 1, nj
               i_brm_write_ = i_brm_write_ + 1
               call write_brm_fptr_(fi, fj, brm(:, i_brm_write_), null())
               fi = fi + dfJi_brm_interp
               fj = fj + dfJj_brm_interp
            enddo
#endif

            ! moving right, shift infl-lines.
            bfm_new_left = bfm_new_right

            i_brm = i_brm_shift

         enddo ! BFM ref points along I dir



         !
         ! Here, right BFM infl line is one where we have old BFM mesh points !
         !

         ! computing BRM offset from pi_prev and pi_curr J infl lines.
         i_brm_shift  = i_brm + i_brm_offsetJ
         i_brm_write_ = i_brm_shift

         ! again, I bases refer to PREV infl line.
         fi_baseptJ = fi_baseptI + dfIi_bfm_ref_
         fj_baseptJ = fj_baseptI + dfIj_bfm_ref_
         fi         = fi_baseptJ
         fj         = fj_baseptJ

         i_bfm_old           = i_bfm_old + 1
         bfmtail             = bfm_undump(:, i_bfm_old)
#ifndef BSA_M3MF_ONLY_PREMESH_
         intg_bfm            = intg_bfm + bfmtail * brd_infl_bfm
#endif
         bfm_new_right(:, 1) = bfmtail
         
         i_brm_shift         = i_brm_shift + 1
         brm(:, i_brm_shift) = getBRM_msh(bfm_new_right(:, 1), fi, fj)
         ! call write_brm_fptr_(fi, fj, brm(:, i_brm_shift), null())
         __FREQ_I_shift_
         __FREQ_J_shift_

         intg = intg + brm(:, i_brm_shift) * brd_infl

         ! compute right infl-line
         i_bfm_interpJ = 1
         do pJhead = 2, this%nj_

            do i_bfm_ref_j = 1, n_pts_bfm_ref_j_ ! loop on all REF BFM pts between 2 old.

               fi_baseptJ = fi_baseptJ + dfJi_bfm_ref_
               fj_baseptJ = fj_baseptJ + dfJj_bfm_ref_
               bfmhead    = getBFM_msh(fi_baseptJ, fj_baseptJ)
#ifndef BSA_M3MF_ONLY_PREMESH_
               intg_bfm   = intg_bfm + bfmhead * ctr_infl_bfm
#endif
      
               ! once we moved head, restore init, distances from head/tail
               dfJhead = dfJ_bfm_ref_
               dfJtail = 0._RDP
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
                  intg  = intg + brm(:, i_brm_shift) * ctr_infl
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
               intg = intg + brm(:, i_brm_shift) * ctr_infl

               bfmtail = bfmhead
            enddo

            ! next head is OLD BFM point
            fi = fi_baseptJ
            fj = fj_baseptJ
            fi_baseptJ = fi_baseptJ + dfJi_bfm_ref_
            fj_baseptJ = fj_baseptJ + dfJj_bfm_ref_
            
            i_bfm_old = i_bfm_old + 1
            bfmhead   = bfm_undump(:, i_bfm_old)
#ifndef BSA_M3MF_ONLY_PREMESH_
            intg_bfm  = intg_bfm + bfmhead * ctr_infl_bfm
#endif

            ! once we moved head, restore init distances from head/tail
            dfJhead = dfJ_bfm_ref_
            dfJtail = 0._RDP
            do pJtail = 1, nipJ
   
               fi = fi + dfJi_brm_interp
               fj = fj + dfJj_brm_interp
   
               ! update actual distances head/tail
               dfJtail = dfJtail + dfJ_brm_interp_
               dfJhead = dfJhead - dfJ_brm_interp_
   
               ! interpolation along J dir (between HEAD-TAIL)
               ! NOTE: save it for later use.
               i_bfm_interpJ = i_bfm_interpJ + 1
               bfm_new_right(:, i_bfm_interpJ) = &
                  (bfmhead * dfJtail + bfmtail * dfJhead) / dfJ_bfm_ref_

               i_brm_shift         = i_brm_shift + 1
               brm(:, i_brm_shift) = getBRM_msh(bfm_new_right(:, i_bfm_interpJ), fi, fj)
               ! call write_brm_fptr_(fi, fj, brm(:, i_brm_shift), null())
               __FREQ_I_shift_
               __FREQ_J_shift_
   
               ! NOTE: it is a center point, except for very last row -> BORDER
               intg  = intg + brm(:, i_brm_shift) * ctr_infl
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
            intg = intg + brm(:, i_brm_shift) * ctr_infl

            ! old head becomes new tail
            bfmtail = bfmhead
         enddo ! pJhead = 2, this%nj_


         ! removing excess of very last HEAD
         ! accounted as center, it is BORDER
         ! NOTE: even worse for very last HEAD which happens to be End point.
         !       there, it is a VERTEX point.
         !       However, it's the very last element in brm, we can remove it after.
         intg     = intg - brm(:, i_brm_shift) * brd_infl
#ifndef BSA_M3MF_ONLY_PREMESH_
         intg_bfm = intg_bfm - bfmtail * brd_infl_bfm
#endif


         ! ok, here we now have BFM values (interpolated along J)
         ! at CURR and PREV (I) index pointers.
         ! We have to interpolate along I between CURR and PREV, i.e.
         ! prev has to start moving toward CURR.


         dfIcurr = dfI_bfm_ref_ ! reset I-dir CURR-PREV distances
         dfIprev = 0._RDP
         do pIprev = 1, nipI ! interpolate along I

            ! bulk I-dir interpolation until pj_head section level.
            ! Then, after treat that triang shaped zone separately.

            dfIprev = dfIprev + dfI_brm_interp_
            dfIcurr = dfIcurr - dfI_brm_interp_
            
            bfm_interp = &
               (  bfm_new_left  * dfIcurr + &
                  bfm_new_right * dfIprev ) / dfI_bfm_ref_


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
            intg = intg + brm(:, i_brm) * brd_infl
            
            do pJtail = 2, nj

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
               intg = intg + brm(:, i_brm) * ctr_infl
            enddo

            ! removing excess from having accounted values at last iter HEAD as
            ! center points (they are BORDER)
            intg = intg - brm(:, i_brm) * brd_infl

         enddo ! pIprev = 1, nipI
         
         ! once finished with this section (CURR-PREV), since we skip J column
         ! at pi_curr infl line, we reset general BRM index to point to
         ! previously shifted one.
         i_brm = i_brm_shift

         ! now update bases along I
         fi_baseptI = fi_baseptI + dfIi_bfm_ref_
         fj_baseptI = fj_baseptI + dfIj_bfm_ref_


#ifndef __BSA_OMP
         ! Now, we can write actual new BFM I-dir infl line.
         fi = fi_baseptI
         fj = fj_baseptI
         do pIprev = 1, nj
            i_brm_write_ = i_brm_write_ + 1
            call write_brm_fptr_(fi, fj, brm(:, i_brm_write_), null())
            fi = fi + dfJi_brm_interp
            fj = fj + dfJj_brm_interp
         enddo
#endif

         
         bfm_new_left = bfm_new_right

      enddo ! pIcurr = 2, this%ni_


! #ifdef __BSA_DEBUG
      if (i_brm /= nPtsPost) then
         print *, ' policy BFM-MLV factors      =  ', pol%interp_bfm_I_fct_, pol%interp_bfm_J_fct_
         print *, ' policy   BRM   factors      =  ', pol%interp_I_fct_, pol%interp_J_fct_
         print *, ' ibrm, nPtsPost  = ', i_brm, nPtsPost
         call bsa_Abort('"i_brm" does not equal rect zone''s n. of (interpolated) points.')
      endif
! #endif


      deallocate(bfm_new_left)
      deallocate(bfm_new_right)
      deallocate(bfm_interp)


      ! removing overestimation for B point
      intg = intg - brm(:, i_brm - nj + 1) * vtx_infl

      ! removing overestimation for last border points
      intg = intg - sum(brm(:, i_brm - nj + 2 : i_brm - 1) * brd_infl, dim=2)

      ! removing overestimation for End point
      intg = intg - brm(:, i_brm) * vtx_infl


      !$omp critical
#ifdef __BSA_OMP
      call write_brm_fptr_(fi_v_, fj_v_, brm, pdata)
#endif

      ! updating main point counters 
      msh_bfmpts_post_     = msh_bfmpts_post_ + ni_bfm_ref_ * nj_bfm_ref_
      msh_brmpts_post_     = msh_brmpts_post_ + nPtsPost

#ifndef BSA_M3MF_ONLY_PREMESH_
      m3mf_msh_ptr_ = m3mf_msh_ptr_ + (intg_bfm * settings%i_bisp_sym_) ! update main BFM integral
#endif
      m3mr_msh_ptr_ = m3mr_msh_ptr_ + (intg     * settings%i_bisp_sym_) ! update main BRM integral
      !$omp end critical
   end subroutine interpolateRZ_HTPC_v3










!    !> Implementation of HTPC interpolation method for a rectangle zone.
!    !> This method IS FASTER (AND MORE ACCURATE!!)
!    subroutine interpolateRZ_HTPC_v2(this)
! #ifdef __BSA_OMP
!       use BsaLib_Data, only: dimM_bisp_, getBRM_msh, m3mr_msh_ptr_
! #else
!       use BsaLib_Data, only: dimM_bisp_, getBRM_msh, bfm_undump, m3mr_msh_ptr_
! #endif
!       use BsaLib_MPolicy, only: MPolicy_t
!       class(MRectZone_t), intent(inout) :: this

!       type(MPolicy_t) :: pol
!       integer   :: ni, nj
!       integer   :: nipI, nipJ, zNintrpPts
!       real(RDP) :: dfIi_old, dfIj_old, dfJi_old, dfJj_old
!       real(RDP) :: dfIi_interp, dfIj_interp, dfJi_interp, dfJj_interp
!       real(RDP) :: dfI_old, dfJ_old, dfI_interp, dfJ_interp, dwI, dwJ
!       real(RDP) :: vtx_infl, brd_infl, ctr_infl

!       ! HTPC indexes
!       integer :: picurr, piprev, pjhead, pjtail

!       !> Pos in general BFM undumped data
!       integer :: i_bfm
!       integer :: i_brm, i_brm_shift, ibrmoffset
!       integer :: i_bfm_interpJ

!       ! freqs
!       real(RDP) :: base_fi, base_fj, fi, fj


!       real(RDP), allocatable :: bfm_new_left(:, :), bfm_new_right(:, :)
!       real(RDP), allocatable :: bfm_interp(:, :)

!       real(RDP) :: dfJtail, dfJhead, dfIcurr, dfIprev
!       real(RDP) :: bfmtail(dimM_bisp_), bfmhead(dimM_bisp_)


!       real(RDP), allocatable :: brm(:, :)
!       real(RDP) :: intg(dimM_bisp_)

!       integer :: iall1, iall2, iall3, iall4


!       ! zone's policy
!       pol = this%policy()


!       ! get original (pre-meshing) deltas
!       ! NOTE: keep them in memory unchanged since they
!       !       might serve later on.
!       call this%getIJfsteps(dfIi_old, dfIj_old, dfJi_old, dfJj_old)


!       ! get interpolated deltas (GRS)
!       dfIi_interp = dfIi_old / pol%interp_I_fct_
!       dfIj_interp = dfIj_old / pol%interp_I_fct_
!       dfJi_interp = dfJi_old / pol%interp_J_fct_
!       dfJj_interp = dfJj_old / pol%interp_J_fct_


!       ! get absolute deltas (in LRS)
!       ! along I and J directions
!       dfI_old    = this%deltaf_I_
!       dfI_interp = dfI_old / pol%interp_I_fct_
      
!       dfJ_old    = this%deltaf_J_
!       dfJ_interp = dfJ_old / pol%interp_J_fct_


!       ! compute influence areas for integration
!       dwI = dfI_interp * CST_PIt2
!       dwJ = dfJ_interp * CST_PIt2
!       ctr_infl = dwI * dwJ
!       brd_infl = ctr_infl / 2._RDP
!       vtx_infl = brd_infl / 2._RDP

!       ! get actualised refinements (along borders)
!       ni = (this%ni_ - 1) * pol%interp_I_fct_ + 1
!       nj = (this%nj_ - 1) * pol%interp_J_fct_ + 1
      
!       ! number of points to interpolate (insert)
!       ! between two know points' direction lines
!       nipI = pol%interp_I_fct_ - 1
!       nipJ = pol%interp_J_fct_ - 1

!       ibrmoffset = nipI * nj

!       ! allocate data
!       zNintrpPts = ni * nj
!       allocate(brm(dimM_bisp_, zNintrpPts), stat=iall1)
!       brm = 0._RDP
!       allocate(bfm_new_left(dimM_bisp_, nj), stat=iall2)
!       bfm_new_left  = 0._RDP
!       allocate(bfm_new_right(dimM_bisp_, nj), stat=iall3)
!       bfm_new_right = 0._RDP
!       allocate(bfm_interp(dimM_bisp_, nj), stat=iall4)
!       bfm_interp = 0._RDP


!       ! before starting, interpolate very first column
!       ! along J dir, to get new mesh from old
!       ! Then for all the others, it will be done inside
!       ! the loop over the NI (old) points of the old mesh.
!       ! NOTE: integrate as well.
      
!       ! init point vertex
!       base_fi = this%Ipt_%freqI()
!       base_fj = this%Ipt_%freqJ()
!       fi      = base_fi
!       fj      = base_fj
!       bfmtail = bfm_undump(:, 1)
!       bfm_new_left(:, 1) = bfmtail
!       brm(:, 1) = getBRM_msh(bfm_new_left(:, 1), fi, fj)
!       intg      = brm(:, 1) * vtx_infl

!       i_brm = 1
!       do pjhead = 2, this%nj_

!          ! OK because we stored it NJ majour.
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
!             intg = intg + brm(:, i_brm) * brd_infl
!          enddo ! pjtail = 1, nipJ

!          ! here, treat head, TAIL==HEAD
!          fi = fi + dfJi_interp
!          fj = fj + dfJj_interp
!          i_brm = i_brm + 1
!          bfm_new_left(:, i_brm) = bfmhead
!          brm(:, i_brm) = getBRM_msh(bfm_new_left(:, i_brm), fi, fj)

!          ! NOTE: it is a border point, except for the very last one (VERTEX)
!          intg = intg + brm(:, i_brm) * brd_infl

!          ! old head becomes new tail
!          ! TODO: we could change bfmhead here, so that
!          !       it might be already ready for next loop.
!          bfmtail = bfmhead
!       enddo

!       ! NOTE: removing excess contribution for last HEAD (VERTEX)
!       intg  = intg - brm(:, i_brm) * vtx_infl


!       i_bfm = this%nj_

!       ! NOTE: starting from 2 because we have to take
!       !       infl lines in consequent groups of two.
!       do picurr = 2, this%ni_


!          ! before doing any computation,
!          ! we need to interpolate BFM along J
!          ! at new CURRENT (I) infl line
!          ! NOTE: once we go through, integrate as well.

!          ! computing BRM offset from pi_prev and pi_curr J infl lines.
!          i_brm_shift = i_brm + ibrmoffset

!          ! reset freqs to point to new base -> prev base moved by old deltas!
!          ! (now pi_prev-pj_tail)
!          ! NOTE: still keep prev base in memory here since they might serve later.
!          fi = base_fi + dfIi_old
!          fj = base_fj + dfIj_old

!          i_bfm   = i_bfm + 1
!          bfmtail = bfm_undump(:, i_bfm)
!          bfm_new_right(:, 1) = bfmtail
         
!          i_brm_shift = i_brm_shift + 1
!          brm(:, i_brm_shift) = getBRM_msh(bfm_new_right(:, 1), fi, fj)
!          intg = intg + brm(:, i_brm_shift) * brd_infl

!          i_bfm_interpJ = 1
!          do pjhead = 2, this%nj_

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
   
!                ! NOTE: it is a center point, except for very last row -> BORDER
!                intg  = intg + brm(:, i_brm_shift) * ctr_infl
!             enddo ! pjtail = 1, nipJ
   
!             ! here, treat head, TAIL==HEAD
!             fi = fi + dfJi_interp
!             fj = fj + dfJj_interp
!             i_bfm_interpJ = i_bfm_interpJ + 1
!             bfm_new_right(:, i_bfm_interpJ) = bfmhead
            
!             i_brm_shift = i_brm_shift + 1
!             brm(:, i_brm_shift) = getBRM_msh(bfm_new_right(:, i_bfm_interpJ), fi, fj)
   
!             ! NOTE: it is a center point, except for the very last one (BORDER)
!             intg = intg + brm(:, i_brm_shift) * ctr_infl
   
!             ! old head becomes new tail
!             bfmtail = bfmhead
            
!          enddo ! pjhead = 2, this%nj_

!          ! removing excess of very last HEAD
!          ! accounted as center, it is BORDER
!          ! NOTE: even worse for very last HEAD which happens to be End point.
!          !       there, it is a VERTEX point.
!          !       However, it's the very last element in brm, we can remove it after.
!          intg = intg - brm(:, i_brm_shift) * brd_infl




!          ! ok, here we now have BFM values (interpolated along J)
!          ! at CURR and PREV (I) index pointers.
!          ! We have to interpolate along I between CURR and PREV, i.e.
!          ! prev has to start moving toward CURR.


!          ! reset I-dir CURR-PREV distances
!          dfIcurr    = dfI_old
!          dfIprev    = 0._RDP

!          ! interpolate along I
!          do piprev = 1, nipI

!             ! bulk I-dir interpolation until pj_head section level.
!             ! Then, after treat that triang shaped zone separately.

!             dfIprev = dfIprev + dfI_interp
!             dfIcurr = dfIcurr - dfI_interp
            
!             bfm_interp = &
!                (  bfm_new_left  * dfIcurr + &
!                   bfm_new_right * dfIprev ) / dfI_old

!             ! once we have the values, go through them to integrate
!             ! NOTE: reset base freqs pointers, this time moving them along INTERP mesh
!             fi = base_fi + (dfIi_interp * piprev)
!             fj = base_fj + (dfIj_interp * piprev)

!             i_brm = i_brm + 1
!             brm(:, i_brm) = getBRM_msh(bfm_interp(:, 1), fi, fj)
!             intg = intg + brm(:, i_brm) * brd_infl
            
!             do pjtail = 2, nj

!                fi = fi + dfJi_interp
!                fj = fj + dfJj_interp
               
!                i_brm = i_brm + 1
!                brm(:, i_brm) = getBRM_msh(bfm_interp(:, pjtail), fi, fj)
!                intg = intg + brm(:, i_brm) * ctr_infl
!             enddo

!             ! removing excess from having accounted values at last iter HEAD as
!             ! center points (they are BORDER)
!             intg = intg - brm(:, i_brm) * brd_infl

!          enddo ! piprev = 1, nipI


!          ! now we can update (PREV) base freqs to match CURR
!          ! NOTE: CURR J infl line has already been computed!
!          base_fi = base_fi + dfIi_old
!          base_fj = base_fj + dfIj_old

         
!          bfm_new_left = bfm_new_right


!          ! once finished with this section (CURR-PREV), since we skip J column
!          ! at pi_curr infl line, we reset general BRM index to point to
!          ! previously shifted one.
!          i_brm = i_brm_shift

!       enddo ! picurr = 2, this%ni_



! #ifdef __BSA_DEBUG
!       if (i_brm /= zNintrpPts) &
!          call bsa_Abort('"i_bfm" does not equal rect zone''s n. of (interpolated) points.')
! #endif


!       ! removing overestimation for B point
!       intg = intg - brm(:, i_brm - nj + 1) * vtx_infl

!       ! removing overestimation for last border points
!       intg = intg - sum(brm(:, i_brm - nj + 2 : i_brm - 1) * brd_infl, dim=2)

!       ! removing overestimation for End point
!       intg = intg - brm(:, i_brm) * vtx_infl

!       ! updating main integral
!       m3mr_msh_ptr_ = m3mr_msh_ptr_ + intg

!    end subroutine interpolateRZ_HTPC_v2




end submodule