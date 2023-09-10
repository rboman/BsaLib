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
module BsaLib_MRectZone

#include "../../../precisions"

   use BsaLib_MPoint
   use BsaLib_M2DPolygZone
   implicit none
   private

   ! imported from <- BsaLib_M2DPolygZone <- BsaLib_MZone
   public :: msh_max_zone_NPts, MZone_ID

   ! Make it visible without needing to import Base class module 
   public :: UndumpZone

   type, public, extends(M2DPolygZone_t) :: MRectZone_t

      !> Init point (from where everything else is made reference)
      type(MPoint_t) :: Ipt_

      !> End point
      type(MPoint_t) :: Ept_

      !> Rect base along I-dir (x axis)
      real(RDP) :: base_I_ = 0._RDP

      !> Rect base along J-dir (y axis)
      real(RDP) :: base_J_ = 0._RDP

      !> Delta freq [Hz] along I-dir (x axis)
      real(RDP) :: deltaf_I_ = 0._RDP

      !> Delta freq [Hz] along J-dir (y axis)
      real(RDP) :: deltaf_J_ = 0._RDP


      logical, private :: refmts_set_ = .false.
      logical, private :: deltas_set_ = .false.

   contains
   
      procedure, pass :: baseI => baseI_rct
      procedure, pass :: baseJ => baseJ_rct
      procedure, pass :: deduceDeltas => deduceDeltas_rect

      procedure, pass :: getAPoint
      procedure, pass :: getBPoint
      procedure, pass :: setDeltas
      procedure, pass :: deduceRefinements
      generic, public :: defineFromDeltas => defineFromDeltas_refinements, defineFromDeltas_maxvalues
      procedure, pass, private :: defineFromDeltas_refinements, defineFromDeltas_maxvalues
      generic :: defineFromEndPtCoordAndBase => defineFromEndPtCoordAndBase_norm, defineFromEndPtCoordAndBase_forceDeltas
      procedure, pass, private :: defineFromEndPtCoordAndBase_norm, defineFromEndPtCoordAndBase_forceDeltas
      procedure, pass :: define
      procedure, pass, private :: setIEpts
      procedure, pass :: validateDeltas
      procedure, pass :: getOtherBase
      procedure, pass :: getIJfsteps
      procedure, pass :: reconstructZoneBaseMesh
      procedure, pass :: compute => compute_s
      procedure, pass, private :: compute_s
      procedure, pass :: getNthQuadVtx
      procedure, pass :: dump => dumpRZ
      procedure, pass :: undump => undumpRZ
      procedure, pass :: interpolate => interpolateRZ
   end type MRectZone_t



   interface MRectZone
      module procedure MRectZone_t_custom_constructor
   end interface
   public :: MRectZone



   ! main interface to module procedures
   interface
      
      module function MRectZone_t_custom_constructor(rot, name) result(this)
         real(RDP), intent(in), optional :: rot
         character(len=*), intent(in), optional :: name
         ! NOTE: compiler uses default initialisation here (built-in)
         type(MRectZone_t) :: this
      end function


      !> Gets rect base along I-dir
      module elemental function baseI_rct(this) result(res)
         class(MRectZone_t), intent(in) :: this
         real(RDP) :: res
      end function


      !> Gets rect base along J-dir
      module elemental function baseJ_rct(this) result(res)
         class(MRectZone_t), intent(in) :: this
         real(RDP) :: res
      end function



      !> Get A point.
      !> A point is the point defined, starting from I point, 
      !> along the J-dir (Y-axis) parallel side.
      module pure function getAPoint(this) result(pt)
         class(MRectZone_t), intent(in) :: this
         type(MPoint_t) :: pt
      end function


      !> Get B point.
      !> B point is the point defined, starting from I point, 
      !> along the I-dir (X-axis) parallel side.
      module pure function getBPoint(this) result(pt)
         class(MRectZone_t), intent(in) :: this
         type(MPoint_t) :: pt
      end function


      !> If refinements are set, deltas are deduced
      module subroutine deduceDeltas_rect(this)
         class(MRectZone_t), intent(inout) :: this
      end subroutine


      !> Set frequency deltas
      module subroutine setDeltas(this, dfi, dfj, adapt)
         class(MRectZone_t), intent(inout) :: this
         real(RDP), intent(in)             :: dfi, dfj
         logical, intent(in), optional     :: adapt
      end subroutine


      !> If deltas are set, deduces refinements
      module subroutine deduceRefinements(this, adapt)
         class(MRectZone_t), intent(inout) :: this
         logical, intent(in) :: adapt
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
      end subroutine





      !> Defines zone bases (sides) from given Deltas, specifying max deltas values.
      !> NOTE: Init and End zone points MUST be known.
      !> If force is present and true, forces deltas to reach max values specified
      !> by 'valI' and 'valJ'. In this case, deltas are then slightly adjusted to fit.
      module subroutine defineFromDeltas_maxvalues(this, pt, loc, dfi, dfj, maxF_i, maxF_j, force, exceed)
         class(MRectZone_t), intent(inout) :: this

         !> Point specification, which location is specified by 'loc'
         class(MPoint_t), intent(in)  :: pt

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
      end subroutine





      module subroutine define(this, pt, loc, base_i, base_j)
         class(MRectZone_t), intent(inout) :: this
         class(MPoint_t), intent(in) :: pt
         character(len=1), optional, intent(in) :: loc
         real(RDP), intent(in), optional  :: base_i, base_j
      end subroutine






      module subroutine defineFromEndPtCoordAndBase_norm(&
         this, Pi, coord_val, coord_ty_ch, baseval, base_dir, called)

         class(MRectZone_t), intent(inout) :: this
         class(MPoint_t), intent(in) :: Pi
         real(RDP), intent(in) :: coord_val
         character(len = 1), intent(in) :: coord_ty_ch
         real(RDP), intent(in) :: baseval
         character(len = 1), intent(in) :: base_dir
         logical, intent(in)            :: called
      end subroutine



      module subroutine defineFromEndPtCoordAndBase_forceDeltas(&
         this, Pi, coord_val, coord_ty_ch, baseval, base_dir, dfi, dfj)

         class(MRectZone_t), intent(inout) :: this
         class(MPoint_t), intent(in) :: Pi
         real(RDP), intent(in) :: coord_val
         character(len = 1), intent(in) :: coord_ty_ch
         real(RDP), intent(in) :: baseval
         character(len = 1), intent(in) :: base_dir
         real(RDP), intent(in)          :: dfi, dfj
      end subroutine







      !> Sets, based on passed point location
      !> Init and End zone points.
      module subroutine setIEpts(this, pt, loc)
         class(MRectZone_t), intent(inout) :: this
         class(MPoint_t), intent(in)  :: pt
         character(len=1), intent(in) :: loc
      end subroutine



      !> Avoid setting a delta smaller than given limit
      module elemental impure subroutine validateDeltas(this, lval)
         class(MRectZone_t), intent(inout) :: this
         real(RDP), intent(in) :: lval
      end subroutine


      !> Automatically computes the second remaining (unknown) rect base
      !> based on the point's coordinates and the base that we already defined.
      !> BUG: cannot see when a base is negative (i.e. when END pt is "behind" INIT one).
      module subroutine getOtherBase(this, pt, base_dir, known_coord, coord_val)
         class(MRectZone_t), intent(inout) :: this
         class(MPoint_t), intent(in)  :: pt
         character(len=1), intent(in) :: base_dir, known_coord
         real(RDP), intent(in) :: coord_val
      end subroutine



      !> Gets actualised frequency deltas along two main
      !> sides directions (I, J), actualised based on *this
      !> zone rotation w.r.t. GRS.
      module subroutine getIJfsteps(this, dfIx, dfIy, dfJx, dfJy)
         class(MRectZone_t), intent(in) :: this
         real(RDP), intent(out) :: dfIx, dfIy, dfJx, dfJy
      end subroutine



      !> Returns the whole zone reconstructed mesh
      !> Maybe used for visualising.
      module function reconstructZoneBaseMesh(this) result(msh)
         class(MRectZone_t), intent(in) :: this
         !> BUG: might be 2-rank array instead of 3!
         real(RDP) :: msh(2, this%nj_, this%ni_)
      end function


      !> Actual zone comutation (pre phase).
      module subroutine compute_s(this)
         class(MRectZone_t), intent(inout) :: this
      end subroutine



      !> Gets vertex point pt coordinates in Nth quadrant
      !> w.r.t. Center point, in zone's LRS.
      module pure function getNthQuadVtx(this, iquad) result(pt)
         class(MRectZone_t), intent(in) :: this
         integer, intent(in) :: iquad
         type(MPoint_t) :: pt
      end function





      !> Dumps a RECTANGULAR zone.
      !>
      !> NOTE: Each specific zone dumping method is called 
      !>       from the STATIC MZone_t procedure DumpZone().
      module subroutine dumpRZ(this)
         class(MRectZone_t), intent(in) :: this
      end subroutine



      !> Undumps a RECTANGULAR zone.
      !>
      !> NOTE: Each specific zone dumping method is called 
      !>       from the STATIC MZone_t procedure UndumpZone().
      !> NOTE: This is the equivalent as reconstruct() of MATLAB.
      module subroutine undumpRZ(this)
         class(MRectZone_t), intent(inout) :: this
      end subroutine



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
#endif
      end subroutine


   end interface


end module BsaLib_MRectZone