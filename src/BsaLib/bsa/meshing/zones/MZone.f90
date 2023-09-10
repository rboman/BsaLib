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
module BsaLib_MZone

#include "../../../precisions"

   use BsaLib_MPolicy
   use BsaLib_IO, only: unit_dump_bfm_, unit_debug_, undebug_fname_
   use BsaLib_CONSTANTS
   use BsaLib_Data, only: bsa_Abort, test_no_bfm_mlr_
   implicit none
   private

   public :: DefaultInitBaseZone, DumpZone, UndumpZone

   ! BUG: make them public to child MZone classes. Not optimal!!
   public :: unit_debug_, unit_dump_bfm_

   ! BUG: same as before, export constants values
   public :: CST_2d3, CST_3d2, CST_PId2
   public :: CST_PIGREC, CST_PIt2, CST_PIt3d2, CST_PIt4


   !> Tracks zone with max N. of points
   integer(kind = 4), public :: msh_max_zone_NPts = 0

   type, public :: MZoneEnum_t
      integer(kind = 4) :: NULL      = 0
      integer(kind = 4) :: RECTANGLE = 1
      integer(kind = 4) :: TRIANGLE  = 2
      integer(kind = 4) :: LINE      = 3
   end type MZoneEnum_t
   type(MZoneEnum_t), public, parameter :: MZone_ID = MZoneEnum_t()
   ! integer, public, parameter :: MZone_RECTANGLE = 1
   ! integer, public, parameter :: MZone_TRIANGLE  = 2
   ! integer, public, parameter :: MZone_LINEAR    = 3

   
   type, abstract, public :: MZone_t

      character(len = 64) :: name_   = ''
      type(MPolicy_t)     :: policy_

      !> Pointer to index of zone's interest modes
      integer(kind = 4), public :: id_im_
   
   contains
      procedure, pass :: zoneName
      procedure, pass :: setPolicy
      procedure, pass :: policy
      procedure, pass :: setInterestModeIndexPtr
      procedure, pass :: disableZonePolicyBfmMLR
      procedure(intf_MZoneIntFct_),    pass, deferred :: zoneTotNPts
      procedure(intf_MZoneGenIn_),     pass, deferred :: dump
      procedure(intf_MZoneGenInOut_),  pass, deferred :: undump
#ifdef __BSA_OMP
      procedure(intf_MZoneInterpOMP_), pass, deferred :: interpolate
#else
      procedure(intf_MZoneGenInOut_),  pass, deferred :: interpolate
#endif
   end type MZone_t



   abstract interface
      pure function intf_MZoneIntFct_(this) result(res)
         import MZone_t
         class(MZone_t), intent(in) :: this
         integer :: res
      end function

      subroutine intf_MZoneGenIn_(this)
         import MZone_t
         class(MZone_t), intent(in) :: this
      end subroutine

      subroutine intf_MZoneGenInOut_(this)
         import MZone_t
         class(MZone_t), intent(inout) :: this
      end subroutine

#ifdef __BSA_OMP
      subroutine intf_MZoneInterpOMP_(this, bfm, pdata)
         import MZone_t, RDP
         class(MZone_t), intent(inout) :: this
         real(RDP), intent(in)         :: bfm(:, :)
         class(*), pointer, intent(in) :: pdata
      end subroutine
#endif
   end interface



contains



   subroutine DefaultInitBaseZone(this)
      class(MZone_t), intent(inout) :: this

      this%policy_ = MPolicy_NULL
   end subroutine


   subroutine zoneName(this, name_in)
      class(MZone_t), intent(inout) :: this
      character(len=*), intent(in)  :: name_in

      this%name_ = name_in(1:len_trim(name_in))
   end subroutine




   subroutine setInterestModeIndexPtr(this, id)
      class(MZone_t), intent(inout) :: this
      integer(kind = 4), intent(in) :: id

      this%id_im_ = id
   end subroutine




   subroutine setPolicy(this, var_in)
      class(MZone_t), intent(inout)  :: this
      class(*), intent(in) :: var_in
      select type (var_in)
         class is (MPolicy_t)
            this%policy_ = var_in
         type is (integer)
            this%policy_ = var_in
         class default
            call bsa_Abort('Unsupported type. Must be either "integer" or "MPolicy_t".')
      end select
   end subroutine

   function policy(this) result(pol_out)
      class(MZone_t), intent(inout)  :: this
      type(MPolicy_t) :: pol_out

      pol_out = this%policy_
   end function






   
   subroutine DumpZone(z, data)
      class(MZone_t), intent(in) :: z
      real(RDP), intent(in)      :: data(:, :)
      ! integer             :: tot

      ! dump specific zone data
      ! NOTE: keep this first since 
      !       we want to read as first parameter,
      !       the actual zone type identifier, so that
      !       we can directly specialise undumping in Mesh()
      !       routine.
      call z%dump()

      ! write common zone data
      write(unit_dump_bfm_) z%name_

      ! policy
      write(unit_dump_bfm_) z%policy_%getID()

      ! zone interest modes index ptr
      write(unit_dump_bfm_) z%id_im_


      ! Dump BFM data.
      !
      ! ! write how many bytes in total to be read
      ! ! afterwards. Then, dimBISP is automatically
      ! ! deferred knowing num of zone's meshing points
      ! tot = size(data)
      ! write(unit_dump_bfm_) tot
      
      write(unit_dump_bfm_) data ! NOTE: dimBISP first, then nj * ni
   end subroutine DumpZone



#ifdef __BSA_OMP
   subroutine UndumpZone(z, bfm_undump)
      use BsaLib_Data, only: dimM_bisp_
      real(RDP), allocatable, intent(inout) :: bfm_undump(:, :)
#else
   subroutine UndumpZone(z)
      use BsaLib_Data, only: bfm_undump
#endif
      class(MZone_t), intent(inout) :: z
      character(len = 64) :: name_hdr
      integer             :: zNp

      call z%undump()  ! read zone's specific data first

      ! read zone common data
      read(unit_dump_bfm_) name_hdr
      call z%zoneName(name_hdr(1:len_trim(name_hdr)))

      ! policy (ID)
      read(unit_dump_bfm_) zNp
      call z%setPolicy(zNp)
      if (test_no_bfm_mlr_) call z%disableZonePolicyBfmMLR()
      
      ! zone interest modes index ptr
      read(unit_dump_bfm_) z%id_im_

      ! once zone is undumped, get its N. of points
      ! NOTE: needed in order to correctly index into bfm_undump !
      zNp = z%zoneTotNPts()

#ifdef __BSA_OMP
      if (.not. allocated(bfm_undump)) then
         allocate(bfm_undump(dimM_bisp_, zNp))
      else
         ! reallocate if more space needed
         if (zNp > size(bfm_undump, 2)) then
            deallocate(bfm_undump)
            allocate(bfm_undump(dimM_bisp_, zNp))
         endif
      endif
#endif

      ! read actual BFM dumped data
      ! NOTE: in second dimension, nj leading over ni
      !       laydown.
      read(unit_dump_bfm_) bfm_undump(:, 1 : zNp)
   end subroutine UndumpZone




   elemental pure subroutine disableZonePolicyBfmMLR(z)
      class(MZone_t), intent(inout) :: z

      z%policy_%interp_bfm_I_fct_ = 1
      z%policy_%interp_bfm_J_fct_ = 1
   end subroutine


end module
