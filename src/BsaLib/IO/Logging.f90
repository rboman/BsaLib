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
module Logging

   use BsaLib_IO, only: unit_dump_bfm_, unit_debug_, undebug_fname_
   use BsaLib_CONSTANTS
   implicit none
   private

! #ifdef _BSA_ALLOC_DEBUG
!    interface allocOKMsg
!       module procedure allocOKMsg_scalar_
!       module procedure allocOKMsg_array_
!    end interface
!    public :: allocOKMsg, deallocOKMsg
! #endif

   public :: allocKOMsg, deallocKOMsg

   type, public :: logger_t

      private
      integer(int32) :: iun_ = 0
      character(len=:), allocatable :: fileName_

   contains

      ! setting up

      procedure, public :: init
      procedure, public :: name
      procedure, public :: setName
      procedure, public :: unit
      procedure, public :: setUnit


      ! logging procedures 

      procedure, public :: logZonePremeshingTotTime
   
      ! procedure, public, pass :: LogWindSpectralAnalysisHeader, LogWindSPeedProfileType
      ! procedure, public, pass :: LogWindZoneData, LogTransFuncType, LogMaxValueType
      ! procedure, public, pass :: LogElementWindLoad, LOgDLMWindCoeffs, LogElementWindIncidenceAngles
      ! procedure, public, pass :: LogElemWindNodalVel
   end type logger_t




   interface

      module subroutine init(this, iun, fname)
         class(logger_t),   intent(inout) :: this
         integer(int32), intent(in)       :: iun
         character(len=*), intent(in), optional :: fname
      end subroutine



      module function name(this) result(nam)
         class(logger_t), intent(in)   :: this
         character(len=:), allocatable :: nam
      end function

      module subroutine setName(this, fname)
         class(logger_t), intent(inout) :: this
         character(len=*), intent(in)   :: fname
      end subroutine



      module function unit(this) result(iun)
         class(logger_t), intent(in) :: this
         integer :: iun
      end function

      module subroutine setUnit(this, iun)
         class(logger_t), intent(inout) :: this
         integer, intent(in), target    :: iun
      end subroutine




      
      module subroutine logZonePremeshingTotTime(this, zname, rtime, npts, print2console)
         class(logger_t), intent(in)   :: this
         character(len=*), intent(in)  :: zname
         real(real64), intent(in)      :: rtime
         integer, intent(in), optional :: npts
         logical, intent(in), optional :: print2console
      end subroutine





!=========================================================================================
!
!     ALLOCATION
!
!=========================================================================================

      module subroutine allocKOMsg(name_, istat, emsg)
         character(len = *), intent(in) :: name_, emsg
         integer, intent(in) :: istat
      end subroutine

      module subroutine deallocKOMsg(name_, istat, emsg)
         character(len = *), intent(in) :: name_, emsg
         integer, intent(in) :: istat
      end subroutine


   end interface


end module Logging