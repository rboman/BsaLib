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
module BsaLib_Timing

   use BsaLib_CONSTANTS, only: real64
   implicit none
   private


   type, public :: timer_t
      private
      real(real64) :: t_init_     = 0._real64
      real(real64) :: t_last_     = 0._real64
      real(real64) :: t_tot_      = 0._real64
      real(real64) :: t_tot_prev_ = 0._real64
   contains
      procedure, public, pass(this) :: init  => InitTimer
      procedure, public, pass(this) :: time  => ClockTimer
      procedure, public, pass(this) :: total => GetTimerTotal
      procedure, public, pass(this) :: reset => ResetTimer
   end type timer_t




contains



   subroutine InitTimer(this)
      class(timer_t) :: this

      ! write(unit_debug_, *), ' @BsaLib_Timing::InitTimer() : Init timer...'
      
      call cpu_time(this%t_init_)
      this%t_last_ = this%t_init_

      ! write(unit_debug_, *), ' @BsaLib_Timing::InitTimer() : Init timer -- ok.'
   end subroutine InitTimer




   function ClockTimer(this) result(dt)
      class(timer_t) :: this
      real(real64) :: dt_tmp
      real(real64) :: dt

! #ifdef __BSA_DEBUG
!       write(unit_debug_, *), ' @BsaLib_Timing::ClockTimer() : save partial time...'
! #endif

      call cpu_time(dt_tmp)
      dt = dt_tmp - this%t_last_

      ! update total
      this%t_tot_ = this%t_tot_ + dt

      ! update last cpu_time call
      this%t_last_ = dt_tmp

! #ifdef __BSA_DEBUG
!       write(unit_debug_, *), ' @BsaLib_Timing::ClockTimer() : save partial time -- ok.'
! #endif
   end function ClockTimer





   pure elemental function GetTimerTotal(this) result(tot)
      class(timer_t), intent(in) :: this
      real(real64) :: tot

! #ifdef __BSA_DEBUG
!       write(unit_debug_, *), ' @BsaLib_Timing::GetTimerTotal() : getting timer total time...'
! #endif

      tot = this%t_tot_

! #ifdef __BSA_DEBUG
!       write(unit_debug_, *), ' @BsaLib_Timing::GetTimerTotal() : getting timer total time -- ok.'
! #endif
   end function GetTimerTotal




   subroutine ResetTimer(this)
      class(timer_t) :: this

! #ifdef __BSA_DEBUG
!       write(unit_debug_, *), ' @BsaLib_Timing::ResetTimer() : Reset timer...'
! #endif

      this%t_init_     = 0._real64
      this%t_last_     = 0._real64
      this%t_tot_      = 0._real64
      this%t_tot_prev_ = 0._real64

! #ifdef __BSA_DEBUG
!       write(unit_debug_, *), ' @BsaLib_Timing::ResetTimer() : Reset timer -- ok.'
! #endif
   end subroutine ResetTimer



end module