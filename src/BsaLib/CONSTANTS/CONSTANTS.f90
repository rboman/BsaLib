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
module BsaLib_CONSTANTS

   use, intrinsic :: iso_fortran_env
#ifndef IK
# define IK int32
#endif
#ifndef RK
# define RK real64
#endif
   implicit none
   public


   !**************************************************************************************
   !   BSA  integer/real types selected kinds
   !**************************************************************************************
   integer, parameter :: bsa_int_t  = IK
   integer, parameter :: bsa_real_t = RK


   !**************************************************************************************
   !   BSA  GENERICS
   !**************************************************************************************
   integer(int32), parameter :: BSA_SPATIAL_SYM_NONE = 1_int32
   integer(int32), parameter :: BSA_SPATIAL_SYM_HALF = 2_int32
   integer(int32), parameter :: BSA_SPATIAL_SYM_FOUR = 4_int32

   integer(int32), parameter :: BSA_PREMESH_MODE_BASE         = 0_int32
   integer(int32), parameter :: BSA_PREMESH_MODE_ZONE_REFINED = 1_int32
   
   integer(int32), parameter :: BSA_PREMESH_TYPE_DIAG_CREST_NO  = 0_int32
   integer(int32), parameter :: BSA_PREMESH_TYPE_DIAG_CREST_YES = 1_int32

   integer(int32), parameter :: BSA_VALIDATE_DELTAS_POLICY_NONE    = 0_int32
   integer(int32), parameter :: BSA_VALIDATE_DELTAS_POLICY_DEFAULT = 1_int32
   integer(int32), parameter :: BSA_VALIDATE_DELTAS_POLICY_LIGHT   = 2_int32
   integer(int32), parameter :: BSA_VALIDATE_DELTAS_POLICY_MEDIUM  = 3_int32
   integer(int32), parameter :: BSA_VALIDATE_DELTAS_POLICY_HIGH    = 4_int32
   integer(int32), parameter :: BSA_VALIDATE_DELTAS_POLICY_STRICT  = 5_int32

   !**************************************************************************************
   !   BSA  I/O  DEFAULTS
   !**************************************************************************************
   character(len = *), parameter :: BSA_OUT_DIRNAME_DEFAULT          = '.\bsaresults\'
   character(len = *), parameter :: BSA_OUT_FILENAME_PREFIX_DEFAULT_ = 'bsaout_def'
   character(len = *), parameter :: BSA_STRUCT_DATA_DUMPFILE         = 'dumpstruct'
   character(len = *), parameter :: BSA_SETTS_DATA_DUMPFILE          = 'dumpsetts'
   character(len = *), parameter :: BSA_WIND_DATA_DUMPFILE           = 'dumpwind'


   !**************************************************************************************
   !   LOGGING TAGs
   !**************************************************************************************
   character(len = *), parameter :: MSGCONT = '            '
   character(len = *), parameter :: INFOMSG = '  --[info]  '
   character(len = *), parameter :: NOTEMSG = '  --[note]  '
   character(len = *), parameter :: WARNMSG = '  --[warn]  '
   character(len = *), parameter :: ERRMSG  = '  --[error] '
   character(len = *), parameter :: DBGMSG  = '  --[debug] '


   !**************************************************************************
   !  I/O CONSTANTs
   !**************************************************************************
   character(len = *), parameter :: IO_ACCESS_DIRECT = 'DIRECT'
   character(len = *), parameter :: IO_ACCESS_SEQUEN = 'SEQUENTIAL'
   character(len = *), parameter :: IO_ACCESS_STREAM = 'STREAM'
   character(len = *), parameter :: IO_ACCESS_APPEND = 'APPEND'

   character(len = *), parameter :: IO_ACTION_WRITE     = 'WRITE'
   character(len = *), parameter :: IO_ACTION_READ      = 'READ'
   character(len = *), parameter :: IO_ACTION_READWRITE = 'READWRITE'

   character(len = *), parameter :: IO_ASYNC_YES = 'YES'
   character(len = *), parameter :: IO_ASYNC_NO  = 'NO'

   character(len = *), parameter :: IO_BUFFERED_YES = 'YES'
   character(len = *), parameter :: IO_BUFFERED_NO  = 'NO'

   character(len = *), parameter :: IO_FORM_FORMATTED   = 'FORMATTED'
   character(len = *), parameter :: IO_FORM_UNFORMATTED = 'UNFORMATTED'
   character(len = *), parameter :: IO_FORM_BINARY      = 'BINARY'

   character(len = *), parameter :: IO_POSITION_ASIS   = 'ASIS'
   character(len = *), parameter :: IO_POSITION_REWIND = 'REWIND'
   character(len = *), parameter :: IO_POSITION_APPEND = 'APPEND'

   character(len = *), parameter :: IO_STATUS_OLD     = 'OLD'
   character(len = *), parameter :: IO_STATUS_NEW     = 'NEW'
   character(len = *), parameter :: IO_STATUS_SCRATCH = 'SCRATCH'
   character(len = *), parameter :: IO_STATUS_REPLACE = 'REPLACE'
   character(len = *), parameter :: IO_STATUS_UNKNOWN = 'UNKNOWN'


   !**************************************************************************
   !  EXPORT CONSTANTs
   !**************************************************************************
   integer(int32), parameter :: BSA_EXPORT_FORMAT_FORMATTED   = 0_int32
   integer(int32), parameter :: BSA_EXPORT_FORMAT_UNFORMATTED = 1_int32
   integer(int32), parameter :: BSA_EXPORT_MODE_APPEND  = 0_int32
   integer(int32), parameter :: BSA_EXPORT_MODE_REPLACE = 1_int32

   integer(int32), parameter :: BSA_EXPORT_BRM_MODE_NONE = 0_int32
   integer(int32), parameter :: BSA_EXPORT_BRM_MODE_BASE = 1_int32
   integer(int32), parameter :: BSA_EXPORT_BRM_MODE_USR  = 9_int32


   character(len = *), parameter :: BSA_EXPORT_M2MF_CLS_FNAME   = "m2mf_cls"
   character(len = *), parameter :: BSA_EXPORT_M2MR_CLS_FNAME   = "m2mr_cls"
   character(len = *), parameter :: BSA_EXPORT_M2O2MR_CLS_FNAME = "m2o2mr_cls"
   character(len = *), parameter :: BSA_EXPORT_M3MF_CLS_FNAME   = "m3mf_cls"
   character(len = *), parameter :: BSA_EXPORT_M3MR_CLS_FNAME   = "m3mr_cls"
   character(len = *), parameter :: BSA_EXPORT_M2MF_MSH_FNAME   = "m2mf_msh"
   character(len = *), parameter :: BSA_EXPORT_M2MR_MSH_FNAME   = "m2mr_msh"
   character(len = *), parameter :: BSA_EXPORT_M2O2MR_MSH_FNAME = "m2o2mr_msh"
   character(len = *), parameter :: BSA_EXPORT_M3MF_MSH_FNAME   = "m3mf_msh"
   character(len = *), parameter :: BSA_EXPORT_M3MR_MSH_FNAME   = "m3mr_msh"


   abstract interface
      subroutine exportBRMinterf_vect_(f1, f2, brm, pdata)
         import :: bsa_real_t
         real(bsa_real_t), intent(in)  :: f1(:), f2(:), brm(:, :)
         class(*), pointer, intent(in) :: pdata
      end subroutine
   end interface




   !**************************************************************************************
   !   NUMERICs
   !**************************************************************************************
   
   !> TO AVOID CRASHING BECAUSE OF MACHINE FLOATING PRECISION ERRORS
   real(real64), parameter :: MACHINE_PRECISION = 1e-12_real64

   real(real64), parameter :: CST_PIGREC = 4.0_real64 * atan(1.0_real64)

   real(real64), parameter :: CST_PIt2   = CST_PIGREC * 2.0_real64
   real(real64), parameter :: CST_PIt4   = CST_PIGREC * 4.0_real64

   real(real64), parameter :: CST_PId2   = CST_PIGREC / 2.0_real64
   real(real64), parameter :: CST_PId4   = CST_PIGREC / 4.0_real64

   
   real(real64), parameter :: CST_2d3    = 2.0_real64 / 3.0_real64
   real(real64), parameter :: CST_3d2    = 3.0_real64 / 2.0_real64
   real(real64), parameter :: CST_PIt3d2 = CST_PIGREC * CST_3d2

   logical, protected :: header_called_ = .false.

contains


   subroutine bsa_printBSAHeader()
      print *
      print *
      print *, '              ____________________________________________ '
      print *, '             |     _____         ____                     |'
      print *, '             |      /   \       /              /\         |'
      print *, '             |     /____/       \___          /  \        |'
      print *, '             |    /    \            \        /____\       |'
      print *, '             |   /_____/  .    _____/  .   _/      \_  .  |'
      print *, '             |____________________________________________|'
      print *
      print *

      header_called_ = .true.
   end subroutine bsa_printBSAHeader


end module