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
submodule(BsaLib_Data) BsaLib_DataImpl

   use BsaLib_CONSTANTS, only: INFOMSG, WARNMSG, ERRMSG, MSGCONT, DBGMSG
   use BsaLib_IO,        only: unit_debug_, unit_dump_bfm_, unit_dump_brm_
   implicit none (type, external)

contains


   module function evaluatePSD(f, nf, itc) result(PSD)
      integer(bsa_int_t), intent(in) :: nf, itc
      real(bsa_real_t), intent(in)   :: f(nf)
      real(bsa_real_t), allocatable, target :: PSD(:, :)

      PSD = wd%evalPSD(nf, f, struct_data%nn_load_, struct_data%n_load_, 1, itc)
   end function




   module subroutine cleanBSAData_()
      integer(int32) :: istat
      character(len = 256) :: emsg

#ifdef BSA_DEBUG
      write(unit_debug_, *) INFOMSG//'@BsaLibData::cleanBSAData_() : cleaning...'
#endif

      if (allocated(wd))          call wd%clean()
      if (allocated(struct_data)) call struct_data%clean()

      ! if (associated(m2mf_cls_ptr_)) nullify(m2mf_cls_ptr_)
      ! if (associated(m2mr_cls_ptr_)) nullify(m2mr_cls_ptr_)
      ! if (associated(m3mf_cls_ptr_)) nullify(m3mf_cls_ptr_)
      ! if (associated(m3mr_cls_ptr_)) nullify(m3mr_cls_ptr_)

      if (associated(m3mf_msh_ptr_)) nullify(m3mf_msh_ptr_)
      if (associated(m3mr_msh_ptr_)) nullify(m3mr_msh_ptr_)


      if (allocated(PHItimesC_local_)) then
         deallocate(PHItimesC_local_, stat=istat, errmsg=emsg)
         if (istat /= 0) call deallocKOMsg('PHItimesC_local_', istat, emsg)
      endif


      if (allocated(peak_exts_)) then
         deallocate(peak_exts_, stat=istat, errmsg=emsg)
         if (istat /= 0) call deallocKOMsg('peak_exts_', istat, emsg)
      endif


      if (allocated(settings)) then
         deallocate(settings, stat=istat, errmsg=emsg)
         if (istat /= 0) call deallocKOMsg('settings', istat, emsg)
      endif

      if (allocated(wd)) then
         deallocate(wd, stat=istat, errmsg=emsg)
         if (istat /= 0) call deallocKOMsg('wd', istat, emsg)
      endif

      if (allocated(struct_data)) then
         deallocate(struct_data, stat=istat, errmsg=emsg)
         if (istat /= 0) call deallocKOMsg('struct_data', istat, emsg)
      endif

      if (allocated(timer)) then
         deallocate(timer, stat=istat, errmsg=emsg)
         if (istat /= 0) call deallocKOMsg('timer', istat, emsg)
      endif

      if (allocated(logger_debug)) then
         deallocate(logger_debug, stat=istat, errmsg=emsg)
         if (istat /= 0) call deallocKOMsg('logger_debug', istat, emsg)
      endif

      is_data_cleaned_ = .true.

#ifdef BSA_DEBUG
      write(unit_debug_, *) INFOMSG//'@BsaLibData::cleanBSAData_() : cleaning -- ok.'
#endif

      call closeBSAUnits_()
   end subroutine



   subroutine closeBSAUnits_()
      logical :: isopn

      ! NOTE: keep conditions since they might be provided from 
      !       host program, so they would not want me to close them.
      !       They'll manage it ;)
      if (close_deb_unit_) then
         inquire(unit = unit_debug_, opened = isopn)
         if (isopn) close(unit_debug_)
      endif

      inquire(unit = unit_dump_bfm_, opened = isopn)
      if (isopn) close(unit_dump_bfm_)

      inquire(unit = unit_dump_brm_, opened = isopn)
      if (isopn) close(unit_dump_brm_)

#ifdef _BSA_EXPORT_POD_TRUNC_INFO
      inquire(unit = iun_POD_trunc_, opened = isopn)
      if (isopn) close(iun_POD_trunc_)
#endif
   end subroutine



   module subroutine bsa_Abort(emsg)
      character(len = *), intent(in), optional :: emsg
      ! external :: abort

      if (present(emsg)) print '(/ 1x, a, a/)', ERRMSG, emsg
      
      call cleanBSAData_() ! free memory before halting
      error stop
   end subroutine

end submodule
