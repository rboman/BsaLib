submodule(BsaLib_IO) BsaLib_IO_Impl

   use BsaLib_Data, only: bsa_Abort
   implicit none

contains

   module subroutine allocKOMsg(name_, istat, emsg)
      character(len = *), intent(in) :: name_, emsg
      integer, intent(in) :: istat

      write(unit_debug_, fmt='(3a)') &
         '[ERROR] variable  "', name_, '"  could not be allocated.'
      write(unit_debug_, fmt='(15x, a, i0, 2a)') &
         'Exit code  ', istat, '. Error message:  ', emsg(1 : len_trim(emsg))

      call bsa_Abort()
   end subroutine



   module subroutine deallocKOMsg(name_, istat, emsg)
      character(len = *), intent(in) :: name_, emsg
      integer, intent(in) :: istat

      write(unit_debug_, fmt='(3a)') &
         '[ERROR] variable  "', name_, '"  could not be de-allocated.'
      write(unit_debug_, fmt='(15x, a, i0, 2a)') &
         'Exit code  ', istat, '. Error message:  ', emsg(1 : len_trim(emsg))

      call bsa_Abort()
   end subroutine
end submodule
