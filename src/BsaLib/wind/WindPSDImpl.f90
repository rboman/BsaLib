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
submodule(BsaLib_WindData) BsaLib_WindPSDImpl

#include "../precisions"

   use BsaLib_IO, only: INFOMSG, WARNMSG, ERRMSG, MSGCONT, DBGMSG
   implicit none

   type :: arr_proc_pointer_t
      procedure(PSDfunc), pointer, nopass :: ptr => null()
   end type arr_proc_pointer_t

   ! TODO: might be statically initialised (parameter)
   type(arr_proc_pointer_t), dimension(5) :: psd_funcs

contains


   module subroutine SetPSDType(this, ipsd)
      class(WindData_t), intent(inout) :: this
      integer(kind = 4), value :: ipsd
      integer :: istat
      character(len = 256) :: emsg

      ! if (ipsd < 1 .or. ipsd > 5) call bsa_Abort('Invalid "ipsd" value.')
      ! if (ipsd < 1 .or. ipsd > 5) ipsd = 1
      if (ipsd < 1 .or. ipsd > 5) error stop ERRMSG//'Invalid "ipsd" value.'

      ! TODO: this might be removed
      this%i_psd_type_ = ipsd

      allocate(this%psd_, stat=istat, errmsg=emsg)
      if (istat == 0) then
#ifdef __BSA_ALLOC_DEBUG
         call allocOKMsg('this % psd_', loc(this%psd_), sizeof(this%psd_))
#endif
      else
         call allocKOMsg('this % psd_', istat, emsg)
      endif

      ! setting function pointer array 
      psd_funcs(1)%ptr => vonKarmanPSD_
      psd_funcs(2)%ptr => kaimalPSD_
      psd_funcs(3)%ptr => davenportPSD_Greisch_
      psd_funcs(5)%ptr => davenportPSD_Uliege_

      call this%psd_%SetPSDFunction(psd_funcs(ipsd)%ptr)

#ifdef __BSA_DEBUG
      print *, INFOMSG, '@WindImpl::SetPSDType() : PSD type set to ', this%i_psd_type_
#endif
   end subroutine SetPSDType





   module function getFullNodalPSD(this, innl, nodesl, PSDvec, f, idir) result(PSDmat)
      use BsaLib_Utility, only: util_getCorrVectIndex
      use BsaLib_Data, only: struct_data
      class(WindData_t), intent(in) :: this
      integer(kind = 4), intent(in) :: innl, idir
      integer(kind = 4), intent(in) :: nodesl(innl)
      real(RDP), intent(in) :: PSDvec(innl)
      real(RDP), intent(in) :: f
      real(RDP) :: PSDmat(innl, innl)
      real(RDP) :: absf
      integer(kind = 4) :: i, j, ni, nj, id

      absf = abs(f)

      ! do concurrent (i = 1 : innl) shared(PSDvec)
      !    PSDmat(:, i) = sqrt(PSDvec * PSDvec(i))
      !    ni = nodesl(i)
      !    do concurrent (j = 1 : innl) shared(this, i)
      !       nj           = nodesl(j)
      !       id           = util_getCorrVectIndex(nj, ni, struct_data%nn_)
      !       PSDmat(j, i) = PSDmat(j, i) * &
      !          this%nod_corr_(id, idir)**(absf)
      !    enddo
      ! enddo
      do i = 1, innl
         PSDmat(:, i) = sqrt(PSDvec * PSDvec(i))
         ni = nodesl(i)
         do j = 1, innl
            nj           = nodesl(j)
            id           = util_getCorrVectIndex(nj, ni, struct_data%nn_)
            PSDmat(j, i) = PSDmat(j, i) * &
               this%nod_corr_(id, idir)**(absf)
         enddo
      enddo
   end function







   module subroutine SetPSDFunction(this, func)
      class(psd_t) :: this
      procedure(PSDfunc), intent(in), pointer :: func

! #ifdef __BSA_DEBUG
!       write(unit_debug_, '(1x, a, a)') INFOMSG, '@WindPSDImpl::SetPSDFunction() : setting PSD function pointer...'
! #endif

      this%psd_fct_ptr => func

#ifdef __BSA_DEBUG
      write(unit_debug_, '(1x, a, a)') INFOMSG, '@WindPSDImpl::SetPSDFunction() : setting PSD function pointer -- ok.'
#endif
   end subroutine SetPSDFunction



   module function evalPSD_(this, nf, f, innl, nnl, idir, itc) result(PSD)
      use BsaLib_Data, only: settings
      class(WindData_t), intent(in) :: this
      integer(kind = 4), intent(in)           :: nf, innl, idir, itc
      integer(kind = 4), intent(in) :: nnl(innl)
      real(RDP), intent(in)  :: f(nf)
      real(RDP) :: PSD(nf, innl)

      if (idir /= 1) then
         print '(/ 1x, a, a, i1, a)', &
            WARNMSG, 'IDIR=  ', idir, ', when usually SHOULD be 1 (X wind direction).'
         print '(1x, a, a/, a, a)', &
            MSGCONT, '(It is uncommon to compute PSDs of wind turbulence assuming', &
            MSGCONT, ' that vortices do not move along X (idir=1), direction of mean wind)'
      endif

      ! invoking internal function pointer
      PSD = this%psd_%psd_fct_ptr(this, nf, f, innl, nnl, idir, itc)

      ! NOTE: PSD scaling (convention based)
      if (settings%i_def_scaling_ == 1) PSD = PSD / CST_PIt4
   end function




   function vonKarmanPSD_(wd, nf, freqs, innl, nnl, idir, itc) result(PSD)
      class(WindData_t), intent(in) :: wd
      integer(kind = 4), intent(in) :: nf       ! n. frequencies
      integer(kind = 4), intent(in) :: innl     ! n. actual nodes loaded
      integer(kind = 4), intent(in) :: idir     ! wind direction
      integer(kind = 4), intent(in) :: itc      ! 
      real(RDP), intent(in)  :: freqs(nf)       ! frequencies
      integer(kind = 4), intent(in) :: nnl(innl)   ! list of actual loaded nodes
      real(RDP) :: PSD(nf, innl)

      ! tmp
      real(RDP), dimension(1, innl) :: L
      real(RDP), allocatable :: rtmp1(:, :)
      integer :: i


! #ifdef __BSA_DEBUG
!       write(unit_debug_, '(1x, a, a)') INFOMSG, '@WindPSDImpl::vonKarmanPSD_() : computing PSD...'
! #endif

      L(1, :) = wd%turb_scales_wz_(itc, idir, wd%wz_node_(1 : innl))

      ! NOTE: how was programmed in FINELG
      ! BUG: check.
      if (idir == 1) then

         ! L/U
         rtmp1 = L / reshape(wd%u_node_(1 : innl), [1, innl])
         PSD   = matmul(reshape(abs(freqs), [nf, 1]), rtmp1)
         PSD   = PSD * PSD ! square
         rtmp1 = rtmp1 * reshape(wd%sigmaUVW_wz_(itc, wd%wz_node_(1 : innl))**2, [1, innl])
         PSD   = (1 + 70.7_RDP * PSD)**(5._RDP/6._RDP)
         

         ! BUG: LOCAL(i) is superflous (should be error..)
         do concurrent (i = 1 : innl) local(i)
            PSD(:, i) = 4._RDP * rtmp1(1, i) / PSD(:, i)
         enddo

      else ! WARNING: should not pass from here

         block
            real(RDP) :: dnlsu(nf, innl), rtmp2(nf, innl), rtmp3(nf, innl)

            dnlsu = 2._RDP * &
               matmul(reshape(freqs, [nf, 1]), &
               reshape(wd%turb_scales_wz_(1, idir, wd%wz_node_(nnl)), [1, innl]) / &
               reshape(wd%u_node_(1 : innl), [1, innl]))

            dnlsu = dnlsu*dnlsu

            rtmp1 = 1._RDP + 70.7_RDP * dnlsu
            rtmp2 = rtmp1 ** (11._RDP / 6._RDP)
            rtmp2 = rtmp2 * reshape(wd%u_node_(1 : innl), [1, innl])
            rtmp3 = reshape(wd%turb_scales_wz_(itc, idir, wd%wz_node_(nnl)), [1, innl]) * &
               (1._RDP + 188.4_RDP * dnlsu) / &
               rtmp2 * reshape(wd%sigmaUVW_wz_(itc, wd%wz_node_(nnl))**2, [1, innl])

            rtmp1 = rtmp3 + rtmp3
            PSD   = rtmp1 + rtmp1
         end block
      endif

! #ifdef __BSA_DEBUG
!       write(unit_debug_, '(1x, a, a)') INFOMSG, '@WindPSDImpl::vonKarmanPSD_() : computing PSD -- ok.'
! #endif
   end function vonKarmanPSD_




   function kaimalPSD_(wd, nf, freqs, innl, nnl, idir, itc) result(PSD)
      class(WindData_t), intent(in) :: wd
      integer(kind = 4), intent(in) :: nf                    ! n. frequencies
      integer(kind = 4), intent(in) :: innl                  ! n. actual nodes loaded
      integer(kind = 4), intent(in) :: idir                  ! wind direction
      integer(kind = 4), intent(in) :: itc                   ! 
      integer(kind = 4), intent(in) :: nnl(innl)     ! list of actual loaded nodes
      real(RDP), intent(in)  :: freqs(nf)   ! frequencies
      real(RDP) :: PSD(nf, innl)

#ifdef __BSA_DEBUG
      write(unit_debug_, '(1x, a, a)') INFOMSG, '@WindPSDImpl::kaimalPSD_() : computing PSD.. [NOT YET IMPLEMENTED]'
#endif

      PSD = 0._RDP

! #ifdef __BSA_DEBUG
!       write(unit_debug_, '(1x, a, a)') &
!             INFOMSG, '@WindPSDImpl::kaimalPSD_() : computing PSD -- ok.'
! #endif
   end function kaimalPSD_





   function davenportPSD_Greisch_(wd, nf, freqs, innl, nnl, idir, itc) result(PSD)
      class(WindData_t), intent(in) :: wd
      integer(kind = 4), intent(in) :: nf          ! n. frequencies
      integer(kind = 4), intent(in) :: innl        ! n. actual nodes loaded
      integer(kind = 4), intent(in) :: idir        ! wind direction
      integer(kind = 4), intent(in) :: itc         ! 
      integer(kind = 4), intent(in) :: nnl(innl)     ! list of actual loaded nodes
      real(RDP), intent(in) :: freqs(nf)   ! frequencies
      real(RDP), parameter :: cst1 = 1200._RDP
      integer   :: i, n
      real(RDP) :: PSD(nf, innl)

! #ifdef __BSA_DEBUG
!       write(unit_debug_, '(1x, a, a)') INFOMSG, '@WindPSDImpl::davenportPSD_Greisch_() : computing PSD...'
! #endif

      do i = 1, innl

         n = nnl(i)

         PSD(:, i) = wd%sigmaUVW_wz_(itc, wd%wz_node_(n)) * wd%sigmaUVW_wz_(itc, wd%wz_node_(n)) * &
            0.65_RDP * cst1 / wd%u_mean_ref_wz_(wd%wz_node_(n)) / &
            (1 + (cst1 * freqs / wd%u_mean_ref_wz_(wd%wz_node_(n)))**2._RDP)**(5._RDP / 6._RDP)
      enddo

! #ifdef __BSA_DEBUG
!       write(unit_debug_, '(1x, a, a)') INFOMSG, '@WindPSDImpl::davenportPSD_Greisch_() : computing PSD -- ok.'
! #endif
   end function davenportPSD_Greisch_





   function davenportPSD_Uliege_(wd, nf, freqs, innl, nnl, idir, itc) result(PSD)
      class(WindData_t), intent(in) :: wd
      integer(kind = 4), intent(in) :: nf                    ! n. frequencies
      integer(kind = 4), intent(in) :: innl                  ! n. actual nodes loaded
      integer(kind = 4), intent(in) :: idir                  ! wind direction
      integer(kind = 4), intent(in) :: itc                   ! 
      integer(kind = 4), intent(in) :: nnl(innl)     ! list of actual loaded nodes
      real(RDP), intent(in)  :: freqs(nf)   ! frequencies
      real(RDP) :: PSD(nf, innl)
      real(RDP) :: cstL_U(1, innl), cstFL_U(nf, innl)
      integer   :: i, n


      cstL_U(1, :) = wd%turb_scales_wz_(itc, idir, wd%wz_node_(nnl)) / wd%u_node_(nnl)
      cstFL_U      = matmul(reshape(abs(freqs), [nf, 1]), cstL_U)


      ! do i = 1, innl
      !    PSD(:, i) = &
      !       (CST_2d3 * cstFL_U(:, i) * cstL_U(1, i) * wd%sigmaUVW_wz_(itc, wd%wz_node_(i))**2) / &
      !       ((1.d0 + (cstFL_U(:, i)**2))**(4.d0 / 3.d0))
      ! enddo
      do concurrent (i = 1 : innl) shared(cstFL_U, cstL_U, itc, wd) local(n)

         n = nnl(i)

         PSD(:, i) = &
            (CST_2d3 * cstFL_U(:, i) * cstL_U(1, i) * wd%sigmaUVW_wz_(itc, wd%wz_node_(n))**2) / &
            ((1.d0 + (cstFL_U(:, i)**2))**(4.d0 / 3.d0))
      enddo

! #ifdef __BSA_DEBUG
!       write(unit_debug_, '(1x, a, a)') INFOMSG, '@WindPSDImpl::davenportPSD_Uliege_() : computing PSD -- ok.'
! #endif
   end function davenportPSD_Uliege_




end submodule