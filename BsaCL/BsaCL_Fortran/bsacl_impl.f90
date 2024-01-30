submodule(BsaCL) BsaCL_impl

#include "_msgtypes.h"

   use, intrinsic :: iso_c_binding
   use BsaCL_c
   implicit none

   procedure(bsacl_EvalFunc_C_to_Fortran__), pointer :: evalfct_ptr_ => null()

contains

   module subroutine bsacl_Init(ierr)
      integer, intent(inout), target :: ierr
      call bsaclInit__(c_loc(ierr))
   end subroutine

   module subroutine bsacl_Run(ierr)
      integer, intent(inout), target :: ierr
      call bsaclRun__(c_loc(ierr))
   end subroutine


   module function bsacl_SetKernelID(kid) result(ierr)
      integer, value, intent(in) :: kid
      integer :: ierr

      ierr = int(bsaclSetKernelID__(int(kid, c_int)))
   end function



   module subroutine bsacl_AcquirePSDId(psdid)
      integer(kind = 4), intent(in) :: psdid
      call bsaclAcquirePSDId__(int(psdid, c_int))
   end subroutine


   module subroutine bsacl_AcquireStructModMat(modmat, natf)
      real(kind = 8), intent(in), target :: modmat(:, :), natf(:)
      integer(c_int) :: ndofs_c_, nmodes_c_
      
      ndofs_c_  = size(modmat, dim=1, kind=c_int)
      nmodes_c_ = size(modmat, dim=2, kind=c_int)
      if (size(natf, kind=c_int) /= nmodes_c_) then
         print '(// 1x, 2a)', INFO_MSG, "N. of modes does not match between modal matrix and nat freqs vector."
         call bsaclAbort__(-200)
      endif
      call bsaclAcquireStructModMat__(c_loc(modmat), c_loc(natf), ndofs_c_, nmodes_c_)
   end subroutine

   module subroutine bsacl_AcquireLoadedNodesList(nodes_load)
      integer(kind = 4), intent(in), target :: nodes_load(:)
      integer(c_int) :: nnodes_l_
      nnodes_l_ = size(nodes_load, kind=c_int)
      call bsaclAcquireLoadedNodesList__(c_loc(nodes_load), nnodes_l_)
   end subroutine

   module subroutine bsacl_AcquireTotalNOfNodes(nn)
      integer(kind = 4), intent(in) :: nn
      call bsaclAcquireTotalNOfNodes__(int(nn, kind=c_int))
   end subroutine

   module subroutine bsacl_AcquireUsedModesList(modes)
      integer(kind = 4), intent(in), target :: modes(:)
      integer(c_int) :: nmodes_eff_
      nmodes_eff_ = size(modes, kind = c_int)
      call bsaclAcquireUsedModesList__(c_loc(modes), nmodes_eff_)
   end subroutine


   module subroutine bsacl_AcquireWindCoeffs(wfc)
      real(kind = 8), intent(in), target :: wfc(:, :, :)
      integer(c_int) :: nlibs_, ndegw_, nnodes_l_

      nlibs_    = size(wfc, 1, c_int)
      ndegw_    = size(wfc, 2, c_int)
      nnodes_l_ = size(wfc, 3, c_int)
      call bsaclAcquireWindCoeffs__(c_loc(wfc), nnodes_l_, nlibs_, ndegw_)
   end subroutine


   module subroutine bsacl_AcquireTurbComponentsList(tc)
      integer(kind = 4), intent(in), target :: tc(:)
      integer(c_int) :: ntc_
      ntc_ = size(tc, kind=c_int)
      call bsaclAcquireTurbComponentsList__(c_loc(tc), ntc_)
   end subroutine


   module subroutine bsacl_AcquirePhiTimesCMat(phi_T_c)
      real(kind = 8), intent(in), target :: phi_T_c(:, :, :)
      integer(c_int) :: nm_eff_, nnodes_l_, ndegw_

      nm_eff_   = size(phi_T_c, 1, c_int)
      nnodes_l_ = size(phi_T_c, 2, c_int)
      ndegw_    = size(phi_T_c, 3, c_int)
      call bsaclAcquirePhiTimesCMat__(c_loc(phi_T_c), nm_eff_, nnodes_l_, ndegw_)
! #ifdef BSA_DEBUG
!       block
!          integer :: m, n, d
!          integer :: id = 0
!          do d = 1, ndegw_
!             do n = 1, nnodes_l_
!                write(4532, '(2x, i5, *(2x, f15.5))') id, phi_T_c(:, n, d)
!                id = id + nm_eff_
!             enddo
!             write(4532, *) ''
!          enddo
!       end block
! #endif
   end subroutine


   module subroutine bsacl_AcquireNodalCorrelation(nod_corr)
      real(kind = 8), intent(in), target :: nod_corr(:, :)
      integer(c_int) :: nnod_corr
      nnod_corr = size(nod_corr, dim=1)
      call bsaclAcquireNodalCorrelation__(c_loc(nod_corr), nnod_corr)
! #ifdef BSA_DEBUG
!       print '(/1x, 2a, i/)', INFO_MSG, 'n. nod_corr  = ', nnod_corr
! #endif
   end subroutine



   module subroutine bsacl_AcquireWindNodalVelocities(nod_vel)
      real(kind = 8), intent(in), target :: nod_vel(:)

      call bsaclAcquireWindNodalVelocities__(c_loc(nod_vel))
   end subroutine
   
   module subroutine bsacl_AcquireWindNodalWindZones(nod_wz)
      integer(kind = 4), intent(in), target :: nod_wz(:)

      call bsaclAcquireWindNodalWindZones__(c_loc(nod_wz))
   end subroutine
   
   module subroutine bsacl_AcquireWindTurbScales(wt_scl, nwz)
      real(kind = 8), intent(in), target :: wt_scl(:, :, :)
      integer(kind = 4), intent(in)      :: nwz

      call bsaclAcquireWindTurbScales__(c_loc(wt_scl), int(nwz, c_int))
   end subroutine
   
   module subroutine bsacl_AcquireWindTurbStd(wt_std, nwz)
      real(kind = 8), intent(in), target :: wt_std(:, :)
      integer(kind = 4), intent(in)      :: nwz

      call bsaclAcquireWindTurbStd__(c_loc(wt_std), int(nwz, c_int))
   end subroutine





   ! This is what gets actually called in the C core
   module subroutine evalFunc_C_to_F_wrapper_(itc, nf, f, innl, res) bind(c)
      integer(c_int), value         :: itc
      integer(c_int), value         :: nf, innl
      real(c_double), intent(in), target :: f
      real(c_double), intent(inout) :: res(*)

      real(kind = 8), allocatable   :: res_(:, :)

      ! locals for Fortran procedure
      integer(kind = 4) :: nf_, itc_, i_
      real(kind = 8), pointer :: f_(:)

      nf_  = int(nf,  kind=4)
      itc_ = int(itc, kind=4)
      call c_f_pointer(c_loc(f), f_, [nf_])

      ! Then, calling internal Fortran function pointer to actual eval function
      res_  = evalfct_ptr_(f_, nf_, itc_)
      do i_ = 1, innl
         res((i_-1)*nf_ + 1 : (i_*nf_)) = real(res_(:, i_), kind=c_double)
      enddo
! #ifdef BSA_DEBUG
!       write(9542, *) "  Evaluation result is (Fortran) :"
!       do itc_ = 1, nf_
!          write(9542, fmt='(2x, i0, 2x, f10.4, *(2x, f12.5))', advance='no') itc, f_(itc_), res(1 : innl)
!       enddo
!       write(9542, *)
! #endif
   end subroutine


   module subroutine bsacl_AcquireEvaluationFunc(evalfct)
      procedure(bsacl_EvalFunc_C_to_Fortran__), pointer, intent(in) :: evalfct

      evalfct_ptr_ => evalfct
      call bsaclAcquireEvaluationFunc__(c_funloc(evalFunc_C_to_F_wrapper_))
   end subroutine


   module subroutine bsacl_AcquireEvalFuncByStrings(strings)
      character(len = *), target, intent(in) :: strings

      call bsaclAcquireEvalFuncByStrings__(c_loc(strings))
   end subroutine


   module subroutine bsacl_AcquireEvalFuncByFile(filename)
      character(len = *), target, intent(in) :: filename

      call bsaclAcquireEvalFuncByFile__(c_loc(filename))
   end subroutine




   module subroutine bsacl_AcquireComputationFreqs(nfi, fi, nfj, fj)
      integer(kind = 4), intent(in)      :: nfi, nfj
      real(kind = 8), intent(in), target :: fi(..), fj(..)
      integer(c_int) :: nfi_, nfj_

      nfi_ = int(nfi, kind=c_int)
      nfj_ = int(nfj, kind=c_int)
      call bsaclAcquireComputationFreqs__(nfi_, c_loc(fi), nfj_, c_loc(fj))
   end subroutine


   module subroutine bsacl_AcquireBaseWindTurbPSD(S_uvw)
      real(kind = 8), intent(in), target :: S_uvw(:, :)
      call bsaclAcquireBaseWindTurbPSD__(c_loc(S_uvw))
   end subroutine



   module subroutine bsacl_AcquireResultBFMVect(m3mf)
      real(kind = 8), intent(in), target :: m3mf(:)
      integer(c_int) :: idim_

      idim_ = size(m3mf, kind=c_int)
      call bsaclAcquireResultBFMVect__(c_loc(m3mf), idim_)
   end subroutine



   module subroutine bsacl_SetDeviceType(itype)
      integer, intent(in) :: itype
      call bsaclSetDeviceType__(int(itype, kind=c_int))
   end subroutine



   module logical function bsacl_VerifyMaxAllocCondition(idim) result(bool)
      integer(kind = 8), intent(in) :: idim
      integer(c_int64_t) :: idim_
      integer(c_int), target :: i_can_ = -1_c_int
      
      idim_ = int(idim, kind=c_int64_t)
      call bsaclVerifyMaxAllocCondition__(idim_, c_loc(i_can_))
      if (i_can_ == 0_c_int) then
         bool = .false.
      else
         bool = .true.
      endif
   end function




   module subroutine bsacl_Abort(ierr)
      integer, intent(in) :: ierr
      call bsaclAbort__(int(ierr, kind=c_int))
   end subroutine


   module subroutine bsacl_Finalise()
      call bsaclFinalise__()
   end subroutine
end submodule
