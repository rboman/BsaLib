module BsaCL

#include "_devtypes.h"
#include "_errtypes.h"

   use iso_c_binding, only: c_int, c_double
   implicit none
   public
   private :: evalFunc_C_to_F_wrapper_
   abstract interface
      function bsacl_EvalFunc_C_to_Fortran__(f_, nf_, itc_) result(res)
         integer(kind = 4), intent(in) :: nf_, itc_
         real(kind = 8), intent(in)    :: f_(nf_)
         real(kind = 8), allocatable, target :: res(:, :)
      end function
   end interface


   integer(kind = 4), parameter :: BSACL_DEVICE_TYPE_CPU              = BSACL_DEVICE_TYPE_CPU_
   integer(kind = 4), parameter :: BSACL_DEVICE_TYPE_GPU              = BSACL_DEVICE_TYPE_GPU_
   integer(kind = 4), parameter :: BSACL_DEVICE_TYPE_ACC              = BSACL_DEVICE_TYPE_ACC_
   integer(kind = 4), parameter :: BSACL_DEVICE_TYPE_DEF              = BSACL_DEVICE_TYPE_DEF_
   integer(kind = 4), parameter :: BSACL_VALUE_MISMATCH_ERROR         = BSACL_VALUE_MISMATCH_ERROR_
   integer(kind = 4), parameter :: BSACL_PLATFORM_FIND_ERROR          = BSACL_PLATFORM_FIND_ERROR_
   integer(kind = 4), parameter :: BSACL_DEF_PLATFORM_ACQUIRE_ERROR   = BSACL_DEF_PLATFORM_ACQUIRE_ERROR_
   integer(kind = 4), parameter :: BSACL_CQUEUES_CREATION_ERROR       = BSACL_CQUEUES_CREATION_ERROR_
   integer(kind = 4), parameter :: BSACL_INVALID_DEVICE_TYPE          = BSACL_INVALID_DEVICE_TYPE_
   integer(kind = 4), parameter :: BSACL_PROBLEM_DIMENSIONS_TOO_SMALL = BSACL_PROBLEM_DIMENSIONS_TOO_SMALL_


   interface
      module subroutine bsacl_Init(ierr)
         integer, intent(inout), target :: ierr
      end subroutine

      module subroutine bsacl_Run(ierr)
         integer, intent(inout), target :: ierr
      end subroutine

      module subroutine bsacl_AcquirePSDId(psdid)
         integer(kind = 4), intent(in) :: psdid
      end subroutine

      module subroutine bsacl_AcquireStructModMat(modmat, natf)
         real(kind = 8), intent(in), target :: modmat(:, :), natf(:)
      end subroutine

      module subroutine bsacl_AcquireLoadedNodesList(nodes_load)
         integer(kind = 4), intent(in), target :: nodes_load(:)
      end subroutine

      module subroutine bsacl_AcquireTotalNOfNodes(nn)
         integer(kind = 4), intent(in) :: nn
      end subroutine

      module subroutine bsacl_AcquireUsedModesList(modes)
         integer(kind = 4), intent(in), target :: modes(:)
      end subroutine

      module subroutine bsacl_AcquireWindCoeffs(wfc)
         real(kind = 8), intent(in), target :: wfc(:, :, :)
      end subroutine

      module subroutine bsacl_AcquireTurbComponentsList(tc)
         integer(kind = 4), intent(in), target :: tc(:)
      end subroutine

      module subroutine bsacl_AcquirePhiTimesCMat(phi_T_c)
         real(kind = 8), intent(in), target :: phi_T_c(:, :, :)
      end subroutine


      module subroutine bsacl_AcquireNodalCorrelation(nod_corr)
         real(kind = 8), intent(in), target :: nod_corr(:, :)
      end subroutine



      module subroutine bsacl_AcquireWindNodalVelocities(nod_vel)
         real(kind = 8), intent(in), target :: nod_vel(:)
      end subroutine

      module subroutine bsacl_AcquireWindNodalWindZones(nod_wz)
         integer(kind = 4), intent(in), target :: nod_wz(:)
      end subroutine

      module subroutine bsacl_AcquireWindTurbScales(wt_scl, nwz)
         real(kind = 8), intent(in), target :: wt_scl(:, :, :)
         integer(kind = 4), intent(in)      :: nwz
      end subroutine

      module subroutine bsacl_AcquireWindTurbStd(wt_std, nwz)
         real(kind = 8), intent(in), target :: wt_std(:, :)
         integer(kind = 4), intent(in)      :: nwz
      end subroutine



      module subroutine evalFunc_C_to_F_wrapper_(itc, nf, f, innl, res) bind(c)
         import c_int, c_double
         integer(c_int), value         :: itc
         integer(c_int), value         :: nf, innl
         real(c_double), intent(in), target :: f
         real(c_double), intent(inout) :: res(*)
      end subroutine

      module subroutine bsacl_AcquireEvaluationFunc(evalfct)
         procedure(bsacl_EvalFunc_C_to_Fortran__), pointer, intent(in) :: evalfct
      end subroutine

      module subroutine bsacl_AcquireEvalFuncByStrings(strings)
         character(len = *), target, intent(in) :: strings
      end subroutine

      module subroutine bsacl_AcquireEvalFuncByFile(filename)
         character(len = *), target, intent(in) :: filename
      end subroutine


      module subroutine bsacl_AcquireComputationFreqs(nfi, fi, nfj, fj)
         integer(kind = 4), intent(in)      :: nfi, nfj
         real(kind = 8), intent(in), target :: fi(..), fj(..)
      end subroutine


      module subroutine bsacl_AcquireBaseWindTurbPSD(S_uvw)
         real(kind = 8), intent(in), target :: S_uvw(:, :)
      end subroutine


      module subroutine bsacl_AcquireResultBFMVect(m3mf)
         real(kind = 8), intent(in), target :: m3mf(:)
      end subroutine


      module subroutine bsacl_SetDeviceType(itype)
         integer, intent(in) :: itype
      end subroutine


      module logical function bsacl_VerifyMaxAllocCondition(idim) result(bool)
         integer(kind = 8), intent(in) :: idim
      end function


      module subroutine bsacl_Abort(ierr)
         integer, intent(in) :: ierr
      end subroutine

      module subroutine bsacl_Finalise()
      end subroutine
   end interface

end module
