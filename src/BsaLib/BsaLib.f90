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
module BsaLib

#include "./precisions"

   use BsaLib_CONSTANTS
   use BsaLib_IO
   implicit none
   public

   private :: mainClassic_, mainMesher_

   ! BUG: these might be moved, and called via an internal function pointer.
   private :: bsa_exportSkewness_compute_, bsa_exportSkewness_nocompute_
   private :: bsa_exportBR_nocompute_



   interface bsa_exportBRdecomp
      module procedure bsa_exportBR_nocompute_
   end interface


   interface bsa_exportSkewness
      module procedure bsa_exportSkewness_compute_
      module procedure bsa_exportSkewness_nocompute_
   end interface



   !**************************************************************************************
   !   INTERAFCE FOR PRIVATE PROCEDURES
   !**************************************************************************************
   interface
      module subroutine mainClassic_(m2mf_cls, m2mr_cls, m2o2mr_cls, m3mf_cls, m3mr_cls)
         real(RDP), allocatable :: m2mf_cls(:), m2mr_cls(:), m2o2mr_cls(:), m3mf_cls(:), m3mr_cls(:)
      end subroutine
      

      module subroutine mainMesher_(m3mf_msh, m3mr_msh)
         real(RDP), target, allocatable :: m3mf_msh(:), m3mr_msh(:)
      end subroutine
   end interface





   !**************************************************************************************
   !   INTERAFCE FOR  PUBLIC  PROCEDURES
   !**************************************************************************************
   interface



      ! --------------------------         GENERAL       ---------------------------------


      module subroutine bsa_setOutputDirectory(dirname)
         character(len=*), intent(in) :: dirname
      end subroutine


      module subroutine bsa_setOutFileName(fname)
         character(len=*), intent(in) :: fname
      end subroutine


      module subroutine bsa_setOutUnit(iunit)
         integer(kind = 4), intent(in) :: iunit
      end subroutine
      

      module subroutine bsa_closeUnitsAtEnd()
      end subroutine


      module subroutine bsa_setExportFileFormat(iform)
         integer(kind = 4), intent(in) :: iform
      end subroutine


      module subroutine bsa_setExportAppendMode(imode)
         integer(kind = 4), intent(in) :: imode
      end subroutine


      module subroutine bsa_setSpatialSymmetry(isym)
         integer(kind = 4), intent(in) :: isym
      end subroutine


      module subroutine bsa_setBfmMLR(bool)
         logical, intent(in) :: bool
      end subroutine


      module subroutine bsa_setPremeshType(itype)
         integer(kind = 4), intent(in) :: itype
      end subroutine


      module subroutine bsa_setPremeshMode(imode)
         integer(kind = 4), intent(in) :: imode
      end subroutine


      module subroutine bsa_doValidateModalData(bool)
         logical, intent(in) :: bool
      end subroutine


      ! module subroutine bsa_doValidateZoneDeltas(bool)
      !    logical, intent(in) :: bool
      ! end subroutine

      module subroutine bsa_setValidateDeltasPolicy(id)
         integer, intent(in) :: id
      end subroutine


      module subroutine bsa_setValidateDeltasValues(ibkg, ires)
         integer, intent(in) :: ibkg, ires
      end subroutine


      module subroutine bsa_Init()
      end subroutine


      module subroutine bsa_forceBsaClsExecution(bool)
         logical, intent(in) :: bool
      end subroutine


      module subroutine bsa_setMaxBkgPeakRestriction(bool)
         logical, intent(in) :: bool
      end subroutine


      module subroutine bsa_setPODTruncationThreshold(rval)
         real(kind = 8), intent(in) :: rval
      end subroutine


      module subroutine bsa_Run(m2mf_cls, m2mr_cls, m2o2mr_cls, m3mf_msh, m3mr_msh, m3mf_cls, m3mr_cls)
         real(RDP), target, allocatable, dimension(:) :: &
            m2mf_cls, m2mr_cls, m2o2mr_cls, m3mf_msh, m3mr_msh, m3mf_cls, m3mr_cls
      end subroutine


      module subroutine bsa_Finalise()
      end subroutine


      logical pure module function bsa_isCleaned()
      end function




      ! --------------------------         SETTINGS       ---------------------------------


      module elemental function bsa_isFullComp() result(bool)
         logical :: bool
      end function


      module subroutine bsa_setSubanType(isuban)
         integer(kind = 4), intent(in) :: isuban
      end subroutine 



      module subroutine bsa_setVersion(ivers)
         integer(kind = 4), intent(in) :: ivers
      end subroutine 



      module subroutine bsa_setScalingConv(iconv)
         integer(kind = 4), intent(in) :: iconv
      end subroutine 



      module subroutine bsa_setSpectraComputation(ipsd, ibisp)
         integer(kind = 4), intent(in), optional :: ipsd, ibisp
      end subroutine 



      module subroutine bsa_setSpectraExtension(ionlydiag)
         integer(kind = 4), intent(in) :: ionlydiag
      end subroutine 
   


      module subroutine bsa_setTestMode(itest)
         integer(kind = 4), intent(in) :: itest
      end subroutine 
   

      module subroutine bsa_setSymmetries(ibispsym, i3dsym)
         integer(kind = 4), intent(in) :: ibispsym, i3dsym
      end subroutine 

   
      module subroutine bsa_setupClassic(nfreqs, df)
         integer(kind = 4), intent(in) :: nfreqs
         real(RDP), intent(in) :: df
      end subroutine 
   


      module subroutine bsa_setupMesher(isvd, bkgrfmt, bkgaext, genpaext, maxaext, ifcov, idumpmod)
         integer(kind = 4), intent(in) :: isvd, bkgrfmt, maxaext
         integer(kind = 4), intent(in) :: bkgaext, genpaext, ifcov, idumpmod
      end subroutine 
   
   
   

      ! --------------------------         WIND       ---------------------------------


      module subroutine bsa_setWindDirections(dirs, ndirs)
         integer(kind = 4), intent(in) :: dirs(:)
         integer(kind = 4), intent(in), optional :: ndirs
      end subroutine


      module subroutine bsa_setWindTurbComps(tc, ntc)
         integer(kind = 4), intent(in) :: tc(:)
         integer(kind = 4), intent(in), optional :: ntc
      end subroutine
   
   
      module subroutine bsa_setWindVertProf(iwprof)
         integer(kind = 4), intent(in) :: iwprof
      end subroutine 
   

   
      module subroutine bsa_setPSDType(ipsd)
         integer(kind = 4), intent(in) :: ipsd
      end subroutine 
   

   
      module subroutine bsa_setWindAltDir(ivert)
         integer(kind = 4), intent(in) :: ivert
      end subroutine 
   

   
      module subroutine bsa_setWindZoneLimits(lim, ilim)
         real(RDP), intent(in) :: lim(..)
         integer(kind = 4), intent(in), optional :: ilim(..)
      end subroutine 
   

   
      module subroutine bsa_setAirDensity(aird)
         real(RDP), intent(in) :: aird
      end subroutine 
   

   
      module subroutine bsa_setGlobalRotMatW2G(rotW2G)
         real(RDP), intent(in) :: rotW2G(3, 3)
      end subroutine 
   

   
      module subroutine bsa_setWZMeanWindVel(mat)
         real(RDP), target, intent(in) :: mat(:)
      end subroutine 
   

   
      module subroutine bsa_setWZRefAlt(Zref)
         real(RDP), target, intent(in) :: Zref(:)
      end subroutine 
   

   
      module subroutine bsa_setTurbWindScales(L)
         real(RDP), target, intent(in) :: L(3, 3, *)
      end subroutine 
   

   
      module subroutine bsa_setTurbWindSDT(sigma)
         real(RDP), target, intent(in) :: sigma(3, *)
      end subroutine 
   

   
      module subroutine bsa_setWindCorrCoeffs(ccoeffs)
         real(RDP), target, intent(in) :: ccoeffs(3, 3, *)
      end subroutine 
   

   
      module subroutine bsa_setWindCorrExpnts(cexpn)
         real(RDP), target, intent(in) :: cexpn(3, 3, *)
      end subroutine 
   

   
      module subroutine bsa_setIncidenceAngles(incang)
         real(RDP), target, intent(in) :: incang(:)
      end subroutine 
   

   
      module subroutine bsa_setWZRotMatW2G(rotW2G_L)
         real(RDP), target, intent(in) :: rotW2G_L(3, 3, *)
      end subroutine 



      module subroutine bsa_setNodalVel(Unod)
         real(RDP), target, intent(in) :: Unod(:)
      end subroutine
   


      module subroutine bsa_setNodalWindZones(NodWZ)
         integer(kind = 4), target, intent(in) :: NodWZ(:)
      end subroutine


      module subroutine bsa_setNodalWindAltitudes(WnodAlt)
         real(RDP), target, intent(in) :: WnodAlt(:)
      end subroutine


      module subroutine bsa_setSpatialNodalCorr(nodCorr)
         real(RDP), target, intent(in) :: nodCorr(:, :)
      end subroutine



      module subroutine bsa_setWindFCoeffs(wfc)
         !> Dimensions should be [nlibs_l, ndegw+3, nnodes_l]
         real(RDP), target, intent(in) :: wfc(:, :, :)
      end subroutine


      module subroutine bsa_setPhitimesC(phiTc)
         real(RDP), target, intent(in) :: phiTc(:, :, :)
      end subroutine

   



      ! --------------------------         STRUCTURAL       ---------------------------------
   
      module subroutine bsa_setNodalCoords(nn, coords)
         integer(kind = 4), intent(in)   :: nn
         real(RDP), target, allocatable  :: coords(:, :)
      end subroutine
   

   
      module subroutine bsa_setNodalNOfDOFs(nlibs)
         integer(kind = 4), intent(in) :: nlibs
      end subroutine 



      module subroutine bsa_setTotalNOfNodes(nn)
         integer(kind = 4), intent(in) :: nn
      end subroutine

   
   
      module subroutine bsa_setLoadedNodalDOFs(libs_l, nlibs_l)
         integer(kind = 4), intent(in), target, allocatable :: libs_l(:)
         integer(kind = 4), intent(in), optional :: nlibs_l
      end subroutine 
   
   
   
      module subroutine bsa_setLoadedNodes(nodes_l, nn_l)
         integer(kind = 4), intent(in), target, allocatable :: nodes_l(:)
         integer(kind = 4), intent(in), optional :: nn_l
      end subroutine
   


      module subroutine bsa_setModalInfo(ndofs, nm, Phi, natf)
         integer(kind = 4), intent(in) :: ndofs, nm
         real(RDP), intent(in), target :: Phi(ndofs, nm), natf(nm)
      end subroutine 



      module subroutine bsa_setKeptModalShapes(modes)
         integer(kind = 4), intent(in) :: modes(:)
      end subroutine
   

   
      module subroutine bsa_setModalMatrices(nm, Mgen, Kgen, Cgen)
         integer(kind = 4), intent(in) :: nm
         real(RDP), intent(in), target, dimension(nm) :: Mgen, Kgen
         real(RDP), intent(in), target :: Cgen(nm, nm)
      end subroutine 
   

   
      module subroutine bsa_setTotDamping(xsi)
         real(RDP), target, intent(in) :: xsi(:)
      end subroutine 



      module pure function bsa_getUsedModeShapes() result(modes)
         integer(kind = 4), allocatable :: modes(:)
      end function





      ! --------------------------         COMPUTING       ---------------------------------

      module subroutine bsa_computeBRdecomp(m2mf, bkg, res)
         real(RDP), intent(in)  :: m2mf(:)
         real(RDP), allocatable, intent(out) :: bkg(:), res(:)
      end subroutine


      module subroutine bsa_computePeakFactors(&
            m2, m2o2, obs_time, peak_g, sk, peak_ng_pos, peak_ng_neg)
         real(kind = 8), intent(in)  :: m2(:), m2o2(:)
         real(kind = 8), intent(in)  :: obs_time
         real(kind = 8), allocatable, intent(inout) :: peak_g(:)
         real(kind = 8), intent(in), allocatable    :: sk(:)
         real(kind = 8), allocatable, intent(inout) :: peak_ng_pos(:)
         real(kind = 8), allocatable, intent(inout), optional :: peak_ng_neg(:)
      end subroutine


      



      ! --------------------------         EXPORTING       ---------------------------------

      module subroutine bsa_setExportDirectory(dirname)
         character(len = *), intent(in) :: dirname
      end subroutine



      module subroutine bsa_setExportInCurrDir()
      end subroutine



      module subroutine bsa_exportBR_nocompute_(fname, bkg, res, xsi)
         character(len = *), intent(in) :: fname
         real(RDP), intent(in) :: bkg(:), res(:), xsi(:)
      end subroutine



      module subroutine bsa_exportMomentToFile(fname, vec)
         character(len = *), intent(in) :: fname
         real(RDP), intent(in)          :: vec(:)
      end subroutine



      module subroutine bsa_exportSkewness_nocompute_(fname, sk)
         character(len = *), intent(in)  :: fname
         real(RDP), intent(in)  :: sk(:)
      end subroutine


      module subroutine bsa_exportSkewness_compute_(fname, m2, m3)
         character(len = *), intent(in)  :: fname
         real(RDP), intent(in)  :: m2(:), m3(:)
      end subroutine



      module subroutine bsa_exportPSDToFile(fname, psd, varname, f)
         character(len = *), intent(in) :: fname
         character(len = *), intent(in), optional :: varname
         real(RDP), intent(in), optional :: f(:)
         real(RDP), intent(in) :: psd(:, :)
      end subroutine



      module subroutine bsa_exportBispToFile(fname, bisp, varname)
         character(len = *), intent(in) :: fname
         character(len = *), intent(in), optional :: varname
         real(RDP), intent(in) :: bisp(:, :, :)
      end subroutine



      module subroutine bsa_saveCoordinatesToFile(fname, coords)
         character(len = *), intent(in)  :: fname
         real(RDP), intent(in), target, optional :: coords(:, :)
      end subroutine



      module subroutine bsa_exportPeakOrExtremesToFile(fname, rvar)
         character(len = *), intent(in) :: fname
         real(RDP), intent(in) :: rvar(:)
      end subroutine



      module subroutine bsa_setBRMExportDefaultMode(imode)
         integer(kind = 4), intent(in) :: imode
      end subroutine


      module subroutine bsa_setBRMExportFunction(fptr)
#ifdef __BSA_OMP
         procedure(exportBRMinterf_vect_all_), pointer, intent(in) :: fptr
#else
         procedure(exportBRMinterf_scalar_),   pointer, intent(in) :: fptr
#endif
      end subroutine


   end interface


end module BsaLib