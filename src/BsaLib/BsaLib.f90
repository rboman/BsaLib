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
         real(bsa_real_t), allocatable :: &
            m2mf_cls(:), m2mr_cls(:), m2o2mr_cls(:), m3mf_cls(:), m3mr_cls(:)
      end subroutine
      

      module subroutine mainMesher_(m3mf_msh, m3mr_msh)
         real(bsa_real_t), target, allocatable :: m3mf_msh(:), m3mr_msh(:)
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
         integer(bsa_int_t), intent(in) :: iunit
      end subroutine
      

      module subroutine bsa_closeUnitsAtEnd()
      end subroutine


      module subroutine bsa_setExportFileFormat(iform)
         integer(bsa_int_t), intent(in) :: iform
      end subroutine


      module subroutine bsa_setExportAppendMode(imode)
         integer(bsa_int_t), intent(in) :: imode
      end subroutine


      module subroutine bsa_setSpatialSymmetry(isym)
         integer(bsa_int_t), intent(in) :: isym
      end subroutine


      module subroutine bsa_setBfmMLR(bool)
         logical, intent(in) :: bool
      end subroutine


      module subroutine bsa_setPremeshType(itype)
         integer(bsa_int_t), intent(in) :: itype
      end subroutine


      module subroutine bsa_setPremeshMode(imode)
         integer(bsa_int_t), intent(in) :: imode
      end subroutine


      module subroutine bsa_doValidateModalData(bool)
         logical, intent(in) :: bool
      end subroutine


      ! module subroutine bsa_doValidateZoneDeltas(bool)
      !    logical, intent(in) :: bool
      ! end subroutine

      module subroutine bsa_setValidateDeltasPolicy(id)
         integer(bsa_int_t), intent(in) :: id
      end subroutine


      module subroutine bsa_setValidateDeltasValues(ibkg, ires)
         integer(bsa_int_t), intent(in) :: ibkg, ires
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
         real(real64), intent(in) :: rval
      end subroutine


      module subroutine bsa_Run(m2mf_cls, m2mr_cls, m2o2mr_cls, m3mf_msh, m3mr_msh, m3mf_cls, m3mr_cls)
         real(bsa_real_t), target, allocatable, dimension(:) :: &
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
         integer(bsa_int_t), intent(in) :: isuban
      end subroutine 



      module subroutine bsa_setVersion(ivers)
         integer(bsa_int_t), intent(in) :: ivers
      end subroutine 



      module subroutine bsa_setScalingConv(iconv)
         integer(bsa_int_t), intent(in) :: iconv
      end subroutine 



      module subroutine bsa_setSpectraComputation(ipsd, ibisp)
         integer(bsa_int_t), intent(in), optional :: ipsd, ibisp
      end subroutine 



      module subroutine bsa_setSpectraExtension(ionlydiag)
         integer(bsa_int_t), intent(in) :: ionlydiag
      end subroutine 
   


      module subroutine bsa_setTestMode(itest)
         integer(bsa_int_t), intent(in) :: itest
      end subroutine 
   

      module subroutine bsa_setSymmetries(ibispsym, i3dsym)
         integer(bsa_int_t), intent(in) :: ibispsym, i3dsym
      end subroutine 

   
      module subroutine bsa_setupClassic(nfreqs, df)
         integer(bsa_int_t), intent(in) :: nfreqs
         real(bsa_real_t), intent(in) :: df
      end subroutine 
   


      module subroutine bsa_setupMesher(isvd, bkgrfmt, bkgaext, genpaext, maxaext, ifcov, idumpmod)
         integer(bsa_int_t), intent(in) :: isvd, bkgrfmt, maxaext
         integer(bsa_int_t), intent(in) :: bkgaext, genpaext, ifcov, idumpmod
      end subroutine 
   
   
   

      ! --------------------------         WIND       ---------------------------------


      module subroutine bsa_setWindDirections(dirs, ndirs)
         integer(bsa_int_t), intent(in) :: dirs(:)
         integer(bsa_int_t), intent(in), optional :: ndirs
      end subroutine


      module subroutine bsa_setWindTurbComps(tc, ntc)
         integer(bsa_int_t), intent(in) :: tc(:)
         integer(bsa_int_t), intent(in), optional :: ntc
      end subroutine
   
   
      module subroutine bsa_setWindVertProf(iwprof)
         integer(bsa_int_t), intent(in) :: iwprof
      end subroutine 
   

   
      module subroutine bsa_setPSDType(ipsd)
         integer(bsa_int_t), intent(in) :: ipsd
      end subroutine 
   

   
      module subroutine bsa_setWindAltDir(ivert)
         integer(bsa_int_t), intent(in) :: ivert
      end subroutine 
   

   
      module subroutine bsa_setWindZoneLimits(lim, ilim)
         real(bsa_real_t), intent(in) :: lim(..)
         integer(bsa_int_t), intent(in), optional :: ilim(..)
      end subroutine 
   

   
      module subroutine bsa_setAirDensity(aird)
         real(bsa_real_t), intent(in) :: aird
      end subroutine 
   

   
      module subroutine bsa_setGlobalRotMatW2G(rotW2G)
         real(bsa_real_t), intent(in) :: rotW2G(3, 3)
      end subroutine 
   

   
      module subroutine bsa_setWZMeanWindVel(mat)
         real(bsa_real_t), target, intent(in) :: mat(:)
      end subroutine 
   

   
      module subroutine bsa_setWZRefAlt(Zref)
         real(bsa_real_t), target, intent(in) :: Zref(:)
      end subroutine 
   

   
      module subroutine bsa_setTurbWindScales(L)
         real(bsa_real_t), target, intent(in) :: L(3, 3, *)
      end subroutine 
   

   
      module subroutine bsa_setTurbWindSDT(sigma)
         real(bsa_real_t), target, intent(in) :: sigma(3, *)
      end subroutine 
   

   
      module subroutine bsa_setWindCorrCoeffs(ccoeffs)
         real(bsa_real_t), target, intent(in) :: ccoeffs(3, 3, *)
      end subroutine 
   

   
      module subroutine bsa_setWindCorrExpnts(cexpn)
         real(bsa_real_t), target, intent(in) :: cexpn(3, 3, *)
      end subroutine 
   

   
      module subroutine bsa_setIncidenceAngles(incang)
         real(bsa_real_t), target, intent(in) :: incang(:)
      end subroutine 
   

   
      module subroutine bsa_setWZRotMatW2G(rotW2G_L)
         real(bsa_real_t), target, intent(in) :: rotW2G_L(3, 3, *)
      end subroutine 



      module subroutine bsa_setNodalVel(Unod)
         real(bsa_real_t), target, intent(in) :: Unod(:)
      end subroutine
   


      module subroutine bsa_setNodalWindZones(NodWZ)
         integer(bsa_int_t), target, intent(in) :: NodWZ(:)
      end subroutine


      module subroutine bsa_setNodalWindAltitudes(WnodAlt)
         real(bsa_real_t), target, intent(in) :: WnodAlt(:)
      end subroutine


      module subroutine bsa_setSpatialNodalCorr(nodCorr)
         real(bsa_real_t), target, intent(in) :: nodCorr(:, :)
      end subroutine



      module subroutine bsa_setWindFCoeffs(wfc)
         !> Dimensions should be [nlibs_l, ndegw+3, nnodes_l]
         real(bsa_real_t), target, intent(in) :: wfc(:, :, :)
      end subroutine


      module subroutine bsa_setPhitimesC(phiTc)
         real(bsa_real_t), target, intent(in) :: phiTc(:, :, :)
      end subroutine

   



      ! --------------------------         STRUCTURAL       ---------------------------------
   
      module subroutine bsa_setNodalCoords(nn, coords)
         integer(bsa_int_t), intent(in)   :: nn
         real(bsa_real_t), target, allocatable  :: coords(:, :)
      end subroutine
   

   
      module subroutine bsa_setNodalNOfDOFs(nlibs)
         integer(bsa_int_t), intent(in) :: nlibs
      end subroutine 



      module subroutine bsa_setTotalNOfNodes(nn)
         integer(bsa_int_t), intent(in) :: nn
      end subroutine

   
   
      module subroutine bsa_setLoadedNodalDOFs(libs_l, nlibs_l)
         integer(bsa_int_t), intent(in), target, allocatable :: libs_l(:)
         integer(bsa_int_t), intent(in), optional :: nlibs_l
      end subroutine 
   
   
   
      module subroutine bsa_setLoadedNodes(nodes_l, nn_l)
         integer(bsa_int_t), intent(in), target, allocatable :: nodes_l(:)
         integer(bsa_int_t), intent(in), optional :: nn_l
      end subroutine
   


      module subroutine bsa_setModalInfo(ndofs, nm, Phi, natf)
         integer(bsa_int_t), intent(in) :: ndofs, nm
         real(bsa_real_t), intent(in), target :: Phi(ndofs, nm), natf(nm)
      end subroutine 



      module subroutine bsa_setKeptModalShapes(modes)
         integer(bsa_int_t), intent(in) :: modes(:)
      end subroutine
   

   
      module subroutine bsa_setModalMatrices(nm, Mgen, Kgen, Cgen)
         integer(bsa_int_t), intent(in) :: nm
         real(bsa_real_t), intent(in), target, dimension(nm) :: Mgen, Kgen
         real(bsa_real_t), intent(in), target :: Cgen(nm, nm)
      end subroutine 
   

   
      module subroutine bsa_setTotDamping(xsi)
         real(bsa_real_t), target, intent(in) :: xsi(:)
      end subroutine 



      module pure function bsa_getUsedModeShapes() result(modes)
         integer(bsa_int_t), allocatable :: modes(:)
      end function





      ! --------------------------         COMPUTING       ---------------------------------

      module subroutine bsa_computeBRdecomp(m2mf, bkg, res)
         real(bsa_real_t), intent(in)  :: m2mf(:)
         real(bsa_real_t), allocatable, intent(out) :: bkg(:), res(:)
      end subroutine


      module subroutine bsa_computePeakFactors(&
            m2, m2o2, obs_time, peak_g, sk, peak_ng_pos, peak_ng_neg)
         real(bsa_real_t), intent(in)  :: m2(:), m2o2(:)
         real(bsa_real_t), intent(in)  :: obs_time
         real(bsa_real_t), allocatable, intent(inout) :: peak_g(:)
         real(bsa_real_t), intent(in), allocatable    :: sk(:)
         real(bsa_real_t), allocatable, intent(inout) :: peak_ng_pos(:)
         real(bsa_real_t), allocatable, intent(inout), optional :: peak_ng_neg(:)
      end subroutine


      



      ! --------------------------         EXPORTING       ---------------------------------

      module subroutine bsa_setExportDirectory(dirname)
         character(len = *), intent(in) :: dirname
      end subroutine



      module subroutine bsa_setExportInCurrDir()
      end subroutine



      module subroutine bsa_exportBR_nocompute_(fname, bkg, res, xsi)
         character(len = *), intent(in) :: fname
         real(bsa_real_t), intent(in) :: bkg(:), res(:), xsi(:)
      end subroutine



      module subroutine bsa_exportMomentToFile(fname, vec)
         character(len = *), intent(in) :: fname
         real(bsa_real_t), intent(in)          :: vec(:)
      end subroutine



      module subroutine bsa_exportSkewness_nocompute_(fname, sk)
         character(len = *), intent(in)  :: fname
         real(bsa_real_t), intent(in)  :: sk(:)
      end subroutine


      module subroutine bsa_exportSkewness_compute_(fname, m2, m3)
         character(len = *), intent(in)  :: fname
         real(bsa_real_t), intent(in)  :: m2(:), m3(:)
      end subroutine



      module subroutine bsa_exportPSDToFile(fname, psd, varname, f)
         character(len = *), intent(in) :: fname
         character(len = *), intent(in), optional :: varname
         real(bsa_real_t), intent(in), optional :: f(:)
         real(bsa_real_t), intent(in) :: psd(:, :)
      end subroutine



      module subroutine bsa_exportBispToFile(fname, bisp, varname)
         character(len = *), intent(in) :: fname
         character(len = *), intent(in), optional :: varname
         real(bsa_real_t), intent(in) :: bisp(:, :, :)
      end subroutine



      module subroutine bsa_saveCoordinatesToFile(fname, coords)
         character(len = *), intent(in)  :: fname
         real(bsa_real_t), intent(in), target, optional :: coords(:, :)
      end subroutine



      module subroutine bsa_exportPeakOrExtremesToFile(fname, rvar)
         character(len = *), intent(in) :: fname
         real(bsa_real_t), intent(in) :: rvar(:)
      end subroutine



      module subroutine bsa_setBRMExportDefaultMode(imode)
         integer(bsa_int_t), intent(in) :: imode
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