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
module BsaLib_Structure

#include "../precisions"
   
   implicit none
   private


   type, public :: StructureModalData_t

      !> n. of (in memory) kept structure vibration modes
      integer(kind = 4) :: nm_ = 0

      !> n. of "usable" structure vibration modes
      !> To account for non 1-normalised modes, if any.
      integer(kind = 4) :: nm_eff_ = 0

      !> List of "usable" structural vibration modes.
      !> They might be less than the actual computed, 
      !> since some of them might refer to torsional
      !> modes, etc..
      integer(kind = 4), allocatable :: modes_(:)

      !> structural natural frequencies
      real(RDP), dimension(:), pointer    :: nat_freqs_;

      !> generalised modal matrix
      real(RDP), dimension(:, :), pointer :: phi_;

      !> modal damping ratios
      real(RDP), dimension(:), pointer :: xsi_;

      !> generalised mass matrix
      real(RDP), dimension(:), pointer :: Mm_;

      !> generalised damping matrix (Rayleigh)
      real(RDP), dimension(:, :), pointer :: Cm_;

      !> generalised stiffness matrix
      real(RDP), dimension(:), pointer :: Km_;
   end type StructureModalData_t



   type, public :: StructureData_t

      !> n. of all nodes
      integer(kind = 4) :: nn_ = 0

      !> n. of all libs (per node)
      integer(kind = 4) :: nlibs_ = 0

      !> n. of total DOFs
      integer(kind = 4) :: ndofs_ = 0
      
      !> n. of actually loaded nodes
      integer(kind = 4) :: nn_load_ = 0

      !> n. of actually loaded DOFs (per loaded node)
      integer(kind = 4) :: nlibs_load_ = 0

      !> list of actual loaded nodes
      integer(kind = 4), pointer :: n_load_(:) => null();

      !> list of actual loaded DOFs (per node)
      integer(kind = 4), pointer :: libs_load_(:) => null();

      ! NOTE: pointer since we do not want to copy.
      !> Nodal coordinates.
      real(RDP), pointer :: coords_(:, :) => null()

      !> modal structure info
      type(StructureModalData_t) :: modal_

      
      ! NOTE: following allocatables since are local instances !

      !> structure time scales
      real(RDP), dimension(:), allocatable :: str_time_scales_

      !> background peak widths
      real(RDP), dimension(:, :), allocatable :: bkg_peak_width_

      !> resonant peak widths
      real(RDP), dimension(:), allocatable :: res_peak_width_

   contains

      procedure, public, pass :: SetNodalCoords
      procedure, public, pass :: SetNOfNodalDOFs
      procedure, public, pass :: SetTotalNOfNodes
      procedure, public, pass :: SetLoadedNodalDOFs
      procedure, public, pass :: SetLoadedNodes
      procedure, public, pass :: SetModalInfo
      procedure, public, pass :: SetKeptModes
      procedure, public, pass :: SetKeptModesDefault
      procedure, public, pass :: SetModalMatrices
      procedure, public, pass :: SetTotDamping
      procedure, public, pass :: ComputeResPeakWidths
      procedure, public, pass :: computeBKGPeakWidths
      procedure, public, pass :: clean
   end type StructureData_t





   interface

      module subroutine SetNodalCoords(this, nn, coords)
         class(StructureData_t), intent(inout) :: this
         integer(kind = 4), intent(in)         :: nn
         real(RDP), target, allocatable        :: coords(:, :)
      end subroutine

      module subroutine SetNOfNodalDOFs(this, nlibs)
         class(StructureData_t), intent(inout) :: this
         integer(kind = 4), intent(in) :: nlibs
      end subroutine


      module subroutine SetTotalNOfNodes(this, nn)
         class(StructureData_t), intent(inout) :: this
         integer(kind = 4), intent(in) :: nn
      end subroutine



      module subroutine SetLoadedNodalDOFs(this, nlib, lib)
         class(StructureData_t), intent(inout) :: this
         integer(kind = 4), intent(in) :: nlib
         integer(kind = 4), target, intent(in) :: lib(:)
      end subroutine



      module subroutine SetLoadedNodes(this, nnl, nl)
         class(StructureData_t), intent(inout) :: this
         integer(kind = 4), intent(in) :: nnl
         integer(kind = 4), target, intent(in) :: nl(:)
      end subroutine






      module subroutine SetModalInfo(this, ndofs, nm, Phi, natf)
         class(StructureData_t), intent(inout) :: this
         integer(kind = 4), intent(in) :: ndofs, nm
         real(RDP), intent(in), target :: Phi(ndofs, nm), natf(nm)
      end subroutine


      module subroutine SetKeptModes(this, modes)
         class(StructureData_t), intent(inout) :: this
         integer(kind = 4), intent(in) :: modes(:)
      end subroutine


      !> NOTE: assumes that modes are not allocated yet.
      module subroutine SetKeptModesDefault(this)
         class(StructureData_t), intent(inout) :: this
      end subroutine



      module subroutine SetModalMatrices(this, nm, Mg, Kg, Cg)
         class(StructureData_t), intent(inout) :: this
         integer(kind = 4), intent(in) :: nm
         real(RDP), intent(in), target, dimension(nm) :: Mg, Kg
         real(RDP), intent(in), target :: Cg(nm, nm)
      end subroutine



      module subroutine SetTotDamping(this, xsi)
         class(StructureData_t), intent(inout) :: this
         real(RDP), target, intent(in) :: xsi(this%modal_%nm_)
      end subroutine






      module subroutine ComputeResPeakWidths(this)
         class(StructureData_t), intent(inout) :: this
      end subroutine



      module subroutine computeBKGPeakWidths(this, wind_scales)
         class(StructureData_t), intent(inout) :: this
         real(RDP), intent(in) :: wind_scales(:, :)
      end subroutine




      module subroutine clean(this)
         class(StructureData_t) :: this
      end subroutine

   end interface


end module BsaLib_Structure