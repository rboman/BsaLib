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
module data

   use BsaLib
   implicit none
   public
   integer(int32), parameter :: IUN_BSADATA = 22222
   integer(int32), parameter :: IUN_FINDATA = 22223
   logical :: l_formmode = .false.
   logical :: fin_data_read_ = .false.
   logical :: bsa_data_read_ = .false.
   character(len = :), allocatable :: FINFILE_FNAME_


   integer(int32)   :: i_suban, i_vers, i_defsc, i_psd, i_bisp, i_onlyd, i_test
   integer(int32)   :: i_bispsym, i_3dsym, i_nfreqs
   real(bsa_real_t) :: r_df
   integer(int32)   :: i_svd, i_bkgrfmt, i_bkgaext, i_genpaext, i_maxaext, i_fcov, i_dumpmod

   integer(int32) :: i_ntc, i_ndirs, tc(3), dirs(3)

   integer(int32) :: i_nnodes, i_nlibs, i_nnodesl, i_nlibsl
   integer(int32), target, allocatable   :: nodesl(:), libsl(:)
   real(bsa_real_t), target, allocatable :: nod_cords(:, :)

   integer(int32)   :: i_varu, i_su, i_vert, i_degw
   integer(int32)   :: i_nzones
   real(bsa_real_t) :: r_aird, r_rotW2G(3, 3)
   real(bsa_real_t), target, allocatable :: r_Zref_z(:), r_UBref_z(:), r_alph_z(:), r_lims_z(:)
   real(bsa_real_t), target, allocatable :: r_L_z(:, :, :), r_std_z(:, :), r_corrC_z(:, :, :)
   real(bsa_real_t), target, allocatable :: r_corrEx_z(:, :, :), r_rotW2G_z(:, :, :), r_incang_z(:)

   integer(int32), target, allocatable   :: i_wzNod(:)
   real(bsa_real_t), target, allocatable :: r_wAltNod(:), r_UBnod(:), r_corrNod(:, :)

   real(bsa_real_t), target, allocatable :: r_wfc(:, :, :)

   integer(int32) :: i_nm, i_ndofs
   real(bsa_real_t), target, allocatable :: r_natf(:), r_modm(:, :)
   real(bsa_real_t), target, allocatable :: r_Mg(:), r_Kg(:), r_Cg(:, :)
   real(bsa_real_t), target, allocatable :: r_xsist(:), r_xsiad(:)

   logical :: use_custom_damping_ = .false.
   real(bsa_real_t), allocatable :: custom_damp_val_(:)

   integer(int32) :: i_exprt_mode_ = BSA_EXPORT_MODE_APPEND
   integer(int32) :: i_exprt_form_ = BSA_EXPORT_FORMAT_FORMATTED
   logical :: export_results_to_files_ = .true.
end module data



!=========================================================================================
!=========================================================================================
! MAIN PROGRAM
!=========================================================================================
!=========================================================================================
program bsa

   use data

   implicit none

   ! local, used to retrieve single run BSA results
   real(bsa_real_t), allocatable :: bkg_(:),  res_(:)
   real(bsa_real_t), allocatable :: m2mf_(:), m2mr_(:), m2o2mr_(:)   ! NOTE: implicitly classic
   real(bsa_real_t), allocatable :: m3mf_cls_(:), m3mr_cls_(:)
   real(bsa_real_t), allocatable :: m3mf_msh_(:), m3mr_msh_(:)


   real(bsa_real_t), allocatable, dimension(:) :: m2_r_diag, m3_r_diag, sk_r_diag, m2o2_r_diag
   real(bsa_real_t), allocatable, dimension(:) :: m2_r_full, m3_r_full, sk_r_full, m2o2_r_full

   real(bsa_real_t), allocatable :: peak_pos_r_diag_g(:), peak_pos_r_diag_ng(:)
   real(bsa_real_t), allocatable :: peak_pos_r_full_g(:), peak_pos_r_full_ng(:)
   real(bsa_real_t), allocatable :: extr_pos_r_diag_g(:), extr_pos_r_diag_ng(:)
   real(bsa_real_t), allocatable :: extr_pos_r_full_g(:), extr_pos_r_full_ng(:)
   
   real(bsa_real_t), allocatable :: peak_neg_r_diag_ng(:)
   real(bsa_real_t), allocatable :: peak_neg_r_full_ng(:)
   real(bsa_real_t), allocatable :: extr_neg_r_diag_ng(:)
   real(bsa_real_t), allocatable :: extr_neg_r_full_ng(:)

   character(len = *), parameter   :: exp_prfx = 'export_'
   character(len = *), parameter   :: cls_sffx = 'cls', msh_sffx = 'msh'
   character(len = *), parameter   :: exp_fext = '.txt', udscr = '_'
   character(len = :), allocatable :: fname, cmb_sffx

   
#ifdef __BSA_CL
#  define BSACL_SFFX_ //'_CL'
#else
# ifdef __BSA_CUDA
#   define BSACL_SFFX_ //'_CUDA'
#  else
#   define BSACL_SFFX_
# endif
#endif


!=========================================================================================
! MAIN BODY
!=========================================================================================
   call bsa_printBSAHeader()
   call printTool()
   call parseArgs()

   ! set defaults for entities not provided by user.
   if (.not. allocated(FINFILE_FNAME_)) FINFILE_FNAME_ = 'bsa.findata'
   call readDataFiles()

   call setup()

   ! NOTE: save in r_xsist the total, before pass it
   if (use_custom_damping_) then
      print '(1x, 2a, g, a)', &
         WARNMSG, 'Using custom damping of  ', custom_damp_val_(1), '  for all modes.'
      r_xsist(:) = custom_damp_val_(1)
   endif
   call bsa_setTotDamping(r_xsist)

   ! BUG: allow bsa_Run to accept already allocated entities
   !      (check for size match).
   call bsa_Run(m2mf_, m2mr_, m2o2mr_, m3mf_msh_, m3mr_msh_, m3mf_cls_, m3mr_cls_)

   ! POST-PROCESSING
   if (export_results_to_files_) then


      if (i_onlyd) then
         cmb_sffx = 'diag'  BSACL_SFFX_
      else
         cmb_sffx = 'full'  BSACL_SFFX_
      endif

      block
         integer(int32), allocatable :: i_modes(:)

         i_modes = bsa_getUsedModeShapes()

         if (allocated(m2mf_)) then
            call bsa_computeBRdecomp(m2mf_, bkg_, res_)
            fname = exp_prfx // 'm2_BR_decomp.txt'
            call bsa_exportBRdecomp(fname, bkg_, res_, r_xsist(i_modes))

            fname = exp_prfx // 'm2_mf_' // cls_sffx // udscr // cmb_sffx // exp_fext
            call bsa_exportMomentToFile(fname, m2mf_)

            if (allocated(m3mf_cls_)) &
               call bsa_exportSkewness(exp_prfx // 'sk_mf_' // cls_sffx // udscr // cmb_sffx // exp_fext, m2mf_, m3mf_cls_)
            if (allocated(m3mf_msh_)) &
               call bsa_exportSkewness(exp_prfx // 'sk_mf_' // msh_sffx // udscr // cmb_sffx // exp_fext, m2mf_, m3mf_msh_)
         endif


         if (allocated(m2mr_)) then
            
            fname = exp_prfx // 'm2_mr_' // cls_sffx // udscr // cmb_sffx // exp_fext
            call bsa_exportMomentToFile(fname, m2mr_)

            if (allocated(m3mr_msh_)) then
               fname = exp_prfx // 'sk_mr_' // msh_sffx // udscr // cmb_sffx // exp_fext
               call bsa_exportSkewness(fname, m2mr_, m3mr_msh_)
               call modalRecombination(r_modm(:, i_modes), m2mr_, m3mr_msh_, m2o2mr_)
            endif
            
            if (allocated(m3mr_cls_)) then
               fname = exp_prfx // 'sk_mr_' // cls_sffx // udscr // cmb_sffx // exp_fext
               call bsa_exportSkewness(fname, m2mr_, m3mr_cls_)
               if (.not. allocated(m3mr_msh_)) &
                  call modalRecombination(r_modm(:, i_modes), m2mr_, m3mr_cls_, m2o2mr_)
            endif


            ! TODO: introduce proper T variable
            call bsa_computePeakFactors(m2_r_diag, m2o2_r_diag, 600.d0, peak_pos_r_diag_g, sk_r_diag, peak_pos_r_diag_ng, peak_neg_r_diag_ng)
            call bsa_computePeakFactors(m2_r_full, m2o2_r_full, 600.d0, peak_pos_r_full_g, sk_r_full, peak_pos_r_full_ng, peak_neg_r_full_ng)


            if (allocated(sk_r_full)) then
               fname = exp_prfx // 'sk_r_full' // exp_fext
               call bsa_exportSkewness(fname, sk_r_full)
            endif


            ! TODO: maybe check if file already exists
            call bsa_saveCoordinatesToFile('coordinates.txt')

            if (i_onlyd == 0) then
               if (allocated(peak_pos_r_diag_g))  then
                  fname = exp_prfx // 'peak_pos_r_' // cmb_sffx // 'D' // '_g' // exp_fext
                  call bsa_exportPeakOrExtremesToFile(fname, peak_pos_r_diag_g)
                  extr_pos_r_diag_g = peak_pos_r_diag_g  * sqrt(m2_r_diag)
                  fname = exp_prfx // 'extr_pos_r_' // cmb_sffx // 'D' // '_g' // exp_fext
                  call bsa_exportPeakOrExtremesToFile(fname, extr_pos_r_diag_g)
               endif

               if (allocated(peak_pos_r_diag_ng)) then
                  fname = exp_prfx // 'peak_pos_r_' // cmb_sffx // 'D' // '_ng' // exp_fext
                  call bsa_exportPeakOrExtremesToFile(fname, peak_pos_r_diag_ng)
                  extr_pos_r_diag_ng = peak_pos_r_diag_ng * sqrt(m2_r_diag)
                  fname = exp_prfx // 'extr_pos_r_' // cmb_sffx // 'D' // '_ng' // exp_fext
                  call bsa_exportPeakOrExtremesToFile(fname, extr_pos_r_diag_ng)
               endif
               if (allocated(peak_neg_r_diag_ng)) then
                  fname = exp_prfx // 'peak_neg_r_' // cmb_sffx // 'D' // '_ng' // exp_fext
                  call bsa_exportPeakOrExtremesToFile(fname, peak_neg_r_diag_ng)
                  extr_neg_r_diag_ng = peak_neg_r_diag_ng * sqrt(m2_r_diag)
                  fname = exp_prfx // 'extr_neg_r_' // cmb_sffx // 'D' // '_ng' // exp_fext
                  call bsa_exportPeakOrExtremesToFile(fname, extr_neg_r_diag_ng)
               endif

               if (allocated(peak_pos_r_full_g))  then
                  fname = exp_prfx // 'peak_pos_r_' // cmb_sffx // 'F' // '_g' // exp_fext
                  call bsa_exportPeakOrExtremesToFile(fname, peak_pos_r_full_g)
                  extr_pos_r_full_g  = peak_pos_r_full_g  * sqrt(m2_r_full)
                  fname = exp_prfx // 'extr_pos_r_' // cmb_sffx // 'F' // '_g' // exp_fext
                  call bsa_exportPeakOrExtremesToFile(fname, extr_pos_r_full_g)
               endif

               if (allocated(peak_pos_r_full_ng)) then
                  fname = exp_prfx // 'peak_pos_r_' // cmb_sffx // 'F' // '_ng' // exp_fext
                  call bsa_exportPeakOrExtremesToFile(fname, peak_pos_r_full_ng)
                  extr_pos_r_full_ng = peak_pos_r_full_ng * sqrt(m2_r_full)
                  fname = exp_prfx // 'extr_pos_r_' // cmb_sffx // 'F' // '_ng' // exp_fext
                  call bsa_exportPeakOrExtremesToFile(fname, extr_pos_r_full_ng)
               endif
               if (allocated(peak_neg_r_full_ng)) then
                  fname = exp_prfx // 'peak_neg_r_' // cmb_sffx // 'F' // '_ng' // exp_fext
                  call bsa_exportPeakOrExtremesToFile(fname, peak_neg_r_full_ng)
                  extr_neg_r_full_ng = peak_neg_r_full_ng * sqrt(m2_r_full)
                  fname = exp_prfx // 'extr_neg_r_' // cmb_sffx // 'F' // '_ng' // exp_fext
                  call bsa_exportPeakOrExtremesToFile(fname, extr_neg_r_full_ng)
               endif

            else

               if (allocated(peak_pos_r_diag_g))  then
                  fname = exp_prfx // 'peak_pos_r_' // cmb_sffx // '_g' // exp_fext
                  call bsa_exportPeakOrExtremesToFile(fname, peak_pos_r_diag_g)
                  extr_pos_r_diag_g  = peak_pos_r_diag_g  * sqrt(m2_r_diag)
                  fname = exp_prfx // 'extr_pos_r_' // cmb_sffx // '_g' // exp_fext
                  call bsa_exportPeakOrExtremesToFile(fname, extr_pos_r_diag_g)
               endif


               if (allocated(peak_pos_r_diag_ng)) then
                  fname = exp_prfx // 'peak_pos_r_' // cmb_sffx // '_ng' // exp_fext
                  call bsa_exportPeakOrExtremesToFile(fname, peak_pos_r_diag_ng)
                  extr_pos_r_diag_ng = peak_pos_r_diag_ng * sqrt(m2_r_diag)
                  fname = exp_prfx // 'extr_pos_r_' // cmb_sffx // '_ng' // exp_fext
                  call bsa_exportPeakOrExtremesToFile(fname, extr_pos_r_diag_ng)
               endif
               if (allocated(peak_neg_r_diag_ng)) then
                  fname = exp_prfx // 'peak_neg_r_' // cmb_sffx // '_ng' // exp_fext
                  call bsa_exportPeakOrExtremesToFile(fname, peak_neg_r_diag_ng)
                  extr_neg_r_diag_ng = peak_neg_r_diag_ng * sqrt(m2_r_diag)
                  fname = exp_prfx // 'extr_neg_r_' // cmb_sffx // '_ng' // exp_fext
                  call bsa_exportPeakOrExtremesToFile(fname, extr_neg_r_diag_ng)
               endif
            endif

         endif ! allocated m2mr
         deallocate(i_modes)

         if (allocated(m3mf_cls_)) then
            fname = exp_prfx // 'm3_mf_' // cls_sffx // udscr // cmb_sffx // exp_fext
            call bsa_exportMomentToFile(fname, m3mf_cls_)
         endif
         if (allocated(m3mr_cls_)) then
            fname = exp_prfx // 'm3_mr_' // cls_sffx // udscr // cmb_sffx // exp_fext
            call bsa_exportMomentToFile(fname, m3mr_cls_)
         endif
         if (allocated(m3mf_msh_)) then
            fname = exp_prfx // 'm3_mf_' // msh_sffx // udscr // cmb_sffx // exp_fext
            call bsa_exportMomentToFile(fname, m3mf_msh_)
         endif
         if (allocated(m3mr_msh_)) then
            fname = exp_prfx // 'm3_mr_' // msh_sffx // udscr // cmb_sffx // exp_fext
            call bsa_exportMomentToFile(fname, m3mr_msh_)
         endif
      end block
   endif

   if (allocated(m2mf_))      deallocate(m2mf_)
   if (allocated(m3mf_cls_))  deallocate(m3mf_cls_)
   if (allocated(m3mf_msh_))  deallocate(m3mf_msh_)
   if (allocated(m2mr_))      deallocate(m2mr_)
   if (allocated(m3mr_cls_))  deallocate(m3mr_cls_)
   if (allocated(m3mr_msh_))  deallocate(m3mr_msh_)
   if (allocated(m3mf_cls_))  deallocate(m3mf_cls_)
   if (allocated(m3mr_cls_))  deallocate(m3mr_cls_)
   if (allocated(m3mf_msh_))  deallocate(m3mf_msh_)
   if (allocated(m3mr_msh_))  deallocate(m3mr_msh_)

   if (allocated(peak_pos_r_diag_g))   deallocate(peak_pos_r_diag_g)
   if (allocated(peak_pos_r_diag_ng))  deallocate(peak_pos_r_diag_ng)
   if (allocated(peak_pos_r_full_g))   deallocate(peak_pos_r_full_g)
   if (allocated(peak_pos_r_full_ng))  deallocate(peak_pos_r_full_ng)
   if (allocated(extr_pos_r_diag_g))   deallocate(extr_pos_r_diag_g)
   if (allocated(extr_pos_r_diag_ng))  deallocate(extr_pos_r_diag_ng)
   if (allocated(extr_pos_r_full_g))   deallocate(extr_pos_r_full_g)
   if (allocated(extr_pos_r_full_ng))  deallocate(extr_pos_r_full_ng)


   call releaseMemory(0)

!=========================================================================================
!=========================================================================================
! END MAIN PROGRAM
!=========================================================================================
!=========================================================================================



contains ! utility procedures




   subroutine printTool()
      print *, '                          BSA  -  main program'
      print *, '   Bispectral Stochastic Analysis of MDOFs systems under stationary actions.'
      print *
   end subroutine


   subroutine usage()
      print *
      print *, ' Syntax:'
      print *, '    bsa.exe  [[options] file]'
      print *
      print *, ' options:'
      print *, '   --readmode, -readmode, /readmode  <val>'
      print *, '        valid  <val>  values:'
      print *, '           formatted, unformatted (default)'
      print *
      print *, '   --force-damping, -force-damping, /force-damping  <val>'
      print *, '        If specified, uses  <val>  damping for all modes.'
      print *, '        For multiple values (parametric analysis), separate'
      print *, '        values using semicolon, ex.  val1;val2;val3;etc..'
      print *, '        In such cases, BSA core is run multiple times.'
      print *
      print *, '   --no-append-exports, -no-append-exports, /no-append-exports'
      print *, '        Overrides instead of appending to existing export files.'
      print *
      print *, '   --export-binary, -export-binary, /export-binary'
      print *, '        Exports results using unformatted (binary) files.'
      print *
      print *, '   --no-export, -no-export, /no-export'
      print *, '        Do not export results to files.'
      print *
      print *, '   --out-file, -out-file, /out-file  <outfilename>'
      print *, '        If specified, uses  <outfilename>  as output file name.'
      print *
      print *, ' file: input file name.'
      print *, '       If none is passed, automatically searches it into'
      print *, '       the current directory (from where the program is invoked)'

      call releaseMemory(1)
   end subroutine




   function removeInputArgPlaceholder_(arg) result(arg_)
!DIR$ ATTRIBUTES FORCEINLINE :: removeInputArgPlaceholder_
      character(len = *), intent(in)  :: arg
      character(len = :), allocatable :: arg_
      character(len = *), parameter   :: MINUS = '-'

      if (arg(1:1) == '/') then

         arg_ = arg(2 : len_trim(arg))

      elseif (arg(1:1) == MINUS) then

         if (arg(2:2) == MINUS) then
            arg_ = arg(3 : len_trim(arg))
         else
            arg_ = arg(2 : len_trim(arg))
         endif

      else ! take it as such

         arg_ = arg(1 : len_trim(arg))
      endif
   end function





   subroutine parseArgs()
      !! Parse input arguments, if any.
      integer :: argc

      argc = command_argument_count()
      if (argc == 0) return

      block
         integer :: i, currargc, istat
         character(len = 64) :: arg
         character(len = :), allocatable :: arg_
      
         currargc = 1

         ! main parsing loop
         do
            call get_command_argument(currargc, arg)
            arg_ = removeInputArgPlaceholder_(arg)

            select case (arg_)

               case ('readmode')
                  currargc = currargc + 1
                  call get_command_argument(currargc, arg, status=istat)
                  if (istat /= 0) call usage()
                  
                  select case (arg(1:len_trim(arg)))
                     case ('formatted')
                        l_formmode = .true.
                     case ('unformatted')
                     case default
                        call usage()
                  end select

               case ('force-damping')
                  currargc = currargc + 1
                  call get_command_argument(currargc, arg, status=istat)
                  if (istat /= 0) call usage()
                  call getCustomXSIValues_(arg)


               case ('no-append-exports')
                  i_exprt_mode_ = BSA_EXPORT_MODE_REPLACE

               case ('export-binary')
                  i_exprt_form_ = BSA_EXPORT_FORMAT_UNFORMATTED

               case ('no-export')
                  export_results_to_files_ = .false.


               case ('out-file')
                  currargc = currargc + 1
                  call get_command_argument(currargc, arg, status=istat)
                  if (istat /= 0) call usage()
                  call bsa_setOutFileName(arg(1 : len_trim(arg)))


               case default ! input file
                  FINFILE_FNAME_ = arg_(1 : len_trim(arg_))

            end select

            if (currargc == argc) exit ! parsing finished
            
            ! go to next input arg
            currargc = currargc + 1
         enddo

      end block

   end subroutine parseArgs





   subroutine getCustomXSIValues_(arg)
      character(len = *), intent(in) :: arg
      integer :: iich, iech, lench, ixsi
      ! BUG: limits to 20 xsi values. Allow infinite.
      real(bsa_real_t) :: xsi_(20)

      lench = len_trim(arg)
      if (lench == 0) goto 99
      iich = 1
      iech = 1
      ixsi = 0
      do while (iech <= lench)
         if (arg(iech:iech) == ';') then
            ixsi = ixsi + 1
            read(unit=arg(iich : iech-1), fmt='(f)', err=99) xsi_(ixsi)
            iech = iech + 1
            iich = iech
         else
            iech = iech + 1
         endif
      enddo
      ! NOTE: treat last one (only one if no semicolon found)
      ixsi = ixsi + 1
      read(unit=arg(iich : iech), fmt='(f)', err=99) xsi_(ixsi)

      ! if we get here, correctly read
      use_custom_damping_ = .true.
      custom_damp_val_    = xsi_(1 : ixsi)
      return

      99 print '(1x, a, a)', &
         ERRMSG, 'Error while parsing damping values. Please check again.'
      call usage()
   end subroutine





   ! subroutine allocGlobFromLocSizes(loc, glob, n)
   !    real(bsa_real_t), allocatable, intent(in)  :: loc(:)
   !    real(bsa_real_t), allocatable, intent(out) :: glob(:, :)
   !    integer, intent(in) :: n
   !    integer :: dim

   !    if (allocated(loc)) then
   !       dim = size(loc)
   !       allocate(glob(dim, n))
   !    endif
   ! end subroutine


   ! subroutine moveAllocToGlob(glob, loc, idx)
   !    real(bsa_real_t), intent(out) :: glob(:, :)
   !    real(bsa_real_t), allocatable :: loc(:)
   !    integer, intent(in) :: idx

   !    if (allocated(loc)) then
   !       ! is a memcpy
   !       glob(:, idx) = loc
   !       deallocate(loc)
   !    endif
   ! end subroutine






   subroutine setup()
      !! This routine collects all the 
      !! calls to the BsaLibrary interface routines,
      !! before launch of the main computing process.


      ! NOTE: everything before bsa_Init() does not affect global (internal) vars.

      call bsa_setExportInCurrDir()  ! NOTE: this to not break plotter functioning..
      call bsa_closeUnitsAtEnd()
      call bsa_setExportFileFormat(i_exprt_form_)
      call bsa_setExportAppendMode(i_exprt_mode_)

      call bsa_setBRMExportDefaultMode(BSA_EXPORT_BRM_MODE_NONE)
      
      call bsa_setBfmMLR(.false.)
      call bsa_forceBsaClsExecution(.true.)

      call bsa_setValidateDeltasPolicy(BSA_VALIDATE_DELTAS_POLICY_NONE)

      call bsa_setMaxBkgPeakRestriction(.true.)

      call bsa_setPODTruncationThreshold(100.d0)

      call bsa_Init()  ! This initialises all necessary instances.

      ! SETTINGS
      if (i_suban /= 0) call bsa_setSubanType(i_suban)
      if (i_vers  /= 0) call bsa_setVersion(i_vers)
      if (i_defsc /= 0) call bsa_setScalingConv(i_defsc)
      call bsa_setSpectraComputation(ipsd=i_psd, ibisp=i_bisp)
      call bsa_setSpectraExtension(i_onlyd)
      call bsa_setSymmetries(i_bispsym, i_3dsym)
      call bsa_setTestMode(i_test)
      call bsa_setupClassic(i_nfreqs, r_df)
      call bsa_setupMesher(&
         i_svd, i_bkgrfmt, i_bkgaext, i_genpaext, i_maxaext, i_fcov, i_dumpmod)
      call bsa_setWindDirections(dirs(1 : i_ndirs), i_ndirs)
      call bsa_setWindTurbComps(tc(1 : i_ntc), i_ntc)

      ! NODAL
      call bsa_setTotalNOfNodes(i_nnodes)
      call bsa_setNodalNOfDOFs(i_nlibs)
      call bsa_setLoadedNodes(nodesl)
      call bsa_setLoadedNodalDOFs(libsl)
      call bsa_setNodalCoords(i_nnodes, nod_cords)

      ! WIND
      call bsa_setWindVertProf(i_varu)
      call bsa_setPSDType(i_su)
      call bsa_setWindAltDir(i_vert)
      call bsa_setWindZoneLimits(r_lims_z)
      call bsa_setAirDensity(r_aird)
      call bsa_setGlobalRotMatW2G(r_rotW2G)
      call bsa_setWZMeanWindVel(r_UBref_z)
      call bsa_setWZRefAlt(r_Zref_z)
      call bsa_setTurbWindScales(r_L_z)
      call bsa_setTurbWindSDT(r_std_z)
      call bsa_setWindCorrCoeffs(r_corrC_z)
      call bsa_setWindCorrExpnts(r_corrEx_z)
      call bsa_setWZRotMatW2G(r_rotW2G_z)
      call bsa_setIncidenceAngles(r_incang_z)
      call bsa_setNodalVel(r_UBnod)
      call bsa_setNodalWindZones(i_wzNod)
      call bsa_setNodalWindAltitudes(r_wAltNod)
      call bsa_setSpatialNodalCorr(r_corrNod)

      ! NOTE: in BSA, we want only LOADED NODES.
      r_wfc = r_wfc(:, :, nodesl)
      call bsa_setWindFCoeffs(r_wfc)
      
      ! MODAL
      call bsa_setModalInfo(i_ndofs, i_nm, r_modm, r_natf)
      call bsa_setModalMatrices(i_nm, r_Mg, r_Kg, r_Cg)

   end subroutine setup


   



   subroutine errAllocVarMsg_(varname, istat, emsg)
      character(len = *), intent(in) :: varname, emsg
      integer, intent(in) :: istat

      print '(/ 1x, 4a, i0)', &
         ERRMSG, 'Cannot allocate  "', varname, '".  Error code  ', istat
      print '(1x, 3a /)', &
         MSGCONT, 'Erorr message:  ', emsg(1 : len_trim(emsg))
      
      call releaseMemory(2)
   end subroutine


   subroutine errDeallocVarMsg_(varname, istat, emsg)
      character(len = *), intent(in) :: varname, emsg
      integer, intent(in) :: istat

      print '(/ 1x, 4a, i0)', &
         ERRMSG, 'Cannot de-allocate  "', varname, '".  Error code  ', istat
      print '(1x, 3a /)', &
         MSGCONT, 'Erorr message:  ', emsg(1 : len_trim(emsg))
      error stop
   end subroutine






   subroutine readDataFiles()

      call getFinData()
      call getBsaData()

      if (.not. (bsa_data_read_ .and. fin_data_read_)) then
         print '(1x, a, a)', &
            ERRMSG, 'Error reading input data from files.'
         call releaseMemory(5)
      endif

#ifdef __BSA_DEBUG
      print '(1x, a, a)', &
         INFOMSG, ' BSA   data read correctly.'
#endif
   end subroutine



   subroutine openFinelgInputFile()
      integer :: istat
      character(len = :), allocatable :: form_, access_

      if (l_formmode) then
         form_   = 'formatted'
         access_ = 'sequential'
      else
         form_   = 'unformatted'
         access_ = 'stream'
      endif

      open(unit=IUN_FINDATA    &
         , file=FINFILE_FNAME_ &
         , iostat=istat   &
         , form=form_     &
         , access=access_ &
         , action=IO_ACTION_READ)

      if (istat == 0) then
#ifdef _DEBUG
         print '(1x, a, a)', DBGMSG, 'Input file correctly opened.'
#endif
         rewind(IUN_FINDATA)
         return
      endif

      print *
      print '(/1x, a, a, i0)', &
         ERRMSG, 'Error opening input file. Error code  ', istat
      error stop
   end subroutine



   subroutine getFinData()
      integer :: istat, itmp
      character(len = 132) :: emsg
      
      call openFinelgInputFile()

      if (l_formmode) then
         read(IUN_FINDATA, *) i_nnodes
         read(IUN_FINDATA, *) i_nlibs
         read(IUN_FINDATA, *) i_nnodesl
         read(IUN_FINDATA, *) i_nlibsl
      else
         read(IUN_FINDATA) i_nnodes
         read(IUN_FINDATA) i_nlibs
         read(IUN_FINDATA) i_nnodesl
         read(IUN_FINDATA) i_nlibsl
      endif

      allocate(nodesl(i_nnodesl), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('nodesl', istat, emsg)

      allocate(libsl(i_nlibsl), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('libsl', istat, emsg)

      allocate(nod_cords(i_nnodes, 3), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('nod_cords', istat, emsg)

      if (l_formmode) then
         read(IUN_FINDATA, *) nodesl
         read(IUN_FINDATA, *) libsl
         read(IUN_FINDATA, *) nod_cords
      else
         read(IUN_FINDATA) nodesl
         read(IUN_FINDATA) libsl
         read(IUN_FINDATA) nod_cords
      endif
      nod_cords = transpose(nod_cords)



      if (l_formmode) then
         read(IUN_FINDATA, *) i_varu
         read(IUN_FINDATA, *) i_su
         read(IUN_FINDATA, *) i_vert
         read(IUN_FINDATA, *) i_degw
         read(IUN_FINDATA, *) r_aird
         read(IUN_FINDATA, *) r_rotW2G
         read(IUN_FINDATA, *) i_nzones
      else
         read(IUN_FINDATA) i_varu
         read(IUN_FINDATA) i_su
         read(IUN_FINDATA) i_vert
         read(IUN_FINDATA) i_degw
         read(IUN_FINDATA) r_aird
         read(IUN_FINDATA) r_rotW2G
         read(IUN_FINDATA) i_nzones
      endif
      
      allocate(r_Zref_z(i_nzones), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('r_Zref_z', istat, emsg)
      allocate(r_UBref_z(i_nzones), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('r_UBref_z', istat, emsg)
      allocate(r_alph_z(i_nzones), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('r_alph_z', istat, emsg)
      allocate(r_L_z(3, 3, i_nzones), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('r_L_z', istat, emsg)
      allocate(r_std_z(3, i_nzones), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('r_std_z', istat, emsg)
      allocate(r_corrC_z(3, 3, i_nzones), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('r_corrC_z', istat, emsg)
      allocate(r_corrEx_z(3, 3, i_nzones), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('r_corrEx_z', istat, emsg)
      allocate(r_lims_z(i_nzones+1), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('r_lims_z', istat, emsg)
      allocate(r_rotW2G_z(3, 3, i_nzones), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('r_rotW2G_z', istat, emsg)
      allocate(r_incang_z(i_nzones), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('r_incang_z', istat, emsg)

      allocate(i_wzNod(i_nnodes), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('i_wzNod', istat, emsg)
      allocate(r_wAltNod(i_nnodes), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('r_wAltNod', istat, emsg)
      allocate(r_UBnod(i_nnodes), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('r_UBnod', istat, emsg)
      itmp = i_nnodes * i_nnodes
      itmp = (itmp + i_nnodes) / 2
      allocate(r_corrNod(itmp, 3), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('r_corrNod', istat, emsg)

      allocate(r_wfc(i_nlibsl, i_degw+3, i_nnodes), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('r_corrNod', istat, emsg)

      if (l_formmode) then
         read(IUN_FINDATA, *) r_Zref_z
         read(IUN_FINDATA, *) r_UBref_z
         read(IUN_FINDATA, *) r_alph_z
         read(IUN_FINDATA, *) r_L_z
         read(IUN_FINDATA, *) r_std_z
         read(IUN_FINDATA, *) r_corrC_z
         read(IUN_FINDATA, *) r_corrEx_z
         read(IUN_FINDATA, *) r_lims_z
         read(IUN_FINDATA, *) r_rotW2G_z
         read(IUN_FINDATA, *) r_incang_z
         read(IUN_FINDATA, *) i_wzNod
         read(IUN_FINDATA, *) r_wAltNod
         read(IUN_FINDATA, *) r_UBnod
         read(IUN_FINDATA, *) r_corrNod

         read(IUN_FINDATA, *) r_wfc
      else
         read(IUN_FINDATA) r_Zref_z
         read(IUN_FINDATA) r_UBref_z
         read(IUN_FINDATA) r_alph_z
         read(IUN_FINDATA) r_L_z
         read(IUN_FINDATA) r_std_z
         read(IUN_FINDATA) r_corrC_z
         read(IUN_FINDATA) r_corrEx_z
         read(IUN_FINDATA) r_lims_z
         read(IUN_FINDATA) r_rotW2G_z
         read(IUN_FINDATA) r_incang_z
         read(IUN_FINDATA) i_wzNod
         read(IUN_FINDATA) r_wAltNod
         read(IUN_FINDATA) r_UBnod
         read(IUN_FINDATA) r_corrNod

         read(IUN_FINDATA) r_wfc
      endif



      if (l_formmode) then
         read(IUN_FINDATA, *) i_nm
         read(IUN_FINDATA, *) i_ndofs
      else
         read(IUN_FINDATA) i_nm
         read(IUN_FINDATA) i_ndofs
      endif
      allocate(r_natf(i_nm), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('r_natf', istat, emsg)
      allocate(r_modm(i_ndofs, i_nm), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('r_modm', istat, emsg)
      allocate(r_Mg(i_nm), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('r_Mg', istat, emsg)
      allocate(r_Kg(i_nm), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('r_Kg', istat, emsg)
      allocate(r_Cg(i_nm, i_nm), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('r_Cg', istat, emsg)
      allocate(r_xsist(i_nm), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('r_xsist', istat, emsg)
      allocate(r_xsiad(i_nm), stat=istat, errmsg=emsg)
      if (istat /= 0) call errAllocVarMsg_('r_xsiad', istat, emsg)
      if (l_formmode) then
         read(IUN_FINDATA, *) r_natf
         read(IUN_FINDATA, *) r_modm
         read(IUN_FINDATA, *) r_Mg
         read(IUN_FINDATA, *) r_Kg
         read(IUN_FINDATA, *) r_Cg
         read(IUN_FINDATA, *) r_xsist
         read(IUN_FINDATA, *) r_xsiad
      else
         read(IUN_FINDATA) r_natf
         read(IUN_FINDATA) r_modm
         read(IUN_FINDATA) r_Mg
         read(IUN_FINDATA) r_Kg
         read(IUN_FINDATA) r_Cg
         read(IUN_FINDATA) r_xsist
         read(IUN_FINDATA) r_xsiad
      endif

      fin_data_read_ = .true.
      close(IUN_FINDATA)
#ifdef __BSA_DEBUG
      print '(1x, a, a)', &
         INFOMSG, 'FINELG data read correctly.'
#endif
   end subroutine





   subroutine getBsaData()
      character(len = 256) :: label
      character(len = *), parameter :: fmt_a = '(a)', fmt_i = '(i)', fmt_f = '(f)'
      integer :: i


      if (.not. fin_data_read_) return

      open(unit=IUN_BSADATA   &
         , file='bsa.bsadata' &
         , form='formatted'   &
         , action=IO_ACTION_READ)

      read(IUN_BSADATA, fmt_a) label
      read(IUN_BSADATA, fmt_i) i_suban
      read(IUN_BSADATA, fmt_i) i_vers
      read(IUN_BSADATA, fmt_i) i_defsc
      read(IUN_BSADATA, fmt_i) i_psd
      read(IUN_BSADATA, fmt_i) i_bisp
      read(IUN_BSADATA, fmt_i) i_onlyd
      read(IUN_BSADATA, fmt_i) i_bispsym
      read(IUN_BSADATA, fmt_i) i_3dsym
      read(IUN_BSADATA, fmt_i) i_test

      read(IUN_BSADATA, fmt_a) label
      read(IUN_BSADATA, fmt_i) i_nfreqs
      read(IUN_BSADATA, fmt_f) r_df

      read(IUN_BSADATA, fmt_a) label
      read(IUN_BSADATA, fmt_i) i_svd
      read(IUN_BSADATA, fmt_i) i_bkgrfmt
      read(IUN_BSADATA, fmt_i) i_bkgaext
      read(IUN_BSADATA, fmt_i) i_genpaext
      read(IUN_BSADATA, fmt_i) i_maxaext
      read(IUN_BSADATA, fmt_i) i_fcov
      read(IUN_BSADATA, fmt_i) i_dumpmod

      ! directions
      read(IUN_BSADATA, fmt_a)   label
      read(IUN_BSADATA, fmt_i)   i_ndirs
      do i = 1, i_ndirs
         read(IUN_BSADATA, fmt_i) dirs(i)
      enddo

      ! turbulence
      read(IUN_BSADATA, fmt_a)   label
      read(IUN_BSADATA, fmt_i)   i_ntc
      do i = 1, i_ntc
         read(IUN_BSADATA, fmt_i) tc(i)
      enddo

      ! ! nodes loaded
      ! read(IUN_BSADATA, fmt_a)   label
      ! read(IUN_BSADATA, fmt_a)   label
      ! call getLoadedNodesFromString(label)

      bsa_data_read_ = .true.
   end subroutine



   ! subroutine getLoadedNodesFromString(label)
   !    !! BUG: for the moment, supports only 1 line (1 range)
   !    character(len = *), intent(in) :: label
   !    character(len = *), parameter  :: col = ':'
   !    integer :: ilen, ibl = 1, i, icount = 0, iini = 1
   !    integer :: vals(3), istat, ival
   !    character(len = 132) :: emsg

   !    do while (label(ibl:ibl) == ' ')
   !       ibl = ibl + 1
   !    enddo

   !    ilen = len_trim(label)
   !    i    = ibl

   !    if ( label(i : ilen) == 'all' ) then
         
   !       i_nnodesl = i_nnodes
   !       nodesl    = [1 : i_nnodesl]
   !       goto 100
   !    end if

   !    do while (i <= ilen)
   !       if (label(i:i) == col) then ! read left-side value
   !          icount = icount + 1
   !          read(label(iini : i-1), fmt='(i)') vals(icount)
   !          iini = i + 1
   !       endif
   !       i = i + 1
   !    enddo
   !    ! treat last value!
   !    icount = icount + 1
   !    read(label(iini : ilen), fmt='(i)') vals(icount)


   !    if (icount == 1) then ! only one node loaded

   !       i_nnodesl = 1
   !       if (vals(1) > i_nnodes) vals(1) = i_nnodes
   !       nodesl = vals(1:1)

   !    elseif (icount == 2) then ! linspace

   !       if (vals(1) > i_nnodes) vals(1) = i_nnodes
   !       if (vals(2) > i_nnodes) vals(2) = i_nnodes

   !       i_nnodesl = vals(2) - vals(1) + 1
   !       allocate(nodesl(i_nnodesl), stat=istat, errmsg=emsg)
   !       if (istat /= 0) call errAllocVarMsg_('nodesl', istat, emsg)
   !       ival = vals(1) - 1
   !       do i = 1, i_nnodesl
   !          nodesl(i) = ival + i
   !       enddo

   !    else ! ==3, range

   !       ! TODO: implement
   !       error stop ERRMSG // ' IMPLEMENT ICOUNT=3'
   !    endif

   !    100 print '(1x, a, a)', &
   !       INFOMSG, 'List of loaded nodes'
   !    print '( 10( "  ", i6) )', &
   !       nodesl
   ! end subroutine







   subroutine modalRecombination(r_Phi, m2mr, m3mr, m2o2mr)
      real(bsa_real_t), intent(in) :: r_Phi(:, :), m2mr(:), m3mr(:)
      real(bsa_real_t), intent(in), allocatable :: m2o2mr(:)
      integer :: ndofs, nmodes, nm3
      integer :: idof, imode
      integer :: imodeM1, posm2, posm3, itmp1
      logical :: is_diag
      real(bsa_real_t) :: phi1, r_tmp


      ndofs  = size(r_Phi, 1)
      nmodes = size(r_Phi, 2)
      nm3    = size(m3mr)

      is_diag = nm3 == nmodes

      allocate(m2_r_diag(ndofs))
      m2_r_diag = 0
      allocate(m3_r_diag(ndofs))
      m3_r_diag = 0
      allocate(m2o2_r_diag(ndofs))
      m2o2_r_diag = 0

      ! UNILATERAL (do it in any case)
      do imode = 1, nmodes

         imodeM1 = imode - 1

         if (is_diag) then
            posm2 = imode
            posm3 = imode
         else
            itmp1 = imodeM1 * nmodes
            posm2 = itmp1 + imode
            posm3 = itmp1*nmodes + posm2
         endif

         ! TODO: use do concurrent
         do idof = 1, ndofs

            phi1  = r_Phi(idof, imode)
            r_tmp = phi1 * phi1


            m2_r_diag(idof) = m2_r_diag(idof) + m2mr(posm2) * r_tmp
            if (allocated(m2o2mr)) m2o2_r_diag(idof) = m2o2_r_diag(idof) + m2o2mr(posm2) * r_tmp

            r_tmp = r_tmp * phi1
            m3_r_diag(idof) = m3_r_diag(idof) + m3mr(posm3) * r_tmp
         enddo ! idof
      enddo ! imode

      sk_r_diag = m3_r_diag / (m2_r_diag)**(CST_3d2)


      if (is_diag) return


      ! CCC (Complete Cubic Combination)
      allocate(m2_r_full(ndofs))
      m2_r_full = 0
      allocate(m3_r_full(ndofs))
      m3_r_full = 0
      allocate(m2o2_r_full(ndofs))
      m2o2_r_full = 0

      block
         integer :: jmode, kmode, jmodeM1, kmodeM1
         integer :: itmp12, itmp13, itmp23
         real(bsa_real_t), dimension(ndofs) :: phi2_, phi3_

         do kmode = 1, nmodes

            kmodeM1 = kmode - 1
            phi3_   = r_Phi(:, kmode)

            itmp12 = kmodeM1 * nmodes
            itmp13 = itmp12  * nmodes
            
            do jmode = 1, nmodes

               jmodeM1 = jmode - 1
               phi2_   = r_Phi(:, jmode)
               itmp23  = jmodeM1 * nmodes
               posm2   = itmp12 + jmode

               ! 2nd order
               ! TODO: use do concurrent
               do idof = 1, ndofs

                  m2_r_full(idof) = m2_r_full(idof) + &
                     m2mr(posm2) * phi2_(idof) * phi3_(idof)

                  if (allocated(m2o2mr)) &
                     m2o2_r_full(idof) = m2o2_r_full(idof) + &
                        m2o2mr(posm2) * phi2_(idof) * phi3_(idof)
               enddo ! idof


               do imode = 1, nmodes

                  posm3 = itmp13 + itmp23 + imode

                  ! TODO: use do concurrent
                  do idof = 1, ndofs

                     phi1 = r_Phi(idof, imode)

                     m3_r_full(idof) = m3_r_full(idof) + &
                        m3mr(posm3) * phi1 * phi2_(idof) * phi3_(idof)
                  enddo ! idof
               enddo ! imode
            enddo ! jmode
         enddo ! kmode

         sk_r_full = m3_r_full / (m2_r_full)**(CST_3d2)

      end block
   end subroutine





   subroutine releaseMemory(iexit)
      integer, intent(in) :: iexit
      integer :: istat
      character(len = 132) :: emsg
      logical :: lflag

      if (.not. bsa_isCleaned()) call bsa_Finalise()

      inquire(unit=IUN_FINDATA, opened=lflag)
      if (lflag) close(IUN_FINDATA)
      inquire(unit=IUN_BSADATA, opened=lflag)
      if (lflag) close(IUN_BSADATA)

      istat = 0

      if (allocated(nodesl)) deallocate(nodesl, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('nodesl', istat, emsg)

      if (allocated(libsl)) deallocate(libsl, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('libsl', istat, emsg)

      if (allocated(nod_cords)) deallocate(nod_cords, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('nod_cords', istat, emsg)


      if (allocated(r_Zref_z)) deallocate(r_Zref_z, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('r_Zref_z', istat, emsg)
      if (allocated(r_UBref_z)) deallocate(r_UBref_z, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('r_UBref_z', istat, emsg)
      if (allocated(r_alph_z)) deallocate(r_alph_z, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('r_alph_z', istat, emsg)
      if (allocated(r_L_z)) deallocate(r_L_z, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('r_L_z', istat, emsg)
      if (allocated(r_std_z)) deallocate(r_std_z, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('r_std_z', istat, emsg)
      if (allocated(r_corrC_z)) deallocate(r_corrC_z, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('r_corrC_z', istat, emsg)
      if (allocated(r_corrEx_z)) deallocate(r_corrEx_z, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('r_corrEx_z', istat, emsg)
      if (allocated(r_lims_z)) deallocate(r_lims_z, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('r_lims_z', istat, emsg)
      if (allocated(r_rotW2G_z)) deallocate(r_rotW2G_z, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('r_rotW2G_z', istat, emsg)
      if (allocated(r_incang_z)) deallocate(r_incang_z, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('r_incang_z', istat, emsg)

      if (allocated(i_wzNod)) deallocate(i_wzNod, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('i_wzNod', istat, emsg)
      if (allocated(r_wAltNod)) deallocate(r_wAltNod, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('r_wAltNod', istat, emsg)
      if (allocated(r_UBnod)) deallocate(r_UBnod, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('r_UBnod', istat, emsg)
      if (allocated(r_corrNod)) deallocate(r_corrNod, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('r_corrNod', istat, emsg)

      if (allocated(r_wfc)) deallocate(r_wfc, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('r_corrNod', istat, emsg)  


      if (allocated(r_natf)) deallocate(r_natf, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('r_natf', istat, emsg)
      if (allocated(r_modm)) deallocate(r_modm, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('r_modm', istat, emsg)
      if (allocated(r_Mg)) deallocate(r_Mg, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('r_Mg', istat, emsg)
      if (allocated(r_Kg)) deallocate(r_Kg, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('r_Kg', istat, emsg)
      if (allocated(r_Cg)) deallocate(r_Cg, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('r_Cg', istat, emsg)
      if (allocated(r_xsist)) deallocate(r_xsist, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('r_xsist', istat, emsg)
      if (allocated(r_xsiad)) deallocate(r_xsiad, stat=istat, errmsg=emsg)
      if (istat /= 0) call errDeallocVarMsg_('r_xsiad', istat, emsg)

      if (iexit == 0) then
         print '(/ 1x, a, a)', INFOMSG, 'BSA terminated correctly.'
      else
         print '(/ 1x, a, a, i0)', &
            ERRMSG, 'BSA terminated with error. Exit status code  ', iexit
      endif
      stop
   end subroutine releaseMemory


end program