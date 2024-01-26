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
submodule(BsaLib) BsaLib_ClassicImpl

   use BsaLib_Data
#ifdef _BSA_DEBUG
   use BsaLib_IO,   only: unit_debug_
#endif
   implicit none (type, external)

contains


   !> BUG: now only supports EVENLY SPACED FREQUENCIES..
   module subroutine mainClassic_(m2mf_cls, m2mr_cls, m2o2mr_cls, m3mf_cls, m3mr_cls)
      use BsaLib_Functions
      real(bsa_real_t), allocatable :: m2mf_cls(:), m2mr_cls(:), m2o2mr_cls(:), m3mf_cls(:), m3mr_cls(:)
      integer :: idim2
      real(bsa_real_t), allocatable :: S_uvw(:, :), f(:)
      integer :: itc_, idir_, idxi, idxe

#ifdef _BSA_DEBUG
      write(unit_debug_, *) INFOMSG, '@BsaClassicImpl::mainClassic_() : Init BSA-Classic main...'
#endif

      call computeFreqsVect_(settings, struct_data, f)

      ! NOTE: reset in case NFREQs has been changed in previous call
      NFREQS = settings%nfreqs_

      ! some shared memory allocation
      if (settings%i_compute_psd_ == 1) then
         allocate(m2mf_cls(dimM_psd_))
         m2mf_cls = 0._bsa_real_t

         allocate(m2mr_cls(dimM_psd_))
         m2mr_cls = 0._bsa_real_t

         allocate(m2o2mr_cls(dimM_psd_))
         m2o2mr_cls = 0._bsa_real_t
      endif
      if (settings%i_compute_bisp_== 1) then
         allocate(m3mf_cls(dimM_bisp_))
         m3mf_cls = 0._bsa_real_t

         allocate(m3mr_cls(dimM_bisp_))
         m3mr_cls = 0._bsa_real_t
      endif


#ifdef _BSA_DEBUG
      print '(1x, a, a, i0)', INFOMSG, 'n. of frequencies to be computed=', settings%nfreqs_
      print '(1x, a, a, i0)', INFOMSG, 'PSD  modal extension=', dimM_psd_
      print '(1x, a, a, i0)', INFOMSG, 'BISP modal extension=', dimM_bisp_
#endif


#ifdef _BSA_USE_GPU
      block
         ! real(bsa_real_t), allocatable :: S_uvw_T_(:, :)

         call bsacl_AcquireResultBFMVect(m3mf_cls)
         call bsacl_AcquireComputationFreqs(NFREQS, f, NFREQS, f)
         ! S_uvw_T_ = transpose(S_uvw)  ! n_dim_ x n_freqs
         ! deallocate(S_uvw)
         ! call bsacl_AcquireBaseWindTurbPSD(S_uvw_T_)
         call bsacl_Run(ierr_cl_)
         if (ierr_cl_ == BSACL_PROBLEM_DIMENSIONS_TOO_SMALL) then
            print '(1x, 2a)', WARNMSG, &
               'Problem dimensions are too small for GPU parallelisation. Using CPU.'
            ! S_uvw = transpose(S_uvw_T_)
            ! deallocate(S_uvw_T_)
            goto 997  ! CPU computation
         endif
      end block
      goto 998
#endif

      997 continue
#ifdef _BSA_DEBUG
      write(unit_debug_, *) INFOMSG, '@BsaClassicImpl::mainClassic_() : computing nodal wind turbulence PSDs...'
#endif

      idim2 = struct_data%nn_load_ * wd%i_ntc_  * wd%i_ndirs_
      allocate(S_uvw(settings%nfreqs_, idim2))
      idxi = 1
      idxe = struct_data%nn_load_
      do itc_ = 1, wd%i_ntc_

         do idir_ = 1, wd%i_ndirs_

            ! BUG:  difference between idir and itc ???
            ! NOTE: done for all the loaded nodes at once !
            S_uvw(:, idxi:idxe) = wd%evalPSD(settings%nfreqs_, f, &
               struct_data%nn_load_, struct_data%n_load_, wd%dirs_(idir_), wd%tc_(itc_))

            idxi = idxe + 1
            idxe = idxe + struct_data%nn_load_
         enddo ! i direction
      enddo ! i turb comp


#ifdef _BSA_DEBUG
      write(unit_debug_, *) INFOMSG, '@BsaClassicImpl::mainClassic_() : computing nodal wind turbulence PSDs -- ok.'

      call bsa_exportPSDToFile('psd_Suvw.txt', S_uvw, f)
#endif

      block
         real(bsa_real_t) :: m2(idim2)

         call intgSpectraVect_(settings%nfreqs_, f, psd=S_uvw, m2=m2)
         call bsa_exportMomentToFile('m2_PSDs.txt', m2)
      end block


      call checkMaxAllocation_()

      ! BUG: set it from .dat file
      settings%i_scalar_vers_ = 0  ! 0=vectorised, 1=scalar


      if (settings%i_scalar_vers_ == 0) then ! VECTORISED

         print '(/1x, 2a)', INFOMSG, 'Using  VECTORISED  version'

         block
            real(bsa_real_t), allocatable :: psd(:, :), bisp(:, :, :)

            !===========================================================
            ! MODAL FORCES
            !
            call getBFM_vect_cls(f, S_uvw, psd, bisp)
            if (allocated(S_uvw)) deallocate(S_uvw)
            call intgSpectraVect_(settings%nfreqs_, f, psd=psd, m2=m2mf_cls, bisp=bisp, m3=m3mf_cls)
            if (settings%i_compute_psd_ == 1) call bsa_exportPSDToFile('psdmf.txt', psd, f)


            !===========================================================
            ! MODAL RESPONSES
            !
            call getBRM_vect_cls(f, psd, bisp)
            call intgSpectraVect_(settings%nfreqs_, f, psd=psd, m2=m2mr_cls, bisp=bisp, m3=m3mr_cls)
            if (settings%i_compute_psd_ == 1) call bsa_exportPSDToFile('psdmr.txt', psd, f)

            block
               real(bsa_real_t), allocatable :: omegas(:)
               integer :: i

               omegas = f * CST_PIt2
               do concurrent (i = 1 : dimM_psd_) shared(omegas, psd)
                  psd(:, i) = psd(:, i) * omegas(:) * omegas(:)
               enddo
               call intgSpectraVect_(settings%nfreqs_, f, psd=psd, m2=m2o2mr_cls)
            end block

#if 0
            block
               real(bsa_real_t), allocatable :: psd_r(:, :)
               integer(int32) :: im, jm, idx

               allocate(psd_r(settings%nfreqs_, dimNr_psd_))
               psd_r = 0._bsa_real_t

               idx = 0
               do jm = 1, struct_data%modal_%nm_eff_
                  do im = 1, struct_data%modal_%nm_eff_

                     idx = idx + 1
                     psd_r = psd_r + &
                        matmul(psd(:, idx:idx), reshape(struct_data%modal_%phi_(:, im:im)*struct_data%modal_%phi_(:, jm:jm), [1, struct_data%ndofs_]))
                  enddo
               enddo

               open(newunit=idx, file='psd_r', action='write', status='replace', form='unformatted')
               rewind(idx)
               write(idx) settings%nfreqs_, struct_data%ndofs_
               do im = 1, settings%nfreqs_
                  write(idx) f(im), psd_r(im, :)
               enddo
               close(idx)

               open(newunit=idx, file='psd_r.txt', action='write', status='replace', form='formatted')
               rewind(idx)
               write(idx, '(*(i))') settings%nfreqs_, struct_data%ndofs_
               do im = 1, settings%nfreqs_
                  write(idx, '(*(g, 2x))') f(im), psd_r(im, :)
               enddo
               close(idx)
            end block
#endif

            if (allocated(psd))  deallocate(psd)
            if (allocated(bisp)) deallocate(bisp)
         end block



            !===============
      else  ! SCALAR VERSION
            !===============


         print '(/1x, 2a)', INFOMSG, 'Using    SCALAR    version'

         print '(/ 1x, a, a /)', &
            WARNMSG, 'For scalar version, computation of m2o2_mr not yet implemented !'

         block
            real(bsa_real_t) :: fi, fj, dw, dw2, omg
            real(bsa_real_t), allocatable :: S_uvw_pad(:, :)

            integer, pointer :: jfr_ext => null()
            integer, target  :: one_ext = 1

            integer :: lpad, indxi, indxe

            real(bsa_real_t), dimension(dimM_psd_)  :: psdfm, psdrm, r_tmp
            real(bsa_real_t), dimension(dimM_bisp_) :: bispfm, bisprm

            psdfm  = 0._bsa_real_t
            psdrm  = 0._bsa_real_t
            bispfm = 0._bsa_real_t
            bisprm = 0._bsa_real_t

            dw  = (f(2) - f(1)) * CST_PIt2 ! [rad/s]
            dw2 = dw*dw
            ! get padded length
            lpad  = (settings%nfreqs_ - 1) / 2
            indxi = lpad + 1
            indxe = lpad + settings%nfreqs_
            lpad  = 2*lpad + settings%nfreqs_

            allocate(S_uvw_pad(idim2, lpad))
            S_uvw_pad(:, indxi : indxe) = transpose(S_uvw)


            if (settings%i_compute_bisp_ == 0) then
               jfr_ext => one_ext
            else
               jfr_ext => settings%nfreqs_
            endif


            do ifr = 1, settings%nfreqs_

               fi  = f(ifr)
               omg = fi * CST_PIt2

               do jfr = 1, jfr_ext

                  fj = f(jfr)

                  call getBFM_scalar_cls(ifr, jfr, fi, fj, S_uvw, S_uvw_pad(:, ifr - 1 + jfr), psdfm, bispfm)

                  ! NOTE: using same infl area for each point in space!
                  m3mf_cls = m3mf_cls + bispfm * dw2


                  call getBRM_scalar_cls(ifr, jfr, fi, fj, psdfm, psdrm, bispfm, bisprm)
                  m3mr_cls = m3mr_cls + bisprm * dw2
               enddo ! i freqs

               if (settings%i_compute_psd_ == 1) then
                  m2mf_cls   = m2mf_cls + psdfm * dw
                  r_tmp      = psdrm * dw
                  m2mr_cls   = m2mr_cls + r_tmp
                  m2o2mr_cls = m2o2mr_cls + r_tmp * omg*omg
               endif

               print '(1x, a, 2(i12, a))', &
                  INFOMSG, ifr*settings%nfreqs_, '  out of  ', settings%nfreqs_**2, '  done...'

            enddo ! j freqs
         end block

      endif ! vect/scalar versions


      998 continue
      if (allocated(f)) deallocate(f)

#ifdef _BSA_DEBUG
      write(unit_debug_, *) &
         INFOMSG, '@BsaClassicImpl::mainClassic_() : Init BSA-Classic main -- ok.'
#endif
   end subroutine mainClassic_







   subroutine checkMaxAllocation_()
      integer(int32) :: itmp

      if (settings%i_test_mode_ == 0) then

         ! Computing max allocation size if it was VECTORISED.
         if (settings%i_compute_bisp_ == 1) then
            itmp = settings%nfreqs_**2  *  dimM_bisp_
         else
            itmp = settings%nfreqs_  *  dimM_bisp_
         endif

         if (settings%i_scalar_vers_ == 0) then

            if (itmp > MAX_VECT_ALLOC_ELEMS) then

               print '( /, 1x, a, a, i0, ")" )', WARNMSG, 'Too high allocation size for VECTORISED BSA version (', itmp
               print '( 1x, a, a, / )', MSGCONT, 'Switching to SCALAR version.'

               settings%i_scalar_vers_ = 1
            endif

         else ! user wants SCALAR.

            ! NOTE: if just PSDs, we can ALWAYS go for vectorised
            if (settings%i_compute_bisp_ == 0) then

               ! Still, check, better.
               if (settings%nfreqs_ * dimM_psd_ < MAX_VECT_ALLOC_ELEMS) then

                  print '( /, 1x, a, a )', NOTEMSG, 'Requested SCALAR BSA version, but for only PSDs computation.' 
                  print '(1x, a, a)', MSGCONT, 'Switching to VECTORISED for perf.'
                  settings%i_scalar_vers_ = 0
               endif

            else ! BISPs as well -> just warn, do not force changing..

               if (itmp < MAX_VECT_ALLOC_ELEMS) then

                  print '( /, 1x, a, a )', NOTEMSG, 'Running SCALAR BSA version, but VECTORISED (preferable) is possible.'
                  print '(1x, a, a)', MSGCONT, 'Consider changing setting.'
               endif
            endif
         endif

      else ! testing mode (==1, yes), keep things as such.

         print '(/1x, a, a)', &
            WARNMSG, 'Frequency definition not being checked for optimal values !'
      endif
   end subroutine








   subroutine computeFreqsVect_(setts, struct, f)
      ! class(bsa_classic_t), intent(inout)   :: this
      class(settings_t), intent(inout)      :: setts
      class(StructureData_t), intent(inout) :: struct
      real(bsa_real_t), allocatable, intent(out) :: f(:)

      logical :: l_df_big = .false.
      integer(int32)   :: nfreqs_1
      real(bsa_real_t) :: df_ref, max_freq, max_freq_ref

      if (setts%nfreqs_ == 0 .or. setts%df_ == 0._bsa_real_t) &
         call bsa_Abort('Either NFREQs or DF are == 0.')


! #ifdef _BSA_DEBUG
!       write(unit_debug_, *) INFOMSG, '@BsaClassicImpl::computeFreqsVect_() : computing frequencies...'
! #endif

      ! NOTE: make sure resonant peak are computed
      call struct%ComputeResPeakWidths()

      ! BUG: let the user choose the dividend
      df_ref = minval(struct%res_peak_width_) / 5

      if (setts%df_ > df_ref) then

         l_df_big = .true.
         print '( /, 1x, a, 2(a, f12.5))', WARNMSG, 'specified  df=', setts%df_, &
            '  is bigger than suggested=', df_ref

      ! BUG: also here, let user choose limit
      elseif (setts%df_ < df_ref / 10) then

         print '( /, 1x, a, 2(a, f12.5))', WARNMSG, 'specified  df=', setts%df_, &
            '  is smaller than  1/10th  of suggested=', df_ref
         print '(1x, a, a /)', MSGCONT, 'Consider increasing it.'
      endif

      nfreqs_1 = setts%nfreqs_ - 1  ! NOTE: do not consider 0 freq

      if (settings%i_test_mode_ == 0) then ! Do actual check only if NO TESTING MODE.

         max_freq     = setts%df_ * nfreqs_1
         max_freq_ref = maxval(struct%modal_%nat_freqs_)

         if (max_freq < max_freq_ref) then

            if (l_df_big) then ! try with suggested one

               max_freq = df_ref * nfreqs_1
               if (max_freq < max_freq_ref) then ! find for the right nfreq value (using df ref)

                  nfreqs_1     = ceiling(max_freq_ref / df_ref)
                  max_freq     = max_freq_ref

                  print '(/ 1x, a, a, i0, a )', &
                     WARNMSG, '"nfreq=', setts%nfreqs_, '"  is too small to reach max frequency (even with suggested "df").'
                  print '( 1x, a, a, i0 /)', &
                     MSGCONT, 'To avoid errors in estimation, it is going to be considered nfreqs= ', nfreqs_1+1
               else

                  print '( /, 1x, a, a, i0, a / )', &
                     WARNMSG, 'with specified "nfreq', setts%nfreqs_, &
                        '", full frequency range coverage is ensured using suggested "df".'
               endif

               setts%df_ = df_ref  ! override df with suggested value


            else ! chosen "df" is OK.


               nfreqs_1 = ceiling(max_freq_ref / setts%df_)
               print '( /, 1x, a, a, i0, a /, 20x, a, i5, / )', &
                  WARNMSG, '"nfreq=', setts%nfreqs_ ,'" is too small to reach max frequency.'

               print '(1x, a, a, i0 /)', &
                  MSGCONT, 'To avoid errors in estimation, it is going to be considered   nfreqs=', nfreqs_1+1

               max_freq = max_freq_ref
            endif

         else

            if (l_df_big) then

               print '( /, 1x, a, a, f12.5, " > ", f12.5, ")")', &
                  WARNMSG, &
                  'chosen "df" is greater than suggested one (', setts%df_, df_ref
            endif
         endif
      endif  ! test mode == 0


      ! Actual freqs computation
      if (setts%i_def_scaling_ == 1) then

         ! NOTE: automatically odd (because of the +1)
         setts%nfreqs_ = (nfreqs_1 * 2) + 1
         if (mod(setts%nfreqs_, 2) == 0) call bsa_Abort('Needing odd n. of frequencies.')
         allocate(f(setts%nfreqs_))
         f = [-nfreqs_1 : nfreqs_1] * setts%df_

      else ! ==2 (frequencies conventions)  [WARNING]

         if (mod(nfreqs_1, 2) /= 0) nfreqs_1 = nfreqs_1 + 1 ! make nfreqs-1 even

         setts%nfreqs_ = nfreqs_1 + 1 ! NOTE: +1 because we consider 0 as well.
         allocate(f(setts%nfreqs_))
         f = [0 : nfreqs_1] * setts%df_

      endif

#ifdef _BSA_DEBUG
      write(unit_debug_, *) INFOMSG, '@BsaClassicImpl::computeFreqsVect_() : computing frequencies -- ok.' 
#endif
   end subroutine computeFreqsVect_






   subroutine intgSpectraVect_(nf, f, psd, m2, bisp, m3)
      integer(bsa_int_t), intent(in) :: nf
      real(bsa_real_t), intent(in)   :: f(nf)
      ! integer(int32), intent(in), optional   :: dimpsd, dimbisp
      real(bsa_real_t), intent(in), optional  :: psd(nf, *), bisp(nf, nf, *)
      real(bsa_real_t), intent(out), optional :: m2(:), m3(:)

      integer(int32)   :: nf_1 = 0, dim = 0, i
      real(bsa_real_t) :: delta
      real(bsa_real_t) :: rtmp, d_2, d2, d2_2


      delta = f(2) - f(1)
      if (settings%i_def_scaling_ == 1) delta = delta * CST_PIt2  ! [rad/s]
      d_2  = delta / 2._bsa_real_t

      !  PSDs
      if (present(psd) .and. present(m2)) then

         ! full integration
         dim = size(m2)
         m2  = sum(psd(:, 1:dim), dim=1) * delta

         ! removing excess from vertexes
         m2(:) = m2(:) - ((psd(1, 1:dim) + psd(nf, 1:dim)) * d_2)
      endif


      ! BISPs
      if (present(bisp) .and. present(m3)) then

         d2   = delta * delta

         ! full integration
         dim = size(m3)
         m3  = sum(sum(bisp(:, :, 1:dim), 1), 1) * d2


         ! removing excess from vertexes/borders

         rtmp = CST_3d2 * d2
         d2_2 = d2 / 2._bsa_real_t
         nf_1 = nf - 1

         ! LEFT
         ! vertex
         m3(:) = m3(:) - (bisp(1, 1, 1:dim) * rtmp)
         ! side
         m3(:) = m3(:) - sum(bisp(2:nf_1, 1, 1:dim) * d2_2, 1)
         ! vertex
         m3(:) = m3(:) - (bisp(nf, 1, 1:dim) * rtmp)

         ! sides (up/down)
         do i = 2, nf_1
            m3(:) = m3(:) - bisp(1, i, 1:dim) * d2_2
            m3(:) = m3(:) - bisp(nf, i, 1:dim) * d2_2
         enddo

         ! RIGHT
         ! vertex
         m3(:) = m3(:) - (bisp(1, nf, 1:dim) * rtmp)
         ! side
         m3(:) = m3(:) - sum(bisp(2:nf_1, nf, 1:dim) * d2_2, 1)
         ! vertex
         m3(:) = m3(:) - (bisp(nf, nf, 1:dim) * rtmp)
      endif
   end subroutine intgSpectraVect_



end submodule
