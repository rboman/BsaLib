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
submodule(BsaLib_Functions) BsaLib_FunctionsImpl

#include "../../precisions"

   use BsaLib_CONSTANTS
   use Logging
   use BsaLib_Utility
   use BsaLib_Data, only: bsa_Abort, do_trunc_POD_, POD_trunc_lim_
   use BsaLib_IO, only: INFOMSG, WARNMSG, ERRMSG, MSGCONT, DBGMSG, NOTEMSG &
                        , unit_dump_bfm_, unit_debug_, undebug_fname_
   implicit none

   
contains


   module subroutine setBsaFunctionLocalVars()

      NFREQS  = settings%nfreqs_
      
      ! nodal
      NNODES  = struct_data%nn_
      NNODESL = struct_data%nn_load_
      NLIBS   = struct_data%nlibs_    ! tot n. of LIBs per node
      NLIBSL  = struct_data%nlibs_load_    ! actual n. of loaded LIBs
      
      ! modal
      NMODES     = struct_data%modal_%nm_
      NMODES_EFF = struct_data%modal_%nm_eff_
      MODES      = struct_data%modal_%modes_
      
      ! wind
      NTCOMPS = wd%i_ntc_
      TCOMPS  = wd%tc_
      NDIRS   = wd%i_ndirs_
      DIRS    = wd%dirs_
      NPSDEL  = NNODESL * NTCOMPS * NDIRS
   end subroutine




   module function getFM_full_tnm_scalar_msh_(fi, fj) result(bfm)
      real(RDP), intent(in) :: fi, fj
      real(RDP) :: bfm(dimM_bisp_)

      real(RDP) :: fiPfj(1), abs_fi, abs_fj, abs_fiPfj

      ! indexes
      integer(kind = 4) :: itc, tc, tcP3, iposM
      integer(kind = 4) :: iposNK, iposNJ, iposNI
      integer(kind = 4) :: ink, inj, ini
      integer(kind = 4) :: nk, nj, ni
      integer(kind = 4) :: posKi, posJi, posIi
      ! integer(kind = 4) :: posKe, posJe, posIe
      integer(kind = 4) :: imk, imj, imi
      integer(kind = 4) :: ilk

      ! modal matrix slices
      real(RDP) :: phiJ(1, NLIBSL), phiI(NLIBSL, 1) !, phiK(1, 1, NLIBSL)
      real(RDP) :: phiK_(NMODES_EFF), phiK

      ! wind forces coeffs
      real(RDP), dimension(1, NLIBSL)    :: ajU, aj
      real(RDP), dimension(NLIBSL, 1)    :: aiU, ai, akU, ak

      ! basic PSDs
      real(RDP), dimension(1, NNODESL) :: S_IJK_fi, S_IJK_fj, S_IJK_fiPfj
      real(RDP) :: S_K_fi, S_K_fj, S_K_fiPfj
      real(RDP) :: S_J_fi, S_J_fj, S_J_fiPfj
      real(RDP) :: S_I_fi, S_I_fj, S_I_fiPfj

      ! nodal spactial correlations
      real(RDP) :: corrIJ, corrIK, corrJK

      ! crossed PSDs
      real(RDP) :: S_JK_fi, S_JK_fj
      real(RDP) :: term1, term2, term3, BF_IJK_ijk(NLIBSL, NLIBSL)
      real(RDP), dimension(NLIBSL, NLIBSL) :: tmp1, tmp2, tmp3

      ! frequencies values
      fiPfj     = fi + fj
      abs_fi    = abs(fi)
      abs_fj    = abs(fj)
      abs_fiPfj = abs(fiPfj(1))


      ! NOTE: preinitalise to 0, to avoid uninitialised precision errors
		!       Like using memset() in C.
      bfm = 0._RDP


      do itc = 1, NTCOMPS

         tc   = wd%tc_(itc)
         tcP3 = tc + 3

         ! prefetch wind turbulence PSDs for all loaded nodes
         S_IJK_fi = wd%evalPSD(1, [fi], NNODESL, struct_data%n_load_, 1, tc)
         
         S_IJK_fj = wd%evalPSD(1, [fj], NNODESL, struct_data%n_load_, 1, tc)
         
         S_IJK_fiPfj = wd%evalPSD(1, fiPfj, NNODESL, struct_data%n_load_, 1, tc)


         iposNK = 1
         do ink = 1, NNODESL

            nk = struct_data%n_load_(ink)
            
            ! BUG: must be this because of anoher bug ahead..
            posKi = (nk - 1) * NLIBS
            ! posKi = (nk - 1) * NLIBSL + 1
            ! posKe = nk * NLIBSL


            S_K_fi    = S_IJK_fi   (1, iposNK)
            S_K_fj    = S_IJK_fj   (1, iposNK)
            S_K_fiPfj = S_IJK_fiPfj(1, iposNK)


            ! NOTE: use node index and NOT values
            !       since we only store actually loaded ones.
            !       This applies also to actually loaded LIBs.
            !       Since we use it all -> : syntax
            ! akU(1, 1, :) = wd%wfc_(tc, :, ink)
            ! ak(1, 1, :)  = wd%wfc_(tcP3, :, ink)


            akU(:, 1) = wd%wfc_(struct_data%libs_load_, tc,   ink)
            ak (:, 1) = wd%wfc_(struct_data%libs_load_, tcP3, ink)


            iposNJ = 1
            do inj = 1, NNODESL

               nj = struct_data%n_load_(inj)
               
               posJi = (nj - 1) * NLIBS
               ! posJi = (nj - 1) * NLIBSL + 1
               ! posJe = nj * NLIBSL

               S_J_fi    = S_IJK_fi   (1, iposNJ)
               S_J_fj    = S_IJK_fj   (1, iposNJ)
               S_J_fiPfj = S_IJK_fiPfj(1, iposNJ)

               ! BUG: should we put what??  DIR??  TC??
               corrJK = wd%nod_corr_(util_getCorrVectIndex(nj, nk, NNODES), tc)

               ! NOTE: we can precompute for perf
               S_JK_fi = corrJK**(abs_fi) * sqrt(S_J_fi * S_K_fi)
               S_JK_fj = corrJK**(abs_fj) * sqrt(S_J_fj * S_K_fj)


               ajU(1, :) = wd%wfc_(struct_data%libs_load_, tc,   inj)
               aj (1, :) = wd%wfc_(struct_data%libs_load_, tcP3, inj)



               iposNI = 1
               do ini = 1, NNODESL

                  ni = struct_data%n_load_(ini)
            
                  posIi = (ni - 1) * NLIBS
                  ! posIi = (ni - 1) * NLIBSL + 1
                  ! posIe = ni * NLIBSL

                  S_I_fi    = S_IJK_fi   (1, iposNI)
                  S_I_fj    = S_IJK_fj   (1, iposNI)
                  S_I_fiPfj = S_IJK_fiPfj(1, iposNI)


                  aiU(:, 1) = wd%wfc_(struct_data%libs_load_, tc,   ini)
                  ai (:, 1) = wd%wfc_(struct_data%libs_load_, tcP3, ini)


                  corrIJ = wd%nod_corr_(util_getCorrVectIndex(ni, nj, NNODES), tc)
                  corrIK = wd%nod_corr_(util_getCorrVectIndex(ni, nk, NNODES), tc)


                  ! term1
                  term1 = corrIJ**(abs_fi) * sqrt(S_I_fi * S_J_fi)
                  term1 = corrIK**(abs_fj) * sqrt(S_I_fj * S_K_fj) * term1

                  tmp1 = matmul(ai, ajU) * term1


                  ! term2
                  term2 = corrIJ**(abs_fiPfj) * sqrt(S_I_fiPfj * S_J_fiPfj)
                  term2 = S_JK_fj * term2

                  tmp2 = matmul(aiU, aj) * term2


                  ! term3
                  term3 = corrIK**(abs_fiPfj) * sqrt(S_I_fiPfj * S_K_fiPfj)
                  term3 = S_JK_fi * term3

                  tmp3 = matmul(aiU, ajU) * term3


                  ! BUG: apparently, cannot make a 3D product..
                  do ilk = 1, NLIBSL


                     phiK_ = struct_data%modal_%phi_(posKi + struct_data%libs_load_(ilk), MODES)


                     ! BUG: this formulation DOES NOT account for
                     !      interaction between turbulent components
                     !      (i.e. uv, uw, vw)
                     BF_IJK_ijk = 2 * &
                        ( &
                           (tmp1 * akU(ilk, 1)) + &
                           (tmp2 * akU(ilk, 1)) + &
                           (tmp3 * ak (ilk, 1))   & 
                        )


                     iposM = 1
                     do imk = 1, NMODES_EFF
                        
                        phiK = phiK_(imk)

                        do imj = 1, NMODES_EFF

                           phiJ(1, :) = struct_data%modal_%phi_(posJi + struct_data%libs_load_, MODES(imj))

                           do imi = 1, NMODES_EFF

                              phiI(:, 1) = struct_data%modal_%phi_(posIi + struct_data%libs_load_, MODES(imi))

                              ! bfm(iposM) = bfm(iposM) + &
                              !    sum(BF_IJK_ijk(:, :, :) * (matmul(phiI, phiJ) * phiK(:, :, :)))
                           
                              bfm(iposM) = bfm(iposM) + &
                                 sum(BF_IJK_ijk(:, :) * (matmul(phiI, phiJ) * phiK))

                              iposM = iposM + 1
                           enddo ! modes I

                        enddo ! modes J

                     enddo ! modes K


                  enddo ! lib K


                  iposNI = iposNI + 1
               enddo ! nodes I

               iposNJ = iposNJ + 1
            enddo ! nodes J

            iposNK = iposNK + 1
         enddo ! nodes K

      enddo ! itc

   end function getFM_full_tnm_scalar_msh_








   module subroutine prefetchSVDWorkDim_()
      double precision :: tmpmat(NNODESL, NNODESL)
      ! double precision, allocatable :: tmpmat(:, :)

      double precision, dimension(NNODESL) :: tmpv

      double precision, dimension(1) :: optWork
      double precision, dimension(1) :: tmp1arr

      integer :: istat
      character(len = 256) :: emsg

      interface
#ifdef BSA_USE_SVD_METHOD__
         SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, &
               VT, LDVT, WORK, LWORK, INFO )
            !    .. Scalar Arguments ..
            CHARACTER          JOBU, JOBVT
            INTEGER            LDA, LDU, LDVT, M, N
            integer            info, lwork
            !     .. Array Arguments ..
            DOUBLE PRECISION   A( LDA, * ), S( * ), U( LDU, * ),&
                                 VT( LDVT, * ), WORK( * )
         end subroutine
#else
         SUBROUTINE dsyev( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
            !     .. Scalar Arguments ..
                  CHARACTER          JOBZ, UPLO
                  INTEGER            INFO, LDA, LWORK, N
            !     ..
            !     .. Array Arguments ..
                  DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * )
         END SUBROUTINE
#endif
      end interface
      
      
      if (.not. allocated(MSHR_SVD_INFO)) then
         allocate(MSHR_SVD_INFO, stat=istat, errmsg=emsg)
         if (istat == 0) then
#ifdef __BSA_ALLOC_DEBUG
            call allocOKMsg('MSHR_SVD_INFO', loc(MSHR_SVD_INFO), sizeof(MSHR_SVD_INFO))
#endif
         else
            call allocKOMsg('MSHR_SVD_INFO', istat, emsg)
         endif
		endif
      
      MSHR_SVD_INFO = 0
#ifdef BSA_USE_SVD_METHOD__
      call dgesvd(&
           'O' &        ! min(M,N) columns of U are returned in array U
         , 'N' &        ! no rows of V are computed
         , NNODESL &    ! n. of rows M
         , NNODESL &    ! n. of cols N
         , tmpmat  &    ! A matrix
         , NNODESL &
         , tmpv    &
         , tmp1arr &    ! U array
         , 1       & 
         , tmp1arr &
         , 1       &
         , optWork &
         , MSHR_SVD_LWORK &
         , MSHR_SVD_INFO  &
      )
#else
      call dsyev('V', 'L', &
         NNODESL, tmpmat, NNODESL, tmp1arr, optWork, MSHR_SVD_LWORK, MSHR_SVD_INFO)
#endif

      if (MSHR_SVD_INFO == 0) then
         
         MSHR_SVD_LWORK = int(optWork(1), kind = 8)
! #ifdef __BSA_DEBUG
         print '(1x, a, a, i0 /)', &
            INFOMSG, 'WORK query ok. Optimal work dimension = ', MSHR_SVD_LWORK
! #endif
         
         if (allocated(MSHR_SVD_WORK)) then
            if (size(MSHR_SVD_WORK) /= MSHR_SVD_LWORK) then
               deallocate(MSHR_SVD_WORK)
               allocate(MSHR_SVD_WORK(MSHR_SVD_LWORK), stat=istat, errmsg=emsg)
               if (istat == 0) then
#ifdef __BSA_ALLOC_DEBUG
                  call allocOKMsg('MSHR_SVD_WORK', &
                     int(MSHR_SVD_LWORK), loc(MSHR_SVD_WORK), sizeof(MSHR_SVD_WORK))
#endif
               else
                  call allocKOMsg('MSHR_SVD_WORK', istat, emsg)
               endif
            endif
         else
            allocate(MSHR_SVD_WORK(MSHR_SVD_LWORK), stat=istat, errmsg=emsg)
            if (istat == 0) then
#ifdef __BSA_ALLOC_DEBUG
               call allocOKMsg('MSHR_SVD_WORK', &
                  int(MSHR_SVD_LWORK), loc(MSHR_SVD_WORK), sizeof(MSHR_SVD_WORK))
#endif
            else
               call allocKOMsg('MSHR_SVD_WORK', istat, emsg)
            endif
         endif
         return ! correct execution flow
      endif


      print '(1x, a, a, i0)', &
         ERRMSG, 'WORK query for SVD decomposition returned code   ', MSHR_SVD_INFO
      print '(1x, a, a)', &
         MSGCONT, 'Please, check again.'
      call bsa_Abort()
   end subroutine




   module subroutine cleanSVDWorkInfo_()
      integer :: istat
      character(len = 256) :: emsg

      ! NOTE: reset to -1 so that next call is going to query again.
      MSHR_SVD_LWORK = - 1

      if (allocated(MSHR_SVD_INFO)) then
         deallocate(MSHR_SVD_INFO, stat=istat, errmsg=emsg)
         if (istat == 0) then
#ifdef __BSA_ALLOC_DEBUG
            call deallocOKMsg('MSHR_SVD_INFO')
#endif
         else
            call deallocKOMsg('MSHR_SVD_INFO', istat, emsg)
         endif
      endif

      if (allocated(MSHR_SVD_WORK)) then
         deallocate(MSHR_SVD_WORK, stat=istat, errmsg=emsg)
         if (istat == 0) then
#ifdef __BSA_ALLOC_DEBUG
            call deallocOKMsg('MSHR_SVD_WORK')
#endif
         else
            call deallocKOMsg('MSHR_SVD_WORK', istat, emsg)
         endif
      endif

#ifdef __BSA_DEBUG
      print '(1x, a, a)', &
         INFOMSG, 'SVD related data cleaned -- ok.'
#endif
   end subroutine






   module function getFM_full_tm_scalar_msh_POD_(fi, fj) result(bfm)
      real(RDP), intent(in) :: fi, fj
      real(RDP) :: bfm(dimM_bisp_)

      real(RDP) :: fiPfj(1), fi_(1), fj_(1)
      
      ! wind turbulent comps indexes
      integer(kind = 4) :: itc, tc, tcP3

      ! n. of kept modes from wind fields decomposition
      integer :: nmw1, nmw2, nmw1w2
      integer :: p, q
      integer :: m, n, o, posm
      
      double precision, allocatable :: S_uvw_w1  (:, :)
      double precision, allocatable :: S_uvw_w2  (:, :)
      double precision, allocatable :: S_uvw_w1w2(:, :)
      
      ! tmp vec for interfacing, BUG: might be avoided?
      double precision :: tmpv(1, NNODESL)

      ! singular values vectors (DECREASING ordering!)
      double precision :: D_S_uvw_w1  (NNODESL)
      double precision :: D_S_uvw_w2  (NNODESL)
      double precision :: D_S_uvw_w1w2(NNODESL)

      double precision, dimension(NNODESL, 1)    :: eigvp, eigvq
      double precision, dimension(NMODES_EFF, 1) :: tmpm1, tmpm2, tmpm3
      double precision :: tmpDp, tmpTq, tmpo, tmpn

      interface
#ifdef BSA_USE_SVD_METHOD__
         SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, &
               VT, LDVT, WORK, LWORK, INFO )
            !    .. Scalar Arguments ..
            CHARACTER          JOBU, JOBVT
            INTEGER            LDA, LDU, LDVT, M, N
            integer            info, lwork
            !     .. Array Arguments ..
            DOUBLE PRECISION   A( LDA, * ), S( * ), U( LDU, * ),&
                                 VT( LDVT, * ), WORK( * )
         end subroutine
#else
         SUBROUTINE dsyev( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
            !     .. Scalar Arguments ..
                  CHARACTER          JOBZ, UPLO
                  INTEGER            INFO, LDA, LWORK, N
            !     ..
            !     .. Array Arguments ..
                  DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * )
         END SUBROUTINE
#endif
      end interface

      bfm      = 0._RDP
      fi_(1)   = fi
      fj_(1)   = fj
      fiPfj(1) = fi + fj


      allocate(S_uvw_w1  (NNODESL, NNODESL))
      allocate(S_uvw_w2  (NNODESL, NNODESL))
      allocate(S_uvw_w1w2(NNODESL, NNODESL))


      do itc = 1, NTCOMPS

         tc   = wd%tc_(itc)
         tcP3 = tc + 3

         !
         ! NODAL WIND TURBULENCEs PSDs (for given tc)
         !
         S_uvw_w1(:, 1:1) = &
            reshape(wd%evalPSD(1, fi_,   NNODESL, struct_data%n_load_, 1, tc), [NNODESL, 1])
         S_uvw_w2(:, 1:1) = &
            reshape(wd%evalPSD(1, fj_,   NNODESL, struct_data%n_load_, 1, tc), [NNODESL, 1])
         S_uvw_w1w2(:, 1:1) = &
            reshape(wd%evalPSD(1, fiPfj, NNODESL, struct_data%n_load_, 1, tc), [NNODESL, 1])

#ifdef __BSA_CHECK_NOD_COH_SVD
         if (itc == 1) then
            write(5482, *) fj_
            do nmw1 = 1, NNODESL
               write(5482, *) S_uvw_w2(nmw1, 1)
            enddo
         endif
#endif

         !
         ! applying spatial nodal coherence
         !
         S_uvw_w1   = wd%getFullNodalPSD(NNODESL, struct_data%n_load_, S_uvw_w1(:, 1),      fi, 1)
         S_uvw_w2   = wd%getFullNodalPSD(NNODESL, struct_data%n_load_, S_uvw_w2(:, 1),      fj, 1)
         S_uvw_w1w2 = wd%getFullNodalPSD(NNODESL, struct_data%n_load_, S_uvw_w1w2(:, 1), fiPfj(1), 1)

#ifdef __BSA_CHECK_NOD_COH_SVD
         if (itc == 1) then
            do nmw1 = 1, NNODESL
               write(5483, *) S_uvw_w2(:, nmw1)
            enddo
         endif
#endif

         !$omp critical
#ifdef BSA_USE_SVD_METHOD__
         call dgesvd(&
              'O' &           ! min(M,N) columns of U are overwritten on array A (saves memory)
            , 'N' &           ! no rows of V are computed
            , NNODESL    &    ! n. of rows M
            , NNODESL    &    ! n. of cols N
            , S_uvw_w1   &    ! A matrix (overwritten with left-singular vectors)
            , NNODESL    &
            , D_S_uvw_w1 &    ! singular values
            , tmpv       &    ! U
            , 1          & 
            , tmpv       &    ! VT
            , 1          &
            , MSHR_SVD_WORK  &
            , MSHR_SVD_LWORK &
            , MSHR_SVD_INFO  &
         )
#else
         call dsyev('V', 'L', &
            NNODESL, S_uvw_w1, NNODESL, D_S_uvw_w1, MSHR_SVD_WORK, MSHR_SVD_LWORK, MSHR_SVD_INFO)
#endif
         if (MSHR_SVD_INFO /= 0) then
            print '(1x, a, a, i0)', &
               ERRMSG, 'Error applying SVD to S_uvw_w1. Exit code  ', MSHR_SVD_INFO
            call bsa_Abort()
         endif


#ifdef BSA_USE_SVD_METHOD__
         call dgesvd(&
              'O' &           ! min(M,N) columns of U are overwritten on array A (saves memory)
            , 'N' &           ! no rows of V are computed
            , NNODESL    &    ! n. of rows M
            , NNODESL    &    ! n. of cols N
            , S_uvw_w2   &    ! A matrix (overwritten with left-singular vectors)
            , NNODESL    &
            , D_S_uvw_w2 &    ! singular values
            , tmpv       &    ! U
            , 1          & 
            , tmpv       &    ! VT
            , 1          &
            , MSHR_SVD_WORK  &
            , MSHR_SVD_LWORK &
            , MSHR_SVD_INFO  &
         )
#else
         call dsyev('V', 'L', &
            NNODESL, S_uvw_w2, NNODESL, D_S_uvw_w2, MSHR_SVD_WORK, MSHR_SVD_LWORK, MSHR_SVD_INFO)
#endif
         if (MSHR_SVD_INFO /= 0) then
            print '(1x, a, a, i0)', &
               ERRMSG, 'Error applying SVD to S_uvw_w2. Exit code  ', MSHR_SVD_INFO
            call bsa_Abort()
         endif

#ifdef __BSA_CHECK_NOD_COH_SVD
         if (itc == 1) then
            write(5484, *) NNODESL
            write(5484, *) D_S_uvw_w2
            do nmw1 = 1, NNODESL
               write(5484, *) S_uvw_w2(:, nmw1)
            enddo
         endif
#endif

#ifdef BSA_USE_SVD_METHOD__
         call dgesvd(&
              'O' &             ! min(M,N) columns of U are overwritten on array A (saves memory)
            , 'N' &             ! no rows of V are computed
            , NNODESL      &    ! n. of rows M
            , NNODESL      &    ! n. of cols N
            , S_uvw_w1w2   &    ! A matrix (overwritten with left-singular vectors)
            , NNODESL      &
            , D_S_uvw_w1w2 &    ! singular values
            , tmpv         &    ! U
            , 1            & 
            , tmpv         &    ! VT
            , 1            &
            , MSHR_SVD_WORK  &
            , MSHR_SVD_LWORK &
            , MSHR_SVD_INFO  &
         )
#else
         call dsyev('V', 'L', &
            NNODESL, S_uvw_w1w2, NNODESL, D_S_uvw_w1w2, MSHR_SVD_WORK, MSHR_SVD_LWORK, MSHR_SVD_INFO)
#endif
         if (MSHR_SVD_INFO /= 0) then
            print '(1x, a, a, i0)', &
               ERRMSG, 'Error applying SVD to S_uvw_w1w2. Exit code  ', MSHR_SVD_INFO
            call bsa_Abort()
         endif
         !$omp end critical

         
#ifdef __BSA_CHECK_NOD_COH_SVD
         return
#endif


         if (do_trunc_POD_) then
            nmw1 = 1
            if (.not. all(D_S_uvw_w1 == D_S_uvw_w1(1))) then
               tmpn = POD_trunc_lim_ * D_S_uvw_w1(1)
               nmw1 = 2
               do while (D_S_uvw_w1(nmw1) >= tmpn)
                  nmw1 = nmw1 + 1
               enddo
               nmw1 = nmw1 - 1
            endif

            nmw2 = 1
            if (.not. all(D_S_uvw_w2 == D_S_uvw_w2(1))) then
               tmpn = POD_trunc_lim_ * D_S_uvw_w2(1)
               nmw2 = 2
               do while (D_S_uvw_w2(nmw2) >= tmpn)
                  nmw2 = nmw2 + 1
               enddo
               nmw2 = nmw2 - 1
            endif

            nmw1w2 = 1
            if (.not. all(D_S_uvw_w1w2 == D_S_uvw_w1w2(1))) then
               tmpn = POD_trunc_lim_ * D_S_uvw_w1w2(1)
               nmw1w2 = 2
               do while (D_S_uvw_w1w2(nmw1w2) >= tmpn)
                  nmw1w2 = nmw1w2 + 1
               enddo
               nmw1w2 = nmw1w2 - 1
            endif

         else
            nmw1   = NNODESL
            nmw2   = nmw1
            nmw1w2 = nmw2
         endif


         ! 5-2-2, 2-2-5 (6-3-3, 3-3-6) (7-4-4, 4-4-7)
         do p = 1, nmw1

            eigvp(:, 1) = S_uvw_w1(:, p)

            ! V_lin_w1
            ! tmpm1 = matmul(wd%phi_times_A_ndegw_(:, :, tc), eigvp)
            do q = 1, NMODES_EFF
               tmpm1(q, 1) = sum(wd%phi_times_A_ndegw_(q, :, tc) * eigvp(:, 1))
            enddo 


            ! D_p_w1
            tmpDp = D_S_uvw_w1(p)


            ! 5-2-2 (6-3-3, 7-4-4)
            do q = 1, nmw2

               eigvq(:, 1) = S_uvw_w2(:, q)

               ! VZ_quad_w1w2
               tmpm2 = matmul(wd%phi_times_A_ndegw_(:, :, tcP3), eigvp * eigvq)

               ! Z_lin_w2
               tmpm3 = matmul(wd%phi_times_A_ndegw_(:, :, tc), eigvq)

               ! T_q_w2
               tmpTq = D_S_uvw_w2(q)


               posm = 1
               do o = 1, NMODES_EFF

                  tmpo = tmpm3(o, 1)

                  do n = 1, NMODES_EFF

                     tmpn = tmpm1(n, 1)

                     do m = 1, NMODES_EFF

                        bfm(posm) = bfm(posm) + &
                           (2 * tmpm2(m, 1) * &
                              tmpn * &
                              tmpo * tmpDp * tmpTq)

                        posm = posm + 1                        
                     enddo ! m modes

                  enddo ! n modes

               enddo ! o modes

            enddo ! q = 1, nmw2



            ! 2-2-5 (3-3-6, 4-4-7)
            do q = 1, nmw1w2

               eigvq(:, 1) = S_uvw_w1w2(:, q)

               ! Z_lin_w1w2
               tmpm2 = matmul(wd%phi_times_A_ndegw_(:, :, tc), eigvq)

               ! VZ_quad_w1w1w2
               tmpm3 = matmul(wd%phi_times_A_ndegw_(:, :, tcP3), eigvp * eigvq)

               ! T_q_w1w2
               tmpTq = D_S_uvw_w1w2(q)

               posm = 1
               do o = 1, NMODES_EFF

                  tmpo = tmpm3(o, 1)

                  do n = 1, NMODES_EFF

                     tmpn = tmpm1(n, 1)

                     do m = 1, NMODES_EFF

                        bfm(posm) = bfm(posm) + &
                           (2 * tmpm2(m, 1) * &
                              tmpn * &
                              tmpo * tmpDp * tmpTq)

                        posm = posm + 1
                     enddo ! m modes
                  enddo ! n modes

               enddo ! o modes
            enddo ! q = 1, nmw1w2

         enddo ! p = 1, nmw1



         ! 2-5-2 (3-6-3, 4-7-4)
         do p = 1, nmw1w2

            eigvp(:, 1) = S_uvw_w1w2(:, p)

            ! V_lin_w1w2
            tmpm1 = matmul(wd%phi_times_A_ndegw_(:, :, tc), eigvp)

            tmpDp = D_S_uvw_w1w2(p)

            do q = 1, nmw2

               eigvq(:, 1) = S_uvw_w2(:, q)

               ! Z_lin_w2
               tmpm2 = matmul(wd%phi_times_A_ndegw_(:, :, tc), eigvq)

               ! VZ_quad_w1w2w2
               tmpm3 = matmul(wd%phi_times_A_ndegw_(:, :, tcP3), eigvp * eigvq)

               tmpTq = D_S_uvw_w2(q)

               posm = 1
               do o = 1, NMODES_EFF

                  tmpo = tmpm2(o, 1)

                  do n = 1, NMODES_EFF

                     tmpn = tmpm3(n, 1)

                     do m = 1, NMODES_EFF

                        bfm(posm) = bfm(posm) + &
                           (2 * tmpm1(m, 1) * &
                              tmpn * &
                              tmpo * tmpDp * tmpTq)

                        posm = posm + 1
                     enddo ! m modes
                  enddo ! n modes
               enddo ! o modes

            enddo ! q = 1, nmw2
         enddo ! p = 1, nmw1w2

      enddo ! itc = 1, NTCOMPS

      ! !$omp critical
      ! !$ write(4382, *) omp_get_thread_num()
      ! write(4383, *) fi, fj
      ! write(4383, '(21g)') S_uvw_w1
      ! write(4383, *) ''
      ! write(4383, '(21g)') S_uvw_w1w2
      ! write(4383, *) MODES
      ! !$ write(4384, *) omp_get_thread_num()
      ! write(4384, '(g)') bfm
      ! write(4384, *) ''
      ! !$omp end critical
   end function getFM_full_tm_scalar_msh_POD_











   module function getRM_full_scalar_msh_(bfm, fi, fj) result(brm)
      real(RDP), intent(in) :: bfm(dimM_bisp_), fi, fj
      real(RDP) :: brm(dimM_bisp_)

      real(RDP) :: wi, wj, wiPwj
      integer(kind = 4)   :: posm, imk, imj, imi

      real(RDP), dimension(NMODES_EFF) :: Cdiag, rpart, ipart, htmp
      real(RDP), dimension(NMODES_EFF) :: H1r, H1i
      real(RDP), dimension(NMODES_EFF) :: H2r, H2i
      real(RDP), dimension(NMODES_EFF) :: H12r, H12i

      real(RDP) :: H12k_r, H12k_i, H2j_r, H2j_i

      wi    = fi * CST_PIt2
      wj    = fj * CST_PIt2
      wiPwj = wi + wj


      ! pre evaluate TFs (per mode)

      ! H1
      rpart = - (wi*wi * struct_data%modal_%Mm_(MODES)) + struct_data%modal_%Km_(MODES)
      do imi = 1, NMODES_EFF
         Cdiag(imi) = struct_data%modal_%Cm_(MODES(imi), MODES(imi))
      enddo
      ipart = Cdiag * wi
      htmp  = rpart*rpart + ipart*ipart
      H1r   =   rpart / htmp
      H1i   = - ipart / htmp

      ! H2
      rpart = - (wj*wj * struct_data%modal_%Mm_(MODES)) + struct_data%modal_%Km_(MODES)
      ipart = Cdiag * wj
      htmp  = rpart*rpart + ipart*ipart
      H2r   =   rpart / htmp
      H2i   = - ipart / htmp


      ! H12
      rpart = - (wiPwj*wiPwj * struct_data%modal_%Mm_(MODES)) + struct_data%modal_%Km_(MODES)
      ipart = Cdiag * wiPwj
      htmp  = rpart*rpart + ipart*ipart
      H12r  =   rpart / htmp
      H12i  = - ipart / htmp

      posm = 1
      do imk = 1, NMODES_EFF

         H12k_r = H12r(imk)
         H12k_i = H12i(imk)

         do imj = 1, NMODES_EFF

            H2j_r = H2r(imj)
            H2j_i = H2i(imj)

            do imi = 1, NMODES_EFF

               brm(posm) = bfm(posm) * &
                  (&
                     H1r(imi) * H2j_r * H12k_r + &
                     H1r(imi) * H2j_i * H12k_i + &
                     H1i(imi) * H2j_r * H12k_i - &
                     H1i(imi) * H2j_i * H12k_r   &
                  )

               posm = posm + 1
            enddo ! imi
         enddo ! imj
      enddo ! imk

   end function getRM_full_scalar_msh_











   module function getFM_diag_tnm_scalar_msh_(fi, fj) result(bfm)
      real(RDP), intent(in) :: fi, fj
      real(RDP) :: bfm(dimM_bisp_)

      real(RDP) :: fiPfj(1)

      integer(kind = 4) :: itc, tc, tcP3, posm, imode
      integer(kind = 4) :: posi, inode, node, ilibk

      real(RDP), dimension(1, NNODESL) :: Suvw_fi, Suvw_fj, Suvw_fiPfj
      real(RDP), dimension(NNODESL) :: Suvw_IJ, Suvw_IJI, Suvw_IJJ

      real(RDP) :: akU, ak, phik(NMODES_EFF)
      real(RDP), dimension(NLIBSL, 1) :: aiU, ai
      real(RDP), dimension(1, NLIBSL) :: ajU, aj
      real(RDP), dimension(NLIBSL, NMODES_EFF) :: phi_

      real(RDP) :: BF_ijk_I(NLIBSL, NLIBSL)
      real(RDP), dimension(NLIBSL, NLIBSL) :: tmp1, tmp2, tmp3

      bfm = 0._RDP

      fiPfj(1) = fi + fj

      do itc = 1, NTCOMPS

         tc   = wd%tc_(itc)
         tcP3 = tc + 3      ! quadratic term coeff

         Suvw_fi = wd%evalPSD(1, [fi], NNODESL, struct_data%n_load_, 1, tc)

         Suvw_fj = wd%evalPSD(1, [fj], NNODESL, struct_data%n_load_, 1, tc)

         Suvw_fiPfj = wd%evalPSD(1, fiPfj, NNODESL, struct_data%n_load_, 1, tc)


         ! precompute for perf
         Suvw_IJ  = Suvw_fi(1, :) * Suvw_fj(1, :)
         Suvw_IJI = Suvw_IJ(:) * Suvw_fi(1, :)
         Suvw_IJJ = Suvw_IJ(:) * Suvw_fj(1, :)


         do inode = 1, NNODESL

            node = int(struct_data%n_load_(inode), 4)
            
            posi = (node - 1) * NLIBS
            phi_ = struct_data%modal_%phi_(posi + struct_data%libs_load_, MODES)


            ajU(1, :) = wd%wfc_(struct_data%libs_load_, tc,   inode)
            aj (1, :) = wd%wfc_(struct_data%libs_load_, tcP3, inode)

            aiU(:, 1) = ajU(1, :)
            ai (:, 1) = aj (1, :)


            ! NOTE: this are tmp values !!!!!
            !       Done for performance.
            tmp1 = Suvw_IJ (inode) * matmul(ai , ajU)
            tmp2 = Suvw_IJJ(inode) * matmul(aiU, aj )
            tmp3 = Suvw_IJI(inode) * matmul(aiU, ajU)


            do ilibk = 1, NLIBSL

               phik = phi_(ilibk, :)

               akU = aiU(ilibk, 1)
               ak  = ai (ilibk, 1)


               BF_ijk_I = 2 * (&
                  tmp1 * akU + &
                  tmp2 * akU + &
                  tmp3 * ak    &
               )


               posm = 1
               do imode = 1, NMODES_EFF

                  bfm(posm) = bfm(posm) + &
                     sum( &
                        BF_ijk_I * &
                           (matmul(phi_(:, imode:imode), transpose(phi_(:, imode:imode))) &
                              * phik(imode)) )

                  posm = posm + 1
               enddo ! modes

            enddo ! libs loaded (k)
         enddo ! nodes loaded
      enddo ! n turb comps
      
   end function getFM_diag_tnm_scalar_msh_




   module function getRM_diag_scalar_msh_(bfm, fi, fj) result(brm)
      real(RDP), intent(in) :: bfm(dimM_bisp_), fi, fj
      real(RDP) :: brm(dimM_bisp_)

      real(RDP) :: wi, wj, wiPwj
      integer(kind = 4)   :: imi

      real(RDP), dimension(NMODES_EFF) :: Cdiag, rpart, ipart, htmp
      real(RDP), dimension(NMODES_EFF) :: H1r, H1i
      real(RDP), dimension(NMODES_EFF) :: H2r, H2i
      real(RDP), dimension(NMODES_EFF) :: H12r, H12i


      wi = fi * CST_PIt2
      wj = fj * CST_PIt2
      wiPwj = wi + wj


      ! pre evaluate TFs (per mode)

      ! H1
      rpart = - (wi*wi * struct_data%modal_%Mm_(MODES)) + struct_data%modal_%Km_(MODES)
      do imi = 1, NMODES_EFF
         Cdiag(imi) = struct_data%modal_%Cm_(MODES(imi), MODES(imi))
      enddo
      ipart = Cdiag * wi
      htmp  = rpart*rpart + ipart*ipart
      H1r   =   rpart / htmp
      H1i   = - ipart / htmp

      ! H2
      rpart = - (wj*wj * struct_data%modal_%Mm_(MODES)) + struct_data%modal_%Km_(MODES)
      ipart = Cdiag * wj
      htmp  = rpart*rpart + ipart*ipart
      H2r   =   rpart / htmp
      H2i   = - ipart / htmp


      ! H12
      rpart  = - (wiPwj*wiPwj * struct_data%modal_%Mm_(MODES)) + struct_data%modal_%Km_(MODES)
      ipart  = Cdiag * wiPwj
      htmp   = rpart*rpart + ipart*ipart
      H12r   =   rpart / htmp
      H12i   = - ipart / htmp

      brm = bfm * (&
         H1r * H2r * H12r + &
         H1r * H2i * H12i + &
         H1i * H2r * H12i - &
         H1i * H2i * H12r   &
      )

   end function getRM_diag_scalar_msh_











!!========================================================================================
!!========================================================================================
!!========================================================================================
!!
!!  classic
!!
!!========================================================================================
!!========================================================================================
!!========================================================================================





   !> BUG: this routine is adapted to the case where we use
   !>      convention on PULSATION.
   !>      Please, adapt it to the case of convention over FREQUENCIES.
   module subroutine getFM_full_tnlm_vect_cls_(f, Suvw, psd, bisp)
      real(RDP), intent(in) :: f(NFREQS)
      real(RDP), intent(in) :: Suvw(NFREQS, NPSDEL)
      real(RDP), allocatable, intent(inout) :: psd(:, :), bisp(:, :, :)

      integer(kind = 4) :: innl3
      integer(kind = 4) :: iin, ien, itmp, ifrj
      integer(kind = 4) :: i_n_pad, i_pad_len

      ! turb components related
      integer(kind = 4) :: itc, tc, tc_posN, tc_pk, tc_pj, tc_pi

      ! nodes indexed values
      integer(kind = 4) :: i_pos_nk, i_pos_nj, i_pos_ni
      integer(kind = 4) :: pos_nk, pos_nj, pos_ni
      integer(kind = 4) :: ink, inj, ini
      integer(kind = 4) :: ni, nj, nk

      ! libs indexed values
      integer(kind = 4) :: ilk, ilj, ili
      ! integer(kind = 4) :: li, lj, lk

      ! modes indexed values
      real(RDP), dimension(NLIBSL, NMODES_EFF) :: phik_, phij_, phii_
      real(RDP) :: phik, phij, phii
      integer(kind = 4) :: posm_
      integer(kind = 4) :: imk, imj, imi
      ! integer(kind = 4) :: mi, mj, mk

      integer(kind = 4) :: i_ncycles = 0

      real(RDP) :: f_abs(NFREQS)

      ! local nodal correlations
      real(RDP) :: corrJK, corrIK, corrIJ

      ! wfc extractions
      integer(kind = 4)   :: tcP3
      real(RDP), dimension(NLIBSL) :: aiU, ai, akU, ak
      real(RDP), dimension(NLIBSL) :: ajU, aj

      ! PSDs local
      real(RDP), allocatable :: S_uvw_i(:), S_uvw_j(:), S_uvw_k(:), PSDF_jk_JK_w(:)
      real(RDP), allocatable :: S_uvw_JK(:), S_uvw_IK(:), S_uvw_IJ(:)
      real(RDP), allocatable :: S_uvw_IK_w1w2(:), S_uvw_IJ_w1w2(:)

      ! BF local
      real(RDP), allocatable :: BF_ijk_IJK_w_w2(:)

      character(len = 256) :: emsg
      !========================================================================                                 


#ifdef __BSA_DEBUG
      write(unit_debug_, '(2a)') &
         INFOMSG, '@BsaClassicImpl::getFM_full_tnlm_vect_cls_() : computing modal forces spectra...'
#endif


      f_abs = abs(f)
      innl3 = NNODESL**3


      ! getting padded length and relative init/end indices (non zero zone)
      itmp      = NFREQS - 1     ! do not consider 0 (point of symmetry)
      i_n_pad   = itmp / 2       ! spread it on the two sides (left / right)
      ! iin      = i_n_pad + 1
      ! ien      = in + itmp
      ien       = i_n_pad + NFREQS
      iin       = i_n_pad + 1
      i_pad_len = itmp + NFREQS


#ifdef __BSA_DEBUG
      print '(1x, a, i5)', '@BsaClassicImpl::getFM_full_tnlm_vect_cls_() : i pad length = ', i_pad_len
      print '(1x, a, i5)', '@BsaClassicImpl::getFM_full_tnlm_vect_cls_() : init index   = ', iin
      print '(1x, a, i5)', '@BsaClassicImpl::getFM_full_tnlm_vect_cls_() : end  index   = ', ien
      print '(1x, a, i5)', '@BsaClassicImpl::getFM_full_tnlm_vect_cls_() : pad range    = ', ien - iin + 1
#endif


      ! these are needed regardlessly of if PSDs or BISPs

      allocate(psd(NFREQS, dimM_psd_), stat=ilk, errmsg=emsg)
      if (ilk == 0) then
#ifdef __BSA_ALLOC_DEBUG
         call allocOKMsg('psd', [NFREQS, dimM_psd_], loc(psd), sizeof(psd))
#endif
      else
         call allocKOMsg('psd', ilk, emsg)
      endif
      psd = 0._RDP

      allocate(S_uvw_k(NFREQS), stat=ilk, errmsg=emsg)
      if (ilk == 0) then
#ifdef __BSA_ALLOC_DEBUG
         call allocOKMsg('S_uvw_k', NFREQS, loc(S_uvw_k), sizeof(S_uvw_k))
#endif
      else
         call allocKOMsg('S_uvw_k', ilk, emsg)
      endif

      allocate(S_uvw_j(NFREQS), stat=ilk, errmsg=emsg)
      if (ilk == 0) then
#ifdef __BSA_ALLOC_DEBUG
         call allocOKMsg('S_uvw_j', NFREQS, loc(S_uvw_j), sizeof(S_uvw_j))
#endif
      else
         call allocKOMsg('S_uvw_j', ilk, emsg)
      endif

      allocate(S_uvw_JK(NFREQS), stat=ilk, errmsg=emsg)
      if (ilk == 0) then
#ifdef __BSA_ALLOC_DEBUG
         call allocOKMsg('S_uvw_JK', NFREQS, loc(S_uvw_JK), sizeof(S_uvw_JK))
#endif
      else
         call allocKOMsg('S_uvw_JK', ilk, emsg)
      endif

      allocate(PSDF_jk_JK_w(NFREQS), stat=ilk, errmsg=emsg)
      if (ilk == 0) then
#ifdef __BSA_ALLOC_DEBUG
         call allocOKMsg('PSDF_jk_JK_w', NFREQS, loc(PSDF_jk_JK_w), sizeof(PSDF_jk_JK_w))
#endif
      else
         call allocKOMsg('PSDF_jk_JK_w', ilk, emsg)
      endif

      if (settings%i_compute_bisp_ == 1) then

         allocate(bisp(NFREQS, NFREQS, dimM_bisp_), stat=ilk, errmsg=emsg)
         if (ilk == 0) then
#ifdef __BSA_ALLOC_DEBUG
            call allocOKMsg('bisp', [NFREQS, NFREQS, dimM_bisp_], loc(bisp), sizeof(bisp))
#endif
         else
            call allocKOMsg('bisp', ilk, emsg)
         endif
         bisp = 0._RDP

         allocate(bf_ijk_IJK_w_w2(NFREQS), stat=ilk, errmsg=emsg)
         if (ilk == 0) then
#ifdef __BSA_ALLOC_DEBUG
            call allocOKMsg(&
               'bf_ijk_IJK_w_w2', NFREQS, loc(bf_ijk_IJK_w_w2), sizeof(bf_ijk_IJK_w_w2))
#endif
         else
            call allocKOMsg('bf_ijk_IJK_w_w2', ilk, emsg)
         endif

         allocate(S_uvw_i(NFREQS), stat=ilk, errmsg=emsg)
         if (ilk == 0) then
#ifdef __BSA_ALLOC_DEBUG
            call allocOKMsg('S_uvw_i', NFREQS, loc(S_uvw_i), sizeof(S_uvw_i))
#endif
         else
            call allocKOMsg('S_uvw_i', ilk, emsg)
         endif

         allocate(S_uvw_IK(NFREQS), stat=ilk, errmsg=emsg)
         if (ilk == 0) then
#ifdef __BSA_ALLOC_DEBUG
            call allocOKMsg('S_uvw_IK', NFREQS, loc(S_uvw_IK), sizeof(S_uvw_IK))
#endif
         else
            call allocKOMsg('S_uvw_IK', ilk, emsg)
         endif

         allocate(S_uvw_IJ(NFREQS), stat=ilk, errmsg=emsg)
         if (ilk == 0) then
#ifdef __BSA_ALLOC_DEBUG
            call allocOKMsg('S_uvw_IJ', NFREQS, loc(S_uvw_IJ), sizeof(S_uvw_IJ))
#endif
         else
            call allocKOMsg('S_uvw_IJ', ilk, emsg)
         endif

         allocate(S_uvw_IK_w1w2(i_pad_len), stat=ilk, errmsg=emsg)
         if (ilk == 0) then
#ifdef __BSA_ALLOC_DEBUG
            call allocOKMsg('S_uvw_IK_w1w2', i_pad_len, loc(S_uvw_IK_w1w2), sizeof(S_uvw_IK_w1w2))
#endif
         else
            call allocKOMsg('S_uvw_IK_w1w2', ilk, emsg)
         endif

         allocate(S_uvw_IJ_w1w2(i_pad_len), stat=ilk, errmsg=emsg)
         if (ilk == 0) then
#ifdef __BSA_ALLOC_DEBUG
            call allocOKMsg('S_uvw_IJ_w1w2', i_pad_len, loc(S_uvw_IJ_w1w2), sizeof(S_uvw_IJ_w1w2))
#endif
         else
            call allocKOMsg('S_uvw_IJ_w1w2', ilk, emsg)
         endif

      endif



      !========================================================================
      ! BUG: for the moment, only considering correlation
      !      between same turbulence component (u, v, w).
      !      No cross-correlation between turbulent components
      !      i.e. E[uv]==E[uw]==E[vw] === 0
      do itc = 1, NTCOMPS

         tc      = wd%tc_(itc) ! get actual turbulent component
         tcP3    = tc + 3   ! quadratic term coeff
         tc_posN = (itc - 1) * NNODESL


         i_pos_nk = 1
         do ink = 1, NNODESL

            nk     = struct_data%n_load_(ink)
            pos_nk = (nk - 1) * NLIBS
            tc_pk  = tc_posN + i_pos_nk

            phik_ = struct_data%modal_%phi_(pos_nk + struct_data%libs_load_, MODES)

            akU(:) = wd%wfc_(struct_data%libs_load_, tc,   ink)
            ak (:) = wd%wfc_(struct_data%libs_load_, tcP3, ink)

            S_uvw_k = Suvw(:, tc_pk)
            ! if (settings%i_only_psd_ == 0) S_uvw_pad_k(iin : ien) = S_uvw_k

            
            i_pos_nj = 1
            do inj = 1, NNODESL

               nj     = struct_data%n_load_(inj)
               pos_nj = (nj - 1) * NLIBS
               tc_pj  = tc_posN + i_pos_nj

               phij_ = struct_data%modal_%phi_(pos_nj + struct_data%libs_load_, MODES)

               ajU(:) = wd%wfc_(struct_data%libs_load_, tc,   inj)
               aj (:) = wd%wfc_(struct_data%libs_load_, tcP3, inj)

               ! BUG: inserted itc, was 1
               corrJK = wd%nod_corr_(util_getCorrVectIndex(nj, nk, NNODES), tc)


               S_uvw_j = Suvw(:, tc_pj)
               ! if (settings%i_only_psd_ == 0) S_uvw_pad_j(iin : ien) = S_uvw_j

               S_uvw_JK = corrJK**(f_abs) * sqrt(S_uvw_k * S_uvw_j)


               !! BISPs
               if (settings%i_compute_bisp_ == 1) then

                  i_pos_ni = 1
                  do ini = 1, NNODESL

                     ni     = struct_data%n_load_(ini)
                     pos_ni = (ni - 1) * NLIBS
                     tc_pi  = tc_posN + i_pos_ni

                     phii_ = struct_data%modal_%phi_(pos_ni + struct_data%libs_load_, MODES)

                     aiU(:) = wd%wfc_(struct_data%libs_load_, tc,   ini)
                     ai (:) = wd%wfc_(struct_data%libs_load_, tcP3, ini)
      
                     corrIK = wd%nod_corr_(util_getCorrVectIndex(ni, nk, NNODES), tc)
                     corrIJ = wd%nod_corr_(util_getCorrVectIndex(ni, nj, NNODES), tc)

                     S_uvw_i = Suvw(:, tc_pi)


                     S_uvw_IK = corrIK**(f_abs) * sqrt(S_uvw_i * S_uvw_k)
                     S_uvw_IK_w1w2(iin : ien) = S_uvw_IK


                     S_uvw_IJ = corrIJ**(f_abs) * sqrt(S_uvw_i * S_uvw_j)
                     S_uvw_IJ_w1w2(iin : ien) = S_uvw_IJ



                     ! loop on frequencies (second dimension, j)
                     itmp = NFREQS
                     do ifrj = 1, NFREQS


                        do ilk = 1, NLIBSL

                           ! lk = struct_data%libs_load_(ilk)

                           do ilj = 1, NLIBSL

                              ! lj   = struct_data%libs_load_(ilj)


                              do ili = 1, NLIBSL

                                 ! li   = struct_data%libs_load_(ili)

                                 
                                 BF_ijk_IJK_w_w2 = 2 * (&
                                    ai (ili) * ajU(ilj) * akU(ilk) * (S_uvw_IJ * S_uvw_IK(ifrj)) + &
                                    aiU(ili) * aj (ilj) * akU(ilk) * (S_uvw_IJ_w1w2(ifrj : itmp) * S_uvw_JK(ifrj)) + &
                                    aiU(ili) * ajU(ilj) * ak (ilk) * (S_uvw_JK * S_uvw_IK_w1w2(ifrj : itmp)) &
                                 &)


                                 ! if (all(BF_ijk_IJK_w_w2 == 0._RDP)) cycle


                                 posm_ = 1
                                 do imk = 1, NMODES_EFF

                                    ! mk   = struct_data%modal_%modes_(imk)
                                    phik = phik_(ilk, imk)

                                    do imj = 1, NMODES_EFF

                                       ! mj   = struct_data%modal_%modes_(imj)
                                       phij = phij_(ilj, imj)


                                       ! TODO: this loop can be suppressed
                                       do imi = 1, NMODES_EFF

                                          ! mi   = struct_data%modal_%modes_(imi)
                                          phii = phii_(ili, imi)

                                          bisp(:, ifrj, posm_) = bisp(:, ifrj, posm_) + &
                                             phik * phij * phii * BF_ijk_IJK_w_w2

                                          posm_ = posm_ + 1
                                       enddo ! i mode
                                    enddo ! j mode
                                 enddo ! k mode                     

                              enddo ! i lib
                           enddo ! j lib
                        enddo ! k lib

                        itmp = itmp + 1
                     enddo ! n freqs j


                     i_pos_ni = i_pos_ni + 1
                  enddo ! i node

#ifdef __BSA_DEBUG
                  i_ncycles = i_ncycles + NNODESL
                  print '(1x, a, a, f10.4, " %")', &
                     INFOMSG, 'getFM_full_tnlm_vect_cls_() :   done  ', &
                        real(i_ncycles, RDP)/innl3*100
#endif

               endif ! bisp computation



               !! PSDs
               do ilk = 1, NLIBSL

                  ! lk   = struct_data%libs_load_(ilk)
                  ! if (akU(ilk) == 0.0_RDP) cycle


                  do ilj = 1, NLIBSL

                     ! lj   = struct_data%libs_load_(ilj)
                     ! if (ajU(ilj) == 0.0_RDP) cycle


                     ! PSD f
                     PSDF_jk_JK_w = ajU(ilj) * akU(ilk) * S_uvw_JK
                     
                     ! if (all(PSDF_jk_JK_w == 0.0_RDP)) cycle


                     posm_ = 1
                     do imk = 1, NMODES_EFF

                        ! mk   = struct_data%modal_%modes_(imk)
                        phik = phik_(ilk, imk)
                        ! if (phik == 0.0_RDP) cycle

                        do imj = 1, NMODES_EFF

                           ! mj   = struct_data%modal_%modes_(imj)
                           phij = phij_(ilj, imj)
                           ! if (phij == 0.0_RDP) cycle

                           psd(:, posm_) = psd(:, posm_) + &
                              phik * phij * PSDF_jk_JK_w

! #ifdef __BSA_DEBUG
!                            write(unit_debug_, &
! 										   '(1x, a, 5(i0, ", "), i0,  "  ; ",  2(2x, g0, " - ", g0) )') &
! 										'  nk, nj, lk, lj, mk, mj :  ', &
!                               nk, nj, &
!                               struct_data%libs_load_(ilk), struct_data%libs_load_(ilj), &
!                               imk, imj, &
! 										akU(ilk), ajU(ilj), &
! 										phik, phij
! #endif


                           posm_ = posm_ + 1
                        enddo ! j mode
                     enddo ! k mode

                  enddo ! j lib
               enddo ! k lib


               i_pos_nj = i_pos_nj + 1
            enddo ! j node


            i_pos_nk = i_pos_nk + 1
         enddo ! k node


      enddo ! itc



      ! deallocation
      if (allocated(S_uvw_i)) deallocate(S_uvw_i)
      if (allocated(S_uvw_j)) deallocate(S_uvw_j)
      if (allocated(S_uvw_k)) deallocate(S_uvw_k)
      if (allocated(PSDF_jk_JK_w)) deallocate(PSDF_jk_JK_w)
      if (allocated(S_uvw_JK)) deallocate(S_uvw_JK)
      if (allocated(S_uvw_IK)) deallocate(S_uvw_IK)
      if (allocated(S_uvw_IJ)) deallocate(S_uvw_IJ)
      if (allocated(S_uvw_IK_w1w2)) deallocate(S_uvw_IK_w1w2)
      if (allocated(S_uvw_IJ_w1w2)) deallocate(S_uvw_IJ_w1w2)
      if (allocated(BF_ijk_IJK_w_w2)) deallocate(BF_ijk_IJK_w_w2)

#ifdef __BSA_DEBUG
      write(unit_debug_, '(2a)') &
         INFOMSG, '@BsaClassicImpl::getFM_full_tnlm_vect_cls_() : computing modal forces spectra -- ok.'
#endif
   end subroutine getFM_full_tnlm_vect_cls_













   !> BUG: this routine is adapted to the case where we use
   !>      convention on PULSATION.
   !>      Please, adapt it to the case of convention over FREQUENCIES.
   module subroutine getFM_full_tnm_vect_cls_(f, Suvw, psd, bisp)
      real(RDP), intent(in) :: f(NFREQS)
      real(RDP), intent(in) :: Suvw(NFREQS, NPSDEL)
      real(RDP), allocatable, intent(inout) :: psd(:, :), bisp(:, :, :)

      integer(kind = 4) :: innl3
      integer(kind = 4) :: iin, ien, itmp, ifrj
      integer(kind = 4) :: i_n_pad, i_pad_len

      ! turb components related
      integer(kind = 4) :: itc, tc, tc_posN, tc_pk, tc_pj, tc_pi

      ! nodes indexed values
      integer(kind = 4) :: i_pos_nk, i_pos_nj, i_pos_ni
      integer(kind = 4) :: pos_nk, pos_nj, pos_ni
      integer(kind = 4) :: ink, inj, ini
      integer(kind = 4) :: ni, nj, nk

      ! modes indexed values
      real(RDP), dimension(NMODES_EFF, 2) :: phik_, phij_, phii_
      real(RDP) :: phij_Ub_, phij_u_, phik_Ub_, phik_u_
      integer(kind = 4) :: posm_
      integer(kind = 4) :: imk, imj, imi

      integer(kind = 4) :: i_ncycles = 0

      real(RDP) :: f_abs(NFREQS)

      ! local nodal correlations
      real(RDP) :: corrJK, corrIK, corrIJ

      ! wfc extractions
      integer(kind = 4)   :: tcP3

      ! PSDs local
      real(RDP), allocatable :: S_uvw_i(:), S_uvw_j(:), S_uvw_k(:), PSDF_jk_JK_w(:)
      real(RDP), allocatable :: S_uvw_JK(:), S_uvw_IK(:), S_uvw_IJ(:)
      real(RDP), allocatable :: S_uvw_IK_w1w2(:), S_uvw_IJ_w1w2(:)

      ! BF local
      real(RDP), allocatable :: BF_ijk_IJK_w_w2(:), tmp1(:), tmp2(:), tmp3(:)

      character(len = 256) :: emsg
      !========================================================================                                 


#ifdef __BSA_DEBUG
      write(unit_debug_, '(2a)') &
         INFOMSG, '@BsaClassicImpl::getFM_full_tnm_vect_cls_() : computing modal forces spectra...'
#endif


      f_abs = abs(f)
      innl3 = NNODESL**3


      ! getting padded length and relative init/end indices (non zero zone)
      itmp      = NFREQS - 1     ! do not consider 0 (point of symmetry)
      i_n_pad   = itmp / 2       ! spread it on the two sides (left / right)
      ! iin      = i_n_pad + 1
      ! ien      = in + itmp
      ien       = i_n_pad + NFREQS
      iin       = i_n_pad + 1
      i_pad_len = itmp + NFREQS


#ifdef __BSA_DEBUG
      print '(1x, a, i5)', '@BsaClassicImpl::getFM_full_tnm_vect_cls_() : i pad length = ', i_pad_len
      print '(1x, a, i5)', '@BsaClassicImpl::getFM_full_tnm_vect_cls_() : init index   = ', iin
      print '(1x, a, i5)', '@BsaClassicImpl::getFM_full_tnm_vect_cls_() : end  index   = ', ien
      print '(1x, a, i5)', '@BsaClassicImpl::getFM_full_tnm_vect_cls_() : pad range    = ', ien - iin + 1
#endif


      ! these are needed regardlessly of if PSDs or BISPs

      allocate(psd(NFREQS, dimM_psd_), stat=itc, errmsg=emsg)
      if (itc == 0) then
#ifdef __BSA_ALLOC_DEBUG
         call allocOKMsg('psd', [NFREQS, dimM_psd_], loc(psd), sizeof(psd))
#endif
      else
         call allocKOMsg('psd', itc, emsg)
      endif
      psd = 0._RDP

      allocate(S_uvw_k(NFREQS), stat=itc, errmsg=emsg)
      if (itc == 0) then
#ifdef __BSA_ALLOC_DEBUG
         call allocOKMsg('S_uvw_k', NFREQS, loc(S_uvw_k), sizeof(S_uvw_k))
#endif
      else
         call allocKOMsg('S_uvw_k', itc, emsg)
      endif

      allocate(S_uvw_j(NFREQS), stat=itc, errmsg=emsg)
      if (itc == 0) then
#ifdef __BSA_ALLOC_DEBUG
         call allocOKMsg('S_uvw_j', NFREQS, loc(S_uvw_j), sizeof(S_uvw_j))
#endif
      else
         call allocKOMsg('S_uvw_j', itc, emsg)
      endif

      allocate(S_uvw_JK(NFREQS), stat=itc, errmsg=emsg)
      if (itc == 0) then
#ifdef __BSA_ALLOC_DEBUG
         call allocOKMsg('S_uvw_JK', NFREQS, loc(S_uvw_JK), sizeof(S_uvw_JK))
#endif
      else
         call allocKOMsg('S_uvw_JK', itc, emsg)
      endif

      allocate(PSDF_jk_JK_w(NFREQS), stat=itc, errmsg=emsg)
      if (itc == 0) then
#ifdef __BSA_ALLOC_DEBUG
         call allocOKMsg('PSDF_jk_JK_w', NFREQS, loc(PSDF_jk_JK_w), sizeof(PSDF_jk_JK_w))
#endif
      else
         call allocKOMsg('PSDF_jk_JK_w', itc, emsg)
      endif

      if (settings%i_compute_bisp_ == 1) then

         allocate(bisp(NFREQS, NFREQS, dimM_bisp_), stat=itc, errmsg=emsg)
         if (itc == 0) then
#ifdef __BSA_ALLOC_DEBUG
            call allocOKMsg('bisp', [NFREQS, NFREQS, dimM_bisp_], loc(bisp), sizeof(bisp))
#endif
         else
            call allocKOMsg('bisp', itc, emsg)
         endif
         bisp = 0._RDP

         allocate(bf_ijk_IJK_w_w2(NFREQS), stat=itc, errmsg=emsg)
         if (itc == 0) then
#ifdef __BSA_ALLOC_DEBUG
            call allocOKMsg(&
               'bf_ijk_IJK_w_w2', NFREQS, loc(bf_ijk_IJK_w_w2), sizeof(bf_ijk_IJK_w_w2))
#endif
         else
            call allocKOMsg('bf_ijk_IJK_w_w2', itc, emsg)
         endif

         allocate(tmp1(NFREQS), stat=itc, errmsg=emsg)
         if (itc == 0) then
#ifdef __BSA_ALLOC_DEBUG
            call allocOKMsg(&
               'tmp1', NFREQS, loc(tmp1), sizeof(tmp1))
#endif
         else
            call allocKOMsg('tmp1', itc, emsg)
         endif

         allocate(tmp2(NFREQS), stat=itc, errmsg=emsg)
         if (itc == 0) then
#ifdef __BSA_ALLOC_DEBUG
            call allocOKMsg(&
               'tmp2', NFREQS, loc(tmp2), sizeof(tmp2))
#endif
         else
            call allocKOMsg('tmp2', itc, emsg)
         endif

         allocate(tmp3(NFREQS), stat=itc, errmsg=emsg)
         if (itc == 0) then
#ifdef __BSA_ALLOC_DEBUG
            call allocOKMsg(&
               'tmp3', NFREQS, loc(tmp3), sizeof(tmp3))
#endif
         else
            call allocKOMsg('tmp3', itc, emsg)
         endif

         allocate(S_uvw_i(NFREQS), stat=itc, errmsg=emsg)
         if (itc == 0) then
#ifdef __BSA_ALLOC_DEBUG
            call allocOKMsg('S_uvw_i', NFREQS, loc(S_uvw_i), sizeof(S_uvw_i))
#endif
         else
            call allocKOMsg('S_uvw_i', itc, emsg)
         endif

         allocate(S_uvw_IK(NFREQS), stat=itc, errmsg=emsg)
         if (itc == 0) then
#ifdef __BSA_ALLOC_DEBUG
            call allocOKMsg('S_uvw_IK', NFREQS, loc(S_uvw_IK), sizeof(S_uvw_IK))
#endif
         else
            call allocKOMsg('S_uvw_IK', itc, emsg)
         endif

         allocate(S_uvw_IJ(NFREQS), stat=itc, errmsg=emsg)
         if (itc == 0) then
#ifdef __BSA_ALLOC_DEBUG
            call allocOKMsg('S_uvw_IJ', NFREQS, loc(S_uvw_IJ), sizeof(S_uvw_IJ))
#endif
         else
            call allocKOMsg('S_uvw_IJ', itc, emsg)
         endif

         allocate(S_uvw_IK_w1w2(i_pad_len), stat=itc, errmsg=emsg)
         if (itc == 0) then
#ifdef __BSA_ALLOC_DEBUG
            call allocOKMsg('S_uvw_IK_w1w2', i_pad_len, loc(S_uvw_IK_w1w2), sizeof(S_uvw_IK_w1w2))
#endif
         else
            call allocKOMsg('S_uvw_IK_w1w2', itc, emsg)
         endif

         allocate(S_uvw_IJ_w1w2(i_pad_len), stat=itc, errmsg=emsg)
         if (itc == 0) then
#ifdef __BSA_ALLOC_DEBUG
            call allocOKMsg('S_uvw_IJ_w1w2', i_pad_len, loc(S_uvw_IJ_w1w2), sizeof(S_uvw_IJ_w1w2))
#endif
         else
            call allocKOMsg('S_uvw_IJ_w1w2', itc, emsg)
         endif

      endif ! i bisp allocation



      !========================================================================
      ! BUG: for the moment, only considering correlation
      !      between same turbulence component (u, v, w).
      !      No cross-correlation between turbulent components
      !      i.e. E[uv]==E[uw]==E[vw] === 0
      do itc = 1, NTCOMPS

         tc      = wd%tc_(itc) ! get actual turbulent component
         tcP3    = tc + 3   ! quadratic term coeff
         tc_posN = (itc - 1) * NNODESL


         i_pos_nk = 1
         do ink = 1, NNODESL

            nk     = struct_data%n_load_(ink)
            pos_nk = (nk - 1) * NLIBS
            tc_pk  = tc_posN + i_pos_nk

            phik_(:, 1) = wd%phi_times_A_ndegw_(:, ink, tc  )
            phik_(:, 2) = wd%phi_times_A_ndegw_(:, ink, tcP3)

            S_uvw_k = Suvw(:, tc_pk)
            
            i_pos_nj = 1
            do inj = 1, NNODESL

               nj     = struct_data%n_load_(inj)
               pos_nj = (nj - 1) * NLIBS
               tc_pj  = tc_posN + i_pos_nj

               phij_(:, 1) = wd%phi_times_A_ndegw_(:, inj, tc  )
               phij_(:, 2) = wd%phi_times_A_ndegw_(:, inj, tcP3)

               corrJK = wd%nod_corr_(util_getCorrVectIndex(nj, nk, NNODES), tc)

               S_uvw_j  = Suvw(:, tc_pj)
               S_uvw_JK = corrJK**(f_abs) * sqrt(S_uvw_k * S_uvw_j)


               !! BISPs
               if (settings%i_compute_bisp_ == 1) then

                  i_pos_ni = 1
                  do ini = 1, NNODESL

                     ni     = struct_data%n_load_(ini)
                     pos_ni = (ni - 1) * NLIBS
                     tc_pi  = tc_posN + i_pos_ni

                     phii_(:, 1) = wd%phi_times_A_ndegw_(:, ini, tc  )
                     phii_(:, 2) = wd%phi_times_A_ndegw_(:, ini, tcP3)
      
                     corrIK = wd%nod_corr_(util_getCorrVectIndex(ni, nk, NNODES), tc)
                     corrIJ = wd%nod_corr_(util_getCorrVectIndex(ni, nj, NNODES), tc)

                     S_uvw_i  = Suvw(:, tc_pi)

                     S_uvw_IK = corrIK**(f_abs) * sqrt(S_uvw_i * S_uvw_k)
                     S_uvw_IK_w1w2(iin : ien) = S_uvw_IK

                     S_uvw_IJ = corrIJ**(f_abs) * sqrt(S_uvw_i * S_uvw_j)
                     S_uvw_IJ_w1w2(iin : ien) = S_uvw_IJ


                     ! loop on frequencies (second dimension, j)
                     itmp = NFREQS
                     do ifrj = 1, NFREQS

                        tmp1 = S_uvw_IJ                   * S_uvw_IK(ifrj)
                        tmp2 = S_uvw_IJ_w1w2(ifrj : itmp) * S_uvw_JK(ifrj)
                        tmp3 = S_uvw_JK                   * S_uvw_IK_w1w2(ifrj : itmp)

                        posm_ = 1
                        do imk = 1, NMODES_EFF

                           phik_Ub_ = phik_(imk, 1)
                           phik_u_  = phik_(imk, 2)

                           do imj = 1, NMODES_EFF

                              phij_Ub_ = phij_(imj, 1)
                              phij_u_  = phij_(imj, 2)

                              ! TODO: this loop can be suppressed
                              do imi = 1, NMODES_EFF

                                 bisp(:, ifrj, posm_) = bisp(:, ifrj, posm_) + &
                                    2 * ( &
                                          (phii_(imi, 2) * phij_Ub_ * phik_Ub_ * tmp1)  &
                                       +  (phii_(imi, 1) * phij_u_  * phik_Ub_ * tmp2)  &
                                       +  (phii_(imi, 1) * phij_Ub_ * phik_u_  * tmp3)  &
                                    )

                                 posm_ = posm_ + 1
                              enddo ! i mode
                           enddo ! j mode
                        enddo ! k mode

                        itmp = itmp + 1
                     enddo ! n freqs j


                     i_pos_ni = i_pos_ni + 1
                  enddo ! i node

! #ifdef __BSA_DEBUG
                  i_ncycles = i_ncycles + NNODESL
                  print '(1x, a, a, f10.4, " %")', &
                     INFOMSG, 'getFM_full_tnm_vect_cls_() :   done  ', &
                        real(i_ncycles, RDP)/innl3*100
! #endif

               endif ! bisp computation



               posm_ = 1
               do imk = 1, NMODES_EFF

                  phik_Ub_ = phik_(imk, 1)

                  do imj = 1, NMODES_EFF

                     psd(:, posm_) = psd(:, posm_) + &
                        (phik_Ub_ * phij_(imj, 1) * S_uvw_JK)

                     posm_ = posm_ + 1
                  enddo ! j mode
               enddo ! k mode


               i_pos_nj = i_pos_nj + 1
            enddo ! j node


            i_pos_nk = i_pos_nk + 1
         enddo ! k node

      enddo ! itc



      ! deallocation
      if (allocated(S_uvw_i)) deallocate(S_uvw_i)
      if (allocated(S_uvw_j)) deallocate(S_uvw_j)
      if (allocated(S_uvw_k)) deallocate(S_uvw_k)
      if (allocated(PSDF_jk_JK_w)) deallocate(PSDF_jk_JK_w)
      if (allocated(S_uvw_JK)) deallocate(S_uvw_JK)
      if (allocated(S_uvw_IK)) deallocate(S_uvw_IK)
      if (allocated(S_uvw_IJ)) deallocate(S_uvw_IJ)
      if (allocated(S_uvw_IK_w1w2)) deallocate(S_uvw_IK_w1w2)
      if (allocated(S_uvw_IJ_w1w2)) deallocate(S_uvw_IJ_w1w2)
      if (allocated(BF_ijk_IJK_w_w2)) deallocate(BF_ijk_IJK_w_w2)

#ifdef __BSA_DEBUG
      write(unit_debug_, '(2a)') &
         INFOMSG, '@BsaClassicImpl::getFM_full_tnm_vect_cls_() : computing modal forces spectra -- ok.'
#endif
   end subroutine getFM_full_tnm_vect_cls_










   module subroutine getRM_full_vect_cls_(f, psd, bisp)
      real(RDP), intent(in)                 :: f(NFREQS)
      real(RDP), allocatable, intent(inout) :: psd(:, :), bisp(:, :, :)

      integer(kind = 4) :: ifrj
      integer(kind = 4) :: posm_psd = 1, posm_bisp = 1
      
      ! modal indexed
      integer(kind = 4) :: imi, imj, imk, mi

      real(RDP) :: omegas(NFREQS, 1)
      real(RDP), allocatable :: r_part(:, :), i_part(:, :), h_tmp(:, :), h_tmp2(:, :)
      real(RDP), allocatable :: Hr_w(:, :), Hi_w(:, :)
      real(RDP), allocatable :: Hr_w1w2(:, :, :), Hi_w1w2(:, :, :)


#ifdef __BSA_DEBUG
      write(unit_debug_, '(2a)') &
         INFOMSG, '@BsaClassicImpl::getRM_full_vect_cls_() : computing modal responses spectra...'
#endif

      ! BUG: check logic if correct
      if (settings%i_compute_bisp_ == 1) then

         allocate(r_part(NFREQS, NFREQS))
         allocate(i_part(NFREQS, NFREQS))
         allocate(h_tmp(NFREQS, NFREQS))
         allocate(h_tmp2(NFREQS, NFREQS))
         allocate(Hr_w1w2(NFREQS, NFREQS, NMODES_EFF))
         allocate(Hi_w1w2(NFREQS, NFREQS, NMODES_EFF))

      elseif (settings%i_compute_psd_ == 1) then

         allocate(r_part(NFREQS, 1))
         allocate(i_part(NFREQS, 1))
         allocate(h_tmp(NFREQS, 1))

      else ! none of them ??
         return
      endif

      
      ! TRANSFER FUNCTION COMPUTATION

      omegas(:, 1) = f * CST_PIt2

      allocate(Hr_w(NFREQS, NMODES_EFF))
      allocate(Hi_w(NFREQS, NMODES_EFF))


      do imi = 1, NMODES_EFF

         mi = MODES(imi)

         ! H1 / H2
         r_part(:, 1) = - (omegas(:, 1)*omegas(:, 1) * struct_data%modal_%Mm_(mi)) + struct_data%modal_%Km_(mi)
         i_part(:, 1) = omegas(:, 1) * struct_data%modal_%Cm_(mi, mi)
         h_tmp(:, 1)  = r_part(:, 1)*r_part(:, 1) + i_part(:, 1)*i_part(:, 1)

         Hr_w(:, imi) =   r_part(:, 1) / h_tmp(:, 1)
         Hi_w(:, imi) = - i_part(:, 1) / h_tmp(:, 1)

         if (settings%i_compute_bisp_ == 1) then
            ! H12
            do ifrj = 1, NFREQS
               h_tmp(:, ifrj) = omegas(:, 1) + omegas(ifrj, 1)
            enddo
            r_part = - (h_tmp*h_tmp * struct_data%modal_%Mm_(mi)) + struct_data%modal_%Km_(mi)
            i_part = h_tmp * struct_data%modal_%Cm_(mi, mi)
            h_tmp  = r_part*r_part + i_part*i_part

            Hr_w1w2(:, :, imi) =   r_part / h_tmp
            Hi_w1w2(:, :, imi) = - i_part / h_tmp
         endif
      enddo


      ! first deallocation
      if (allocated(r_part)) deallocate(r_part)
      if (allocated(i_part)) deallocate(i_part)


      ! BUG: optimise memory accesses !!!
      do imk = 1, NMODES_EFF

         if (settings%i_compute_bisp_ == 1) then
            h_tmp = Hr_w1w2(:, :, imk)
            h_tmp2= Hi_w1w2(:, :, imk)
         endif 

         do imj = 1, NMODES_EFF

            psd(:, posm_psd) = &
               psd(:, posm_psd) * &
               (  Hr_w(:, imk) * Hr_w(:, imj) + &
                  Hi_w(:, imk) * Hi_w(:, imj) )

            ! BISPs
            if (settings%i_compute_bisp_ == 1) then

               do imi = 1, NMODES_EFF

                  do ifrj = 1, NFREQS

                     bisp(:, ifrj, posm_bisp) = &
                        bisp(:, ifrj, posm_bisp) * &
                        (  Hr_w(:, imi) * Hr_w(ifrj, imj) * h_tmp(:, ifrj) + &
                           Hr_w(:, imi) * Hi_w(ifrj, imj) * h_tmp2(:, ifrj) + &
                           Hi_w(:, imi) * Hr_w(ifrj, imj) * h_tmp2(:, ifrj) - &
                           Hi_w(:, imi) * Hi_w(ifrj, imj) * h_tmp(:, ifrj) )
                  enddo

                  posm_bisp = posm_bisp + 1
               enddo ! i mode
            endif ! i bisp

            posm_psd = posm_psd + 1
         enddo ! j mode
      enddo ! k mode


      ! deallocate what s remained
      if (allocated(Hr_w)) deallocate(Hr_w)
      if (allocated(Hi_w)) deallocate(Hi_w)
      if (allocated(h_tmp)) deallocate(h_tmp)

      if (settings%i_compute_bisp_ == 1) then
         if (allocated(h_tmp2)) deallocate(h_tmp2)
         if (allocated(Hr_w1w2)) deallocate(Hr_w1w2)
         if (allocated(Hi_w1w2)) deallocate(Hi_w1w2)
      endif


#ifdef __BSA_DEBUG
      write(unit_debug_, '(2a)') &
         INFOMSG, '@BsaClassicImpl::getRM_full_vect_cls_() : computing modal responses spectra -- ok.'
#endif
   end subroutine getRM_full_vect_cls_








   module subroutine getFM_diag_tnlm_vect_cls_(f, Suvw, psd, bisp)
      real(RDP), intent(in) :: f(NFREQS)
      real(RDP), intent(in) :: Suvw(NFREQS, NPSDEL)
      real(RDP), intent(inout), allocatable :: psd(:, :), bisp(:, :, :)

      integer(kind = 4) :: iin = 0, ien = 0, itmp = 0, ifrj = 0
      integer(kind = 4) :: i_n_pad = 0, i_pad_len = 0

      ! turb components related
      integer(kind = 4) :: itc   = 0, tc_pos = 0, tc_posN = 0
      integer(kind = 4) :: tc, tcP3

      ! node related indexes
      integer(kind = 4) :: in, n, posN

      ! libs related vars
      integer(kind = 4) :: ilk, ilj, ili, lk, lj, li

      ! modal amtrix related
      real(RDP), dimension(NMODES_EFF) :: phik, phij, phii
      integer(kind = 4) :: posk, posj, posi

      ! wind forces coeffs
      real(RDP) :: ai, aiU, aj, ajU, ak, akU

      ! modes related
      integer(kind = 4) :: im

      ! 
      ! real(RDP) :: S_N_curr(settings%nfreqs_)
      real(RDP), allocatable :: Suvw_N_T(:, :), BF_ijk_III_w1w2(:, :)
      real(RDP), allocatable :: tmp1(:, :), tmp2(:, :), tmp3(:, :)
      real(RDP), allocatable :: Suvw_N_w1w2(:, :), Suvw_N_pad(:)
      real(RDP), allocatable :: PSDF_jk_JJ_w(:)


#ifdef __BSA_DEBUG
      write(unit_debug_, '(2a)') &
         INFOMSG, '@BsaClassicImpl::getFM_diag_tnlm_vect_cls_() : computing modal forces spectra...'
#endif

      ! getting padded length and relative init/end indices (non zero zone)
      itmp      = NFREQS - 1    ! do not consider 0 (point of symmetry)
      i_n_pad   = itmp / 2                ! spread it on the two sides (left / right)
      ien       = i_n_pad + NFREQS
      iin       = i_n_pad + 1
      i_pad_len = itmp + NFREQS


      if (settings%i_compute_psd_ == 1) then

         allocate(PSDF_jk_JJ_w(NFREQS))
         allocate(psd(NFREQS, dimM_psd_))
         psd = 0._RDP
      endif

      if (settings%i_compute_bisp_ == 1) then

         allocate(Suvw_N_T(1, NFREQS))
         allocate(Suvw_N_pad(i_pad_len))
         allocate(Suvw_N_w1w2(NFREQS, NFREQS))

         allocate(tmp1(NFREQS, NFREQS))
         allocate(tmp2(NFREQS, NFREQS))
         allocate(tmp3(NFREQS, NFREQS))

         allocate(BF_ijk_III_w1w2(NFREQS, NFREQS))
         allocate(bisp(NFREQS, NFREQS, dimM_bisp_))
         bisp = 0._RDP
      endif


      !========================================================================
      ! BUG: for the moment, only considering correlation
      !      between same turbulence component (u, v, w).
      !      No cross-correlation between turbulent components
      !      i.e. E[uv]==E[uw]==E[vw] === 0
      do itc = 1, NTCOMPS

         tc      = wd%tc_(itc)
         tcP3    = tc + 3
         tc_posN = (itc - 1) * NNODESL

         do in = 1, NNODESL

            n      = struct_data%n_load_(in)
            posN   = (n - 1) * NLIBS
            tc_pos = tc_posN + in


            ! BISP
            if (settings%i_compute_bisp_ == 1) then

               Suvw_N_pad(iin : ien) = Suvw(:, tc_pos)
               Suvw_N_T(1, :)        = Suvw(:, tc_pos)  ! storing transpose once for all

               itmp = 0
!DIR$ UNROLL= 10
               do ifrj = 1, NFREQS

                  Suvw_N_w1w2(:, ifrj) = Suvw_N_pad(ifrj : NFREQS+itmp)
                  
                  tmp1(:, ifrj) = Suvw_N_w1w2(:, ifrj) * Suvw_N_T(1, ifrj)
                  tmp2(:, ifrj) = Suvw_N_w1w2(:, ifrj) * Suvw(:, tc_pos)
                  ! tmp3(:, ifrj) = Suvw(:, tc_pos) * Suvw_N_T(1, ifrj)

                  itmp = itmp + 1
               enddo ! omega j

               ! original version
               ! NOTE: had stack overflow not using /heap-arrays0 compile option
               tmp3 = matmul(Suvw(:, tc_pos:tc_pos), Suvw_N_T)


               do ilk = 1, NLIBSL

                  lk   = struct_data%libs_load_(ilk)
                  posk = posN + lk
                  phik = struct_data%modal_%phi_(posk, MODES)

                  akU = wd%wfc_(lk, tc,   in)
                  ak  = wd%wfc_(lk, tcP3, in)


                  do ilj = 1, NLIBSL

                     lj   = struct_data%libs_load_(ilj)
                     posj = posN + lj
                     phij = struct_data%modal_%phi_(posj, MODES)

                     ajU = wd%wfc_(lj, tc,   in)
                     aj  = wd%wfc_(lj, tcP3, in)

!DIR$ UNROLL= 6
                     do ili = 1, NLIBSL

                        li   = struct_data%libs_load_(ili)
                        posi = posN + li
                        phii = struct_data%modal_%phi_(posi, MODES)

                        aiU = wd%wfc_(li, tc,   in)
                        ai  = wd%wfc_(li, tcP3, in)


                        BF_ijk_III_w1w2 = 2 * (&
                           ai  * ajU * akU * tmp3 + &
                           aiU * aj  * akU * tmp1 + &
                           aiU * ajU * ak  * tmp2 &
                        &)


                        ! NOTE: no need to retrieve actual mode
                        !       since we are supposed to store only
                        !       kept mode for all related variables.
!DIR$ UNROLL= 10
                        do im = 1, NMODES_EFF

                           bisp(:, :, im) = bisp(:, :, im) + &
                              phik(im) * phij(im) * phii(im) * BF_ijk_III_w1w2
                        enddo

                     enddo ! i lib
                  enddo ! j lib
               enddo ! k lib

            endif ! do bisp computation


            ! PSDs
            do ilk = 1, NLIBSL

               lk   = struct_data%libs_load_(ilk)
               posk = posN + lk
               phik = struct_data%modal_%phi_(posk, MODES)

               akU = wd%wfc_(lk, tc, in)

!DIR$ UNROLL= 6
               do ilj = 1, NLIBSL

                  lj   = struct_data%libs_load_(ilj)
                  posj = posN + lj
                  phij = struct_data%modal_%phi_(posj, MODES)

                  ajU = wd%wfc_(lj, tc, in)


                  PSDF_jk_JJ_w = ajU * akU * Suvw(:, tc_pos)

!DIR$ UNROLL= 10
                  do im = 1, NMODES_EFF

                     psd(:, im) = psd(:, im) + &
                        phij(im) * phik(im) * PSDF_jk_JJ_w
                  enddo
               enddo ! j lib
            enddo ! k lib


#ifdef __BSA_DEBUG
            print '(1x, 2a, f10.4, " %")', &
               INFOMSG, 'getFM_diag_tnlm_vect_cls_() :   done  ', real(in, RDP) / NNODESL * 100
#endif

         enddo ! nodes
      enddo ! turb comps

      if (allocated(PSDF_jk_JJ_w))    deallocate(PSDF_jk_JJ_w)
      if (allocated(Suvw_N_T))        deallocate(Suvw_N_T)
      if (allocated(Suvw_N_pad))      deallocate(Suvw_N_pad)
      if (allocated(Suvw_N_w1w2))     deallocate(Suvw_N_w1w2)
      if (allocated(tmp1))            deallocate(tmp1)
      if (allocated(tmp2))            deallocate(tmp2)
      if (allocated(tmp3))            deallocate(tmp3)
      if (allocated(BF_ijk_III_w1w2)) deallocate(BF_ijk_III_w1w2)

#ifdef __BSA_DEBUG
      write(unit_debug_, '(2a)') &
         INFOMSG, '@BsaClassicImpl::getFM_diag_tnlm_vect_cls_() : computing modal forces spectra -- ok.'
#endif
   end subroutine



   module subroutine getRM_diag_vect_cls_(f, psd, bisp)
      real(RDP), intent(in)                 :: f(NFREQS)
      real(RDP), allocatable, intent(inout) :: psd(:, :), bisp(:, :, :)

      integer(kind = 4) :: ifrj
      integer(kind = 4) :: pos = 1
      
      ! modal indexed
      integer(kind = 4) :: imi, mi

      real(RDP) :: omegas(NFREQS, 1)
      real(RDP), allocatable :: r_part(:, :), i_part(:, :), h_tmp(:, :)
      real(RDP) :: Mgi, Kgi, Cgi
      real(RDP), allocatable :: Hr_w(:), Hi_w(:)
      real(RDP), allocatable :: Hr_w1w2(:, :), Hi_w1w2(:, :)


#ifdef __BSA_DEBUG
      write(unit_debug_, '(2a)') &
         INFOMSG, '@BsaClassicImpl::getRM_diag_vect_cls_() : computing modal responses spectra...'
#endif


      ! BUG: check logic if correct
      if (settings%i_compute_bisp_ == 1) then

         allocate(r_part(NFREQS, NFREQS))
         allocate(i_part(NFREQS, NFREQS))
         allocate(h_tmp(NFREQS, NFREQS))
         allocate(Hr_w1w2(NFREQS, NFREQS))
         allocate(Hi_w1w2(NFREQS, NFREQS))

      elseif (settings%i_compute_psd_ == 1) then

         allocate(r_part(NFREQS, 1))
         allocate(i_part(NFREQS, 1))
         allocate(h_tmp(NFREQS, 1))

      else ! none of them. NOTE: if so, we don't ever get here..
         return
      endif

      
      ! TRANSFER FUNCTION COMPUTATION

      ! pulsations
      omegas(:, 1) = f * CST_PIt2


      allocate(Hr_w(NFREQS))
      allocate(Hi_w(NFREQS))


      do imi = 1, NMODES_EFF

         mi = MODES(imi)

         ! NOTE: assumes we are storing only kept modes
         !       related modal info
         Mgi = struct_data%modal_%Mm_(mi)
         Kgi = struct_data%modal_%Km_(mi)
         Cgi = struct_data%modal_%Cm_(mi, mi)

         ! H1 / H2
         r_part(:, 1) = - (omegas(:, 1)*omegas(:, 1) * Mgi) + Kgi
         i_part(:, 1) = omegas(:, 1) * Cgi
         h_tmp(:, 1)  = r_part(:, 1)*r_part(:, 1) + i_part(:, 1)*i_part(:, 1)

         Hr_w =   r_part(:, 1) / h_tmp(:, 1)
         Hi_w = - i_part(:, 1) / h_tmp(:, 1)


         ! PSD
         psd(:, pos) = &
            psd(:, pos) * &
            (  Hr_w * Hr_w + &
               Hi_w * Hi_w )


         if (settings%i_compute_bisp_ == 1) then

            ! w1+w2^T
            do ifrj = 1, NFREQS
               h_tmp(:, ifrj) = omegas(:, 1) + omegas(ifrj, 1)
            enddo

            r_part = - (h_tmp*h_tmp * Mgi) + Kgi
            i_part = h_tmp * Cgi
            h_tmp  = r_part*r_part + i_part*i_part

            Hr_w1w2 =   r_part / h_tmp
            Hi_w1w2 = - i_part / h_tmp

            ! BISP
            do ifrj = 1, NFREQS

               bisp(:, ifrj, pos) = &
                  bisp(:, ifrj, pos) * &
                  (  Hr_w * Hr_w(ifrj) * Hr_w1w2(:, ifrj) + &
                     Hr_w * Hi_w(ifrj) * Hi_w1w2(:, ifrj) + &
                     Hi_w * Hr_w(ifrj) * Hi_w1w2(:, ifrj) - &
                     Hi_w * Hi_w(ifrj) * Hr_w1w2(:, ifrj) )
            enddo
         endif

         pos = pos + 1
      enddo ! i mode


      ! deallocate
      if (allocated(r_part)) deallocate(r_part)
      if (allocated(i_part)) deallocate(i_part)
      if (allocated(Hr_w))   deallocate(Hr_w)
      if (allocated(Hi_w))   deallocate(Hi_w)
      if (allocated(h_tmp))  deallocate(h_tmp)

      if (settings%i_compute_bisp_ == 1) then
         if (allocated(Hr_w1w2)) deallocate(Hr_w1w2)
         if (allocated(Hi_w1w2)) deallocate(Hi_w1w2)
      endif


#ifdef __BSA_DEBUG
      write(unit_debug_, '(2a)') &
         INFOMSG, '@BsaClassicImpl::getRM_diag_vect_cls_() : computing modal responses spectra -- ok.'
#endif
   end subroutine getRM_diag_vect_cls_












   !> BUG: this routine is adapted to the case where we use
   !>      convention on PULSATION.
   !>      Please, adapt it to the case of convention over FREQUENCIES.
   module pure subroutine getFM_full_tnlm_scalar_cls_(ii, ij, fi, fj, Suvw, Suvw_pad, psd, bisp)
      integer, intent(in)      :: ii, ij
      real(RDP), intent(in)    :: fi, fj
      real(RDP), intent(in)    :: Suvw(NFREQS, NPSDEL)
      real(RDP), intent(in)    :: Suvw_pad(NPSDEL)
      real(RDP), intent(inout) :: psd(dimM_psd_), bisp(dimM_bisp_)

      ! turb components related
      integer(kind = 4) :: itc, tc_posN, tc_pk, tc_pj, tc_pi

      ! freqs related
      real(RDP) :: fiabs, fjabs, fiPfj, fiPfjabs

      ! nodes indexed values
      integer(kind = 4) :: i_pos_nk, i_pos_nj, i_pos_ni
      integer(kind = 4) :: pos_nk, pos_nj, pos_ni
      integer(kind = 4) :: ink, inj, ini
      integer(kind = 4) :: ni, nj, nk

      ! libs indexed values
      integer(kind = 4) :: ilk !, ilj, ili

      ! modes indexed values
      real(RDP), dimension(NLIBSL, NMODES_EFF) :: phik_, phij_, phii_
      real(RDP) :: phik(NMODES_EFF), phikk, phij(1, NLIBSL), phii(NLIBSL, 1)
      integer   :: posm
      integer   :: imk, imj, imi 
      ! integer(kind = 4) :: mi, mj, mk

      ! local nodal correlations
      real(RDP) :: corrJK, corrIK, corrIJ

      ! wfc extractions
      integer(kind = 4) :: tc, tcP3
      real(RDP), dimension(NLIBSL, 1) :: aiU, ai, akU, ak
      real(RDP), dimension(1, NLIBSL) :: ajU, aj

      ! PSDs local
      real(RDP) :: S_uvw_i_i, S_uvw_i_j, S_uvw_i_ij
      real(RDP) :: S_uvw_j_i, S_uvw_j_j, S_uvw_j_ij
      real(RDP) :: S_uvw_k_i, S_uvw_k_j, S_uvw_k_ij
      real(RDP) :: S_uvw_JK_i, S_uvw_JK_j
      real(RDP) :: S_uvw_IK_j, S_uvw_IK_ij
      real(RDP) :: S_uvw_IJ_i, S_uvw_IJ_ij

      ! BF local
      real(RDP), dimension(NLIBSL, NLIBSL) :: BF_ijk_IJK_w_w2, PSD_jk_JK_w
      real(RDP), dimension(NLIBSL, NLIBSL) :: tmp1, tmp2, tmp3
      !========================================================================

      psd = 0._RDP
      bisp= 0._RDP

      fiabs = abs(fi)
      fjabs = abs(fj)

      fiPfj    = fi + fj
      fiPfjabs = abs(fiPfj)


      !========================================================================
      ! BUG: for the moment, only considering correlation
      !      between same turbulence component (u, v, w).
      !      No cross-correlation between turbulent components
      !      i.e. E[uv]==E[uw]==E[vw] === 0
      do itc = 1, wd%i_ntc_

         tc_posN = (itc - 1) * NNODESL
         tc      = wd%tc_(itc) ! get actual turbulent component
         tcP3    = tc + 3      ! quadratic term coeff


         i_pos_nk = 1
         do ink = 1, NNODESL

            nk     = struct_data%n_load_(ink)
            pos_nk = (nk - 1) * NLIBS
            tc_pk  = tc_posN + i_pos_nk

            phik_  = struct_data%modal_%phi_(pos_nk + struct_data%libs_load_, MODES)

            S_uvw_k_i  = Suvw(ii, tc_pk)
            S_uvw_k_j  = Suvw(ij, tc_pk)
            S_uvw_k_ij = Suvw_pad(tc_pk)


            akU(:, 1) = wd%wfc_(struct_data%libs_load_, tc,   ink)
            ak (:, 1) = wd%wfc_(struct_data%libs_load_, tcP3, ink)

            
            i_pos_nj = 1
            do inj = 1, NNODESL

               nj     = struct_data%n_load_(inj)
               pos_nj = (nj - 1) * NLIBS
               tc_pj  = tc_posN + i_pos_nj

               phij_  = struct_data%modal_%phi_(pos_nj + struct_data%libs_load_, MODES)

               corrJK = wd%nod_corr_(util_getCorrVectIndex(nj, nk, NNODES), tc)

               S_uvw_j_i  = Suvw(ii, tc_pj)
               S_uvw_j_j  = Suvw(ij, tc_pj)
               S_uvw_j_ij = Suvw_pad(tc_pj)


               S_uvw_JK_i = corrJK**(fiabs) * sqrt(S_uvw_j_i * S_uvw_k_i)
               S_uvw_JK_j = corrJK**(fjabs) * sqrt(S_uvw_j_j * S_uvw_k_j)


               ajU(1, :) = wd%wfc_(struct_data%libs_load_, tc,   inj)
               aj (1, :) = wd%wfc_(struct_data%libs_load_, tcP3, inj)



               !! BISPs

               if (settings%i_compute_bisp_ == 1) then

                  i_pos_ni = 1
                  do ini = 1, NNODESL

                     ni     = struct_data%n_load_(ini)
                     pos_ni = (ni - 1) * NLIBS
                     tc_pi  = tc_posN + i_pos_ni

                     phii_  = struct_data%modal_%phi_(pos_ni + struct_data%libs_load_, MODES)
      
                     corrIK = wd%nod_corr_(util_getCorrVectIndex(ni, nk, NNODES), tc)
                     corrIJ = wd%nod_corr_(util_getCorrVectIndex(ni, nj, NNODES), tc)

                     S_uvw_i_i  = Suvw(ii, tc_pi)
                     S_uvw_i_j  = Suvw(ij, tc_pi)
                     S_uvw_i_ij = Suvw_pad(tc_pi)


                     S_uvw_IK_j  = corrIK**(fjabs)    * sqrt(S_uvw_i_j  * S_uvw_k_j )
                     S_uvw_IK_ij = corrIK**(fiPfjabs) * sqrt(S_uvw_i_ij * S_uvw_k_ij)


                     S_uvw_IJ_i  = corrIJ**(fiabs)    * sqrt(S_uvw_i_i  * S_uvw_j_i )
                     S_uvw_IJ_ij = corrIJ**(fiPfjabs) * sqrt(S_uvw_i_ij * S_uvw_j_ij)


                     aiU(:, 1) = wd%wfc_(struct_data%libs_load_, tc,   ini)
                     ai (:, 1) = wd%wfc_(struct_data%libs_load_, tcP3, ini)


                     tmp1 = matmul(ai, ajU ) * S_uvw_IJ_i  * S_uvw_IK_j
                     tmp2 = matmul(aiU, aj ) * S_uvw_IJ_ij * S_uvw_JK_j
                     tmp3 = matmul(aiU, ajU) * S_uvw_JK_i  * S_uvw_IK_ij


                     do ilk = 1, NLIBSL

                        phik = phik_(ilk, :)

                              
                        BF_ijk_IJK_w_w2 = 2 * (&
                           tmp1 * akU(ilk, 1) + &
                           tmp2 * akU(ilk, 1) + &
                           tmp3 * ak (ilk, 1)   &
                        &)


                        ! if (all(BF_ijk_IJK_w_w2 == 0._RDP)) cycle


                        posm = 1
                        do imk = 1, NMODES_EFF

                           phikk = phik(imk)

                           do imj = 1, NMODES_EFF

                              phij(1, :) = phij_(:, imj)

                              do imi = 1, NMODES_EFF

                                 phii(:, 1) = phii_(:, imi)

                                 bisp(posm) = bisp(posm) + &
                                    sum(matmul(phii, phij) * BF_ijk_IJK_w_w2 * phikk)

                                 posm = posm + 1
                              enddo ! i mode
                           enddo ! j mode
                        enddo ! k mode                     


                     enddo ! k lib


                     i_pos_ni = i_pos_ni + 1
                  enddo ! i node

               endif ! bisp computation


               !! PSDs
               if (ij == 1) then


                  ! PSD f
                  PSD_jk_JK_w = matmul(akU, ajU) * S_uvw_JK_i


                  posm = 1
                  do imk = 1, NMODES_EFF

                     ! NOTE: reusing variable, but naming is wrong!!!
                     phii = phik_(:, imk:imk)

                     do imj = 1, NMODES_EFF

                        phij(1, :) = phij_(:, imj)

                        psd(posm) = psd(posm) + &
                           sum(matmul(phii, phij) * PSD_jk_JK_w)
!                                     phik, phij

                        posm = posm + 1
                     enddo ! j mode
                  enddo ! k mode

               endif ! PSD computation


               i_pos_nj = i_pos_nj + 1
            enddo ! j node

            i_pos_nk = i_pos_nk + 1
         enddo ! k node

      enddo ! itc

   end subroutine getFM_full_tnlm_scalar_cls_











   !> BUG: this routine is adapted to the case where we use
   !>      convention on PULSATION.
   !>      Please, adapt it to the case of convention over FREQUENCIES.
   module pure subroutine getFM_full_tnm_scalar_cls_(ii, ij, fi, fj, Suvw, Suvw_pad, psd, bisp)
      integer, intent(in)      :: ii, ij
      real(RDP), intent(in)    :: fi, fj
      real(RDP), intent(in)    :: Suvw(NFREQS, NPSDEL)
      real(RDP), intent(in)    :: Suvw_pad(NPSDEL)
      real(RDP), intent(inout) :: psd(dimM_psd_), bisp(dimM_bisp_)

      ! turb components related
      integer(kind = 4) :: itc, tc_posN, tc_pk, tc_pj

      ! freqs related
      real(RDP) :: fiabs, fjabs, fiPfj, fiPfjabs

      ! nodes indexed values
      integer(kind = 4) :: i_pos_nk, i_pos_nj
      integer(kind = 4) :: pos_nk, pos_nj
      integer(kind = 4) :: ink, inj
      integer(kind = 4) :: nj, nk

      ! libs indexed values
      integer(kind = 4) :: ilk !, ilj, ili

      ! modes indexed values
      real(RDP), dimension(NMODES_EFF, 2) :: phik_, phij_
      real(RDP) :: phij_Ub_, phij_u_, phik_Ub_, phik_u_
      integer   :: posm
      integer   :: imk, imj, imi

      ! local nodal correlations
      real(RDP) :: corrJK

      ! wfc extractions
      integer(kind = 4) :: tc, tcP3

      ! PSDs local
      real(RDP) :: S_uvw_j_i,  S_uvw_j_j,  S_uvw_j_ij
      real(RDP) :: S_uvw_k_i,  S_uvw_k_j,  S_uvw_k_ij
      real(RDP) :: S_uvw_JK_i, S_uvw_JK_j
      !========================================================================

      psd  = 0._RDP
      bisp = 0._RDP

      fiabs = abs(fi)
      fjabs = abs(fj)

      fiPfj    = fi + fj
      fiPfjabs = abs(fiPfj)


      !========================================================================
      ! BUG: for the moment, only considering correlation
      !      between same turbulence component (u, v, w).
      !      No cross-correlation between turbulent components
      !      i.e. E[uv]==E[uw]==E[vw] === 0
      do itc = 1, wd%i_ntc_

         tc_posN = (itc - 1) * NNODESL
         tc      = wd%tc_(itc) ! get actual turbulent component
         tcP3    = tc + 3      ! quadratic term coeff


         i_pos_nk = 1
         do ink = 1, NNODESL

            nk     = struct_data%n_load_(ink)
            pos_nk = (nk - 1) * NLIBS
            tc_pk  = tc_posN + i_pos_nk

            phik_(:, 1) = wd%phi_times_A_ndegw_(:, ink, tc  )
            phik_(:, 2) = wd%phi_times_A_ndegw_(:, ink, tcP3)

            S_uvw_k_i  = Suvw(ii, tc_pk)
            S_uvw_k_j  = Suvw(ij, tc_pk)
            S_uvw_k_ij = Suvw_pad(tc_pk)

            
            i_pos_nj = 1
            do inj = 1, NNODESL

               nj     = struct_data%n_load_(inj)
               pos_nj = (nj - 1) * NLIBS
               tc_pj  = tc_posN + i_pos_nj

               phij_(:, 1) = wd%phi_times_A_ndegw_(:, inj, tc  )
               phij_(:, 2) = wd%phi_times_A_ndegw_(:, inj, tcP3)

               corrJK = wd%nod_corr_(util_getCorrVectIndex(nj, nk, NNODES), tc)

               S_uvw_j_i  = Suvw(ii, tc_pj)
               S_uvw_j_j  = Suvw(ij, tc_pj)
               S_uvw_j_ij = Suvw_pad(tc_pj)


               S_uvw_JK_i = corrJK**(fiabs) * sqrt(S_uvw_j_i * S_uvw_k_i)
               S_uvw_JK_j = corrJK**(fjabs) * sqrt(S_uvw_j_j * S_uvw_k_j)


               !! BISPs

               if (settings%i_compute_bisp_ == 1) then

                  block
                     integer(kind = 4) :: i_pos_ni, ini, ni, pos_ni
                     integer(kind = 4) :: tc_pi

                     real(RDP) :: corrIK, corrIJ
                     real(RDP) :: S_uvw_i_i,  S_uvw_i_j,  S_uvw_i_ij
                     real(RDP) :: S_uvw_IK_j, S_uvw_IK_ij
                     real(RDP) :: S_uvw_IJ_i, S_uvw_IJ_ij
                     real(RDP) :: tmp1, tmp2, tmp3

                     real(RDP), dimension(NMODES_EFF, 2) :: phii_


                     i_pos_ni = 1
                     do ini = 1, NNODESL

                        ni     = struct_data%n_load_(ini)
                        pos_ni = (ni - 1) * NLIBS
                        tc_pi  = tc_posN + i_pos_ni

                        phii_(:, 1) = wd%phi_times_A_ndegw_(:, ini, tc  )
                        phii_(:, 2) = wd%phi_times_A_ndegw_(:, ini, tcP3)
         
                        corrIK = wd%nod_corr_(util_getCorrVectIndex(ni, nk, NNODES), tc)
                        corrIJ = wd%nod_corr_(util_getCorrVectIndex(ni, nj, NNODES), tc)

                        S_uvw_i_i  = Suvw(ii, tc_pi)
                        S_uvw_i_j  = Suvw(ij, tc_pi)
                        S_uvw_i_ij = Suvw_pad(tc_pi)


                        S_uvw_IK_j  = corrIK**(fjabs)    * sqrt(S_uvw_i_j  * S_uvw_k_j )
                        S_uvw_IK_ij = corrIK**(fiPfjabs) * sqrt(S_uvw_i_ij * S_uvw_k_ij)

                        S_uvw_IJ_i  = corrIJ**(fiabs)    * sqrt(S_uvw_i_i  * S_uvw_j_i )
                        S_uvw_IJ_ij = corrIJ**(fiPfjabs) * sqrt(S_uvw_i_ij * S_uvw_j_ij)


                        tmp1 = S_uvw_IJ_i  * S_uvw_IK_j
                        tmp2 = S_uvw_IJ_ij * S_uvw_JK_j
                        tmp3 = S_uvw_JK_i  * S_uvw_IK_ij


                        posm = 1
                        do imk = 1, NMODES_EFF

                           phik_Ub_ = phik_(imk, 1)
                           phik_u_  = phik_(imk, 2)

                           do imj = 1, NMODES_EFF

                              phij_Ub_ = phij_(imj, 1)
                              phij_u_  = phij_(imj, 2)

                              do imi = 1, NMODES_EFF

                                 bisp(posm) = bisp(posm) +  &
                                    2 * ( &
                                          (phii_(imi, 2) * phij_Ub_ * phik_Ub_ * tmp1)  &
                                       +  (phii_(imi, 1) * phij_u_  * phik_Ub_ * tmp2)  &
                                       +  (phii_(imi, 1) * phij_Ub_ * phik_u_  * tmp3)  &
                                    )

                                 posm = posm + 1
                              enddo ! i mode
                           enddo ! j mode
                        enddo ! k mode


                        i_pos_ni = i_pos_ni + 1
                     enddo ! i node

                  end block

               endif ! bisp computation


               !! PSDs
               if (ij == 1) then

                  posm = 1
                  do imk = 1, NMODES_EFF

                     phik_Ub_ = phik_(imk, 1)

                     do imj = 1, NMODES_EFF

                        psd(posm) = psd(posm) + &
                           (phik_Ub_ * phij_(imj, 1) * S_uvw_JK_i)

                        posm = posm + 1
                     enddo ! j mode
                  enddo ! k mode

               endif ! PSD computation


               i_pos_nj = i_pos_nj + 1
            enddo ! j node

            i_pos_nk = i_pos_nk + 1
         enddo ! k node

      enddo ! itc

   end subroutine getFM_full_tnm_scalar_cls_










   module subroutine getRM_full_scalar_cls_(ii, ij, fi, fj, psdin, psdout, bispin, bispout)
      integer, intent(in)      :: ii, ij
      real(RDP), intent(in)    :: fi, fj
      real(RDP), intent(in)    :: psdin(dimM_psd_), bispin(dimM_bisp_)
      real(RDP), intent(out)   :: psdout(dimM_psd_), bispout(dimM_bisp_)

      real(RDP) :: wi
      integer(kind = 4) :: posm, imk, imj

      real(RDP), dimension(NMODES_EFF) :: Cdiag, rpart, ipart, htmp
      real(RDP), dimension(NMODES_EFF) :: H1r, H1i

      real(RDP) :: H2j_r, H2j_i

      wi = fi * CST_PIt2


      ! pre evaluate TFs (per mode)

      ! H1
      rpart = - (wi*wi * struct_data%modal_%Mm_(MODES)) + struct_data%modal_%Km_(MODES)
      do imk = 1, NMODES_EFF
         Cdiag(imk) = struct_data%modal_%Cm_(MODES(imk), MODES(imk))
      enddo
      ipart = Cdiag * wi
      htmp  = rpart*rpart + ipart*ipart
      H1r   =   rpart / htmp
      H1i   = - ipart / htmp


      if (ij == 1) then ! also PSDr

         posm = 1
         do imk = 1, NMODES_EFF

            H2j_r = H1r(imk)
            H2j_i = H1i(imk)

            do imj = 1, NMODES_EFF

               psdout(posm) = psdin(posm) * (&
                  H1r(imj) * H2j_r + &
                  H1i(imj) * H2j_i   &
               )

               posm = posm + 1
            enddo ! j modes
         enddo ! k modes

         if (settings%i_compute_bisp_ == 0) return
      endif


      ! BISPr

      block
         real(RDP), dimension(NMODES_EFF) :: H2r, H2i
         real(RDP), dimension(NMODES_EFF) :: H12r, H12i
         real(RDP) :: H12k_r, H12k_i
         real(RDP) :: wj, wiPwj, tmp1, tmp2, tmp3, tmp4
         integer   :: imi

         wj    = fj * CST_PIt2
         wiPwj = wi + wj

         ! H2
         rpart = - (wj*wj * struct_data%modal_%Mm_(MODES)) + struct_data%modal_%Km_(MODES)
         ipart = Cdiag * wj
         htmp  = rpart*rpart + ipart*ipart
         H2r   =   rpart / htmp
         H2i   = - ipart / htmp


         ! H12
         rpart  = - (wiPwj*wiPwj * struct_data%modal_%Mm_(MODES)) + struct_data%modal_%Km_(MODES)
         ipart  = Cdiag * wiPwj
         htmp   = rpart*rpart + ipart*ipart
         H12r   =   rpart / htmp
         H12i   = - ipart / htmp


         posm = 1
         do imk = 1, NMODES_EFF

            H12k_r = H12r(imk)
            H12k_i = H12i(imk)

            do imj = 1, NMODES_EFF

               H2j_r = H2r(imj)
               H2j_i = H2i(imj)

               tmp1 = H2j_r * H12k_r  ! real
               tmp2 = H2j_i * H12k_i  ! real
               tmp3 = H2j_r * H12k_i  ! imag
               tmp4 = H2j_i * H12k_r  ! imag

               do imi = 1, NMODES_EFF

                  bispout(posm) = bispin(posm) * (&
                     H1r(imi) * tmp1 + &
                     H1r(imi) * tmp2 + &
                     H1i(imi) * tmp3 - &
                     H1i(imi) * tmp4   &
                  )

                  posm = posm + 1
               enddo ! i modes
            enddo ! j modes
         enddo ! k modes

      end block

   end subroutine getRM_full_scalar_cls_














   !> BUG: this routine is adapted to the case where we use
   !>      convention on PULSATION.
   !>      Please, adapt it to the case of convention over FREQUENCIES.
   module pure subroutine getFM_diag_tnlm_scalar_cls_(ii, ij, fi, fj, Suvw, Suvw_pad, psd, bisp)
      integer, intent(in)      :: ii, ij
      real(RDP), intent(in)    :: fi, fj
      real(RDP), intent(in)    :: Suvw(NFREQS, NPSDEL)
      real(RDP), intent(in)    :: Suvw_pad(NPSDEL)
      real(RDP), intent(inout) :: psd(dimM_psd_), bisp(dimM_bisp_)

      integer(kind = 4) :: itc, tc, tcP3, iposN, inode, n, imode, posNi

      real(RDP) :: Suvw_i

      real(RDP), dimension(NLIBSL, NMODES_EFF) :: phiN_
      real(RDP), dimension(NLIBSL, NLIBSL) :: PSDF_jk_JJ_w, tmp3, phiIJ_
      real(RDP), dimension(NLIBSL) :: aNU, aN

      ! logical :: lflag = .false.


      psd  = 0._RDP
      bisp = 0._RDP


      iposN = 1
      do itc = 1, NTCOMPS

         tc   = wd%tc_(itc)
         tcP3 = tc + 3


         do inode = 1, NNODESL

            n     = struct_data%n_load_(inode)
            posNi = (n - 1) * NLIBS

            phiN_ = struct_data%modal_%phi_(posNi + struct_data%libs_load_, MODES)
            
            aNU   = wd%wfc_(struct_data%libs_load_, tc,   inode)
            an    = wd%wfc_(struct_data%libs_load_, tcP3, inode)

            Suvw_i = Suvw(ii, iposN)

            ! PSD comp
            if (ij == 1) then

               do imode = 1, NLIBSL
                  tmp3(:, imode) = aNU * aNU(imode)
               enddo
               PSDF_jk_JJ_w = tmp3 * Suvw_i

               do imode = 1, NMODES_EFF

                  phiIJ_ = matmul(phiN_(:, imode:imode), transpose(phiN_(:, imode:imode)))

                  psd(imode) = psd(imode) + &
                     sum(phiIJ_ * PSDF_jk_JJ_w)

               enddo ! n modes

            endif ! ij == 1



            if (settings%i_compute_bisp_ == 1) then

               block
                  real(RDP) :: Suvw_j, Suvw_ij
                  integer   :: ilk

                  real(RDP) :: BF_ijk_III_wiwj(NLIBSL, NLIBSL)
                  real(RDP), dimension(NLIBSL, NLIBSL) :: tmp2, tmp1

                  Suvw_j  = Suvw(ij, iposN)
                  Suvw_ij = Suvw_pad(iposN)

                  tmp1 = matmul(reshape(aN , [NLIBSL, 1]), reshape(aNU, [1, NLIBSL]))
                  tmp2 = matmul(reshape(aNU, [NLIBSL, 1]), reshape(aN , [1, NLIBSL]))
                  tmp3 = matmul(reshape(aNU, [NLIBSL, 1]), reshape(aNU, [1, NLIBSL]))

                  tmp1 = tmp1 * (Suvw_i  * Suvw_j)
                  tmp2 = tmp2 * (Suvw_ij * Suvw_j)
                  tmp3 = tmp3 * (Suvw_ij * Suvw_i)

                  do ilk = 1, NLIBSL

                     BF_ijk_III_wiwj = 2 * (&
                        tmp1 * aNU(ilk) + &
                        tmp2 * aNU(ilk) + &
                        tmp3 * aN (ilk)   &
                     )

                     do imode = 1, NMODES_EFF

                        phiIJ_ = matmul(phiN_(:, imode:imode), transpose(phiN_(:, imode:imode)))

                        bisp(imode) = bisp(imode) + &
                           sum(phiIJ_ * BF_ijk_III_wiwj * phiN_(ilk, imode))
                     enddo ! n modes

                  enddo ! k libs

               end block

            endif ! compute BISP

            iposN = iposN + 1
         enddo ! nodes loaded

      enddo ! turb comps

   end subroutine getFM_diag_tnlm_scalar_cls_






   module subroutine getRM_diag_scalar_cls_(ii, ij, fi, fj, psdin, psdout, bispin, bispout)
      integer, intent(in)      :: ii, ij  ! freqs indexes
      real(RDP), intent(in)    :: fi, fj
      real(RDP), intent(in)    :: psdin(dimM_psd_), bispin(dimM_bisp_)
      real(RDP), intent(out)   :: psdout(dimM_psd_), bispout(dimM_bisp_)

      real(RDP) :: wi
      integer(kind = 4)   :: imi

      real(RDP), dimension(NMODES_EFF) :: Cdiag, rpart, ipart, htmp
      real(RDP), dimension(NMODES_EFF) :: H1r, H1i


      wi = fi * CST_PIt2


      ! pre evaluate TFs (per mode)

      ! H1
      rpart = - (wi*wi * struct_data%modal_%Mm_(MODES)) + struct_data%modal_%Km_(MODES)
      do imi = 1, NMODES_EFF
         Cdiag(imi) = struct_data%modal_%Cm_(MODES(imi), MODES(imi))
      enddo
      ipart = Cdiag * wi
      htmp  = rpart*rpart + ipart*ipart
      H1r   =   rpart / htmp
      H1i   = - ipart / htmp

      
      if (ij == 1) psdout = psdin * (H1r * H1r + H1i * H1i)


      block
         real(RDP) :: wj, wiPwj
         real(RDP), dimension(NMODES_EFF) :: H2r, H2i
         real(RDP), dimension(NMODES_EFF) :: H12r, H12i
         real(RDP) :: H12k_r, H12k_i, H2j_r, H2j_i
         

         wj    = fj * CST_PIt2
         wiPwj = wi + wj

         ! H2
         rpart = - (wj*wj * struct_data%modal_%Mm_(MODES)) + struct_data%modal_%Km_(MODES)
         ipart = Cdiag * wj
         htmp  = rpart*rpart + ipart*ipart
         H2r   =   rpart / htmp
         H2i   = - ipart / htmp


         ! H12
         rpart = - (wiPwj*wiPwj * struct_data%modal_%Mm_(MODES)) + struct_data%modal_%Km_(MODES)
         ipart = Cdiag * wiPwj
         htmp  = rpart*rpart + ipart*ipart
         H12r  =   rpart / htmp
         H12i  = - ipart / htmp


         bispout = bispin * (&
            H1r * H2r * H12r + &
            H1r * H2i * H12i + &
            H1i * H2r * H12i - &
            H1i * H2i * H12r   &
         )

      end block
      
   end subroutine getRM_diag_scalar_cls_












   module pure subroutine getBR_SFm_val_(nm, Suvw, fnat, im, m, psd)
      !! BUG: very unoptimised..
      !!      is basically a small copy of "getFM_full_tnlm_scalar_cls_"
      integer, intent(in)      :: im, m, nm
      real(RDP), intent(in)    :: Suvw(nm, NPSDEL), fnat
      real(RDP), intent(inout) :: psd

      ! turb components related
      integer(kind = 4) :: itc, tc_posN, tc_pk, tc_pj

      ! freqs related
      real(RDP) :: fiabs

      ! nodes indexed values
      integer(kind = 4) :: i_pos_nk, i_pos_nj
      integer(kind = 4) :: pos_nk, pos_nj
      integer(kind = 4) :: ink, inj
      integer(kind = 4) :: nj, nk, ilk

      ! local nodal correlations
      real(RDP) :: corrJK
      
      integer(kind = 4) :: tc, imk, imj
      real(RDP), dimension(NLIBSL, 1) :: akU, phij_
      real(RDP), dimension(1, NLIBSL) :: ajU, phik_

      ! PSDs local
      real(RDP) :: S_uvw_j_i
      real(RDP) :: S_uvw_k_i
      real(RDP) :: S_uvw_JK_i

      real(RDP), dimension(NLIBSL, NLIBSL) :: PSD_jk_JK_w
      !========================================================================

      psd   = 0._RDP
      fiabs = abs(fnat)


      !========================================================================
      ! BUG: for the moment, only considering correlation
      !      between same turbulence component (u, v, w).
      !      No cross-correlation between turbulent components
      !      i.e. E[uv]==E[uw]==E[vw] === 0
      do itc = 1, wd%i_ntc_

         tc_posN = (itc - 1) * NNODESL
         tc      = wd%tc_(itc) ! get actual turbulent component

         i_pos_nk = 1
         do ink = 1, NNODESL

            nk     = struct_data%n_load_(ink)
            pos_nk = (nk - 1) * NLIBS
            tc_pk  = tc_posN + i_pos_nk

            phik_  = transpose(struct_data%modal_%phi_(pos_nk + struct_data%libs_load_, m:m))

            akU(:, 1) = wd%wfc_(struct_data%libs_load_, tc,   ink)
            
            S_uvw_k_i = Suvw(im, tc_pk)


            
            i_pos_nj = 1
            do inj = 1, NNODESL

               nj     = struct_data%n_load_(inj)
               pos_nj = (nj - 1) * NLIBS
               tc_pj  = tc_posN + i_pos_nj

               phij_  = struct_data%modal_%phi_(pos_nj + struct_data%libs_load_, m:m)

               corrJK = wd%nod_corr_(util_getCorrVectIndex(nj, nk, NNODES), tc)

               S_uvw_j_i  = Suvw(im, tc_pj)
               S_uvw_JK_i = corrJK**(fiabs) * sqrt(S_uvw_j_i * S_uvw_k_i)


               ajU(1, :) = wd%wfc_(struct_data%libs_load_, tc,   inj)


               ! PSD f
               PSD_jk_JK_w = matmul(akU, ajU) * S_uvw_JK_i

               psd = psd + sum(matmul(phij_, phik_) * PSD_jk_JK_w)

               i_pos_nj = i_pos_nj + 1
            enddo ! j node

            i_pos_nk = i_pos_nk + 1
         enddo ! k node

      enddo ! itc
   end subroutine












end submodule BsaLib_FunctionsImpl