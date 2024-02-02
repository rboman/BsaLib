/**
 * This file is part of BSA Library.
 * Copyright (C) 2023  Michele Esposito Marzino 
 *
 * BSA Library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BSA Library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with BSA Library.  If not, see <https://www.gnu.org/licenses/>.
 * */
#ifdef BSACL_USE_CUDA__
#  include "_base.h"
#else
#  ifdef BSACL_PI
#    undef  BSACL_PI
#    define BSACL_PI (BSACL_REAL)M_PI
#  endif
#  ifndef BSACL_BASE_DIR
#    define BSACL_BASE_DIR ./
#  endif
#  define STRINGIFYMACRO_LITERAL(X) #X
#  define xstr(s) STRINGIFYMACRO_LITERAL(s)
#  define STRINGIFYMACRO_VALUE(X) xstr(X)
#  define concatenate(X, Y) X/Y
#  define INC concatenate(BSACL_BASE_DIR, _base.h)
#  include STRINGIFYMACRO_VALUE(INC)
#endif


#ifdef BSACL_USE_CUDA__
# include <cuda_runtime.h>
#else
# ifdef BSACL_REAL_IS_DOUBLE
#  pragma OPENCL EXTENSION cl_khr_fp64 : enable
#  undef BSACL_REAL_IS_DOUBLE
# endif
#endif


#define _use_fast_reduction_scheme_


DEVICE UINT getCorrId(
      PRIVATE UINT ni,
      PRIVATE UINT nj,
      const   UINT NTOT
)
{
   BOOL invert_ = (BOOL)(nj < ni);
   UINT res_;

   if (invert_) {
      UINT itemp_ = ni;
      ni = nj;
      nj = itemp_;
   }
   res_ = ni * (NTOT-1) + nj;
   res_ = res_ - (UINT)((ni*ni - ni) / 2.f);
   return res_;
}



#ifdef BSACL_USE_CUDA__
# ifdef BSACL_WIND_PSD_ID
#  undef BSACL_WIND_PSD_ID
# endif
#else
# ifndef BSACL_WIND_PSD_ID
#  define BSACL_WIND_PSD_ID 1
# endif
#endif




DEVICE BSACL_REAL evalFct(
      const BSACL_REAL f,
#ifdef BSACL_USE_CUDA__
      const UINT BSACL_WIND_PSD_ID,
#endif
      const BSACL_REAL w_scl,
      const BSACL_REAL w_std,
      const BSACL_REAL w_nodvel
)
{
   BSACL_REAL rtmp, res = (BSACL_REAL)0.f;

   const BSACL_REAL cL_U  = w_scl / w_nodvel;
   const BSACL_REAL cFL_U = f * cL_U;

#if (defined(BSACL_USE_CUDA__)) || (BSACL_WIND_PSD_ID==BSACL_PSD_TYPE_VONKARMAN)
# ifdef BSACL_USE_CUDA__
   if (BSACL_WIND_PSD_ID==BSACL_PSD_TYPE_VONKARMAN) {
# endif
      rtmp = cFL_U*cFL_U;
      rtmp = (rtmp * (BSACL_REAL)70.7f) + (BSACL_REAL)1.f;
      rtmp = POW(rtmp, ((BSACL_REAL)5.f/(BSACL_REAL)6.f));
      rtmp = 1.f / rtmp;

      res  = (4.f * cL_U * w_std*w_std) * rtmp;
# ifdef BSACL_USE_CUDA__
   }
# endif
#endif // BSACL_WIND_PSD_ID==BSACL_PSD_TYPE_VONKARMAN


#if (defined BSACL_USE_CUDA__) || (BSACL_WIND_PSD_ID==BSACL_PSD_TYPE_DAVENPORT)
# ifdef BSACL_USE_CUDA__
   if (BSACL_WIND_PSD_ID==BSACL_PSD_TYPE_DAVENPORT) {
# endif
      rtmp = cFL_U*cFL_U + 1.f;
      rtmp = POW(rtmp, ((BSACL_REAL)4.f/(BSACL_REAL)3.f));
      rtmp = 1.f / rtmp;

      res  = (2.f/3.f * cFL_U * cL_U * w_std*w_std) * rtmp;
# ifdef BSACL_USE_CUDA__
   }
# endif
#endif // BSACL_WIND_PSD_ID==BSACL_PSD_TYPE_DAVENPORT

#ifdef BSACL_CONV_PULSATION
   res =  res / (4 * BSACL_PI);
#endif

   return res;
} // evalFct




#ifdef BSACL_PASS_PARAMS_BY_MACRO__
#  ifndef NTC__
#   error NTC__  is not defined!
#  endif
#  ifndef NNL__
#   error NNL__  is not defined!
#  endif
#  ifndef NN__
#   error NN__  is not defined!
#  endif
#  ifndef NM_EFF__
#   error NM_EFF__  is not defined!
#  endif
#endif // BSACL_PASS_PARAMS_BY_MACRO__


#ifdef BSACL_USE_CUDA__
# define PSD_ID_ARG BSACL_WIND_PSD_ID,
#else 
# define PSD_ID_ARG
#endif





#if (BSACL_KERNEL_ID==3)

/**
 * BFM kernel using a total of NN^3 WI organised into NWGs.
 * Also, second dimension holds NM^3 WG, each one of which holds a unique 
 *   modal indexes combination.
 * Hence, each WI will uniquely perform the kernel for a single and 
 *   unique combination of NODAL indexes (I,J,K). 
 *
 * BSACL_KERNEL_ID==3:
 *  This version loads all the 2D frequency vectors, and compute loops internally
 *    to avoid multiple kernel enqueueing.
*/
KERNEL void bfm_kernel(
#ifdef BSACL_USE_CUDA__
      const    UINT          BSACL_WIND_PSD_ID,
#endif
#ifndef BSACL_PASS_PARAMS_BY_MACRO__
      const    UINT          NTC__,
#endif
      CONSTANT UINT          *tc,
#ifndef BSACL_PASS_PARAMS_BY_MACRO__
      const    UINT          NNL__,
#endif
      GLOBAL   UINT          *nodes_load,
      GLOBAL   BSACL_REAL          *fi,
      const    UINT          NFI__,
      GLOBAL   BSACL_REAL          *fj,
      const    UINT          NFJ__,
      const    BSACL_REAL          dInfl,
#ifndef BSACL_PASS_PARAMS_BY_MACRO__
      const    UINT          NM_EFF__,
      const    UINT          NDEGW__,
#endif
      const    GLOBAL   BSACL_REAL *phiTc,
#ifndef BSACL_PASS_PARAMS_BY_MACRO__
      const    UINT          NN__,
      const    UINT          NNOD_CORR__,   // BUG: NOT used
#endif
      const    GLOBAL   BSACL_REAL *nod_corr,
      const    GLOBAL   BSACL_REAL *wind_nod_vel,
      const    GLOBAL   BSACL_REAL *wind_turb_scl,
      const    GLOBAL   BSACL_REAL *wind_turb_std,
      const    GLOBAL   int  *wind_nod_winz,
      GLOBAL BSACL_REAL *m3mf
)
{
   const size_t gid0_  = GLOBAL_ID_X_DIM0;
   const size_t lid0_  = LOCAL_ID_X_DIM0;

   UINT itmp_ = (NNL__ * NNL__ * NNL__) - 1;
   if (gid0_ > itmp_) return;

   LOCAL  BSACL_REAL  m3mf_wg_x_[BSACL_WIpWG];
   m3mf_wg_x_[lid0_] = 0.f;

   /** get UNIQUE nodal indexes for this WI */
   itmp_ = NNL__ * NNL__;
   const UINT ink_ = (UINT)gid0_ / itmp_;
   itmp_           = (UINT)gid0_ - (ink_ * itmp_);
   const UINT inj_ = itmp_ / NNL__;
   const UINT ini_ = itmp_ - (inj_ * NNL__);
   const UINT ni_  = nodes_load[ini_]-1;
   const UINT nj_  = nodes_load[inj_]-1; 
   const UINT nk_  = nodes_load[ink_]-1;

   const BSACL_REAL ubni_ = wind_nod_vel[ni_];
   const BSACL_REAL ubnj_ = wind_nod_vel[nj_];
   const BSACL_REAL ubnk_ = wind_nod_vel[nk_];

   /** 
      BUG: consider all 3 spatial direction 
      BUG: we should do it on the frequency (exponent).
           However, we do it to the base cause we do it only once (faster).
           Otherwise, we should split the loops..
   */
   const BSACL_REAL corrIJ_ = nod_corr[getCorrId(ni_, nj_, NN__)] < BSACL_REAL_MIN ? BSACL_REAL_MIN : nod_corr[getCorrId(ni_, nj_, NN__)];
   const BSACL_REAL corrIK_ = nod_corr[getCorrId(ni_, nk_, NN__)] < BSACL_REAL_MIN ? BSACL_REAL_MIN : nod_corr[getCorrId(ni_, nk_, NN__)];
   const BSACL_REAL corrJK_ = nod_corr[getCorrId(nj_, nk_, NN__)] < BSACL_REAL_MIN ? BSACL_REAL_MIN : nod_corr[getCorrId(nj_, nk_, NN__)];


   /** determines which combination (M, N, O) of modal indexes apply to this WG. */
   const size_t wgid1_ = BLOCK_ID_Y_DIM1;
   itmp_ = NM_EFF__ * NM_EFF__;
   const UINT mk_ = (UINT)wgid1_ / itmp_;
   itmp_          = (UINT)wgid1_ - (mk_ * itmp_);
   const UINT mj_ = itmp_ / NM_EFF__;
   const UINT mi_ = itmp_ - (mj_ * NM_EFF__);


   /** NOTE: each WI has its own unique NODAL indexes!
             3 turb components.
             2 coeffs (Uu and u^2, for each t.c.)
             3 modal indexes
   */
   LOCAL  BSACL_REAL  phiTc_[6 * 3 * BSACL_WIpWG];

   const UINT phiTc_offst_ = 18U * (UINT)lid0_;

   itmp_ = NM_EFF__ * NNL__;
   phiTc_[0 + phiTc_offst_] = phiTc[mi_ + (ini_*NM_EFF__)];                // Uu
   phiTc_[1 + phiTc_offst_] = phiTc[mi_ + (ini_*NM_EFF__) + (1 * itmp_)];  // Uv
   phiTc_[2 + phiTc_offst_] = phiTc[mi_ + (ini_*NM_EFF__) + (2 * itmp_)];  // Uw

   phiTc_[3 + phiTc_offst_] = phiTc[mj_ + (inj_*NM_EFF__)];                // Uu
   phiTc_[4 + phiTc_offst_] = phiTc[mj_ + (inj_*NM_EFF__) + (1 * itmp_)];  // Uv
   phiTc_[5 + phiTc_offst_] = phiTc[mj_ + (inj_*NM_EFF__) + (2 * itmp_)];  // Uw

   phiTc_[6 + phiTc_offst_] = phiTc[mk_ + (ink_*NM_EFF__)];                // Uu
   phiTc_[7 + phiTc_offst_] = phiTc[mk_ + (ink_*NM_EFF__) + (1 * itmp_)];  // Uv
   phiTc_[8 + phiTc_offst_] = phiTc[mk_ + (ink_*NM_EFF__) + (2 * itmp_)];  // Uw


   phiTc_[9  + phiTc_offst_] = phiTc[mi_ + (ini_*NM_EFF__) + (3 * itmp_)]; // u^2
   phiTc_[10 + phiTc_offst_] = phiTc[mi_ + (ini_*NM_EFF__) + (4 * itmp_)]; // v^2
   phiTc_[11 + phiTc_offst_] = phiTc[mi_ + (ini_*NM_EFF__) + (5 * itmp_)]; // w^2

   phiTc_[12 + phiTc_offst_] = phiTc[mj_ + (inj_*NM_EFF__) + (3 * itmp_)]; // u^2
   phiTc_[13 + phiTc_offst_] = phiTc[mj_ + (inj_*NM_EFF__) + (4 * itmp_)]; // v^2
   phiTc_[14 + phiTc_offst_] = phiTc[mj_ + (inj_*NM_EFF__) + (5 * itmp_)]; // w^2

   phiTc_[15 + phiTc_offst_] = phiTc[mk_ + (ink_*NM_EFF__) + (3 * itmp_)]; // u^2
   phiTc_[16 + phiTc_offst_] = phiTc[mk_ + (ink_*NM_EFF__) + (4 * itmp_)]; // v^2
   phiTc_[17 + phiTc_offst_] = phiTc[mk_ + (ink_*NM_EFF__) + (5 * itmp_)]; // w^2


   BSACL_REAL S_uvw_IJ_i, S_uvw_IJ_ij;
   BSACL_REAL S_uvw_IK_j, S_uvw_IK_ij;
   BSACL_REAL S_uvw_JK_i, S_uvw_JK_j;

   for (UINT itc_ = 0; itc_ < NTC__; ++itc_) {

      UINT tc_   = tc[itc_] - 1;

      BSACL_REAL wstd_ = wind_turb_std[tc_];  // BUG: account for multiple wind zones!!
      BSACL_REAL wscl_ = wind_turb_scl[tc_];  // BUG: account for multiple wind zones!!

      for (UINT ifj_=0; ifj_ < NFJ__; ++ifj_) {

         BSACL_REAL fj_ = fj[ifj_];

         S_uvw_IK_j   = evalFct(fj_,  PSD_ID_ARG   wscl_, wstd_, ubni_);
         S_uvw_IK_j  *= evalFct(fj_,  PSD_ID_ARG   wscl_, wstd_, ubnk_);
         S_uvw_IK_j   = sqrt(S_uvw_IK_j);
         S_uvw_IK_j  *= POW(corrIK_, (FABS((BSACL_REAL)fj_)));

         S_uvw_JK_j   = evalFct(fj_,  PSD_ID_ARG   wscl_, wstd_, ubnj_);
         S_uvw_JK_j  *= evalFct(fj_,  PSD_ID_ARG   wscl_, wstd_, ubnk_);
         S_uvw_JK_j   = sqrt(S_uvw_JK_j);
         S_uvw_JK_j  *= POW(corrJK_, (FABS((BSACL_REAL)fj_)));

         for (UINT ifi_=0; ifi_ < NFI__; ++ifi_) {

            BSACL_REAL fi_    = fi[ifi_];
            BSACL_REAL fiPfj_ = fi_ + fj_;

            S_uvw_IJ_i   = evalFct(fi_,  PSD_ID_ARG   wscl_, wstd_, ubni_);
            S_uvw_IJ_i  *= evalFct(fi_,  PSD_ID_ARG   wscl_, wstd_, ubnj_);
            S_uvw_IJ_i   = sqrt(S_uvw_IJ_i);
            S_uvw_IJ_i  *= POW(corrIJ_, (FABS((BSACL_REAL)fi_)));

            S_uvw_IJ_ij  = evalFct(fiPfj_,  PSD_ID_ARG   wscl_, wstd_, ubni_);
            S_uvw_IJ_ij *= evalFct(fiPfj_,  PSD_ID_ARG   wscl_, wstd_, ubnj_);
            S_uvw_IJ_ij  = sqrt(S_uvw_IJ_ij);
            S_uvw_IJ_ij *= POW(corrIJ_, (FABS((BSACL_REAL)fiPfj_)));


            S_uvw_IK_ij  = evalFct(fiPfj_,  PSD_ID_ARG   wscl_, wstd_, ubni_);
            S_uvw_IK_ij *= evalFct(fiPfj_,  PSD_ID_ARG   wscl_, wstd_, ubnk_);
            S_uvw_IK_ij  = sqrt(S_uvw_IK_ij);
            S_uvw_IK_ij *= POW(corrIK_, (FABS((BSACL_REAL)fiPfj_)));


            S_uvw_JK_i   = evalFct(fi_,  PSD_ID_ARG   wscl_, wstd_, ubnj_);
            S_uvw_JK_i  *= evalFct(fi_,  PSD_ID_ARG   wscl_, wstd_, ubnk_);
            S_uvw_JK_i   = sqrt(S_uvw_JK_i);
            S_uvw_JK_i  *= POW(corrJK_, (FABS((BSACL_REAL)fi_)));


            m3mf_wg_x_[lid0_] += 2.f * (
                 phiTc_[9+phiTc_offst_+tc_] * phiTc_[3 +phiTc_offst_+tc_] * phiTc_[6 +phiTc_offst_+tc_] 
                  * (S_uvw_IJ_i  * S_uvw_IK_j )
               + phiTc_[0+phiTc_offst_+tc_] * phiTc_[12+phiTc_offst_+tc_] * phiTc_[6 +phiTc_offst_+tc_] 
                  * (S_uvw_IJ_ij * S_uvw_JK_j )
               + phiTc_[0+phiTc_offst_+tc_] * phiTc_[3 +phiTc_offst_+tc_] * phiTc_[15+phiTc_offst_+tc_] 
                  * (S_uvw_JK_i  * S_uvw_IK_ij)
            );

         } // fi
      } // fj
   } // NTC_

   // Multiply by reference area
   m3mf_wg_x_[lid0_] *= dInfl;

   // BUG: apparently, removing this barrier leads to wrong results..
   LOCAL_WORKGROUP_BARRIER;

#ifdef _use_fast_reduction_scheme_
   UINT active = BSACL_WIpWG;
   while (active > 1) {
      LOCAL_WORKGROUP_BARRIER;
      active /= 2;
      if (lid0_ < active) {
         m3mf_wg_x_[lid0_] += m3mf_wg_x_[lid0_+active];
      }
   }
   // Then store sum into global variable
   if (0 == lid0_) {
      itmp_ = (UINT)BLOCK_ID_X_DIM0;
      m3mf[(itmp_*NM_EFF__*NM_EFF__*NM_EFF__) + wgid1_] += m3mf_wg_x_[0];
   }
#else
   if (0 == lid0_) {

      /** Reduce among all WI of current WG. */
      for (itmp_ = 1; itmp_ < BSACL_WIpWG; itmp_++)
         m3mf_wg_x_[0] += m3mf_wg_x_[itmp_];

      itmp_ = (UINT)BLOCK_ID_X_DIM0;

      /** Reduce on global result array */
      m3mf[(itmp_*NM_EFF__*NM_EFF__*NM_EFF__) + wgid1_] += m3mf_wg_x_[0];
   }
#endif
   LOCAL_WORKGROUP_BARRIER;
}

#endif // (BSACL_KERNEL_ID==3)




#if (BSACL_KERNEL_ID==4)

/**
 * BFM kernel using a total of NF^2 x NM^3 work items.
 * Each work group is a 1D column-like vector, i.e. takes only 
 * one column. So that, each work group has a unique set of modes.
 * Each work group will iterate through all combinations of loaded nodes.
 */
KERNEL void bfm_kernel(
#ifdef BSACL_USE_CUDA__
      const    UINT          BSACL_WIND_PSD_ID,
#endif
#ifndef BSACL_PASS_PARAMS_BY_MACRO__
      const    UINT          NTC__,
#endif
      CONSTANT UINT          *tc,
#ifndef BSACL_PASS_PARAMS_BY_MACRO__
      const    UINT          NNL__,
#endif
      GLOBAL   UINT          *nodes_load,
      GLOBAL   __real        *fi,
#ifndef BSACL_PASS_PARAMS_BY_MACRO__
      const    UINT          NFI__,
#endif
      GLOBAL   __real        *fj,
#ifndef BSACL_PASS_PARAMS_BY_MACRO__
      const    UINT          NFJ__,
      const    UINT          NM_EFF__,
      const    UINT          NDEGW__,
#endif
      const    GLOBAL  __real  *phiTc,
#ifndef BSACL_PASS_PARAMS_BY_MACRO__
      const    UINT            NN__,
      const    UINT            NNOD_CORR__,   // BUG: NOT used
#endif
      const    GLOBAL  __real  *nod_corr,
      const    GLOBAL  __real  *wind_nod_vel,
      const    GLOBAL  __real  *wind_turb_scl,
      const    GLOBAL  __real  *wind_turb_std,
      const    GLOBAL   int    *wind_nod_winz,
      GLOBAL __real *m3mf
)
{
   const size_t gid0_ = GLOBAL_ID_X_DIM0; // this determines the pair of freqs
   const size_t lid0_ = LOCAL_ID_X_DIM0;

   /* NOTE: do set local variable to 0 for all WIs, avoids getting garbage in reduction op. */
   LOCAL BSACL_REAL m3mf_loc_[BSACL_WIpWG];
   m3mf_loc_[lid0_] = (BSACL_REAL)0.f;

   UINT itmp_ = (NFI__ * NFJ__) - 1;
   if (gid0_ > itmp_) return;


   /** determine which combination (M, N, O) of modal indexes apply to this WG. */
   const size_t wgid1_ = BLOCK_ID_Y_DIM1;
   itmp_ = NM_EFF__ * NM_EFF__;
   const UINT mk_ = (UINT)wgid1_ / itmp_;
   itmp_          = (UINT)wgid1_ - (mk_ * itmp_);
   const UINT mj_ = itmp_ / NM_EFF__;
   const UINT mi_ = itmp_ - (mj_ * NM_EFF__);


   /**
    * Prefetch phiTc_ for this WG m-n-o modal indexes.
    * phiTc_mno_ will be grouped by Nodal info.
    * [ [], [] ] where each inner [] contains everything at the nodal level.
    * Then, each inner [] will be itself divided as:
    *     [ [6 coeffs mode m], [6 coeffs mode n], [6 coeffs mode o]  ]
    * */
#ifdef BSACL_USE_CUDA__
   /* BUG: we can't use a RT parameter for automatic allocation.. */
   LOCAL BSACL_REAL phiTc_mno_[3 * 6 * 4]; // 6 : Uu, Uv, Uw, u2, v2, w2; 3 modes
#else
   LOCAL BSACL_REAL phiTc_mno_[3 * 6 * NNL__]; // 6 : Uu, Uv, Uw, u2, v2, w2; 3 modes
#endif
   const UINT phiTc_offst_ = NM_EFF__ * NNL__;
   if (BSACL_WIpWG < NNL__) {
      const BSACL_USHORT n_reps = NNL__ / BSACL_WIpWG;
      // full BSACL_WIpWG batches
      for (BSACL_USHORT r=0; r < n_reps; ++r) {
         for (BSACL_USHORT d=0; d < 6; ++d) {
            UINT nid = lid0_ + r*BSACL_WIpWG;
            phiTc_mno_[18*nid +      d] = (BSACL_REAL)phiTc[mi_ + (nid * NM_EFF__) + (d * phiTc_offst_)];
            phiTc_mno_[18*nid +  6 + d] = (BSACL_REAL)phiTc[mj_ + (nid * NM_EFF__) + (d * phiTc_offst_)];
            phiTc_mno_[18*nid + 12 + d] = (BSACL_REAL)phiTc[mk_ + (nid * NM_EFF__) + (d * phiTc_offst_)];
         }
      }
      // last batch covering (NNL__-(BSACL_WIpWG * n_reps))
      if (lid0_ < (NNL__ - (BSACL_WIpWG * n_reps))) {
         for (BSACL_USHORT d=0; d < 6; ++d) {
            UINT nid = lid0_ + n_reps*BSACL_WIpWG;
            phiTc_mno_[18*nid +      d] = (BSACL_REAL)phiTc[mi_ + (nid * NM_EFF__) + (d * phiTc_offst_)];
            phiTc_mno_[18*nid +  6 + d] = (BSACL_REAL)phiTc[mj_ + (nid * NM_EFF__) + (d * phiTc_offst_)];
            phiTc_mno_[18*nid + 12 + d] = (BSACL_REAL)phiTc[mk_ + (nid * NM_EFF__) + (d * phiTc_offst_)];
         }
      }
   } else { // NWI > NNL
      if (lid0_ < NNL__) {
         UINT nid = 18*lid0_;
         for (BSACL_USHORT d=0; d < 6; ++d) {
            phiTc_mno_[nid +      d] = (BSACL_REAL)phiTc[mi_ + (lid0_ * NM_EFF__) + (d * phiTc_offst_)];
            phiTc_mno_[nid +  6 + d] = (BSACL_REAL)phiTc[mj_ + (lid0_ * NM_EFF__) + (d * phiTc_offst_)];
            phiTc_mno_[nid + 12 + d] = (BSACL_REAL)phiTc[mk_ + (lid0_ * NM_EFF__) + (d * phiTc_offst_)];
         }
      }
   }
   LOCAL_WORKGROUP_BARRIER;


   // Compute pair of frequencies
   itmp_ = gid0_ / NFI__;
   const BSACL_REAL fj_    = (BSACL_REAL)fj[itmp_];
   itmp_ = gid0_ - (itmp_ * NFI__);
   const BSACL_REAL fi_    = (BSACL_REAL)fi[itmp_];
   const BSACL_REAL fiPfj_ = fi_ + fj_;


   BSACL_REAL nod_corr_;
   BSACL_REAL S_uvw_IJ_i, S_uvw_IJ_ij;
   BSACL_REAL S_uvw_IK_j, S_uvw_IK_ij;
   BSACL_REAL S_uvw_JK_i, S_uvw_JK_j;

   for (UINT itc_=0; itc_ < NTC__; ++itc_) {

      UINT tc_   = tc[itc_] - 1;

      BSACL_REAL wstd_ = (BSACL_REAL)wind_turb_std[tc_];  // BUG: account for multiple wind zones!!
      BSACL_REAL wscl_ = (BSACL_REAL)wind_turb_scl[tc_];  // BUG: account for multiple wind zones!!

      for (UINT ink_=0; ink_ < NNL__; ++ink_) {

         UINT nk_offs_ = 18*ink_;
         UINT nk_      = nodes_load[ink_]-1;
         BSACL_REAL ubnk_ = (BSACL_REAL)wind_nod_vel[nk_];

         BSACL_REAL S_uvw_K_i_  = evalFct(fi_,     PSD_ID_ARG   wscl_, wstd_, ubnk_);
         BSACL_REAL S_uvw_K_j_  = evalFct(fj_,     PSD_ID_ARG   wscl_, wstd_, ubnk_);
         BSACL_REAL S_uvw_K_ij_ = evalFct(fiPfj_,  PSD_ID_ARG   wscl_, wstd_, ubnk_);

         for (UINT inj_=0; inj_ < NNL__; ++inj_) {

            UINT nj_offs_ = 18*inj_;
            UINT nj_      = nodes_load[inj_]-1;
            BSACL_REAL ubnj_ = (BSACL_REAL)wind_nod_vel[nj_];

            nod_corr_ = (BSACL_REAL)nod_corr[getCorrId(nj_, nk_, NN__)];
            BSACL_REAL corrJK_   = nod_corr_ < BSACL_REAL_MIN ? BSACL_REAL_MIN : nod_corr_;

            BSACL_REAL S_uvw_J_i_ = evalFct(fi_,  PSD_ID_ARG   wscl_, wstd_, ubnj_);
            S_uvw_JK_i  = S_uvw_J_i_ * S_uvw_K_i_;
            S_uvw_JK_i  = sqrt(S_uvw_JK_i);
            S_uvw_JK_i *= POW(corrJK_, (FABS((BSACL_REAL)fi_)));

            S_uvw_JK_j  = evalFct(fj_,  PSD_ID_ARG   wscl_, wstd_, ubnj_);
            S_uvw_JK_j *= S_uvw_K_j_;
            S_uvw_JK_j  = sqrt(S_uvw_JK_j);
            S_uvw_JK_j *= POW(corrJK_, (FABS((BSACL_REAL)fj_)));

            BSACL_REAL S_uvw_J_ij_  = evalFct(fiPfj_,  PSD_ID_ARG   wscl_, wstd_, ubnj_);

            for (UINT ini_=0; ini_ < NNL__; ++ini_) {

               UINT ni_offs_ = 18*ini_;
               UINT ni_      = nodes_load[ini_]-1;
               BSACL_REAL ubni_ = (BSACL_REAL)wind_nod_vel[ni_];

               nod_corr_ = (BSACL_REAL)nod_corr[getCorrId(ni_, nk_, NN__)];
               BSACL_REAL corrIK_ = nod_corr_ < BSACL_REAL_MIN ? BSACL_REAL_MIN : nod_corr_;
               nod_corr_ = (BSACL_REAL)nod_corr[getCorrId(ni_, nj_, NN__)];
               BSACL_REAL corrIJ_ = nod_corr_ < BSACL_REAL_MIN ? BSACL_REAL_MIN : nod_corr_;

               S_uvw_IK_j   = evalFct(fj_,  PSD_ID_ARG   wscl_, wstd_, ubni_);
               S_uvw_IK_j  *= S_uvw_K_j_;
               S_uvw_IK_j   = sqrt(S_uvw_IK_j);
               S_uvw_IK_j  *= POW(corrIK_, (FABS((BSACL_REAL)fj_)));

               S_uvw_IK_ij  = evalFct(fiPfj_,  PSD_ID_ARG   wscl_, wstd_, ubni_);
               S_uvw_IK_ij *= S_uvw_K_ij_;
               S_uvw_IK_ij  = sqrt(S_uvw_IK_ij);
               S_uvw_IK_ij *= POW(corrIK_, (FABS((BSACL_REAL)fiPfj_)));

               S_uvw_IJ_i   = evalFct(fi_,  PSD_ID_ARG   wscl_, wstd_, ubni_);
               S_uvw_IJ_i  *= S_uvw_J_i_;
               S_uvw_IJ_i   = sqrt(S_uvw_IJ_i);
               S_uvw_IJ_i  *= POW(corrIJ_, (FABS((BSACL_REAL)fi_)));

               S_uvw_IJ_ij  = evalFct(fiPfj_,  PSD_ID_ARG   wscl_, wstd_, ubni_);
               S_uvw_IJ_ij *= S_uvw_J_ij_;
               S_uvw_IJ_ij  = sqrt(S_uvw_IJ_ij);
               S_uvw_IJ_ij *= POW(corrIJ_, (FABS((BSACL_REAL)fiPfj_)));

               m3mf_loc_[lid0_] += 2.f * (
                    phiTc_mno_[ni_offs_ + 3 + tc_] * phiTc_mno_[nj_offs_ + 6 +     tc_] * phiTc_mno_[nk_offs_ + 12 +     tc_] * (S_uvw_IJ_i  * S_uvw_IK_j )
                  + phiTc_mno_[ni_offs_ +     tc_] * phiTc_mno_[nj_offs_ + 6 + 3 + tc_] * phiTc_mno_[nk_offs_ + 12 +     tc_] * (S_uvw_IJ_ij * S_uvw_JK_j )
                  + phiTc_mno_[ni_offs_ +     tc_] * phiTc_mno_[nj_offs_ + 6 +     tc_] * phiTc_mno_[nk_offs_ + 12 + 3 + tc_] * (S_uvw_JK_i  * S_uvw_IK_ij)
               );
            }
         }
      }
   }

#ifdef _use_fast_reduction_scheme_
   UINT active = BSACL_WIpWG;
   while (active > 1) {
      LOCAL_WORKGROUP_BARRIER;
      active /= 2;
      if (lid0_ < active) {
         m3mf_loc_[lid0_] += m3mf_loc_[lid0_+active];
      }
   }
   // Then store sum into global variable
   const size_t wgid0_ = BLOCK_ID_X_DIM0;
   const size_t nwgd0_ = N_BLOCKS_WGROUPS_X_DIM0;
   if (0 == lid0_) {
      m3mf[wgid1_*nwgd0_ + wgid0_] = (__real)m3mf_loc_[0];
   }
#else
   if (0 == lid0_) {

      /** Reduce among all WI of current WG. */
      for (itmp_ = 1; itmp_ < BSACL_WIpWG; ++itmp_)
         m3mf_loc_[0] += m3mf_loc_[itmp_];

      // Then store sum into global variable
      const size_t wgid0_ = BLOCK_ID_X_DIM0;
      const size_t nwgd0_ = N_BLOCKS_WGROUPS_X_DIM0;
      m3mf[wgid1_*nwgd0_ + wgid0_] = (__real)m3mf_loc_[0];
   }
#endif
   LOCAL_WORKGROUP_BARRIER;
}

#endif // (BSACL_KERNEL_ID==4)

