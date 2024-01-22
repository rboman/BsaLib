#define STRINGIFYMACRO_LITERAL(X) #X
#define xstr(s) STRINGIFYMACRO_LITERAL(s)
#define STRINGIFYMACRO_VALUE(X) xstr(X)

#ifdef BSACL_USE_CUDA__
#  include "_base.h"
#else
#  ifdef BSACL_PI
#    undef  BSACL_PI
#    define BSACL_PI M_PI
#  endif
#  ifndef BSACL_BASE_DIR
#    define BSACL_BASE_DIR ./
#  endif
#  define concatenate(X, Y) X/Y
#  define INC concatenate(BSACL_BASE_DIR, _base.h)
#  include STRINGIFYMACRO_VALUE(INC)
#endif


#ifndef REAL
# define REAL double
# define REAL_IS_DOUBLE
#endif


#ifdef BSACL_USE_CUDA__
# include <cuda_runtime.h>
# pragma message("   --- [NOTE]:  Compiling kernel using  CUDA  runtime")
#else
# ifdef REAL_IS_DOUBLE
#  pragma OPENCL EXTENSION cl_khr_fp64 : enable
# endif
#endif


// NOTE: default define BSACL_WIpWG if not passed as argument when 
//       compiling on-the-fly this CL source.
#ifndef BSACL_WIpWG
#  define BSACL_WIpWG 64
#endif

#ifndef BSACL_KERNEL_ID
#  define BSACL_KERNEL_ID 3
#endif


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









#if (BSACL_KERNEL_ID==1)

/**
 * BFM kernel using NWG work-groups (1D)
 * Each work-group takes a slice NM_EFF^3 / NWG of global result array.
 *
 * I.e. each work item will loop on:
 *   - n. of turb. components  (NTC)
 *   - 3-loops on n. of loaded nodes (NNL)
 *
 * This kernel uses pre-computed Phi x C matrix product, 
 *  for avoiding explicit loops on nodal degrees of freedom (LIBs)
*/
KERNEL void bfm_kernel(
      const    UINT       NTC,
      CONSTANT UINT      *tc,
      const    UINT       NNL,
      GLOBAL   UINT      *nodes_load,
      const    GLOBAL   REAL  *S_uvw_fi,
      CONSTANT REAL           *fi_abs,
      const    GLOBAL   REAL  *S_uvw_fj,
      CONSTANT REAL           *fj_abs,
      const    GLOBAL   REAL  *S_uvw_fiPfj,
      CONSTANT REAL           *fiPfj_abs,
      const    REAL            dInfl,
      const    UINT              NM_EFF,
      const    UINT              NDEGW,
      const    GLOBAL   REAL    *phiTc,
      const    UINT              NN,
      const    UINT              NNOD_CORR,
      const    GLOBAL   REAL    *nod_corr,
      const    GLOBAL   REAL  *wind_nod_vel,
      const    GLOBAL   REAL  *wind_turb_scl,
      const    GLOBAL   REAL  *wind_turb_std,
      const    GLOBAL   REAL  *wind_nod_winz,
      GLOBAL REAL *m3mf
) 
{

   const size_t id_ = GLOBAL_ID_X_DIM0;
   UINT itmp_ = (NM_EFF * NM_EFF * NM_EFF) - 1;
   if (id_ > itmp_) return;

#ifndef TEST

   LOCAL REAL  S_uvw_K_i;
   LOCAL REAL  S_uvw_K_j;
   LOCAL REAL  S_uvw_K_ij;
   LOCAL REAL  coorJK;
   LOCAL REAL  S_uvw_J_i;
   LOCAL REAL  S_uvw_J_j;
   LOCAL REAL  S_uvw_J_ij;
   LOCAL REAL  S_uvw_JK_i;
   LOCAL REAL  S_uvw_JK_j;
   LOCAL REAL  coorIJ;
   LOCAL REAL  coorIK;
   LOCAL REAL  S_uvw_I_i;
   LOCAL REAL  S_uvw_I_j;
   LOCAL REAL  S_uvw_I_ij;
   LOCAL REAL  S_uvw_IJ_i;
   LOCAL REAL  S_uvw_IJ_ij;
   LOCAL REAL  S_uvw_IK_j;
   LOCAL REAL  S_uvw_IK_ij;
   LOCAL REAL  BF_ijk_IJK_;

   itmp_ = NM_EFF * NM_EFF;

   UINT imk_ = (UINT)id_ / itmp_;
   itmp_ = (UINT)id_ - (imk_ * itmp_);
   UINT imj_ = itmp_ / NM_EFF;
   UINT imi_ = itmp_ - (imj_ * NM_EFF);


   // if (id_ != 0) return;
   // if (id_ != (NM_EFF*NM_EFF - 1)) return;

   // printf("\n imk - imj - imi =   %d - %d - %d", imk_, imj_, imi_);
   // printf("\n NTC       = %d", NTC);
   // printf("\n NNL       = %d", NNL);
   // printf("\n NN        = %d", NN);
   // printf("\n NM_EFF    = %d", NM_EFF);
   // printf("\n NDEGW     = %d", NDEGW);
   // printf("\n NNOD_CORR = %d", NNOD_CORR);
   // printf("\n fi_abs    = %f", (float)*fi_abs);
   // printf("\n fj_abs    = %f", (float)*fj_abs);
   // printf("\n fiPfj_abs = %f", (float)*fiPfj_abs);
   // printf("\n dInfl     = %f", dInfl);

   // printf("\n\n TCs:\n");
   // for (BSACL_USHORT i = 0; i < NTC; ++i)
   //    printf("  %d", tc[i]);

   // printf("\n\n NODES LOAD:\n");
   // for (BSACL_USHORT i = 0; i < NNL; ++i)
   //    printf("  %d", nodes_load[i]);


   UINT nodcorrID_IJ_;
   UINT nodcorrID_IK_;
   UINT nodcorrID_JK_;


   for (UINT itc_ = 0; itc_ < NTC; ++itc_) {

      UINT tc_offs_ = itc_ * NNL;

      itmp_ = tc[itc_] - 1;
      UINT nod_corr_offs_ = itmp_ * NNOD_CORR;
      UINT tc_offs_phiTc_ = itmp_ * (NM_EFF * NNL);

      itmp_ = itmp_ + 3U;
      UINT tcP3_offs_ = itmp_ * (NM_EFF * NNL);

      for (UINT ink_ = 0; ink_ < NNL; ++ink_) {

         UINT ink_offs_ = ink_ * NM_EFF;
         UINT       nk_ = nodes_load[ink_] - 1;

         S_uvw_K_i  = S_uvw_fi[   ink_ + tc_offs_];
         S_uvw_K_j  = S_uvw_fj[   ink_ + tc_offs_];
         S_uvw_K_ij = S_uvw_fiPfj[ink_ + tc_offs_];

         // printf("\n Accessing   S_uvw_fiPfj[%d]  =  %f .", ink_ + tc_offs_, S_uvw_fiPfj[ink_ + tc_offs_]);

         for (UINT inj_ = 0; inj_ < NNL; ++inj_) {

            UINT inj_offs_ = inj_ * NM_EFF;
            UINT       nj_ = nodes_load[inj_] - 1;

            nodcorrID_JK_  = getCorrId(nj_, nk_, NN);
            coorJK         = nod_corr[nodcorrID_JK_ + nod_corr_offs_];


            // if (itc_ == 0)
            //    printf("\n Accessing   nod_corr[%d]  (inj - ink) = (%d-%d)  =  %f", 
            //       nodcorrID_JK_ + itmp_, inj_+1, ink_+1, nod_corr[nodcorrID_JK_ + itmp_]);


            S_uvw_J_i  = S_uvw_fi[   inj_ + tc_offs_];
            S_uvw_J_j  = S_uvw_fj[   inj_ + tc_offs_];
            S_uvw_J_ij = S_uvw_fiPfj[inj_ + tc_offs_];

            S_uvw_JK_i = pow(coorJK, *fi_abs) * sqrt(S_uvw_J_i * S_uvw_K_i);
            S_uvw_JK_j = pow(coorJK, *fj_abs) * sqrt(S_uvw_J_j * S_uvw_K_j);


            // if (inj_ == 0 && ink_ == 0)
            //    printf("\n  S_uvw_K_ij  =  %f", S_uvw_JK_i);


            for (UINT ini_ = 0; ini_ < NNL; ++ini_) {

               UINT ini_offs_ = ini_ * NM_EFF;
               UINT       ni_ = nodes_load[ini_] - 1;

               nodcorrID_IJ_ = getCorrId(ni_, nj_, NN);
               nodcorrID_IK_ = getCorrId(ni_, nk_, NN);
               coorIJ        = nod_corr[nodcorrID_IJ_ + nod_corr_offs_];
               coorIK        = nod_corr[nodcorrID_IK_ + nod_corr_offs_];

               S_uvw_I_i  = S_uvw_fi[   ini_ + tc_offs_];
               S_uvw_I_j  = S_uvw_fj[   ini_ + tc_offs_];
               S_uvw_I_ij = S_uvw_fiPfj[ini_ + tc_offs_];

               S_uvw_IJ_i  = pow(coorIJ, *fi_abs   ) * sqrt(S_uvw_I_i  * S_uvw_J_i );
               S_uvw_IJ_ij = pow(coorIJ, *fiPfj_abs) * sqrt(S_uvw_I_ij * S_uvw_J_ij);

               S_uvw_IK_j  = pow(coorIK, *fj_abs   ) * sqrt(S_uvw_I_j  * S_uvw_K_j );
               S_uvw_IK_ij = pow(coorIK, *fiPfj_abs) * sqrt(S_uvw_I_ij * S_uvw_K_ij);


               // if (ink_ == 0 && inj_ == 0 && ini_ == 0) {
               //    printf("\n phiTc_U[%d]  =  %f", imk_ + ink_offs_ + tc_offs_phiTc_, phiTc[imk_ + ink_offs_ + tc_offs_phiTc_]);
               //    printf("\n phiTc_u[%d]  =  %f", imk_ + ink_offs_ + tcP3_offs_, phiTc[imk_ + ink_offs_ + tcP3_offs_]);
               // }



               // BUG: prefetch phiTc in local work-group memory
               BF_ijk_IJK_ = 
                  2. * (
                     (phiTc[imi_ + ini_offs_ + tcP3_offs_] * 
                        phiTc[imj_ + inj_offs_ + tc_offs_phiTc_  ] * 
                           phiTc[imk_ + ink_offs_ + tc_offs_phiTc_  ] * (S_uvw_IJ_i  * S_uvw_IK_j ))
                  +  (phiTc[imi_ + ini_offs_ + tc_offs_phiTc_  ] * 
                        phiTc[imj_ + inj_offs_ + tcP3_offs_] * 
                           phiTc[imk_ + ink_offs_ + tc_offs_phiTc_  ] * (S_uvw_IJ_ij * S_uvw_JK_j ))
                  +  (phiTc[imi_ + ini_offs_ + tc_offs_phiTc_  ] * 
                        phiTc[imj_ + inj_offs_ + tc_offs_phiTc_  ] * 
                           phiTc[imk_ + ink_offs_ + tcP3_offs_] * (S_uvw_JK_i  * S_uvw_IK_ij))
                  );


               // printf("\n  BF_ijk_IJK_  =  %f", BF_ijk_IJK_);


               // Integral update (global)
               m3mf[id_] += BF_ijk_IJK_ * dInfl;


            } // ni
         } // nj

      } // nk


      // printf("\n");


   } // ntc


#else

   // printf("\n  id  =  %llu", id_);
   m3mf[id_] += (REAL)(id_ + 1);

#endif // not defined TEST
}

#endif // BSACL_KERNEL_ID==1









#if (BSACL_KERNEL_ID==2) || (BSACL_KERNEL_ID==3)



#ifdef BSACL_USE_CUDA__
# ifdef BSACL_WIND_PSD_ID
#  undef BSACL_WIND_PSD_ID
# endif
#else
# ifndef BSACL_WIND_PSD_ID
#  define BSACL_WIND_PSD_ID 1
# endif
#endif




DEVICE REAL evalFct(
      const REAL f,
# ifdef BSACL_USE_CUDA__
      const UINT BSACL_WIND_PSD_ID,
# endif
      const REAL w_scl,
      const REAL w_std,
      const REAL w_nodvel
) 
{
   REAL rtmp, res = (REAL)0;

   const REAL cL_U  = w_scl / w_nodvel;
   const REAL cFL_U = f * cL_U;

#if (defined BSACL_USE_CUDA__) || (BSACL_WIND_PSD_ID==BSACL_PSD_TYPE_VONKARMAN)
# ifdef BSACL_USE_CUDA__
   if (BSACL_WIND_PSD_ID==BSACL_PSD_TYPE_VONKARMAN) {
# endif
      rtmp = cFL_U*cFL_U;
      rtmp = rtmp * 70.7f + 1;
      rtmp = POWR(rtmp, (REAL)(5.f/6.f));
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
      rtmp = POWR(rtmp, (REAL)(4.f/3.f));
      rtmp = 1.f / rtmp;

      res  = (2.f/3.f * cFL_U * cL_U * w_std*w_std) * rtmp;
# ifdef BSACL_USE_CUDA__
   }
# endif
#endif // BSACL_WIND_PSD_ID==BSACL_PSD_TYPE_DAVENPORT

#ifdef BSACL_CONV_PULSATION
   res =  res / (4*BSACL_PI);
#endif

   return res;
} // evalFct




# ifdef BSACL_PASS_PARAMS_BY_MACRO
#   ifndef NTC__
#    error NTC__  is not defined!
#   endif
#   ifndef NNL__
#    error NNL__  is not defined!
#   endif
#   ifndef NN__
#    error NN__  is not defined!
#   endif
#   ifndef NM_EFF__
#    error NM_EFF__  is not defined!
#   endif
# endif // BSACL_PASS_PARAMS_BY_MACRO

// # if (BSACL_KERNEL_ID==3)
// #  ifndef NFI__
// #   error NFI__  is not defined!
// #  endif
// #  ifndef NFJ__
// #   error NFJ__  is not defined!
// #  endif
// # endif


# ifdef BSACL_USE_CUDA__
#  define PSD_ID_ARG BSACL_WIND_PSD_ID,
# else 
#  define PSD_ID_ARG
# endif


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
      const    UINT              BSACL_WIND_PSD_ID,
#endif
#ifndef BSACL_PASS_PARAMS_BY_MACRO
      const    UINT              NTC__, 
#endif
      CONSTANT UINT             *tc,
#ifndef BSACL_PASS_PARAMS_BY_MACRO
      const    UINT              NNL__, 
#endif
      GLOBAL   UINT             *nodes_load,
# if (BSACL_KERNEL_ID==2)
      const    REAL          fi_,
      const    REAL          fj_,
# else
      GLOBAL   REAL         *fi,
      const    UINT          NFI__,
      GLOBAL   REAL         *fj,
      const    UINT          NFJ__,
# endif
      const    REAL          dInfl,
#ifndef BSACL_PASS_PARAMS_BY_MACRO
      const    UINT              NM_EFF__,
      const    UINT              NDEGW__, 
#endif
      const    GLOBAL   REAL    *phiTc,
#ifndef BSACL_PASS_PARAMS_BY_MACRO
      const    UINT              NN__,
      const    UINT              NNOD_CORR__,   // BUG: NOT used
#endif
      const    GLOBAL   REAL    *nod_corr,
      const    GLOBAL   REAL  *wind_nod_vel,
      const    GLOBAL   REAL  *wind_turb_scl,
      const    GLOBAL   REAL  *wind_turb_std,
      const    GLOBAL   int   *wind_nod_winz,
      GLOBAL REAL *m3mf
) 
{
   const size_t gid0_  = GLOBAL_ID_X_DIM0;
   const size_t lid0_  = LOCAL_ID_X_DIM0;

   UINT itmp_ = (NNL__ * NNL__ * NNL__) - 1;
   if (gid0_ > itmp_) return;

   LOCAL  REAL  m3mf_wg_x_[BSACL_WIpWG];
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

   const REAL ubni_ = wind_nod_vel[ni_];
   const REAL ubnj_ = wind_nod_vel[nj_];
   const REAL ubnk_ = wind_nod_vel[nk_];

   /** 
      BUG: consider all 3 spatial direction 
      BUG: we should do it on the frequency (exponent).
           However, we do it to the base cause we do it only once (faster).
           Otherwise, we should split the loops..
   */
   const REAL corrIJ_ = nod_corr[getCorrId(ni_, nj_, NN__)] < REAL_MIN ? REAL_MIN : nod_corr[getCorrId(ni_, nj_, NN__)];
   const REAL corrIK_ = nod_corr[getCorrId(ni_, nk_, NN__)] < REAL_MIN ? REAL_MIN : nod_corr[getCorrId(ni_, nk_, NN__)];
   const REAL corrJK_ = nod_corr[getCorrId(nj_, nk_, NN__)] < REAL_MIN ? REAL_MIN : nod_corr[getCorrId(nj_, nk_, NN__)];


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
   LOCAL  REAL  phiTc_[6 * 3 * BSACL_WIpWG];

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


   REAL S_uvw_IJ_i, S_uvw_IJ_ij;
   REAL S_uvw_IK_j, S_uvw_IK_ij;
   REAL S_uvw_JK_i, S_uvw_JK_j;

   for (UINT itc_ = 0; itc_ < NTC__; ++itc_) {

      UINT tc_   = tc[itc_] - 1;

      REAL wstd_ = wind_turb_std[tc_];  // BUG: account for multiple wind zones!!
      REAL wscl_ = wind_turb_scl[tc_];  // BUG: account for multiple wind zones!!

#if (BSACL_KERNEL_ID==3)
      for (UINT ifj_=0; ifj_ < NFJ__; ++ifj_) {

         REAL fj_ = fj[ifj_];
#endif

         S_uvw_IK_j   = evalFct(fj_,  PSD_ID_ARG   wscl_, wstd_, ubni_);
         S_uvw_IK_j  *= evalFct(fj_,  PSD_ID_ARG   wscl_, wstd_, ubnk_);
         S_uvw_IK_j   = sqrt(S_uvw_IK_j);
         S_uvw_IK_j  *= POWR(corrIK_, (REAL)(fabs(fj_)));

         S_uvw_JK_j   = evalFct(fj_,  PSD_ID_ARG   wscl_, wstd_, ubnj_);
         S_uvw_JK_j  *= evalFct(fj_,  PSD_ID_ARG   wscl_, wstd_, ubnk_);
         S_uvw_JK_j   = sqrt(S_uvw_JK_j);
         S_uvw_JK_j  *= POWR(corrJK_, (REAL)(fabs(fj_)));

#if (BSACL_KERNEL_ID==3)
         for (UINT ifi_=0; ifi_ < NFI__; ++ifi_) {

            REAL fi_    = fi[ifi_];
#endif
            REAL fiPfj_ = fi_ + fj_;

            S_uvw_IJ_i   = evalFct(fi_,  PSD_ID_ARG   wscl_, wstd_, ubni_);
            S_uvw_IJ_i  *= evalFct(fi_,  PSD_ID_ARG   wscl_, wstd_, ubnj_);
            S_uvw_IJ_i   = sqrt(S_uvw_IJ_i);
            S_uvw_IJ_i  *= POWR(corrIJ_, (REAL)(fabs(fi_)));

            S_uvw_IJ_ij  = evalFct(fiPfj_,  PSD_ID_ARG   wscl_, wstd_, ubni_);
            S_uvw_IJ_ij *= evalFct(fiPfj_,  PSD_ID_ARG   wscl_, wstd_, ubnj_);
            S_uvw_IJ_ij  = sqrt(S_uvw_IJ_ij);
            S_uvw_IJ_ij *= POWR(corrIJ_, (REAL)(fabs(fiPfj_)));


            S_uvw_IK_ij  = evalFct(fiPfj_,  PSD_ID_ARG   wscl_, wstd_, ubni_);
            S_uvw_IK_ij *= evalFct(fiPfj_,  PSD_ID_ARG   wscl_, wstd_, ubnk_);
            S_uvw_IK_ij  = sqrt(S_uvw_IK_ij);
            S_uvw_IK_ij *= POWR(corrIK_, (REAL)(fabs(fiPfj_)));


            S_uvw_JK_i   = evalFct(fi_,  PSD_ID_ARG   wscl_, wstd_, ubnj_);
            S_uvw_JK_i  *= evalFct(fi_,  PSD_ID_ARG   wscl_, wstd_, ubnk_);
            S_uvw_JK_i   = sqrt(S_uvw_JK_i);
            S_uvw_JK_i  *= POWR(corrJK_, (REAL)(fabs(fi_)));


            m3mf_wg_x_[lid0_] += 2.f * (
                 phiTc_[9+phiTc_offst_+tc_] * phiTc_[3 +phiTc_offst_+tc_] * phiTc_[6 +phiTc_offst_+tc_] 
                  * (S_uvw_IJ_i  * S_uvw_IK_j )
               + phiTc_[0+phiTc_offst_+tc_] * phiTc_[12+phiTc_offst_+tc_] * phiTc_[6 +phiTc_offst_+tc_] 
                  * (S_uvw_IJ_ij * S_uvw_JK_j )
               + phiTc_[0+phiTc_offst_+tc_] * phiTc_[3 +phiTc_offst_+tc_] * phiTc_[15+phiTc_offst_+tc_] 
                  * (S_uvw_JK_i  * S_uvw_IK_ij)
            );

#if (BSACL_KERNEL_ID==3)
         } // fi
      } // fj
#endif

   } // NTC_

   // Multiply by reference area
   m3mf_wg_x_[lid0_] *= dInfl;

   // NOTE: apparently, removing this barrier leads to wrong results..
   LOCAL_WORKGROUP_BARRIER;


   // BUG: unoptimal reduction scheme !!
   if (0 == lid0_) {

      /** Reduce among all WI of current WG. */
      for (itmp_ = 1; itmp_ < BSACL_WIpWG; itmp_++)
         m3mf_wg_x_[0] += m3mf_wg_x_[itmp_];

      itmp_ = (UINT)BLOCK_ID_X_DIM0;

      /** Reduce on global result array */
      m3mf[(itmp_*NM_EFF__*NM_EFF__*NM_EFF__) + wgid1_] += m3mf_wg_x_[0];
   }
   LOCAL_WORKGROUP_BARRIER;
}

#endif // (BSACL_KERNEL_ID==2) || (BSACL_KERNEL_ID==3)

