#ifndef BSACL_EXTDATA__H_
#define BSACL_EXTDATA__H_

#ifndef NULL
#  define NULL 0
#endif

typedef void (*evalFct_t)(int, int, const double*, int, double*);

typedef struct extdata_t {
   
   int          *nodes_load__;  // list of loaded nodes in the structural model
   unsigned int NNODES_LOAD__;
   unsigned int NN__;
   
   unsigned int PSD_ID__;

   unsigned int NLIBS__;
   unsigned int NDOFS__;

   unsigned int NMODES_EFF__;
   unsigned int NDEGW__;

   unsigned int NTC__;
   unsigned int NNOD_CORR__; // n. of nodal correlation entries (per direction!)

   unsigned int DIM_M3_M__;
   unsigned int NMODES__;

   unsigned int NWZ__;  // N. of wind zones

   double       *modmat__;
   double       *natfreqs__;
   
   int *modes_eff__;    // list of effective modes used
   int *tc__;           // list of effective turbulent components
   
   double *wfc__;       // wind force coefficients

   double *phi_T_c__;   // Phi x C matrix

   double *nod_corr__;  // nodal correlation array (no duplicates)

   double *wind_nodal_vel__;
   double *wind_turb_scales__;
   double *wind_turb_std__;
   int    *wind_nodal_windz__;

   double *m3mf__;
} extdata_t;


static inline void freeExtData(extdata_t *extdata) {
   extdata->nodes_load__  = NULL;
   extdata->modmat__      = NULL;
   extdata->natfreqs__    = NULL;
   extdata->modes_eff__   = NULL;
   extdata->wfc__         = NULL;
   extdata->nod_corr__    = NULL;
   extdata->phi_T_c__     = NULL;
   extdata->wind_nodal_vel__   = NULL;
   extdata->wind_nodal_windz__ = NULL;
   extdata->wind_turb_scales__ = NULL;
   extdata->wind_turb_std__    = NULL;
   extdata->NWZ__         = 0;
   extdata->NNODES_LOAD__ = 0;
   extdata->NLIBS__       = 0;
   extdata->NDOFS__       = 0;
   extdata->NMODES__      = 0;
   extdata->NMODES_EFF__  = 0;
   extdata->NDEGW__       = 0;
   extdata->NTC__         = 0;
   extdata->NNOD_CORR__   = 0;
   extdata->PSD_ID__      = 0;
}



#endif // BSACL_EXTDATA__H_