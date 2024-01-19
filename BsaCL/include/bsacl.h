#pragma once

#include <stdint.h>

#include "_macros.h"
#include "_extdata.h"
#include "_devtypes.h"
#include "_msgtypes.h"
#include "_errtypes.h"
#include "_base.h"

#ifdef __cplusplus
extern "C" {
#endif

void bsaclInit(int *__EXT_PTR_CONST ierr);
void bsaclRun(int *__EXT_PTR_CONST ierr);
void bsaclAcquirePSDId(const uint32_t psdid);
void bsaclAcquireStructModMat(REAL *__EXT_PTR_CONST modmat, REAL *__EXT_PTR_CONST natf, const uint32_t ndofs, const uint32_t nmodes);
void bsaclAcquireLoadedNodesList(int *__EXT_PTR_CONST nodes_load, const uint32_t nnodes_l);
void bsaclAcquireUsedModesList(int *__EXT_PTR_CONST modes, const uint32_t nmodes_eff);
void bsaclAcquireWindCoeffs(REAL *__EXT_PTR_CONST wfc, const uint32_t nnodes_l, const uint32_t nlibs, const uint32_t ndegw);
void bsaclAcquireTurbComponentsList(int *__EXT_PTR_CONST tc, const uint32_t ntc);
void bsaclAcquirePhiTimesCMat(REAL *__EXT_PTR_CONST phi_T_c, const uint32_t nmodes_eff, const uint32_t nnodes_l, const uint32_t ndegw);
void bsaclAcquireNodalCorrelation(REAL *__EXT_PTR_CONST nod_corr, const uint32_t nnod_corr);
void bsaclAcquireWindNodalVelocities(REAL *__EXT_PTR_CONST nod_vel);
void bsaclAcquireWindNodalWindZones(int *__EXT_PTR_CONST nod_wz);
void bsaclAcquireWindTurbScales(REAL *__EXT_PTR_CONST wt_scl, const uint32_t nwz);
void bsaclAcquireWindTurbStd(REAL *__EXT_PTR_CONST wt_std, const uint32_t nwz);
void bsaclAcquireEvalFunc(evalFct_t fct);
void bsaclAcquireEvalFuncByStrings(char **__EXT_PTR_CONST strings);
void bsaclAcquireEvalFuncByFile(char *__EXT_PTR_CONST filename);
void bsaclAcquireComputationFreqs(const uint32_t nfi, REAL *__EXT_PTR_CONST fi, const uint32_t nfj, REAL *__EXT_PTR_CONST fj);
void bsaclAcquireBaseWindTurbPSD(REAL *__EXT_PTR_CONST S_uvw);
void bsaclAcquireResultBFMVect(REAL *__EXT_PTR_CONST m3mf, const uint32_t idim);
void bsaclSetDeviceType(const uint32_t itype);
void bsaclVerifyMaxAllocCondition(size_t idim, uint32_t *__EXT_PTR_CONST ican);
void bsaclAbort(const int ierr);
void bsaclFinalise(void);

#ifdef __cplusplus
}
#endif