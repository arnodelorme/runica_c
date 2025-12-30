/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_runica_simple_api.h
 *
 * Code generation for function 'runica_simple'
 *
 */

#ifndef _CODER_RUNICA_SIMPLE_API_H
#define _CODER_RUNICA_SIMPLE_API_H

/* Include files */
#include "emlrt.h"
#include "mex.h"
#include "tmwtypes.h"
#include <string.h>

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void runica_simple(real_T data[976128], boolean_T extended,
                   real_T weights[1024], real_T sphere[1024]);

void runica_simple_api(const mxArray *const prhs[2], int32_T nlhs,
                       const mxArray *plhs[2]);

void runica_simple_atexit(void);

void runica_simple_initialize(void);

void runica_simple_terminate(void);

void runica_simple_xil_shutdown(void);

void runica_simple_xil_terminate(void);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (_coder_runica_simple_api.h) */
