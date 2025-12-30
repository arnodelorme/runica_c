/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * runica_simple.h
 *
 * Code generation for function 'runica_simple'
 *
 */

#ifndef RUNICA_SIMPLE_H
#define RUNICA_SIMPLE_H

/* Include files */
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
extern void runica_simple(double data[976128], boolean_T extended,
                          double weights[1024], double sphere[1024]);

extern void runica_simple_initialize(void);

extern void runica_simple_terminate(void);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (runica_simple.h) */
