/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * runica_simple.c
 *
 * Code generation for function 'runica_simple'
 *
 */

/* Include files */
#include "runica_simple.h"
#include "rt_nonfinite.h"
#include "rt_defines.h"
#include "rt_nonfinite.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

/* Variable Definitions */
static unsigned int state[625];

static boolean_T isInitialized_runica_simple = false;

/* Function Declarations */
static void b_rand(double *r, int samples);

static void b_sqrt(creal_T *x);

static double b_xnrm2(int n, const double x[3]);

static void c_eml_rand_mt19937ar_stateful_i(void);

static void inv(const creal_T *x, creal_T *y, int nchan);

static void mean(const double *x, double *y, int nchan, int samples);

static void randperm(double *p, int samples);

static double rt_atan2d_snf(double u0, double u1);

static double rt_hypotd_snf(double u0, double u1);

static void sort(creal_T x[32], int idx[32]);

static boolean_T sortLE(const creal_T v[32], int idx1, int idx2);

static void sqrtm(const double *A, creal_T *X, int nchan);

static double xdlanv2(double *a, double *b, double *c, double *d, double *rt1i,
                      double *rt2r, double *rt2i, double *cs, double *sn);

static int xhseqr(double h[1024], double z[1024]);

static double xnrm2(int n, const double x[1024], int ix0);

static void xrot(int n, double x[1024], int ix0, int iy0, double c, double s);

static void xzlarf(int m, int n, int iv0, double tau, double C[1024], int ic0,
                   double work[32]);

/* Function Definitions */
static void b_rand(double *r, int samples)
{
  int j;
  int k;
  int kk;
  for (k = 0; k < samples; k++) {
    double b_r;
    /* ========================= COPYRIGHT NOTICE ============================
     */
    /*  This is a uniform (0,1) pseudorandom number generator based on: */
    /*                                                                         */
    /*  A C-program for MT19937, with initialization improved 2002/1/26. */
    /*  Coded by Takuji Nishimura and Makoto Matsumoto. */
    /*                                                                         */
    /*  Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura, */
    /*  All rights reserved. */
    /*                                                                         */
    /*  Redistribution and use in source and binary forms, with or without */
    /*  modification, are permitted provided that the following conditions */
    /*  are met: */
    /*                                                                         */
    /*    1. Redistributions of source code must retain the above copyright */
    /*       notice, this list of conditions and the following disclaimer. */
    /*                                                                         */
    /*    2. Redistributions in binary form must reproduce the above copyright
     */
    /*       notice, this list of conditions and the following disclaimer */
    /*       in the documentation and/or other materials provided with the */
    /*       distribution. */
    /*                                                                         */
    /*    3. The names of its contributors may not be used to endorse or */
    /*       promote products derived from this software without specific */
    /*       prior written permission. */
    /*                                                                         */
    /*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS */
    /*  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT */
    /*  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR */
    /*  A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT */
    /*  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, */
    /*  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT */
    /*  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, */
    /*  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY */
    /*  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT */
    /*  (INCLUDING  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
     */
    /*  OF THIS  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. */
    /*                                                                         */
    /* =============================   END   =================================
     */
    unsigned int u[2];
    do {
      for (j = 0; j < 2; j++) {
        unsigned int mti;
        unsigned int y;
        mti = state[624] + 1U;
        if (state[624] + 1U >= 625U) {
          for (kk = 0; kk < 227; kk++) {
            y = (state[kk] & 2147483648U) | (state[kk + 1] & 2147483647U);
            if ((y & 1U) == 0U) {
              y >>= 1U;
            } else {
              y = y >> 1U ^ 2567483615U;
            }
            state[kk] = state[kk + 397] ^ y;
          }
          for (kk = 0; kk < 396; kk++) {
            y = (state[kk + 227] & 2147483648U) |
                (state[kk + 228] & 2147483647U);
            if ((y & 1U) == 0U) {
              y >>= 1U;
            } else {
              y = y >> 1U ^ 2567483615U;
            }
            state[kk + 227] = state[kk] ^ y;
          }
          y = (state[623] & 2147483648U) | (state[0] & 2147483647U);
          if ((y & 1U) == 0U) {
            y >>= 1U;
          } else {
            y = y >> 1U ^ 2567483615U;
          }
          state[623] = state[396] ^ y;
          mti = 1U;
        }
        y = state[(int)mti - 1];
        state[624] = mti;
        y ^= y >> 11U;
        y ^= y << 7U & 2636928640U;
        y ^= y << 15U & 4022730752U;
        u[j] = y ^ y >> 18U;
      }
      u[0] >>= 5U;
      u[1] >>= 6U;
      b_r =
          1.1102230246251565E-16 * ((double)u[0] * 6.7108864E+7 + (double)u[1]);
    } while (b_r == 0.0);
    r[k] = b_r;
  }
}

static void b_sqrt(creal_T *x)
{
  double absxi;
  double absxr;
  double xi;
  double xr;
  xr = x->re;
  xi = x->im;
  if (xi == 0.0) {
    if (xr < 0.0) {
      absxr = 0.0;
      absxi = sqrt(-xr);
    } else {
      absxr = sqrt(xr);
      absxi = 0.0;
    }
  } else if (xr == 0.0) {
    if (xi < 0.0) {
      absxr = sqrt(-xi / 2.0);
      absxi = -absxr;
    } else {
      absxr = sqrt(xi / 2.0);
      absxi = absxr;
    }
  } else if (rtIsNaN(xr)) {
    absxr = rtNaN;
    absxi = rtNaN;
  } else if (rtIsNaN(xi)) {
    absxr = rtNaN;
    absxi = rtNaN;
  } else if (rtIsInf(xi)) {
    absxr = fabs(xi);
    absxi = xi;
  } else if (rtIsInf(xr)) {
    if (xr < 0.0) {
      absxr = 0.0;
      absxi = xi * -xr;
    } else {
      absxr = xr;
      absxi = 0.0;
    }
  } else {
    absxr = fabs(xr);
    absxi = fabs(xi);
    if ((absxr > 4.4942328371557893E+307) ||
        (absxi > 4.4942328371557893E+307)) {
      absxr *= 0.5;
      absxi = rt_hypotd_snf(absxr, absxi * 0.5);
      if (absxi > absxr) {
        absxr = sqrt(absxi) * sqrt(absxr / absxi + 1.0);
      } else {
        absxr = sqrt(absxi) * 1.4142135623730951;
      }
    } else {
      absxr = sqrt((rt_hypotd_snf(absxr, absxi) + absxr) * 0.5);
    }
    if (xr > 0.0) {
      absxi = 0.5 * (xi / absxr);
    } else {
      if (xi < 0.0) {
        absxi = -absxr;
      } else {
        absxi = absxr;
      }
      absxr = 0.5 * (xi / absxi);
    }
  }
  x->re = absxr;
  x->im = absxi;
}

static double b_xnrm2(int n, const double x[3])
{
  double y;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = fabs(x[1]);
    } else {
      double absxk;
      double scale;
      double t;
      scale = 3.3121686421112381E-170;
      absxk = fabs(x[1]);
      if (absxk > 3.3121686421112381E-170) {
        y = 1.0;
        scale = absxk;
      } else {
        t = absxk / 3.3121686421112381E-170;
        y = t * t;
      }
      absxk = fabs(x[2]);
      if (absxk > scale) {
        t = scale / absxk;
        y = y * t * t + 1.0;
        scale = absxk;
      } else {
        t = absxk / scale;
        y += t * t;
      }
      y = scale * sqrt(y);
    }
  }
  return y;
}

static void c_eml_rand_mt19937ar_stateful_i(void)
{
  int mti;
  unsigned int r;
  memset(&state[0], 0, 625U * sizeof(unsigned int));
  r = 5489U;
  state[0] = 5489U;
  for (mti = 0; mti < 623; mti++) {
    r = ((r ^ r >> 30U) * 1812433253U + (unsigned int)mti) + 1U;
    state[mti + 1] = r;
  }
  state[624] = 624U;
}

static void inv(const creal_T *x, creal_T *y, int nchan)
{
  /* TODO: Generalize to variable nchan - currently only supports nchan=32 */
  if (nchan != 32) {
    fprintf(stderr, "Error: inv() currently only supports nchan=32, got %d\n", nchan);
    return;
  }
  creal_T b_x[1024];
  double ai;
  double ar;
  double bi;
  double re;
  double s;
  double sgnbr;
  double smax;
  double temp_re;
  int b_i;
  int i;
  int i2;
  int j;
  int jA;
  int jBcol;
  int jp1j;
  int k;
  int pipk;
  signed char ipiv[32];
  signed char p[32];
  for (i = 0; i < 1024; i++) {
    y[i].re = 0.0;
    y[i].im = 0.0;
    b_x[i] = x[i];
  }
  for (i = 0; i < 32; i++) {
    ipiv[i] = (signed char)(i + 1);
  }
  for (j = 0; j < 31; j++) {
    int b_tmp;
    int mmj_tmp;
    mmj_tmp = 30 - j;
    b_tmp = j * 33;
    jp1j = b_tmp + 2;
    pipk = 32 - j;
    jA = 0;
    smax = fabs(b_x[b_tmp].re) + fabs(b_x[b_tmp].im);
    for (k = 2; k <= pipk; k++) {
      jBcol = (b_tmp + k) - 1;
      s = fabs(b_x[jBcol].re) + fabs(b_x[jBcol].im);
      if (s > smax) {
        jA = k - 1;
        smax = s;
      }
    }
    pipk = b_tmp + jA;
    if ((b_x[pipk].re != 0.0) || (b_x[pipk].im != 0.0)) {
      if (jA != 0) {
        jA += j;
        ipiv[j] = (signed char)(jA + 1);
        for (k = 0; k < 32; k++) {
          pipk = k << 5;
          jBcol = j + pipk;
          temp_re = b_x[jBcol].re;
          smax = b_x[jBcol].im;
          b_x[jBcol] = b_x[jA + pipk];
          i = jA + (k << 5);
          b_x[i].re = temp_re;
          b_x[i].im = smax;
        }
      }
      i = (b_tmp - j) + 32;
      for (b_i = jp1j; b_i <= i; b_i++) {
        ar = b_x[b_i - 1].re;
        ai = b_x[b_i - 1].im;
        temp_re = b_x[b_tmp].re;
        bi = b_x[b_tmp].im;
        if (bi == 0.0) {
          if (ai == 0.0) {
            re = ar / temp_re;
            smax = 0.0;
          } else if (ar == 0.0) {
            re = 0.0;
            smax = ai / temp_re;
          } else {
            re = ar / temp_re;
            smax = ai / temp_re;
          }
        } else if (temp_re == 0.0) {
          if (ar == 0.0) {
            re = ai / bi;
            smax = 0.0;
          } else if (ai == 0.0) {
            re = 0.0;
            smax = -(ar / bi);
          } else {
            re = ai / bi;
            smax = -(ar / bi);
          }
        } else {
          s = fabs(temp_re);
          smax = fabs(bi);
          if (s > smax) {
            s = bi / temp_re;
            smax = temp_re + s * bi;
            re = (ar + s * ai) / smax;
            smax = (ai - s * ar) / smax;
          } else if (smax == s) {
            if (temp_re > 0.0) {
              sgnbr = 0.5;
            } else {
              sgnbr = -0.5;
            }
            if (bi > 0.0) {
              smax = 0.5;
            } else {
              smax = -0.5;
            }
            re = (ar * sgnbr + ai * smax) / s;
            smax = (ai * sgnbr - ar * smax) / s;
          } else {
            s = temp_re / bi;
            smax = bi + s * temp_re;
            re = (s * ar + ai) / smax;
            smax = (s * ai - ar) / smax;
          }
        }
        b_x[b_i - 1].re = re;
        b_x[b_i - 1].im = smax;
      }
    }
    jA = b_tmp;
    for (jp1j = 0; jp1j <= mmj_tmp; jp1j++) {
      pipk = (b_tmp + (jp1j << 5)) + 32;
      smax = b_x[pipk].re;
      sgnbr = b_x[pipk].im;
      if ((smax != 0.0) || (sgnbr != 0.0)) {
        temp_re = -smax - sgnbr * 0.0;
        smax = smax * 0.0 - sgnbr;
        i = jA + 34;
        i2 = (jA - j) + 64;
        for (pipk = i; pipk <= i2; pipk++) {
          jBcol = ((b_tmp + pipk) - jA) - 33;
          sgnbr = b_x[jBcol].re;
          s = b_x[jBcol].im;
          b_x[pipk - 1].re += sgnbr * temp_re - s * smax;
          b_x[pipk - 1].im += sgnbr * smax + s * temp_re;
        }
      }
      jA += 32;
    }
  }
  for (i = 0; i < 32; i++) {
    p[i] = (signed char)(i + 1);
  }
  for (k = 0; k < 31; k++) {
    signed char i1;
    i1 = ipiv[k];
    if (i1 > k + 1) {
      pipk = p[i1 - 1];
      p[i1 - 1] = p[k];
      p[k] = (signed char)pipk;
    }
  }
  for (k = 0; k < 32; k++) {
    i = (p[k] - 1) << 5;
    i2 = k + i;
    y[i2].re = 1.0;
    y[i2].im = 0.0;
    for (j = k + 1; j < 33; j++) {
      jBcol = (j + i) - 1;
      if ((y[jBcol].re != 0.0) || (y[jBcol].im != 0.0)) {
        i2 = j + 1;
        for (b_i = i2; b_i < 33; b_i++) {
          smax = y[jBcol].re;
          pipk = (b_i + ((j - 1) << 5)) - 1;
          sgnbr = b_x[pipk].im;
          s = y[jBcol].im;
          temp_re = b_x[pipk].re;
          pipk = (b_i + i) - 1;
          y[pipk].re -= smax * temp_re - s * sgnbr;
          y[pipk].im -= smax * sgnbr + s * temp_re;
        }
      }
    }
  }
  for (j = 0; j < 32; j++) {
    jBcol = j << 5;
    for (k = 31; k >= 0; k--) {
      jA = k << 5;
      jp1j = k + jBcol;
      ar = y[jp1j].re;
      ai = y[jp1j].im;
      if ((ar != 0.0) || (ai != 0.0)) {
        pipk = k + jA;
        temp_re = b_x[pipk].re;
        bi = b_x[pipk].im;
        if (bi == 0.0) {
          if (ai == 0.0) {
            re = ar / temp_re;
            smax = 0.0;
          } else if (ar == 0.0) {
            re = 0.0;
            smax = ai / temp_re;
          } else {
            re = ar / temp_re;
            smax = ai / temp_re;
          }
        } else if (temp_re == 0.0) {
          if (ar == 0.0) {
            re = ai / bi;
            smax = 0.0;
          } else if (ai == 0.0) {
            re = 0.0;
            smax = -(ar / bi);
          } else {
            re = ai / bi;
            smax = -(ar / bi);
          }
        } else {
          s = fabs(temp_re);
          smax = fabs(bi);
          if (s > smax) {
            s = bi / temp_re;
            smax = temp_re + s * bi;
            re = (ar + s * ai) / smax;
            smax = (ai - s * ar) / smax;
          } else if (smax == s) {
            if (temp_re > 0.0) {
              sgnbr = 0.5;
            } else {
              sgnbr = -0.5;
            }
            if (bi > 0.0) {
              smax = 0.5;
            } else {
              smax = -0.5;
            }
            re = (ar * sgnbr + ai * smax) / s;
            smax = (ai * sgnbr - ar * smax) / s;
          } else {
            s = temp_re / bi;
            smax = bi + s * temp_re;
            re = (s * ar + ai) / smax;
            smax = (s * ai - ar) / smax;
          }
        }
        y[jp1j].re = re;
        y[jp1j].im = smax;
        for (b_i = 0; b_i < k; b_i++) {
          smax = y[jp1j].re;
          pipk = b_i + jA;
          sgnbr = b_x[pipk].im;
          s = y[jp1j].im;
          temp_re = b_x[pipk].re;
          i = b_i + jBcol;
          y[i].re -= smax * temp_re - s * sgnbr;
          y[i].im -= smax * sgnbr + s * temp_re;
        }
      }
    }
  }
}

static void mean(const double *x, double *y, int nchan, int samples)
{
  double *bsum = (double *)malloc(nchan * sizeof(double));
  int ib;
  int k;
  int xj;
  int xoffset;

  if (!bsum) {
    fprintf(stderr, "Error: Failed to allocate memory in mean\n");
    return;
  }

  /* Simple implementation: sum across samples dimension */
  for (xj = 0; xj < nchan; xj++) {
    y[xj] = 0.0;
  }

  for (k = 0; k < samples; k++) {
    xoffset = k * nchan;
    for (xj = 0; xj < nchan; xj++) {
      y[xj] += x[xoffset + xj];
    }
  }

  for (k = 0; k < nchan; k++) {
    y[k] /= (double)samples;
  }

  free(bsum);
}

static void randperm(double *p, int samples)
{
  int *idx = (int *)malloc(samples * sizeof(int));
  int *iwork = (int *)malloc(samples * sizeof(int));
  double d;
  int b_i;
  int i;
  int k;
  int pEnd;
  int qEnd;

  if (!idx || !iwork) {
    free(idx);
    free(iwork);
    fprintf(stderr, "Error: Failed to allocate memory in randperm\n");
    return;
  }

  b_rand(p, samples);

  /* Handle odd samples case */
  int samples_even = (samples / 2) * 2;
  for (k = 0; k < samples_even; k += 2) {
    d = p[k + 1];
    if ((p[k] <= d) || rtIsNaN(d)) {
      idx[k] = k + 1;
      idx[k + 1] = k + 2;
    } else {
      idx[k] = k + 2;
      idx[k + 1] = k + 1;
    }
  }
  /* Handle remaining odd sample */
  if (samples % 2 == 1) {
    idx[samples - 1] = samples;
  }

  i = 2;
  while (i < samples) {
    int i2;
    int j;
    i2 = i << 1;
    j = 1;
    for (pEnd = i + 1; pEnd <= samples + 1; pEnd = qEnd + i) {
      int b_p;
      int kEnd;
      int q;
      b_p = j;
      q = pEnd - 1;
      qEnd = j + i2;
      if (qEnd > samples + 1) {
        qEnd = samples + 1;
      }
      k = 0;
      kEnd = qEnd - j;
      while (k + 1 <= kEnd) {
        d = p[idx[q] - 1];
        b_i = idx[b_p - 1];
        if ((p[b_i - 1] <= d) || rtIsNaN(d)) {
          iwork[k] = b_i;
          b_p++;
          if (b_p == pEnd) {
            while (q + 1 < qEnd) {
              k++;
              iwork[k] = idx[q];
              q++;
            }
          }
        } else {
          iwork[k] = idx[q];
          q++;
          if (q + 1 == qEnd) {
            while (b_p < pEnd) {
              k++;
              iwork[k] = idx[b_p - 1];
              b_p++;
            }
          }
        }
        k++;
      }
      for (k = 0; k < kEnd; k++) {
        idx[(j + k) - 1] = iwork[k];
      }
      j = qEnd;
    }
    i = i2;
  }
  for (b_i = 0; b_i < samples; b_i++) {
    p[b_i] = idx[b_i];
  }

  free(idx);
  free(iwork);
}

static double rt_atan2d_snf(double u0, double u1)
{
  double y;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else if (rtIsInf(u0) && rtIsInf(u1)) {
    int i;
    int i1;
    if (u0 > 0.0) {
      i = 1;
    } else {
      i = -1;
    }
    if (u1 > 0.0) {
      i1 = 1;
    } else {
      i1 = -1;
    }
    y = atan2(i, i1);
  } else if (u1 == 0.0) {
    if (u0 > 0.0) {
      y = RT_PI / 2.0;
    } else if (u0 < 0.0) {
      y = -(RT_PI / 2.0);
    } else {
      y = 0.0;
    }
  } else {
    y = atan2(u0, u1);
  }
  return y;
}

static double rt_hypotd_snf(double u0, double u1)
{
  double a;
  double b;
  double y;
  a = fabs(u0);
  b = fabs(u1);
  if (a < b) {
    a /= b;
    y = b * sqrt(a * a + 1.0);
  } else if (a > b) {
    b /= a;
    y = a * sqrt(b * b + 1.0);
  } else if (rtIsNaN(b)) {
    y = rtNaN;
  } else {
    y = a * 1.4142135623730951;
  }
  return y;
}

static void sort(creal_T x[32], int idx[32])
{
  creal_T xwork[32];
  int iwork[32];
  int i;
  int k;
  int pEnd;
  int qEnd;
  for (k = 0; k <= 30; k += 2) {
    if (sortLE(x, k + 1, k + 2)) {
      idx[k] = k + 1;
      idx[k + 1] = k + 2;
    } else {
      idx[k] = k + 2;
      idx[k + 1] = k + 1;
    }
  }
  i = 2;
  while (i < 32) {
    int i2;
    int j;
    i2 = i << 1;
    j = 1;
    for (pEnd = i + 1; pEnd < 33; pEnd = qEnd + i) {
      int kEnd;
      int p;
      int q;
      p = j;
      q = pEnd;
      qEnd = j + i2;
      k = 0;
      kEnd = qEnd - j;
      while (k + 1 <= kEnd) {
        int b_i;
        int i1;
        b_i = idx[p - 1];
        i1 = idx[q - 1];
        if (sortLE(x, b_i, i1)) {
          iwork[k] = b_i;
          p++;
          if (p == pEnd) {
            while (q < qEnd) {
              k++;
              iwork[k] = idx[q - 1];
              q++;
            }
          }
        } else {
          iwork[k] = i1;
          q++;
          if (q == qEnd) {
            while (p < pEnd) {
              k++;
              iwork[k] = idx[p - 1];
              p++;
            }
          }
        }
        k++;
      }
      for (k = 0; k < kEnd; k++) {
        idx[(j + k) - 1] = iwork[k];
      }
      j = qEnd;
    }
    i = i2;
  }
  memcpy(&xwork[0], &x[0], 32U * sizeof(creal_T));
  for (k = 0; k < 32; k++) {
    x[k] = xwork[idx[k] - 1];
  }
}

static boolean_T sortLE(const creal_T v[32], int idx1, int idx2)
{
  boolean_T p;
  if (rtIsNaN(v[idx2 - 1].re) || rtIsNaN(v[idx2 - 1].im)) {
    p = (rtIsNaN(v[idx1 - 1].re) || rtIsNaN(v[idx1 - 1].im));
  } else if (rtIsNaN(v[idx1 - 1].re) || rtIsNaN(v[idx1 - 1].im)) {
    p = true;
  } else {
    double x;
    double y;
    boolean_T SCALEA;
    boolean_T SCALEB;
    if ((fabs(v[idx1 - 1].re) > 8.9884656743115785E+307) ||
        (fabs(v[idx1 - 1].im) > 8.9884656743115785E+307)) {
      SCALEA = true;
    } else {
      SCALEA = false;
    }
    if ((fabs(v[idx2 - 1].re) > 8.9884656743115785E+307) ||
        (fabs(v[idx2 - 1].im) > 8.9884656743115785E+307)) {
      SCALEB = true;
    } else {
      SCALEB = false;
    }
    if (SCALEA || SCALEB) {
      x = rt_hypotd_snf(v[idx1 - 1].re / 2.0, v[idx1 - 1].im / 2.0);
      y = rt_hypotd_snf(v[idx2 - 1].re / 2.0, v[idx2 - 1].im / 2.0);
    } else {
      x = rt_hypotd_snf(v[idx1 - 1].re, v[idx1 - 1].im);
      y = rt_hypotd_snf(v[idx2 - 1].re, v[idx2 - 1].im);
    }
    if (x == y) {
      double b_x_tmp;
      double b_y_tmp;
      double x_tmp;
      double y_tmp;
      x_tmp = v[idx1 - 1].re;
      b_x_tmp = v[idx1 - 1].im;
      x = rt_atan2d_snf(b_x_tmp, x_tmp);
      y_tmp = v[idx2 - 1].re;
      b_y_tmp = v[idx2 - 1].im;
      y = rt_atan2d_snf(b_y_tmp, y_tmp);
      if (x == y) {
        if (x_tmp != y_tmp) {
          if (x >= 0.0) {
            x = y_tmp;
            b_y_tmp = x_tmp;
          } else {
            x = x_tmp;
            b_y_tmp = y_tmp;
          }
        } else if (x_tmp < 0.0) {
          x = b_y_tmp;
          b_y_tmp = b_x_tmp;
        } else {
          x = b_x_tmp;
        }
        y = b_y_tmp;
        if (x == b_y_tmp) {
          x = 0.0;
          y = 0.0;
        }
      }
    }
    p = (x >= y);
  }
  return p;
}

static void sqrtm(const double *A, creal_T *X, int nchan)
{
  /* TODO: Generalize to variable nchan - currently only supports nchan=32 */
  if (nchan != 32) {
    fprintf(stderr, "Error: sqrtm() currently only supports nchan=32, got %d\n", nchan);
    return;
  }
  creal_T Q[1024];
  creal_T R[1024];
  creal_T T[1024];
  double Vr[1024];
  double b_A[1024];
  double work[32];
  double a;
  double b_d;
  double d;
  double d_tmp;
  double mu1_im;
  double mu1_re;
  double rt1i;
  double rt1r;
  double s;
  double sn;
  double t1_im;
  double xnorm;
  int b_i;
  int exitg1;
  int i;
  int i1;
  int ic0;
  int ix0;
  int j;
  int knt;
  int m;
  int mm1;
  boolean_T exitg2;
  boolean_T p;
  p = true;
  for (m = 0; m < 1024; m++) {
    if (p) {
      d = A[m];
      if (rtIsInf(d) || rtIsNaN(d)) {
        p = false;
      }
    } else {
      p = false;
    }
  }
  if (!p) {
    for (i = 0; i < 1024; i++) {
      Q[i].re = rtNaN;
      Q[i].im = 0.0;
    }
    knt = 3;
    for (j = 0; j < 30; j++) {
      if (knt <= 32) {
        memset(&Q[(j * 32 + knt) + -1], 0,
               (unsigned int)(-knt + 33) * sizeof(creal_T));
      }
      knt++;
    }
    for (i = 0; i < 1024; i++) {
      T[i].re = rtNaN;
      T[i].im = 0.0;
    }
  } else {
    double tau[31];
    memcpy(&b_A[0], &A[0], 1024U * sizeof(double));
    memset(&work[0], 0, 32U * sizeof(double));
    for (b_i = 0; b_i < 31; b_i++) {
      int alpha1_tmp_tmp_tmp;
      int in;
      int lastc;
      int lastv;
      knt = b_i << 5;
      in = (b_i + 1) << 5;
      alpha1_tmp_tmp_tmp = b_i + knt;
      t1_im = b_A[alpha1_tmp_tmp_tmp + 1];
      if (b_i + 3 <= 32) {
        i = b_i + 1;
      } else {
        i = 30;
      }
      ix0 = (i + knt) + 2;
      tau[b_i] = 0.0;
      xnorm = xnrm2(30 - b_i, b_A, ix0);
      if (xnorm != 0.0) {
        xnorm = rt_hypotd_snf(t1_im, xnorm);
        if (t1_im >= 0.0) {
          xnorm = -xnorm;
        }
        if (fabs(xnorm) < 1.0020841800044864E-292) {
          knt = 0;
          i = (ix0 - b_i) + 29;
          do {
            knt++;
            for (m = ix0; m <= i; m++) {
              b_A[m - 1] *= 9.9792015476736E+291;
            }
            xnorm *= 9.9792015476736E+291;
            t1_im *= 9.9792015476736E+291;
          } while ((fabs(xnorm) < 1.0020841800044864E-292) && (knt < 20));
          xnorm = rt_hypotd_snf(t1_im, xnrm2(30 - b_i, b_A, ix0));
          if (t1_im >= 0.0) {
            xnorm = -xnorm;
          }
          tau[b_i] = (xnorm - t1_im) / xnorm;
          a = 1.0 / (t1_im - xnorm);
          for (m = ix0; m <= i; m++) {
            b_A[m - 1] *= a;
          }
          for (m = 0; m < knt; m++) {
            xnorm *= 1.0020841800044864E-292;
          }
          t1_im = xnorm;
        } else {
          tau[b_i] = (xnorm - t1_im) / xnorm;
          a = 1.0 / (t1_im - xnorm);
          i = (ix0 - b_i) + 29;
          for (m = ix0; m <= i; m++) {
            b_A[m - 1] *= a;
          }
          t1_im = xnorm;
        }
      }
      b_A[alpha1_tmp_tmp_tmp + 1] = 1.0;
      ic0 = in + 1;
      if (tau[b_i] != 0.0) {
        lastv = 30 - b_i;
        knt = (alpha1_tmp_tmp_tmp - b_i) + 31;
        while ((lastv + 1 > 0) && (b_A[knt] == 0.0)) {
          lastv--;
          knt--;
        }
        lastc = 32;
        exitg2 = false;
        while ((!exitg2) && (lastc > 0)) {
          knt = in + lastc;
          m = knt;
          do {
            exitg1 = 0;
            if (m <= knt + (lastv << 5)) {
              if (b_A[m - 1] != 0.0) {
                exitg1 = 1;
              } else {
                m += 32;
              }
            } else {
              lastc--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);
          if (exitg1 == 1) {
            exitg2 = true;
          }
        }
      } else {
        lastv = -1;
        lastc = 0;
      }
      if (lastv + 1 > 0) {
        if (lastc != 0) {
          memset(&work[0], 0, (unsigned int)lastc * sizeof(double));
          knt = alpha1_tmp_tmp_tmp + 1;
          i = (in + (lastv << 5)) + 1;
          for (ix0 = ic0; ix0 <= i; ix0 += 32) {
            i1 = (ix0 + lastc) - 1;
            for (m = ix0; m <= i1; m++) {
              mm1 = m - ix0;
              work[mm1] += b_A[m - 1] * b_A[knt];
            }
            knt++;
          }
        }
        if (!(-tau[b_i] == 0.0)) {
          knt = in;
          for (j = 0; j <= lastv; j++) {
            d = b_A[(alpha1_tmp_tmp_tmp + j) + 1];
            if (d != 0.0) {
              xnorm = d * -tau[b_i];
              i = knt + 1;
              i1 = lastc + knt;
              for (ix0 = i; ix0 <= i1; ix0++) {
                b_A[ix0 - 1] += work[(ix0 - knt) - 1] * xnorm;
              }
            }
            knt += 32;
          }
        }
      }
      xzlarf(31 - b_i, 31 - b_i, alpha1_tmp_tmp_tmp + 2, tau[b_i], b_A,
             (b_i + in) + 2, work);
      b_A[alpha1_tmp_tmp_tmp + 1] = t1_im;
    }
    memcpy(&Vr[0], &b_A[0], 1024U * sizeof(double));
    for (j = 30; j >= 0; j--) {
      m = (j + 1) << 5;
      for (b_i = 0; b_i <= j; b_i++) {
        Vr[m + b_i] = 0.0;
      }
      i = j + 3;
      for (b_i = i; b_i < 33; b_i++) {
        knt = m + b_i;
        Vr[knt - 1] = Vr[knt - 33];
      }
    }
    memset(&Vr[0], 0, 32U * sizeof(double));
    Vr[0] = 1.0;
    memset(&work[0], 0, 32U * sizeof(double));
    for (b_i = 30; b_i >= 0; b_i--) {
      knt = (b_i + (b_i << 5)) + 33;
      if (b_i + 1 < 31) {
        Vr[knt] = 1.0;
        xzlarf(31 - b_i, 30 - b_i, knt + 1, tau[b_i], Vr, knt + 33, work);
        ix0 = knt + 2;
        i = (knt - b_i) + 31;
        for (m = ix0; m <= i; m++) {
          Vr[m - 1] *= -tau[b_i];
        }
      }
      Vr[knt] = 1.0 - tau[b_i];
      for (j = 0; j < b_i; j++) {
        Vr[(knt - j) - 1] = 0.0;
      }
    }
    xhseqr(b_A, Vr);
    for (i = 0; i < 1024; i++) {
      T[i].re = b_A[i];
      T[i].im = 0.0;
      Q[i].re = Vr[i];
      Q[i].im = 0.0;
    }
    for (m = 30; m >= 0; m--) {
      mm1 = m + 1;
      i = m << 5;
      i1 = m + i;
      d = b_A[i1 + 1];
      if (d != 0.0) {
        a = b_A[i1];
        ic0 = (m + 1) << 5;
        knt = m + ic0;
        xnorm = b_A[knt];
        t1_im = d;
        d_tmp = b_A[knt + 1];
        b_d = d_tmp;
        rt1r =
            xdlanv2(&a, &xnorm, &t1_im, &b_d, &rt1i, &mu1_im, &mu1_re, &s, &sn);
        mu1_re = rt1r - d_tmp;
        xnorm = rt_hypotd_snf(rt_hypotd_snf(mu1_re, rt1i), d);
        if (rt1i == 0.0) {
          a = mu1_re / xnorm;
          sn = 0.0;
        } else if (mu1_re == 0.0) {
          a = 0.0;
          sn = rt1i / xnorm;
        } else {
          a = mu1_re / xnorm;
          sn = rt1i / xnorm;
        }
        s = d / xnorm;
        for (j = mm1; j < 33; j++) {
          ix0 = m + ((j - 1) << 5);
          xnorm = T[ix0].re;
          t1_im = T[ix0].im;
          d_tmp = T[ix0 + 1].re;
          rt1r = T[ix0 + 1].im;
          T[ix0].re = (a * xnorm + sn * t1_im) + s * d_tmp;
          T[ix0].im = (a * t1_im - sn * xnorm) + s * rt1r;
          mu1_re = a * d_tmp - sn * rt1r;
          mu1_im = a * rt1r + sn * d_tmp;
          T[ix0 + 1].re = mu1_re - s * xnorm;
          T[ix0 + 1].im = mu1_im - s * t1_im;
        }
        for (b_i = 0; b_i <= m + 1; b_i++) {
          ix0 = b_i + i;
          xnorm = T[ix0].re;
          t1_im = T[ix0].im;
          knt = b_i + ic0;
          d_tmp = T[knt].re;
          rt1r = T[knt].im;
          mu1_re = a * xnorm - sn * t1_im;
          mu1_im = a * t1_im + sn * xnorm;
          T[ix0].re = mu1_re + s * d_tmp;
          T[ix0].im = mu1_im + s * rt1r;
          T[knt].re = (a * d_tmp + sn * rt1r) - s * xnorm;
          T[knt].im = (a * rt1r - sn * d_tmp) - s * t1_im;
        }
        for (b_i = 0; b_i < 32; b_i++) {
          ix0 = b_i + i;
          xnorm = Q[ix0].re;
          t1_im = Q[ix0].im;
          knt = b_i + ic0;
          d_tmp = Q[knt].re;
          rt1r = Q[knt].im;
          mu1_re = a * xnorm - sn * t1_im;
          mu1_im = a * t1_im + sn * xnorm;
          Q[ix0].re = mu1_re + s * d_tmp;
          Q[ix0].im = mu1_im + s * rt1r;
          Q[knt].re = (a * d_tmp + sn * rt1r) - s * xnorm;
          Q[knt].im = (a * rt1r - sn * d_tmp) - s * t1_im;
        }
        T[(m + (m << 5)) + 1].re = 0.0;
        T[(m + (m << 5)) + 1].im = 0.0;
      }
    }
  }
  j = -1;
  int exitg3;
  do {
    exitg3 = 0;
    if (j + 1 < 32) {
      b_i = 0;
      do {
        exitg1 = 0;
        if (b_i <= j) {
          knt = b_i + ((j + 1) << 5);
          if ((T[knt].re != 0.0) || (T[knt].im != 0.0)) {
            memset(&R[0], 0, 1024U * sizeof(creal_T));
            for (j = 0; j < 32; j++) {
              mm1 = j << 5;
              ic0 = j + mm1;
              R[ic0] = T[ic0];
              b_sqrt(&R[ic0]);
              for (b_i = j; b_i >= 1; b_i--) {
                mu1_re = 0.0;
                mu1_im = 0.0;
                i = b_i + 1;
                for (m = i; m <= j; m++) {
                  knt = (b_i + ((m - 1) << 5)) - 1;
                  xnorm = R[knt].re;
                  ix0 = (m + mm1) - 1;
                  t1_im = R[ix0].im;
                  d_tmp = R[knt].im;
                  rt1r = R[ix0].re;
                  mu1_re += xnorm * rt1r - d_tmp * t1_im;
                  mu1_im += xnorm * t1_im + d_tmp * rt1r;
                }
                knt = (b_i + mm1) - 1;
                mu1_re = T[knt].re - mu1_re;
                mu1_im = T[knt].im - mu1_im;
                ix0 = (b_i + ((b_i - 1) << 5)) - 1;
                t1_im = R[ix0].re + R[ic0].re;
                d_tmp = R[ix0].im + R[ic0].im;
                if (d_tmp == 0.0) {
                  if (mu1_im == 0.0) {
                    a = mu1_re / t1_im;
                    sn = 0.0;
                  } else if (mu1_re == 0.0) {
                    a = 0.0;
                    sn = mu1_im / t1_im;
                  } else {
                    a = mu1_re / t1_im;
                    sn = mu1_im / t1_im;
                  }
                } else if (t1_im == 0.0) {
                  if (mu1_re == 0.0) {
                    a = mu1_im / d_tmp;
                    sn = 0.0;
                  } else if (mu1_im == 0.0) {
                    a = 0.0;
                    sn = -(mu1_re / d_tmp);
                  } else {
                    a = mu1_im / d_tmp;
                    sn = -(mu1_re / d_tmp);
                  }
                } else {
                  rt1r = fabs(t1_im);
                  xnorm = fabs(d_tmp);
                  if (rt1r > xnorm) {
                    s = d_tmp / t1_im;
                    b_d = t1_im + s * d_tmp;
                    a = (mu1_re + s * mu1_im) / b_d;
                    sn = (mu1_im - s * mu1_re) / b_d;
                  } else if (xnorm == rt1r) {
                    if (t1_im > 0.0) {
                      t1_im = 0.5;
                    } else {
                      t1_im = -0.5;
                    }
                    if (d_tmp > 0.0) {
                      xnorm = 0.5;
                    } else {
                      xnorm = -0.5;
                    }
                    a = (mu1_re * t1_im + mu1_im * xnorm) / rt1r;
                    sn = (mu1_im * t1_im - mu1_re * xnorm) / rt1r;
                  } else {
                    s = t1_im / d_tmp;
                    b_d = d_tmp + s * t1_im;
                    a = (s * mu1_re + mu1_im) / b_d;
                    sn = (s * mu1_im - mu1_re) / b_d;
                  }
                }
                R[knt].re = a;
                R[knt].im = sn;
              }
            }
            exitg1 = 1;
          } else {
            b_i++;
          }
        } else {
          j++;
          exitg1 = 2;
        }
      } while (exitg1 == 0);
      if (exitg1 == 1) {
        exitg3 = 1;
      }
    } else {
      memset(&R[0], 0, 1024U * sizeof(creal_T));
      for (j = 0; j < 32; j++) {
        knt = j + (j << 5);
        R[knt] = T[knt];
        b_sqrt(&R[knt]);
      }
      exitg3 = 1;
    }
  } while (exitg3 == 0);
  for (i = 0; i < 32; i++) {
    for (i1 = 0; i1 < 32; i1++) {
      a = 0.0;
      sn = 0.0;
      for (mm1 = 0; mm1 < 32; mm1++) {
        ix0 = i + (mm1 << 5);
        xnorm = Q[ix0].re;
        knt = mm1 + (i1 << 5);
        t1_im = R[knt].im;
        d_tmp = Q[ix0].im;
        rt1r = R[knt].re;
        a += xnorm * rt1r - d_tmp * t1_im;
        sn += xnorm * t1_im + d_tmp * rt1r;
      }
      mm1 = i + (i1 << 5);
      T[mm1].re = a;
      T[mm1].im = sn;
    }
    for (i1 = 0; i1 < 32; i1++) {
      a = 0.0;
      sn = 0.0;
      for (mm1 = 0; mm1 < 32; mm1++) {
        ix0 = mm1 << 5;
        knt = i1 + ix0;
        xnorm = Q[knt].re;
        t1_im = -Q[knt].im;
        knt = i + ix0;
        d_tmp = T[knt].re;
        rt1r = T[knt].im;
        a += d_tmp * xnorm - rt1r * t1_im;
        sn += d_tmp * t1_im + rt1r * xnorm;
      }
      mm1 = i + (i1 << 5);
      X[mm1].re = a;
      X[mm1].im = sn;
    }
  }
  for (i = 0; i < 1024; i++) {
    b_A[i] = X[i].im;
  }
  xnorm = 0.0;
  j = 0;
  exitg2 = false;
  while ((!exitg2) && (j < 32)) {
    s = 0.0;
    for (b_i = 0; b_i < 32; b_i++) {
      s += fabs(b_A[b_i + (j << 5)]);
    }
    if (rtIsNaN(s)) {
      xnorm = rtNaN;
      exitg2 = true;
    } else {
      if (s > xnorm) {
        xnorm = s;
      }
      j++;
    }
  }
  t1_im = 0.0;
  j = 0;
  exitg2 = false;
  while ((!exitg2) && (j < 32)) {
    s = 0.0;
    for (b_i = 0; b_i < 32; b_i++) {
      knt = b_i + (j << 5);
      s += rt_hypotd_snf(X[knt].re, X[knt].im);
    }
    if (rtIsNaN(s)) {
      t1_im = rtNaN;
      exitg2 = true;
    } else {
      if (s > t1_im) {
        t1_im = s;
      }
      j++;
    }
  }
  if (xnorm <= 7.1054273576010019E-14 * t1_im) {
    for (i = 0; i < 1024; i++) {
      X[i].im = 0.0;
    }
  }
}

static double xdlanv2(double *a, double *b, double *c, double *d, double *rt1i,
                      double *rt2r, double *rt2i, double *cs, double *sn)
{
  double rt1r;
  if (*c == 0.0) {
    *cs = 1.0;
    *sn = 0.0;
  } else if (*b == 0.0) {
    double temp;
    *cs = 0.0;
    *sn = 1.0;
    temp = *d;
    *d = *a;
    *a = temp;
    *b = -*c;
    *c = 0.0;
  } else {
    double temp;
    temp = *a - *d;
    if ((temp == 0.0) && ((*b < 0.0) != (*c < 0.0))) {
      *cs = 1.0;
      *sn = 0.0;
    } else {
      double bcmax;
      double bcmis;
      double p;
      double scale;
      double z;
      int count;
      int i;
      p = 0.5 * temp;
      bcmis = fabs(*b);
      scale = fabs(*c);
      bcmax = fmax(bcmis, scale);
      if (!(*b < 0.0)) {
        count = 1;
      } else {
        count = -1;
      }
      if (!(*c < 0.0)) {
        i = 1;
      } else {
        i = -1;
      }
      bcmis = fmin(bcmis, scale) * (double)count * (double)i;
      scale = fmax(fabs(p), bcmax);
      z = p / scale * p + bcmax / scale * bcmis;
      if (z >= 8.8817841970012523E-16) {
        double tau;
        *a = sqrt(scale) * sqrt(z);
        if (p < 0.0) {
          *a = -*a;
        }
        z = p + *a;
        *a = *d + z;
        *d -= bcmax / z * bcmis;
        tau = rt_hypotd_snf(*c, z);
        *cs = z / tau;
        *sn = *c / tau;
        *b -= *c;
        *c = 0.0;
      } else {
        double tau;
        bcmis = *b + *c;
        scale = fmax(fabs(temp), fabs(bcmis));
        count = 0;
        while ((scale >= 7.4428285367870146E+137) && (count <= 20)) {
          bcmis *= 1.3435752215134178E-138;
          temp *= 1.3435752215134178E-138;
          scale = fmax(fabs(temp), fabs(bcmis));
          count++;
        }
        while ((scale <= 1.3435752215134178E-138) && (count <= 20)) {
          bcmis *= 7.4428285367870146E+137;
          temp *= 7.4428285367870146E+137;
          scale = fmax(fabs(temp), fabs(bcmis));
          count++;
        }
        tau = rt_hypotd_snf(bcmis, temp);
        *cs = sqrt(0.5 * (fabs(bcmis) / tau + 1.0));
        if (!(bcmis < 0.0)) {
          count = 1;
        } else {
          count = -1;
        }
        *sn = -(0.5 * temp / (tau * *cs)) * (double)count;
        bcmax = *a * *cs + *b * *sn;
        scale = -*a * *sn + *b * *cs;
        z = *c * *cs + *d * *sn;
        bcmis = -*c * *sn + *d * *cs;
        *b = scale * *cs + bcmis * *sn;
        *c = -bcmax * *sn + z * *cs;
        temp = 0.5 * ((bcmax * *cs + z * *sn) + (-scale * *sn + bcmis * *cs));
        *a = temp;
        *d = temp;
        if (*c != 0.0) {
          if (*b != 0.0) {
            if ((*b < 0.0) == (*c < 0.0)) {
              bcmis = sqrt(fabs(*b));
              scale = sqrt(fabs(*c));
              *a = bcmis * scale;
              if (!(*c < 0.0)) {
                p = *a;
              } else {
                p = -*a;
              }
              tau = 1.0 / sqrt(fabs(*b + *c));
              *a = temp + p;
              *d = temp - p;
              *b -= *c;
              *c = 0.0;
              bcmax = bcmis * tau;
              bcmis = scale * tau;
              temp = *cs * bcmax - *sn * bcmis;
              *sn = *cs * bcmis + *sn * bcmax;
              *cs = temp;
            }
          } else {
            *b = -*c;
            *c = 0.0;
            temp = *cs;
            *cs = -*sn;
            *sn = temp;
          }
        }
      }
    }
  }
  rt1r = *a;
  *rt2r = *d;
  if (*c == 0.0) {
    *rt1i = 0.0;
    *rt2i = 0.0;
  } else {
    *rt1i = sqrt(fabs(*b)) * sqrt(fabs(*c));
    *rt2i = -*rt1i;
  }
  return rt1r;
}

static int xhseqr(double h[1024], double z[1024])
{
  double bb;
  double d;
  double h22;
  double rt2r;
  double s;
  double s_tmp_tmp;
  double tr;
  double tst;
  int b_i;
  int b_k;
  int c_k;
  int i;
  int info;
  int k;
  int kdefl;
  int nr;
  boolean_T exitg1;
  info = 0;
  for (nr = 0; nr < 29; nr++) {
    i = nr + (nr << 5);
    h[i + 2] = 0.0;
    h[i + 3] = 0.0;
  }
  h[959] = 0.0;
  kdefl = 0;
  b_i = 31;
  exitg1 = false;
  while ((!exitg1) && (b_i + 1 >= 1)) {
    double h12;
    int i1;
    int its;
    int knt;
    int l;
    boolean_T converged;
    boolean_T exitg2;
    l = 1;
    converged = false;
    its = 0;
    exitg2 = false;
    while ((!exitg2) && (its < 961)) {
      boolean_T exitg3;
      k = b_i;
      exitg3 = false;
      while ((!exitg3) && (k + 1 > l)) {
        i = k + ((k - 1) << 5);
        d = fabs(h[i]);
        if (d <= 3.2066693760143564E-291) {
          exitg3 = true;
        } else {
          knt = k + (k << 5);
          h12 = h[knt];
          tr = fabs(h12);
          bb = h[i - 1];
          tst = fabs(bb) + tr;
          if (tst == 0.0) {
            if (k - 1 >= 1) {
              tst = fabs(h[(k + ((k - 2) << 5)) - 1]);
            }
            if (k + 2 <= 32) {
              tst += fabs(h[knt + 1]);
            }
          }
          if (d <= 2.2204460492503131E-16 * tst) {
            h22 = fabs(h[knt - 1]);
            h12 = fabs(bb - h12);
            tst = fmax(tr, h12);
            bb = fmin(tr, h12);
            s = tst + bb;
            if (fmin(d, h22) * (fmax(d, h22) / s) <=
                fmax(3.2066693760143564E-291,
                     2.2204460492503131E-16 * (bb * (tst / s)))) {
              exitg3 = true;
            } else {
              k--;
            }
          } else {
            k--;
          }
        }
      }
      l = k + 1;
      if (k + 1 > 1) {
        h[k + ((k - 1) << 5)] = 0.0;
      }
      if (k + 1 >= b_i) {
        converged = true;
        exitg2 = true;
      } else {
        double v[3];
        double rt1r;
        int m;
        kdefl++;
        if (kdefl - kdefl / 20 * 20 == 0) {
          s = fabs(h[b_i + ((b_i - 1) << 5)]) +
              fabs(h[(b_i + ((b_i - 2) << 5)) - 1]);
          tst = 0.75 * s + h[b_i + (b_i << 5)];
          h12 = -0.4375 * s;
          bb = s;
          h22 = tst;
        } else if (kdefl - kdefl / 10 * 10 == 0) {
          knt = k + (k << 5);
          s = fabs(h[knt + 1]) + fabs(h[(k + ((k + 1) << 5)) + 2]);
          tst = 0.75 * s + h[knt];
          h12 = -0.4375 * s;
          bb = s;
          h22 = tst;
        } else {
          knt = b_i + ((b_i - 1) << 5);
          tst = h[knt - 1];
          bb = h[knt];
          knt = b_i + (b_i << 5);
          h12 = h[knt - 1];
          h22 = h[knt];
        }
        s = ((fabs(tst) + fabs(h12)) + fabs(bb)) + fabs(h22);
        if (s == 0.0) {
          rt1r = 0.0;
          tr = 0.0;
          rt2r = 0.0;
          h22 = 0.0;
        } else {
          tst /= s;
          bb /= s;
          h12 /= s;
          h22 /= s;
          tr = (tst + h22) / 2.0;
          tst = (tst - tr) * (h22 - tr) - h12 * bb;
          bb = sqrt(fabs(tst));
          if (tst >= 0.0) {
            rt1r = tr * s;
            rt2r = rt1r;
            tr = bb * s;
            h22 = -tr;
          } else {
            rt1r = tr + bb;
            rt2r = tr - bb;
            if (fabs(rt1r - h22) <= fabs(rt2r - h22)) {
              rt1r *= s;
              rt2r = rt1r;
            } else {
              rt2r *= s;
              rt1r = rt2r;
            }
            tr = 0.0;
            h22 = 0.0;
          }
        }
        m = b_i - 1;
        exitg3 = false;
        while ((!exitg3) && (m >= k + 1)) {
          knt = m + ((m - 1) << 5);
          tst = h[knt];
          s_tmp_tmp = h[knt - 1];
          bb = s_tmp_tmp - rt2r;
          s = (fabs(bb) + fabs(h22)) + fabs(tst);
          h12 = tst / s;
          knt = m + (m << 5);
          v[0] = (h12 * h[knt - 1] + bb * (bb / s)) - tr * (h22 / s);
          tst = h[knt];
          v[1] = h12 * (((s_tmp_tmp + tst) - rt1r) - rt2r);
          v[2] = h12 * h[knt + 1];
          s = (fabs(v[0]) + fabs(v[1])) + fabs(v[2]);
          v[0] /= s;
          v[1] /= s;
          v[2] /= s;
          if (m == k + 1) {
            exitg3 = true;
          } else {
            i = m + ((m - 2) << 5);
            if (fabs(h[i - 1]) * (fabs(v[1]) + fabs(v[2])) <=
                2.2204460492503131E-16 * fabs(v[0]) *
                    ((fabs(h[i - 2]) + fabs(s_tmp_tmp)) + fabs(tst))) {
              exitg3 = true;
            } else {
              m--;
            }
          }
        }
        for (c_k = m; c_k <= b_i; c_k++) {
          knt = (b_i - c_k) + 2;
          if (knt >= 3) {
            nr = 3;
          } else {
            nr = knt;
          }
          if (c_k > m) {
            knt = (((c_k - 2) << 5) + c_k) - 1;
            for (b_k = 0; b_k < nr; b_k++) {
              v[b_k] = h[knt + b_k];
            }
          }
          h12 = v[0];
          tr = 0.0;
          if (nr > 0) {
            tst = b_xnrm2(nr - 1, v);
            if (tst != 0.0) {
              bb = rt_hypotd_snf(v[0], tst);
              if (v[0] >= 0.0) {
                bb = -bb;
              }
              if (fabs(bb) < 1.0020841800044864E-292) {
                knt = 0;
                do {
                  knt++;
                  for (b_k = 2; b_k <= nr; b_k++) {
                    v[b_k - 1] *= 9.9792015476736E+291;
                  }
                  bb *= 9.9792015476736E+291;
                  h12 *= 9.9792015476736E+291;
                } while ((fabs(bb) < 1.0020841800044864E-292) && (knt < 20));
                bb = rt_hypotd_snf(h12, b_xnrm2(nr - 1, v));
                if (h12 >= 0.0) {
                  bb = -bb;
                }
                tr = (bb - h12) / bb;
                tst = 1.0 / (h12 - bb);
                for (b_k = 2; b_k <= nr; b_k++) {
                  v[b_k - 1] *= tst;
                }
                for (b_k = 0; b_k < knt; b_k++) {
                  bb *= 1.0020841800044864E-292;
                }
                h12 = bb;
              } else {
                tr = (bb - v[0]) / bb;
                tst = 1.0 / (v[0] - bb);
                for (b_k = 2; b_k <= nr; b_k++) {
                  v[b_k - 1] *= tst;
                }
                h12 = bb;
              }
            }
          }
          if (c_k > m) {
            i = c_k + ((c_k - 2) << 5);
            h[i - 1] = h12;
            h[i] = 0.0;
            if (c_k < b_i) {
              h[i + 1] = 0.0;
            }
          } else if (m > k + 1) {
            i = (c_k + ((c_k - 2) << 5)) - 1;
            h[i] *= 1.0 - tr;
          }
          d = v[1];
          tst = tr * v[1];
          if (nr == 3) {
            s_tmp_tmp = v[2];
            h12 = tr * v[2];
            for (nr = c_k; nr < 33; nr++) {
              i = c_k + ((nr - 1) << 5);
              rt2r = h[i - 1];
              rt1r = h[i];
              s = h[i + 1];
              bb = (rt2r + d * rt1r) + s_tmp_tmp * s;
              rt2r -= bb * tr;
              h[i - 1] = rt2r;
              rt1r -= bb * tst;
              h[i] = rt1r;
              s -= bb * h12;
              h[i + 1] = s;
            }
            if (c_k + 3 <= b_i + 1) {
              i = c_k + 2;
            } else {
              i = b_i;
            }
            for (nr = 0; nr <= i; nr++) {
              b_k = nr + ((c_k - 1) << 5);
              rt2r = h[b_k];
              i1 = nr + (c_k << 5);
              rt1r = h[i1];
              knt = nr + ((c_k + 1) << 5);
              s = h[knt];
              bb = (rt2r + d * rt1r) + s_tmp_tmp * s;
              rt2r -= bb * tr;
              h[b_k] = rt2r;
              rt1r -= bb * tst;
              h[i1] = rt1r;
              s -= bb * h12;
              h[knt] = s;
            }
            for (nr = 0; nr < 32; nr++) {
              i = nr + ((c_k - 1) << 5);
              rt2r = z[i];
              b_k = nr + (c_k << 5);
              rt1r = z[b_k];
              i1 = nr + ((c_k + 1) << 5);
              s = z[i1];
              bb = (rt2r + d * rt1r) + s_tmp_tmp * s;
              rt2r -= bb * tr;
              z[i] = rt2r;
              rt1r -= bb * tst;
              z[b_k] = rt1r;
              s -= bb * h12;
              z[i1] = s;
            }
          } else if (nr == 2) {
            for (nr = c_k; nr < 33; nr++) {
              i = c_k + ((nr - 1) << 5);
              s_tmp_tmp = h[i - 1];
              rt2r = h[i];
              bb = s_tmp_tmp + d * rt2r;
              s_tmp_tmp -= bb * tr;
              h[i - 1] = s_tmp_tmp;
              rt2r -= bb * tst;
              h[i] = rt2r;
            }
            for (nr = 0; nr <= b_i; nr++) {
              i = nr + ((c_k - 1) << 5);
              s_tmp_tmp = h[i];
              b_k = nr + (c_k << 5);
              rt2r = h[b_k];
              bb = s_tmp_tmp + d * rt2r;
              s_tmp_tmp -= bb * tr;
              h[i] = s_tmp_tmp;
              rt2r -= bb * tst;
              h[b_k] = rt2r;
            }
            for (nr = 0; nr < 32; nr++) {
              i = nr + ((c_k - 1) << 5);
              s_tmp_tmp = z[i];
              b_k = nr + (c_k << 5);
              rt2r = z[b_k];
              bb = s_tmp_tmp + d * rt2r;
              s_tmp_tmp -= bb * tr;
              z[i] = s_tmp_tmp;
              rt2r -= bb * tst;
              z[b_k] = rt2r;
            }
          }
        }
        its++;
      }
    }
    if (!converged) {
      info = b_i + 1;
      exitg1 = true;
    } else {
      if ((l != b_i + 1) && (l == b_i)) {
        i = b_i << 5;
        b_k = b_i + i;
        d = h[b_k - 1];
        i1 = (b_i - 1) << 5;
        knt = b_i + i1;
        s_tmp_tmp = h[knt];
        rt2r = h[b_k];
        xdlanv2(&h[knt - 1], &d, &s_tmp_tmp, &rt2r, &s, &tst, &bb, &tr, &h22);
        h[b_k - 1] = d;
        h[knt] = s_tmp_tmp;
        h[b_k] = rt2r;
        if (b_i + 1 < 32) {
          knt = ((b_i + 1) << 5) + b_i;
          b_k = 30 - b_i;
          for (k = 0; k <= b_k; k++) {
            nr = knt + (k << 5);
            tst = h[nr];
            h12 = h[nr - 1];
            h[nr] = tr * tst - h22 * h12;
            h[nr - 1] = tr * h12 + h22 * tst;
          }
        }
        xrot(b_i - 1, h, i1 + 1, i + 1, tr, h22);
        xrot(32, z, i1 + 1, i + 1, tr, h22);
      }
      kdefl = 0;
      b_i = l - 2;
    }
  }
  for (nr = 0; nr < 30; nr++) {
    for (b_i = nr + 3; b_i < 33; b_i++) {
      h[(b_i + (nr << 5)) - 1] = 0.0;
    }
  }
  return info;
}

static double xnrm2(int n, const double x[1024], int ix0)
{
  double y;
  int k;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = fabs(x[ix0 - 1]);
    } else {
      double scale;
      int kend;
      scale = 3.3121686421112381E-170;
      kend = (ix0 + n) - 1;
      for (k = ix0; k <= kend; k++) {
        double absxk;
        absxk = fabs(x[k - 1]);
        if (absxk > scale) {
          double t;
          t = scale / absxk;
          y = y * t * t + 1.0;
          scale = absxk;
        } else {
          double t;
          t = absxk / scale;
          y += t * t;
        }
      }
      y = scale * sqrt(y);
    }
  }
  return y;
}

static void xrot(int n, double x[1024], int ix0, int iy0, double c, double s)
{
  int k;
  if (n >= 1) {
    int i;
    i = (unsigned char)n;
    for (k = 0; k < i; k++) {
      double b_temp_tmp;
      double temp_tmp;
      int b_temp_tmp_tmp;
      int temp_tmp_tmp;
      temp_tmp_tmp = (iy0 + k) - 1;
      temp_tmp = x[temp_tmp_tmp];
      b_temp_tmp_tmp = (ix0 + k) - 1;
      b_temp_tmp = x[b_temp_tmp_tmp];
      x[temp_tmp_tmp] = c * temp_tmp - s * b_temp_tmp;
      x[b_temp_tmp_tmp] = c * b_temp_tmp + s * temp_tmp;
    }
  }
}

static void xzlarf(int m, int n, int iv0, double tau, double C[1024], int ic0,
                   double work[32])
{
  int i;
  int ia;
  int lastc;
  int lastv;
  if (tau != 0.0) {
    boolean_T exitg2;
    lastv = m;
    i = iv0 + m;
    while ((lastv > 0) && (C[i - 2] == 0.0)) {
      lastv--;
      i--;
    }
    lastc = n;
    exitg2 = false;
    while ((!exitg2) && (lastc > 0)) {
      int exitg1;
      i = ic0 + ((lastc - 1) << 5);
      ia = i;
      do {
        exitg1 = 0;
        if (ia <= (i + lastv) - 1) {
          if (C[ia - 1] != 0.0) {
            exitg1 = 1;
          } else {
            ia++;
          }
        } else {
          lastc--;
          exitg1 = 2;
        }
      } while (exitg1 == 0);
      if (exitg1 == 1) {
        exitg2 = true;
      }
    }
  } else {
    lastv = 0;
    lastc = 0;
  }
  if (lastv > 0) {
    double c;
    int b_i;
    int i1;
    if (lastc != 0) {
      b_i = (unsigned char)lastc;
      memset(&work[0], 0, (unsigned int)b_i * sizeof(double));
      b_i = ic0 + ((lastc - 1) << 5);
      for (i = ic0; i <= b_i; i += 32) {
        c = 0.0;
        i1 = (i + lastv) - 1;
        for (ia = i; ia <= i1; ia++) {
          c += C[ia - 1] * C[((iv0 + ia) - i) - 1];
        }
        i1 = (i - ic0) >> 5;
        work[i1] += c;
      }
    }
    if (!(-tau == 0.0)) {
      i = ic0;
      b_i = (unsigned char)lastc;
      for (lastc = 0; lastc < b_i; lastc++) {
        c = work[lastc];
        if (c != 0.0) {
          c *= -tau;
          i1 = lastv + i;
          for (ia = i; ia < i1; ia++) {
            C[ia - 1] += C[((iv0 + ia) - i) - 1] * c;
          }
        }
        i += 32;
      }
    }
  }
}

void runica_simple(double *data, double *weights, double *sphere,
                   int nchan, int samples)
{
  const int ndata = nchan * samples;
  const int nmatrix = nchan * nchan;
  const int block = 52;
  const int nblocks = (samples + block - 1) / block;
  const int xb_size = nchan * (samples - nchan);

  /* Large arrays - dynamically allocated */
  double *b_y = (double *)malloc(ndata * sizeof(double));
  double *x = (double *)malloc(ndata * sizeof(double));
  double *Xb_data = (double *)malloc(xb_size * sizeof(double));
  double *tmp_data = (double *)malloc(xb_size * sizeof(double));
  double *timeperm = (double *)malloc(samples * sizeof(double));

  if (!b_y || !x || !Xb_data || !tmp_data || !timeperm) {
    free(b_y); free(x); free(Xb_data); free(tmp_data); free(timeperm);
    fprintf(stderr, "Error: Failed to allocate memory in runica_simple\n");
    return;
  }

  /* Stack arrays for nchan x nchan matrices */
  double BI_d[nmatrix];  /* Will be initialized as block * I */
  signed char BI[nmatrix];

  /* Initialize BI as block * I (identity scaled by block size) */
  memset(BI, 0, nmatrix * sizeof(signed char));
  for (int ii = 0; ii < nchan; ii++) {
    BI[ii + ii * nchan] = block;
  }

  /* Keep this initialization for compatibility - remove if code works without it */
  static const signed char BI_UNUSED[1024] = {
      52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 52};
  (void)BI_UNUSED;  /* Suppress unused warning */

  creal_T b_weights[nmatrix];
  creal_T c_weights[nmatrix];
  creal_T sphere_complex[nmatrix];
  creal_T winv[nmatrix];
  creal_T y[nchan];
  double dW[nmatrix];
  double d_weights[nmatrix];
  double olddelta[nmatrix];
  double oldweights[nmatrix];
  double rowmeans[nchan];
  double b_dW;
  double lrate;
  double muj;
  double oldchange;
  double re;
  double weights_re_tmp;
  int iidx[nchan];
  int ar;
  int br;
  int i;
  int idx;
  int j;
  int k;
  int step;
  int xpageoffset;
  signed char startweights[nmatrix];
  if (!isInitialized_runica_simple) {
    runica_simple_initialize();
  }
  /*  RUNICA_SIMPLE  Infomax ICA (logistic) with optional extended ICA (tanh)
   * and optional PCA. */
  /*  */
  /*  Usage */
  /*    [W,S] = runica_simple(X) */
  /*    [W,S,meanvar,bias,signs,lrates,acts] = runica_simple(X,'key',value,...)
   */
  /*  */
  /*  Input */
  /*    data: (chans x samples) */
  /*  */
  /*  Keys kept (others removed: specgram, posact, interrupt) */
  /*    'extended'   : 0 off, >0 update signs every N blocks, <0 fix #subgauss =
   * -N */
  /*    'pca'        : 0 off, N keep N PCs */
  /*    'sphering'   : 'on' or 'off' */
  /*    'weights'    : initial W (ncomps x chans) or 0 to use identity */
  /*    'lrate'      : learning rate */
  /*    'block'      : block size */
  /*    'anneal'     : anneal factor in (0,1] */
  /*    'annealdeg'  : anneal threshold in degrees */
  /*    'stop'       : stopping threshold on weight change */
  /*    'maxsteps'   : max passes */
  /*    'bias'       : 'on' or 'off' */
  /*    'momentum'   : in [0,1) */
  /*    'verbose'    : 'on' or 'off' */
  /*    'logfile'    : filename or [] */
  /*    'rndreset'   : 'yes' or 'no' or logical */
  /*  */
  /*  Notes */
  /*    Unmixing is W*S (if sphering on). If pca on, returned W already includes
   * PCA, and S = eye(original_chans). */
  /*  Minimal defaults */
  /*  Extended ICA internals */
  /*  Validate */
  /*  Logging */
  /*  Random seed */
  /* if resetSeed, rand('state',sum(100*clock)); else, rand('state',0); end
   * %#ok<RAND> */
  /*  Remove channel means */
  mean(data, rowmeans, nchan, samples);
  for (i = 0; i < samples; i++) {
    for (k = 0; k < nchan; k++) {
      idx = k + i * nchan;
      data[idx] -= rowmeans[k];
    }
  }
  /*  Optional PCA */
  /*  Sphering */
  /*  fprintf('Computing sphering matrix\n'); */
  for (j = 0; j < nchan; j++) {
    muj = 0.0;
    for (xpageoffset = 0; xpageoffset < samples; xpageoffset++) {
      weights_re_tmp = data[j + (xpageoffset  * nchan)];
      x[xpageoffset + samples * j] = weights_re_tmp;
      muj += weights_re_tmp;
    }
    muj /= (double)samples;
    for (xpageoffset = 0; xpageoffset < samples; xpageoffset++) {
      idx = xpageoffset + samples * j;
      x[idx] -= muj;
    }
  }
  for (idx = 0; idx <= (nchan-1)*nchan; idx += nchan) {
    i = idx + 1;
    k = idx + nchan;
    if (i <= k) {
      memset(&dW[i + -1], 0, (unsigned int)((k - i) + 1) * sizeof(double));
    }
  }
  br = -1;
  for (idx = 0; idx <= (nchan-1)*nchan; idx += nchan) {
    ar = -1;
    i = idx + 1;
    k = idx + nchan;
    for (j = i; j <= k; j++) {
      muj = 0.0;
      for (xpageoffset = 0; xpageoffset < samples; xpageoffset++) {
        muj += x[(xpageoffset + ar) + 1] * x[(xpageoffset + br) + 1];
      }
      dW[j - 1] += 3.2783660623545228E-5 * muj;
      ar += samples;
    }
    br += samples;
  }
  sqrtm(dW, b_weights, nchan);
  inv(b_weights, sphere_complex, nchan);
  for (i = 0; i < nmatrix; i++) {
    sphere_complex[i].re *= 2.0;
    sphere_complex[i].im *= 2.0;
  }
  for (i = 0; i < nchan; i++) {
    for (k = 0; k < samples; k++) {
      weights_re_tmp = 0.0;
      for (ar = 0; ar < nchan; ar++) {
        weights_re_tmp += sphere_complex[i + (ar  * nchan)].re * data[ar + (k  * nchan)];
      }
      x[i + (k  * nchan)] = weights_re_tmp;
    }
  }
  memcpy(&data[0], &x[0], ndata * sizeof(double));
  /*  Init weights */
  memset(&weights[0], 0, nmatrix * sizeof(double));
  /*  Extended ICA configuration */
  /*  Training constants */
  lrate = 0.00018755035531556522;
  for (k = 0; k < nchan; k++) {
    weights[k + (k  * nchan)] = 1.0;
    rowmeans[k] = 0.0;
  }
  /*  Signs init for extended ICA */
  for (i = 0; i < nmatrix; i++) {
    k = (int)weights[i];
    startweights[i] = (signed char)k;
    oldweights[i] = k;
    olddelta[i] = 0.0;
  }
  oldchange = rtInf;
  /*  fprintf('Beginning ICA training: ncomps=%f, extended=%f, bias=%f\n',
   * ncomps, extended, useBias); */
  step = 0;
  int exitg1;
  do {
    exitg1 = 0;
    if (step < 512) {
      int t;
      boolean_T exitg2;
      boolean_T wts_blowup;
      randperm(timeperm, samples);
      wts_blowup = false;
      t = 0;
      exitg2 = false;
      while ((!exitg2) && (t < 586)) {
        double U[1664];
        double b_x[1664];
        idx = t * 52;
        for (j = 0; j < 52; j++) {
          br = j  * nchan;
          for (xpageoffset = 0; xpageoffset < nchan; xpageoffset++) {
            Xb_data[xpageoffset + nchan * j] =
                data[xpageoffset + (((int)timeperm[idx + j] - 1)  * nchan)];
            tmp_data[br + xpageoffset] = 0.0;
          }
          for (k = 0; k < nchan; k++) {
            ar = k  * nchan;
            muj = Xb_data[br + k];
            for (xpageoffset = 0; xpageoffset < nchan; xpageoffset++) {
              i = br + xpageoffset;
              tmp_data[i] += weights[ar + xpageoffset] * muj;
            }
          }
        }
        for (i = 0; i < 52; i++) {
          for (k = 0; k < nchan; k++) {
            U[k + (i  * nchan)] = tmp_data[k + nchan * i] + rowmeans[k];
          }
        }
        for (k = 0; k < 1664; k++) {
          b_x[k] = 1.0 - 2.0 * (1.0 / (exp(-U[k]) + 1.0));
        }
        for (i = 0; i < nchan; i++) {
          for (k = 0; k < nchan; k++) {
            weights_re_tmp = 0.0;
            for (ar = 0; ar < 52; ar++) {
              idx = ar  * nchan;
              weights_re_tmp += b_x[i + idx] * U[k + idx];
            }
            idx = i + (k  * nchan);
            dW[idx] = weights_re_tmp + (double)BI[idx];
          }
        }
        for (i = 0; i < nmatrix; i++) {
          dW[i] *= lrate;
        }
        for (i = 0; i < nchan; i++) {
          for (k = 0; k < nchan; k++) {
            weights_re_tmp = 0.0;
            for (ar = 0; ar < nchan; ar++) {
              weights_re_tmp += dW[i + (ar  * nchan)] * weights[ar + (k  * nchan)];
            }
            idx = i + (k  * nchan);
            d_weights[idx] = weights[idx] + weights_re_tmp;
          }
        }
        memcpy(&weights[0], &d_weights[0], nmatrix * sizeof(double));
        for (i = 0; i < nchan; i++) {
          for (k = 0; k < 52; k++) {
            U[k + 52 * i] = b_x[i + (k  * nchan)];
          }
        }
        for (j = 0; j < nchan; j++) {
          xpageoffset = j * 52;
          weights_re_tmp = U[xpageoffset];
          for (k = 0; k < 51; k++) {
            weights_re_tmp += U[(xpageoffset + k) + 1];
          }
          rowmeans[j] += lrate * weights_re_tmp;
        }
        for (k = 0; k < nmatrix; k++) {
          dW[k] = fabs(weights[k]);
        }
        if (!rtIsNaN(dW[0])) {
          idx = 1;
        } else {
          boolean_T exitg3;
          idx = 0;
          k = 2;
          exitg3 = false;
          while ((!exitg3) && (k < 1025)) {
            if (!rtIsNaN(dW[k - 1])) {
              idx = k;
              exitg3 = true;
            } else {
              k++;
            }
          }
        }
        if (idx == 0) {
          muj = dW[0];
        } else {
          muj = dW[idx - 1];
          i = idx + 1;
          for (k = i; k < 1025; k++) {
            weights_re_tmp = dW[k - 1];
            if (muj < weights_re_tmp) {
              muj = weights_re_tmp;
            }
          }
        }
        if (muj > 1.0E+8) {
          wts_blowup = true;
          exitg2 = true;
        } else {
          /*  Extended sign updates by kurtosis */
          t++;
        }
      }
      if (wts_blowup) {
        lrate *= 0.9;
        /*  fprintf('Weights blew up, restarting with lrate=%g\n', lrate); */
        for (i = 0; i < nmatrix; i++) {
          idx = startweights[i];
          weights[i] = idx;
          oldweights[i] = idx;
        }
        memset(&rowmeans[0], 0, nchan * sizeof(double));
        step = 0;
      } else {
        /*  Step metrics */
        step++;
        b_dW = 0.0;
        for (i = 0; i < nmatrix; i++) {
          weights_re_tmp = weights[i] - oldweights[i];
          dW[i] = weights_re_tmp;
          b_dW += weights_re_tmp * weights_re_tmp;
        }
        if ((step > 2) && ((!rtIsInf(oldchange)) && (!rtIsNaN(oldchange))) &&
            (oldchange > 0.0)) {
          muj = 0.0;
          for (i = 0; i < nmatrix; i++) {
            weights_re_tmp = weights[i] - oldweights[i];
            oldweights[i] = weights_re_tmp;
            muj += weights_re_tmp * olddelta[i];
          }
          muj = acos(muj / sqrt(b_dW * oldchange));
        } else {
          muj = 0.0;
        }
        muj *= 57.295779513082323;
        printf("step %g  lrate %g  wchange %.9g  angledelta %.3g deg\n",
               (double)step, lrate, b_dW, muj);
        fflush(stdout);
        /*  Anneal */
        if (muj > 60.0) {
          lrate *= 0.9;
          memcpy(&olddelta[0], &dW[0], nmatrix * sizeof(double));
          oldchange = b_dW;
        } else if (step == 1) {
          memcpy(&olddelta[0], &dW[0], nmatrix * sizeof(double));
          oldchange = b_dW;
        }
        /*  Stop and blowup guard */
        if ((step > 2) && (b_dW < 1.0E-6)) {
          exitg1 = 1;
        } else {
          if (b_dW > 1.0E+9) {
            lrate *= 0.8;
          }
          memcpy(&oldweights[0], &weights[0], nmatrix * sizeof(double));
        }
      }
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);
  /*  Activations for output */
  /*  If PCA was used, compose matrices back to original channel space */
  /*  Component variance ranking */
  for (i = 0; i < nmatrix; i++) {
    b_weights[i].re = weights[i];
    b_weights[i].im = 0.0;
  }
  for (i = 0; i < nchan; i++) {
    for (k = 0; k < nchan; k++) {
      re = 0.0;
      muj = 0.0;
      for (ar = 0; ar < nchan; ar++) {
        idx = i + (ar  * nchan);
        b_dW = b_weights[idx].re;
        br = ar + (k  * nchan);
        weights_re_tmp = sphere_complex[br].im;
        lrate = b_weights[idx].im;
        oldchange = sphere_complex[br].re;
        re += b_dW * oldchange - lrate * weights_re_tmp;
        muj += b_dW * weights_re_tmp + lrate * oldchange;
      }
      ar = i + (k  * nchan);
      c_weights[ar].re = re;
      c_weights[ar].im = muj;
    }
  }
  inv(c_weights, winv, nchan);
  for (k = 0; k < nmatrix; k++) {
    weights_re_tmp = winv[k].re;
    muj = winv[k].im;
    b_weights[k].re = weights_re_tmp * weights_re_tmp - muj * muj;
    weights_re_tmp *= muj;
    b_weights[k].im = weights_re_tmp + weights_re_tmp;
  }
  for (j = 0; j < nchan; j++) {
    xpageoffset = j  * nchan;
    re = b_weights[xpageoffset].re;
    muj = b_weights[xpageoffset].im;
    for (k = 0; k < 31; k++) {
      idx = (xpageoffset + k) + 1;
      re += b_weights[idx].re;
      muj += b_weights[idx].im;
    }
    y[j].re = re;
    y[j].im = muj;
  }
  for (i = 0; i < samples; i++) {
    for (k = 0; k < nchan; k++) {
      weights_re_tmp = 0.0;
      for (ar = 0; ar < nchan; ar++) {
        weights_re_tmp += weights[k + (ar  * nchan)] * data[ar + (i  * nchan)];
      }
      x[i + samples * k] = weights_re_tmp;
    }
  }
  for (k = 0; k < ndata; k++) {
    weights_re_tmp = x[k];
    b_y[k] = weights_re_tmp * weights_re_tmp;
  }
  /*  Sort descending */
  for (j = 0; j < nchan; j++) {
    xpageoffset = j * samples;
    weights_re_tmp = b_y[xpageoffset];
    for (k = 0; k < samples - 1; k++) {
      weights_re_tmp += b_y[(xpageoffset + k) + 1];
    }
    rowmeans[j] = weights_re_tmp;
    for (idx = 0; idx < 29; idx++) {
      br = xpageoffset + ((idx + 1)  * nmatrix);
      muj = b_y[br];
      if (idx + 2 == 30) {
        ar = 808;
      } else {
        ar = nmatrix;
      }
      for (k = 2; k <= ar; k++) {
        muj += b_y[(br + k) - 1];
      }
      rowmeans[j] += muj;
    }
    weights_re_tmp = rowmeans[j];
    muj = y[j].re * weights_re_tmp;
    b_dW = y[j].im * weights_re_tmp;
    if (b_dW == 0.0) {
      re = muj / ndata - 1.0;
      muj = 0.0;
    } else if (muj == 0.0) {
      re = 0.0;
      muj = b_dW / ndata - 1.0;
    } else {
      re = muj / ndata - 1.0;
      muj = b_dW / ndata - 1.0;
    }
    y[j].re = re;
    y[j].im = muj;
  }
  sort(y, iidx);
  for (i = 0; i < nchan; i++) {
    idx = i  * nchan;
    for (k = 0; k < nchan; k++) {
      d_weights[k + idx] = weights[(iidx[k] + idx) - 1];
    }
  }
  memcpy(&weights[0], &d_weights[0], nmatrix * sizeof(double));
  /* Copy real part of sphere_complex to sphere output */
  for (i = 0; i < nmatrix; i++) {
    sphere[i] = sphere_complex[i].re;
  }

  /* Free dynamically allocated arrays */
  free(b_y);
  free(x);
  free(Xb_data);
  free(tmp_data);
  free(timeperm);
}

void runica_simple_initialize(void)
{
  c_eml_rand_mt19937ar_stateful_i();
  isInitialized_runica_simple = true;
}

void runica_simple_terminate(void)
{
  isInitialized_runica_simple = false;
}

/* End of code generation (runica_simple.c) */
