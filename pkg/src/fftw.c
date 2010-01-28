/*
 * fftw.c - wrapper around fftw3
 *
 * Wraps the plan creation and execution phases of fftw. Currently
 * both 1D c2c FFTs and r2r DCTs are supported.
 * 
 * Author:
 *  Sebastian Krey <skrey@statistik.tu-dortmund.de>
 *  Olaf Mersmann  <olafm@statistik.tu-dortmund.de>
 */

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Applic.h>
#include <fftw3.h>

#include <assert.h>

#define ALLOC_COMPLEX_VECTOR(S, D, N)		\
  SEXP S;					\
  PROTECT(S = allocVector(CPLXSXP, N));		\
  Rcomplex *D = COMPLEX(S);

#define ALLOC_REAL_VECTOR(S, D, N)		\
  SEXP S;					\
  PROTECT(S = allocVector(REALSXP, N));		\
  double *D = REAL(S);

/* Is fftw initialized? */
static int initialized = FALSE;

/* Holds a plan for the forward and reverse FFT of size 'size'. */
typedef struct {
  R_len_t size;
  fftw_complex *in;
  fftw_complex *out;
  fftw_plan forward;
  fftw_plan backward;
} fft_plan;

/* fft_plan destructor */
void fft_plan_finalizer(SEXP s_plan) {
  fft_plan *plan = R_ExternalPtrAddr(s_plan);

  if (NULL != plan->in)
    fftw_free(plan->in);
  if (NULL != plan->out)
    fftw_free(plan->out);
  if (NULL != plan->forward)
    fftw_destroy_plan(plan->forward);
  if (NULL != plan->backward)
    fftw_destroy_plan(plan->backward);
  Free(plan);
}

/* Holds a plan for the forward and reverse type 'type DCT of size 'size'. */
typedef struct {
  R_len_t size;
  int type;
  double *in;
  double *out;
  fftw_plan forward;
  fftw_plan backward;
} dct_plan;

void dct_plan_finalizer(SEXP s_plan) {
  dct_plan *plan = R_ExternalPtrAddr(s_plan);

  if (NULL != plan->in)
    fftw_free(plan->in);
  if (NULL != plan->out)
    fftw_free(plan->out);
  if (NULL != plan->forward)
    fftw_destroy_plan(plan->forward);
  /* Check that forward and backward are not the same plan (type I DCT): */
  if (NULL != plan->backward && plan->backward != plan->forward)
    fftw_destroy_plan(plan->backward);
  Free(plan);
}

/* Helper routines: */
static int choose_effort(SEXP s_effort) {
  int effort = INTEGER(s_effort)[0];
  
  if (effort <= 0) {
    return FFTW_ESTIMATE;
  } else if (effort == 1) {
    return FFTW_MEASURE;
  } else if (effort == 2) {
    return FFTW_PATIENT;
  } else {
    return FFTW_EXHAUSTIVE;
  }
}

/* General FFT case: */
SEXP FFT_print_plan(SEXP s_plan) {
  fft_plan *plan = R_ExternalPtrAddr(s_plan);
  assert(plan != NULL);

  Rprintf("plan->size     : %i\n", plan->size);
  Rprintf("plan->in       : 0x%08x\n", plan->in);
  Rprintf("plan->out      : 0x%08x\n", plan->out);
  Rprintf("plan->forward  : 0x%08x\n", plan->forward);
  Rprintf("plan->backward : 0x%08x\n", plan->backward);
  return R_NilValue;
}

SEXP FFT_plan(SEXP s_n, SEXP s_effort) {
  fft_plan *plan;
  int effort = choose_effort(s_effort);

  /* If s_n is a single integer, assume it is the length */
  R_len_t n = length(s_n);
  if (n == 1)
    n = INTEGER(s_n)[0];

  /* Possibly initalize fftw: */
  if (!initialized) {
    fftw_import_system_wisdom();
    initialized = TRUE;
  }

  plan           = Calloc(1, fft_plan);
  plan->size     = n;
  plan->in       = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
  plan->out      = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
  plan->forward  = fftw_plan_dft_1d(plan->size, plan->in, plan->out, 
				    FFTW_FORWARD, 
				    FFTW_DESTROY_INPUT | effort);
  plan->backward = fftw_plan_dft_1d(plan->size, plan->in, plan->out, 
				    FFTW_BACKWARD, 
				    FFTW_DESTROY_INPUT | effort);

  SEXP s_ret = R_MakeExternalPtr((void *)plan, R_NilValue, R_NilValue);
  assert(s_ret != R_NilValue);
  R_RegisterCFinalizer(s_ret, fft_plan_finalizer);
  return s_ret;
}

SEXP FFT_execute(SEXP s_plan, SEXP s_x, SEXP s_inv) {
  R_len_t i, n;
  fft_plan *plan;
  fftw_plan p;

  plan = R_ExternalPtrAddr(s_plan);
  assert(plan != NULL);

  /* Extract fftw plan: */
  if (INTEGER(s_inv)[0] == FALSE) {    
    p = plan->forward;
  } else {
    p = plan->backward;
  }

  /* Extract input vector: */
  n = length(s_x);
  if (n != plan->size) {
    error("Input and plan size differ.");
    return R_NilValue;
  }
  if (CPLXSXP == TYPEOF(s_x)) {
    Rcomplex *x = COMPLEX(s_x);    
    for (i = 0; i < n; ++i) {
      plan->in[i][0] = x[i].r;
      plan->in[i][1] = x[i].i;
    }  
  } else if (REALSXP == TYPEOF(s_x)) {
    double *x = REAL(s_x);
    for (i = 0; i < n; ++i) {
      plan->in[i][0] = x[i];
      plan->in[i][1] = 0.0;
    }  
  } else {
    error("'s_x' must be real or complex.");
  }

  /* Execute plan: */
  fftw_execute(p);

  /* Copy output: */
  ALLOC_COMPLEX_VECTOR(s_ret, ret, n);
  for (i = 0; i < n; ++i) {
    ret[i].r = plan->out[i][0];
    ret[i].i = plan->out[i][1];
  }

  UNPROTECT(1); /* s_ret */
  return s_ret;
}

/* real to real DCT case: */
SEXP DCT_plan(SEXP s_n, SEXP s_type, SEXP s_effort) {
  dct_plan *plan;
  R_len_t n = length(s_n);
  int type = INTEGER(s_type)[0];  
  fftw_r2r_kind fw_type, bw_type;

  int effort = choose_effort(s_effort);

  /* Decide on type of DCT */
  if (type == 1) {
    fw_type = FFTW_REDFT00;
    bw_type = fw_type;
  } else if (type == 2) {
    fw_type = FFTW_REDFT10;
    bw_type = FFTW_REDFT01;
  } else if (type == 3) {
    fw_type = FFTW_REDFT01;
    bw_type = FFTW_REDFT10;
  } else if (type == 4) {
    fw_type = FFTW_REDFT11;
    bw_type = fw_type;
  } else {
    error("Unknown type specified.");
    return NULL;
  }
    
  /* If s_n is a single integer, assume it is the length */
  if (n == 1)
    n = INTEGER(s_n)[0];
  
  if (!initialized) {
    fftw_import_system_wisdom();
    initialized = TRUE;
  }
  
  plan           = Calloc(1, dct_plan);
  plan->size     = n;
  plan->in       = (double*) fftw_malloc(sizeof(double) * n);
  plan->out      = (double*) fftw_malloc(sizeof(double) * n);
  plan->forward  = fftw_plan_r2r_1d(plan->size, plan->in, plan->out,
				    fw_type, 
				    FFTW_DESTROY_INPUT | effort);
  if (bw_type == fw_type) {
    plan->backward = plan->forward;
  } else {
    plan->backward = fftw_plan_r2r_1d(plan->size, plan->in, plan->out,
				      bw_type, 
				      FFTW_DESTROY_INPUT | effort);
  }

  SEXP s_ret = R_MakeExternalPtr((void *)plan, R_NilValue, R_NilValue);
  assert(s_ret != R_NilValue);
  R_RegisterCFinalizer(s_ret, fft_plan_finalizer);
  return s_ret;
}

SEXP DCT_execute(SEXP s_plan, SEXP s_x, SEXP s_inv) {
  R_len_t i, n;
  dct_plan *plan;
  fftw_plan p;

  plan = R_ExternalPtrAddr(s_plan);
  assert(plan != NULL);

  /* Extract fftw plan: */
  if (INTEGER(s_inv)[0] == FALSE) {    
    p = plan->forward;
  } else {
    p = plan->backward;
  }

  /* Extract input vector: */
  n = length(s_x);
  if (n != plan->size) {
    error("Input and plan size differ.");
    return R_NilValue;
  }
  if (REALSXP == TYPEOF(s_x)) {
    double *x = REAL(s_x);
    for (i = 0; i < n; ++i) {
      plan->in[i] = x[i];
    }  
  } else {
    error("'s_x' must be real.");
  }

  /* Execute plan: */
  fftw_execute(p);

  /* Copy output: */
  ALLOC_REAL_VECTOR(s_ret, ret, n);
  for (i = 0; i < n; ++i) {
    ret[i] = plan->out[i];
  }

  UNPROTECT(1); /* s_ret */
  return s_ret;
}
