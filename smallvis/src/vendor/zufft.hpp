/* zufft's C interface for other packages (design §7).
 *
 * Use it with `LinkingTo: zufft` and `Imports: zufft` (import something
 * from zufft in NAMESPACE, so that zufft is loaded first). Resolve the table
 * once, on the R thread, during your package's initialisation:
 *
 *     static const zufft_api *fft;
 *     void R_init_mypkg(DllInfo *dll) {
 *       fft = zufft_api_get();
 *       if (fft == NULL) ...;  // zufft too old, or not loaded
 *     }
 *
 * The table's members never call R and may be called from any thread;
 * zufft serialises top-level transforms internally. After a failure,
 * fft->message() returns this thread's diagnostic until the next call.
 *
 * Complex buffers hold pairs of doubles (re, im), such as R's Rcomplex,
 * C99 double _Complex or C++ std::complex<double>. Axes are 0-based.
 * Strides count elements: one double for real buffers, one (re, im) pair
 * for complex ones. A NULL stride means dense column-major (R's layout).
 *
 * The table is append-only: new members go at the end and struct_size tells
 * a consumer which exist. A layout change would register "zufft_api_v2".
 * Define ZUFFT_NO_R before including this header to use the types and
 * zufft_scale() without R's headers (zufft's own fuzzing does); it leaves
 * out zufft_api_get().
 *
 * This header is MIT-licensed, like zufft's own code. */
#ifndef ZUFFT_H
#define ZUFFT_H

#include <stddef.h>
#include <math.h>
#ifndef ZUFFT_NO_R
#include <R_ext/Rdynload.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define ZUFFT_ABI_VERSION 1

typedef struct {
  int ndim;                /* number of dimensions, >= 1 */
  const size_t *shape;     /* extents, each >= 1 */
  const ptrdiff_t *stride; /* elements per step, >= 0; NULL = column-major */
} zufft_layout;

typedef enum {
  ZUFFT_OK = 0,
  ZUFFT_EINVAL = 1,    /* bad argument: NULL, axes, type, overlap */
  ZUFFT_ESHAPE = 2,    /* inconsistent or empty shapes */
  ZUFFT_EOVERFLOW = 3, /* sizes or addresses overflow */
  ZUFFT_ENOMEM = 4,    /* out of memory */
  ZUFFT_ETHREAD = 5,   /* a worker thread could not start */
  ZUFFT_EENGINE = 6,   /* the transform engine failed; see message() */
  ZUFFT_EINTERNAL = 7  /* a bug in zufft */
} zufft_status;

typedef struct {
  size_t struct_size;       /* sizeof(zufft_api) in the provider */
  unsigned int abi_version; /* ZUFFT_ABI_VERSION */
  unsigned int flags;       /* reserved, zero */

  /* Complex-to-complex over `axes`; in == out with identical layouts runs
   * in place. forward != 0 uses exp(-2 pi i jk / n). The result is scaled
   * by fct (1 for FFTW's unnormalised convention). nthreads: a count, -k
   * for all cores but k - 1, or 0 for zufft's current fft_workers(). */
  int (*c2c)(const zufft_layout *lin, const double *in,
             const zufft_layout *lout, double *out,
             int naxes, const int *axes, int forward,
             double fct, int nthreads);
  /* Real-to-complex: the output's extent on the last listed axis is
   * n / 2 + 1. in and out must not overlap. */
  int (*r2c)(const zufft_layout *lin, const double *in,
             const zufft_layout *lout, double *out,
             int naxes, const int *axes, int forward,
             double fct, int nthreads);
  /* Complex-to-real: the input's extent on the last listed axis is
   * n / 2 + 1 for an output extent n. The input is not modified (FFTW's
   * c2r destroys it). in and out must not overlap. */
  int (*c2r)(const zufft_layout *lin, const double *in,
             const zufft_layout *lout, double *out,
             int naxes, const int *axes, int forward,
             double fct, int nthreads);
  /* DCT and DST of type 1-4 (SciPy's definitions, unnormalised unless fct
   * says otherwise); ortho != 0 applies SciPy's orthogonalize correction.
   * For an inverse, swap types 2 and 3 yourself, as SciPy does. */
  int (*dct)(const zufft_layout *lin, const double *in,
             const zufft_layout *lout, double *out,
             int naxes, const int *axes, int type,
             double fct, int ortho, int nthreads);
  int (*dst)(const zufft_layout *lin, const double *in,
             const zufft_layout *lout, double *out,
             int naxes, const int *axes, int type,
             double fct, int ortho, int nthreads);
  /* Smallest fast length >= n (real != 0: factors 2, 3, 5; otherwise 2, 3,
   * 5, 7, 11). Returns 0 for n == 0 or n > 2^52. */
  size_t (*good_size)(size_t n, int real);
  /* This thread's last failure message, "" after a success. Copy it before
   * the next call on this thread. */
  const char *(*message)(void);
  /* ... appended over time ... */
} zufft_api;

/* SciPy's normalisations. */
typedef enum {
  ZUFFT_NORM_NONE = 0,  /* no scaling */
  ZUFFT_NORM_ORTHO = 1, /* 1 / sqrt(N) */
  ZUFFT_NORM_BY_N = 2   /* 1 / N */
} zufft_norm;

/* The scale factor for a transform over `naxes` axes of logical lengths
 * len[], computed as SciPy does: N is the product over axes of
 * mult * (len + delta), in long double, rounded once to double. Use
 * mult = 1, delta = 0 for Fourier transforms (with the output lengths for
 * c2r), mult = 2 and delta = -1 for DCT-I, mult = 2 and delta = 1 for
 * DST-I, mult = 2 and delta = 0 for the other DCT and DST types. Returns
 * ZUFFT_OK and sets *fct, or ZUFFT_EINVAL for a bad argument. */
static inline int zufft_scale(zufft_norm norm, int naxes, const size_t *len,
                              size_t mult, int delta, double *fct) {
  long double n = 1.0L;
  int i;
  if (fct == NULL || naxes < 1 || len == NULL || mult < 1) return ZUFFT_EINVAL;
  for (i = 0; i < naxes; ++i) {
    long double f = (long double) mult * ((long double) len[i] + (long double) delta);
    if (!(f >= 1.0L)) return ZUFFT_EINVAL;
    n *= f;
  }
  switch (norm) {
    case ZUFFT_NORM_NONE: *fct = 1.0; return ZUFFT_OK;
    case ZUFFT_NORM_BY_N: *fct = (double) (1.0L / n); return ZUFFT_OK;
    case ZUFFT_NORM_ORTHO: *fct = (double) (1.0L / sqrtl(n)); return ZUFFT_OK;
  }
  return ZUFFT_EINVAL;
}

/* Resolves zufft's table through R_GetCCallable. Call on the R thread, with
 * zufft loaded (R_GetCCallable signals an R error otherwise). Returns NULL
 * if the provider's table is older than this header or of another ABI. */
#ifndef ZUFFT_NO_R
static inline const zufft_api *zufft_api_get(void) {
  typedef const zufft_api *(*zufft_getter)(void);
  zufft_getter get = (zufft_getter) R_GetCCallable("zufft", "zufft_api_v1");
  const zufft_api *api = get ? get() : NULL;
  if (api == NULL || api->abi_version != ZUFFT_ABI_VERSION ||
      api->struct_size < sizeof(zufft_api))
    return NULL;
  return api;
}
#endif

#ifdef __cplusplus
}
#endif

#endif
