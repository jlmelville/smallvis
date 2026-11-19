# Smallvis correctness and maintenance

The user asked to record the findings from the ChatGPT Pro review and the independent
follow-up, then start working through them in order. The review baseline is
`426fe13d574dd601da316079f3a6a3d2da34b3ef`. Targeted probes on an R 4.6.2 installation
reproduced all eight numbered findings and the three smaller numerical failures.
This plan lives at the repository root, outside the nested R package.
The supplied review also reported clean Ubuntu and Windows package checks for
source commit `0e06962` and sanitizer runs for the native tree and threading
harnesses. The local reproduction and the checks below are a separate record.

## Current state

All ten items are complete. The integrated package checks passed. Deferred
proposals remain below with their revisit conditions.

## Plan of work and acceptance

1. [x] **Stabilize Gaussian perplexity calibration.** For
   `X <- diag(sqrt(1 + 1e-4 * (0:23)))`, dense R and C++ each fell back on all
   24 rows for perplexity 4. The kNN example with squared distances
   `2 + 1e-4 * (0:11)` produced perplexity 12, versus 4.00003 after a common
   distance shift. Solve for entropy using distances relative to the nearest
   non-self neighbor. Check R, native dense, and saved-neighbor Gaussian
   paths against the independent conditional probabilities. Preserve raw `V`
   and degree diagnostics, raw weighting when row normalization is disabled
   or follows symmetrization, and route-specific fallbacks for true ties.
   Verify ordinary raw weights against `exp(-beta * D^2)` and avoid extra
   dense matrices on the default route. Start with `R/perplexity.R`,
   `src/perplexity.cpp`, and `tests/testthat/test_native_perplexity.R` under
   `smallvis/`.
2. [x] **Implement the ABSSNE limits and close the ABSNE threshold hole.**
   Equal `P` and `Q` gave ABSSNE costs of 5.025271 at the default parameters
   and about 5.6e11 at `alpha = lambda = 0`. The reverse-KL limit had the wrong
   directional derivative. Use explicit limiting formulas and compare scores
   and derivatives to independent references, including the default,
   `alpha = 0`, `beta = 0`, `lambda = 0`, the joint zero limit, and
   `alpha = lambda = 1e-5` for ABSNE. Factor mathematical logic only where
   both kernels can reuse it directly.
3. [x] **Guard function-based optimization for unnormalized t-SNE.** For a
   non-unit-mass `P`, the reported KL differentiates to a force using
   `P - sum(P) * Q`, while the experimental force uses `P - Q`. L-BFGS
   accepted the mismatch. Preserve the force, use the existing guard, and
   document the reported KL's diagnostic role. Check rejection and continued
   support for normalized L-BFGS and unnormalized gradient-only fits.
4. [x] **Make the FFT grid translation invariant.** On the 24-point unit-circle
   fixture, translating only x by 1,000 increased the relative maximum force
   error from about 0.00028% to 6.39%. Set spacing from axis extents, use
   per-axis origins, and compare a translated grid and one-step fit against
   direct forces. Preserve the rectangular padding and shared kernel spacing.
5. [x] **Honor dense HLSNE and CETSNE symmetrization.** Both accepted
   `"mutual"` and returned the same `P` as `"symmetric"`. Forward the selected
   option and compare affinities on a fixture where the options differ.
6. [x] **Repair the maintained examples.** Replace `lv_eps` with `gr_eps`;
   both UMAP/L-BFGS examples ran with `gr_eps = 0`. Explain that feature PCA
   needs features and that PCA initialization accepts `dist` via classical
   MDS. Execute the corrected examples and regenerate affected help.
7. [x] **Keep normalized spectral initialization matrix-shaped.** Both
   eigensolvers dropped the column for `k = 1` and a one-step fit errored.
   Keep the dimension, apply the sibling degree and dimension checks, and
   exercise both solver branches plus the public zero- and one-step routes.
8. [x] **Respect the perplexity-stepping budget and preserve the handoff.**
   With `max_iter = 12`, preliminary budgets of 1, 5, and 6 ran 14, 10, and
   12 total iterations with early stopping disabled. An explicit
   `Y_init_sdev` also reapplied at every stage. Allocate the remainder,
   scale only the first start, and report a clear error when a positive
   budget cannot fund the schedule. Verify exact budgets and transferred
   coordinates, including the final stage and the empty-schedule path.
9. [x] **Handle zero-gradient NDBD and the HSSNE Gaussian limit.** Scaled NDBD
   divides by zero at a stationary start. Leave zero gradients unscaled and
   preserve optimizer momentum. HSSNE at `alpha = 0` underflows for separated
   points. On that fixture SSNE's probability floor gave 30.68 versus 2.25e6
   for a direct log-space Gaussian KL calculation. Validate HSSNE against a
   log-space reference and specify the role of any clipping.
10. [x] **Remove the duplicate dense distance buffer.** The two distance
    kernels each allocate an extra `8N^2` bytes, or 200 MB at 5,000
    observations. Follow `tweight_cpp()`'s main-thread allocation and
    disjoint worker writes. Test values, shape, diagonal, and one versus
    multiple threads. The review did not measure runtime.

## Deferred proposals

- Keep failed calibration row IDs and a termination reason as follow-ups to
  the repairs. Report achieved perplexity only with a defined normalization
  and stage. FFT already exposes `repulsion` and `fft_backend` through
  `ret_extra`. Revisit a wider resolved-settings API when needed.
- Profile the repeated Barnes–Hut sparse buffer copies before changing them.
  The FFT attraction path offers an existing buffer access pattern.
- Consider optional `vizier` and `Rfast` only if installation friction
  warrants changing default plotting behavior or guarding geodesic MDS.
  A base partial sort might replace `Rfast::nth()` in `ballmmds`.
- A short NEWS entry and tag make sense for a deliberate release checkpoint.
  Keep the existing cost object, dense reference implementations, and
  historical articles. File splitting and broader infrastructure are outside
  this cleanup. The three exports, shared repeated-start input preparation,
  final-coordinate scoring, and native thread ownership are useful boundaries.
  Extend the existing mathematical and public-route tests for each repair.

## Verification and recovery

Use focused testthat files from the package directory, then the full suite and
configured Air and lintr checks. For native changes, build a copy of
`smallvis/` in a temporary directory so compiled objects stay out of the
checkout. Regenerate Rcpp wrappers if registered signatures change. Record
completed checks and open decisions here before moving to the next item.

## Completed checks

- **Calibration.** The near-tie public regression failed against the original
  installation and passed against the rebuilt package for dense R and C++,
  saved neighbors with dense output, and saved neighbors with sparse output.
  Compare against a per-row entropy root found independently with
  `uniroot()`. Conditional paths use shifted weights. Optional raw diagnostics
  recover their row scale (and can report zero raw weights), while raw
  weighting modes retain their fallback when no raw weight is representable.
  The `tsneu` and `pstsne` variants reuse calibrated affinities for their
  row-normalized kernels and retain their separate smooth-kNN normalization.
  On R 4.6.2 the targeted tests passed 999 assertions, and the full suite
  passed 2,183 assertions with zero failures, warnings, or skips. Rcpp
  regeneration was idempotent. Air, clang-format on the changed native
  source, and lintr 3.4.0 passed.
- **Alpha-beta divergences.** ABSNE and ABSSNE share the inclusive threshold
  dispatch. ABSSNE now evaluates the analytic limits and differentiates the
  Gaussian normalizer. The independent-reference checks failed on the
  original build and passed after the fix. They exercise equal distributions,
  forward and reverse KL, other exponents in the zero-parameter limits, and an
  interior case. The focused suite passed 953 assertions and the full suite
  passed 2,203 assertions with zero failures, warnings, or skips. Air,
  clang-format on the changed native source, and lintr 3.4.0 passed.
- **Unnormalized t-SNE.** L-BFGS is rejected before calibration. A normalized
  fit still runs L-BFGS, and a gradient-only unnormalized fit matches an
  independent one-step update with the original attraction scaling. The
  optimizer, dispatch, and FFT t-SNE tests passed 405 assertions. The
  regenerated function help explains how to interpret the reported KL value.
- **FFT grid.** The longest axis extent determines the shared spacing, and
  coordinate offsets and force moments use separate axis origins. The
  translated unit-circle test preserved the grid dimensions and one-step
  coordinates, with force error below 0.002% against direct evaluation.
  The FFT and native Barnes–Hut tests passed 274 assertions, including the
  existing rectangular-padding and wide-grid cases.
- **Graph symmetrization.** Dense HLSNE and CETSNE match the t-SNE
  reference affinities for both symmetric and mutual weighting and reject
  invalid options. Input equivalence, gradient, and API tests passed
  615 assertions.
- **Help and input guide.** Regenerated `smallvis.Rd` and executed the three
  corrected UMAP expressions extracted from that help. All returned finite
  coordinates. Rendered the inputs guide against the installed package,
  including its new `dist` and `"spca"` example.
- **Spectral initialization.** Both the dense and RSpectra branches retain
  the matrix dimension and solve the expected normalized eigenproblem for
  a graph with unequal degrees. Public zero- and one-step fits produce
  finite `N × 1` matrices. Spectral, API, and dispatch tests passed
  511 assertions.
- **Perplexity stepping.** The callback trace records 2, 2, 1, and 7
  iterations for a five-iteration preliminary budget and twelve total
  iterations. The preceding stage's last coordinates match the next
  stage's initial coordinates, including the final run. The wrapper
  rejects fractional or insufficient budgets. The stepping, repeat,
  and saved-neighbor tests passed 257 assertions. Regenerated the stepping
  help with the budget and initialization rules.
- **NDBD and HSSNE.** Scaled NDBD leaves a zero gradient unscaled and applies
  its usual momentum recurrence. Exact HSSNE at `alpha = 0` computes the
  Gaussian score from unfloored log probabilities. Independent log-space
  scores and directional derivatives match in both the R and native
  distance routes, even when individual output probabilities underflow.
  Optimizer, diagnostic, gradient, FFT SNE, and API tests passed
  700 assertions.
- **Native distances.** Both workers fill disjoint columns of the protected
  R output directly. Tests matched ordinary and squared distances against
  base R for five-dimensional data with duplicate points and verified
  identical results with one or three threads. The focused native,
  diagnostic, and input-equivalence tests passed 226 assertions.
- **Integrated validation.** The local suite and the suite run by
  `R CMD check` each passed 2,268 assertions with zero failures, warnings,
  or skips. On macOS with R 4.6.2, a fresh source build followed by
  `R CMD check --no-manual --no-build-vignettes` reported `Status: OK`.
  The check used `_R_CHECK_FORCE_SUGGESTS_=false` because the local libraries
  lack `fftwtools` and `zufft`. FFT tests exercised the default `stats`
  backend. Air, R Markdown formatting, lintr 3.4.0, generated Rd checks,
  and strict C++17 compilation and formatting checks on the changed native
  files passed.
